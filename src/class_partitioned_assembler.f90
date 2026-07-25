#include "scalar.fpp"

!=====================================================================!
! The partitioned system: a spatial assembler whose vectors are
! DISTRIBUTED across coarray images. The solver is untouched - it
! calls the same inner_product and jacobian-vector product as in
! serial - but after setup_partition a vector here means this image's
! owned slab, living in the graph's local frame:
!
!    global picture                each image after setup
!    ┌──────────────────┐          ┌───────────┐   halo   ┌───────────┐
!    │  all n dofs,     │   ──▶    │ own | gh  │ <------> │ own | gh  │
!    │  every image     │          └───────────┘          └───────────┘
!    │  a photocopy     │          vectors are the owned prefix only;
!    └──────────────────┘          the ghost tail exists just inside
!                                  the product's exchange buffer
!
! The three distributed questions, and what crosses the wire:
!
!    inner_product   owned slabs dot locally ... ONE SCALAR crosses
!    product         exchange the halo, then the local block's rows
!                    dot their edges ........... ONE VALUE PER CUT EDGE
!    residual        owned rows of b + skew - A x (the skew still
!                    reads the whole picture - the stage-4 tail)
!
! The words moved per product ARE the edge cut: partition quality and
! communication volume are the same number, and rcb has been
! minimizing both all along.
!
! block_preconditioner shrank to almost nothing here: the frame did
! the plumbing, so additive Schwarz is just each image's local block
! solve - no gather, no scatter, no collective.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_partitioned_assembler

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix
  use class_mesh             , only : mesh
  use class_stored_graph     , only : stored_graph
  use class_assembler        , only : spatial_assembler => assembler
  use interface_linear_solver, only : preconditioner
  use module_solve_mode      , only : FORWARD, REVERSE, WHOLE, &
       &                              is_valid_mode, is_valid_part

  implicit none

  private
  public :: partitioned_assembler
  public :: block_preconditioner

  ! the exchange wires. module variables, not object state: the
  ! standard forbids a coarray component on an allocatable object,
  ! and every driver holds its assembler allocatable. nothing is
  ! smuggled through them - each is written and read only inside its
  ! exchange, between two syncs.
  !   post  - owners post the owned slab, ghosts PULL   (forward)
  !   gpost - ghosts post their contributions, owners PULL and sum
  !                                                     (reverse)
  real(dp), allocatable :: post(:)[:]
  real(dp), allocatable :: gpost(:)[:]

  !===================================================================!
  ! The partitioned system
  !===================================================================!

  type, extends(spatial_assembler) :: partitioned_assembler

     ! set by setup_partition (after boundary conditions are applied);
     ! from then on every solver-facing vector is the owned slab
     logical :: partitioned = .false.

     ! the FACE frame:  [ owned dofs | face-ghost dofs ]
     !                    1 .. nown    nown+1 .. nloc
     ! the product's halo - one cell across each cut face
     integer              :: nown = 0, ngh = 0, nloc = 0
     integer, allocatable :: own(:)         ! global dofs of the owned prefix
     integer, allocatable :: gh_owner(:)    ! ghost j's owning image
     integer, allocatable :: gh_slot(:)     ! ...and its slot in that owner's frame

     ! the REVERSE book: the transpose scatters owned-row values into
     ! frame columns, and the ghost columns belong to other images -
     ! so each owned row COLLECTS contributions from the images that
     ! ghost it. grouped by owned slot: contributor image + its ghost
     ! index there. the transpose of the forward pull, zero messages
     ! to build (the partition is replicated).
     integer, allocatable :: rev_ptr(:)     ! (nown+1) owned slot -> its contributors
     integer, allocatable :: rev_img(:)     ! contributor image
     integer, allocatable :: rev_slot(:)    ! ...and its ghost index there

     ! the NODE frame:  [ owned dofs | node-ghost dofs ]
     !                    1 .. nown    nown+1 .. nwide
     ! the skew's halo - cells sharing a mesh POINT reach wider than
     ! the face halo (see node_graph on the mesh)
     integer              :: nwide = 0
     integer, allocatable :: wide_owner(:)  ! node-ghost j's owning image
     integer, allocatable :: wide_slot(:)   ! ...and its slot in that owner's frame
     integer, allocatable :: wide_dofs(:)   ! the wide frame's global dofs (own ++ node-ghosts)

     ! the owned rows, read in the frame (columns 1..nloc); the source's
     ! owned slab (state-independent, cached); and the lumped mass -
     ! diag(cell volume) - the owned slab of dR/dudot, purely local
     type(csr_matrix)      :: A_local
     real(dp), allocatable :: b_loc(:)
     real(dp), allocatable :: m_own(:)

   contains

     procedure :: setup_partition

     ! the distributed system queries
     procedure :: inner_product
     procedure :: get_jacobian_residual_product
     procedure :: state_residual

     ! the frozen seats, halo-aware (the base loops all cells with
     ! global dof indexing - impossible on a slab, so we override)
     procedure :: add_residual
     procedure :: add_jacobian_vector_product
     procedure :: add_jacobian_vector_product_transpose

     ! the adjoint's design accumulation, owned cells only
     procedure :: add_design_residual_transpose_product

     ! the transpose is verified whole-only in parallel (operator
     ! parts are a serial-only diagnostic; the reverse halo is genuine)
     procedure :: verify_transpose_consistency

     ! the door between the frame and the world outside it
     procedure :: replicate

     procedure, private :: exchange_halo
     procedure, private :: exchange_halo_reverse
     procedure, private :: spatial_local
     procedure, private :: spatial_local_transpose

  end type partitioned_assembler

  interface partitioned_assembler
     module procedure construct_partitioned
  end interface partitioned_assembler

  !===================================================================!
  ! Per-image block preconditioner (additive Schwarz without overlap).
  ! The residual arriving here IS this image's owned slab - the frame
  ! already did the plumbing - so the apply is the local block solve,
  ! nothing else.
  !===================================================================!

  type, extends(preconditioner) :: block_preconditioner

     class(preconditioner), allocatable :: block   ! per-image, owned-block sized

   contains

     procedure :: apply

  end type block_preconditioner

  interface block_preconditioner
     module procedure construct_block_preconditioner
  end interface block_preconditioner

contains

  !===================================================================!
  ! Construct from a mesh, exactly like the serial spatial assembler.
  ! Apply boundary conditions and the equation as usual, then call
  ! setup_partition.
  !===================================================================!

  impure type(partitioned_assembler) function construct_partitioned(grid) result(this)

    class(mesh), intent(in) :: grid

    this % spatial_assembler = spatial_assembler(grid)

  end function construct_partitioned

  !===================================================================!
  ! Partition the system's graph across the images and move into the
  ! local frame:
  !
  !    1. rcb stamps the parts (deterministic, replicated - every
  !       image computes the identical partition)
  !    2. the frame: owned dofs, then ghosts; and the address book -
  !       each ghost's owner and its slot in the owner's frame,
  !       computed locally from the replicated lists, no messages
  !    3. the owned rows of the assembled operator, far ends
  !       renumbered into the frame (local_block)
  !    4. the source's owned slab, cached
  !    5. the exchange wire, sized to the largest owned slab
  !    6. the system reports its OWNED length - from here on,
  !       vectors are slabs
  !===================================================================!

  impure subroutine setup_partition(this)

    class(partitioned_assembler), intent(inout) :: this

    type(csr_matrix)      :: A_global
    integer , allocatable :: gh_dofs(:), floc(:)
    real(dp), allocatable :: bfull(:)
    integer :: me, j, g, v, p, maxown

    me = this_image()

    ! the grid IS the graph - partition it in place
    call this % grid % partition_rcb(this % grid % cell_centers, num_images())

    ! ---- the frame ----
    this % own  = this % grid % dofs_of(this % grid % owned(me))
    gh_dofs     = this % grid % dofs_of(this % grid % ghosts(me))
    this % nown = size(this % own)
    this % ngh  = size(gh_dofs)
    this % nloc = this % nown + this % ngh
    floc        = this % grid % frame_inverse(me)

    ! ---- the exchange tables come from the graph: it holds the
    ! replicated partition, so it can answer where every ghost lives
    ! with no communication. the assembler only keeps the wires. ----
    call this % grid % ghost_owners(me, this % gh_owner, this % gh_slot)

    ! ---- the WIDE (node-ring) halo: the same partition, read on the
    ! mesh's node-adjacency graph (cells sharing a mesh point). its
    ! ghosts reach wider than the face halo - exactly the cells the
    ! skew correction interpolates through. the node graph answers its
    ! own table; the slots agree with the face frames because both
    ! graphs sort the SAME part stamps into the SAME owned lists. ----
    build_wide: block
      type(stored_graph) :: ng
      integer :: vw
      ng = this % grid % node_graph()
      ng % num_variables = this % grid % num_variables   ! same dof arithmetic
      call ng % set_partition([(this % grid % part_of(vw), vw = 1, this % grid % num_vertices)])
      this % wide_dofs = [this % own, this % grid % dofs_of(ng % ghosts(me))]
      this % nwide     = size(this % wide_dofs)
      call ng % ghost_owners(me, this % wide_owner, this % wide_slot)
    end block build_wide

    ! ---- the owned rows, in the frame ----
    call this % get_operator_csr(A_global)
    this % A_local = A_global % local_block(this % own, floc, this % nloc)

    ! ---- the source's owned slab (state-independent) ----
    allocate(bfull(this % grid % num_dofs()))
    call this % get_source(bfull)
    this % b_loc = bfull(this % own)

    ! ---- the lumped mass on the owned rows: each dof's cell volume.
    ! the mass is diagonal, so dR/dudot needs no halo at all ----
    allocate(this % m_own(this % nown))
    do j = 1, this % nown
       g = this % own(j)
       v = (g - 1)/this % grid % num_variables + 1
       this % m_own(j) = this % grid % cell_volumes(v)
    end do

    ! ---- the wire: coarray allocation must agree across images, so
    ! size it to the largest owned slab of any part ----
    maxown = 0
    do p = 1, this % grid % nparts
       maxown = max(maxown, (this % grid % own_ptr(p+1) - this % grid % own_ptr(p)) &
            &               * this % grid % num_variables)
    end do
    if (allocated(post)) deallocate(post)
    allocate(post(maxown)[*])

    ! ---- the reverse table comes from the graph too: who keeps
    ! copies of my owned dofs, and where in their ghost lists - the
    ! transpose exchange sends along exactly these entries ----
    call this % grid % ghost_copies(me, this % rev_ptr, this % rev_img, this % rev_slot)

    ! ---- the reverse wire: sized to the largest ghost count of any
    ! part (computable locally - the lists are replicated) ----
    reverse_wire: block
      integer :: img, maxgh
      maxgh = 0
      do img = 1, this % grid % nparts
         maxgh = max(maxgh, &
              & (this % grid % gh_ptr(img+1) - this % grid % gh_ptr(img)) &
              & * this % grid % num_variables)
      end do
      if (allocated(gpost)) deallocate(gpost)
      allocate(gpost(max(maxgh, 1))[*])
    end block reverse_wire

    ! ---- the vectors shrink: the state length becomes the owned slab,
    ! and the initial-condition field follows it (a transient march
    ! seeds U(:,1) from phi, so phi must be the owned slab too) ----
    this % phi = this % phi(this % own)
    this % num_state_vars = this % nown

    this % partitioned = .true.

  end subroutine setup_partition

  !===================================================================!
  ! The halo exchange: one slab out, one frame in. Generic over WHICH
  ! halo - hand it a book (owner list, slot list) and it fills that
  ! ring. The wire always posts the owned slab, so the face halo and
  ! the wider node halo share it; only the book differs.
  !
  !    post(1:nown) = my owned values          every image posts its
  !            sync                            slab, then pulls its
  !    buf(nown+j) = post(slot_j)[owner_j]     ghosts straight from
  !            sync                            the owners' frames
  !
  !    ┌─ image 1 ──┐          ┌─ image 2 ──┐
  !    │ own: a b c │ ──a───▶  │ gh:  a     │    each arrow is one
  !    │ gh:  x     │ ◀───x──  │ own: x y z │    cut edge's value -
  !    └────────────┘          └────────────┘    the traffic IS the cut
  !===================================================================!

  subroutine exchange_halo(this, v, buf, owner, slot)

    class(partitioned_assembler), intent(in)  :: this
    real(dp)                    , intent(in)  :: v(:)        ! owned slab
    real(dp)                    , intent(out) :: buf(:)      ! frame length
    integer                     , intent(in)  :: owner(:)    ! ghost -> its owner image
    integer                     , intent(in)  :: slot(:)     ! ghost -> its slot there

    integer :: j

    buf(1:this % nown)  = v(1:this % nown)
    post(1:this % nown) = v(1:this % nown)
    sync all
    do j = 1, size(owner)
       buf(this % nown + j) = post(slot(j))[owner(j)]
    end do
    sync all

  end subroutine exchange_halo

  !===================================================================!
  ! Distributed inner product: the owned slabs dot locally and one
  ! scalar crosses the wire. Each dof is owned exactly once, so the
  ! sum of the slabs' dots is the whole graph's dot, identically on
  ! every image.
  !===================================================================!

  impure real(dp) function inner_product(this, a, b)

    class(partitioned_assembler), intent(in) :: this
    real(dp)                    , intent(in) :: a(:)
    real(dp)                    , intent(in) :: b(:)

    inner_product = dot_product(a, b)
    if (this % partitioned) call co_sum(inner_product)

  end function inner_product

  !===================================================================!
  ! Distributed product: exchange the halo, then the local block's
  ! rows dot their edges -
  !
  !    v (owned slab) ──exchange──▶ [ v | ghosts ] ──A_local──▶ w
  !                                                  (owned slab)
  !
  ! No whole-vector collective anywhere: the input's ghosts are
  ! pulled fresh from their owners every product, so the only
  ! invariant a vector must keep is a valid owned slab - which
  ! elementwise arithmetic cannot break.
  !===================================================================!

  impure subroutine get_jacobian_residual_product(this, w, v, mode, part)

    class(partitioned_assembler), intent(in)           :: this
    real(dp)                    , intent(out)          :: w(:)
    real(dp)                    , intent(in)           :: v(:)
    integer                     , intent(in), optional :: mode
    integer                     , intent(in), optional :: part

    integer :: dir, sub

    dir = FORWARD
    if (present(mode)) dir = mode
    sub = WHOLE
    if (present(part)) sub = part

    ! a wrong tag dies at the door with its name
    if (.not. is_valid_mode(dir)) then
       write(*,'(1x,a,i0)') "partitioned_assembler: invalid mode tag ", dir
       error stop "partitioned_assembler: mode must be FORWARD or REVERSE"
    end if
    if (.not. is_valid_part(sub)) then
       write(*,'(1x,a,i0)') "partitioned_assembler: invalid part tag ", sub
       error stop "partitioned_assembler: part must be WHOLE, DIAGONAL, " // &
            & "LOWER_TRIANGLE or UPPER_TRIANGLE"
    end if

    if (this % partitioned) then
       if (sub .ne. WHOLE) then
          error stop "partitioned_assembler: distributed operator parts are a tracked deferral"
       end if

       ! do we need the spatial operator's TRANSPOSE? a plain REVERSE
       ! asks for it; a frozen transpose freeze marched forward also
       ! does - and a REVERSE of a transposed freeze is forward again
       ! (the same XOR the serial seat composes).
       spatial: block
         logical :: trans
         real(dp), allocatable :: Av(:)
         trans = (dir .eq. REVERSE)
         if (allocated(this % lin_coeff)) trans = this % lin_transpose .neqv. (dir .eq. REVERSE)

         if (trans) then
            Av = this % spatial_local_transpose(v)
         else
            Av = this % spatial_local(v)
         end if

         ! the frozen linearization J = beta*M - alpha*A rides the SAME
         ! halo as the steady operator - the mass is diagonal (so M = Mᵀ,
         ! no halo), only the spatial term crosses the wire:
         !
         !    J v = beta*(m_own * v_own)  -  alpha*(A[ᵀ] v)
         !          └── diagonal, local ─┘     └─ halo (fwd or rev) ─┘
         if (allocated(this % lin_coeff)) then
            w = real(this % lin_coeff(2), dp)*this % m_own*v(1:this % nown) &
                 & - real(this % lin_coeff(1), dp)*Av
         else
            w = Av
         end if
       end block spatial
       return
    end if

    ! a frozen linearization has no REPLICATED path either - refuse
    ! loudly rather than march the wrong operator (before setup only)
    if (allocated(this % lin_coeff)) then
       error stop "partitioned_assembler: a frozen linearization needs setup_partition first"
    end if

    ! before setup: the plain replicated composition
    if (dir .eq. REVERSE) then
       call this % transpose_product(w, v, sub)
    else if (sub .eq. WHOLE) then
       call this % get_jacobian_vector_product(w, v)
    else
       call this % get_jacobian_vector_product(w, v, filter = sub)
    end if

  end subroutine get_jacobian_residual_product

  !===================================================================!
  ! The steady spatial product on the owned rows: exchange the face
  ! halo, then the local block dots its edges. The one place the
  ! spatial operator crosses the wire - every frozen product and the
  ! nonlinear residual route through it.
  !
  !    v (owned slab) ──exchange──▶ [ v | ghosts ] ──A_local──▶ Av
  !                                                  (owned slab)
  !===================================================================!

  impure function spatial_local(this, v) result(Av)

    class(partitioned_assembler), intent(in) :: this
    real(dp)                    , intent(in) :: v(:)

    real(dp), allocatable :: Av(:), buf(:)

    allocate(Av(this % nown), buf(this % nloc))
    call this % exchange_halo(v, buf, this % gh_owner, this % gh_slot)
    call this % A_local % matvec(buf, Av)

  end function spatial_local

  !===================================================================!
  ! The reverse exchange: the transpose scatters owned-row values
  ! into frame columns, and the ghost tail belongs to other images.
  ! So each owner PULLS the contributions aimed at its rows and sums
  ! them - the mirror of the forward pull:
  !
  !    forward:   me  ◀──reads owned── owner        (owner → ghost)
  !    reverse:   owner ◀──reads ghost── me          (ghost → owner)
  !               and each owner SUMS what lands on its rows
  !
  !    ┌─ image 1 ──┐          ┌─ image 2 ──┐
  !    │ own: +=g   │ ◀──g───  │ gh: g      │   image 2 posts its ghost
  !    │            │          │            │   contribution; image 1
  !    └────────────┘          └────────────┘   pulls and adds it in
  !===================================================================!

  subroutine exchange_halo_reverse(this, yframe)

    class(partitioned_assembler), intent(in)    :: this
    real(dp)                    , intent(inout) :: yframe(:)   ! frame length; own rows updated

    integer :: s, k

    ! post my ghost-column contributions, tagged by ghost index
    gpost(1:this % ngh) = yframe(this % nown + 1 : this % nloc)
    sync all
    ! each owned row sums the contributions the ghosting images posted
    do s = 1, this % nown
       do k = this % rev_ptr(s), this % rev_ptr(s+1) - 1
          yframe(s) = yframe(s) + gpost(this % rev_slot(k))[this % rev_img(k)]
       end do
    end do
    sync all

  end subroutine exchange_halo_reverse

  !===================================================================!
  ! The transpose spatial product on the owned rows: A_local walked
  ! against its arrows (matvec_transpose) fills a frame vector, then
  ! the reverse exchange sums the ghost-column contributions back to
  ! their owners.
  !
  !    v (owned slab) ──A_localᵀ──▶ [ own | ghost contribs ]
  !                                  ── reverse exchange ──▶ owned slab
  !===================================================================!

  impure function spatial_local_transpose(this, v) result(Atv)

    class(partitioned_assembler), intent(in) :: this
    real(dp)                    , intent(in) :: v(:)

    real(dp), allocatable :: Atv(:), yframe(:)

    allocate(yframe(this % nloc))
    call this % A_local % matvec_transpose(v(1:this % nown), yframe)
    call this % exchange_halo_reverse(yframe)
    Atv = yframe(1:this % nown)

  end function spatial_local_transpose

  !===================================================================!
  ! The semi-discrete residual on the owned rows, halo-aware:
  !
  !    R = M*udot  -  A*u  +  b        (S(:,1)=u, S(:,2)=udot)
  !        │           │        │
  !        diagonal    spatial, halo  cached owned source
  !        m_own*udot  spatial_local
  !
  ! The base loops every cell with global dof indexing - impossible
  ! on a slab - so the frame answers it here. Newton reads its
  ! convergence residual through this seat, and so does the frozen
  ! linearized residual.
  !===================================================================!

  impure subroutine add_residual(this, residual, filter)

    class(partitioned_assembler), intent(in)           :: this
    type(scalar)                , intent(inout)        :: residual(:)
    integer                     , intent(in), optional :: filter

    if (present(filter)) then
       error stop "partitioned_assembler: a filtered residual (an operator part) " // &
            & "is a tracked deferral"
    end if

    residual(1:this % nown) = residual(1:this % nown) &
         & + this % m_own * this % S(1:this % nown, 2) &
         & - this % spatial_local(this % S(1:this % nown, 1)) &
         & + this % b_loc

  end subroutine add_residual

  !===================================================================!
  ! The jacobian-vector product on the owned rows, halo-aware:
  !
  !    pdt += [ beta*M - alpha*A ] vec        scalars = [alpha, beta]
  !            │          │
  !            diagonal   spatial, halo (spatial_local)
  !
  ! This is the frozen operator's action; the base loops all cells,
  ! so the frame answers it. The linearized residual b - Jx routes
  ! its J here.
  !===================================================================!

  impure subroutine add_jacobian_vector_product(this, pdt, vec, scalars, filter)

    class(partitioned_assembler), intent(in)           :: this
    type(scalar)                , intent(inout)        :: pdt(:)
    type(scalar)                , intent(in)           :: vec(:)
    type(scalar)                , intent(in)           :: scalars(:)
    integer                     , intent(in), optional :: filter

    if (present(filter)) then
       error stop "partitioned_assembler: a filtered product (an operator part) " // &
            & "is a tracked deferral"
    end if

    pdt(1:this % nown) = pdt(1:this % nown) &
         & + real(scalars(2), dp)*this % m_own*real(vec(1:this % nown), dp) &
         & - real(scalars(1), dp)*this % spatial_local(real(vec(1:this % nown), dp))

  end subroutine add_jacobian_vector_product

  !===================================================================!
  ! The transpose jacobian-vector product on the owned rows:
  !
  !    pdt += [ beta*Mᵀ - alpha*Aᵀ ] vec = [ beta*M - alpha*Aᵀ ] vec
  !            │           │
  !            diagonal    spatial transpose (reverse halo)
  !
  ! The adjoint marches a transposed freeze FORWARD, so the linearized
  ! residual b - Jᵀx routes its Jᵀ here. The mass is diagonal, so only
  ! the spatial transpose crosses the (reverse) wire.
  !===================================================================!

  impure subroutine add_jacobian_vector_product_transpose(this, pdt, vec, scalars, filter)

    class(partitioned_assembler), intent(in)           :: this
    type(scalar)                , intent(inout)        :: pdt(:)
    type(scalar)                , intent(in)           :: vec(:)
    type(scalar)                , intent(in)           :: scalars(:)
    integer                     , intent(in), optional :: filter

    if (present(filter)) then
       error stop "partitioned_assembler: a filtered transpose (an operator part) " // &
            & "is a tracked deferral"
    end if

    ! a pure mass action (alpha = 0, the bdf's step couplings) is
    ! diagonal - no wire at all. every image holds the same scalars,
    ! so every image takes the same branch and the syncs stay lined up.
    if (real(scalars(1), dp) .eq. 0.0_dp) then
       pdt(1:this % nown) = pdt(1:this % nown) &
            & + real(scalars(2), dp)*this % m_own*real(vec(1:this % nown), dp)
       return
    end if

    pdt(1:this % nown) = pdt(1:this % nown) &
         & + real(scalars(2), dp)*this % m_own*real(vec(1:this % nown), dp) &
         & - real(scalars(1), dp)*this % spatial_local_transpose(real(vec(1:this % nown), dp))

  end subroutine add_jacobian_vector_product_transpose

  !===================================================================!
  ! The adjoint's design accumulation, on the slab:
  !
  !    state scratch:  S(:,1) slab ──face exchange──▶ [own | ghosts]
  !                    ── dropped in at the frame's global dofs ──▶ x
  !    psi scratch:    zeros everywhere except my owned rows
  !
  !    my owned cells ──design_residual_rows──▶ my rows' share of
  !                                             psi' dR/dx
  !    one small sum over the design space  ──▶ dfdx (identical on
  !                                             every image)
  !
  ! Each row is owned exactly once, so the images' shares add up to
  ! the serial answer. The reduction carries ndv numbers - the size
  ! of the design space, not the mesh.
  !===================================================================!

  impure subroutine add_design_residual_transpose_product(this, dfdx, psi)

    class(partitioned_assembler), intent(inout) :: this
    real(dp)                    , intent(inout) :: dfdx(:)
    type(scalar)                , intent(in)    :: psi(:)

    real(dp)    , allocatable :: buf(:), xfull(:), delta(:)
    type(scalar), allocatable :: psifull(:)
    integer :: n

    n = this % grid % num_dofs()
    allocate(buf(this % nloc), xfull(n), psifull(n), delta(size(dfdx)))

    call this % exchange_halo(this % S(:,1), buf, this % gh_owner, this % gh_slot)
    xfull = 0.0_dp
    xfull(this % grid % frame(this_image())) = buf

    psifull = 0.0_dp
    psifull(this % own) = psi(1:this % nown)

    delta = 0.0_dp
    call this % design_residual_rows(delta, psifull, xfull, &
         & this % grid % owned(this_image()))
    call co_sum(delta)
    dfdx = dfdx + delta

  end subroutine add_design_residual_transpose_product

  !===================================================================!
  ! Verify the transpose whole-only: the reverse halo is a GENUINE
  ! transpose, not a symmetry claim, so the identity <w, A v> =
  ! <Aᵀ w, v> is a real self-check of the reverse exchange. Operator
  ! parts (the triangles) are a serial-only diagnostic - the colored
  ! smoothers replaced the triangle sweeps - so they are not checked
  ! here. Returns the relative defect over the owned rows, reduced.
  !===================================================================!

  impure real(dp) function verify_transpose_consistency(this) result(defect)

    class(partitioned_assembler), intent(in) :: this

    real(dp), allocatable :: v(:), w(:), Av(:), Atw(:)
    real(dp) :: lhs, rhs, scale
    integer  :: i

    allocate(v(this % nown), w(this % nown))
    do i = 1, this % nown
       v(i) = sin(real(this % own(i), dp)*0.7_dp) + 0.3_dp
       w(i) = cos(real(this % own(i), dp)*0.5_dp) + 0.2_dp
    end do

    Av  = this % spatial_local(v)
    Atw = this % spatial_local_transpose(w)

    lhs   = this % inner_product(w, Av)
    rhs   = this % inner_product(Atw, v)
    scale = max(abs(lhs), abs(rhs), 1.0_dp)
    defect = abs(lhs - rhs)/scale

  end function verify_transpose_consistency

  !===================================================================!
  ! The steady residual, owned rows only:
  !
  !    r  =  b     +  s(own)  -  A x
  !          cached   skew       the local block's product,
  !          slab     term       halo exchanged inside
  !
  ! The skew interpolates VERTEX values, whose ring of cells reaches
  ! wider than the face halo - so it rides the NODE halo, not the
  ! face one. Exchange the wide ring, drop it into a scratch by its
  ! global dofs, and let the skew read locally:
  !
  !     x (own slab) ──wide exchange──▶ [ own | node-ghosts ]
  !                                            │ scatter by global dof
  !                                            ▼
  !                          xwide = 0 everywhere except own+node-ghost
  !                                            │ get_skew_source
  !                                            ▼   (owned rows are exact:
  !                          keep sfull(own)       every point an owned
  !                                                cell touches has its
  !                                                whole ring in the halo)
  !
  ! No whole-vector co_sum: the wide exchange moves one value per node
  ! cut, not the entire field. The scratch xwide is still full length
  ! (get_skew_source loops all cells and indexes globally) - a memory
  ! O(n), not a collective; making the skew loop itself frame-local is
  ! a later refinement.
  !===================================================================!

  impure subroutine state_residual(this, r, x)

    class(partitioned_assembler), intent(in)  :: this
    real(dp)                    , intent(out) :: r(:)
    real(dp)                    , intent(in)  :: x(:)

    real(dp), allocatable :: Ax(:), xwide(:), buf(:), sfull(:)
    integer :: n

    if (.not. this % partitioned) then
       error stop "partitioned_assembler: a solve before setup_partition - " // &
            & "the frame does not exist yet"
    end if

    allocate(Ax(this % nown))
    call this % get_jacobian_residual_product(Ax, x)

    ! the node halo into a zeroed scratch, by global dof - no collective
    n = this % grid % num_dofs()
    allocate(buf(this % nwide), xwide(n), sfull(n))
    call this % exchange_halo(x, buf, this % wide_owner, this % wide_slot)
    xwide = 0.0_dp
    xwide(this % wide_dofs) = buf
    call this % get_skew_source(sfull, xwide)

    r = this % b_loc + sfull(this % own) - Ax

  end subroutine state_residual

  !===================================================================!
  ! The door: replicate an owned slab into a full global vector -
  ! scatter my slab into zeros, sum the images' contributions (each
  ! dof owned exactly once, so the sum IS the assembly):
  !
  !    image 1:  [ a b . . . ]      +
  !    image 2:  [ . . c d . ]      +      =   [ a b c d e ]
  !    image 3:  [ . . . . e ]                  everywhere
  !
  ! For the world outside the frame - writers, comparisons, the
  ! wide-reaching skew. The solve's hot path never calls it.
  !===================================================================!

  impure subroutine replicate(this, xloc, xfull)

    class(partitioned_assembler), intent(in)  :: this
    real(dp)                    , intent(in)  :: xloc(:)
    real(dp)                    , intent(out) :: xfull(:)

    xfull = 0.0_dp
    call this % grid % scatter(this_image(), xloc, xfull)
    call co_sum(xfull)

  end subroutine replicate

  !===================================================================!
  ! Block preconditioner constructor: just the per-image block,
  ! sized to the owned slab. No graph, no lists - the frame already
  ! did the plumbing.
  !===================================================================!

  impure type(block_preconditioner) function construct_block_preconditioner(block) &
       & result(this)

    class(preconditioner), intent(in) :: block

    allocate(this % block, source = block)

  end function construct_block_preconditioner

  !===================================================================!
  ! z = M^-1 r, additive over the images: r IS this image's owned
  ! slab, so apply the local block and be done. The images' slabs
  ! never meet - the next product's exchange carries whatever a
  ! neighbour needs to see.
  !===================================================================!

  subroutine apply(this, r, z)

    class(block_preconditioner), intent(in)  :: this
    real(dp)                   , intent(in)  :: r(:)
    real(dp)                   , intent(out) :: z(:)

    call this % block % apply(r, z)

  end subroutine apply

end module class_partitioned_assembler
