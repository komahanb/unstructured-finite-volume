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
  use class_assembler        , only : spatial_assembler => assembler
  use interface_linear_solver, only : preconditioner
  use module_solve_mode      , only : FORWARD, REVERSE, WHOLE, &
       &                              is_valid_mode, is_valid_part

  implicit none

  private
  public :: partitioned_assembler
  public :: block_preconditioner

  ! the exchange wire. a module variable, not object state: the
  ! standard forbids a coarray component on an allocatable object,
  ! and every driver holds its assembler allocatable. nothing is
  ! smuggled through it - it is written and read only inside
  ! exchange_halo, between two syncs.
  real(dp), allocatable :: post(:)[:]

  !===================================================================!
  ! The partitioned system
  !===================================================================!

  type, extends(spatial_assembler) :: partitioned_assembler

     ! set by setup_partition (after boundary conditions are applied);
     ! from then on every solver-facing vector is the owned slab
     logical :: partitioned = .false.

     ! the frame:  [ owned dofs | ghost dofs ]
     !               1 .. nown    nown+1 .. nloc
     integer              :: nown = 0, ngh = 0, nloc = 0
     integer, allocatable :: own(:)         ! global dofs of the owned prefix
     integer, allocatable :: gh_owner(:)    ! ghost j's owning image
     integer, allocatable :: gh_slot(:)     ! ...and its slot in that owner's frame

     ! the owned rows, read in the frame (columns 1..nloc), and the
     ! source's owned slab (state-independent, cached once)
     type(csr_matrix)      :: A_local
     real(dp), allocatable :: b_loc(:)

   contains

     procedure :: setup_partition

     ! the distributed system queries
     procedure :: inner_product
     procedure :: get_jacobian_residual_product
     procedure :: state_residual

     procedure, private :: exchange_halo

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
    integer , allocatable :: gh_dofs(:), floc(:), owner_inv(:)
    real(dp), allocatable :: bfull(:)
    integer :: me, j, g, v, p, prev_p, maxown

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

    ! ---- the address book: owner and slot of every ghost, computed
    ! locally - the partition bookkeeping is replicated, so no
    ! communication at all ----
    allocate(this % gh_owner(this % ngh), this % gh_slot(this % ngh))
    prev_p = 0
    do j = 1, this % ngh
       g = gh_dofs(j)
       v = (g - 1)/this % grid % num_variables + 1     ! the dof's vertex
       p = this % grid % part_of(v)
       if (p .ne. prev_p) then
          owner_inv = this % grid % frame_inverse(p)
          prev_p    = p
       end if
       this % gh_owner(j) = p
       this % gh_slot(j)  = owner_inv(g)
    end do

    ! ---- the owned rows, in the frame ----
    call this % get_operator_csr(A_global)
    this % A_local = A_global % local_block(this % own, floc, this % nloc)

    ! ---- the source's owned slab (state-independent) ----
    allocate(bfull(this % grid % num_dofs()))
    call this % get_source(bfull)
    this % b_loc = bfull(this % own)

    ! ---- the wire: coarray allocation must agree across images, so
    ! size it to the largest owned slab of any part ----
    maxown = 0
    do p = 1, this % grid % nparts
       maxown = max(maxown, (this % grid % own_ptr(p+1) - this % grid % own_ptr(p)) &
            &               * this % grid % num_variables)
    end do
    if (allocated(post)) deallocate(post)
    allocate(post(maxown)[*])

    ! ---- the vectors shrink ----
    this % num_state_vars = this % nown

    this % partitioned = .true.

  end subroutine setup_partition

  !===================================================================!
  ! The halo exchange: one slab out, one slab in.
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

  subroutine exchange_halo(this, v, buf)

    class(partitioned_assembler), intent(in)  :: this
    real(dp)                    , intent(in)  :: v(:)     ! owned slab
    real(dp)                    , intent(out) :: buf(:)   ! frame length

    integer :: j

    buf(1:this % nown)  = v(1:this % nown)
    post(1:this % nown) = v(1:this % nown)
    sync all
    do j = 1, this % ngh
       buf(this % nown + j) = post(this % gh_slot(j))[this % gh_owner(j)]
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

    real(dp), allocatable :: buf(:)
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

    ! a frozen linearization has no distributed path yet - refuse
    ! loudly rather than march the wrong operator
    if (allocated(this % lin_coeff)) then
       error stop "partitioned_assembler: a frozen linearization is a tracked deferral"
    end if

    if (this % partitioned) then
       if (dir .eq. REVERSE) then
          error stop "partitioned_assembler: a distributed transpose march is a tracked deferral"
       end if
       if (sub .ne. WHOLE) then
          error stop "partitioned_assembler: distributed operator parts are a tracked deferral"
       end if
       allocate(buf(this % nloc))
       call this % exchange_halo(v, buf)
       call this % A_local % matvec(buf, w)
       return
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
  ! The steady residual, owned rows only:
  !
  !    r  =  b     +  s(own)  -  A x
  !          cached   skew       the local block's product,
  !          slab     term       halo exchanged inside
  !
  ! The skew correction still reads the whole solution (its flux
  ! loops have not learned the frame yet), so ONE whole-vector
  ! assembly survives here, per outer pass, off the hot path - the
  ! stage-4 tail of the dereplication plan.
  !===================================================================!

  impure subroutine state_residual(this, r, x)

    class(partitioned_assembler), intent(in)  :: this
    real(dp)                    , intent(out) :: r(:)
    real(dp)                    , intent(in)  :: x(:)

    real(dp), allocatable :: Ax(:), xfull(:), sfull(:)
    integer :: n

    if (.not. this % partitioned) then
       error stop "partitioned_assembler: a solve before setup_partition - " // &
            & "the frame does not exist yet"
    end if

    allocate(Ax(this % nown))
    call this % get_jacobian_residual_product(Ax, x)

    n = this % grid % num_dofs()
    allocate(xfull(n), sfull(n))
    xfull = 0.0_dp
    call this % grid % scatter(this_image(), x, xfull)
    call co_sum(xfull)
    call this % get_skew_source(sfull, xfull)

    r = this % b_loc + sfull(this % own) - Ax

  end subroutine state_residual

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
