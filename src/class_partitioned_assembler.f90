!=====================================================================!
! The partitioned system: a spatial assembler whose vectors are
! distributed across coarray images. The solver is untouched - it calls
! the same inner_product and jacobian-vector product as in serial, and
! parallelism is provided here, on the system side:
!
!   inner_product - each image sums its owned entries, then reduces
!                   across images (co_sum)
!   product       - each image computes only its owned rows of A v,
!                   then the result is assembled across images
!
! The partition comes from the assembler's own graph (recursive
! coordinate bisection over the cell centroids); each image's owned dof
! list is derived from the owned cells. Every image runs the identical
! deterministic iteration, so the collective calls always match up.
! On one image both operations reduce exactly to their serial forms.
!
! The module also provides block_preconditioner: a per-image
! preconditioner (e.g. an algebraic multigrid built on the owned
! diagonal block) applied to the owned entries and assembled across
! images - additive Schwarz without overlap. It keeps the vectors
! consistent across images, so it composes with the solver's
! pre_conditioner slot like any other preconditioner.
!
! Operator parts (diagonal/triangles) fall back to the replicated
! serial implementation - correct on every image, just not distributed.
! REVERSE products route through the one transpose seat, so the
! override, the declared-symmetric claim, and the refusal behave
! exactly as on the serial system. The distributed fast path covers the
! forward whole-operator product, which is the kernel's inner loop.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_partitioned_assembler

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix
  use class_mesh             , only : mesh
  use class_graph            , only : mesh_graph
  use class_assembler        , only : spatial_assembler => assembler
  use interface_linear_solver, only : preconditioner
  use module_solve_mode      , only : FORWARD, REVERSE, WHOLE, &
       &                              is_valid_mode, is_valid_part

  implicit none

  private
  public :: partitioned_assembler
  public :: block_preconditioner

  !===================================================================!
  ! The partitioned system
  !===================================================================!

  type, extends(spatial_assembler) :: partitioned_assembler

     ! set by setup_partition (after boundary conditions are applied)
     logical              :: partitioned = .false.
     integer, allocatable :: own(:)         ! this image's owned dofs
     type(csr_matrix)     :: A              ! assembled operator for the owned rows

   contains

     procedure :: setup_partition
     procedure :: owned_dofs

     ! the distributed system queries
     procedure :: inner_product
     procedure :: get_jacobian_residual_product

  end type partitioned_assembler

  interface partitioned_assembler
     module procedure construct_partitioned
  end interface partitioned_assembler

  !===================================================================!
  ! Per-image block preconditioner (additive Schwarz without overlap):
  ! applies a preconditioner built on this image's owned block to the
  ! owned entries of the residual, then assembles across images.
  !===================================================================!

  type, extends(preconditioner) :: block_preconditioner

     class(preconditioner), allocatable :: block   ! per-image, owned-block sized
     integer, allocatable :: own(:)

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
  ! Partition the system's graph across the images (recursive
  ! coordinate bisection), derive this image's owned dof list, and
  ! assemble the operator for the owned-row products. Deterministic and
  ! replicated, so every image computes the identical partition.
  !===================================================================!

  impure subroutine setup_partition(this)

    class(partitioned_assembler), intent(inout) :: this

    type(mesh_graph) :: gp

    gp = this % g
    call gp % partition_rcb(this % grid % cell_centers, num_images())

    this % own = this % owned_dofs(gp, this_image())

    call this % get_operator_csr(this % A)

    this % partitioned = .true.

  end subroutine setup_partition

  !===================================================================!
  ! The dofs of the cells part k owns (variable-fastest interleaving)
  !===================================================================!

  pure function owned_dofs(this, gp, k) result(dofs)

    class(partitioned_assembler), intent(in) :: this
    type(mesh_graph)            , intent(in) :: gp
    integer                     , intent(in) :: k

    integer, allocatable :: cells(:), dofs(:)
    integer :: i, ivar, pos

    cells = gp % owned(k)
    allocate(dofs(size(cells) * gp % num_variables))

    pos = 0
    do i = 1, size(cells)
       do ivar = 1, gp % num_variables
          pos = pos + 1
          dofs(pos) = gp % dof(cells(i), ivar)
       end do
    end do

  end function owned_dofs

  !===================================================================!
  ! Distributed inner product: sum over this image's owned entries,
  ! reduced across images. Each dof is owned exactly once, so the
  ! result equals the serial dot product, identically on every image.
  !===================================================================!

  impure real(dp) function inner_product(this, a, b)

    class(partitioned_assembler), intent(in) :: this
    real(dp)                    , intent(in) :: a(:)
    real(dp)                    , intent(in) :: b(:)

    integer :: i

    if (.not. this % partitioned) then
       inner_product = dot_product(a, b)
       return
    end if

    inner_product = 0.0_dp
    do i = 1, size(this % own)
       inner_product = inner_product + a(this % own(i))*b(this % own(i))
    end do
    call co_sum(inner_product)

  end function inner_product

  !===================================================================!
  ! Distributed product: this image computes only its owned rows of
  ! A v, and the full result is assembled across images. Requires v to
  ! be consistent across images, which the solver maintains because
  ! every vector it forms descends from assembled products and the
  ! replicated residual. Parts and REVERSE fall back to the replicated
  ! serial composition.
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

    ! a wrong tag dies at the door with its name - this override is
    ! reached by dynamic dispatch, so it carries its own door (the
    ! reserved stored-CSR fast path below is untouched)
    if (.not. is_valid_mode(dir)) then
       write(*,'(1x,a,i0)') "partitioned_assembler: invalid mode tag ", dir
       error stop "partitioned_assembler: mode must be FORWARD or REVERSE"
    end if
    if (.not. is_valid_part(sub)) then
       write(*,'(1x,a,i0)') "partitioned_assembler: invalid part tag ", sub
       error stop "partitioned_assembler: part must be WHOLE, DIAGONAL, " // &
            & "LOWER_TRIANGLE or UPPER_TRIANGLE"
    end if

    if (this % partitioned .and. dir .eq. FORWARD .and. sub .eq. WHOLE) then
       w = 0.0_dp
       call this % A % matvec_rows(v, w, this % own)
       call co_sum(w)
       return
    end if

    ! the REVERSE fallback routes through the one transpose seat, so a
    ! genuine override, a declared-symmetric claim, and the refusal all
    ! behave exactly as on the serial system (replicated on every image)
    if (dir .eq. REVERSE) then
       call this % transpose_product(w, v, sub)
       return
    end if

    ! replicated serial forward fallback (correct on every image)
    if (sub .eq. WHOLE) then
       call this % get_jacobian_vector_product(w, v)
    else
       call this % get_jacobian_vector_product(w, v, filter = sub)
    end if

  end subroutine get_jacobian_residual_product

  !===================================================================!
  ! Block preconditioner constructor: the per-image preconditioner
  ! (sized to the owned block) and the owned dof list it applies to.
  !===================================================================!

  impure type(block_preconditioner) function construct_block_preconditioner(block, own) &
       & result(this)

    class(preconditioner), intent(in) :: block
    integer              , intent(in) :: own(:)

    allocate(this % block, source = block)
    this % own = own

  end function construct_block_preconditioner

  !===================================================================!
  ! z = M^-1 r, additive over the images: each image applies its block
  ! preconditioner to its owned entries, the rest of z is zero, and the
  ! assembled sum is the block-diagonal preconditioned residual -
  ! identical on every image.
  !===================================================================!

  subroutine apply(this, r, z)

    class(block_preconditioner), intent(in)  :: this
    real(dp)                   , intent(in)  :: r(:)
    real(dp)                   , intent(out) :: z(:)

    real(dp), allocatable :: r_loc(:), z_loc(:)
    integer :: i, nown

    nown = size(this % own)
    allocate(r_loc(nown), z_loc(nown))

    do i = 1, nown
       r_loc(i) = r(this % own(i))
    end do

    call this % block % apply(r_loc, z_loc)

    z = 0.0_dp
    do i = 1, nown
       z(this % own(i)) = z_loc(i)
    end do
    call co_sum(z)

  end subroutine apply

end module class_partitioned_assembler
