#include "scalar.fpp"

!=====================================================================!
! Test fixture: an assembled operator wrapped as a system, so kernel-
! level suites (krylov, advection) can drive solvers through the same
! contract production code uses - the residual, the jacobian-vector
! product, and the inner product. The transpose is genuine
! (matvec_transpose on the stored operator), so REVERSE products are
! legitimate on non-symmetric operators without a symmetry claim.
!
! Operator parts (diagonal/triangles) are refused with a clear error:
! the kernels this fixture serves run on whole-operator products only.
!=====================================================================!

module class_csr_system

  use iso_fortran_env    , only : dp => REAL64
  use class_csr          , only : csr_matrix
  use interface_assembler, only : assembler, WHOLE
  use module_solve_mode  , only : FORWARD, REVERSE, is_valid_mode

  implicit none

  private
  public :: csr_system

  type, extends(assembler) :: csr_system

     type(csr_matrix)      :: A
     real(dp), allocatable :: b(:)

   contains

     ! the system contract on the stored operator
     procedure :: state_residual                => csr_residual
     procedure :: get_jacobian_residual_product => csr_product
     procedure :: transpose_product             => csr_transpose_product

     ! the deferred assembler contract
     procedure :: add_residual                => csr_add_residual
     procedure :: add_jacobian_vector_product => csr_add_product
     procedure :: add_initial_condition       => csr_add_initial_condition

  end type csr_system

  interface csr_system
     module procedure create
  end interface csr_system

contains

  !===================================================================!
  ! Wrap an assembled operator (and optionally its right-hand side)
  !===================================================================!

  impure type(csr_system) function create(A, b) result(this)

    type(csr_matrix), intent(in)           :: A
    real(dp)        , intent(in), optional :: b(:)

    this % A = A
    this % num_state_vars     = A % nrows
    this % differential_order = 0

    allocate(this % b(A % nrows))
    this % b = 0.0_dp
    if (present(b)) this % b = b

  end function create

  !===================================================================!
  ! The residual of the stored system at x: r = b - A x
  !===================================================================!

  impure subroutine csr_residual(this, r, x)

    class(csr_system), intent(in)  :: this
    real(dp)         , intent(out) :: r(:)
    real(dp)         , intent(in)  :: x(:)

    real(dp), allocatable :: w(:)

    allocate(w(this % num_state_vars))
    call this % A % matvec(x, w)
    r = this % b - w

  end subroutine csr_residual

  !===================================================================!
  ! The product on the stored operator: forward and genuine transpose.
  ! Parts are refused - the kernels this fixture serves use the whole
  ! operator only.
  !===================================================================!

  impure subroutine csr_product(this, w, v, mode, part)

    class(csr_system), intent(in)           :: this
    real(dp)         , intent(out)          :: w(:)
    real(dp)         , intent(in)           :: v(:)
    integer          , intent(in), optional :: mode
    integer          , intent(in), optional :: part

    integer :: dir, sub

    dir = FORWARD
    if (present(mode)) dir = mode
    sub = WHOLE
    if (present(part)) sub = part

    ! a wrong tag dies at the door - dynamic dispatch bypasses the base
    ! door, so this override carries its own
    if (.not. is_valid_mode(dir)) then
       write(*,'(1x,a,i0)') "csr_system: invalid mode tag ", dir
       error stop "csr_system: mode must be FORWARD or REVERSE"
    end if

    if (sub .ne. WHOLE) then
       error stop "csr_system: operator parts are not provided - " // &
            & "this fixture serves whole-operator kernels"
    end if

    if (dir .eq. REVERSE) then
       call this % A % matvec_transpose(v, w)
    else
       call this % A % matvec(v, w)
    end if

  end subroutine csr_product

  !===================================================================!
  ! The genuine transpose on the stored operator (no symmetry claim)
  !===================================================================!

  impure subroutine csr_transpose_product(this, w, v, sub)

    class(csr_system), intent(in)  :: this
    real(dp)         , intent(out) :: w(:)
    real(dp)         , intent(in)  :: v(:)
    integer          , intent(in)  :: sub

    if (sub .ne. WHOLE) then
       error stop "csr_system: operator parts are not provided - " // &
            & "this fixture serves whole-operator kernels"
    end if

    call this % A % matvec_transpose(v, w)

  end subroutine csr_transpose_product

  !===================================================================!
  ! The steady residual of the stored system: res += A S(:,1) - b
  !===================================================================!

  impure subroutine csr_add_residual(this, residual, filter)

    class(csr_system), intent(in)           :: this
    type(scalar)     , intent(inout)        :: residual(:)
    type(integer)    , intent(in), optional :: filter

    real(dp), allocatable :: w(:), x(:)

    allocate(w(this % num_state_vars), x(this % num_state_vars))
    x = real(this % S(:,1), dp)
    call this % A % matvec(x, w)
    residual = residual + (w - this % b)

  end subroutine csr_add_residual

  !===================================================================!
  ! pdt += scalars(1) * A vec on the stored operator
  !===================================================================!

  impure subroutine csr_add_product(this, pdt, vec, scalars, filter)

    class(csr_system), intent(in)           :: this
    type(scalar)     , intent(inout)        :: pdt(:)
    type(scalar)     , intent(in)           :: vec(:)
    type(scalar)     , intent(in)           :: scalars(:)
    type(integer)    , intent(in), optional :: filter

    real(dp), allocatable :: w(:), v(:)

    allocate(w(this % num_state_vars), v(this % num_state_vars))
    v = real(vec, dp)
    call this % A % matvec(v, w)
    pdt = pdt + scalars(1)*w

  end subroutine csr_add_product

  !===================================================================!
  ! The stored system starts from rest
  !===================================================================!

  impure subroutine csr_add_initial_condition(this, U)

    class(csr_system), intent(in)    :: this
    type(scalar)     , intent(inout) :: U(:,:)

    U = 0.0_dp

  end subroutine csr_add_initial_condition

end module class_csr_system
