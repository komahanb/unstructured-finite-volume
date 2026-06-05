#include "scalar.fpp"

!=====================================================================!
! Backward Difference Formula time integrator (orders 1-6). Extends the
! abstract integrator: get_bandwidth gives how many past steps the
! formula uses at a step, step predicts the new state from the BDF
! stencil and, when implicit, drives the residual to zero via newton.
! The linearization coefficients coeff(n+1) = (A(p,1)/h)^n are the
! [alpha, beta, gamma, ...] handed to the jacobian-vector product.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_bdf

  use iso_fortran_env     , only : dp => REAL64
  use interface_integrator, only : integrator
  use interface_assembler , only : assembler
  use class_newton_solver , only : newton

  implicit none

  private
  public :: bdf

  !-------------------------------------------------------------------!
  ! BDF integrator type
  !-------------------------------------------------------------------!

  type, extends(integrator) :: bdf

     integer                   :: max_order = 6
     type(scalar), allocatable :: A(:,:)        ! BDF coefficient table

   contains

     procedure :: step
     procedure :: get_bandwidth
     procedure :: get_linearization_coeff
     procedure :: get_stencil_coeff

  end type bdf

  interface bdf
     module procedure create
  end interface bdf

contains

  !===================================================================!
  ! Construct a BDF integrator of the requested accuracy order
  !===================================================================!

  type(bdf) function create(system, tinit, tfinal, h, max_order) result(this)

    class(assembler), intent(in) :: system
    real(dp)        , intent(in) :: tinit, tfinal, h
    integer         , intent(in) :: max_order

    call this % construct(system, tinit, tfinal, h, implicit = .true.)

    this % max_order = min(max_order, 6)

    ! Coefficient table - http://www.scholarpedia.org/article/Backward_differentiation_formulas
    allocate(this % A(this % max_order, this % max_order + 1))
    this % A = 0.0d0

    if (this % max_order .ge. 1) this % A(1,1:2) = [1.0d0, -1.0d0]
    if (this % max_order .ge. 2) this % A(2,1:3) = [3.0d0, -4.0d0, 1.0d0]/2.0d0
    if (this % max_order .ge. 3) this % A(3,1:4) = [11.0d0, -18.0d0, 9.0d0, -2.0d0]/6.0d0
    if (this % max_order .ge. 4) this % A(4,1:5) = [25.0d0, -48.0d0, 36.0d0, -16.0d0, 3.0d0]/12.0d0
    if (this % max_order .ge. 5) this % A(5,1:6) = [137.0d0, -300.0d0, 300.0d0, -200.0d0, 75.0d0, -12.0d0]/60.0d0
    if (this % max_order .ge. 6) this % A(6,1:7) = [147.0d0, -360.0d0, 450.0d0, -400.0d0, 225.0d0, -72.0d0, 10.0d0]/60.0d0

  end function create

  !===================================================================!
  ! Number of past steps the formula uses (ramps up to max_order)
  !===================================================================!

  pure integer function get_bandwidth(this, step_index) result(width)

    class(bdf), intent(in) :: this
    integer   , intent(in) :: step_index

    width = step_index - 1

    if (width .gt. this % max_order) width = this % max_order

  end function get_bandwidth

  !===================================================================!
  ! Linearization coefficients  coeff(n+1) = (A(p,1)/h)^n,  n = 0..order
  !===================================================================!

  subroutine get_linearization_coeff(this, p, h, coeff)

    class(bdf)  , intent(in)    :: this
    integer     , intent(in)    :: p
    real(dp)    , intent(in)    :: h
    type(scalar), intent(inout) :: coeff(:)

    integer :: n

    do n = 0, this % system % get_differential_order()
       coeff(n+1) = (this % A(p,1)/h)**n
    end do

  end subroutine get_linearization_coeff

  !===================================================================!
  ! First-derivative stencil at bandwidth p:  scoeff(i+1) = A(p,i+1)/h,
  ! i = 0..p - the same coefficients that form udot from past states in
  ! step(). The adjoint backward sweep uses the whole stencil.
  !===================================================================!

  subroutine get_stencil_coeff(this, p, h, scoeff)

    class(bdf)  , intent(in)               :: this
    integer     , intent(in)               :: p
    real(dp)    , intent(in)               :: h
    type(scalar), intent(out), allocatable :: scoeff(:)

    integer :: i

    allocate(scoeff(p+1))

    do i = 0, p
       scoeff(i+1) = this % A(p, i+1)/h
    end do

  end subroutine get_stencil_coeff

  !===================================================================!
  ! Advance one step: predict the new state from the BDF stencil, then
  ! (implicit) drive the residual to zero with a newton solve.
  !===================================================================!

  impure subroutine step(this, t, U, h, p, ierr)

    class(bdf)  , intent(inout) :: this
    real(dp)    , intent(inout) :: t(:)
    type(scalar), intent(inout) :: U(:,:,:)       ! (window, nvars, order+1)
    integer     , intent(in)    :: p
    real(dp)    , intent(in)    :: h
    integer     , intent(out)   :: ierr

    type(scalar), allocatable :: coeff(:)
    type(newton)              :: nlsolver
    integer                   :: kk, torder, n, i

    ierr   = 0
    kk     = size(U, dim=1)                       ! newest state at window end
    torder = this % system % get_differential_order()

    ! Advance the time
    t(kk) = t(kk-1) + h

    ! Predictor: carry the last value for the lowest order
    U(kk,:,1) = U(kk-1,:,1)

    ! Higher-order states from the BDF stencil
    do n = 1, torder

       U(kk,:,n+1) = 0.0d0

       do i = 0, p
          U(kk,:,n+1) = U(kk,:,n+1) + (this % A(p,i+1)/h)*U(kk-i,:,n)
       end do

    end do

    ! Implicit correction driving the residual to zero
    if (this % implicit) then

       allocate(coeff(torder+1))

       call this % get_linearization_coeff(p, h, coeff)
       call nlsolver % solve(this % system, coeff, U(kk,:,:))

       deallocate(coeff)

    end if

  end subroutine step

end module class_bdf
