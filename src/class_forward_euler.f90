#include "scalar.fpp"

!=====================================================================!
! Forward euler: the plainest march there is.
!
!     ( 1 )──▶( 2 )──▶( 3 )──▶ ... ──▶( n )      the step chain at
!                                                power 1 - each step
!     u_new = u_old - h * R(u_old, 0) / m        needs only the one
!                                                before it
!
! R is the residual at the old state with the time derivative set to
! zero, and m is the lumped mass, so -R/m is the state's velocity;
! the step walks down it. Nothing is implicit and nothing is solved.
!
! Why it exists here: at h = 1 on a pure reaction law (no flux) the
! step IS the law's own map, exactly - marching the mandelbrot source
! this way makes the trajectory the textbook iteration z -> z^2 + c.
! Other marchers on the same law paint different fractals (an
! implicit step is itself a newton solve - newton-fractal country);
! this one paints the classical set.
!=====================================================================!

module class_forward_euler

  use iso_fortran_env     , only : dp => REAL64
  use interface_integrator, only : integrator
  use class_assembler     , only : assembler

  implicit none

  private
  public :: forward_euler

  type, extends(integrator) :: forward_euler

   contains

     procedure :: step
     procedure :: get_bandwidth
     procedure :: get_stencil_coeff

  end type forward_euler

  interface forward_euler
     module procedure create
  end interface forward_euler

contains

  !===================================================================!
  ! Construct over [tinit, tfinal] at step h. The scheme is explicit,
  ! so no solver rides along.
  !===================================================================!

  impure type(forward_euler) function create(system, tinit, tfinal, h) result(this)

    class(assembler), intent(in) :: system
    real(dp)        , intent(in) :: tinit, tfinal, h

    call this % construct(system, tinit, tfinal, h, implicit = .false.)

  end function create

  !===================================================================!
  ! Advance one step: the residual at the old state, with the time
  ! derivative zero, is minus the mass times the velocity, so
  !
  !     u_new  =  u_old - h * R / m
  !===================================================================!

  impure subroutine step(this, t, U, h, p, ierr)

    class(forward_euler), intent(inout) :: this
    real(dp)            , intent(inout) :: t(:)
    type(scalar)        , intent(inout) :: U(:,:,:)
    integer             , intent(in)    :: p
    real(dp)            , intent(in)    :: h
    integer             , intent(out)   :: ierr

    type(scalar), allocatable :: res(:)
    real(dp)    , allocatable :: m(:)
    integer                   :: n, kk

    ierr = 0
    kk   = size(U, dim = 1)
    n    = this % system % get_num_state_vars()

    t(kk) = t(kk-1) + h

    allocate(res(n), m(n))

    ! Evaluate the residual at the old state with no time derivative.
    this % system % S(:,1) = U(kk-1,:,1)
    this % system % S(:,2) = 0.0_dp
    res = 0.0_dp
    call this % system % add_residual(res)
    call this % system % get_lumped_mass(m)

    U(kk,:,1) = U(kk-1,:,1) - h*real(res, dp)/m
    U(kk,:,2) = -real(res, dp)/m           ! The velocity, for the record.

  end subroutine step

  !===================================================================!
  ! One step back is all this scheme ever reads.
  !===================================================================!

  pure integer function get_bandwidth(this, step_index) result(width)

    class(forward_euler), intent(in) :: this
    integer             , intent(in) :: step_index

    width = 1

  end function get_bandwidth

  !===================================================================!
  ! Return the two-point difference (u_k - u_{k-1})/h for anything
  ! that asks for the stencil.
  !===================================================================!

  impure subroutine get_stencil_coeff(this, p, h, scoeff)

    class(forward_euler), intent(in)               :: this
    integer             , intent(in)               :: p
    real(dp)            , intent(in)               :: h
    type(scalar)        , intent(out), allocatable :: scoeff(:)

    allocate(scoeff(p+1))
    scoeff(1) =  1.0_dp/h
    scoeff(2) = -1.0_dp/h

  end subroutine get_stencil_coeff

end module class_forward_euler
