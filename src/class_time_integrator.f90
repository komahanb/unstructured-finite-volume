!=====================================================================!
! Backward-euler time marching for the finite volume system. Each step
! solves the transient operator (M/dt - A) phi^{n+1} = M/dt phi^n - b
! that the assembler builds when transient is on; the linear solve is a
! conjugate gradient, and phi_old is advanced between steps.
!
! Higher-order bdf would carry more past states - this is bdf1.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_time_integrator

  use iso_fortran_env          , only : dp => REAL64
  use class_assembler          , only : assembler
  use class_conjugate_gradient , only : conjugate_gradient

  implicit none

  private
  public :: time_integrator

  type :: time_integrator

     type(conjugate_gradient) :: lsolver   ! inner linear solve each step
     real(dp)                 :: tinit
     real(dp)                 :: tfinal
     real(dp)                 :: dt

   contains

     procedure :: integrate

  end type time_integrator

  interface time_integrator
     module procedure create
  end interface time_integrator

contains

  type(time_integrator) function create(system, tinit, tfinal, dt, max_it, max_tol) &
       & result(this)

    type(assembler), intent(inout) :: system
    real(dp)       , intent(in)    :: tinit, tfinal, dt
    integer        , intent(in)    :: max_it
    real(dp)       , intent(in)    :: max_tol

    ! Turn on the transient operator and build the inner solver around it
    call system % set_transient(dt)
    this % tinit  = tinit
    this % tfinal = tfinal
    this % dt     = dt
    this % lsolver = conjugate_gradient(system, max_it, max_tol, 0)

  end function create

  !===================================================================!
  ! March from tinit to tfinal starting from phi0; return the final
  ! state in phi.
  !===================================================================!

  subroutine integrate(this, phi0, phi)

    class(time_integrator), intent(inout) :: this
    real(dp)              , intent(in)    :: phi0(:)
    real(dp), allocatable , intent(out)   :: phi(:)

    real(dp), allocatable :: xnew(:)
    integer  :: k, nsteps
    real(dp) :: t

    allocate(phi, source = phi0)
    nsteps = nint((this % tfinal - this % tinit)/this % dt)
    t = this % tinit

    write(*,'(1x,a,i0,a,es10.3)') "marching ", nsteps, " backward-euler steps, dt = ", this % dt

    do k = 1, nsteps
       this % lsolver % FVAssembler % phi_old = phi   ! previous state
       call this % lsolver % solve(xnew)              ! (M/dt - A) phi = M/dt phi_old - b
       phi = xnew
       t = t + this % dt
    end do

  end subroutine integrate

end module class_time_integrator
