#include "scalar.fpp"

!=====================================================================!
! The mandelbrot law: the escape-time fractal recast as physics.
!
! The iteration z -> z^2 + c is forward euler with step one on
!
!     dz/dt = z^2 + c - z
!
! (write euler's step z + h*(z^2 + c - z) and set h = 1: the z
! cancels and the map remains - an identity, not an approximation).
! split z = (u, v) and c = (c1, c2): a two-variable REACTING source
! with no flux, and the first law in this codebase whose coupling
! graph has an edge -
!
!     (u) ───────── (v)         the complex square mixes the two
!                               parts, both ways:
!     dS_u/dv = -2v             the edge, written as numbers,
!     dS_v/du =  2v             lives in dsource_dq below
!
! Two laws in one type: julia (c held fixed, the cell is z0) and
! mandelbrot proper (the cell IS c, z0 = 0) - the c this source
! reads either comes from its stored constant or from the point
! st % x it is asked at.
!
! One sign to know about: the assembler's semi-discrete residual is
! M*udot - A*u + b = 0, so a fluxless law marches du/dt = -S. This
! source therefore returns MINUS the map's velocity,
!
!     S = z - z^2 - c
!
! and the marched map comes out z -> z^2 + c on the nose.
!=====================================================================!

module class_mandelbrot

  use iso_fortran_env  , only : dp => REAL64
  use interface_physics, only : source, point_state

  implicit none

  private
  public :: mandelbrot_source

  type, extends(source) :: mandelbrot_source

     real(dp) :: c(2)   = 0.0_dp    ! the constant (julia mode)
     logical  :: c_is_x = .false.   ! mandelbrot mode: c = the point

   contains

     procedure :: value
     procedure :: dsource_dq

  end type mandelbrot_source

  interface mandelbrot_source
     module procedure create
  end interface mandelbrot_source

contains

  !===================================================================!
  ! julia:       mandelbrot_source(c)   z0 = the cell (prime phi)
  ! mandelbrot:  mandelbrot_source()    c  = the cell, z0 = 0
  !===================================================================!

  pure type(mandelbrot_source) function create(c) result(this)

    real(dp), intent(in), optional :: c(2)

    if (present(c)) then
       this % c      = c
       this % c_is_x = .false.
    else
       this % c_is_x = .true.
    end if

    ! two variables (Re z, Im z), one coupling edge: the square
    ! mixes them - the first non-empty coupling graph in the house
    call this % form_coupling(2, tails = [1], heads = [2])

  end function create

  !===================================================================!
  ! S = z - z^2 - c  (minus the map's velocity; see the banner) -
  ! and ZERO outside the circle of no return:
  !
  !        │z│ <= 2 :  the map runs
  !        │z│ >  2 :  S = 0, the point holds still
  !
  ! once │z│ passes 2 the orbit provably never comes back (the
  ! classical escape bound, │c│ <= 2), so freezing it changes no
  ! escape count - and it keeps the arithmetic finite: an escaped
  ! orbit squares its way to overflow in a handful of steps if the
  ! law keeps pushing it.
  !===================================================================!

  pure function value(this, st) result(S)

    class(mandelbrot_source), intent(in) :: this
    type(point_state)       , intent(in) :: st

    type(scalar) :: S(this % num_vertices)
    real(dp)     :: c(2)

    if (st % q(1)*st % q(1) + st % q(2)*st % q(2) .gt. 4.0_dp) then
       S = 0.0_dp
       return
    end if

    c = this % c
    if (this % c_is_x) c = st % x(1:2)

    S(1) = st % q(1) - (st % q(1)*st % q(1) - st % q(2)*st % q(2)) - c(1)
    S(2) = st % q(2) -  2.0_dp*st % q(1)*st % q(2)                 - c(2)

  end function value

  !===================================================================!
  ! The coupling edge as numbers (the jacobian of S):
  !
  !     dS/dq = [ 1 - 2u     2v   ]      and ZERO in the frozen
  !             [  -2v     1 - 2u ]      zone, like S itself - which
  !                                      is what will keep the
  !                                      adjoint finite there
  !===================================================================!

  pure function dsource_dq(this, st) result(dS)

    class(mandelbrot_source), intent(in) :: this
    type(point_state)       , intent(in) :: st

    type(scalar) :: dS(this % num_vertices, this % num_vertices)

    if (st % q(1)*st % q(1) + st % q(2)*st % q(2) .gt. 4.0_dp) then
       dS = 0.0_dp
       return
    end if

    dS(1,1) = 1.0_dp - 2.0_dp*st % q(1)
    dS(1,2) =          2.0_dp*st % q(2)
    dS(2,1) =        - 2.0_dp*st % q(2)
    dS(2,2) = 1.0_dp - 2.0_dp*st % q(1)

  end function dsource_dq

end module class_mandelbrot
