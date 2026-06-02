!=====================================================================!
! A boundary condition on a named physical group.
!
! Everything is a robin condition  a*phi + b*dphi/dn = c  on the
! boundary face. dirichlet and neumann just fall out:
!
!   dirichlet (phi = g)        : a=1, b=0, c=g
!   neumann   (dphi/dn = g)    : a=0, b=1, c=g
!
! With a one sided gradient  dphi/dn ~ (phi_b - phi_p)/delta  the face
! contributes a diffusive flux  area*(phi_b - phi_p)/delta  to cell p.
! Eliminating the (unknown) face value phi_b leaves a piece that depends
! on phi_p (goes to the diagonal) and a constant (goes to the rhs):
!
!   flux = lhs_coeff(area,delta)*phi_p  -  rhs_coeff(area,delta)
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_boundary_condition

  use iso_fortran_env, only : dp => REAL64
  use class_string   , only : string

  implicit none

  private
  public :: boundary_condition
  public :: BC_DIRICHLET, BC_NEUMANN, BC_ROBIN
  public :: dirichlet, neumann, robin

  ! Kinds (for printing/identity; the maths is always robin)
  integer, parameter :: BC_DIRICHLET = 1
  integer, parameter :: BC_NEUMANN   = 2
  integer, parameter :: BC_ROBIN     = 3

  type :: boundary_condition

     type(string) :: name             ! physical group name
     integer      :: tag  = -1        ! resolved tag number
     integer      :: kind = BC_DIRICHLET
     real(dp)     :: a = 1.0_dp        ! a*phi + b*dphi/dn = c
     real(dp)     :: b = 0.0_dp
     real(dp)     :: c = 0.0_dp        ! default: homogeneous dirichlet phi=0

   contains

     procedure :: lhs_coeff           ! diagonal contribution (multiplies phi_p)
     procedure :: rhs_coeff           ! constant contribution to the rhs
     procedure :: print

  end type boundary_condition

contains

  !===================================================================!
  ! Named constructors - set the (a,b,c) triad and the kind tag
  !===================================================================!

  pure type(boundary_condition) function dirichlet(value) result(bc)
    real(dp), intent(in) :: value
    bc % kind = BC_DIRICHLET
    bc % a = 1.0_dp; bc % b = 0.0_dp; bc % c = value
  end function dirichlet

  pure type(boundary_condition) function neumann(flux) result(bc)
    real(dp), intent(in) :: flux
    bc % kind = BC_NEUMANN
    bc % a = 0.0_dp; bc % b = 1.0_dp; bc % c = flux
  end function neumann

  pure type(boundary_condition) function robin(a, b, c) result(bc)
    real(dp), intent(in) :: a, b, c
    bc % kind = BC_ROBIN
    bc % a = a; bc % b = b; bc % c = c
  end function robin

  !===================================================================!
  ! Coefficient on phi_p contributed to the diagonal of the operator
  !===================================================================!

  pure real(dp) function lhs_coeff(this, area, delta)
    class(boundary_condition), intent(in) :: this
    real(dp), intent(in) :: area, delta
    real(dp) :: denom
    denom = this % a + this % b/delta
    lhs_coeff = - area*this % a/(delta*denom)
  end function lhs_coeff

  !===================================================================!
  ! Constant the face contributes to the right hand side (already
  ! carries the sign for moving to the rhs)
  !===================================================================!

  pure real(dp) function rhs_coeff(this, area, delta)
    class(boundary_condition), intent(in) :: this
    real(dp), intent(in) :: area, delta
    real(dp) :: denom
    denom = this % a + this % b/delta
    rhs_coeff = - area*this % c/(delta*denom)
  end function rhs_coeff

  !===================================================================!
  ! Pretty print
  !===================================================================!

  subroutine print(this)
    class(boundary_condition), intent(in) :: this
    character(len=9) :: kindname
    select case (this % kind)
    case (BC_DIRICHLET); kindname = "dirichlet"
    case (BC_NEUMANN);   kindname = "neumann"
    case (BC_ROBIN);     kindname = "robin"
    case default;        kindname = "unknown"
    end select
    write(*,'(1x,a,1x,a,1x,a,i0,1x,3(a,es12.4))') &
         & "bc", trim(kindname), "tag=", this % tag, &
         & "a=", this % a, " b=", this % b, " c=", this % c
  end subroutine print

end module class_boundary_condition
