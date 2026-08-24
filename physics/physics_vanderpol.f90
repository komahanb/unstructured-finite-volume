!=====================================================================!
! The van der Pol oscillator as a governing constraint, and an energy
! as a functional integrand.
!
! The equation, at degree N,
!
!      R  =  q^(N)  -  nu (1 - q^2) q^(N-1)  +  q  =  0
!
! is the ordinary oscillator at N = 2 and its higher-degree
! continuation above that: the damping always acts on the derivative
! one below the highest, and the restoring term always acts on the
! value. The functional integrand supplied beside it is the energy
!
!      F  =  ( q^2 + (q')^2 ) / 2 .
!
! Each is stated once, at one instant. Everything else is written in
! physics_integrand, which both extend.
!
!             THE PARTIALS OF THE RESIDUAL
!
!      dR/dq^(N)     =  1
!      dR/dq^(N-1)   =  -nu (1 - q^2)
!      dR/dq         =  2 nu q q^(N-1)  +  1
!      dR/dq^(d)     =  0     for every other d
!      dR/dnu        =  -(1 - q^2) q^(N-1)
!
! and at N > 2 the degrees between one and N-2 are exactly zero,
! which is what a partial indexed one place out would not be. They
! are not written out in code: the rule is differentiated as
! written.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_vanderpol

  use iso_fortran_env       , only : dp => REAL64
  use physics_integrand     , only : nodal_integrand
  use util_derivative_terms , only : derivative_terms, &
       & operator(+), operator(-), operator(*)

  implicit none

  private
  public :: van_der_pol, van_der_pol_energy

  type, extends(nodal_integrand) :: van_der_pol
   contains
     procedure :: name       => residual_name
     procedure :: at_instant => residual_at_instant
  end type van_der_pol

  type, extends(nodal_integrand) :: van_der_pol_energy
   contains
     procedure :: name       => energy_name
     procedure :: at_instant => energy_at_instant
  end type van_der_pol_energy

  interface van_der_pol
     module procedure create_residual
  end interface van_der_pol

  interface van_der_pol_energy
     module procedure create_energy
  end interface van_der_pol_energy

contains

  function create_residual(degree) result(this)

    integer, intent(in) :: degree
    type(van_der_pol) :: this

    call this % declare_degree(degree)

  end function create_residual

  !===================================================================!
  ! The energy reads a velocity, so a degree-zero problem has none.
  !===================================================================!

  function create_energy(degree) result(this)

    integer, intent(in) :: degree
    type(van_der_pol_energy) :: this

    if (degree < 1) then
       error stop 'physics_vanderpol: the energy needs a velocity'
    end if

    call this % declare_degree(degree)

  end function create_energy

  pure function residual_name(this) result(name)

    class(van_der_pol), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'van der pol residual'

  end function residual_name

  pure function energy_name(this) result(name)

    class(van_der_pol_energy), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'van der pol energy'

  end function energy_name

  !===================================================================!
  ! THE TWO RULES.
  !===================================================================!

  pure function residual_at_instant(this, q, nu) result(r)

    class(van_der_pol)    , intent(in) :: this
    type(derivative_terms), intent(in) :: q(0:)
    type(derivative_terms), intent(in) :: nu
    type(derivative_terms) :: r

    type(derivative_terms) :: one
    integer :: n

    n   = this % equation_degree()
    one = derivative_terms(1.0_dp, nu)

    r = q(n) - nu * (one - q(0) * q(0)) * q(n - 1) + q(0)

  end function residual_at_instant

  pure function energy_at_instant(this, q, nu) result(r)

    class(van_der_pol_energy), intent(in) :: this
    type(derivative_terms)   , intent(in) :: q(0:)
    type(derivative_terms)   , intent(in) :: nu
    type(derivative_terms) :: r

    associate (u1 => this, u2 => nu); end associate

    r = 0.5_dp * (q(0) * q(0) + q(1) * q(1))

  end function energy_at_instant

end module physics_vanderpol
