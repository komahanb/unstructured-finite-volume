!=====================================================================!
! The van der Pol oscillator as a governing constraint, and two
! functional integrands beside it, each stated once at one instant.
!
! The equation, at degree N,
!
!      R  =  q^(N)  -  nu (1 - q^2) q^(N-1)  +  q  =  0
!
! is the ordinary oscillator at N = 2 and its higher-degree
! continuation above that: the damping acts on the derivative one
! below the highest, and the restoring term on the value. The energy
! is
!
!      F  =  ( q^2 + (q')^2 ) / 2
!
! and the dissipation, the power the damping term draws,
!
!      F  =  nu (1 - q^2) (q')^2 ,
!
! a functional that reads the design itself, so its own partial in
! the design is not zero.
!
! Each is an expression over the unknown and the design, and its
! partials in either are taken by evaluating it; none is written out
! here. A functional reads the velocity, so a degree-zero problem has
! none: stating it at degree zero stops the program.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_vanderpol

  use util_precision    , only : dp
  use physics_expression, only : expression, unknown, design, derivative, stated, &
       & operator(+), operator(-), operator(*), operator(**)

  implicit none

  private
  public :: van_der_pol, van_der_pol_energy, van_der_pol_dissipation

contains

  function van_der_pol(degree) result(r)

    integer, intent(in) :: degree
    type(expression) :: r

    type(expression) :: q, nu

    q  = unknown()
    nu = design()

    r = stated(derivative(q, degree) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, degree - 1) &
         & + derivative(q, 0), degree, 'van der pol residual')

  end function van_der_pol

  function van_der_pol_energy(degree) result(f)

    integer, intent(in) :: degree
    type(expression) :: f

    type(expression) :: q

    q = unknown()

    f = stated(0.5_dp * (derivative(q, 0)**2 + derivative(q, 1)**2), degree, 'van der pol energy')

  end function van_der_pol_energy

  function van_der_pol_dissipation(degree) result(f)

    integer, intent(in) :: degree
    type(expression) :: f

    type(expression) :: q, nu

    q  = unknown()
    nu = design()

    f = stated(nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, 1) * derivative(q, 1), &
         & degree, 'van der pol dissipation')

  end function van_der_pol_dissipation

end module physics_vanderpol
