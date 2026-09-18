!=====================================================================!
! THE TWO-STAGE, THIRD-ORDER SDIRK TABLEAU.
!
! Let Omega be a point. Let (W11, W21, W22, B1, B2, delta1, delta2)
! : Omega -> R^7 be the entries of the tableau, scaled by the step h,
!
!      A = | W11   0  |      b = (B1, B2)      c = (delta1, delta2)
!          | W21  W22 |
!
! satisfying the seven order conditions
!
!      W11 - delta1                                    = 0
!      W21 + W22 - delta2                              = 0
!      B1 + B2 - h                                     = 0
!      B1 delta1 + B2 delta2 - h^2 / 2                 = 0
!      B1 delta1^2 + B2 delta2^2 - h^3 / 3             = 0
!      B1 W11 delta1 + B2 (W21 delta1 + W22 delta2) - h^3 / 6 = 0
!      W11 - W22                                       = 0
!
! Omega has no coordinate, so Omega_h = Omega, no derivative is
! approximated, and the discrete residual R_h : R^7 -> R^7 is the
! seven conditions themselves. The discrete solution is the zero of
! R_h near the Crouzeix root, W11 = W22 = h (3 + sqrt 3) / 6.
!=====================================================================!

program dirk_coefficients

  use util_precision    , only : dp
  use operation_manifold, only : continuous_manifold, discrete_manifold
  use operation_field   , only : continuous_field, discrete_field
  use operation_field   , only : operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual
  use util_verbosity    , only : set_verbosity

  implicit none

  real(dp), parameter :: h = 0.1_dp
  character(len=6), parameter :: names(7) = &
       & [character(len=6) :: 'W11', 'W21', 'W22', 'B1', 'B2', 'delta1', 'delta2']

  type(continuous_manifold) :: omega
  type(discrete_manifold)   :: omega_h

  type(continuous_field)    :: W11, W21, W22, B1, B2, delta1, delta2, equation(7), crouzeix
  type(discrete_field)      :: estimate, solution

  type(continuous_residual) :: r
  type(discrete_residual)   :: r_h

  real(dp) :: values(7)
  integer  :: i

  ! Omega a point; Omega_h = Omega
  omega   = continuous_manifold()
  omega_h = omega % discretize()

  ! the unknown functions Omega -> R: the entries of the tableau
  W11    = omega % unknown('W11')
  W21    = omega % unknown('W21')
  W22    = omega % unknown('W22')
  B1     = omega % unknown('B1')
  B2     = omega % unknown('B2')
  delta1 = omega % unknown('delta1')
  delta2 = omega % unknown('delta2')

  ! the seven order conditions as functions Omega -> R of the unknowns
  equation(1) = W11 - delta1
  equation(2) = W21 + W22 - delta2
  equation(3) = B1 + B2 - h
  equation(4) = B1*delta1 + B2*delta2 - h**2/2.0_dp
  equation(5) = B1*delta1**2 + B2*delta2**2 - h**3/3.0_dp
  equation(6) = B1*W11*delta1 + B2*(W21*delta1 + W22*delta2) - h**3/6.0_dp
  equation(7) = W11 - W22

  ! R on Omega, and R_h = R on Omega_h
  r   = continuous_residual(omega, equation)
  r_h = r % discretize(omega_h)

  ! the initial estimate near the Crouzeix root, a field Omega -> R^7,
  ! and the zero of R_h from it
  crouzeix = continuous_field(omega, h * [0.8_dp, -0.6_dp, 0.8_dp, 0.5_dp, 0.5_dp, 0.8_dp, 0.2_dp])
  estimate = crouzeix % discretize(omega_h)
  ! the convergence of every Newton step and linear solve is printed
  ! by image 1
  if (this_image() == 1) call set_verbosity(1)
  call r_h % minimize(estimate, solution)

  call solution % values(values)
  print '(a6, " = ", es24.16)', (names(i), values(i), i = 1, 7)

end program dirk_coefficients
