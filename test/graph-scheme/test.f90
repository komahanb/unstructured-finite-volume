!=====================================================================!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_algebra

  use util_precision
  use operation_expression

  implicit none

  type(expression) :: W11, W21, W22, B1, B2, delta1, delta2
  type(expression) :: h, C(7)

  W11    = unknown(1)
  W21    = unknown(2)
  W22    = unknown(3)
  B1     = unknown(4)
  B2     = unknown(5)
  delta1 = unknown(6)
  delta2 = unknown(7)

  h = design() ! supplied physical step size

  C(1) = W11 - delta1
  C(2) = W21 + W22 - delta2
  C(3) = B1 + B2 - h
  C(4) = B1*delta1 + B2*delta2 - h**2/2.0_dp
  C(5) = B1*delta1**2 + B2*delta2**2 - h**3/3.0_dp
  C(6) = B1*W11*delta1 + B2*(W21*delta1 + W22*delta2) - h**3/6.0_dp
  C(7) = W11 - W22

end program test_graph_algebra
