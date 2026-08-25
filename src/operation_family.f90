!=====================================================================!
! The time-marching family: an edge function on a scheme's coupling
! that also declares how the scheme marches.
!
! Beyond the edge rule it inherits, a family declares how far back
! its widest constraint reads, how many stages one instant holds, and
! which derivative degree it solves for. The coefficients it produces
! are dimensionless: the weight an edge finally carries is the
! coefficient times a power of the step, and that power is fixed by
! the two degrees the edge joins, which operation_step_scaling
! supplies and operation_weight multiplies in.
!
!=====================================================================!
!
!                    THE TWO LAGRANGE FUNCTIONALS
!
! A multistep family interpolates through the instants behind k at
! scaled offsets
!
!      theta_j  =  (t_k - t_(k-j)) / dt_k ,     theta_0 = 0 ,
!
! and every coefficient is either the slope at zero of a Lagrange
! basis function through those nodes (a difference) or its integral
! over the last step (a quadrature). Both functionals are supplied
! here for the families that extend this type. On a uniform grid
! every theta_j is j and the tabulated coefficients come out; they
! are the uniform value of the formula, not a separate case.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family

  use iso_fortran_env         , only : dp => REAL64
  use operation_edge_function , only : edge_function
  use util_derivative_terms   , only : derivative_terms, value, &
       & operator(+), operator(-), operator(*), operator(/)

  implicit none

  private
  public :: family
  public :: offsets, slope_at_zero, integral_over_step

  type, abstract, extends(edge_function) :: family

   contains

     procedure(family_degree_interface) , deferred :: history_depth
     procedure(family_count_interface)  , deferred :: num_stages
     procedure(family_degree_interface) , deferred :: primary_degree
     procedure(family_pattern_interface), deferred :: row_pattern

  end type family

  abstract interface

     pure integer function family_count_interface(this)
       import :: family
       class(family), intent(in) :: this
     end function family_count_interface

     !----------------------------------------------------------------!
     ! Two counts that read the degree of the equation, because the
     ! rows a family makes depend on how many derivatives there are
     ! to determine: how far back the widest of those rows reaches,
     ! and which degree the governing constraint determines. The
     ! second is the value for a difference family and the highest
     ! degree for a quadrature or a stage family; every other degree
     ! is determined by a derived row.
     !----------------------------------------------------------------!

     pure integer function family_degree_interface(this, equation_degree)
       import :: family
       class(family), intent(in) :: this
       integer      , intent(in) :: equation_degree
     end function family_degree_interface

     !----------------------------------------------------------------!
     ! The sources of the derived row that determines one degree: how
     ! many instants back each one lies, and what degree it holds. A
     ! degree the family determines by no derived row gives an empty
     ! pattern rather than stopping, so a caller may ask about every
     ! degree.
     !----------------------------------------------------------------!

     pure subroutine family_pattern_interface(this, determines, equation_degree, &
          & offset, source_degree)
       import :: family
       class(family), intent(in) :: this
       integer      , intent(in) :: determines, equation_degree
       integer, allocatable, intent(out) :: offset(:), source_degree(:)
     end subroutine family_pattern_interface

  end interface

contains

  !===================================================================!
  ! THE LAGRANGE FUNCTIONALS.
  !
  ! The scaled offsets at instant k of the n instants k, k-1, ...,
  ! k-n+1: theta_j is the distance back to instant k - j in units of
  ! the step ending at k. Every step read must be positive; a zero
  ! or negative one stops the program, the ratio being undefined.
  !===================================================================!

  pure function offsets(dt, k, n) result(theta)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, n
    type(derivative_terms) :: theta(0:n-1)

    integer :: j

    do j = k - n + 2, k
       if (value(dt(j)) <= 0.0_dp) then
          error stop 'operation_family: every time step read is positive'
       end if
    end do

    theta(0) = derivative_terms(0.0_dp, dt(k))
    do j = 1, n - 1
       theta(j) = theta(j - 1) + dt(k - j + 1) / dt(k)
    end do

  end function offsets

  !===================================================================!
  ! The slope at zero of the j-th Lagrange basis function through the
  ! nodes u: the sum over the other nodes m of one over (u_j - u_m)
  ! times the product over the remaining nodes of (0 - u_i)/(u_j - u_i).
  !===================================================================!

  pure function slope_at_zero(u, j) result(s)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: s

    type(derivative_terms) :: term
    integer :: m, i, n

    n = ubound(u, 1)
    s = derivative_terms(0.0_dp, u(j))

    do m = 0, n
       if (m == j) cycle
       term = derivative_terms(1.0_dp, u(j)) / (u(j) - u(m))
       do i = 0, n
          if (i == j .or. i == m) cycle
          term = term * (-u(i)) / (u(j) - u(i))
       end do
       s = s + term
    end do

  end function slope_at_zero

  !===================================================================!
  ! The coefficients of the j-th Lagrange basis function through the
  ! nodes u as a polynomial in u, lowest power first: the product of
  ! (u - u_m)/(u_j - u_m) over the other nodes, multiplied out one
  ! factor at a time.
  !===================================================================!

  pure function basis_polynomial(u, j) result(c)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: c(0:ubound(u, 1))

    type(derivative_terms) :: shifted(0:ubound(u, 1))
    integer :: m, l, n

    n = ubound(u, 1)
    do l = 0, n
       c(l) = derivative_terms(0.0_dp, u(j))
    end do
    c(0) = derivative_terms(1.0_dp, u(j))

    do m = 0, n
       if (m == j) cycle
       shifted(0) = derivative_terms(0.0_dp, u(j))
       do l = 1, n
          shifted(l) = c(l - 1)
       end do
       do l = 0, n
          c(l) = (shifted(l) - u(m) * c(l)) / (u(j) - u(m))
       end do
    end do

  end function basis_polynomial

  !===================================================================!
  ! The integral from -1 to 0 of the j-th Lagrange basis function
  ! through the nodes u: each power u^n contributes (-1)^n / (n + 1).
  ! That is the last step in scaled units, so the result is the
  ! quadrature weight of node j over the step ending at zero.
  !===================================================================!

  pure function integral_over_step(u, j) result(w)

    type(derivative_terms), intent(in) :: u(0:)
    integer               , intent(in) :: j
    type(derivative_terms) :: w

    type(derivative_terms) :: c(0:ubound(u, 1))
    integer :: n

    c = basis_polynomial(u, j)
    w = derivative_terms(0.0_dp, u(j))

    do n = 0, ubound(u, 1)
       w = w + (real((-1)**n, dp) / real(n + 1, dp)) * c(n)
    end do

  end function integral_over_step
end module operation_family
