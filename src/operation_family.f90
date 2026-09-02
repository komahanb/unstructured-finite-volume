!=====================================================================!
! The time-marching family: an edge function on a scheme's coupling
! that also declares how the scheme marches.
!
! Beyond the edge rule it inherits, a family declares how far back
! its widest constraint reads, how many stages one instant contains, and
! which derivative degree it solves for. The coefficients it produces
! are dimensionless: the weight an edge finally stores is the
! coefficient times a power of the step, and that power is fixed by
! the two degrees the edge joins, which operation_weight multiplies
! in.
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
! every theta_j is j and the tabulated coefficients result; they
! are the uniform value of the formula, not a separate case.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family

  use util_precision  , only : dp
  use operation_edge_function , only : edge_function
  use util_derivative_terms   , only : derivative_terms, value, &
       & operator(+), operator(-), operator(*), operator(/)

  implicit none

  private
  public :: family
  public :: offsets, slope_at_zero, integral_over_step, negated

  type, abstract, extends(edge_function) :: family

   contains

     procedure :: history_depth  => family_history_depth
     procedure :: num_stages     => family_num_stages
     procedure :: stage_weight   => family_stage_weight
     procedure :: step_quadrature => family_step_quadrature
     procedure :: primary_degree => family_primary_degree
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
     ! many instants back each one lies, and what degree it has. A
     ! degree the family determines by no derived row gives an empty
     ! pattern rather than stopping, so a caller may query every
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
  !===================================================================!
  ! The defaults of a family unless it overrides them: one stage, one
  ! instant of history, and the primary unknown the equation's own
  ! highest derivative. bdf reads further back and solves for the value;
  ! dirk has as many stages as its tableau has weights.
  !===================================================================!

  pure integer function family_num_stages(this)

    class(family), intent(in) :: this

    associate (u1 => this); end associate
    family_num_stages = 1

  end function family_num_stages

  !===================================================================!
  ! THE QUADRATURE OVER ONE STEP, as the weights of the instants the
  ! step reaches back over: weight(j) belongs to the instant j - 1
  ! back from k, so a rule of m nodes returns m weights.
  !
  ! A rule on m nodes integrates the degree m - 1 interpolant through
  ! them exactly. Over a step of width h that leaves a local error of
  ! order h**(m+1), and over the T/h steps of the horizon an error of
  ! order h**m. So m nodes give order m.
  !
  ! WHAT A FAMILY RETURNS UNLESS IT OVERRIDES THIS: one node at
  ! weight one, which is the rectangle rule and first order. A family
  ! whose stencil already contains the instants of an interpolatory
  ! rule returns that rule instead, and a stage family is never
  ! called - its quadrature is the tableau, read through stage_weight.
  !
  ! MATCHING THE RULE TO THE STATES. The values integrated are
  ! themselves accurate to order p, so a rule finer than p gains
  ! nothing: the error is of order h**min(m,p) either way. Taking m
  ! as the family's own order is therefore exactly enough, and the
  ! weights grow and alternate in sign beyond it.
  !===================================================================!

  pure subroutine family_step_quadrature(this, dt, k, weight)

    class(family)         , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k
    type(derivative_terms), allocatable, intent(out) :: weight(:)

    associate (u1 => this); end associate
    if (k < 1 .or. k > size(dt)) then
       error stop 'operation_family: a quadrature is evaluated at an instant of the block'
    end if
    allocate(weight(1))
    weight(1) = derivative_terms(1.0_dp, dt(k))

  end subroutine family_step_quadrature

  !===================================================================!
  ! The quadrature weight of one stage: the tableau's b for a stage
  ! family, one for a multistep family, whose single quadrature point
  ! is the instant itself. An index outside the stages stops the
  ! program.
  !===================================================================!

  pure real(dp) function family_stage_weight(this, i)

    class(family), intent(in) :: this
    integer      , intent(in) :: i

    if (i /= 1) then
       error stop 'operation_family: a multistep family has one quadrature point per instant'
    end if
    family_stage_weight = 1.0_dp

  end function family_stage_weight

  pure integer function family_history_depth(this, equation_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: equation_degree

    associate (u1 => this, u2 => equation_degree); end associate
    family_history_depth = 1

  end function family_history_depth

  pure integer function family_primary_degree(this, equation_degree)

    class(family), intent(in) :: this
    integer      , intent(in) :: equation_degree

    associate (u1 => this); end associate
    family_primary_degree = equation_degree

  end function family_primary_degree

  !===================================================================!
  ! The offsets negated: the nodes a family reads its coefficients at
  ! lie on the opposite side of the instant they are measured from.
  !===================================================================!

  pure function negated(u) result(minus_u)

    type(derivative_terms), intent(in) :: u(0:)
    type(derivative_terms) :: minus_u(0:ubound(u, 1))

    integer :: i

    do i = 0, ubound(u, 1)
       minus_u(i) = (-1.0_dp) * u(i)
    end do

  end function negated

end module operation_family
