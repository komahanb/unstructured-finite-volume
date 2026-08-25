!=====================================================================!
! The backward-difference family of order p.
!
! The primary unknown is the value q_k. Every source is a value at an
! earlier instant, and two constraints read them:
!
!      determines 1   the velocity constraint, coefficient on
!                     q_(k-j) the slope at zero of the j-th basis
!                     function through instants k .. k-p
!      determines d   the constraint on the d-th derivative: the
!                     velocity constraint applied d times, each inner
!                     one at instant k - i with its own offsets and
!                     the step ratio carrying its 1/dt onto instant
!                     k's,
!
!         c(d, j at k) = sum over i of
!                        alpha_i(at k) c(d-1, j-i at k-i) dt_k/dt_(k-i)
!
!                     which reaches d p instants back, so the widest
!                     row of an equation of degree N reaches N p, and
!                     the history a block needs grows with the
!                     equation as well as with the order. At d = 2 on a
!                     uniform grid this is the convolution of alpha
!                     with itself, and the tabulated second-difference
!                     coefficients come out of it.
!
! The primary unknown is the value, so the governing constraint takes
! the degree-zero row and every degree above it is a derived row.
!
! An edge from a source that is not a value, an edge into a row on
! the d-th derivative that reaches past d p instants, a row of degree
! below one, or an edge running from a later instant: each stops the
! program, because the family defines no coefficient for it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family_bdf

  use util_precision  , only : dp
  use operation_family     , only : family, offsets, slope_at_zero, negated
  use util_derivative_terms, only : derivative_terms, &
       & operator(+), operator(*), operator(/)

  implicit none

  private
  public :: bdf_family

  type, extends(family) :: bdf_family

     integer, private :: order = 1

   contains

     procedure :: name             => bdf_name
     procedure :: history_depth    => bdf_history_depth
     procedure :: num_stages       => bdf_num_stages
     procedure :: primary_degree   => bdf_primary_degree
     procedure :: row_pattern      => bdf_row_pattern
     procedure :: edge_coefficient => bdf_edge_coefficient

  end type bdf_family

  interface bdf_family
     module procedure create
  end interface bdf_family

contains

  !===================================================================!
  ! An order below one stops the program.
  !===================================================================!

  function create(order) result(this)

    integer, intent(in) :: order
    type(bdf_family) :: this

    if (order < 1) then
       error stop 'operation_family_bdf: the order is positive'
    end if

    this % order = order
    call this % declare_arguments(3)

  end function create

  pure function bdf_name(this) result(name)

    class(bdf_family), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'bdf'

  end function bdf_name

  pure integer function bdf_history_depth(this, equation_degree)

    class(bdf_family), intent(in) :: this
    integer          , intent(in) :: equation_degree

    bdf_history_depth = equation_degree * this % order

  end function bdf_history_depth

  pure integer function bdf_num_stages(this)

    class(bdf_family), intent(in) :: this

    associate (u1 => this); end associate
    bdf_num_stages = 1

  end function bdf_num_stages

  pure integer function bdf_primary_degree(this, equation_degree)

    class(bdf_family), intent(in) :: this
    integer          , intent(in) :: equation_degree

    associate (u1 => this, u2 => equation_degree); end associate
    bdf_primary_degree = 0

  end function bdf_primary_degree

  !===================================================================!
  ! The row on the d-th derivative reads the values over the d p
  ! instants it reaches back. The degree-zero row belongs to the
  ! governing constraint and has no pattern here.
  !===================================================================!

  pure subroutine bdf_row_pattern(this, determines, equation_degree, &
       & offset, source_degree)

    class(bdf_family), intent(in) :: this
    integer          , intent(in) :: determines, equation_degree
    integer, allocatable, intent(out) :: offset(:), source_degree(:)

    integer :: j, reach

    if (determines < 1 .or. determines > equation_degree) then
       allocate(offset(0), source_degree(0))
       return
    end if

    reach  = determines * this % order
    offset = [(j, j = 0, reach)]
    allocate(source_degree(reach + 1), source=0)

  end subroutine bdf_row_pattern

  !===================================================================!
  ! alpha_j at instant k. The nodes are the offsets negated, since
  ! instant k - j lies at time t_k - theta_j dt_k.
  !===================================================================!

  pure function velocity_coefficient(dt, k, j, p) result(alpha)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, j, p
    type(derivative_terms) :: alpha

    alpha = slope_at_zero(negated(offsets(dt, k, p + 1)), j)

  end function velocity_coefficient

  !===================================================================!
  ! The coefficient of the row on the d-th derivative, at offset j
  ! from instant k: the velocity coefficients composed d times, each
  ! inner one read at the instant its outer factor points to.
  !===================================================================!

  pure recursive function derivative_coefficient(dt, k, j, p, d) result(c)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, j, p, d
    type(derivative_terms) :: c

    integer :: i

    if (d == 1) then
       c = velocity_coefficient(dt, k, j, p)
       return
    end if

    c = derivative_terms(0.0_dp, dt(k))

    do i = max(0, j - (d - 1) * p), min(j, p)
       c = c + velocity_coefficient(dt, k, i, p) &
            & * derivative_coefficient(dt, k - i, j - i, p, d - 1) &
            & * dt(k) / dt(k - i)
    end do

  end function derivative_coefficient

  pure function bdf_edge_coefficient(this, dt, tail, head, &
       & source_degree, determines) result(c)

    class(bdf_family)     , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, source_degree, determines
    type(derivative_terms) :: c

    integer :: j

    j = head - tail

    if (j < 0) then
       error stop 'operation_family_bdf: an edge runs from an earlier instant'
    end if
    if (source_degree /= 0) then
       error stop 'operation_family_bdf: every source is a value'
    end if

    if (determines < 1) then
       error stop 'operation_family_bdf: a derived row determines a derivative'
    end if
    if (j > determines * this % order) then
       error stop 'operation_family_bdf: the row on the d-th derivative reaches d p instants'
    end if

    c = derivative_coefficient(dt, head, j, this % order, determines)

  end function bdf_edge_coefficient

end module operation_family_bdf
