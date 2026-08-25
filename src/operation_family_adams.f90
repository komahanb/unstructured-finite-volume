!=====================================================================!
! The implicit Adams-Moulton family of order p.
!
! The primary unknown is the acceleration q"_k. Each lower component
! at k is the one at k - 1 plus the integral over the last step of
! the component above it, that integrand interpolated through the
! instants k .. k-p+1:
!
!      determines 1   q'_k = q'_(k-1) + dt sum_i alpha_i q"_(k-i)
!      determines 0   q_k  = q_(k-1)  + dt sum_i alpha_i q'_(k-i)
!
! with alpha_i the integral over the last step, in scaled units, of
! the i-th basis function through those instants. So one edge rule
! serves both rows: a source one degree above the constraint's, at
! offset i < p, carries alpha_i; the source of the same degree at
! offset one carries one. The row determining q reads q'_k itself
! at offset zero, and that value is determined at the same instant
! by the row above - the coupled block solves both, and no
! elimination is written here.
!
! Any other edge stops the program: a same-degree source not one
! instant back, a source more than one degree above, an offset at or
! past p, a constraint determining any other degree, or an edge from
! a later instant.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family_adams

  use util_precision  , only : dp
  use operation_family     , only : family, offsets, integral_over_step
  use util_derivative_terms, only : derivative_terms, operator(*)

  implicit none

  private
  public :: adams_family

  type, extends(family) :: adams_family

     integer, private :: order = 1

   contains

     procedure :: name             => adams_name
     procedure :: history_depth    => adams_history_depth
     procedure :: num_stages       => adams_num_stages
     procedure :: primary_degree   => adams_primary_degree
     procedure :: row_pattern      => adams_row_pattern
     procedure :: edge_coefficient => adams_edge_coefficient

  end type adams_family

  interface adams_family
     module procedure create
  end interface adams_family

contains

  function create(order) result(this)

    integer, intent(in) :: order
    type(adams_family) :: this

    if (order < 1) then
       error stop 'operation_family_adams: the order is positive'
    end if

    this % order = order
    call this % declare_arguments(3)

  end function create

  pure function adams_name(this) result(name)

    class(adams_family), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'adams-moulton'

  end function adams_name

  pure integer function adams_history_depth(this, equation_degree)

    class(adams_family), intent(in) :: this
    integer            , intent(in) :: equation_degree

    associate (u1 => equation_degree); end associate
    adams_history_depth = max(this % order - 1, 1)

  end function adams_history_depth

  pure integer function adams_num_stages(this)

    class(adams_family), intent(in) :: this

    associate (u1 => this); end associate
    adams_num_stages = 1

  end function adams_num_stages

  pure integer function adams_primary_degree(this, equation_degree)

    class(adams_family), intent(in) :: this
    integer            , intent(in) :: equation_degree

    associate (u1 => this); end associate
    adams_primary_degree = equation_degree

  end function adams_primary_degree

  !===================================================================!
  ! The row on degree d carries the same degree one instant back and
  ! quadratures the degree above it over the last p instants. The
  ! highest degree belongs to the governing constraint and has no
  ! pattern here.
  !===================================================================!

  pure subroutine adams_row_pattern(this, determines, equation_degree, &
       & offset, source_degree)

    class(adams_family), intent(in) :: this
    integer            , intent(in) :: determines, equation_degree
    integer, allocatable, intent(out) :: offset(:), source_degree(:)

    integer :: i

    if (determines < 0 .or. determines >= equation_degree) then
       allocate(offset(0), source_degree(0))
       return
    end if

    offset        = [1, (i, i = 0, this % order - 1)]
    source_degree = [determines, (determines + 1, i = 0, this % order - 1)]

  end subroutine adams_row_pattern

  !===================================================================!
  ! alpha_i at instant k: the quadrature weight over the last step of
  ! the i-th basis function through instants k .. k-p+1.
  !===================================================================!

  pure function quadrature_weight(dt, k, i, p) result(alpha)

    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: k, i, p
    type(derivative_terms) :: alpha

    alpha = integral_over_step(negated(offsets(dt, k, p)), i)

  end function quadrature_weight

  pure function negated(u) result(minus_u)

    type(derivative_terms), intent(in) :: u(0:)
    type(derivative_terms) :: minus_u(0:ubound(u, 1))

    integer :: i

    do i = 0, ubound(u, 1)
       minus_u(i) = (-1.0_dp) * u(i)
    end do

  end function negated

  pure function adams_edge_coefficient(this, dt, tail, head, &
       & source_degree, determines) result(c)

    class(adams_family)   , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, source_degree, determines
    type(derivative_terms) :: c

    integer :: i

    i = head - tail

    if (i < 0) then
       error stop 'operation_family_adams: an edge runs from an earlier instant'
    end if
    if (determines < 0) then
       error stop 'operation_family_adams: a constraint determines a degree at or above the value'
    end if

    if (source_degree == determines) then
       if (i /= 1) error stop 'operation_family_adams: the same degree is carried one instant'
       c = derivative_terms(1.0_dp, dt(head))
    else if (source_degree == determines + 1) then
       if (i >= this % order) error stop 'operation_family_adams: the quadrature reaches p instants'
       c = quadrature_weight(dt, head, i, this % order)
    else
       error stop 'operation_family_adams: a source is the constraint''s degree or one above'
    end if

  end function adams_edge_coefficient

end module operation_family_adams
