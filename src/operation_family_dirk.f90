!=====================================================================!
! The diagonally implicit Runge-Kutta family, from a Butcher tableau.
!
! The coupling here is within one step. Its vertices are
!
!      1              the instant the step leaves from
!      1 + i          stage i, i = 1 .. s
!      2 + s          the instant the step arrives at
!
! and the primary unknown at a stage is its acceleration. Each lower
! component at a stage is the incoming one plus the step times the
! tableau row applied to the component above it, over the stages at
! or before it:
!
!      u'_i = q'_in + dt sum_j a_ij u"_j        determines 1
!      u_i  = q_in  + dt sum_j a_ij u'_j        determines 0
!
! and the arriving instant reads the stages through the weights b:
!
!      q"_out = sum_j b_j u"_j                   determines 2
!      q'_out = q'_in + dt sum_j b_j u"_j        determines 1
!      q_out  = q_in  + dt sum_j b_j u'_j        determines 0
!
! So the edge rule is: from the incoming instant, one; from stage j
! into stage i, a_ij with j at or before i; from stage j into the
! arriving instant, b_j. No coefficient depends on dt, so the
! tableau serves any step.
!
! The incoming instant reaches the arriving one directly: every row
! below the highest degree carries its own degree across the step
! unchanged, which is the one in q'_out = q'_in + dt sum_j b_j u"_j.
! The highest degree has no such term, and a caller assembling the
! step simply does not make that edge.
!
! An edge from a stage after its head, an edge into the incoming
! instant, an edge out of the arriving instant, or a source degree
! that is neither the constraint's nor one above it: each stops the
! program. A tableau with an entry above its diagonal is refused at
! construction, since the stages could not then be ordered.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_family_dirk

  use iso_fortran_env , only : dp => REAL64
  use operation_family     , only : family
  use util_derivative_terms, only : derivative_terms

  implicit none

  private
  public :: dirk_family
  public :: implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage, &
       & hairer_wanner_five_stage

  type, extends(family) :: dirk_family

     real(dp), private, allocatable :: a(:,:)
     real(dp), private, allocatable :: b(:)

   contains

     procedure :: name             => dirk_name
     procedure :: history_depth    => dirk_history_depth
     procedure :: num_stages       => dirk_num_stages
     procedure :: primary_degree   => dirk_primary_degree
     procedure :: row_pattern      => dirk_row_pattern
     procedure :: edge_coefficient => dirk_edge_coefficient

  end type dirk_family

  interface dirk_family
     module procedure create
  end interface dirk_family

contains

  !===================================================================!
  ! A tableau: a square matrix a and weights b of the same extent.
  ! An entry above the diagonal, or extents that disagree, stops the
  ! program.
  !===================================================================!

  function create(a, b) result(this)

    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in) :: b(:)
    type(dirk_family) :: this

    integer :: i, j

    if (size(a, 1) /= size(a, 2) .or. size(b) /= size(a, 1)) then
       error stop 'operation_family_dirk: the tableau is square with one weight per stage'
    end if

    do i = 1, size(a, 1)
       do j = i + 1, size(a, 2)
          if (a(i, j) /= 0.0_dp) then
             error stop 'operation_family_dirk: a diagonally implicit tableau has no entry above its diagonal'
          end if
       end do
    end do

    this % a = a
    this % b = b
    call this % declare_arguments(3)

  end function create

  pure function dirk_name(this) result(name)

    class(dirk_family), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'dirk'

  end function dirk_name

  pure integer function dirk_history_depth(this)

    class(dirk_family), intent(in) :: this

    associate (u1 => this); end associate
    dirk_history_depth = 1

  end function dirk_history_depth

  pure integer function dirk_num_stages(this)

    class(dirk_family), intent(in) :: this

    dirk_num_stages = size(this % b)

  end function dirk_num_stages

  pure integer function dirk_primary_degree(this, equation_degree)

    class(dirk_family), intent(in) :: this
    integer           , intent(in) :: equation_degree

    associate (u1 => this); end associate
    dirk_primary_degree = equation_degree

  end function dirk_primary_degree

  !===================================================================!
  ! A stage family's rows run between the stages of one step, not
  ! between instants, so they have no pattern in instant offsets. A
  ! caller assembling a stage block reads the tableau through
  ! num_stages and edge_coefficient instead, and an empty pattern is
  ! how this family says so.
  !===================================================================!

  pure subroutine dirk_row_pattern(this, determines, equation_degree, &
       & offset, source_degree)

    class(dirk_family), intent(in) :: this
    integer           , intent(in) :: determines, equation_degree
    integer, allocatable, intent(out) :: offset(:), source_degree(:)

    associate (u1 => this, u2 => determines, u3 => equation_degree); end associate

    allocate(offset(0), source_degree(0))

  end subroutine dirk_row_pattern

  pure function dirk_edge_coefficient(this, dt, tail, head, &
       & source_degree, determines) result(c)

    class(dirk_family)    , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, source_degree, determines
    type(derivative_terms) :: c

    integer :: s, i, j

    s = size(this % b)

    if (head == 1 .or. tail == 2 + s) then
       error stop 'operation_family_dirk: an edge runs from the incoming instant or a stage into a later vertex'
    end if
    if (source_degree /= determines .and. source_degree /= determines + 1) then
       error stop 'operation_family_dirk: a source is the constraint''s degree or one above'
    end if

    if (tail == 1) then
       c = derivative_terms(1.0_dp, dt(head))
       return
    end if

    j = tail - 1
    if (head == 2 + s) then
       c = derivative_terms(this % b(j), dt(head))
       return
    end if

    i = head - 1
    if (j > i) then
       error stop 'operation_family_dirk: a stage reads stages at or before it'
    end if
    c = derivative_terms(this % a(i, j), dt(head))

  end function dirk_edge_coefficient

  !===================================================================!
  ! Four tableaux. The implicit midpoint rule, order two; the
  ! Crouzeix two-stage, order three; the Crouzeix three-stage, order
  ! four; and the Hairer-Wanner five-stage L-stable method, order
  ! four.
  !===================================================================!

  function implicit_midpoint() result(this)

    type(dirk_family) :: this

    this = dirk_family(reshape([0.5_dp], [1, 1]), [1.0_dp])

  end function implicit_midpoint

  function crouzeix_two_stage() result(this)

    type(dirk_family) :: this

    real(dp) :: g

    g = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp

    this = dirk_family(reshape([g, 1.0_dp - 2.0_dp * g, 0.0_dp, g], [2, 2]), &
         & [0.5_dp, 0.5_dp])

  end function crouzeix_two_stage

  function crouzeix_three_stage() result(this)

    type(dirk_family) :: this

    real(dp) :: g, w, pi

    pi = acos(-1.0_dp)
    g  = cos(pi / 18.0_dp) / sqrt(3.0_dp) + 0.5_dp
    w  = 1.0_dp / (6.0_dp * (1.0_dp - 2.0_dp * g)**2)

    this = dirk_family(reshape( &
         & [g,             0.5_dp - g, 2.0_dp * g,      &
         &  0.0_dp,        g,          1.0_dp - 4.0_dp * g, &
         &  0.0_dp,        0.0_dp,     g], [3, 3]),    &
         & [w, 1.0_dp - 2.0_dp * w, w])

  end function crouzeix_three_stage

  function hairer_wanner_five_stage() result(this)

    type(dirk_family) :: this

    real(dp) :: a(5, 5)

    a = 0.0_dp
    a(1, 1:1) = [1.0_dp / 4.0_dp]
    a(2, 1:2) = [1.0_dp / 2.0_dp, 1.0_dp / 4.0_dp]
    a(3, 1:3) = [17.0_dp / 50.0_dp, -1.0_dp / 25.0_dp, 1.0_dp / 4.0_dp]
    a(4, 1:4) = [371.0_dp / 1360.0_dp, -137.0_dp / 2720.0_dp, 15.0_dp / 544.0_dp, 1.0_dp / 4.0_dp]
    a(5, 1:5) = [25.0_dp / 24.0_dp, -49.0_dp / 48.0_dp, 125.0_dp / 16.0_dp, -85.0_dp / 12.0_dp, &
         &       1.0_dp / 4.0_dp]

    this = dirk_family(a, a(5, :))

  end function hairer_wanner_five_stage

end module operation_family_dirk
