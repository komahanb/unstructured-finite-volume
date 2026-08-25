!=====================================================================!
! The derived constraints of a scheme, assembled as one stencil.
!
!             THE SIGN CONVENTION, OWNED HERE
!
! A constraint that determines a component reads
!
!      residual  =  (the component it determines)
!                     -  sum over its edges of weight * (its source)
!
! so the determined component's column carries plus one and every
! source column carries minus its weight. The constant is zero: a
! derived constraint is homogeneous, and the affine part belongs to
! the initial conditions, which are not assembled here.
!
! That convention is written once, in this module, so that no caller
! chooses a sign. A caller supplies, for each edge, the row it enters
! and the column it reads, together with the weight; the diagonal is
! added here.
!
!             THE SQUARE BLOCK
!
! Rows and columns share one index space, so the result is square in
! the unknowns. One row per component: at each instant the family's
! primary degree is the row the governing constraint occupies and
! every other degree is a row assembled here. A component that no
! constraint determines leaves an empty row, which is what the
! initial conditions fill.
!
!             THE JACOBIAN IS THE SAME OBJECT
!
! A derived constraint is linear in the state, so the stencil built
! here IS its Jacobian: the stencil's partial action applies the same
! weights to a direction. There is no second object to build and no
! second path to keep in step.
!
!             WHAT IS REFUSED
!
! An edge whose source is the component its own constraint
! determines: the plus one and the minus weight would then fall in
! the same entry, which hides a family that has declared a circular
! row. Two edges entering the same row from the same column: the
! caller has supplied one dependency twice, and combining them
! silently would hide it. An index outside the unknowns. Arrays of
! disagreeing extent.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_scheme_stencil

  use util_precision  , only : dp
  use operation_stencil, only : stencil

  implicit none

  private
  public :: derived_constraints

contains

  !===================================================================!
  ! Refuse indices outside the unknowns, extents that disagree, and
  ! an edge that reads what its own row determines.
  !===================================================================!

  subroutine require_within(determined, source, weight, num_unknowns)

    integer , intent(in) :: determined(:), source(:)
    real(dp), intent(in) :: weight(:)
    integer , intent(in) :: num_unknowns

    if (size(source) /= size(determined) .or. size(weight) /= size(determined)) then
       error stop 'operation_scheme_stencil: one row, one column and one weight per edge'
    end if

    if (any(determined < 1) .or. any(determined > num_unknowns) .or. &
         & any(source < 1) .or. any(source > num_unknowns)) then
       error stop 'operation_scheme_stencil: every index names an unknown'
    end if

    if (any(determined == source)) then
       error stop 'operation_scheme_stencil: a constraint reads no source it determines itself'
    end if

  end subroutine require_within

  !===================================================================!
  ! Refuse two edges that enter the same row from the same column.
  ! The edges are counted into their rows first, so each row's
  ! columns are compared only among themselves: a row holds as many
  ! columns as the scheme reaches, which is small, and the whole
  ! check costs one pass over the edges and one over the unknowns.
  !===================================================================!

  subroutine require_distinct(determined, source, num_unknowns)

    integer, intent(in) :: determined(:), source(:), num_unknowns

    integer, allocatable :: start(:), order(:)
    integer :: r, i, j

    call grouped_by_row(determined, num_unknowns, start, order)

    do r = 1, num_unknowns
       do i = start(r + 1), start(r + 2) - 1
          do j = start(r + 1), i - 1
             if (source(order(i)) == source(order(j))) then
                error stop 'operation_scheme_stencil: a row reads each column once'
             end if
          end do
       end do
    end do

  end subroutine require_distinct

  !===================================================================!
  ! The edges counted into their rows: start(r+1) is where row r
  ! begins in order, and order lists the edges row by row.
  !===================================================================!

  subroutine grouped_by_row(determined, num_unknowns, start, order)

    integer, intent(in) :: determined(:), num_unknowns
    integer, allocatable, intent(out) :: start(:), order(:)

    integer, allocatable :: place(:)
    integer :: e, r

    allocate(start(num_unknowns + 2), source=0)

    do e = 1, size(determined)
       start(determined(e) + 2) = start(determined(e) + 2) + 1
    end do

    start(1:2) = 1
    do r = 2, num_unknowns + 1
       start(r + 1) = start(r + 1) + start(r)
    end do

    allocate(order(size(determined)))
    place = start(1:num_unknowns + 1)

    do e = 1, size(determined)
       r = determined(e)
       order(place(r + 1)) = e
       place(r + 1) = place(r + 1) + 1
    end do

  end subroutine grouped_by_row

  !===================================================================!
  ! The rows a scheme's weights make. The distinct rows entered by
  ! the edges are the constraints, and each takes a plus one on its
  ! own column.
  !===================================================================!

  function derived_constraints(determined, source, weight, num_unknowns, label) &
       & result(rows)

    integer         , intent(in)           :: determined(:)
    integer         , intent(in)           :: source(:)
    real(dp)        , intent(in)           :: weight(:)
    integer         , intent(in)           :: num_unknowns
    character(len=*), intent(in), optional :: label
    type(stencil) :: rows

    integer , allocatable :: diagonal(:)
    real(dp), allocatable :: zero(:)

    call require_within(determined, source, weight, num_unknowns)
    call require_distinct(determined, source, num_unknowns)

    diagonal = rows_entered(determined, num_unknowns)
    allocate(zero(num_unknowns), source=0.0_dp)

    rows = stencil([determined, diagonal], [source, diagonal], &
         & [-weight, spread(1.0_dp, 1, size(diagonal))], zero, label)

  end function derived_constraints

  !===================================================================!
  ! The distinct rows the edges enter, in increasing order. Each one
  ! is a constraint and takes the plus one on its own column.
  !===================================================================!

  function rows_entered(determined, num_unknowns) result(diagonal)

    integer, intent(in) :: determined(:), num_unknowns
    integer, allocatable :: diagonal(:)

    logical, allocatable :: entered(:)
    integer :: e, r, n

    allocate(entered(num_unknowns), source=.false.)

    do e = 1, size(determined)
       entered(determined(e)) = .true.
    end do

    allocate(diagonal(count(entered)))
    n = 0

    do r = 1, num_unknowns
       if (entered(r)) then
          n = n + 1
          diagonal(n) = r
       end if
    end do

  end function rows_entered

end module operation_scheme_stencil
