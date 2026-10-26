!=====================================================================!
! Module that assembles the derived constraints of a scheme as one
! stencil.
!
! The sign convention, owned here. A constraint that determines a
! component reads
!
!      residual  =  (the component it determines)
!                     -  sum over its edges of weight * (its source)
!
! so the determined component's column has plus one and every source
! column has minus its weight:
!
!                  q(k-2)    q(k-1)    q(k)      q'(k)      the columns
!   row q'(k)      -w_2      -w_1      -w_0      +1
!
! The constant is zero: a derived constraint is homogeneous, and the
! affine part belongs to the initial conditions, which are not
! assembled here.
!
! The convention is written once, in this module, so that no caller
! chooses a sign. A caller supplies, for each edge, the row it enters
! and the column it reads, together with the weight; the diagonal is
! added here.
!
! The square block. Rows and columns share one index space, so the
! result is square in the unknowns. There is one row per component:
! at each instant the family's primary degree is the row the governing
! constraint occupies, and every other degree is a row assembled here.
! A component that no constraint determines leaves an empty row, which
! the initial conditions fill.
!
! The Jacobian is the same object. A derived constraint is linear in
! the state, so the stencil built here is its Jacobian: the stencil's
! partial action applies the same weights to a direction. There is no
! second object to build and no second code path to keep consistent.
!
! What is refused:
!
!   - an edge whose source is the component its own constraint
!     determines: the plus one and the minus weight would then fall in
!     the same entry, which hides a family that has declared a
!     circular row;
!   - two edges entering the same row from the same column: the caller
!     has supplied one dependency twice, and combining them silently
!     would hide it;
!   - an index outside the unknowns;
!   - arrays of disagreeing extent.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_scheme_stencil

  use util_precision   , only : dp
  use operation_stencil, only : stencil
  use relation_binary  , only : group_by_key

  implicit none

  private

  public :: derived_constraints

contains

  !===================================================================!
  ! Refuses indices outside the unknowns, extents that disagree, and an
  ! edge that reads what its own row determines
  !===================================================================!

  subroutine require_within(determined, source, weight, num_unknowns)

    integer , intent(in) :: determined(:), source(:)
    real(dp), intent(in) :: weight(:)
    integer , intent(in) :: num_unknowns

    character(len=200) :: message
    integer            :: bad

    if (size(source) /= size(determined) .or. size(weight) /= size(determined)) then
       write(message,'(a,i0,a,i0,a,i0)') 'operation_scheme_stencil: require_within requires &
            &one row, one column and one weight per edge; size(determined) = ', size(determined), &
            & ', size(source) = ', size(source), ', size(weight) = ', size(weight)
       error stop trim(message)
    end if

    if (any(determined < 1) .or. any(determined > num_unknowns) .or. &
         & any(source < 1) .or. any(source > num_unknowns)) then
       write(message,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'operation_scheme_stencil: require_within &
            &requires every index in 1..num_unknowns (', num_unknowns, '); determined ranges &
            &from ', minval(determined), ' to ', maxval(determined), ', source ranges from ', &
            & minval(source), ' to ', maxval(source)
       error stop trim(message)
    end if

    if (any(determined == source)) then
       bad = findloc(determined == source, .true., dim=1)
       write(message,'(a,i0,a,i0)') 'operation_scheme_stencil: require_within requires a &
            &constraint to read no source it determines itself; edge ', bad, ' has &
            &determined = source = ', determined(bad)
       error stop trim(message)
    end if

  end subroutine require_within

  !===================================================================!
  ! Refuses two edges that enter the same row from the same column.
  !
  ! The edges are counted into their rows first, so each row's columns
  ! are compared only among themselves. A row has as many columns as
  ! the scheme reaches, which is few, and the whole check costs one
  ! pass over the edges and one over the unknowns:
  !
  !   edges grouped by row      row 1: e3 e7         row 2: e1 e4 e9
  !   compared within a row     e3 / e7              e1 / e4, e1 / e9, e4 / e9
  !===================================================================!

  subroutine require_distinct(determined, source, num_unknowns)

    integer, intent(in) :: determined(:), source(:), num_unknowns

    integer, allocatable :: start(:), order(:)
    integer              :: r, i, j, e
    character(len=200)   :: message

    call group_by_key(num_unknowns, determined, [(e, e = 1, size(determined))], start, order)

    do r = 1, num_unknowns
       do i = start(r), start(r + 1) - 1
          do j = start(r), i - 1
             if (source(order(i)) == source(order(j))) then
                write(message,'(a,i0,a,i0)') 'operation_scheme_stencil: require_distinct &
                     &requires a row to read each column once; row ', r, ' reads column ', &
                     & source(order(i))
                error stop trim(message)
             end if
          end do
       end do
    end do

  end subroutine require_distinct

  !===================================================================!
  ! Returns the rows a scheme's weights make. The distinct rows the
  ! edges enter are the constraints, and each takes a plus one on its
  ! own column.
  !===================================================================!

  function derived_constraints(determined, source, weight, num_unknowns, label) &
       & result(rows)

    integer         , intent(in)           :: determined(:)
    integer         , intent(in)           :: source(:)
    real(dp)        , intent(in)           :: weight(:)
    integer         , intent(in)           :: num_unknowns
    character(len=*), intent(in), optional :: label
    type(stencil)                          :: rows

    integer , allocatable :: diagonal(:)
    real(dp), allocatable :: zero(:)

    call require_within(determined, source, weight, num_unknowns)
    call require_distinct(determined, source, num_unknowns)

    diagonal = rows_entered(determined, num_unknowns)
    allocate(zero(num_unknowns), source=0.0_dp)

    ! The sources' minus weights, then the plus ones on the diagonal
    rows = stencil([determined, diagonal], [source, diagonal], &
         & [-weight, spread(1.0_dp, 1, size(diagonal))], zero, label)

  end function derived_constraints

  !===================================================================!
  ! Returns the distinct rows the edges enter, in increasing order.
  ! Each one is a constraint and takes the plus one on its own column.
  !===================================================================!

  function rows_entered(determined, num_unknowns) result(diagonal)

    integer, intent(in)  :: determined(:), num_unknowns
    integer, allocatable :: diagonal(:)

    logical, allocatable :: entered(:)
    integer              :: e, r, n

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
