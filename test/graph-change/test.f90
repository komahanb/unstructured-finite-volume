!=====================================================================!
! Tests for graph_change_protocol, graph_value_map, and
! graph_value_change:
!
!  - the controller lifecycle apply -> check -> keep | revert on
!    accepted, rejected, vetoed, failing, and mixed changes, and
!    the change_result flags after each outcome
!  - the value map status transitions UNATTACHED -> UNKNOWN ->
!    KNOWN and back, and its storage rules: rows are keyed on
!    copied identity tokens rather than position, values are
!    copied on write, and rows outlive the variables that created
!    them
!  - value_change run through the same controller against every
!    possible prior status, checking that revert restores the
!    prior state exactly
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_change

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use graph_change_protocol, only : change_controller, change_result
  use graph_value_map      , only : value_map, &
       & VALUE_UNATTACHED, VALUE_UNKNOWN, VALUE_KNOWN
  use graph_value_change   , only : value_change
  use toy_changes          , only : counting_change, mixed_change

  implicit none

  type(change_controller) :: controller
  type(change_result)     :: result

  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph reversible change suite"
  write(*,'(1x,a)') "============================================="

  call check_the_lifecycle(nfail)
  call check_the_map_ontology(nfail)
  call check_the_storage_laws(nfail)
  call check_the_value_change(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all reversible change checks passed"
  else
     error stop
  end if

contains

  !===================================================================!
  ! Run the counting toy through the controller four ways -
  ! accepted, rejected, vetoed by check, and failing in apply -
  ! and verify the counter and the change_result flags after each
  ! run. Then run the mixed toy accepted and rejected, and check
  ! that reset clears every flag.
  !===================================================================!

  subroutine check_the_lifecycle(nfail)

    integer, intent(inout) :: nfail

    type(counting_change) :: change
    type(mixed_change)    :: both

    call controller % run(change, .true., result)
    call report(change % rooms == 1 .and. &
         & result % attempted .and. result % applied .and. &
         & result % checked .and. result % check_passed .and. &
         & result % accepted .and. result % kept .and. &
         & .not. result % reverted .and. .not. result % failed .and. &
         & result % touches_structure .and. .not. result % touches_value, &
         & "an accepted change is kept, and the record says so", nfail)

    call controller % run(change, .false., result)
    call report(change % rooms == 1 .and. &
         & result % reverted .and. .not. result % kept .and. &
         & .not. result % accepted, &
         & "a rejected change is reverted: the counter is restored", nfail)

    change % check_passes = .false.
    call controller % run(change, .true., result)
    call report(change % rooms == 1 .and. &
         & result % checked .and. .not. result % check_passed .and. &
         & result % reverted .and. .not. result % accepted, &
         & "a vetoed change is reverted before acceptance is asked", nfail)

    change % check_passes = .true.
    change % fail_apply   = .true.
    call controller % run(change, .true., result)
    call report(change % rooms == 1 .and. &
         & result % failed .and. result % reverted .and. &
         & .not. result % applied .and. .not. result % kept, &
         & "a failed apply is reverted and never kept", nfail)

    call controller % run(both, .true., result)
    call report(both % rooms == 1 .and. both % value == 2.0_dp .and. &
         & result % touches_structure .and. result % touches_value .and. &
         & result % kept, &
         & "an accepted mixed change keeps both mutations", nfail)

    call controller % run(both, .false., result)
    call report(both % rooms == 1 .and. both % value == 2.0_dp .and. &
         & result % reverted, &
         & "a rejected mixed change restores structure and value", nfail)

    call result % reset()
    call report(.not. (result % attempted .or. result % applied .or. &
         & result % checked .or. result % check_passed .or. &
         & result % accepted .or. result % kept .or. result % reverted &
         & .or. result % failed .or. result % touches_structure .or. &
         & result % touches_value), &
         & "reset returns the whole record to false", nfail)

  end subroutine check_the_lifecycle

  !===================================================================!
  ! Take one map entry through every status transition and verify
  ! attached / status_of / value_of after each step:
  ! UNATTACHED -> attach -> UNKNOWN -> mark_known -> KNOWN ->
  ! update -> KNOWN -> mark_unknown -> UNKNOWN -> detach ->
  ! UNATTACHED.
  !===================================================================!

  subroutine check_the_map_ontology(nfail)

    integer, intent(inout) :: nfail

    type(value_map) :: map
    type(graph)     :: a
    real(dp), allocatable :: rv(:)

    call a % declare()

    call report(.not. map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNATTACHED, &
         & "a fresh map reports UNATTACHED for an unmapped graph", nfail)

    call map % attach_unknown(a)
    call report(map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNKNOWN, &
         & "attach_unknown: attached is true and the status is UNKNOWN", nfail)

    call map % mark_known(a, [2.5_dp])
    call map % value_of(a, rv)
    call report(map % status_of(a) == VALUE_KNOWN .and. &
         & size(rv) == 1 .and. rv(1) == 2.5_dp, &
         & "mark_known: the status is KNOWN and the value reads back", nfail)

    call map % mark_known(a, [3.5_dp, 4.5_dp])
    call map % value_of(a, rv)
    call report(size(rv) == 2 .and. rv(1) == 3.5_dp .and. rv(2) == 4.5_dp, &
         & "mark_known on a KNOWN entry updates the stored values", nfail)

    call map % mark_unknown(a)
    call report(map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNKNOWN, &
         & "mark_unknown: the status returns to UNKNOWN, the entry remains", &
         & nfail)

    call map % detach(a)
    call report(.not. map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNATTACHED, &
         & "detach removes the entry: UNATTACHED again", nfail)

  end subroutine check_the_map_ontology

  !===================================================================!
  ! Storage checks: (1) detach the middle of three rows and verify
  ! the other two are still found, because lookup is by identity
  ! token and not by row position; (2) an assignment copy of a
  ! graph carries the same token and finds the same row; (3)
  ! mark_known copies the caller's array, so mutating that array
  ! afterwards must not change the stored value; (4) a row created
  ! from a subroutine-local graph must still be readable after
  ! that variable is gone.
  !===================================================================!

  subroutine check_the_storage_laws(nfail)

    integer, intent(inout) :: nfail

    type(value_map) :: map
    type(graph)     :: a, b, c, keeper
    type(graph)     :: b_copy
    real(dp), allocatable :: rv(:)
    real(dp) :: mine(2)

    call a % declare()
    call b % declare()
    call c % declare()

    call map % attach_unknown(a)
    call map % attach_unknown(b)
    call map % attach_unknown(c)
    call map % mark_known(a, [1.0_dp])
    call map % mark_known(c, [3.0_dp])

    call map % detach(b)
    call map % value_of(a, rv)
    call report(rv(1) == 1.0_dp .and. map % status_of(c) == VALUE_KNOWN &
         & .and. .not. map % attached(b), &
         & "after detaching a middle row the others are found by identity", &
         & nfail)

    b_copy = b
    call report(.not. map % attached(b_copy), &
         & "an assignment copy of a graph finds the same row", nfail)

    mine = [5.0_dp, 6.0_dp]
    call map % mark_known(a, mine)
    mine = [-9.0_dp, -9.0_dp]
    call map % value_of(a, rv)
    call report(rv(1) == 5.0_dp .and. rv(2) == 6.0_dp, &
         & "mark_known copies values: mutating the source array changes &
         &nothing", nfail)

    call attach_from_a_life(map, keeper)
    call map % value_of(keeper, rv)
    call report(map % status_of(keeper) == VALUE_KNOWN .and. &
         & rv(1) == 7.0_dp, &
         & "the map outlives the variable that built its row", nfail)

  end subroutine check_the_storage_laws

  !-------------------------------------------------------------------!
  ! Creates a map row keyed on a subroutine-local graph and
  ! returns a copy of that graph; the caller reads the row after
  ! the local variable no longer exists.
  !-------------------------------------------------------------------!

  subroutine attach_from_a_life(map, keeper)

    type(value_map), intent(inout) :: map
    type(graph)    , intent(out)   :: keeper

    type(graph) :: temporary

    call temporary % declare()
    call map % attach_unknown(temporary)
    call map % mark_known(temporary, [7.0_dp])

    keeper = temporary

  end subroutine attach_from_a_life

  !===================================================================!
  ! value_change through the controller, one case per prior
  ! status: an accepted update on an unattached graph must leave a
  ! KNOWN entry; a rejected update on a KNOWN entry must restore
  ! the old values; a vetoed update on an UNKNOWN entry must leave
  ! it UNKNOWN; a rejected update on an unattached graph must
  ! remove the entry its apply created.
  !===================================================================!

  subroutine check_the_value_change(nfail)

    integer, intent(inout) :: nfail

    type(value_map)    :: map
    type(value_change) :: update
    type(graph)        :: d, e, f
    real(dp), allocatable :: rv(:)

    call d % declare()
    call e % declare()
    call f % declare()

    call update % bind(map, d, [1.5_dp, 2.5_dp])
    call controller % run(update, .true., result)
    call map % value_of(d, rv)
    call report(result % kept .and. result % touches_value .and. &
         & .not. result % touches_structure .and. &
         & map % status_of(d) == VALUE_KNOWN .and. &
         & rv(1) == 1.5_dp .and. rv(2) == 2.5_dp, &
         & "an accepted update on an unattached graph: KNOWN and kept", nfail)

    call update % bind(map, d, [9.0_dp, 9.0_dp])
    call controller % run(update, .false., result)
    call map % value_of(d, rv)
    call report(result % reverted .and. &
         & rv(1) == 1.5_dp .and. rv(2) == 2.5_dp, &
         & "a rejected update restores the old values exactly", nfail)

    call map % attach_unknown(e)
    call update % bind(map, e, [3.0_dp], check_passes=.false.)
    call controller % run(update, .true., result)
    call report(result % reverted .and. map % attached(e) .and. &
         & map % status_of(e) == VALUE_UNKNOWN, &
         & "a vetoed update on an UNKNOWN entry leaves it UNKNOWN", nfail)

    call update % bind(map, f, [4.0_dp])
    call controller % run(update, .false., result)
    call report(result % reverted .and. .not. map % attached(f) .and. &
         & map % status_of(f) == VALUE_UNATTACHED, &
         & "a rejected update on an unattached graph removes its entry", nfail)

  end subroutine check_the_value_change

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_graph_change
