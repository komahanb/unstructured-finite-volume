!=====================================================================!
! The reversible change suite: one program over the whole change
! seam -
!
!      protocol    apply -> check -> keep | revert, owned by the
!                  one controller for structural, value, and mixed
!                  changes alike; failure is a lawful answer at
!                  every step
!      value map   graph identity -> value status x field, keyed
!                  on copied tokens and never on position, with
!                  the closed vocabulary UNATTACHED / UNKNOWN /
!                  KNOWN
!      value change  the pure attached-value member of the family:
!                  the map updated through the same controller,
!                  and revert restoring the seat EXACTLY -
!                  unattached seats leave, unknown seats stand
!                  untrusted again, known seats get their old
!                  values back
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
  ! THE LIFECYCLE, WALKED FOUR WAYS on the structural toy and once
  ! on the mixed one: acceptance keeps, rejection reverts, a veto
  ! reverts before acceptance is asked, and a failed apply reverts
  ! without ever having mutated - and every terminal record says
  ! exactly what happened.
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
         & "a rejected change is reverted: the room ungrows", nfail)

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
         & "a failed apply is a lawful answer: reverted, never kept", nfail)

    call controller % run(both, .true., result)
    call report(both % rooms == 1 .and. both % value == 2.0_dp .and. &
         & result % touches_structure .and. result % touches_value .and. &
         & result % kept, &
         & "a mixed change rides the same lifecycle: both kept", nfail)

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
  ! THE CLOSED STATUS VOCABULARY, walked whole on one seat:
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
         & "a fresh map holds nothing: unattached is a lawful answer", nfail)

    call map % attach_unknown(a)
    call report(map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNKNOWN, &
         & "attach opens the seat: attached and UNKNOWN", nfail)

    call map % mark_known(a, [2.5_dp])
    call map % value_of(a, rv)
    call report(map % status_of(a) == VALUE_KNOWN .and. &
         & size(rv) == 1 .and. rv(1) == 2.5_dp, &
         & "mark_known vouches: KNOWN, and the value reads back", nfail)

    call map % mark_known(a, [3.5_dp, 4.5_dp])
    call map % value_of(a, rv)
    call report(size(rv) == 2 .and. rv(1) == 3.5_dp .and. rv(2) == 4.5_dp, &
         & "updating a KNOWN seat is lawful: values move", nfail)

    call map % mark_unknown(a)
    call report(map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNKNOWN, &
         & "mark_unknown withdraws trust, the seat stands", nfail)

    call map % detach(a)
    call report(.not. map % attached(a) .and. &
         & map % status_of(a) == VALUE_UNATTACHED, &
         & "detach closes the seat: unattached again", nfail)

  end subroutine check_the_map_ontology

  !===================================================================!
  ! THE STORAGE LAWS. Rows key on copied identity, never position:
  ! detaching a middle row re-aims nothing, a copied graph IS its
  ! identity, values are copied at attach so the caller's array is
  ! not borrowed, and the map outlives every variable that built
  ! it.
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
         & "detaching a middle row re-aims nothing: identity, not position", &
         & nfail)

    b_copy = b
    call report(.not. map % attached(b_copy), &
         & "a copied graph IS its identity: the copy answers as the &
         &original", nfail)

    mine = [5.0_dp, 6.0_dp]
    call map % mark_known(a, mine)
    mine = [-9.0_dp, -9.0_dp]
    call map % value_of(a, rv)
    call report(rv(1) == 5.0_dp .and. rv(2) == 6.0_dp, &
         & "values are copied at attach: mutating the source moves &
         &nothing", nfail)

    call attach_from_a_life(map, keeper)
    call map % value_of(keeper, rv)
    call report(map % status_of(keeper) == VALUE_KNOWN .and. &
         & rv(1) == 7.0_dp, &
         & "the map outlives the variable that built its row", nfail)

  end subroutine check_the_storage_laws

  !-------------------------------------------------------------------!
  ! One row built from a life that ends here: the local graph dies
  ! at return, its copied token and minted field do not.
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
  ! THE VALUE CHANGE THROUGH THE ONE CONTROLLER. Every prior state
  ! of the seat is restored exactly on revert: a fresh seat leaves,
  ! an untrusted seat stands untrusted again, a known seat gets its
  ! old values back - and an accepted update simply stays.
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
         & "an accepted update vouches a fresh seat: KNOWN and kept", nfail)

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
         & "a vetoed update on an untrusted seat leaves it untrusted", nfail)

    call update % bind(map, f, [4.0_dp])
    call controller % run(update, .false., result)
    call report(result % reverted .and. .not. map % attached(f) .and. &
         & map % status_of(f) == VALUE_UNATTACHED, &
         & "a rejected update on a fresh seat detaches it whole", nfail)

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
