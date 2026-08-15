!=====================================================================!
! The computational graph suite: the laws of the epistemic pair
! (COMPUTATIONAL-GRAPH.md).
!
! G = (Q, R): two seats, four states, one type. The checks below
! hold the classifier to exactly-one, keep bottom distinct from
! empty - an allocated zero-length payload is REALIZED data - keep
! void distinct from structurally empty - a void graph rides a
! whole GAMMA = (S, P) untroubled - and keep realized distinct
! from solved: an inconsistent pair is a realized graph, and this
! module cannot tell, by design.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_state

  use graph_carrier  , only : counted_set
  use graph_relation , only : stored_relation, slot
  use graph_structure, only : relational_graph, held_set, &
       &                      held_relation
  use graph_state    , only : computational_graph, &
       &                      GRAPH_STATE_VOID, GRAPH_STATE_DATA, &
       &                      GRAPH_STATE_OPERATOR, &
       &                      GRAPH_STATE_REALIZED, state_name, &
       &                      void_graph, data_graph, operator_graph, &
       &                      realized_graph

  implicit none

  !===================================================================!
  ! The payloads. Q is not one field and R is not graph_operation:
  ! any realized knowledge may take a seat, so the fixtures are
  ! deliberately nobody's production types.
  !===================================================================!

  type :: observed_values

     real, allocatable :: entries(:)

  end type observed_values

  type :: residual_stub

     character(len=:), allocatable :: law

  end type residual_stub

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "computational graph suite (COMPUTATIONAL-GRAPH.md)"
  write(*,'(1x,a)') "============================================="

  call check_four_states(nfail)
  call check_exactly_one(nfail)
  call check_bottom_is_not_empty(nfail)
  call check_void_is_not_topological(nfail)
  call check_realized_is_not_solved(nfail)
  call check_canonical_names(nfail)
  call check_seats_answer_what_was_handed(nfail)
  call check_graph_identity(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all computational-state checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " state check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Two seats, four states: each constructor lands its graph in the
  ! state its name promises, and the two queries agree with it.
  !===================================================================!

  subroutine check_four_states(nfail)

    integer, intent(inout) :: nfail

    type(observed_values)     :: q
    type(residual_stub)       :: r
    type(computational_graph) :: g00, g10, g01, g11

    allocate(q % entries, source=[2.0, 0.0])
    r % law = 'qdot + S(q) = 0'

    g00 = void_graph('nothing yet')
    g10 = data_graph('measurements', q)
    g01 = operator_graph('governing law', r)
    g11 = realized_graph('candidate pair', q, r)

    call report(.not. g00 % has_data() .and. &
         &      .not. g00 % has_operator() .and. &
         &      g00 % state() .eq. GRAPH_STATE_VOID, &
         & "the void graph: neither seat, and it says so", nfail)

    call report(g10 % has_data() .and. &
         &      .not. g10 % has_operator() .and. &
         &      g10 % state() .eq. GRAPH_STATE_DATA, &
         & "the data graph: Q realized, R bottom", nfail)

    call report(.not. g01 % has_data() .and. &
         &      g01 % has_operator() .and. &
         &      g01 % state() .eq. GRAPH_STATE_OPERATOR, &
         & "the operator graph: R realized, Q bottom", nfail)

    call report(g11 % has_data() .and. &
         &      g11 % has_operator() .and. &
         &      g11 % state() .eq. GRAPH_STATE_REALIZED, &
         & "the realized graph: both seats occupied", nfail)

  end subroutine check_four_states

  !===================================================================!
  ! The classifier is a partition: for every construction, exactly
  ! one of the four constants answers, and each state constant is
  ! equivalent to its row of the truth table.
  !===================================================================!

  subroutine check_exactly_one(nfail)

    integer, intent(inout) :: nfail

    type(observed_values)     :: q
    type(residual_stub)       :: r
    type(computational_graph) :: gs(4)
    integer                   :: states(4)
    integer                   :: k
    logical                   :: one, rows

    allocate(q % entries, source=[1.0])
    r % law = 'R(q) = 0'

    gs(1) = void_graph('g00')
    gs(2) = data_graph('g10', q)
    gs(3) = operator_graph('g01', r)
    gs(4) = realized_graph('g11', q, r)

    states = [GRAPH_STATE_VOID, GRAPH_STATE_DATA, &
         &    GRAPH_STATE_OPERATOR, GRAPH_STATE_REALIZED]

    one  = .true.
    rows = .true.
    do k = 1, 4
       one = one .and. (count(states .eq. gs(k) % state()) .eq. 1)
       rows = rows .and. &
            & ((gs(k) % state() .eq. GRAPH_STATE_VOID) .eqv. &
            &  (.not. gs(k) % has_data() .and. &
            &   .not. gs(k) % has_operator())) .and. &
            & ((gs(k) % state() .eq. GRAPH_STATE_DATA) .eqv. &
            &  (gs(k) % has_data() .and. &
            &   .not. gs(k) % has_operator())) .and. &
            & ((gs(k) % state() .eq. GRAPH_STATE_OPERATOR) .eqv. &
            &  (.not. gs(k) % has_data() .and. &
            &   gs(k) % has_operator())) .and. &
            & ((gs(k) % state() .eq. GRAPH_STATE_REALIZED) .eqv. &
            &  (gs(k) % has_data() .and. gs(k) % has_operator()))
    end do

    call report(one, &
         & "exactly one state holds, for every construction", nfail)
    call report(rows, &
         & "and each state is its row of the truth table, no other", &
         & nfail)

  end subroutine check_exactly_one

  !===================================================================!
  ! Bottom is not empty. An allocated payload with zero entries is
  ! realized knowledge whose content happens to be small; an
  ! unallocated seat is no knowledge at all.
  !===================================================================!

  subroutine check_bottom_is_not_empty(nfail)

    integer, intent(inout) :: nfail

    type(observed_values)     :: none_measured
    type(computational_graph) :: honest, ignorant

    allocate(none_measured % entries(0))

    honest   = data_graph('an empty run, honestly recorded', &
         &                none_measured)
    ignorant = void_graph('no run at all')

    call report(honest % has_data() .and. &
         &      honest % state() .eq. GRAPH_STATE_DATA, &
         & "an allocated zero-length payload is realized data", nfail)

    call report(.not. ignorant % has_data() .and. &
         &      ignorant % state() .eq. GRAPH_STATE_VOID, &
         & "an unallocated seat is bottom, not a small value", nfail)

  end subroutine check_bottom_is_not_empty

  !===================================================================!
  ! Void speaks of knowledge, never topology: a void computational
  ! graph rides a whole, materialized GAMMA = (S, P) - the
  ! calculator of AGENTS.md 8.2 - and stays void.
  !===================================================================!

  subroutine check_void_is_not_topological(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)                :: ops, vals, ports
    type(stored_relation)            :: flow
    type(relational_graph), target   :: gamma
    type(computational_graph)        :: g
    class(relational_graph), pointer :: host

    ops   = counted_set('operations', 3)
    vals  = counted_set('values'    , 5)
    ports = counted_set('ports'     , 2)

    flow = stored_relation('flow', &
         & [slot(ops), slot(vals), slot(ports)], &
         & reshape([1,2,1,  1,3,1,  1,5,2,  2,5,1,  2,4,2], [3, 5]))

    gamma = relational_graph('calculator', &
         & [held_set(ops), held_set(vals), held_set(ports)], &
         & [held_relation(flow)])

    g = void_graph('all unknowns', structure=gamma)

    host => g % structure()
    call report(g % state() .eq. GRAPH_STATE_VOID .and. &
         &      host % same_as(gamma), &
         & "a void graph rides a whole structure, and stays void", &
         & nfail)

  end subroutine check_void_is_not_topological

  !===================================================================!
  ! Realized is not solved: the constructor asserts occupancy and
  ! nothing else, so a deliberately inconsistent pair seats itself
  ! without complaint. Satisfied, consistent, converged are other
  ! words, asserted elsewhere.
  !===================================================================!

  subroutine check_realized_is_not_solved(nfail)

    integer, intent(inout) :: nfail

    type(observed_values)     :: q
    type(residual_stub)       :: r
    type(computational_graph) :: g

    allocate(q % entries, source=[1.0, 2.0])
    r % law = 'q = 0'

    g = realized_graph('deliberately inconsistent', q, r)

    call report(g % state() .eq. GRAPH_STATE_REALIZED, &
         & "an inconsistent pair is realized - solved is another word", &
         & nfail)

  end subroutine check_realized_is_not_solved

  !===================================================================!
  ! The canonical names, and no synonyms: void, data, operator,
  ! realized.
  !===================================================================!

  subroutine check_canonical_names(nfail)

    integer, intent(inout) :: nfail

    call report(state_name(GRAPH_STATE_VOID)     == 'void graph' &
         & .and. state_name(GRAPH_STATE_DATA)     == 'data graph' &
         & .and. state_name(GRAPH_STATE_OPERATOR) == 'operator graph' &
         & .and. state_name(GRAPH_STATE_REALIZED) == 'realized graph', &
         & "the four states answer their canonical names", nfail)

  end subroutine check_canonical_names

  !===================================================================!
  ! A seat answers what it was handed: the accessors reference the
  ! owned occupants, recovered whole across the polymorphic seat.
  !===================================================================!

  subroutine check_seats_answer_what_was_handed(nfail)

    integer, intent(inout) :: nfail

    type(observed_values)             :: q
    type(residual_stub)               :: r
    type(computational_graph), target :: g
    class(*), pointer                 :: seat
    logical                           :: ok

    allocate(q % entries, source=[2.0, 0.0])
    r % law = 'qdot + S(q) = 0'

    g = realized_graph('candidate pair', q, r)

    ok = .false.
    seat => g % data()
    select type (seat)
    type is (observed_values)
       ok = size(seat % entries) .eq. 2 .and. &
            & all(seat % entries .eq. [2.0, 0.0])
    end select
    call report(ok, "the data seat answers the payload it was handed", &
         & nfail)

    ok = .false.
    seat => g % residual()
    select type (seat)
    type is (residual_stub)
       ok = seat % law == 'qdot + S(q) = 0'
    end select
    call report(ok, "and the residual seat answers its law", nfail)

  end subroutine check_seats_answer_what_was_handed

  !===================================================================!
  ! The fourth citizen on the identity roll.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(computational_graph) :: g, h, copy

    g = void_graph('one')
    h = void_graph('two')

    call report(.not. g % same_as(h), &
         & "two declarations are two graphs", nfail)
    copy = g
    call report(copy % same_as(g), &
         & "a copy is the same declared graph", nfail)
    call report(g % name() == 'one', &
         & "and the name is the reader's, as everywhere", nfail)

  end subroutine check_graph_identity

end program test_graph_state
