!=====================================================================!
! CALCULATOR TOWER . LEVEL 3 . THE RELATIONAL GRAPH
!
! The level answers one question: HOW CARRIERS AND RELATIONS COEXIST
! AS ONE STRUCTURE. The persistent calculator becomes
!
!      GAMMA = ( { X, O, P }, { T_flow, D } )
!
! where D is ADMITTED, not redefined: it is derived here once more
! by the approved Level-2 road - restrict, project, compose - and
! the graph receives the materialized result. The pair (+, x) is
! never written on the way in. No provenance machinery rides the
! graph; the derivation stands in this file, as the level contract
! asks.
!
! The strongest negative truth of the level is the import list: the
! calculator compiles as a relational graph WITHOUT the ordinary
! profile - no vertex, no edge, no tail, no head anywhere.
!
! A second, smaller graph witnesses the coexistence law: D and its
! own restriction to the times-only subset - the empty relation
! over O x O - share one signature and remain two citizens.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_3

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : relational_graph, held_set, held_relation

  implicit none

  type(counted_set)              :: x, o, p
  type(subset_set)               :: p_out, p_in
  type(stored_relation)          :: flow, t_out3, t_in3
  type(stored_relation)          :: produces, consumes
  class(relation), allocatable   :: d
  type(relational_graph), target :: g
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 3 . relational graph"
  write(*,'(1x,a)') "============================================="

  x = counted_set('value-slots', 5)
  o = counted_set('operations' , 2)
  p = counted_set('ports'      , 3)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]

  flow = stored_relation('flow', [o, x, p], table)

  ! The approved Level-2 road, walked once more; the graph admits
  ! what the algebra derived.
  p_out    = subset_set('output-port', p, [PORT_OUT])
  p_in     = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  t_out3   = restrict_slot(flow, 3, p_out)
  t_in3    = restrict_slot(flow, 3, p_in)
  produces = project_slots(t_out3, [1, 2])
  consumes = project_slots(t_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

  g = relational_graph('calculator', &
       & [held_set(x), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])

  call check_ownership(nfail)
  call check_signature_closure(nfail)
  call check_ternary_stands(nfail)
  call check_coexistence(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! The graph owns exactly the three carriers and the two relations
  ! it was handed - by structural identity, through the generators
  ! alone: holds_set for the carriers, and a scan of relation_at
  ! identities for the relations, composed here rather than added
  ! to the contract.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(g % num_member_sets() .eq. 3, &
         & "the graph owns three member sets", nfail)
    call report(g % num_relations() .eq. 2, &
         & "and two relations", nfail)

    call report(g % holds_set(x) .and. g % holds_set(o) .and. &
         &      g % holds_set(p), &
         & "X, O and P are its own, by identity", nfail)

    call report(graph_holds_relation(g, flow), &
         & "the flow is owned, the same declared relation", nfail)
    call report(graph_holds_relation(g, d), &
         & "and so is the derived dependency", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! Signature closure: every slot of every owned relation resolves
  ! to a carrier the graph holds.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    class(member_set), allocatable :: dom
    integer                        :: k, s
    logical                        :: ok

    ok = .true.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. g % holds_set(dom)
       end do
    end do
    call report(ok, &
         & "every relation slot resolves to an owned carrier", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! The ternary relation remains first-class inside the graph: the
  ! owned flow answers arity three and its six tuples; the owned
  ! dependency answers its one. Nothing asked it to become binary.
  !===================================================================!

  subroutine check_ternary_stands(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer :: rp
    integer                  :: k

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % same_as(flow)) then
          call report(rp % arity() .eq. 3 .and. &
               &      rp % num_tuples() .eq. 6 .and. &
               &      rp % has([OP_TIMES, SLOT_E, PORT_OUT]), &
               & "the owned flow is still ternary, six tuples strong", &
               & nfail)
       else if (rp % same_as(d)) then
          call report(rp % num_tuples() .eq. 1 .and. &
               &      rp % has([OP_PLUS, OP_TIMES]), &
               & "the owned dependency still holds its one derived pair", &
               & nfail)
       end if
    end do

  end subroutine check_ternary_stands

  !===================================================================!
  ! Two relations of one signature are two citizens: D, and D
  ! restricted to the times-only subset - the empty relation over
  ! O x O - coexist in one graph with their own identities and
  ! their own extensions.
  !===================================================================!

  subroutine check_coexistence(nfail)

    integer, intent(inout) :: nfail

    type(subset_set)               :: times_only
    type(stored_relation)          :: d_empty
    type(relational_graph), target :: aux
    class(relation), pointer       :: r1, r2
    class(member_set), allocatable :: da, db
    logical                        :: ok

    times_only = subset_set('times-only', o, [OP_TIMES])
    d_empty    = restrict_slot(d, 1, times_only)

    call report(d_empty % num_tuples() .eq. 0, &
         & "restricting D to times-only leaves the empty relation", nfail)

    aux = relational_graph('coexistence', [held_set(o)], &
         & [held_relation(d), held_relation(d_empty)])

    r1 => aux % relation_at(1)
    r2 => aux % relation_at(2)

    call report(.not. r1 % same_as(r2), &
         & "one signature, two citizens - no collision", nfail)

    da = r1 % domain(1)
    db = r2 % domain(1)
    ok = da % same_as(db)
    da = r1 % domain(2)
    db = r2 % domain(2)
    ok = ok .and. da % same_as(db)
    call report(ok, &
         & "over the very same signature, slot for slot", nfail)

    call report(r1 % has([OP_PLUS, OP_TIMES]) .and. &
         &      r2 % num_tuples() .eq. 0, &
         & "and each keeps its own extension", nfail)

  end subroutine check_coexistence

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! the generators, as the generation rule asks.
  !===================================================================!

  logical function graph_holds_relation(g, r)

    type(relational_graph), target, intent(in) :: g
    class(relation)               , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % same_as(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program calculator_level_3
