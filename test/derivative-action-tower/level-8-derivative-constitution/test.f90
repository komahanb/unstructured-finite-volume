!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 8 . DERIVATIVE CONSTITUTION
!
! Gate B. The level answers one question: CAN NUMERICAL TANGENT AND
! REVERSE ACTIONS TRAVEL OVER THE SAME CONSTITUTED COMPUTATION
! STRUCTURE. Only here do the symbols mean something - product
! multiplies, sum adds - and each operation owes exactly one more
! thing: its local linearization at the primal inputs. ONE shadow,
! TWO actions:
!
!      Jv        forward, along the derived order
!      J^T zbar  backward, along the same order reversed - the SAME
!                coefficients, ends swapped, accumulated += per
!                operation/input-port incidence
!
! sealed by the global duality law
!
!      < zbar, Jv >_Z  =  < J^T zbar, v >_X
!
! at the primary base (18 = 18), at an asymmetric secondary base
! (-38 = -38), under unit seeds, and under reversed tuple order.
! Seeds and results are ORDINARY fields on ordinary domains: no
! tangent_field, no cotangent_field, no dual number. The evaluators
! consume R_flow incidences, the derived order, primal values and
! the law table; the structural J_ZX pattern is consulted by the
! TEST as support metadata - never fed to the numerical path. And
! the accumulation is certified, not assumed: y receives exactly
! two incidence contributions - sum.in2 and product.in2 - counted
! by generic machinery that names no slot.
!
! No adjoint solve, no minimization, no full Jacobian, no finite
! difference, no path counting - anywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_8

  use iso_fortran_env  , only : dp => REAL64
  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : relational_graph, held_set, held_relation
  use graph_profile    , only : directed_adjacency_view
  use graph_algorithms , only : reachable, topological_order
  use class_graph_field, only : field
  use derivative_constitution_fixture, only : apply_law, &
       & local_linearization, slot_for_port, primal_execution, &
       & tangent_action, reverse_action

  implicit none

  type(counted_set)              :: v, o, p
  type(subset_set)               :: x_dom, c, z_dom, p_in, p_out
  type(stored_relation)          :: flow, backwards
  class(relation), allocatable   :: d, d2, av
  type(relational_graph), target :: g, g2, g_av
  type(directed_adjacency_view)  :: view, view2, dep_view
  type(field)                    :: qx, vx, zbar_f, jv_f, xbar_f
  integer, allocatable           :: order(:), order2(:), hits(:)
  real(dp), allocatable          :: obs(:), seedv(:), zseed(:), xbar(:)
  real(dp), allocatable          :: base(:), dot(:)
  logical, allocatable           :: avail(:), davail(:)
  real(dp)                       :: jv
  integer                        :: table(3, 6)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 8 . action"
  write(*,'(1x,a)') "============================================="

  v     = counted_set('value-slots', 4)
  o     = counted_set('operations' , 2)
  p     = counted_set('ports'      , 3)
  x_dom = subset_set('independent', v, [SLOT_Y, SLOT_X])
  c     = subset_set('computed'   , v, [SLOT_U, SLOT_Z])
  z_dom = subset_set('response'   , c, [SLOT_Z])

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  flow = stored_relation('flow', [o, v, p], table)

  p_in  = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  p_out = subset_set('output-port', p, [PORT_OUT])

  ! The lower road, walked again: the structure supplies the order.
  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out), [1, 2]), &
       & project_slots(restrict_slot(flow, 3, p_in ), [2, 1]))
  g = relational_graph('derivative specimen', &
       & [held_set(v), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])
  view = directed_adjacency_view(g, d)
  call topological_order(view, order)

  allocate(base(v % size()), avail(v % size()))
  allocate(dot(v % size()), davail(v % size()))
  allocate(xbar(x_dom % size()))

  call check_derived_order(nfail)
  call check_laws(nfail)
  call check_linearization(nfail)
  call check_local_pairing(nfail)
  call check_primal(nfail)
  call check_forward_action(nfail)
  call check_reverse_action(nfail)
  call check_accumulation(nfail)
  call check_duality_primary(nfail)
  call check_support_versus_action(nfail)
  call check_secondary_base(nfail)
  call check_order_invariance(nfail)

  call verdict(nfail, "level 8")

contains

  !===================================================================!
  ! The order is DERIVED; the evaluators receive what the sort
  ! answered - [product, sum] appears as expectation only.
  !===================================================================!

  subroutine check_derived_order(nfail)

    integer, intent(inout) :: nfail

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PRODUCT .and. order(2) .eq. OP_SUM, &
         & "the derived order is [product, sum], exactly", nfail)

  end subroutine check_derived_order

  !===================================================================!
  ! The primal law table, certified independently.
  !===================================================================!

  subroutine check_laws(nfail)

    integer, intent(inout) :: nfail

    call report(abs(apply_law(OP_PRODUCT, 2.0_dp, 3.0_dp) - 6.0_dp) &
         &      < 1.0d-14, &
         & "product is multiplication: product(2, 3) = 6", nfail)
    call report(abs(apply_law(OP_SUM, 6.0_dp, 3.0_dp) - 9.0_dp) &
         &      < 1.0d-14, &
         & "sum is addition: sum(6, 3) = 9", nfail)

  end subroutine check_laws

  !===================================================================!
  ! The local linear shadow: port-relative coefficients at the
  ! primal inputs - and a nonsymmetric point so a swapped port
  ! would be caught.
  !===================================================================!

  subroutine check_linearization(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: cf(2)

    cf = local_linearization(OP_PRODUCT, 2.0_dp, 3.0_dp)
    call report(abs(cf(1) - 3.0_dp) < 1.0d-14 .and. &
         &      abs(cf(2) - 2.0_dp) < 1.0d-14, &
         & "L_product(2, 3) = [3, 2]: each port carries the OTHER " // &
         & "input", nfail)
    cf = local_linearization(OP_PRODUCT, 5.0_dp, 7.0_dp)
    call report(abs(cf(1) - 7.0_dp) < 1.0d-14 .and. &
         &      abs(cf(2) - 5.0_dp) < 1.0d-14, &
         & "L_product(5, 7) = [7, 5]: nonsymmetric, port-true", nfail)
    cf = local_linearization(OP_SUM, 5.0_dp, 7.0_dp)
    call report(abs(cf(1) - 1.0_dp) < 1.0d-14 .and. &
         &      abs(cf(2) - 1.0_dp) < 1.0d-14, &
         & "L_sum = [1, 1] wherever it stands", nfail)

  end subroutine check_linearization

  !===================================================================!
  ! The local pairing law, per operation at a nonsymmetric point:
  ! < obar, L v > = < L^T obar, v >. Both actions read ONE row -
  ! there is no separate reverse formula to disagree with.
  !===================================================================!

  subroutine check_local_pairing(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: cf(2), obar, v1, v2, lhs, rhs

    obar = 1.3_dp
    v1   = 0.2_dp
    v2   = -0.7_dp

    cf  = local_linearization(OP_PRODUCT, 5.0_dp, 7.0_dp)
    lhs = obar * (cf(1) * v1 + cf(2) * v2)
    rhs = (cf(1) * obar) * v1 + (cf(2) * obar) * v2
    call report(abs(lhs - rhs) < 1.0d-13, &
         & "product: <obar, Lv> = <L^T obar, v> at (5, 7)", nfail)

    cf  = local_linearization(OP_SUM, 5.0_dp, 7.0_dp)
    lhs = obar * (cf(1) * v1 + cf(2) * v2)
    rhs = (cf(1) * obar) * v1 + (cf(2) * obar) * v2
    call report(abs(lhs - rhs) < 1.0d-13, &
         & "sum: the same pairing, the same one shadow", nfail)

  end subroutine check_local_pairing

  !===================================================================!
  ! The primal at the primary base: seed [3, 2] on X = { y, x },
  ! execute the derived order, and the computed slots earn their
  ! numbers - u = 6, z = 9.
  !===================================================================!

  subroutine check_primal(nfail)

    integer, intent(inout) :: nfail

    qx = field('base point', x_dom)
    call qx % set_real_vector([3.0_dp, 2.0_dp])
    call qx % get_real_vector(obs)

    call primal_execution(flow, v, x_dom, obs, c, order, base, avail)

    call report(abs(base(v % local_index(SLOT_U)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "u = 6, computed into the workspace", nfail)
    call report(abs(base(v % local_index(SLOT_Z)) - 9.0_dp) &
         &      < 1.0d-12, &
         & "z = 9: the base point stands complete", nfail)
    call report(count(avail) .eq. v % size(), &
         & "every slot is available after execution", nfail)

  end subroutine check_primal

  !===================================================================!
  ! The tangent action at the primary base: seed dy = -1, dx = 4 as
  ! an ORDINARY field on X, ride the derived order, and Jv = 9
  ! comes back as an ordinary field on Z.
  !===================================================================!

  subroutine check_forward_action(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: dom
    real(dp), allocatable          :: out_val(:)

    vx = field('tangent seed', x_dom)
    call vx % set_real_vector([-1.0_dp, 4.0_dp])
    call vx % get_real_vector(seedv)

    call tangent_action(flow, v, x_dom, seedv, c, order, base, &
         & dot, davail)

    call report(abs(dot(v % local_index(SLOT_U)) - 10.0_dp) &
         &      < 1.0d-12, &
         & "du = 3(4) + 2(-1) = 10, by the local shadow", nfail)

    jv = dot(v % local_index(SLOT_Z))
    jv_f = field('tangent response', z_dom)
    call jv_f % set_real_vector([jv])

    call jv_f % domain(dom)
    call jv_f % get_real_vector(out_val)
    call report(dom % same_as(z_dom) .and. &
         &      abs(out_val(z_dom % local_index(SLOT_Z)) - 9.0_dp) &
         &      < 1.0d-12, &
         & "Jv = 9, an ordinary field on Z", nfail)

  end subroutine check_forward_action

  !===================================================================!
  ! The reverse action at the primary base: seed zbar = 2 as an
  ! ordinary field on Z, ride the same order backwards with the
  ! same coefficients, and J^T zbar = [6, 6] lands on X = { y, x }.
  !===================================================================!

  subroutine check_reverse_action(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: dom

    zbar_f = field('reverse seed', z_dom)
    call zbar_f % set_real_vector([2.0_dp])
    call zbar_f % get_real_vector(zseed)

    call reverse_action(flow, v, x_dom, order, base, z_dom, &
         & zseed, xbar, hits)

    call report(abs(xbar(x_dom % local_index(SLOT_Y)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "ybar = 2 + 2(2) = 6: two contributions, one slot", nfail)
    call report(abs(xbar(x_dom % local_index(SLOT_X)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "xbar = 3(2) = 6", nfail)

    xbar_f = field('reverse result', x_dom)
    call xbar_f % set_real_vector(xbar)
    call xbar_f % domain(dom)
    call report(dom % same_as(x_dom), &
         & "J^T zbar is an ordinary field on X, by identity", nfail)

  end subroutine check_reverse_action

  !===================================================================!
  ! Accumulation, certified - not assumed: the generic counter
  ! (which names no slot) shows y received exactly TWO INPUT-PORT
  ! INCIDENCE ACCUMULATION EVENTS - sum.in2 and product.in2 - while
  ! x and u received one each and the response none.
  !
  ! An event is a traversal fact, counted whenever bar(in) += ...
  ! executes; whether the added value is nonzero is the business of
  ! the coefficients and the seed. Incidence multiplicity is NOT
  ! numerical nonzero multiplicity. (At this base both of y's events
  ! do carry nonzero value - 2 and 4 - which is why ybar reads 6;
  ! that follows from L and zbar, not from the count.)
  !===================================================================!

  subroutine check_accumulation(nfail)

    integer, intent(inout) :: nfail

    call report(hits(v % local_index(SLOT_Y)) .eq. 2, &
         & "y took two input-port incidence accumulation events", &
         & nfail)
    call report(hits(v % local_index(SLOT_X)) .eq. 1 .and. &
         &      hits(v % local_index(SLOT_U)) .eq. 1 .and. &
         &      hits(v % local_index(SLOT_Z)) .eq. 0, &
         & "x and u one event each, the response none: " // &
         & "incidence-local +=, never path counting", nfail)

  end subroutine check_accumulation

  !===================================================================!
  ! The global duality law at the primary base: both pairings
  ! computed from the returned actions, equal to tight tolerance -
  ! and equal to 18, pinned in the assertion alone.
  !===================================================================!

  subroutine check_duality_primary(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: lhs, rhs

    lhs = zseed(1) * jv
    rhs = dot_product(xbar, seedv)

    call report(abs(lhs - rhs) < 1.0d-12, &
         & "< zbar, Jv >_Z = < J^T zbar, v >_X", nfail)
    call report(abs(lhs - 18.0_dp) < 1.0d-12, &
         & "and both say 18", nfail)

  end subroutine check_duality_primary

  !===================================================================!
  ! Support and action, complementary: the reachability that built
  ! Gate A's J-pattern agrees with where the reverse action landed
  ! nonzero incidence - yet the evaluators never consumed A_V, J_ZX
  ! or its transpose: their signatures take R_flow, order, values
  ! and laws, nothing more. Support says WHO; incidence says HOW.
  !===================================================================!

  subroutine check_support_versus_action(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m
    logical :: ok

    av = compose_binary( &
         & project_slots(restrict_slot(flow, 3, p_in ), [2, 1]), &
         & project_slots(restrict_slot(flow, 3, p_out), [1, 2]))
    g_av = relational_graph('value dependency', [held_set(v)], &
         & [held_relation(av)])
    dep_view = directed_adjacency_view(g_av, av)

    ok = .true.
    do i = 1, x_dom % size()
       m = x_dom % member(i)
       ok = ok .and. (reachable(dep_view, m, SLOT_Z) .eqv. &
            &         (hits(v % local_index(m)) .gt. 0))
    end do
    call report(ok, &
         & "the structural support agrees with the action's " // &
         & "incidence - and was never fed to it", nfail)

  end subroutine check_support_versus_action

  !===================================================================!
  ! The asymmetric second base, same everything: x = 4, y = 2 makes
  ! dz/dx = 2 and dz/dy = 5 differ, so a swapped slot cannot hide.
  ! Same flow, same order, same laws, same evaluators.
  !===================================================================!

  subroutine check_secondary_base(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: obs2(:), base2(:), dot2(:), xbar2(:)
    logical , allocatable :: avail2(:), davail2(:)
    real(dp)              :: jv2, lhs, rhs
    type(field)           :: qx2

    allocate(base2(v % size()), avail2(v % size()))
    allocate(dot2(v % size()), davail2(v % size()))
    allocate(xbar2(x_dom % size()))

    qx2 = field('base point two', x_dom)
    call qx2 % set_real_vector([2.0_dp, 4.0_dp])
    call qx2 % get_real_vector(obs2)

    call primal_execution(flow, v, x_dom, obs2, c, order, &
         & base2, avail2)
    call report(abs(base2(v % local_index(SLOT_U)) - 8.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(base2(v % local_index(SLOT_Z)) - 10.0_dp) &
         &      < 1.0d-12, &
         & "second base: u = 8, z = 10", nfail)

    ! dy = 5, dx = -3 on X = { y, x }
    call tangent_action(flow, v, x_dom, [5.0_dp, -3.0_dp], c, &
         & order, base2, dot2, davail2)
    jv2 = dot2(v % local_index(SLOT_Z))
    call report(abs(dot2(v % local_index(SLOT_U)) - 14.0_dp) &
         &      < 1.0d-12 .and. abs(jv2 - 19.0_dp) < 1.0d-12, &
         & "du = 2(-3) + 4(5) = 14 and Jv = 19", nfail)

    call reverse_action(flow, v, x_dom, order, base2, z_dom, &
         & [-2.0_dp], xbar2)
    call report(abs(xbar2(x_dom % local_index(SLOT_Y)) + 10.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(xbar2(x_dom % local_index(SLOT_X)) + 4.0_dp) &
         &      < 1.0d-12, &
         & "J^T(-2) = [-10, -4] on { y, x }", nfail)

    lhs = (-2.0_dp) * jv2
    rhs = dot_product(xbar2, [5.0_dp, -3.0_dp])
    call report(abs(lhs - rhs) < 1.0d-12 .and. &
         &      abs(lhs + 38.0_dp) < 1.0d-12, &
         & "duality at the asymmetric base: -38 = -38", nfail)

    ! Unit tangent seeds: independent action probes, never an
    ! assembled Jacobian.
    call tangent_action(flow, v, x_dom, [0.0_dp, 1.0_dp], c, &
         & order, base2, dot2, davail2)
    call report(abs(dot2(v % local_index(SLOT_Z)) - 2.0_dp) &
         &      < 1.0d-12, &
         & "x-basis seed: Jv = 2 = dz/dx", nfail)

    call tangent_action(flow, v, x_dom, [1.0_dp, 0.0_dp], c, &
         & order, base2, dot2, davail2)
    call report(abs(dot2(v % local_index(SLOT_Z)) - 5.0_dp) &
         &      < 1.0d-12, &
         & "y-basis seed: Jv = 5 = dz/dy", nfail)

    call reverse_action(flow, v, x_dom, order, base2, z_dom, &
         & [1.0_dp], xbar2)
    call report(abs(xbar2(x_dom % local_index(SLOT_Y)) - 5.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(xbar2(x_dom % local_index(SLOT_X)) - 2.0_dp) &
         &      < 1.0d-12, &
         & "unit reverse seed: J^T(1) = [5, 2] - the gradient row, " // &
         & "one backward pass", nfail)

  end subroutine check_secondary_base

  !===================================================================!
  ! The six facts handed backwards: their OWN derivation - D2,
  ! order2, ports rediscovered - drives the same primal, the same
  ! actions, the same duality. Meaning is order-independent.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    integer               :: rev(3, 6), k
    real(dp), allocatable :: baseb(:), dotb(:), xbarb(:)
    logical , allocatable :: availb(:), davailb(:)
    real(dp)              :: jvb, lhs, rhs

    allocate(baseb(v % size()), availb(v % size()))
    allocate(dotb(v % size()), davailb(v % size()))
    allocate(xbarb(x_dom % size()))

    do k = 1, 6
       rev(:, k) = table(:, 7 - k)
    end do
    backwards = stored_relation('flow backwards', [o, v, p], rev)

    d2 = compose_binary( &
         & project_slots(restrict_slot(backwards, 3, p_out), [1, 2]), &
         & project_slots(restrict_slot(backwards, 3, p_in ), [2, 1]))
    g2 = relational_graph('derivative specimen backwards', &
         & [held_set(v), held_set(o), held_set(p)], &
         & [held_relation(backwards), held_relation(d2)])
    view2 = directed_adjacency_view(g2, d2)
    call topological_order(view2, order2)

    call primal_execution(backwards, v, x_dom, obs, c, order2, &
         & baseb, availb)
    call tangent_action(backwards, v, x_dom, seedv, c, order2, &
         & baseb, dotb, davailb)
    jvb = dotb(v % local_index(SLOT_Z))
    call reverse_action(backwards, v, x_dom, order2, baseb, &
         & z_dom, zseed, xbarb)

    lhs = zseed(1) * jvb
    rhs = dot_product(xbarb, seedv)

    call report(abs(baseb(v % local_index(SLOT_U)) - 6.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(baseb(v % local_index(SLOT_Z)) - 9.0_dp) &
         &      < 1.0d-12, &
         & "backwards facts, same primal: u = 6, z = 9", nfail)
    call report(abs(jvb - 9.0_dp) < 1.0d-12, &
         & "same tangent action: Jv = 9", nfail)
    call report(abs(xbarb(x_dom % local_index(SLOT_Y)) - 6.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(xbarb(x_dom % local_index(SLOT_X)) - 6.0_dp) &
         &      < 1.0d-12, &
         & "same reverse action: [6, 6]", nfail)
    call report(abs(lhs - rhs) < 1.0d-12, &
         & "and the duality holds, either way written", nfail)

  end subroutine check_order_invariance

end program derivative_level_8
