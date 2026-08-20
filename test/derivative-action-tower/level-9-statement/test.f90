!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 9 . THE STATEMENT
!
! Gate C. The level asks the only question left: WHAT COMPLETE
! DIFFERENTIATION PROBLEM IS BEING ASKED. The statement is: given
! the relational computation u = product(x, y), z = sum(u, y),
! constituted by product(a,b) = ab and sum(a,b) = a + b, evaluate
! the primal and the derivative of the response z at x = 2, y = 3.
! It SELECTS - it invents nothing:
!
!      structure       T_flow, owned by GAMMA; D derived; order derived
!      base point      x = 2, y = 3 on X = { y, x }
!      constitution    the Level-8 laws and the ONE local
!                      linearization, reused, never redone
!      response        Z = { z }, a subdomain - no location relation
!      requested       the derivative of z with respect to X
!
! The answer is not a scalar. It is a FIELD on the independent
! domain,
!
!      Dz^T(1) = [3, 3]  on  X = { y, x }
!
! obtained by ONE reverse traversal with zbar = 1, and certified
! forward one action at a time: the y-basis seed and the x-basis
! seed each answer 3, the nontrivial direction v = [-1, 4] answers
! Jv = 9, the returned field acts as the derivative linear
! functional on that same v - and the pairing closes at 18 = 18
! through this whole path. The primal z = 9 stands beside it as a
! secondary truth.
!
! The statement evaluates GRAPH-OWNED structure: the external flow
! selector is located in GAMMA by identity and then DESTROYED, and
! every number below is computed through the graph's own relation.
! Level 7 - minimization - is not on this road at all: this orbit
! does not inhabit that radial contract, and no solver, no host and
! no residual appear anywhere. Nothing here is an adjoint: reverse
! derivative action is not the adjoint method, and that distinction
! is the next tower's to earn.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_9

  use iso_fortran_env  , only : dp => REAL64
  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary   , only : stored_relation, relation
  use relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use relation_algorithms , only : topological_order
  use field_stored, only : stored_field
  use derivative_constitution_fixture, only : primal_execution, &
       &                                      tangent_action, reverse_action
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & has_set

  implicit none


  type(graph)                  :: v, o, p
  type(graph)                   :: x_dom, c, z_dom, p_in, p_out
  type(stored_relation), allocatable :: flow
  class(relation), allocatable       :: d
  class(relation), pointer           :: gflow => null()
  class(relation), pointer           :: rp    => null()
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(stored_field)                        :: qx, zbar_f, grad_f, vseed_f
  type(graph)     :: dom
  integer, allocatable               :: order(:)
  real(dp), allocatable              :: obs(:), gradient(:), gvals(:)
  real(dp), allocatable              :: seed1(:), vdir(:), base(:), dot(:)
  logical, allocatable               :: avail(:), davail(:)
  real(dp)                           :: jv, lhs, rhs
  integer                            :: table(3, 6)
  integer                            :: nfail, i
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 9 . statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure the statement selects
  call v % declare()
  call sets % bind(v, counted_set_representation(4))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))
  call x_dom % declare()
  call sets       % bind(x_dom, listed_set_representation([SLOT_Y, SLOT_X]))
  call inclusions % include_in(x_dom, v)
  call c % declare()
  call sets       % bind(c, listed_set_representation([SLOT_U, SLOT_Z]))
  call inclusions % include_in(c, v)
  call z_dom % declare()
  call sets       % bind(z_dom, listed_set_representation([SLOT_Z]))
  call inclusions % include_in(z_dom, c)

  table(:, 1) = [OP_PRODUCT, SLOT_X, PORT_IN1]
  table(:, 2) = [OP_PRODUCT, SLOT_Y, PORT_IN2]
  table(:, 3) = [OP_PRODUCT, SLOT_U, PORT_OUT]
  table(:, 4) = [OP_SUM    , SLOT_U, PORT_IN1]
  table(:, 5) = [OP_SUM    , SLOT_Y, PORT_IN2]
  table(:, 6) = [OP_SUM    , SLOT_Z, PORT_OUT]
  allocate(flow)
  flow = stored_relation('flow', [o, v, p], table, sets)

  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)

  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out, sets, inclusions), [1, 2], sets), &
       & project_slots(restrict_slot(flow, 3, p_in , sets, inclusions), [2, 1], sets), sets)

  ! 'derivative specimen': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 2
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), v)
  call bnd % bind_set(selem(2), o)
  call bnd % bind_set(selem(3), p)
  call bnd % bind_relation(relem(1), flow)
  call bnd % bind_relation(relem(2), d)

  do kcell = 1, 3
     scell(kcell) % branch(1) = known_branch(selem(kcell))
     if (kcell .lt. 3) scell(kcell) % branch(2) = &
          & known_branch(scell(kcell + 1))
  end do
  do kcell = 1, 2
     rcell(kcell) % branch(1) = known_branch(relem(kcell))
     if (kcell .lt. 2) rcell(kcell) % branch(2) = &
          & known_branch(rcell(kcell + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  ! -- the order comes from the graph's own dependency: the relation is
  !    made, the dependency selector dies, and the walk still answers
  call topological_order(d, sets, order)

  ! -- the statement keeps the GRAPH-OWNED flow, located by identity;
  !    then the external selector dies too
  do i = 1, num_relations(g)
     rp => relation_at(g, bnd, i)
     if (rp % same_as(flow)) gflow => rp
  end do
  if (.not. associated(gflow)) then
     error stop 'statement: the graph does not own the selected flow'
  end if
  deallocate(flow)

  allocate(base(sets % num_members_of(v)), avail(sets % num_members_of(v)))
  allocate(dot(sets % num_members_of(v)), davail(sets % num_members_of(v)))
  allocate(gradient(sets % num_members_of(x_dom)))

  call check_selection(nfail)
  call check_primal(nfail)
  call check_derivative_field(nfail)
  call check_forward_agreement(nfail)
  call check_functional_action(nfail)
  call check_duality(nfail)

  ! -- the tower's one result: a FIELD, not a scalar - and a
  !    REAL-valued field, serialized as it stands. One token per
  !    member, in X's own enumeration order, at full round-trip
  !    precision: the marker carries the actual field, never a
  !    rounded image of it. No integer conversion lives on this path.
  call grad_f % real_vector(gvals)
  write(*,'(1x,a)', advance='no') "DERIVATIVE_RESULT ="
  do i = 1, sets % num_members_of(x_dom)
     write(*,'(es24.16)', advance='no') gvals(i)
  end do
  write(*,'(a)') ""

  call verdict(nfail, "level 9")

contains

  !===================================================================!
  ! What the statement selected: a derived order, a graph-owned
  ! flow, and no surviving selector to have leaned on.
  !===================================================================!

  subroutine check_selection(nfail)

    integer, intent(inout) :: nfail

    call report(size(order) .eq. 2 .and. &
         &      order(1) .eq. OP_PRODUCT .and. order(2) .eq. OP_SUM, &
         & "the order is derived from the graph's own dependency, " // &
         & "the selector already dead", nfail)
    call report(.not. allocated(flow) .and. associated(gflow) .and. &
         &      gflow % arity() .eq. 3 .and. &
         &      gflow % num_tuples() .eq. 6, &
         & "the statement holds the GRAPH-OWNED flow: six ternary " // &
         & "facts, no external selector left alive", nfail)

  end subroutine check_selection

  !===================================================================!
  ! The primal at the statement's base point - through the
  ! graph-owned structure and the Level-8 constitution.
  !===================================================================!

  subroutine check_primal(nfail)

    integer, intent(inout) :: nfail

    qx = stored_field('base point', x_dom, sets % num_members_of(x_dom))
    call qx % set_real_vector([3.0_dp, 2.0_dp])
    call qx % real_vector(obs)

    call primal_execution(gflow, v, sets, x_dom, obs, c, order, base, avail)

    call report(abs(base(sets % index_in(v, SLOT_U)) - 6.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(base(sets % index_in(v, SLOT_Z)) - 9.0_dp) &
         &      < 1.0d-12, &
         & "the primal answers u = 6 and z = 9", nfail)

  end subroutine check_primal

  !===================================================================!
  ! THE statement's answer: one reverse traversal with zbar = 1
  ! returns the whole derivative as an ordinary FIELD on X - the
  ! natural first derivative object of this framework.
  !===================================================================!

  subroutine check_derivative_field(nfail)

    integer, intent(inout) :: nfail

    zbar_f = stored_field('reverse seed', z_dom, sets % num_members_of(z_dom))
    call zbar_f % set_real_vector([1.0_dp])
    call zbar_f % real_vector(seed1)

    call reverse_action(gflow, v, sets, x_dom, order, base, z_dom, &
         & seed1, gradient)

    grad_f = stored_field('derivative of z on X', x_dom, sets % num_members_of(x_dom))
    call grad_f % set_real_vector(gradient)

    dom = grad_f % domain()
    call report(dom % same_as(x_dom), &
         & "the derivative is a field on X, by identity - not a " // &
         & "raw array, not a matrix", nfail)
    call report(abs(gradient(sets % index_in(x_dom, SLOT_Y)) - 3.0_dp) &
         &      < 1.0d-12, &
         & "dz/dy = 3, read through the domain map", nfail)
    call report(abs(gradient(sets % index_in(x_dom, SLOT_X)) - 3.0_dp) &
         &      < 1.0d-12, &
         & "dz/dx = 3", nfail)
    call report(grad_f % num_entries() .eq. 2, &
         & "one reverse traversal, the complete row", nfail)

  end subroutine check_derivative_field

  !===================================================================!
  ! Forward certification of the same linear map, one action per
  ! seed: the y basis and the x basis each answer 3, agreeing with
  ! the components the single reverse traversal returned. Two
  ! forward actions were needed; one reverse action sufficed - for
  ! THIS specimen, and no complexity claim beyond it.
  !===================================================================!

  subroutine check_forward_agreement(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: jy, jx

    ! X = { y, x }: the y basis is [1, 0], the x basis is [0, 1]
    call tangent_action(gflow, v, sets, x_dom, [1.0_dp, 0.0_dp], c, order, &
         & base, dot, davail)
    jy = dot(sets % index_in(v, SLOT_Z))

    call tangent_action(gflow, v, sets, x_dom, [0.0_dp, 1.0_dp], c, order, &
         & base, dot, davail)
    jx = dot(sets % index_in(v, SLOT_Z))

    call report(abs(jy - 3.0_dp) < 1.0d-12 .and. &
         &      abs(jx - 3.0_dp) < 1.0d-12, &
         & "forward probes: J e_y = 3 and J e_x = 3", nfail)
    call report(abs(jy - gradient(sets % index_in(x_dom, SLOT_Y))) &
         &      < 1.0d-12 .and. &
         &      abs(jx - gradient(sets % index_in(x_dom, SLOT_X))) &
         &      < 1.0d-12, &
         & "and each agrees with the component the one reverse " // &
         & "pass returned", nfail)

  end subroutine check_forward_agreement

  !===================================================================!
  ! The returned field ACTS: on the nontrivial direction
  ! v = [-1, 4] the tangent action answers 9, and the derivative
  ! field paired with the same v answers 9 - the field is the
  ! derivative linear functional, and no matrix was built to say so.
  !===================================================================!

  subroutine check_functional_action(nfail)

    integer, intent(inout) :: nfail

    vseed_f = stored_field('tangent direction', x_dom, sets % num_members_of(x_dom))
    call vseed_f % set_real_vector([-1.0_dp, 4.0_dp])
    call vseed_f % real_vector(vdir)

    call tangent_action(gflow, v, sets, x_dom, vdir, c, order, base, &
         & dot, davail)
    jv = dot(sets % index_in(v, SLOT_Z))

    call report(abs(jv - 9.0_dp) < 1.0d-12, &
         & "Jv = 9 for v = [-1, 4]", nfail)
    call report(abs(dot_product(gradient, vdir) - jv) < 1.0d-12, &
         & "and the derivative field paired with v says the same: " // &
         & "it is the linear functional", nfail)

  end subroutine check_functional_action

  !===================================================================!
  ! The pairing law, once more and through the WHOLE statement
  ! path: graph-owned structure, derived order, reused constitution.
  !===================================================================!

  subroutine check_duality(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: xbar2(:), seed2(:)
    type(stored_field)           :: zbar2_f

    allocate(xbar2(sets % num_members_of(x_dom)))

    zbar2_f = stored_field('reverse seed two', z_dom, sets % num_members_of(z_dom))
    call zbar2_f % set_real_vector([2.0_dp])
    call zbar2_f % real_vector(seed2)

    call reverse_action(gflow, v, sets, x_dom, order, base, z_dom, &
         & seed2, xbar2)

    lhs = seed2(1) * jv
    rhs = dot_product(xbar2, vdir)

    call report(abs(lhs - rhs) < 1.0d-12, &
         & "< zbar, Jv >_Z = < J^T zbar, v >_X, statement-wide", &
         & nfail)
    call report(abs(lhs - 18.0_dp) < 1.0d-12, &
         & "and both say 18", nfail)

  end subroutine check_duality

end program derivative_level_9
