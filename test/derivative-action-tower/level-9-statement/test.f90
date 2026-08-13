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
!      structure       R_flow, owned by G; D derived; order derived
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
! selector is located in G by identity and then DESTROYED, and
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
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : relational_graph, held_set, held_relation
  use graph_profile    , only : directed_adjacency_view
  use graph_algorithms , only : topological_order
  use class_graph_field, only : field
  use derivative_constitution_fixture, only : primal_execution, &
       &                                      tangent_action, reverse_action

  implicit none

  type(counted_set)                  :: v, o, p
  type(subset_set)                   :: x_dom, c, z_dom, p_in, p_out
  type(stored_relation), allocatable :: flow
  class(relation), allocatable       :: d
  class(relation), pointer           :: gflow => null()
  class(relation), pointer           :: rp    => null()
  type(relational_graph), target     :: g
  type(directed_adjacency_view)      :: view
  type(field)                        :: qx, zbar_f, grad_f, vseed_f
  class(member_set), allocatable     :: dom
  integer, allocatable               :: order(:)
  real(dp), allocatable              :: obs(:), gradient(:), gvals(:)
  real(dp), allocatable              :: seed1(:), vdir(:), base(:), dot(:)
  logical, allocatable               :: avail(:), davail(:)
  real(dp)                           :: jv, lhs, rhs
  integer                            :: table(3, 6)
  integer                            :: nfail, i

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 9 . statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure the statement selects
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
  allocate(flow)
  flow = stored_relation('flow', [o, v, p], table)

  p_in  = subset_set('input-ports', p, [PORT_IN1, PORT_IN2])
  p_out = subset_set('output-port', p, [PORT_OUT])

  d = compose_binary( &
       & project_slots(restrict_slot(flow, 3, p_out), [1, 2]), &
       & project_slots(restrict_slot(flow, 3, p_in ), [2, 1]))

  g = relational_graph('derivative specimen', &
       & [held_set(v), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])

  ! -- the order comes from the graph's own dependency: the view is
  !    made, the dependency selector dies, and the walk still answers
  view = directed_adjacency_view(g, d)
  deallocate(d)
  call topological_order(view, order)

  ! -- the statement keeps the GRAPH-OWNED flow, located by identity;
  !    then the external selector dies too
  do i = 1, g % num_relations()
     rp => g % relation_at(i)
     if (rp % same_as(flow)) gflow => rp
  end do
  if (.not. associated(gflow)) then
     error stop 'statement: the graph does not own the selected flow'
  end if
  deallocate(flow)

  allocate(base(v % size()), avail(v % size()))
  allocate(dot(v % size()), davail(v % size()))
  allocate(gradient(x_dom % size()))

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
  call grad_f % get_real_vector(gvals)
  write(*,'(1x,a)', advance='no') "DERIVATIVE_RESULT ="
  do i = 1, x_dom % size()
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

    qx = field('base point', x_dom)
    call qx % set_real_vector([3.0_dp, 2.0_dp])
    call qx % get_real_vector(obs)

    call primal_execution(gflow, v, x_dom, obs, c, order, base, avail)

    call report(abs(base(v % local_index(SLOT_U)) - 6.0_dp) &
         &      < 1.0d-12 .and. &
         &      abs(base(v % local_index(SLOT_Z)) - 9.0_dp) &
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

    zbar_f = field('reverse seed', z_dom)
    call zbar_f % set_real_vector([1.0_dp])
    call zbar_f % get_real_vector(seed1)

    call reverse_action(gflow, v, x_dom, order, base, z_dom, &
         & seed1, gradient)

    grad_f = field('derivative of z on X', x_dom)
    call grad_f % set_real_vector(gradient)

    call grad_f % domain(dom)
    call report(dom % same_as(x_dom), &
         & "the derivative is a field on X, by identity - not a " // &
         & "raw array, not a matrix", nfail)
    call report(abs(gradient(x_dom % local_index(SLOT_Y)) - 3.0_dp) &
         &      < 1.0d-12, &
         & "dz/dy = 3, read through the domain map", nfail)
    call report(abs(gradient(x_dom % local_index(SLOT_X)) - 3.0_dp) &
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
    call tangent_action(gflow, v, x_dom, [1.0_dp, 0.0_dp], c, order, &
         & base, dot, davail)
    jy = dot(v % local_index(SLOT_Z))

    call tangent_action(gflow, v, x_dom, [0.0_dp, 1.0_dp], c, order, &
         & base, dot, davail)
    jx = dot(v % local_index(SLOT_Z))

    call report(abs(jy - 3.0_dp) < 1.0d-12 .and. &
         &      abs(jx - 3.0_dp) < 1.0d-12, &
         & "forward probes: J e_y = 3 and J e_x = 3", nfail)
    call report(abs(jy - gradient(x_dom % local_index(SLOT_Y))) &
         &      < 1.0d-12 .and. &
         &      abs(jx - gradient(x_dom % local_index(SLOT_X))) &
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

    vseed_f = field('tangent direction', x_dom)
    call vseed_f % set_real_vector([-1.0_dp, 4.0_dp])
    call vseed_f % get_real_vector(vdir)

    call tangent_action(gflow, v, x_dom, vdir, c, order, base, &
         & dot, davail)
    jv = dot(v % local_index(SLOT_Z))

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
    type(field)           :: zbar2_f

    allocate(xbar2(x_dom % size()))

    zbar2_f = field('reverse seed two', z_dom)
    call zbar2_f % set_real_vector([2.0_dp])
    call zbar2_f % get_real_vector(seed2)

    call reverse_action(gflow, v, x_dom, order, base, z_dom, &
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
