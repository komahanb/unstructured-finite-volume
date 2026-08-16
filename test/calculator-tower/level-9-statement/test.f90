!=====================================================================!
! CALCULATOR TOWER . LEVEL 9 . THE STATEMENT
!
! The level asks the only question left: WHAT PROBLEM IS BEING
! ASKED. The statement is: evaluate (2 + 3) x 4. It SELECTS - it
! invents nothing:
!
!      structure       the calculator relational graph GAMMA
!      inputs          a = 2, b = 3, d = 4 on K = { d, a, b }
!      constitution    the Level-8 law table, reused, never redone
!      discretization  the Level-6 location relation L
!      requested       the output slot e
!
! The unknowns are declared U = { c, e } - the OPPOSITE enumeration
! of Levels 7 and 8 - so the requested e sits second in storage and
! any raw solution(1) habit fails. The residual operation is a
! test-local adapter that holds the GRAPH-OWNED flow (the external
! selector dies before the solve) and delegates every number to the
! Level-8 fixture; ordinary GMRES solves through its own operation
! face, rhs on Y in, solution on U out; and the answer is read
! through U's local_index alone. The literal 20 stands only in the
! final assertion. The seven-vertex host is compatibility scenery -
! not part of the mathematical statement.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module constituted_residual_fixture

  ! The Level-9 adapter: Level-8 semantics wearing the legacy
  ! graph_operation face so the ordinary solver can drive it. Its
  ! entire numerical act is delegation.

  use iso_fortran_env  , only : dp => REAL64
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_set_representation, only : set_representation
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relation   , only : relation
  ! The kernel's graph and the ordinary view's ordinary_graph are two
  ! types with two names now, so nothing is renamed at this door.
  use fractal_graph        , only : graph, known_branch, null_branch
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set
  use graph_operation_view, only : graph_operation
  use graph_ordinary_view, only : ordinary_graph
  use graph_field_calculus, only : graph_field
  use class_graph_field, only : field
  use arithmetic_constitution_fixture, only : generated_residual

  implicit none

  private
  public :: constituted_residual

  type, extends(graph_operation) :: constituted_residual
     class(relation), allocatable :: flow      ! the GRAPH-OWNED copy
     class(relation), allocatable :: located
     type(set_graph)            :: xs, os, ys
     type(set_graph)             :: known, unknown
     ! Identities, counts, and the action's OWN coordinates - copied
     ! in at construction, the way a CSR relation copies the numbering
     ! it executes against. A representation carries no identity and
     ! answers only where a member stands, so holding one keeps this
     ! type free of any caller's map.
     integer                      :: n_xs = 0, n_os = 0, n_ys = 0
     integer                      :: n_known = 0, n_unknown = 0
     class(set_representation), allocatable :: c_xs, c_os, c_ys
     class(set_representation), allocatable :: c_known, c_unknown
     real(dp), allocatable        :: known_values(:)
   contains
     procedure :: name   => cr_name
     procedure :: domain => cr_domain
     procedure :: apply  => cr_apply
  end type constituted_residual

  interface constituted_residual
     module procedure create_adapter
  end interface constituted_residual

contains

  !===================================================================!
  ! The selector only NAMES the relation; what the adapter keeps is
  ! the graph-owned citizen, found by identity and refused if GAMMA
  ! does not own it.
  !===================================================================!

  function create_adapter(g, b, selector, located, xs, os, ys, sets, &
       &                  known, known_values, unknown) result(this)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(relation)  , intent(in) :: selector, located
    type(set_graph), intent(in) :: xs, os, ys
    type(set_map)  , intent(in) :: sets
    type(set_graph) , intent(in) :: known, unknown
    real(dp)         , intent(in) :: known_values(:)
    type(constituted_residual)    :: this

    class(relation), pointer :: rp
    integer :: kk
    logical :: found
    integer         :: n_known
    integer         :: n_os
    integer         :: n_unknown
    integer         :: n_xs
    integer         :: n_ys

    found = .false.
    do kk = 1, num_relations(g)
       rp => relation_at(g, b, kk)
       if (rp % same_as(selector)) then
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       error stop 'statement: the graph does not own the selected flow'
    end if
    allocate(this % flow, source=rp)

    allocate(this % located, source=located)
    this % xs = xs ; this % n_xs = sets % size_of(xs)
    this % os = os ; this % n_os = sets % size_of(os)
    this % ys = ys ; this % n_ys = sets % size_of(ys)
    this % known   = known   ; this % n_known   = sets % size_of(known)
    this % unknown = unknown ; this % n_unknown = sets % size_of(unknown)

    call sets % extent_of(xs,      this % c_xs)
    call sets % extent_of(os,      this % c_os)
    call sets % extent_of(ys,      this % c_ys)
    call sets % extent_of(known,   this % c_known)
    call sets % extent_of(unknown, this % c_unknown)
    this % known_values = known_values

  end function create_adapter

  pure function cr_name(this) result(name)
    class(constituted_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted residual'
  end function cr_name

  subroutine cr_domain(this, input_graph, domain, nentries)
    class(constituted_residual), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    integer         :: n_ys
    associate (u1 => input_graph); end associate
    domain   = this % ys
    nentries = this % n_ys
  end subroutine cr_domain

  subroutine cr_apply(this, input_graph, input_data, output)

    class(constituted_residual), intent(in)        :: this
    class(ordinary_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_map) :: mine

    type(field)                    :: out
    type(set_graph) :: dom
    real(dp), allocatable          :: ustate(:), r(:)
    integer         :: n_os
    integer         :: n_xs
    integer         :: n_ys

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'statement: the residual needs a state to judge'
    end if
    if (size(input_data) < 1) then
       error stop 'statement: the residual needs a state to judge'
    end if

    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % unknown)) then
       error stop 'statement: the state must live on the unknown domain'
    end if
    call input_data(1) % get_real_vector(ustate)

    !----------------------------------------------------------------!
    ! The action's own coordinates, rebound into a LOCAL map for the
    ! duration of the call. The representations are its own; the map
    ! is a temporary that never leaves this scope.
    !----------------------------------------------------------------!

    call mine % bind(this % xs,      this % c_xs)
    call mine % bind(this % os,      this % c_os)
    call mine % bind(this % ys,      this % c_ys)
    call mine % bind(this % known,   this % c_known)
    call mine % bind(this % unknown, this % c_unknown)

    call generated_residual(this % flow, this % located, &
         & this % xs, this % os, this % ys, &
         & mine, this % known, this % known_values, &
         & this % unknown, ustate, r)

    out = field('residual', this % ys, this % n_ys)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine cr_apply

end module constituted_residual_fixture

program calculator_level_9

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_field_calculus, only : graph_field
  use graph_relation   , only : stored_relation, relation
  use fractal_graph        , only : graph, set_graph => graph, &
       & known_branch, null_branch
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_inclusion_map  , only : inclusion_map, declared_subobject
  use graph_relational_view, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use class_graph      , only : ordinary_stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use constituted_residual_fixture, only : constituted_residual

  implicit none

  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type(set_graph)     :: x, o, p, y
  type(set_graph)      :: k, u, p_out, p_in
  type(stored_relation), allocatable :: flow
  type(stored_relation) :: located, t_out3, t_in3, produces, consumes
  class(relation), allocatable       :: d
  type(graph)             , target :: g
  type(graph)             , target :: scell(3), selem(3)
  type(graph)             , target :: rcell(2), relem(2)
  type(relational_binding)         :: bnd
  integer                          :: kcell
  type(ordinary_stored_graph)    :: host
  type(constituted_residual) :: residual_op
  type(gmres)           :: solver
  type(field)           :: rhsf
  type(set_graph)  :: dom
  class(graph_field), allocatable :: sol, rf
  real(dp), allocatable :: gv(:), solval(:), rv(:)
  real(dp)              :: answer
  integer               :: table(3, 6)
  integer               :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions
  integer         :: n_dom

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure, and its graph
  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call o % declare()
  call sets % bind(o, counted_set_representation(2))
  call p % declare()
  call sets % bind(p, counted_set_representation(3))
  call y % declare()
  call sets % bind(y, counted_set_representation(2))

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  allocate(flow)
  flow = stored_relation('flow', [o, x, p], table, sets)

  call p_out % declare()
  call sets       % bind(p_out, listed_set_representation([PORT_OUT]))
  call inclusions % include_in(p_out, p)
  t_out3   = restrict_slot(flow, 3, p_out, sets, inclusions)
  call p_in % declare()
  call sets       % bind(p_in, listed_set_representation([PORT_IN1, PORT_IN2]))
  call inclusions % include_in(p_in, p)
  t_in3    = restrict_slot(flow, 3, p_in, sets, inclusions)
  produces = project_slots(t_out3, [1, 2], sets)
  consumes = project_slots(t_in3 , [2, 1], sets)
  d        = compose_binary(produces, consumes, sets)

  ! 'calculator': (S, P) as one sequence on each branch.
  call g % declare()
  do kcell = 1, 3
     call scell(kcell) % declare()
     call selem(kcell) % declare()
  end do
  do kcell = 1, 2
     call rcell(kcell) % declare()
     call relem(kcell) % declare()
  end do

  call bnd % bind_set(selem(1), x)
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

  ! -- the discretization, distinct from the graph on purpose
  located = stored_relation('located', [y, x], &
       & reshape([ROW_C, SLOT_C,  ROW_E, SLOT_E], [2, 2]), sets)

  ! -- the inputs, and the unknowns in a NEW enumeration {c, e}
  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_D, SLOT_A, SLOT_B]))
  call inclusions % include_in(k, x)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_C, SLOT_E]))
  call inclusions % include_in(u, x)

  ! -- the adapter keeps the GRAPH-OWNED flow; the selector dies.
  residual_op = constituted_residual(g, bnd, flow, located, x, o, y, sets, &
       &                             k, [4.0_dp, 2.0_dp, 3.0_dp], u)
  deallocate(flow)

  ! -- the compatibility scenery: not part of the statement.
  host = ordinary_stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  ! -- the ordinary solver, through its own operation face
  call solver % attach(residual_op, host, u, sets % size_of(u))
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom, n_dom)
  call report(dom % same_as(u), &
       & "the solver answers on U = { c, e }", nfail)

  call solver % constant(gv)
  rhsf = field('rhs', y, sets % size_of(y))
  call rhsf % set_real_vector(-gv)

  call solver % apply(host, [rhsf], sol)

  dom = sol % domain()
  call report(dom % same_as(u), &
       & "and the solution field lives there, by identity", nfail)

  ! -- verify the constituted residual vanishes at the solution
  select type (sol)
  type is (field)
     call residual_op % apply(host, [sol], rf)
  end select
  call rf % get_real_vector(rv)
  call report(sqrt(sum(rv * rv)) < 1.0d-9, &
       & "the constituted residual vanishes at the solution", nfail)

  ! -- the statement's one question
  call sol % get_real_vector(solval)
  answer = solval(sets % index_in(u, SLOT_E))

  call report(abs(answer - 20.0_dp) < 1.0d-9, &
       & "evaluate (2 + 3) x 4: the tower answers 20", nfail)

  write(*,'(1x,a,i0)') "CALCULATOR_RESULT = ", nint(answer)

  call verdict(nfail, "level 9")

end program calculator_level_9
