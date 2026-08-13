!=====================================================================!
! CALCULATOR TOWER . LEVEL 9 . THE STATEMENT
!
! The level asks the only question left: WHAT PROBLEM IS BEING
! ASKED. The statement is: evaluate (2 + 3) x 4. It SELECTS - it
! invents nothing:
!
!      structure       the calculator relational graph G
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
  use graph_carrier    , only : member_set, counted_set, subset_set
  use graph_relation   , only : relation
  use graph_structure  , only : relational_graph
  use graph_grammar    , only : graph, graph_field, graph_operation
  use class_graph_field, only : field
  use arithmetic_constitution_fixture, only : generated_residual

  implicit none

  private
  public :: constituted_residual

  type, extends(graph_operation) :: constituted_residual
     class(relation), allocatable :: flow      ! the GRAPH-OWNED copy
     class(relation), allocatable :: located
     type(counted_set)            :: xs, os, ys
     type(subset_set)             :: known, unknown
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
  ! the graph-owned citizen, found by identity and refused if G
  ! does not own it.
  !===================================================================!

  function create_adapter(g, selector, located, xs, os, ys, &
       &                  known, known_values, unknown) result(this)

    class(relational_graph), target, intent(in) :: g
    class(relation)  , intent(in) :: selector, located
    type(counted_set), intent(in) :: xs, os, ys
    type(subset_set) , intent(in) :: known, unknown
    real(dp)         , intent(in) :: known_values(:)
    type(constituted_residual)    :: this

    class(relation), pointer :: rp
    integer :: kk
    logical :: found

    found = .false.
    do kk = 1, g % num_relations()
       rp => g % relation_at(kk)
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
    this % xs = xs
    this % os = os
    this % ys = ys
    this % known   = known
    this % unknown = unknown
    this % known_values = known_values

  end function create_adapter

  pure function cr_name(this) result(name)
    class(constituted_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted residual'
  end function cr_name

  subroutine cr_domain(this, input_graph, domain)
    class(constituted_residual), intent(in) :: this
    class(graph), intent(in) :: input_graph
    class(member_set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % ys)
  end subroutine cr_domain

  subroutine cr_apply(this, input_graph, input_data, output)

    class(constituted_residual), intent(in)        :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)                    :: out
    class(member_set), allocatable :: dom
    real(dp), allocatable          :: ustate(:), r(:)

    associate (u1 => input_graph); end associate

    call input_data(1) % domain(dom)
    if (.not. dom % same_as(this % unknown)) then
       error stop 'statement: the state must live on the unknown domain'
    end if
    call input_data(1) % get_real_vector(ustate)

    call generated_residual(this % flow, this % located, &
         & this % xs, this % os, this % ys, &
         & this % known, this % known_values, this % unknown, ustate, r)

    out = field('residual', this % ys)
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
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_grammar    , only : graph_field
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_structure  , only : relational_graph, held_set, held_relation
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use constituted_residual_fixture, only : constituted_residual

  implicit none

  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type(counted_set)     :: x, o, p, y
  type(subset_set)      :: k, u, p_out
  type(stored_relation), allocatable :: flow
  type(stored_relation) :: located, r_out3, r_in3, produces, consumes
  class(relation), allocatable       :: d
  type(relational_graph), target     :: g
  type(stored_graph)    :: host
  type(constituted_residual) :: residual_op
  type(gmres)           :: solver
  type(field)           :: rhsf
  class(member_set), allocatable  :: dom
  class(graph_field), allocatable :: sol, rf
  real(dp), allocatable :: gv(:), solval(:), rv(:)
  real(dp)              :: answer
  integer               :: table(3, 6)
  integer               :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 9 . the statement"
  write(*,'(1x,a)') "============================================="

  ! -- the structure, and its graph
  x = counted_set('value-slots'  , 5)
  o = counted_set('operations'   , 2)
  p = counted_set('ports'        , 3)
  y = counted_set('residual-rows', 2)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  allocate(flow)
  flow = stored_relation('flow', [o, x, p], table)

  p_out    = subset_set('output-port', p, [PORT_OUT])
  r_out3   = restrict_slot(flow, 3, p_out)
  r_in3    = restrict_slot(flow, 3, &
       &       subset_set('input-ports', p, [PORT_IN1, PORT_IN2]))
  produces = project_slots(r_out3, [1, 2])
  consumes = project_slots(r_in3 , [2, 1])
  d        = compose_binary(produces, consumes)

  g = relational_graph('calculator', &
       & [held_set(x), held_set(o), held_set(p)], &
       & [held_relation(flow), held_relation(d)])

  ! -- the discretization, distinct from the graph on purpose
  located = stored_relation('located', [y, x], &
       & reshape([ROW_C, SLOT_C,  ROW_E, SLOT_E], [2, 2]))

  ! -- the inputs, and the unknowns in a NEW enumeration {c, e}
  k = subset_set('known'  , x, [SLOT_D, SLOT_A, SLOT_B])
  u = subset_set('unknown', x, [SLOT_C, SLOT_E])

  ! -- the adapter keeps the GRAPH-OWNED flow; the selector dies.
  residual_op = constituted_residual(g, flow, located, x, o, y, &
       &                             k, [4.0_dp, 2.0_dp, 3.0_dp], u)
  deallocate(flow)

  ! -- the compatibility scenery: not part of the statement.
  host = stored_graph(7, tails=[1,2,3,4,5,6], heads=[2,3,4,5,6,7])

  ! -- the ordinary solver, through its own operation face
  call solver % attach(residual_op, host, u)
  solver % tolerance      = 1.0d-12
  solver % max_iterations = 50

  call solver % domain(host, dom)
  call report(dom % same_as(u), &
       & "the solver answers on U = { c, e }", nfail)

  call solver % constant(gv)
  rhsf = field('rhs', y)
  call rhsf % set_real_vector(-gv)

  call solver % apply(host, [rhsf], sol)

  call sol % domain(dom)
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
  answer = solval(u % local_index(SLOT_E))

  call report(abs(answer - 20.0_dp) < 1.0d-9, &
       & "evaluate (2 + 3) x 4: the tower answers 20", nfail)

  write(*,'(1x,a,i0)') "CALCULATOR_RESULT = ", nint(answer)

  call verdict(nfail, "level 9")

end program calculator_level_9
