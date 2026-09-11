!=====================================================================!
! The minimization refusal, EXPECTED TO TERMINATE: the current solver
! family is square - an action whose residual dimension differs
! from the unknown dimension is refused when stated.
!=====================================================================!
module rectangular_statement_fixture
  use iso_fortran_env  , only : dp => REAL64
  ! An action names a domain and counts it: the identity and the
  ! count are the whole of it.
  use graph_fractal    , only : graph
  use operation_action, only : operation, binding
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  implicit none
  private
  public :: rectangular_statement
  type, extends(operation) :: rectangular_statement
     type(graph) :: y
     integer         :: n_y = 0
   contains
     procedure :: name => rectangular_name
     procedure :: domain => rectangular_domain
     procedure :: apply => rectangular_apply
  end type rectangular_statement
  interface rectangular_statement
     module procedure create_rectangular_statement
  end interface rectangular_statement
contains
  ! The constructor declares the one argument, the state.
  function create_rectangular_statement() result(this)
    type(rectangular_statement) :: this
    call this % declare_arguments(1)
  end function create_rectangular_statement
  pure function rectangular_name(this) result(name)
    class(rectangular_statement), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'rectangular statement'
  end function rectangular_name
  subroutine rectangular_domain(this, input_graph, domain, num_entries)
    class(rectangular_statement), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    associate (u => input_graph); end associate
    domain   = this % y
    num_entries = this % n_y
  end subroutine rectangular_domain
  subroutine rectangular_apply(this, input_graph, inputs, output)
    class(rectangular_statement), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(binding), intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: out
    associate (u => input_graph, u2 => present(inputs)); end associate
    out = stored_field('r', this % y, this % n_y)
    call out % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine rectangular_apply
end module rectangular_statement_fixture

program minimization_refusal
  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_inclusion     , only : inclusion_map
  use view_directed_stored      , only : stored_directed_graph
  use operation_gmres, only : gmres
  use rectangular_statement_fixture , only : rectangular_statement
  implicit none
  type(stored_directed_graph)  :: host
  type(graph)     :: x, u
  type(rectangular_statement)      :: action
  type(gmres)         :: solver
  type(set_map)       :: sets
  type(inclusion_map) :: inclusions
  host = stored_directed_graph(4, tails=[1,2,3], heads=[2,3,4])
  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call u % declare()
  call sets       % bind(u, listed_set_representation([5, 1])) ! two unknowns
  call inclusions % include_in(u, x)
  action = rectangular_statement()
  call action % y % declare()
  call sets % bind(action % y, counted_set_representation(3))  ! three residuals
  action % n_y = 3
  call solver % state(action, host, u, sets % num_members_of(u))
  write(*,'(1x,a)') "REACHED PAST THE REFUSAL"
end program minimization_refusal
