!=====================================================================!
! The minimization refusal, EXPECTED TO DIE: the current solver
! family is square - an action whose residual dimension differs
! from the unknown dimension is refused at attach.
!=====================================================================!
module lopsided_fixture
  use iso_fortran_env  , only : dp => REAL64
  ! An action names a domain and counts it. It holds no map: the
  ! identity and the count are the whole of what it is entitled to.
  use graph_fractal    , only : graph
  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  implicit none
  private
  public :: lopsided
  type, extends(operation) :: lopsided
     type(graph) :: y
     integer         :: n_y = 0
   contains
     procedure :: name => l_name
     procedure :: domain => l_domain
     procedure :: apply => l_apply
  end type lopsided
  interface lopsided
     module procedure create_lopsided
  end interface lopsided
contains
  ! The constructor declares the one argument, the state.
  function create_lopsided() result(this)
    type(lopsided) :: this
    call this % declare_arguments(1)
  end function create_lopsided
  pure function l_name(this) result(name)
    class(lopsided), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'lopsided'
  end function l_name
  subroutine l_domain(this, input_graph, domain, num_entries)
    class(lopsided), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries
    associate (u => input_graph); end associate
    domain   = this % y
    num_entries = this % n_y
  end subroutine l_domain
  subroutine l_apply(this, input_graph, input_data, output)
    class(lopsided), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field), intent(in), optional :: input_data(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: out
    associate (u => input_graph, u2 => input_data); end associate
    out = stored_field('r', this % y, this % n_y)
    call out % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine l_apply
end module lopsided_fixture

program minimization_refusal
  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_inclusion     , only : inclusion_map
  use view_directed_stored      , only : stored_directed_graph
  use operation_gmres, only : gmres
  use lopsided_fixture , only : lopsided
  implicit none
  type(stored_directed_graph)  :: host
  type(graph)     :: x, u
  type(lopsided)      :: action
  type(gmres)         :: solver
  type(set_map)       :: sets
  type(inclusion_map) :: inclusions
  host = stored_directed_graph(4, tails=[1,2,3], heads=[2,3,4])
  call x % declare()
  call sets % bind(x, counted_set_representation(5))
  call u % declare()
  call sets       % bind(u, listed_set_representation([5, 1])) ! two unknowns
  call inclusions % include_in(u, x)
  action = lopsided()
  call action % y % declare()
  call sets % bind(action % y, counted_set_representation(3))  ! three residuals
  action % n_y = 3
  call solver % attach(action, host, u, sets % num_members_of(u))
  write(*,'(1x,a)') "REACHED PAST THE REFUSAL"
end program minimization_refusal
