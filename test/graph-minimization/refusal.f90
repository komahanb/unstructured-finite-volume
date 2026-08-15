!=====================================================================!
! The minimization refusal, EXPECTED TO DIE: the current solver
! family is square - an action whose residual dimension differs
! from the unknown dimension is refused at attach.
!=====================================================================!
module lopsided_fixture
  use iso_fortran_env  , only : dp => REAL64
  use graph_set    , only : set, index_set
  use graph_grammar    , only : ordinary_graph, graph_field, graph_operation
  use class_graph_field, only : field
  implicit none
  private
  public :: lopsided
  type, extends(graph_operation) :: lopsided
     type(index_set) :: y
   contains
     procedure :: name => l_name
     procedure :: domain => l_domain
     procedure :: apply => l_apply
  end type lopsided
contains
  pure function l_name(this) result(name)
    class(lopsided), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'lopsided'
  end function l_name
  subroutine l_domain(this, input_graph, domain)
    class(lopsided), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    class(set), allocatable, intent(out) :: domain
    associate (u => input_graph); end associate
    allocate(domain, source=this % y)
  end subroutine l_domain
  subroutine l_apply(this, input_graph, input_data, output)
    class(lopsided), intent(in) :: this
    class(ordinary_graph), intent(in) :: input_graph
    class(graph_field), intent(in), optional :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output
    type(field) :: out
    associate (u => input_graph, u2 => input_data); end associate
    out = field('r', this % y)
    call out % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine l_apply
end module lopsided_fixture

program minimization_refusal
  use graph_set    , only : index_set, subset
  use class_graph      , only : stored_graph
  use class_graph_gmres, only : gmres
  use lopsided_fixture , only : lopsided
  implicit none
  type(stored_graph) :: host
  type(index_set)  :: x
  type(subset)   :: u
  type(lopsided)     :: action
  type(gmres)        :: solver
  host = stored_graph(4, tails=[1,2,3], heads=[2,3,4])
  x = index_set('slots', 5)
  u = subset('unknowns', x, [5, 1])       ! two unknowns
  action % y = index_set('rows', 3)          ! three residuals
  call solver % attach(action, host, u)
  write(*,'(1x,a)') "REACHED PAST THE REFUSAL"
end program minimization_refusal
