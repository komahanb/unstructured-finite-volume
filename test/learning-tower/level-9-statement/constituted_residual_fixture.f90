!=====================================================================!
! THE CONSTITUTED LEARNING RESIDUAL - the Level-9 adapter, test-
! local: Level-8 semantics wearing the legacy graph_operation face
! so the ordinary solver can drive them. Its entire numerical act
! is delegation to learning_constitution_fixture::generated_residual
! - no multiplication here, no subtraction, no slot name, no order
! literal, no 2*w - 6 anywhere. What it holds is the STATEMENT'S
! furniture: the GRAPH-OWNED flow (found by identity, refused if
! the graph does not own it - the external selector may die), the
! location relation L, the carriers, the role subdomains, the
! observed values and the DERIVED execution order handed down by
! the caller's own topological sort.
!
! The adapter answers on Y; the minimizer answers on Theta. The
! trainable state arrives as a field on Theta, is judged by the
! complete constituted model - laws INTO the computed slots, the
! residual read at L's home - and leaves as a field on Y.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module constituted_residual_fixture

  use iso_fortran_env  , only : dp => REAL64
  use graph_carrier    , only : member_set, counted_set, subset_set
  use graph_relation   , only : relation
  ! graph_grammar exports a type named graph too; the kernel keeps
  ! the name and the grammar's is renamed at the door.
  use fractal_graph        , only : graph
  use graph_relational_view, only : relational_binding, &
       & num_relations, relation_at
  use graph_grammar    , only : grammar_graph => graph, graph_field, graph_operation
  use class_graph_field, only : field
  use learning_constitution_fixture, only : generated_residual

  implicit none

  private
  public :: constituted_learning_residual

  type, extends(graph_operation) :: constituted_learning_residual
     class(relation), allocatable :: flow      ! the GRAPH-OWNED copy
     class(relation), allocatable :: located
     type(counted_set)            :: slots, rows
     type(subset_set)             :: observed, trainable, computed
     real(dp), allocatable        :: observed_values(:)
     integer , allocatable        :: order(:)  ! derived by the caller
   contains
     procedure :: name   => clr_name
     procedure :: domain => clr_domain
     procedure :: apply  => clr_apply
  end type constituted_learning_residual

  interface constituted_learning_residual
     module procedure create_adapter
  end interface constituted_learning_residual

contains

  !===================================================================!
  ! The selector only NAMES the flow; what the adapter keeps is the
  ! graph-owned citizen, found by identity and refused if GAMMA does
  ! not own it.
  !===================================================================!

  function create_adapter(g, b, selector, located, slots, rows, &
       & observed, observed_values, trainable, computed, order) &
       & result(this)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(relation)  , intent(in) :: selector, located
    type(counted_set), intent(in) :: slots, rows
    type(subset_set) , intent(in) :: observed, trainable, computed
    real(dp)         , intent(in) :: observed_values(:)
    integer          , intent(in) :: order(:)
    type(constituted_learning_residual) :: this

    class(relation), pointer :: rp
    integer :: kk
    logical :: found

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
    this % slots     = slots
    this % rows      = rows
    this % observed  = observed
    this % trainable = trainable
    this % computed  = computed
    this % observed_values = observed_values
    this % order           = order

  end function create_adapter

  pure function clr_name(this) result(name)
    class(constituted_learning_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted learning residual'
  end function clr_name

  subroutine clr_domain(this, input_graph, domain)
    class(constituted_learning_residual), intent(in) :: this
    class(grammar_graph), intent(in) :: input_graph
    class(member_set), allocatable, intent(out) :: domain
    associate (u1 => input_graph); end associate
    allocate(domain, source=this % rows)
  end subroutine clr_domain

  subroutine clr_apply(this, input_graph, input_data, output)

    class(constituted_learning_residual), intent(in) :: this
    class(grammar_graph), intent(in)                 :: input_graph
    class(graph_field), intent(in), optional         :: input_data(:)
    class(graph_field), allocatable, intent(inout)   :: output

    type(field)                    :: out
    class(member_set), allocatable :: dom
    real(dp), allocatable          :: tstate(:), r(:)

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'statement: the residual needs a state to judge'
    end if
    if (size(input_data) < 1) then
       error stop 'statement: the residual needs a state to judge'
    end if

    call input_data(1) % domain(dom)
    if (.not. dom % same_as(this % trainable)) then
       error stop 'statement: the state must live on the trainable domain'
    end if
    call input_data(1) % get_real_vector(tstate)

    allocate(r(this % rows % size()))
    call generated_residual(this % flow, this % located, &
         & this % slots, this % rows, &
         & this % observed, this % observed_values, &
         & this % trainable, tstate, &
         & this % computed, this % order, r)

    out = field('residual', this % rows)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine clr_apply

end module constituted_residual_fixture
