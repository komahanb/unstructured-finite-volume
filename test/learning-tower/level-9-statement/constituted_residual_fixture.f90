!=====================================================================!
! THE CONSTITUTED LEARNING RESIDUAL - the Level-9 adapter, test-
! local: Level-8 semantics wearing the legacy operation face
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
  use graph_fractal        , only : set_graph => graph
  use map_set        , only : set_map
  use map_set_representation, only : set_representation
  use relation_finitary   , only : relation
  ! The kernel's graph and the ordinary view's directed_graph are two
  ! types with two names now, so nothing is renamed at this door.
  use graph_fractal        , only : graph
  use view_relational, only : relational_binding, &
       & num_relations, relation_at
  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  use learning_constitution_fixture, only : generated_residual

  implicit none

  private
  public :: constituted_learning_residual

  type, extends(operation) :: constituted_learning_residual
     class(relation), allocatable :: flow      ! the GRAPH-OWNED copy
     class(relation), allocatable :: located
     ! Identities, counts, and the action's OWN coordinates. A
     ! representation carries no identity, so holding one keeps this
     ! type free of any caller's map.
     type(set_graph)            :: slots, rows
     type(set_graph)             :: observed, trainable, computed
     integer :: n_slots = 0, n_rows = 0
     integer :: n_observed = 0, n_trainable = 0, n_computed = 0
     class(set_representation), allocatable :: c_slots, c_rows
     class(set_representation), allocatable :: c_observed, c_trainable, c_computed
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

  function create_adapter(g, b, selector, located, slots, rows, sets, &
       & observed, observed_values, trainable, computed, order) &
       & result(this)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(relation)  , intent(in) :: selector, located
    type(set_graph), intent(in) :: slots, rows
    type(set_map)  , intent(in) :: sets
    type(set_graph) , intent(in) :: observed, trainable, computed
    real(dp)         , intent(in) :: observed_values(:)
    integer          , intent(in) :: order(:)
    type(constituted_learning_residual) :: this

    class(relation), pointer :: rp
    integer :: kk
    logical :: found
    integer         :: n_computed
    integer         :: n_observed
    integer         :: n_rows
    integer         :: n_slots
    integer         :: n_trainable

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

    this % n_slots     = sets % size_of(slots)
    this % n_rows      = sets % size_of(rows)
    this % n_observed  = sets % size_of(observed)
    this % n_trainable = sets % size_of(trainable)
    this % n_computed  = sets % size_of(computed)

    call sets % extent_of(slots,     this % c_slots)
    call sets % extent_of(rows,      this % c_rows)
    call sets % extent_of(observed,  this % c_observed)
    call sets % extent_of(trainable, this % c_trainable)
    call sets % extent_of(computed,  this % c_computed)

  end function create_adapter

  pure function clr_name(this) result(name)
    class(constituted_learning_residual), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'constituted learning residual'
  end function clr_name

  subroutine clr_domain(this, input_graph, domain, nentries)
    class(constituted_learning_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    integer         :: n_rows
    associate (u1 => input_graph); end associate
    domain   = this % rows
    nentries = this % n_rows
  end subroutine clr_domain

  subroutine clr_apply(this, input_graph, input_data, output)

    class(constituted_learning_residual), intent(in) :: this
    class(directed_graph), intent(in)                 :: input_graph
    class(field), intent(in), optional         :: input_data(:)
    class(field), allocatable, intent(inout)   :: output

    type(stored_field)                    :: out
    type(set_graph) :: dom
    real(dp), allocatable          :: tstate(:), r(:)

    !----------------------------------------------------------------!
    ! The action's own coordinates, rebound into a LOCAL map for the
    ! duration of the call. The representations belong to the action;
    ! the map is a temporary and never leaves this scope.
    !----------------------------------------------------------------!

    type(set_map) :: mine
    integer         :: n_rows

    associate (u1 => input_graph); end associate

    call mine % bind(this % slots,     this % c_slots)
    call mine % bind(this % rows,      this % c_rows)
    call mine % bind(this % observed,  this % c_observed)
    call mine % bind(this % trainable, this % c_trainable)
    call mine % bind(this % computed,  this % c_computed)

    if (.not. present(input_data)) then
       error stop 'statement: the residual needs a state to judge'
    end if
    if (size(input_data) < 1) then
       error stop 'statement: the residual needs a state to judge'
    end if

    dom = input_data(1) % domain()
    if (.not. dom % same_as(this % trainable)) then
       error stop 'statement: the state must live on the trainable domain'
    end if
    call input_data(1) % get_real_vector(tstate)

    allocate(r(this % n_rows))
    call generated_residual(this % flow, this % located, &
         & this % slots, mine, this % rows, &
         & this % observed, this % observed_values, &
         & this % trainable, tstate, &
         & this % computed, this % order, r)

    out = stored_field('residual', this % rows, this % n_rows)
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine clr_apply

end module constituted_residual_fixture
