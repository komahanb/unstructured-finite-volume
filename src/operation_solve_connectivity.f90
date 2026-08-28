!=====================================================================!
! Solver assembly connectivity.
!
! A mesh connectivity names which point indices make one element.
! This module does the same thing for solver assembly: it names which
! operation-input slots and derivative-degree graphs make one residual,
! one Jacobian, one forward tangent block, or one transposed adjoint
! block. It carries indices only. The state, residual values, operation,
! graph host, and linear solver arrive later in the driver.
!
! The point is to keep the repeated derivative-degree partition
! structure as data. Newton, Halley corrections, forward tangents, and
! reverse adjoints then reuse the same integer stencil instead of each
! spelling the same slot laws inline.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_solve_connectivity

  use operation_chain_rule, only : partial_connectivity, partial_connectivity_graph

  implicit none

  private
  public :: assembly_connectivity
  public :: newton_residual_connectivity
  public :: newton_jacobian_connectivity
  public :: halley_connectivity
  public :: tangent_residual_connectivity
  public :: tangent_jacobian_connectivity
  public :: adjoint_residual_connectivity
  public :: adjoint_jacobian_connectivity

  !===================================================================!
  ! One algorithm-level assembly stencil.
  !
  ! input_slot names the tuple/path positions read by the assembly.
  ! state_slot_value is the position occupied by the unknown/current
  ! state. differentiated_slot_value is the operation argument whose
  ! first partial supplies a Jacobian-like block. history_slot names
  ! the time-read slots used by tangent and adjoint propagation.
  !
  ! partial_graph(k) is the lower-level chain-rule connectivity for
  ! D^k R(Q, nu)[s_1, ..., s_k]. That graph says which derivative
  ! degrees open the multilinear slots; this type says where the
  ! whole solver assembly uses those graphs.
  !===================================================================!

  type :: assembly_connectivity

     private

     integer :: state_slot_value = 1
     integer :: differentiated_slot_value = 0
     logical :: transposed_value = .false.
     integer, allocatable :: input_slot(:)
     integer, allocatable :: history_slot(:)
     integer, allocatable :: parameter_slot(:)
     type(partial_connectivity_graph), allocatable :: partial_graph(:)

   contains

     procedure :: num_inputs      => connectivity_num_inputs
     procedure :: input           => connectivity_input
     procedure :: state           => connectivity_state
     procedure :: differentiated  => connectivity_differentiated
     procedure :: transposed      => connectivity_transposed
     procedure :: num_history     => connectivity_num_history
     procedure :: history         => connectivity_history
     procedure :: num_parameters  => connectivity_num_parameters
     procedure :: parameter       => connectivity_parameter
     procedure :: max_degree      => connectivity_max_degree
     procedure :: partial         => connectivity_partial

  end type assembly_connectivity

  interface newton_residual_connectivity
     module procedure create_newton_residual_connectivity
  end interface newton_residual_connectivity

  interface newton_jacobian_connectivity
     module procedure create_newton_jacobian_connectivity
  end interface newton_jacobian_connectivity

  interface halley_connectivity
     module procedure create_halley_connectivity
  end interface halley_connectivity

  interface tangent_residual_connectivity
     module procedure create_tangent_residual_connectivity
  end interface tangent_residual_connectivity

  interface tangent_jacobian_connectivity
     module procedure create_tangent_jacobian_connectivity
  end interface tangent_jacobian_connectivity

  interface adjoint_residual_connectivity
     module procedure create_adjoint_residual_connectivity
  end interface adjoint_residual_connectivity

  interface adjoint_jacobian_connectivity
     module procedure create_adjoint_jacobian_connectivity
  end interface adjoint_jacobian_connectivity

contains

  !===================================================================!
  ! Newton residual assembly:
  !
  !      D^0 R(U, held) = R(U, held)
  !
  ! There is no differentiated slot and no derivative-degree graph.
  ! The connectivity only says that the state is input slot 1 and the
  ! rest of the input tuple is read unchanged.
  !===================================================================!

  function create_newton_residual_connectivity(num_inputs) result(this)

    integer, intent(in) :: num_inputs
    type(assembly_connectivity) :: this

    call set_inputs(this, num_inputs)
    this % state_slot_value = 1

  end function create_newton_residual_connectivity

  !===================================================================!
  ! Newton Jacobian assembly:
  !
  !      D_U R(U, held)[dU]
  !
  ! It reads the same tuple as the residual, but now slot 1 is opened
  ! as the differentiated slot. The driver decides whether this is
  ! compiled to a stencil or kept as a frozen linearization operation.
  !===================================================================!

  function create_newton_jacobian_connectivity(num_inputs) result(this)

    integer, intent(in) :: num_inputs
    type(assembly_connectivity) :: this

    call set_inputs(this, num_inputs)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    call set_partials(this, 1)

  end function create_newton_jacobian_connectivity

  !===================================================================!
  ! Higher Newton/Halley residual correction assembly:
  !
  !      D^s R(U)[delta_1, ..., delta_{s-1}, 0]
  !
  ! The unknown delta_s is deliberately absent from the right hand
  ! side; the missing term is the Jacobian block solved separately.
  !===================================================================!

  function create_halley_connectivity(max_degree) result(this)

    integer, intent(in) :: max_degree
    type(assembly_connectivity) :: this

    call set_inputs(this, 1)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    call set_partials(this, max_degree)

  end function create_halley_connectivity

  !===================================================================!
  ! Forward tangent residual assembly:
  !
  !      known contractions from the state, history, and parameter
  !      paths form the right hand side for the current tangent solve.
  !
  ! The reach-dependent history slots are data because BDF startup and
  ! BDF2 do not have the same stencil.
  !===================================================================!

  function create_tangent_residual_connectivity(max_degree, reach, &
       & num_parameters) result(this)

    integer, intent(in) :: max_degree
    integer, intent(in) :: reach
    integer, intent(in) :: num_parameters
    type(assembly_connectivity) :: this

    integer :: j

    call require_nonnegative(reach)
    call require_nonnegative(num_parameters)

    call set_inputs(this, 1 + reach + num_parameters)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    call set_partials(this, max_degree)

    allocate(this % history_slot(reach))
    do j = 1, reach
       this % history_slot(j) = j
    end do

    allocate(this % parameter_slot(num_parameters))
    do j = 1, num_parameters
       this % parameter_slot(j) = 1 + reach + j
    end do

  end function create_tangent_residual_connectivity

  !===================================================================!
  ! Forward tangent Jacobian assembly:
  !
  !      D_state R_k[ V_k^(s) ]
  !
  ! This is the diagonal block reused across derivative orders at the
  ! configured instant.
  !===================================================================!

  function create_tangent_jacobian_connectivity(num_inputs) result(this)

    integer, intent(in) :: num_inputs
    type(assembly_connectivity) :: this

    call set_inputs(this, num_inputs)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    call set_partials(this, 1)

  end function create_tangent_jacobian_connectivity

  !===================================================================!
  ! Adjoint residual-block assembly:
  !
  !      (D_history R_k)^T lambda_k
  !
  ! These are the off-diagonal residual partials read in reverse by
  ! the adjoint system. The connectivity is marked transposed because
  ! every block is consumed through its dual action.
  !===================================================================!

  function create_adjoint_residual_connectivity(reach) result(this)

    integer, intent(in) :: reach
    type(assembly_connectivity) :: this

    integer :: j

    call require_nonnegative(reach)

    call set_inputs(this, 1 + reach)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    this % transposed_value = .true.

    allocate(this % history_slot(reach))
    do j = 1, reach
       this % history_slot(j) = j
    end do

  end function create_adjoint_residual_connectivity

  !===================================================================!
  ! Adjoint Jacobian assembly:
  !
  !      (D_state R_k)^T W_k
  !
  ! This is the diagonal state block of the adjoint solve, the
  ! transpose of the tangent Jacobian block.
  !===================================================================!

  function create_adjoint_jacobian_connectivity(num_inputs) result(this)

    integer, intent(in) :: num_inputs
    type(assembly_connectivity) :: this

    call set_inputs(this, num_inputs)
    this % state_slot_value = 1
    this % differentiated_slot_value = 1
    this % transposed_value = .true.
    call set_partials(this, 1)

  end function create_adjoint_jacobian_connectivity

  subroutine set_inputs(this, num_inputs)

    type(assembly_connectivity), intent(inout) :: this
    integer                    , intent(in)    :: num_inputs

    integer :: j

    if (num_inputs < 1) then
       error stop 'connectivity: an assembly reads at least one input slot'
    end if

    allocate(this % input_slot(num_inputs))
    do j = 1, num_inputs
       this % input_slot(j) = j
    end do

  end subroutine set_inputs

  subroutine set_partials(this, max_degree)

    type(assembly_connectivity), intent(inout) :: this
    integer                    , intent(in)    :: max_degree

    integer :: degree

    call require_nonnegative(max_degree)

    allocate(this % partial_graph(max_degree))
    do degree = 1, max_degree
       this % partial_graph(degree) = partial_connectivity(degree)
    end do

  end subroutine set_partials

  pure subroutine require_nonnegative(value)

    integer, intent(in) :: value

    if (value < 0) then
       error stop 'connectivity: a count is nonnegative'
    end if

  end subroutine require_nonnegative

  pure subroutine require_index(count, index)

    integer, intent(in) :: count
    integer, intent(in) :: index

    if (index < 1 .or. index > count) then
       error stop 'connectivity: the requested slot exists'
    end if

  end subroutine require_index

  pure integer function connectivity_num_inputs(this) result(num_inputs)

    class(assembly_connectivity), intent(in) :: this

    if (allocated(this % input_slot)) then
       num_inputs = size(this % input_slot)
    else
       num_inputs = 0
    end if

  end function connectivity_num_inputs

  pure integer function connectivity_input(this, slot) result(input)

    class(assembly_connectivity), intent(in) :: this
    integer                     , intent(in) :: slot

    call require_index(this % num_inputs(), slot)

    input = this % input_slot(slot)

  end function connectivity_input

  pure integer function connectivity_state(this) result(state)

    class(assembly_connectivity), intent(in) :: this

    state = this % state_slot_value

  end function connectivity_state

  pure integer function connectivity_differentiated(this) result(slot)

    class(assembly_connectivity), intent(in) :: this

    if (this % differentiated_slot_value < 1) then
       error stop 'connectivity: the assembly has a differentiated slot'
    end if

    slot = this % differentiated_slot_value

  end function connectivity_differentiated

  pure logical function connectivity_transposed(this) result(transposed)

    class(assembly_connectivity), intent(in) :: this

    transposed = this % transposed_value

  end function connectivity_transposed

  pure integer function connectivity_num_history(this) result(num_history)

    class(assembly_connectivity), intent(in) :: this

    if (allocated(this % history_slot)) then
       num_history = size(this % history_slot)
    else
       num_history = 0
    end if

  end function connectivity_num_history

  pure integer function connectivity_history(this, slot) result(history)

    class(assembly_connectivity), intent(in) :: this
    integer                     , intent(in) :: slot

    call require_index(this % num_history(), slot)

    history = this % history_slot(slot)

  end function connectivity_history

  pure integer function connectivity_num_parameters(this) result(num_parameters)

    class(assembly_connectivity), intent(in) :: this

    if (allocated(this % parameter_slot)) then
       num_parameters = size(this % parameter_slot)
    else
       num_parameters = 0
    end if

  end function connectivity_num_parameters

  pure integer function connectivity_parameter(this, slot) result(parameter)

    class(assembly_connectivity), intent(in) :: this
    integer                     , intent(in) :: slot

    call require_index(this % num_parameters(), slot)

    parameter = this % parameter_slot(slot)

  end function connectivity_parameter

  pure integer function connectivity_max_degree(this) result(max_degree)

    class(assembly_connectivity), intent(in) :: this

    if (allocated(this % partial_graph)) then
       max_degree = size(this % partial_graph)
    else
       max_degree = 0
    end if

  end function connectivity_max_degree

  pure function connectivity_partial(this, degree) result(connectivity)

    class(assembly_connectivity), intent(in) :: this
    integer                     , intent(in) :: degree
    type(partial_connectivity_graph) :: connectivity

    call require_index(this % max_degree(), degree)

    connectivity = this % partial_graph(degree)

  end function connectivity_partial

end module operation_solve_connectivity
