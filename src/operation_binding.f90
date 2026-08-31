!=====================================================================!
! ONE OPERATION CALL, BOUND TO THE DATA GRAPH.
!
! An operation declares arguments and their contracts. The driver
! follows the graph, obtains one field for every argument and names
! the data vertices at which the result is to be placed. This module
! joins those facts into one value:
!
!       bound_operation = (input bindings, output destinations)
!
! An input binding owns the field snapshot and checks the argument's
! contract. An output destination is not an argument and carries no
! value before the operation runs; it names where the driver places
! the result afterwards.
!
! The destination ordinal is meaningful only in the pairing held by
! the driver that created it. Giving data graphs explicit identity
! would let that identity travel here as well; until then this type
! prevents inputs and outputs from being confused but cannot compare
! destinations belonging to different pairings.
!=====================================================================!

module operation_binding

  use util_precision  , only : dp
  use operation_action, only : operation, argument, contract
  use field_calculus  , only : field

  implicit none

  private
  public :: binding, output_destination, bound_operation
  public :: bind_operation
  public :: bound_inputs, bound_available_inputs, bound_value
  public :: bound_integer_vector, bound_real_vector

  type :: binding

     type(argument)            , private :: to
     class(field), allocatable , private :: value

   contains

     procedure :: argument_is => binding_argument_is
     procedure :: argument    => binding_argument
     procedure :: place       => binding_place

  end type binding

  interface binding
     module procedure create_binding
  end interface binding

  type :: output_destination

     integer, private :: vertex = 0

   contains

     procedure :: data_vertex => destination_data_vertex

  end type output_destination

  type :: bound_operation

     type(binding), allocatable            , private :: inputs(:)
     type(output_destination), allocatable , private :: outputs(:)

   contains

     procedure :: num_inputs
     procedure :: input
     procedure :: num_outputs
     procedure :: output

  end type bound_operation

  interface bind_operation
     module procedure bind_existing
     module procedure bind_with_inputs
     module procedure bind_without_inputs
  end interface bind_operation

contains

  function create_binding(to, value) result(this)

    type(argument), intent(in) :: to
    class(field)  , intent(in) :: value
    type(binding) :: this
    type(contract) :: required

    if (.not. to % is_named()) then
       error stop 'operation_binding: a binding names an argument'
    end if
    required = to % contract()
    if (.not. required % accepts(value)) then
       error stop 'operation_binding: a bound field satisfies its argument contract'
    end if

    this % to = to
    allocate(this % value, source=value)

  end function create_binding

  pure logical function binding_argument_is(this, a)

    class(binding), intent(in) :: this
    type(argument), intent(in) :: a

    binding_argument_is = this % to % matches(a)

  end function binding_argument_is

  pure function binding_argument(this) result(a)

    class(binding), intent(in) :: this
    type(argument) :: a

    a = this % to

  end function binding_argument

  subroutine binding_place(this, value)

    class(binding), intent(in) :: this
    class(field), allocatable, intent(inout) :: value

    if (.not. allocated(this % value)) then
       error stop 'operation_binding: a binding carries a field'
    end if

    if (allocated(value)) deallocate(value)
    allocate(value, source=this % value)

  end subroutine binding_place

  subroutine bound_inputs(action, input_data, bound)

    class(operation), intent(in) :: action
    class(field)    , intent(in) :: input_data(:)
    type(binding), allocatable, intent(out) :: bound(:)

    integer :: k

    if (size(input_data) /= action % num_arguments()) then
       error stop 'operation_binding: every declared argument is bound exactly once'
    end if

    allocate(bound(size(input_data)))
    do k = 1, size(input_data)
       bound(k) = binding(action % argument(k), input_data(k))
    end do

  end subroutine bound_inputs

  subroutine bound_available_inputs(action, input_data, bound)

    class(operation), intent(in) :: action
    class(field)    , intent(in) :: input_data(:)
    type(binding), allocatable, intent(out) :: bound(:)

    integer :: k

    if (size(input_data) > action % num_arguments()) then
       error stop 'operation_binding: supplied inputs name declared arguments'
    end if

    allocate(bound(size(input_data)))
    do k = 1, size(input_data)
       bound(k) = binding(action % argument(k), input_data(k))
    end do

  end subroutine bound_available_inputs

  subroutine bound_value(bound, a, value)

    type(binding), intent(in) :: bound(:)
    type(argument), intent(in) :: a
    class(field), allocatable, intent(inout) :: value

    integer :: k, found

    if (.not. a % is_named()) then
       error stop 'operation_binding: a bound value is named by an argument'
    end if

    found = 0
    do k = 1, size(bound)
       if (bound(k) % argument_is(a)) then
          if (found /= 0) then
             error stop 'operation_binding: an argument is bound once'
          end if
          found = k
       end if
    end do

    if (found == 0) then
       error stop 'operation_binding: the argument is bound'
    end if

    call bound(found) % place(value)

  end subroutine bound_value

  subroutine bound_integer_vector(bound, a, values)

    type(binding), intent(in) :: bound(:)
    type(argument), intent(in) :: a
    integer, allocatable, intent(out) :: values(:)

    class(field), allocatable :: value

    call bound_value(bound, a, value)
    call value % integer_vector(values)

  end subroutine bound_integer_vector

  subroutine bound_real_vector(bound, a, values)

    type(binding), intent(in) :: bound(:)
    type(argument), intent(in) :: a
    real(dp), allocatable, intent(out) :: values(:)

    class(field), allocatable :: value

    call bound_value(bound, a, value)
    call value % real_vector(values)

  end subroutine bound_real_vector

  function bind_existing(action, inputs, output_vertices) result(this)

    class(operation), intent(in) :: action
    type(binding)   , intent(in) :: inputs(:)
    integer         , intent(in) :: output_vertices(:)
    type(bound_operation) :: this

    integer :: i, k, found

    if (size(inputs) /= action % num_arguments()) then
       error stop 'operation_binding: every declared argument is bound exactly once'
    end if

    do i = 1, size(inputs)
       if (.not. action % owns(inputs(i) % argument())) then
          error stop 'operation_binding: every input belongs to the bound operation'
       end if
    end do

    do k = 1, action % num_arguments()
       found = 0
       do i = 1, size(inputs)
          if (inputs(i) % argument_is(action % argument(k))) found = found + 1
       end do
       if (found /= 1) then
          error stop 'operation_binding: every declared argument is bound exactly once'
       end if
    end do

    this % inputs = inputs
    call bind_outputs(output_vertices, this % outputs)

  end function bind_existing

  function bind_with_inputs(action, input_data, output_vertices) result(this)

    class(operation), intent(in) :: action
    class(field)    , intent(in) :: input_data(:)
    integer         , intent(in) :: output_vertices(:)
    type(bound_operation) :: this

    call bound_inputs(action, input_data, this % inputs)
    call bind_outputs(output_vertices, this % outputs)

  end function bind_with_inputs

  function bind_without_inputs(action, output_vertices) result(this)

    class(operation), intent(in) :: action
    integer         , intent(in) :: output_vertices(:)
    type(bound_operation) :: this

    if (action % num_arguments() /= 0) then
       error stop 'operation_binding: every declared argument is bound exactly once'
    end if

    allocate(this % inputs(0))
    call bind_outputs(output_vertices, this % outputs)

  end function bind_without_inputs

  subroutine bind_outputs(vertices, outputs)

    integer, intent(in) :: vertices(:)
    type(output_destination), allocatable, intent(out) :: outputs(:)

    integer :: i

    if (size(vertices) < 1 .or. any(vertices < 1)) then
       error stop 'operation_binding: an operation result has a data destination'
    end if

    do i = 1, size(vertices)
       if (count(vertices == vertices(i)) /= 1) then
          error stop 'operation_binding: each output destination is named once'
       end if
    end do

    allocate(outputs(size(vertices)))
    do i = 1, size(vertices)
       outputs(i) % vertex = vertices(i)
    end do

  end subroutine bind_outputs

  pure integer function destination_data_vertex(this) result(vertex)

    class(output_destination), intent(in) :: this

    if (this % vertex < 1) then
       error stop 'operation_binding: an output destination names a data vertex'
    end if
    vertex = this % vertex

  end function destination_data_vertex

  pure integer function num_inputs(this)

    class(bound_operation), intent(in) :: this

    num_inputs = 0
    if (allocated(this % inputs)) num_inputs = size(this % inputs)

  end function num_inputs

  pure function input(this, k) result(one)

    class(bound_operation), intent(in) :: this
    integer               , intent(in) :: k
    type(binding) :: one

    if (.not. allocated(this % inputs) .or. k < 1 .or. k > size(this % inputs)) then
       error stop 'operation_binding: an input is one the operation binds'
    end if
    one = this % inputs(k)

  end function input

  pure integer function num_outputs(this)

    class(bound_operation), intent(in) :: this

    num_outputs = 0
    if (allocated(this % outputs)) num_outputs = size(this % outputs)

  end function num_outputs

  pure function output(this, k) result(one)

    class(bound_operation), intent(in) :: this
    integer               , intent(in) :: k
    type(output_destination) :: one

    if (.not. allocated(this % outputs) .or. k < 1 .or. k > size(this % outputs)) then
       error stop 'operation_binding: an output is one the operation binds'
    end if
    one = this % outputs(k)

  end function output

end module operation_binding
