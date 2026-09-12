module linear_rules
  use util_precision, only : dp
  use operation_action, only : operation, binding, is_bound, bound_real_vector, emit_real
  use operation_driver, only : vertex_rule
  use view_directed, only : directed_graph
  use field_calculus, only : field
  implicit none
  private
  public :: linear_rule, placed_rule

  type, extends(operation) :: linear_rule
     real(dp), allocatable :: coefficients(:)
     real(dp) :: constant = 0.0_dp
     logical :: defined = .true.
     integer, pointer :: evaluations => null()
   contains
     procedure :: apply
  end type linear_rule

  interface linear_rule
     module procedure create
  end interface linear_rule

  !===================================================================!
  ! A rule told its place by the driver: its value is the linear
  ! combination of its bound inputs with one coefficient per read,
  ! emitted with the vertex it was applied at and the sum of the
  ! vertices it read, so a caller can check what the driver stated.
  !===================================================================!

  type, extends(vertex_rule) :: placed_rule
     real(dp) :: coefficient = 1.0_dp
     real(dp) :: constant = 0.0_dp
   contains
     procedure :: apply => placed_apply
  end type placed_rule

contains

  function create(coefficients, constant, evaluations, defined) result(this)
    real(dp), intent(in) :: coefficients(:), constant
    integer, target, intent(inout), optional :: evaluations
    logical, intent(in), optional :: defined
    type(linear_rule) :: this
    this % coefficients = coefficients
    this % constant = constant
    if (present(evaluations)) this % evaluations => evaluations
    if (present(defined)) this % defined = defined
    call this % declare_arguments(size(coefficients), label='linear combination')
  end function create

  subroutine apply(this, input_graph, inputs, output)
    class(linear_rule), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(binding), intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output
    real(dp), allocatable :: values(:)
    real(dp) :: total
    integer :: argument

    if (associated(this % evaluations)) this % evaluations = this % evaluations + 1
    if (.not. this % defined) return
    total = this % constant
    if (present(inputs)) then
       do argument = 1, size(this % coefficients)
          if (.not. is_bound(inputs, this % argument(argument))) cycle
          call bound_real_vector(inputs, this % argument(argument), values)
          total = total + this % coefficients(argument) * sum(values)
       end do
    end if
    call emit_real('linear combination', input_graph % vertex_set(), 1, [total], output)
  end subroutine apply

  subroutine placed_apply(this, input_graph, inputs, output)
    class(placed_rule), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(binding), intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output
    real(dp), allocatable :: values(:)
    real(dp) :: total
    integer :: argument

    if (.not. allocated(this % reads)) error stop 'placed_rule: the driver states the vertices read'
    if (this % num_arguments() /= size(this % reads)) error stop 'placed_rule: one argument per vertex read'
    total = this % constant
    if (present(inputs)) then
       do argument = 1, size(this % reads)
          if (.not. is_bound(inputs, this % argument(argument))) cycle
          call bound_real_vector(inputs, this % argument(argument), values)
          total = total + this % coefficient * values(1)
       end do
    end if
    call emit_real('placed combination', input_graph % vertex_set(), 1, &
         & [total, real(this % vertex, dp), real(sum(this % reads), dp)], output, num_components=3)
  end subroutine placed_apply

end module linear_rules
