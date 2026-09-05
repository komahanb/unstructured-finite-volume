!=====================================================================!
! The step rules a driver evaluates at each instant. A march is a
! digraph of instants: one rule per step, one datum per instant, the
! arcs stating which instant each step reads. Every rule here reads
! its history from the bindings the driver passes and returns the
! next state; nothing here records an instant's index.
!
!      explicit_step      q_n = q_(n-1) - h S(q_(n-1))
!      implicit_step      sum_j a_j q_(n-j) + h S(q_n) = 0, solved by
!                         newton; the coefficients a_j are the slopes
!                         at zero of the Lagrange basis through the
!                         uniform nodes 0, -1, ..., -k, with k the
!                         number of history states bound - backward
!                         euler at the first step, bdf k after
!      adjoint_step       lambda_(n-1) = lambda_n - h J^T lambda_n, the
!                         transpose derived from the law's tangent at
!                         the recorded state and compiled as a stencil
!
! Every law returns MINUS its velocity, so the state moves against
! the statement.
!=====================================================================!

module march_rules_fixture

  use iso_fortran_env       , only : dp => REAL64
  use operation_action      , only : operation, binding, bound_real_vector, emit_real
  use view_directed         , only : directed_graph
  use field_calculus        , only : field
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
  use operation_newton      , only : newton
  use operation_gmres       , only : gmres
  use operation_linearization, only : linearization, tangent_of
  use operation_stencil     , only : stencil
  use operation_family      , only : slope_at_zero
  use util_derivative_terms , only : derivative_terms, value

  implicit none

  private
  public :: explicit_step, implicit_step, adjoint_step

  type, extends(operation) :: explicit_step
     class(operation), allocatable :: law
     real(dp) :: h = 1.0_dp
   contains
     procedure :: apply => explicit_apply
  end type explicit_step

  interface explicit_step
     module procedure create_explicit
  end interface explicit_step

  !===================================================================!
  ! The residual of one implicit step in the unknown q_n, argument 1,
  ! with the history q_(n-1), ..., q_(n-k) as arguments 2 .. k+1.
  !===================================================================!

  type, extends(operation) :: implicit_residual
     class(operation), allocatable :: law
     real(dp) :: h = 1.0_dp
     real(dp), allocatable :: a(:)
   contains
     procedure :: apply => residual_apply
  end type implicit_residual

  type, extends(operation) :: implicit_step
     class(operation), allocatable :: law
     real(dp) :: h = 1.0_dp
     real(dp) :: tolerance        = 1.0d-12
     real(dp) :: linear_tolerance = 1.0d-14
   contains
     procedure :: apply => implicit_apply
  end type implicit_step

  interface implicit_step
     module procedure create_implicit
  end interface implicit_step

  type, extends(operation) :: adjoint_step
     class(operation), allocatable :: law
     real(dp) :: h = 1.0_dp
   contains
     procedure :: apply => adjoint_apply
  end type adjoint_step

  interface adjoint_step
     module procedure create_adjoint
  end interface adjoint_step

contains

  !===================================================================!
  ! The constructors declare the arguments: the previous state for
  ! the explicit step, `order` history states for the implicit step,
  ! the costate then the recorded state for the adjoint step.
  !===================================================================!

  function create_explicit(law, h) result(this)
    class(operation), intent(in) :: law
    real(dp)        , intent(in) :: h
    type(explicit_step) :: this
    allocate(this % law, source=law)
    this % h = h
    call this % declare_arguments(1, label='explicit step')
  end function create_explicit

  function create_implicit(law, h, order) result(this)
    class(operation), intent(in) :: law
    real(dp)        , intent(in) :: h
    integer         , intent(in) :: order
    type(implicit_step) :: this
    if (order < 1) error stop 'implicit_step: the order is positive'
    allocate(this % law, source=law)
    this % h = h
    call this % declare_arguments(order, label='implicit step')
  end function create_implicit

  function create_residual(law, h, a) result(this)
    class(operation), intent(in) :: law
    real(dp)        , intent(in) :: h
    real(dp)        , intent(in) :: a(0:)
    type(implicit_residual) :: this
    allocate(this % law, source=law)
    this % h = h
    this % a = a
    call this % declare_arguments(size(a), label='implicit residual')
  end function create_residual

  function create_adjoint(law, h) result(this)
    class(operation), intent(in) :: law
    real(dp)        , intent(in) :: h
    type(adjoint_step) :: this
    allocate(this % law, source=law)
    this % h = h
    call this % declare_arguments(2, label='adjoint step')
  end function create_adjoint

  !===================================================================!
  ! The state as a field on the vertex set: one entry per vertex,
  ! the component count read from the vector length.
  !===================================================================!

  function state_field(input_graph, label, q) result(state)
    class(directed_graph), intent(in) :: input_graph
    character(len=*)     , intent(in) :: label
    real(dp)             , intent(in) :: q(:)
    type(stored_field) :: state
    integer :: nv
    nv = input_graph % num_vertices()
    if (mod(size(q), nv) /= 0) error stop 'march_rules: the state has one block per vertex'
    state = stored_field(label, input_graph % vertex_set(), nv, num_components=size(q) / nv)
    call state % set_real_vector(q)
  end function state_field

  subroutine emit_state(input_graph, q, output)
    class(directed_graph), intent(in) :: input_graph
    real(dp)             , intent(in) :: q(:)
    class(field), allocatable, intent(inout) :: output
    integer :: nv
    nv = input_graph % num_vertices()
    call emit_real('state', input_graph % vertex_set(), nv, q, output, &
         & num_components=size(q) / nv)
  end subroutine emit_state

  !===================================================================!
  ! Minus the velocity, S(q), from the law.
  !===================================================================!

  subroutine velocity(law, input_graph, q, s)
    class(operation)     , intent(in) :: law
    class(directed_graph), intent(in) :: input_graph
    real(dp)             , intent(in) :: q(:)
    real(dp), allocatable, intent(out) :: s(:)
    class(field), allocatable :: out
    call law % apply(input_graph, law % bind([state_field(input_graph, 'state', q)]), out)
    call out % real_vector(s)
    if (size(s) /= size(q)) error stop 'march_rules: the law returns one value per state entry'
  end subroutine velocity

  !===================================================================!
  ! The explicit step. Without a bound previous state there is no
  ! step, and the output is left unallocated.
  !===================================================================!

  subroutine explicit_apply(this, input_graph, inputs, output)

    class(explicit_step) , intent(in)           :: this
    class(directed_graph), intent(in)           :: input_graph
    type(binding)        , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)    :: output

    real(dp), allocatable :: q(:), s(:)

    if (allocated(output)) deallocate(output)
    if (.not. present(inputs)) return

    call bound_real_vector(inputs, this % argument(1), q)
    call velocity(this % law, input_graph, q, s)
    call emit_state(input_graph, q - this % h * s, output)

  end subroutine explicit_apply

  !===================================================================!
  ! The residual sum_j a_j q_(n-j) + h S(q_n).
  !===================================================================!

  subroutine residual_apply(this, input_graph, inputs, output)

    class(implicit_residual), intent(in)        :: this
    class(directed_graph), intent(in)           :: input_graph
    type(binding)        , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)    :: output

    real(dp), allocatable :: q(:), s(:), r(:), history(:)
    integer :: j

    if (.not. present(inputs)) error stop 'implicit_residual: the unknown is bound'

    call bound_real_vector(inputs, this % argument(1), q)
    call velocity(this % law, input_graph, q, s)
    r = this % a(0) * q + this % h * s
    do j = 1, size(this % a) - 1
       call bound_real_vector(inputs, this % argument(j + 1), history)
       r = r + this % a(j) * history
    end do

    call emit_real('residual', input_graph % vertex_set(), input_graph % num_vertices(), &
         & r, output, num_components=size(r) / input_graph % num_vertices())

  end subroutine residual_apply

  !===================================================================!
  ! The implicit step: the coefficients from the uniform nodes over
  ! the bound history, the residual stated to newton with the history
  ! as its stored inputs, the previous state as the initial iterate.
  !===================================================================!

  subroutine implicit_apply(this, input_graph, inputs, output)

    class(implicit_step) , intent(in)           :: this
    class(directed_graph), intent(in)           :: input_graph
    type(binding)        , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)    :: output

    type(derivative_terms), allocatable :: nodes(:)
    type(stored_field), allocatable :: history(:)
    real(dp), allocatable :: a(:), q(:), zero(:), previous(:)
    type(newton) :: solver
    real(dp) :: achieved
    integer :: k, j, nv

    if (allocated(output)) deallocate(output)
    if (.not. present(inputs)) return

    k  = size(inputs)
    nv = input_graph % num_vertices()

    allocate(nodes(0:k), a(0:k))
    do j = 0, k
       nodes(j) = derivative_terms(-real(j, dp), 0)
    end do
    do j = 0, k
       a(j) = value(slope_at_zero(nodes, j))
    end do

    allocate(history(k))
    do j = 1, k
       call bound_real_vector(inputs, this % argument(j), previous)
       history(j) = state_field(input_graph, 'history', previous)
    end do
    call bound_real_vector(inputs, this % argument(1), q)

    allocate(solver % inner, source=gmres())
    solver % inner % tolerance = this % linear_tolerance
    solver % tolerance         = this % tolerance
    call solver % state(create_residual(this % law, this % h, a), input_graph, &
         & input_graph % vertex_set(), nv, num_components=size(q) / nv, &
         & stored_inputs=history)

    allocate(zero(size(q)))
    zero = 0.0_dp
    call solver % solve(zero, q, achieved)

    call emit_state(input_graph, q, output)

  end subroutine implicit_apply

  !===================================================================!
  ! The adjoint step. The tangent of the law is frozen at the recorded
  ! state, compiled onto a stencil over the state's width, transposed,
  ! and applied to the costate over the stencil's own pattern.
  !===================================================================!

  subroutine adjoint_apply(this, input_graph, inputs, output)

    class(adjoint_step)  , intent(in)           :: this
    class(directed_graph), intent(in)           :: input_graph
    type(binding)        , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout)    :: output

    type(linearization) :: tangent
    type(stencil) :: compiled, transposed
    type(stored_field) :: costate
    class(field), allocatable :: pushed
    real(dp), allocatable :: lambda(:), q(:), jt_lambda(:)

    if (allocated(output)) deallocate(output)
    if (.not. present(inputs)) return

    call bound_real_vector(inputs, this % argument(1), lambda)
    call bound_real_vector(inputs, this % argument(2), q)

    tangent = tangent_of(this % law)
    call tangent % freeze([state_field(input_graph, 'state', q)])
    compiled   = stencil(tangent, input_graph, size(q))
    transposed = compiled % transpose()

    costate = stored_field('costate', transposed % pattern % vertex_set(), size(lambda))
    call costate % set_real_vector(lambda)
    call transposed % apply(transposed % pattern, transposed % bind([costate]), pushed)
    call pushed % real_vector(jt_lambda)

    call emit_state(input_graph, lambda - this % h * jt_lambda, output)

  end subroutine adjoint_apply

end module march_rules_fixture
