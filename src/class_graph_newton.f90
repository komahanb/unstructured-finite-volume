!=====================================================================!
! Newton's iteration on the tower.
!
! For a nonlinear statement action(q) = 0: linearize where you
! stand, solve the linear question, step, repeat,
!
!      J(q) dq = -action(q)          q <- q + dq
!
! The Jacobian is never formed. Its action on a direction is a
! difference of two residuals,
!
!      J v  ~  ( action(q + eps v) - action(q) ) / eps
!
! wrapped as an operation - so the linear solver inside sees an
! ordinary operation and asks no questions. The difference buys
! generality and pays in precision: the residual floor is the
! machine epsilon over the step times the residual scale, about
! eight digits. A statement that can linearize itself exactly may
! be attached in place of the difference when those digits matter. The inner solver is
! HELD, not inherited: newton governs a linear solver of its own
! choosing, one level down, and hands it one linear question per
! step. That is the governance discipline, in one component.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_newton

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph, graph_field, graph_operation
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use graph_minimization , only : linear_solver

  implicit none

  private
  public :: newton

  !===================================================================!
  ! The derivative of an operation along a direction, as an
  ! operation. Private machinery: newton builds one per step.
  !===================================================================!

  type, extends(graph_operation) :: directional_derivative

     class(graph_operation), allocatable :: of
     real(dp), allocatable :: at(:)
     real(dp), allocatable :: base(:)
     real(dp) :: step = 1.0d-7

   contains

     procedure :: name   => derivative_name
     procedure :: domain => derivative_domain
     procedure :: apply  => derivative_apply

  end type directional_derivative

  !===================================================================!
  ! Newton itself: an iteration budget, a tolerance, and the linear
  ! solver it governs.
  !===================================================================!

  type :: newton

     class(linear_solver), allocatable :: inner

     integer  :: max_iterations = 20
     real(dp) :: tolerance      = 1.0d-10

   contains

     procedure :: solve

  end type newton

  interface directional_derivative
     module procedure directional_derivative_at
  end interface directional_derivative

contains

  !===================================================================!
  ! Drive action(q) toward zero from the given q.
  !===================================================================!

  subroutine solve(this, action, on, q, achieved)

    class(newton), intent(inout)          :: this
    class(graph_operation), intent(in)    :: action
    class(graph)          , intent(in)    :: on
    real(dp), intent(inout)               :: q(:)
    real(dp), intent(out)                 :: achieved

    type(directional_derivative) :: jacobian
    real(dp), allocatable :: residual(:), g(:), dq(:)
    real(dp) :: linear_achieved
    integer :: it

    allocate(dq(size(q)))

    do it = 1, this % max_iterations

       ! Where we stand: the residual, through the inner solver's
       ! own plumbing.
       call this % inner % attach(action, on)
       call this % inner % constant(g)
       call this % inner % matvec(q, residual)
       residual = residual + g

       achieved = this % inner % norm(residual)
       if (achieved < this % tolerance) return

       ! The linear question at this point, and its answer.
       jacobian = directional_derivative(action, q, residual)

       call this % inner % attach(jacobian, on)
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)

       q = q + dq

    end do

    call this % inner % attach(action, on)
    call this % inner % constant(g)
    call this % inner % matvec(q, residual)
    achieved = this % inner % norm(residual + g)

  end subroutine solve

  !===================================================================!
  ! The derivative operation, built where newton stands.
  !===================================================================!

  type(directional_derivative) function directional_derivative_at(of, at, base) &
       & result(this)

    class(graph_operation), intent(in) :: of
    real(dp), intent(in) :: at(:)
    real(dp), intent(in) :: base(:)

    allocate(this % of, source=of)
    allocate(this % at  , source=at)
    allocate(this % base, source=base)

  end function directional_derivative_at

  pure function derivative_name(this) result(name)

    class(directional_derivative), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'derivative of ' // this % of % name()

  end function derivative_name

  subroutine derivative_domain(this, input_graph, domain)

    class(directional_derivative), intent(in) :: this
    class(graph), intent(in)                  :: input_graph
    class(graph), allocatable, intent(out)    :: domain

    call this % of % domain(input_graph, domain)

  end subroutine derivative_domain

  !===================================================================!
  ! J v as two residuals and a difference. The direction arrives as
  ! the input field; the answer leaves on the same cells.
  !===================================================================!

  subroutine derivative_apply(this, input_graph, input_data, output)

    class(directional_derivative), intent(in)      :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(support) :: cells
    type(field)   :: state
    class(graph_field), allocatable :: pushed
    real(dp), allocatable :: v(:), y(:)
    integer :: nv, i

    nv = input_graph % num_vertices()

    if (present(input_data)) then
       call input_data(1) % get_real_vector(v)
    else
       allocate(v(nv))
       v = 0.0_dp
    end if

    cells = support(GRAPH_SIDE_VERTEX, [(i, i = 1, nv)])
    state = field('state', cells)
    call state % set_real_vector(this % at + this % step * v)

    call this % of % apply(input_graph, [state], pushed)
    call pushed % get_real_vector(y)

    y = (y - this % base) / this % step

    state = field('J v', cells)
    call state % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=state)

  end subroutine derivative_apply

end module class_graph_newton
