!=====================================================================!
! The linearization: the tangent of a statement S in one of its
! arguments, at a frozen input tuple, behind the operation interface,
! so a minimizer reads an ordinary linear operation. The primal S is
! written once; its tangent is this derived operation, evaluated by
! one of two modes chosen from S's max_degree:
!
!      exact       D_a S(x) [v]                          max_degree >= 1,
!                  one partial action, one variation (a, v)
!      difference  ( S(x + eps v e_a) - S(x) ) / eps     otherwise,
!                  two residuals, about eight digits
!
! The argument a is any argument S owns; it defaults to S's first.
! The frozen point is the whole input tuple [x_1, ..., x_m], moved by
! freeze. A base residual passed to freeze is used by the difference
! mode and ignored by the exact mode.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_linearization

  use util_precision  , only : dp, half_digits
  use operation_action, only : operation, argument, variation, contract
  use operation_action, only : binding, bound_value
  use operation_action, only : emit
  use view_directed, only : directed_graph
  use field_calculus, only : field, FIELD_REAL
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field

  implicit none

  private
  public :: linearization
  public :: tangent_of

  type, extends(operation) :: linearization

     class(operation), allocatable :: of

     type(argument) :: wrt

     type(stored_field), allocatable :: at(:)
     real(dp), allocatable :: base(:)
     real(dp) :: step = half_digits

   contains

     procedure :: name   => linearization_name
     procedure :: domain => linearization_domain
     procedure :: apply  => linearization_apply
     procedure :: exact  => linearization_exact

     procedure, private :: freeze_inputs
     generic :: freeze => freeze_inputs

  end type linearization

contains

  !===================================================================!
  ! The tangent of a statement in one argument (the first unless
  ! named); the frozen point and base arrive through freeze. An
  ! argument the statement does not own stops the program.
  !===================================================================!

  function tangent_of(of, wrt) result(this)

    class(operation), intent(in)             :: of
    type(argument), intent(in), optional     :: wrt
    type(linearization)                      :: this

    allocate(this % of, source=of)

    if (present(wrt)) then
       if (.not. of % owns(wrt)) then
          error stop 'linearization: the argument is one the statement owns'
       end if
       this % wrt = wrt
    else
       this % wrt = of % argument(1)
    end if

    ! the tangent reads one direction, shaped as the argument differentiated
    call this % declare_arguments(1, [this % wrt % contract()])

  end function tangent_of

  !===================================================================!
  ! The exact mode is available when the statement computes at least
  ! a first partial action.
  !===================================================================!

  pure function linearization_exact(this) result(exact)

    class(linearization), intent(in) :: this
    logical :: exact

    exact = this % of % max_degree() >= 1

  end function linearization_exact

  !===================================================================!
  ! Move the frozen point to the whole input tuple. An empty tuple
  ! stops the program. The base residual is stored when passed, and
  ! cleared otherwise so the difference mode recomputes it.
  !===================================================================!

  subroutine freeze_inputs(this, at_inputs, base)

    class(linearization), intent(inout)  :: this
    type(stored_field), intent(in)       :: at_inputs(:)
    real(dp), intent(in), optional       :: base(:)

    if (size(at_inputs) < 1) then
       error stop 'linearization: the frozen tuple contains the statement''s inputs'
    end if

    this % at = at_inputs

    if (present(base)) then
       this % base = base
    else
       if (allocated(this % base)) deallocate(this % base)
    end if

  end subroutine freeze_inputs

  pure function linearization_name(this) result(name)

    class(linearization), intent(in) :: this
    character(len=:), allocatable :: name

    if (this % exact()) then
       name = 'exact derivative of ' // this % of % name()
    else
       name = 'derivative of ' // this % of % name()
    end if

  end function linearization_name

  subroutine linearization_domain(this, input_graph, domain, num_entries)

    class(linearization), intent(in)  :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer    , intent(out) :: num_entries

    call this % of % domain(input_graph, domain, num_entries)

  end subroutine linearization_domain

  !===================================================================!
  ! The frozen tuple and the position of the differentiated
  ! argument in it. Checks, each stopping the program: a point must
  ! have been frozen; the tuple must include the differentiated
  ! argument.
  !===================================================================!

  subroutine frozen_tuple(this, tuple, position)

    class(linearization), intent(in) :: this
    type(stored_field), allocatable, intent(out) :: tuple(:)
    integer             , intent(out) :: position

    integer :: k

    if (.not. allocated(this % at)) then
       error stop 'linearization: the tangent is taken at a frozen state'
    end if
    tuple = this % at

    position = 0
    do k = 1, this % of % num_arguments()
       if (this % wrt % matches(this % of % argument(k))) position = k
    end do
    if (position < 1 .or. position > size(tuple)) then
       error stop 'linearization: the frozen tuple includes the differentiated argument'
    end if

  end subroutine frozen_tuple

  !===================================================================!
  ! D_a S(x)[v] at the frozen tuple x, on the statement's own domain.
  ! Checks, each stopping the program: the domain must be nonempty;
  ! a point must have been frozen; the direction must be defined on
  ! the differentiated argument's domain and match its width; every
  ! result of the statement must be defined on the statement's domain,
  ! because a field of equal length from another domain would pass
  ! otherwise. Without input data the direction is zero.
  !===================================================================!

  subroutine linearization_apply(this, input_graph, inputs, output)

    class(linearization), intent(in)         :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field), allocatable :: tuple(:)
    type(stored_field)   :: direction, out
    class(field), allocatable :: pushed, bound_direction
    type(graph) :: on, along
    real(dp), allocatable :: v(:), y(:), base(:), x(:)
    integer :: n_on, p, width

    call this % of % domain(input_graph, on, n_on)

    if (n_on <= 0) then
       error stop 'linearization: the operation''s domain is empty'
    end if

    call frozen_tuple(this, tuple, p)

    along = tuple(p) % domain()
    call tuple(p) % real_vector(x)
    width = size(x)

    if (present(inputs)) then
       call bound_value(inputs, this % argument(1), bound_direction)
       if (.not. bound_direction % defined_on(along)) then
          error stop 'linearization: the direction must be defined on the differentiated argument''s domain'
       end if
       call bound_direction % real_vector(v)
       if (size(v) /= width) then
          error stop 'linearization: the direction must match the frozen state''s width'
       end if
    else
       allocate(v(width))
       v = 0.0_dp
    end if

    direction = stored_field('direction', along, tuple(p) % num_entries(), &
         & num_components=tuple(p) % num_components())
    call direction % set_real_vector(v)

    if (this % exact()) then

       call this % of % partial_action(input_graph, this % of % bind(tuple), &
            & [variation(this % wrt, direction)], pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

    else

       ! the base residual: from freeze when passed, computed here
       ! once otherwise
       if (allocated(this % base)) then
          base = this % base
       else
          call this % of % apply(input_graph, this % of % bind(tuple), pushed)
          call require_domain(pushed, on)
          call pushed % real_vector(base)
       end if

       call tuple(p) % set_real_vector(x + this % step * v)
       call this % of % apply(input_graph, this % of % bind(tuple), pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

       y = (y - base) / this % step

    end if

    out = stored_field('J v', on, n_on, num_components=max(size(y) / n_on, 1))
    call out % set_real_vector(y)
    call emit(out, output)

  end subroutine linearization_apply

  !===================================================================!
  ! A same-domain tangent subtracts or contracts results, so each
  ! must have come from the statement's domain; equal length is not
  ! that claim.
  !===================================================================!

  subroutine require_domain(result, expected)

    class(field), intent(in) :: result
    type(graph) , intent(in) :: expected

    if (.not. result % defined_on(expected)) then
       error stop 'linearization: the operation result is defined on its stated domain'
    end if

  end subroutine require_domain

end module operation_linearization
