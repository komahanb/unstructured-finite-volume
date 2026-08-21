!=====================================================================!
! The linearization: the tangent of a statement S in one of its
! arguments, at a frozen input tuple, behind the operation interface,
! so a minimizer sees an ordinary linear operation. The primal S is
! written once; its tangent is this derived operation, evaluated by
! one of two roads chosen from S's max_degree:
!
!      exact       D_a S(x) [v]                          max_degree >= 1,
!                  one partial action, one variation (a, v)
!      difference  ( S(x + eps v e_a) - S(x) ) / eps     otherwise,
!                  two residuals, about eight digits
!
! The argument a is any argument S owns; it defaults to S's first.
! The frozen point is the whole input tuple [x_1, ..., x_m]; freeze
! also accepts the first argument's values alone, for a statement of
! one argument, and builds the tuple on S's domain when applied. A
! base residual handed to freeze is used by the difference road and
! ignored by the exact one.
!
! dual_by_basis forms (D_a S)^T lambda under the Euclidean pairing
! on stored values, one application per basis vector of argument a:
! it serves any block, square or not, where the compiled transpose
! of a stencil serves only a square one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_linearization

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation, argument, variation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field

  implicit none

  private
  public :: linearization
  public :: tangent_of
  public :: dual_by_basis

  type, extends(operation) :: linearization

     class(operation), allocatable :: of

     type(argument) :: wrt

     type(stored_field), allocatable :: at(:)
     real(dp), allocatable :: at_values(:)
     real(dp), allocatable :: base(:)
     real(dp) :: step = 1.0d-7

   contains

     procedure :: name   => linearization_name
     procedure :: domain => linearization_domain
     procedure :: apply  => linearization_apply
     procedure :: exact  => linearization_exact

     procedure, private :: freeze_values
     procedure, private :: freeze_inputs
     generic :: freeze => freeze_values, freeze_inputs

  end type linearization

contains

  !===================================================================!
  ! The tangent of a statement in one argument (the first unless
  ! named); the frozen point and base may arrive now or through
  ! freeze. An argument the statement does not own stops the
  ! program.
  !===================================================================!

  function tangent_of(of, wrt, at_inputs, at, base) result(this)

    class(operation), intent(in)             :: of
    type(argument), intent(in), optional     :: wrt
    type(stored_field), intent(in), optional :: at_inputs(:)
    real(dp), intent(in), optional           :: at(:)
    real(dp), intent(in), optional           :: base(:)
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

    if (present(at_inputs)) call this % freeze(at_inputs, base)
    if (present(at))        call this % freeze(at, base)

    ! the tangent reads one input, the direction
    call this % declare_arguments(1)

  end function tangent_of

  !===================================================================!
  ! The exact road is open when the statement computes at least a
  ! first partial action.
  !===================================================================!

  pure function linearization_exact(this) result(exact)

    class(linearization), intent(in) :: this
    logical :: exact

    exact = this % of % max_degree() >= 1

  end function linearization_exact

  !===================================================================!
  ! Move the frozen point: the whole input tuple, or the first
  ! argument's values alone. The base residual is stored when handed
  ! over, and forgotten otherwise so the difference road measures it
  ! fresh.
  !===================================================================!

  subroutine freeze_inputs(this, at_inputs, base)

    class(linearization), intent(inout)  :: this
    type(stored_field), intent(in)       :: at_inputs(:)
    real(dp), intent(in), optional       :: base(:)

    if (size(at_inputs) < 1) then
       error stop 'linearization: the frozen tuple holds the statement''s inputs'
    end if

    if (allocated(this % at_values)) deallocate(this % at_values)
    this % at = at_inputs

    call keep_base(this, base)

  end subroutine freeze_inputs

  subroutine freeze_values(this, at, base)

    class(linearization), intent(inout) :: this
    real(dp), intent(in)           :: at(:)
    real(dp), intent(in), optional :: base(:)

    if (.not. this % wrt % matches(this % of % argument(1))) then
       error stop 'linearization: values alone freeze the first argument; &
            &freeze the input tuple for another'
    end if

    if (allocated(this % at)) deallocate(this % at)
    this % at_values = at

    call keep_base(this, base)

  end subroutine freeze_values

  subroutine keep_base(this, base)

    class(linearization), intent(inout) :: this
    real(dp), intent(in), optional :: base(:)

    if (present(base)) then
       this % base = base
    else
       if (allocated(this % base)) deallocate(this % base)
    end if

  end subroutine keep_base

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
  ! argument in it. Values frozen alone become the state on the
  ! statement's domain. Checks, each stopping the program: a point
  ! must have been frozen; values frozen alone must hold a whole
  ! number of components per domain member; the tuple must reach the
  ! differentiated argument.
  !===================================================================!

  subroutine frozen_tuple(this, on, n_on, tuple, position)

    class(linearization), intent(in) :: this
    type(graph)         , intent(in) :: on
    integer             , intent(in) :: n_on
    type(stored_field), allocatable, intent(out) :: tuple(:)
    integer             , intent(out) :: position

    integer :: k

    if (allocated(this % at)) then
       tuple = this % at
    else if (allocated(this % at_values)) then
       if (mod(size(this % at_values), n_on) /= 0) then
          error stop 'linearization: the frozen state must carry a whole number &
               &per member of the operation''s domain'
       end if
       allocate(tuple(1))
       tuple(1) = stored_field('state', on, n_on, &
            & num_components=max(size(this % at_values) / n_on, 1))
       call tuple(1) % set_real_vector(this % at_values)
    else
       error stop 'linearization: the tangent is taken at a frozen state'
    end if

    position = 0
    do k = 1, this % of % num_arguments()
       if (this % wrt % matches(this % of % argument(k))) position = k
    end do
    if (position < 1 .or. position > size(tuple)) then
       error stop 'linearization: the frozen tuple reaches the differentiated argument'
    end if

  end subroutine frozen_tuple

  !===================================================================!
  ! D_a S(x)[v] at the frozen tuple x, on the statement's own domain.
  ! Checks, each stopping the program: the domain must be nonempty;
  ! a point must have been frozen; the direction must live on the
  ! differentiated argument's domain and match its width; every
  ! result of the statement must live on the statement's domain,
  ! because a field of equal length from another domain would pass
  ! otherwise. Without input data the direction is zero.
  !===================================================================!

  subroutine linearization_apply(this, input_graph, input_data, output)

    class(linearization), intent(in)         :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field), allocatable :: tuple(:)
    type(stored_field)   :: direction, out
    class(field), allocatable :: pushed
    type(graph) :: on, along
    real(dp), allocatable :: v(:), y(:), base(:), x(:)
    integer :: n_on, p, width

    call this % of % domain(input_graph, on, n_on)

    if (n_on <= 0) then
       error stop 'linearization: the operation''s domain is empty'
    end if

    call frozen_tuple(this, on, n_on, tuple, p)

    along = tuple(p) % domain()
    call tuple(p) % real_vector(x)
    width = size(x)

    if (present(input_data)) then
       if (.not. input_data(1) % defined_on(along)) then
          error stop 'linearization: the direction must live on the differentiated argument''s domain'
       end if
       call input_data(1) % real_vector(v)
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

       call this % of % partial_action(input_graph, tuple, &
            & [variation(this % wrt, direction)], pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

    else

       ! the base residual: from freeze when handed over, measured
       ! here once when not
       if (allocated(this % base)) then
          base = this % base
       else
          call this % of % apply(input_graph, tuple, pushed)
          call require_domain(pushed, on)
          call pushed % real_vector(base)
       end if

       call tuple(p) % set_real_vector(x + this % step * v)
       call this % of % apply(input_graph, tuple, pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

       y = (y - base) / this % step

    end if

    out = stored_field('J v', on, n_on, num_components=max(size(y) / n_on, 1))
    call out % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine linearization_apply

  !===================================================================!
  ! (D_a S)^T lambda under the Euclidean pairing on stored values:
  ! entry j is < D_a S e_j, lambda >, one application of the tangent
  ! per basis vector of the differentiated argument. lambda must
  ! match the tangent's result width; a mismatch stops the program.
  ! The tangent must be frozen.
  !===================================================================!

  subroutine dual_by_basis(tangent, input_graph, lambda, g)

    type(linearization)  , intent(in)  :: tangent
    class(directed_graph), intent(in)  :: input_graph
    real(dp)             , intent(in)  :: lambda(:)
    real(dp), allocatable, intent(out) :: g(:)

    type(stored_field), allocatable :: tuple(:)
    type(stored_field) :: basis
    class(field), allocatable :: pushed
    type(graph) :: on, along
    real(dp), allocatable :: e(:), y(:)
    integer :: n_on, p, width, j

    call tangent % of % domain(input_graph, on, n_on)
    call frozen_tuple(tangent, on, n_on, tuple, p)

    along = tuple(p) % domain()
    call tuple(p) % real_vector(e)
    width = size(e)

    allocate(g(width))

    do j = 1, width
       e    = 0.0_dp
       e(j) = 1.0_dp
       basis = stored_field('basis', along, tuple(p) % num_entries(), &
            & num_components=tuple(p) % num_components())
       call basis % set_real_vector(e)
       call tangent % apply(input_graph, [basis], pushed)
       call pushed % real_vector(y)
       if (size(y) /= size(lambda)) then
          error stop 'linearization: the dual pairs the tangent''s result with lambda'
       end if
       g(j) = dot_product(y, lambda)
    end do

  end subroutine dual_by_basis

  !===================================================================!
  ! A same-domain tangent subtracts or contracts results, so each
  ! must have come from the statement's domain; equal length is not
  ! that claim.
  !===================================================================!

  subroutine require_domain(result, expected)

    class(field), intent(in) :: result
    type(graph) , intent(in) :: expected

    if (.not. result % defined_on(expected)) then
       error stop 'linearization: the operation result lives on its stated domain'
    end if

  end subroutine require_domain

end module operation_linearization
