!=====================================================================!
! The linearization: the tangent J v of a statement S at a frozen
! state q, behind the operation interface, so a minimizer sees an
! ordinary linear operation. The primal S is written once; its
! tangent is this derived operation, evaluated by one of two roads
! chosen from S's max_degree:
!
!      exact       J v = D S(q) [v]                    max_degree >= 1,
!                  one partial action in input slot 1
!      difference  J v ~ ( S(q + eps v) - S(q) ) / eps  otherwise,
!                  two residuals, about eight digits
!
! freeze moves the frozen state between a governor's steps; a base
! residual handed to freeze is used by the difference road and
! ignored by the exact one. The tangent linearizes S : U -> U in its
! first input slot only.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_linearization

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field

  implicit none

  private
  public :: linearization
  public :: tangent_of

  type, extends(operation) :: linearization

     class(operation), allocatable :: of

     real(dp), allocatable :: at(:)
     real(dp), allocatable :: base(:)
     real(dp) :: step = 1.0d-7

   contains

     procedure :: name   => linearization_name
     procedure :: domain => linearization_domain
     procedure :: apply  => linearization_apply
     procedure :: freeze => linearization_freeze
     procedure :: exact  => linearization_exact

  end type linearization

contains

  !===================================================================!
  ! The tangent of a statement; the state and base may arrive now or
  ! through freeze.
  !===================================================================!

  function tangent_of(of, at, base) result(this)

    class(operation), intent(in)   :: of
    real(dp), intent(in), optional :: at(:)
    real(dp), intent(in), optional :: base(:)
    type(linearization)            :: this

    allocate(this % of, source=of)
    if (present(at))   allocate(this % at  , source=at)
    if (present(base)) allocate(this % base, source=base)

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
  ! Move the frozen state. The base residual is stored when handed
  ! over, and forgotten otherwise so the difference road measures it
  ! fresh.
  !===================================================================!

  subroutine linearization_freeze(this, at, base)

    class(linearization), intent(inout) :: this
    real(dp), intent(in)           :: at(:)
    real(dp), intent(in), optional :: base(:)

    this % at = at

    if (present(base)) then
       this % base = base
    else
       if (allocated(this % base)) deallocate(this % base)
    end if

  end subroutine linearization_freeze

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
  ! J v at the frozen state, on the statement's own domain. Checks,
  ! each stopping the program: the domain must be nonempty; a state
  ! must have been frozen; the frozen state must hold a whole number
  ! of components per domain member; the direction must live on the
  ! statement's domain and match the frozen state's size; every
  ! result of the statement must live on that same domain, because
  ! a field of equal length from another domain would pass
  ! otherwise. Without input data the direction is zero.
  !===================================================================!

  subroutine linearization_apply(this, input_graph, input_data, output)

    class(linearization), intent(in)         :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)   :: state, direction, out
    class(field), allocatable :: pushed
    type(graph) :: on, given
    real(dp), allocatable :: v(:), y(:), base(:)
    integer :: n_on, num_components

    call this % of % domain(input_graph, on, n_on)

    if (n_on <= 0) then
       error stop 'linearization: the operation''s domain is empty'
    end if
    if (.not. allocated(this % at)) then
       error stop 'linearization: the tangent is taken at a frozen state'
    end if
    if (mod(size(this % at), n_on) /= 0) then
       error stop 'linearization: the frozen state must carry a whole number &
            &per member of the operation''s domain'
    end if

    num_components = max(size(this % at) / n_on, 1)

    if (present(input_data)) then
       given = input_data(1) % domain()
       if (.not. given % same_as(on)) then
          error stop 'linearization: the direction must live on the operation''s domain'
       end if
       call input_data(1) % real_vector(v)
       if (size(v) /= size(this % at)) then
          error stop 'linearization: the direction must match the frozen state''s width'
       end if
    else
       allocate(v(size(this % at)))
       v = 0.0_dp
    end if

    state = stored_field('state', on, n_on, num_components=num_components)
    call state % set_real_vector(this % at)

    if (this % exact()) then

       direction = stored_field('direction', on, n_on, num_components=num_components)
       call direction % set_real_vector(v)
       call this % of % partial_action(input_graph, [state], [1], [direction], pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

    else

       ! the base residual: from freeze when handed over, measured
       ! here once when not
       if (allocated(this % base)) then
          base = this % base
       else
          call this % of % apply(input_graph, [state], pushed)
          call require_domain(pushed, on)
          call pushed % real_vector(base)
       end if

       call state % set_real_vector(this % at + this % step * v)
       call this % of % apply(input_graph, [state], pushed)
       call require_domain(pushed, on)
       call pushed % real_vector(y)

       y = (y - base) / this % step

    end if

    out = stored_field('J v', on, n_on, num_components=num_components)
    call out % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine linearization_apply

  !===================================================================!
  ! A same-domain tangent subtracts or contracts results, so each
  ! must have come from the statement's domain; equal length is not
  ! that claim.
  !===================================================================!

  subroutine require_domain(result, expected)

    class(field), intent(in) :: result
    type(graph) , intent(in) :: expected

    type(graph) :: got

    got = result % domain()
    if (.not. got % same_as(expected)) then
       error stop 'linearization: the operation result lives on its stated domain'
    end if

  end subroutine require_domain

end module operation_linearization
