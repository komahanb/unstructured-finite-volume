!=====================================================================!
! The residual of one block, as one operation.
!
! A block's rows are of three kinds and they are added here into a
! single statement, because a minimizer drives one statement to zero:
!
!      derived     the scheme's own rows, already assembled as a
!                  stencil; linear in the state, so that stencil is
!                  also their jacobian
!      governing   the physics, at each slice's primary degree - the
!                  one degree no derived row determines
!      carried     the instants a block reaches back over, whose
!                  components are known before it starts; their rows
!                  are the identity less what they hold, so the block
!                  is square and nonsingular
!
! Its two arguments are the state and the design, in that order,
! which is what a minimizer supplies when the design is handed to it
! as a held input.
!
!             THE JACOBIAN
!
! Both halves carry exact partials - the stencil by being linear, the
! physics by differentiating its own rule - so the tangent is exact
! and nothing is differenced. A variation arrives named for this
! statement's argument and is renamed for each half before it is
! passed on, since each half checks the variation against its own.
!
!             WHAT IS REFUSED
!
! A state that is not one component per degree per slice; a missing
! argument; a carried row outside the unknowns.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_block

  use iso_fortran_env      , only : dp => REAL64
  use operation_action     , only : operation, variation
  use view_directed        , only : directed_graph
  use view_directed_stored , only : stored_directed_graph
  use field_calculus       , only : field
  use field_stored         , only : stored_field
  use graph_fractal        , only : graph
  use operation_stencil    , only : stencil
  use physics_integrand    , only : nodal_integrand

  implicit none

  private
  public :: block_residual

  type, extends(operation) :: block_residual

     type(stencil)                      , private :: derived
     class(nodal_integrand), allocatable , private :: physics
     type(stored_directed_graph)        , private :: instants
     integer , allocatable              , private :: carried(:)
     real(dp), allocatable              , private :: held(:)
     integer                            , private :: degrees = 0
     integer                            , private :: slices  = 0
     integer                            , private :: primary = 0

   contains

     procedure :: name           => block_name
     procedure :: domain         => block_domain
     procedure :: apply          => block_apply
     procedure :: max_degree     => block_max_degree
     procedure :: partial_action => block_partial_action
     procedure :: num_unknowns

  end type block_residual

  interface block_residual
     module procedure create
  end interface block_residual

contains

  function create(derived, physics, slices, degrees, primary, carried, held) result(this)

    type(stencil)         , intent(in) :: derived
    class(nodal_integrand), intent(in) :: physics
    integer               , intent(in) :: slices, degrees, primary
    integer               , intent(in) :: carried(:)
    real(dp)              , intent(in) :: held(:)
    type(block_residual) :: this

    if (size(carried) /= size(held)) then
       error stop 'gti_block: one value per carried component'
    end if
    if (any(carried < 1) .or. any(carried > slices * degrees)) then
       error stop 'gti_block: every carried row names an unknown'
    end if

    this % derived = derived
    allocate(this % physics, source=physics)
    this % slices  = slices
    this % degrees = degrees
    this % primary = primary
    this % carried = carried
    this % held    = held

    this % instants = stored_directed_graph(slices, tails=[integer ::], heads=[integer ::])
    call this % declare_arguments(2)

  end function create

  pure integer function num_unknowns(this)

    class(block_residual), intent(in) :: this

    num_unknowns = this % slices * this % degrees

  end function num_unknowns

  pure function block_name(this) result(name)

    class(block_residual), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'block residual'

  end function block_name

  subroutine block_domain(this, input_graph, domain, num_entries)

    class(block_residual), intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()

  end subroutine block_domain

  !===================================================================!
  ! The stencil supplies one order and the physics many, so the
  ! statement supplies one.
  !===================================================================!

  pure integer function block_max_degree(this)

    class(block_residual), intent(in) :: this

    associate (u1 => this); end associate
    block_max_degree = 1

  end function block_max_degree

  !===================================================================!
  ! The governing value at each slice, placed at the row its primary
  ! degree holds.
  !===================================================================!

  pure subroutine placed(this, governing, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: governing(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: k

    do k = 1, this % slices
       r((k - 1) * this % degrees + this % primary + 1) = &
            & r((k - 1) * this % degrees + this % primary + 1) + governing(k)
    end do

  end subroutine placed

  !===================================================================!
  ! The rows of the instants a block reaches back over: what they
  ! hold, less what is proposed for them.
  !===================================================================!

  pure subroutine carry(this, x, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: x(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = x(this % carried(i)) - this % held(i)
    end do

  end subroutine carry

  subroutine require_state(this, input_data)

    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)

    real(dp), allocatable :: x(:)

    if (size(input_data) < 2) then
       error stop 'gti_block: the state and the design are given'
    end if

    call input_data(1) % real_vector(x)
    if (size(x) /= this % num_unknowns()) then
       error stop 'gti_block: the state holds one component per degree per slice'
    end if

  end subroutine require_state

  subroutine block_apply(this, input_graph, input_data, output)

    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out, state
    class(field), allocatable :: half
    real(dp), allocatable :: r(:), governing(:), x(:)

    if (.not. present(input_data)) then
       error stop 'gti_block: the state and the design are given'
    end if
    call require_state(this, input_data)
    call input_data(1) % real_vector(x)

    state = stored_field('state', input_graph % vertex_set(), size(x))
    call state % set_real_vector(x)

    call this % derived % apply(input_graph, [state], half)
    call half % real_vector(r)

    call this % physics % apply(this % instants, input_data, half)
    call half % real_vector(governing)

    call placed(this, governing, r)
    call carry(this, x, r)

    out = stored_field(this % name(), input_graph % vertex_set(), size(r))
    call out % set_real_vector(r)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine block_apply

  !===================================================================!
  ! The tangent: each half differentiated in its own argument, the
  ! variation renamed for it, and the carried rows contributing the
  ! direction itself.
  !===================================================================!

  subroutine block_partial_action(this, input_graph, input_data, variations, output)

    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field)         , intent(in)        :: input_data(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out, state
    class(field), allocatable :: half
    real(dp), allocatable :: r(:), governing(:), v(:), x(:)

    call this % require_owned(variations)

    if (size(variations) /= 1) then
       error stop 'gti_block: the requested order is within max_degree'
    end if
    call require_state(this, input_data)
    call variations(1) % direction(v)
    call input_data(1) % real_vector(x)

    state = stored_field('state', input_graph % vertex_set(), size(x))
    call state % set_real_vector(x)

    call tangent_halves(this, input_graph, input_data, variations, state, r, governing)

    call placed(this, governing, r)
    call carry_direction(this, v, r)

    out = stored_field(this % name(), input_graph % vertex_set(), size(r))
    call out % set_real_vector(r)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine block_partial_action

  !===================================================================!
  ! Each half differentiated in its own argument, the variation
  ! renamed for it, since each half checks a variation against the
  ! argument it declared itself.
  !===================================================================!

  subroutine tangent_halves(this, input_graph, input_data, variations, state, r, governing)

    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    type(stored_field)   , intent(in) :: state
    real(dp), allocatable, intent(out) :: r(:), governing(:)

    class(field), allocatable :: half

    call this % derived % partial_action(input_graph, [state], &
         & [variations(1) % with_argument(this % derived % argument(1))], half)
    call half % real_vector(r)

    call this % physics % partial_action(this % instants, input_data, &
         & [variations(1) % with_argument(this % physics % argument(1))], half)
    call half % real_vector(governing)

  end subroutine tangent_halves

  pure subroutine carry_direction(this, v, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: v(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = v(this % carried(i))
    end do

  end subroutine carry_direction

end module gti_block
