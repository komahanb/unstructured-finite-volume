!=====================================================================!
! A rule stated at one instant of a time hierarchy, over the
! components held there and the design.
!
!      input graph    the instants
!      input 1        the state, one flat field over instant and
!                     degree: the component of degree d at instant k
!                     is held at (k-1)(N+1) + d + 1
!      input 2        the design, one value per instant, so a design
!                     that varies in time needs no other shape
!      output         one value per instant
!
! A concretion supplies the rule at one instant and nothing else. The
! instants are traversed, the state is read, and the partials in
! either argument are taken here once, for every concretion.
!
!             THE PARTIALS
!
! Each rule is evaluated over derivative_terms, so the value and
! every mixed partial in the directions asked for are carried
! together and are exact. A variation may name either argument, so
! the Q-partials, the X-partials and the mixed Q-X partials all come
! from one path and to any degree. A Newton block wants the
! Q-partials as a row of numbers, which is one call per degree with a
! direction that is one at that degree and zero elsewhere.
!
!             WHAT IS REFUSED
!
! A degree below one; a state whose extent is not the instants times
! N+1; a design that is not one value per instant; a missing
! argument; a variation naming neither argument.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_integrand

  use util_precision  , only : dp
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use field_calculus        , only : field
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
  use util_derivative_terms , only : derivative_terms, mixed_partial, &
       & max_subset_width, operator(*)

  implicit none

  private
  public :: nodal_integrand
  public :: zero_integrand

  type, abstract, extends(operation) :: nodal_integrand

     integer, private :: degree = 2

   contains

     procedure(integrand_rule_interface), deferred :: at_instant

     procedure :: equation_degree
     procedure :: declare_degree
     procedure :: domain         => integrand_domain
     procedure :: apply          => integrand_apply
     procedure :: max_degree     => integrand_max_degree
     procedure :: partial_action => integrand_partial_action

  end type nodal_integrand

  !-------------------------------------------------------------------!
  ! The integrand that is zero at every point: the physics of a
  ! linear statement, whose whole content is its stencil.
  !-------------------------------------------------------------------!

  type, extends(nodal_integrand) :: zero_integrand
   contains
     procedure :: name       => zero_name
     procedure :: at_instant => zero_rule
  end type zero_integrand

  interface zero_integrand
     module procedure create_zero
  end interface zero_integrand

  abstract interface

     pure function integrand_rule_interface(this, q, nu) result(r)
       import :: nodal_integrand, derivative_terms
       class(nodal_integrand), intent(in) :: this
       type(derivative_terms), intent(in) :: q(0:)
       type(derivative_terms), intent(in) :: nu
       type(derivative_terms) :: r
     end function integrand_rule_interface

  end interface

contains

  !===================================================================!
  ! The degree of the equation, declared once by a concretion. A
  ! degree below one stops the program: there is no highest
  ! derivative then.
  !===================================================================!

  subroutine declare_degree(this, degree)

    class(nodal_integrand), intent(inout) :: this
    integer               , intent(in)    :: degree

    if (degree < 1) then
       error stop 'physics_integrand: the degree of the equation is positive'
    end if

    this % degree = degree
    call this % declare_arguments(2)

  end subroutine declare_degree

  pure integer function equation_degree(this)

    class(nodal_integrand), intent(in) :: this

    equation_degree = this % degree

  end function equation_degree

  !===================================================================!
  ! THE MACHINERY.
  !
  ! One value per instant.
  !===================================================================!

  subroutine integrand_domain(this, input_graph, domain, num_entries)

    class(nodal_integrand), intent(in)  :: this
    class(directed_graph) , intent(in)  :: input_graph
    type(graph)           , intent(out) :: domain
    integer               , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()

  end subroutine integrand_domain

  pure integer function integrand_max_degree(this)

    class(nodal_integrand), intent(in) :: this

    associate (u1 => this); end associate
    integrand_max_degree = max_subset_width()

  end function integrand_max_degree

  !===================================================================!
  ! The state and the design as derivative terms, each direction
  ! seeded on the argument its variation names.
  !===================================================================!

  subroutine seeded(this, input_data, variations, q, nu)

    class(nodal_integrand), intent(in) :: this
    class(field)          , intent(in) :: input_data(:)
    type(variation)       , intent(in) :: variations(:)
    type(derivative_terms), allocatable, intent(out) :: q(:), nu(:)

    real(dp), allocatable :: state(:), design(:), v(:)
    integer :: n, i

    if (size(input_data) < 2) then
       error stop 'physics_vanderpol: the state and the design are given'
    end if

    call input_data(1) % real_vector(state)
    call input_data(2) % real_vector(design)

    n = size(variations)
    call constants(state, n, q)
    call constants(design, n, nu)

    do i = 1, n
       call variations(i) % direction(v)
       if (variations(i) % argument_is(this % argument(1))) then
          call seed(q, i, v)
       else if (variations(i) % argument_is(this % argument(2))) then
          call seed(nu, i, v)
       else
          error stop 'physics_vanderpol: a variation names the state or the design'
       end if
    end do

  end subroutine seeded

  !===================================================================!
  ! A real vector as terms that carry no derivative yet.
  !===================================================================!

  subroutine constants(x, n, terms)

    real(dp), intent(in) :: x(:)
    integer , intent(in) :: n
    type(derivative_terms), allocatable, intent(out) :: terms(:)

    integer :: j

    allocate(terms(size(x)))

    do j = 1, size(x)
       terms(j) = derivative_terms(x(j), n)
    end do

  end subroutine constants

  !===================================================================!
  ! One direction laid into every entry of a quantity.
  !===================================================================!

  subroutine seed(x, i, v)

    type(derivative_terms), intent(inout) :: x(:)
    integer               , intent(in)    :: i
    real(dp)              , intent(in)    :: v(:)

    integer :: j

    if (size(v) /= size(x)) then
       error stop 'physics_vanderpol: a direction has one entry per unknown it varies'
    end if

    do j = 1, size(x)
       call x(j) % set_direction(i, v(j))
    end do

  end subroutine seed

  !===================================================================!
  ! The rule at every instant, answering the coefficient of the full
  ! subset: the value with no directions, the mixed partial with n.
  !===================================================================!

  subroutine evaluated(this, input_graph, q, nu, output)

    class(nodal_integrand), intent(in) :: this
    class(directed_graph) , intent(in) :: input_graph
    type(derivative_terms), intent(in) :: q(:), nu(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out
    real(dp), allocatable :: values(:)
    integer :: k, nd, base

    nd = this % degree + 1

    if (size(q) /= input_graph % num_vertices() * nd) then
       error stop 'physics_vanderpol: the state holds one component per degree per instant'
    end if
    if (size(nu) /= input_graph % num_vertices()) then
       error stop 'physics_vanderpol: the design holds one value per instant'
    end if

    allocate(values(input_graph % num_vertices()))

    do k = 1, input_graph % num_vertices()
       base = (k - 1) * nd
       values(k) = mixed_partial(this % at_instant(q(base + 1:base + nd), nu(k)))
    end do

    out = stored_field(this % name(), input_graph % vertex_set(), &
         & input_graph % num_vertices())
    call out % set_real_vector(values)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine evaluated

  subroutine integrand_apply(this, input_graph, input_data, output)

    class(nodal_integrand), intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(variation), allocatable :: none(:)
    type(derivative_terms), allocatable :: q(:), nu(:)

    if (.not. present(input_data)) then
       error stop 'physics_vanderpol: the state and the design are given'
    end if

    allocate(none(0))
    call seeded(this, input_data, none, q, nu)
    call evaluated(this, input_graph, q, nu, output)

  end subroutine integrand_apply

  subroutine integrand_partial_action(this, input_graph, input_data, variations, output)

    class(nodal_integrand), intent(in)       :: this
    class(directed_graph) , intent(in)       :: input_graph
    class(field)          , intent(in)       :: input_data(:)
    type(variation)       , intent(in)       :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: q(:), nu(:)

    call this % require_owned(variations)

    if (size(variations) > this % max_degree()) then
       error stop 'physics_vanderpol: the requested order is within max_degree'
    end if

    call seeded(this, input_data, variations, q, nu)
    call evaluated(this, input_graph, q, nu, output)

  end subroutine integrand_partial_action

  function create_zero(degree) result(this)

    integer, intent(in) :: degree
    type(zero_integrand) :: this

    call this % declare_degree(degree)

  end function create_zero

  pure function zero_name(this) result(name)

    class(zero_integrand), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'zero'

  end function zero_name

  pure function zero_rule(this, q, nu) result(r)

    class(zero_integrand) , intent(in) :: this
    type(derivative_terms), intent(in) :: q(0:)
    type(derivative_terms), intent(in) :: nu
    type(derivative_terms) :: r

    associate (u1 => this, u2 => nu); end associate
    r = 0.0_dp * q(0)

  end function zero_rule

end module physics_integrand
