!=====================================================================!
! Test fixtures shared by test.f90 and refusal.f90: four concrete
! operations used to exercise the chain rule, the exact
! linearization, the step table, and the marcher's derivative
! routines.
!
!      quartic_form     Phi(q, xi) = q^4 + q^3 xi + q^2 xi^2
!                                  + q xi^3 + xi^4,  max_degree 4
!      power8_form      Phi(q, xi) = (q + xi)^8,     max_degree 8
!      equilibrium_law  S_i(q, xi) = q_i^2 - xi,     max_degree 8;
!                       every partial of order 3 or higher is zero
!      linear_law       S(q) = q; implements graph_operation only,
!                       so tangent_of must select the difference
!                       linearization for it
!
! Partial derivatives are computed rather than tabulated: the
! quartic differentiates each of its five monomials with falling
! factorials, and every order-k mixed partial of power8 equals
! 8!/(8-k)! (q + xi)^(8-k).
!
! xi is read from input slot 2 when the caller passes two input
! fields, and from the xi_default component otherwise. The chain
! rule tests pass both inputs; exact_linearization passes only the
! state, so the default keeps one operation usable by both.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module toy_differentiable_forms

  use iso_fortran_env     , only : dp => REAL64
  use operation_action, only : graph_operation
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_discretization      , only : differentiable_operation
  use graph_fractal       , only : set_graph => graph
  use class_graph_field   , only : field
  use class_graph_chain_rule, only : argument_path

  implicit none

  private
  public :: quartic_form, power8_form, equilibrium_law, linear_law
  public :: scalar_pair, fill_path

  !===================================================================!
  ! Phi = q^4 + q^3 xi + q^2 xi^2 + q xi^3 + xi^4, scalar q in slot
  ! one and scalar xi in slot two, exact to order 4.
  !===================================================================!

  type, extends(differentiable_operation) :: quartic_form
     real(dp) :: xi_default = 2.0_dp
   contains
     procedure :: name           => quartic_name
     procedure :: domain         => quartic_domain
     procedure :: apply          => quartic_apply
     procedure :: max_degree     => quartic_max_degree
     procedure :: partial_action => quartic_partial_action
  end type quartic_form

  !===================================================================!
  ! Phi = (q + xi)^8. With u = q + xi, every mixed partial of
  ! order k, in any mix of the two slots, equals 8!/(8-k)! u^(8-k).
  ! Used by the degree-8 chain-rule test, whose expected values
  ! are recomputed independently by Taylor convolution.
  !===================================================================!

  type, extends(differentiable_operation) :: power8_form
     real(dp) :: xi_default = 2.0_dp
   contains
     procedure :: name           => power8_name
     procedure :: domain         => power8_domain
     procedure :: apply          => power8_apply
     procedure :: max_degree     => power8_max_degree
     procedure :: partial_action => power8_partial_action
  end type power8_form

  !===================================================================!
  ! S_i(q, xi) = q_i^2 - xi, elementwise: the operation the
  ! marcher tests integrate. At xi = 1 the fixed point is q = 1,
  ! so a march started there produces a constant trajectory while
  ! its parameter sensitivities are nonzero, which keeps the
  ! expected values in closed form. max_degree is 8 because the
  ! tests request compositions up to degree 8; every partial of
  ! order 3 or higher is returned as exactly zero.
  !===================================================================!

  type, extends(differentiable_operation) :: equilibrium_law
     real(dp) :: xi_default = 1.0_dp
   contains
     procedure :: name           => equilibrium_name
     procedure :: domain         => equilibrium_domain
     procedure :: apply          => equilibrium_apply
     procedure :: max_degree     => equilibrium_max_degree
     procedure :: partial_action => equilibrium_partial_action
  end type equilibrium_law

  !===================================================================!
  ! S(q) = q, implementing graph_operation only. Used to check
  ! that tangent_of falls back to the difference linearization for
  ! a non-differentiable operation, and as the statement of the
  ! nonuniform implicit march test.
  !===================================================================!

  type, extends(graph_operation) :: linear_law
   contains
     procedure :: name   => linear_name
     procedure :: domain => linear_domain
     procedure :: apply  => linear_apply
  end type linear_law

contains

  !===================================================================!
  ! Domain used by every toy: the input graph's vertex set.
  !===================================================================!

  subroutine vertex_domain(input_graph, domain, nentries)

    class(directed_graph), intent(in)  :: input_graph
    type(set_graph)      , intent(out) :: domain
    integer              , intent(out) :: nentries

    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()

  end subroutine vertex_domain

  !===================================================================!
  ! Read the state vector from input slot 1, and xi from input
  ! slot 2 when two inputs were passed or from xi_default when
  ! only the state was passed.
  !===================================================================!

  subroutine read_arguments(input_data, xi_default, q, xi)

    class(graph_field), intent(in)     :: input_data(:)
    real(dp)          , intent(in)     :: xi_default
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)          , intent(out)    :: xi

    real(dp), allocatable :: xv(:)

    call input_data(1) % get_real_vector(q)

    if (size(input_data) >= 2) then
       call input_data(2) % get_real_vector(xv)
       xi = xv(1)
    else
       xi = xi_default
    end if

  end subroutine read_arguments

  !===================================================================!
  ! Copy the given values into a field on the input graph's vertex
  ! set and place it in output.
  !===================================================================!

  subroutine pack_output(input_graph, values, output)

    class(directed_graph), intent(in) :: input_graph
    real(dp), intent(in) :: values(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph) :: cells
    type(field)     :: out
    integer         :: nentries

    cells    = input_graph % vertex_set()
    nentries = input_graph % num_vertices()

    out = field('toy value', cells, nentries, &
         & ncomp = size(values) / nentries)
    call out % set_real_vector(values)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine pack_output

  !===================================================================!
  ! The falling factorial n (n-1) ... (n-k+1): the coefficient of
  ! the k-th derivative of an n-th power. It is zero when k > n,
  ! which makes derivatives past an exponent vanish.
  !===================================================================!

  pure function falling(n, k) result(f)

    integer, intent(in) :: n, k
    real(dp) :: f

    integer :: j

    f = 1.0_dp
    do j = 0, k - 1
       f = f * real(n - j, dp)
    end do

  end function falling

  !===================================================================!
  ! Count how many entries of slots(:) name input 1 (a) and input
  ! 2 (b), and accumulate the product of the first value of every
  ! direction field.
  !===================================================================!

  subroutine count_seats(slots, directions, a, b, product)

    integer           , intent(in)  :: slots(:)
    class(graph_field), intent(in)  :: directions(:)
    integer           , intent(out) :: a, b
    real(dp)          , intent(out) :: product

    real(dp), allocatable :: v(:)
    integer :: j

    a = count(slots == 1)
    b = count(slots == 2)

    product = 1.0_dp
    do j = 1, size(slots)
       call directions(j) % get_real_vector(v)
       product = product * v(1)
    end do

  end subroutine count_seats

  !===================================================================!
  ! Build the two scalar input fields (q, xi) on the given vertex
  ! set.
  !===================================================================!

  subroutine scalar_pair(q, xi, cells, inputs)

    real(dp)       , intent(in)  :: q, xi
    type(set_graph), intent(in)  :: cells
    type(field)    , intent(out) :: inputs(2)

    inputs(1) = field('q', cells, 1, ncomp=1)
    call inputs(1) % set_real_vector([q])
    inputs(2) = field('xi', cells, 1, ncomp=1)
    call inputs(2) % set_real_vector([xi])

  end subroutine scalar_pair

  !===================================================================!
  ! Build an argument_path for the given input slot whose
  ! derivative(k) holds the k-th scalar, every entry marked
  ! occupied.
  !===================================================================!

  subroutine fill_path(path, slot, derivatives, cells)

    type(argument_path), intent(inout) :: path
    integer            , intent(in)    :: slot
    real(dp)           , intent(in)    :: derivatives(:)
    type(set_graph)    , intent(in)    :: cells

    integer :: k

    path % slot = slot
    if (allocated(path % derivative)) deallocate(path % derivative)
    allocate(path % derivative(size(derivatives)))

    do k = 1, size(derivatives)
       path % derivative(k) % occupied  = .true.
       path % derivative(k) % direction = field('path', cells, 1, ncomp=1)
       call path % derivative(k) % direction % &
            & set_real_vector([derivatives(k)])
    end do

  end subroutine fill_path

  !===================================================================!
  ! The quartic.
  !===================================================================!

  pure function quartic_name(this) result(name)
    class(quartic_form), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'toy quartic'
  end function quartic_name

  subroutine quartic_domain(this, input_graph, domain, nentries)
    class(quartic_form), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    call vertex_domain(input_graph, domain, nentries)
  end subroutine quartic_domain

  pure function quartic_max_degree(this) result(degree)
    class(quartic_form), intent(in) :: this
    integer :: degree
    associate (u1 => this); end associate
    degree = 4
  end function quartic_max_degree

  subroutine quartic_apply(this, input_graph, input_data, output)

    class(quartic_form), intent(in)                 :: this
    class(directed_graph), intent(in)               :: input_graph
    class(graph_field), intent(in), optional        :: input_data(:)
    class(graph_field), allocatable, intent(inout)  :: output

    real(dp), allocatable :: q(:)
    real(dp) :: xi, phi

    phi = 0.0_dp
    if (present(input_data)) then
       call read_arguments(input_data, this % xi_default, q, xi)
       phi = quartic_mixed(0, 0, q(1), xi)
    end if

    call pack_output(input_graph, [phi], output)

  end subroutine quartic_apply

  !-------------------------------------------------------------------!
  ! The mixed partial d^a/dq^a d^b/dxi^b of Phi, computed per
  ! monomial q^i xi^j with j = 4 - i: each term contributes
  ! falling(i, a) falling(j, b) q^(i-a) xi^(j-b) and is skipped
  ! when a > i or b > j. At (a, b) = (0, 0) this returns Phi.
  !-------------------------------------------------------------------!

  pure function quartic_mixed(a, b, q, xi) result(m)

    integer , intent(in) :: a, b
    real(dp), intent(in) :: q, xi
    real(dp) :: m

    integer :: i, j

    m = 0.0_dp
    do i = 0, 4
       j = 4 - i
       if (a > i .or. b > j) cycle
       m = m + falling(i, a) * falling(j, b) * q ** (i - a) * xi ** (j - b)
    end do

  end function quartic_mixed

  subroutine quartic_partial_action(this, input_graph, input_data, slots, &
       & directions, output)

    class(quartic_form), intent(in)                :: this
    class(directed_graph), intent(in)              :: input_graph
    class(graph_field), intent(in)                 :: input_data(:)
    integer           , intent(in)                 :: slots(:)
    class(graph_field), intent(in)                 :: directions(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: xi, product
    integer  :: a, b

    call read_arguments(input_data, this % xi_default, q, xi)
    call count_seats(slots, directions, a, b, product)

    call pack_output(input_graph, &
         & [quartic_mixed(a, b, q(1), xi) * product], output)

  end subroutine quartic_partial_action

  !===================================================================!
  ! The power eight.
  !===================================================================!

  pure function power8_name(this) result(name)
    class(power8_form), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'toy power eight'
  end function power8_name

  subroutine power8_domain(this, input_graph, domain, nentries)
    class(power8_form), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    call vertex_domain(input_graph, domain, nentries)
  end subroutine power8_domain

  pure function power8_max_degree(this) result(degree)
    class(power8_form), intent(in) :: this
    integer :: degree
    associate (u1 => this); end associate
    degree = 8
  end function power8_max_degree

  subroutine power8_apply(this, input_graph, input_data, output)

    class(power8_form), intent(in)                  :: this
    class(directed_graph), intent(in)               :: input_graph
    class(graph_field), intent(in), optional        :: input_data(:)
    class(graph_field), allocatable, intent(inout)  :: output

    real(dp), allocatable :: q(:)
    real(dp) :: xi, phi

    phi = 0.0_dp
    if (present(input_data)) then
       call read_arguments(input_data, this % xi_default, q, xi)
       phi = (q(1) + xi) ** 8
    end if

    call pack_output(input_graph, [phi], output)

  end subroutine power8_apply

  subroutine power8_partial_action(this, input_graph, input_data, slots, &
       & directions, output)

    class(power8_form), intent(in)                 :: this
    class(directed_graph), intent(in)              :: input_graph
    class(graph_field), intent(in)                 :: input_data(:)
    integer           , intent(in)                 :: slots(:)
    class(graph_field), intent(in)                 :: directions(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:)
    real(dp) :: xi, product, term
    integer  :: a, b, order

    call read_arguments(input_data, this % xi_default, q, xi)
    call count_seats(slots, directions, a, b, product)

    ! d^k (q + xi)^8 = 8!/(8-k)! (q + xi)^(8-k) for every mix of
    ! q and xi slots of total order k
    order = a + b
    term  = falling(8, order) * (q(1) + xi) ** (8 - order)

    call pack_output(input_graph, [term * product], output)

  end subroutine power8_partial_action

  !===================================================================!
  ! The equilibrium law.
  !===================================================================!

  pure function equilibrium_name(this) result(name)
    class(equilibrium_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'toy equilibrium'
  end function equilibrium_name

  subroutine equilibrium_domain(this, input_graph, domain, nentries)
    class(equilibrium_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    call vertex_domain(input_graph, domain, nentries)
  end subroutine equilibrium_domain

  pure function equilibrium_max_degree(this) result(degree)
    class(equilibrium_law), intent(in) :: this
    integer :: degree
    associate (u1 => this); end associate
    degree = 8
  end function equilibrium_max_degree

  subroutine equilibrium_apply(this, input_graph, input_data, output)

    class(equilibrium_law), intent(in)              :: this
    class(directed_graph), intent(in)               :: input_graph
    class(graph_field), intent(in), optional        :: input_data(:)
    class(graph_field), allocatable, intent(inout)  :: output

    real(dp), allocatable :: q(:), s(:)
    real(dp) :: xi

    if (present(input_data)) then
       call read_arguments(input_data, this % xi_default, q, xi)
       s = q * q - xi
    else
       allocate(s(input_graph % num_vertices()))
       s = 0.0_dp
    end if

    call pack_output(input_graph, s, output)

  end subroutine equilibrium_apply

  subroutine equilibrium_partial_action(this, input_graph, input_data, &
       & slots, directions, output)

    class(equilibrium_law), intent(in)             :: this
    class(directed_graph), intent(in)              :: input_graph
    class(graph_field), intent(in)                 :: input_data(:)
    integer           , intent(in)                 :: slots(:)
    class(graph_field), intent(in)                 :: directions(:)
    class(graph_field), allocatable, intent(inout) :: output

    real(dp), allocatable :: q(:), s(:), v1(:), v2(:), w(:)
    real(dp) :: xi
    integer  :: a, b, j

    call read_arguments(input_data, this % xi_default, q, xi)

    a = count(slots == 1)
    b = count(slots == 2)

    ! collect up to two state direction vectors whole; only the
    ! first value of a xi direction is used
    do j = 1, size(slots)
       if (slots(j) == 1) then
          if (.not. allocated(v1)) then
             call directions(j) % get_real_vector(v1)
          else
             call directions(j) % get_real_vector(v2)
          end if
       else
          call directions(j) % get_real_vector(w)
       end if
    end do

    ! the three nonzero partials of q^2 - xi; every other slot
    ! combination is zero
    allocate(s(size(q)))
    if (a == 1 .and. b == 0) then
       s = 2.0_dp * q * v1
    else if (a == 2 .and. b == 0) then
       s = 2.0_dp * v1 * v2
    else if (a == 0 .and. b == 1) then
       s = -w(1)
    else
       s = 0.0_dp
    end if

    call pack_output(input_graph, s, output)

  end subroutine equilibrium_partial_action

  !===================================================================!
  ! The linear law.
  !===================================================================!

  pure function linear_name(this) result(name)
    class(linear_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u1 => this); end associate
    name = 'toy linear'
  end function linear_name

  subroutine linear_domain(this, input_graph, domain, nentries)
    class(linear_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries
    associate (u1 => this); end associate
    call vertex_domain(input_graph, domain, nentries)
  end subroutine linear_domain

  subroutine linear_apply(this, input_graph, input_data, output)

    class(linear_law), intent(in)                   :: this
    class(directed_graph), intent(in)               :: input_graph
    class(graph_field), intent(in), optional        :: input_data(:)
    class(graph_field), allocatable, intent(inout)  :: output

    real(dp), allocatable :: s(:)

    associate (u1 => this); end associate

    if (present(input_data)) then
       call input_data(1) % get_real_vector(s)
    else
       allocate(s(input_graph % num_vertices()))
       s = 0.0_dp
    end if

    call pack_output(input_graph, s, output)

  end subroutine linear_apply

end module toy_differentiable_forms
