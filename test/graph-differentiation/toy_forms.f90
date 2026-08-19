!=====================================================================!
! THE DIFFERENTIABLE TOYS: the first concretions of the calculus
! contract, one fixture set serving the whole differentiation
! tower - the chain rule, the exact linearization, the step table,
! and the marcher's three derivative walks.
!
!      quartic      Phi(q, xi) = q^4 + q^3 xi + q^2 xi^2
!                              + q xi^3 + xi^4          degree 4
!      power8       Phi(q, xi) = (q + xi)^8             degree 8
!      equilibrium  S_i(q, xi) = q_i^2 - xi             degree 8,
!                   every partial past the second exactly zero
!      linear       S(q) = q                            not
!                   differentiable - the difference road's witness
!
! Every mixed partial is COMPUTED, never tabulated: the quartic
! walks its five monomials with falling factorials, the power8
! collapses to one falling factorial of u = q + xi, and the
! equilibrium law's calculus is three live cases and zero. The
! second input slot carries xi; a statement reached with the state
! alone reads its own default, so the same toy serves the chain
! rule's two-channel compositions and the linearization family's
! same-domain tangent.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module toy_differentiable_forms

  use iso_fortran_env     , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_calculus      , only : differentiable_operation
  use fractal_graph       , only : set_graph => graph
  use class_graph_field   , only : field
  use class_graph_chain_rule, only : chain_channel

  implicit none

  private
  public :: quartic_form, power8_form, equilibrium_law, linear_law
  public :: scalar_pair, fill_scalar_channel

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
  ! Phi = (q + xi)^8: with u = q + xi, every mixed partial of order
  ! k in any mix of the two slots is the one number 8!/(8-k)!
  ! u^(8-k) - the degree-8 witness whose path derivatives a test
  ! can recompute independently by Taylor convolution.
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
  ! S_i(q, xi) = q_i^2 - xi, elementwise: the marched statement.
  ! At xi = 1 its equilibrium is q = 1, so the walk stands still
  ! while every derivative of the walk stays alive - the chain the
  ! directional and adjoint oracles are priced on. Declared to
  ! degree 8 because the declaration is a promise to ANSWER, and
  ! zero is an exact answer for every partial past the second.
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
  ! S(q) = q, an ordinary operation and deliberately NOT a
  ! differentiable one: tangent_of must hand it the difference
  ! road, and the nonuniform governed march prices 1/3 on it.
  !===================================================================!

  type, extends(graph_operation) :: linear_law
   contains
     procedure :: name   => linear_name
     procedure :: domain => linear_domain
     procedure :: apply  => linear_apply
  end type linear_law

contains

  !===================================================================!
  ! The one domain every toy answers: the graph's own vertices.
  !===================================================================!

  subroutine vertex_domain(input_graph, domain, nentries)

    class(directed_graph), intent(in)  :: input_graph
    type(set_graph)      , intent(out) :: domain
    integer              , intent(out) :: nentries

    domain   = input_graph % all_vertices()
    nentries = input_graph % num_vertices()

  end subroutine vertex_domain

  !===================================================================!
  ! The two arguments as the toys read them: the state whole from
  ! slot one, xi from slot two when the caller supplies it and the
  ! statement's own default otherwise - one reading discipline for
  ! the chain rule's two-channel inputs and the linearization
  ! family's state-only reach.
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
  ! Land an answer on the graph's own vertices.
  !===================================================================!

  subroutine pack_answer(input_graph, values, output)

    class(directed_graph), intent(in) :: input_graph
    real(dp), intent(in) :: values(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph) :: cells
    type(field)     :: out
    integer         :: nentries

    cells    = input_graph % vertex_set()
    nentries = input_graph % num_vertices()

    out = field('toy answer', cells, nentries, &
         & ncomp = size(values) / nentries)
    call out % set_real_vector(values)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine pack_answer

  !===================================================================!
  ! The falling factorial n (n-1) ... (n-k+1): the k-th derivative
  ! of a power is one of these times the lowered power, and it is
  ! zero exactly when the derivative outruns the exponent.
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
  ! How many of the requested seats perturb each slot, and the
  ! product of the scalar direction values - the whole bookkeeping
  ! a scalar toy's mixed partial needs.
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
  ! One scalar input pair (q, xi) as fields on the given cells -
  ! the composition point both test programs stand at.
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
  ! One derivative channel: the named slot, seat k carrying the
  ! k-th scalar, every seat occupied.
  !===================================================================!

  subroutine fill_scalar_channel(channel, slot, seats, cells)

    type(chain_channel), intent(inout) :: channel
    integer            , intent(in)    :: slot
    real(dp)           , intent(in)    :: seats(:)
    type(set_graph)    , intent(in)    :: cells

    integer :: k

    channel % slot = slot
    if (allocated(channel % derivative)) deallocate(channel % derivative)
    allocate(channel % derivative(size(seats)))

    do k = 1, size(seats)
       channel % derivative(k) % occupied  = .true.
       channel % derivative(k) % direction = field('path', cells, 1, ncomp=1)
       call channel % derivative(k) % direction % set_real_vector([seats(k)])
    end do

  end subroutine fill_scalar_channel

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

    call pack_answer(input_graph, [phi], output)

  end subroutine quartic_apply

  !-------------------------------------------------------------------!
  ! The mixed partial M(a, b) of Phi, computed by walking the five
  ! monomials q^i xi^(4-i): each contributes its falling factorials
  ! times the lowered powers, and dies exactly when a derivative
  ! outruns an exponent. At (a, b) = (0, 0) this IS Phi.
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

    call pack_answer(input_graph, &
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

    call pack_answer(input_graph, [phi], output)

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

    ! every slot mix of order k answers the same one number:
    ! d^k u^8 = 8!/(8-k)! u^(8-k), u = q + xi
    order = a + b
    term  = falling(8, order) * (q(1) + xi) ** (8 - order)

    call pack_answer(input_graph, [term * product], output)

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

    call pack_answer(input_graph, s, output)

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

    ! the state directions ride whole, xi's rides by its one value
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

    ! the whole calculus of q^2 - xi: three live partials, the
    ! rest exactly zero
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

    call pack_answer(input_graph, s, output)

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

    call pack_answer(input_graph, s, output)

  end subroutine linear_apply

end module toy_differentiable_forms
