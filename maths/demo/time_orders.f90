! Time integration on the chain graph: orders of accuracy of the
! three marching rules on dq/dt = -xi q, the adjoint against the
! closed form and against a finite difference, and forward
! directional derivatives in the rate xi to order 2.
module decay_fixture
  use iso_fortran_env, only : dp => REAL64
  use operation_action, only : operation, variation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal, only : graph
  use field_stored, only : stored_field
  implicit none
  private
  public :: decay_law

  ! S(q, xi) = xi q : minus the velocity, so dq/dt = -xi q
  type, extends(operation) :: decay_law
   contains
     procedure :: name           => law_name
     procedure :: domain         => law_domain
     procedure :: apply          => law_apply
     procedure :: max_degree     => law_max_degree
     procedure :: partial_action => law_partial_action
  end type decay_law

  interface decay_law
     module procedure create
  end interface decay_law

contains

  function create() result(this)
    type(decay_law) :: this
    call this % declare_arguments(2)
  end function create

  pure function law_name(this) result(name)
    class(decay_law), intent(in) :: this
    character(len=:), allocatable :: name
    associate (u => this); end associate
    name = 'decay'
  end function law_name

  subroutine law_domain(this, input_graph, domain, num_entries)
    class(decay_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(graph), intent(out) :: domain
    integer, intent(out) :: num_entries
    associate (u => this); end associate
    domain = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()
  end subroutine law_domain

  pure function law_max_degree(this) result(degree)
    class(decay_law), intent(in) :: this
    integer :: degree
    associate (u => this); end associate
    degree = 8
  end function law_max_degree

  subroutine law_apply(this, input_graph, input_data, output)
    class(decay_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field), intent(in), optional :: input_data(:)
    class(field), allocatable, intent(inout) :: output
    real(dp), allocatable :: q(:), xi(:)
    type(stored_field) :: out
    associate (u => this); end associate
    call input_data(1) % real_vector(q)
    call input_data(2) % real_vector(xi)
    out = stored_field('decay', input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(xi(1) * q)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine law_apply

  ! D_q S [v] = xi v ; D_xi S [w] = w q ; D_q D_xi S [v, w] = w v ;
  ! every other partial is zero.
  subroutine law_partial_action(this, input_graph, input_data, variations, output)
    class(decay_law), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field), intent(in) :: input_data(:)
    type(variation), intent(in) :: variations(:)
    class(field), allocatable, intent(inout) :: output
    real(dp), allocatable :: q(:), xi(:), v(:), w(:), y(:)
    type(stored_field) :: out
    integer :: n_q, n_xi, k

    call this % require_owned(variations)
    call input_data(1) % real_vector(q)
    call input_data(2) % real_vector(xi)
    allocate(y(size(q)))
    y = 0.0_dp

    n_q = 0; n_xi = 0
    do k = 1, size(variations)
       if (variations(k) % argument_is(this % argument(1))) then
          n_q = n_q + 1;  call variations(k) % direction(v)
       else
          n_xi = n_xi + 1; call variations(k) % direction(w)
       end if
    end do

    if (n_q == 1 .and. n_xi == 0) y = xi(1) * v
    if (n_q == 0 .and. n_xi == 1) y = w(1) * q
    if (n_q == 1 .and. n_xi == 1) y = w(1) * v

    out = stored_field('decay partial', input_graph % vertex_set(), input_graph % num_vertices())
    call out % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine law_partial_action

end module decay_fixture

program time_orders
  use iso_fortran_env, only : dp => REAL64
  use view_directed_stored, only : stored_directed_graph
  use graph_fractal, only : graph
  use field_stored, only : stored_field
  use operation_chain_rule, only : argument_path
  use operation_marching, only : marcher, MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2
  use operation_newton, only : newton
  use operation_gmres, only : gmres
  use decay_fixture, only : decay_law
  implicit none

  type(stored_directed_graph) :: lone
  type(graph) :: cells
  type(decay_law) :: law
  type(marcher) :: clock
  type(stored_field) :: xif
  type(argument_path) :: xipath(1)
  real(dp), allocatable :: trajectory(:,:), sens(:,:,:)
  real(dp) :: q(1), lambda(1), xi, h, err(3, 6), e_plus(1), e_minus(1), fd, exact1, exact2
  integer :: rule, k, n, nsteps(6), rules(3)
  character(len=16) :: names(3)

  nsteps = [10, 20, 40, 80, 160, 320]
  names  = ['forward Euler ', 'backward Euler', 'BDF2          ']
  rules  = [MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2]
  xi     = 1.0_dp

  lone  = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
  cells = lone % vertex_set()
  law   = decay_law()
  xif   = stored_field('xi', cells, 1); call xif % set_real_vector([xi])

  allocate(clock % inner, source=newton())
  select type (nw => clock % inner)
  type is (newton)
     allocate(nw % inner, source=gmres())
     nw % inner % tolerance = 1.0d-14
     nw % tolerance = 1.0d-13
  end select

  ! ---- orders of accuracy on q(1) = exp(-1)
  do rule = 1, 3
     clock % rule = rules(rule)
     do k = 1, 6
        clock % step = 1.0_dp / nsteps(k)
        q = [1.0_dp]
        call clock % march(law, lone, q, nsteps(k), parameters=[xif])
        err(rule, k) = abs(q(1) - exp(-xi))
     end do
  end do

  print '(a)', 'TABLE orders'
  print '(a)', 'rule & N & error & order'
  do rule = 1, 3
     do k = 1, 6
        if (k == 1) then
           print '(a,a,i4,a,es10.2,a)', trim(names(rule)), ' & ', nsteps(k), ' & ', err(rule, k), ' & -- \\'
        else
           print '(a,a,i4,a,es10.2,a,f5.2,a)', trim(names(rule)), ' & ', nsteps(k), ' & ', err(rule, k), &
                & ' & ', log(err(rule, k - 1) / err(rule, k)) / log(2.0_dp), ' \\'
        end if
     end do
  end do

  ! ---- adjoint: seed 1 at the last instant, lambda at the first = dq_N/dq_0
  print '(a)', 'TABLE adjoint'
  do rule = 2, 3
     clock % rule = rules(rule)
     n = 16; h = 1.0_dp / n; clock % step = h
     q = [1.0_dp]
     call clock % march(law, lone, q, n, trajectory=trajectory, parameters=[xif])
     lambda = [1.0_dp]
     call clock % march_adjoint(law, lone, lambda, n, trajectory, parameters=[xif])
     ! finite difference of the march in q0
     q = [1.0_dp + 1.0d-6]; call clock % march(law, lone, q, n, parameters=[xif]); e_plus = q
     q = [1.0_dp - 1.0d-6]; call clock % march(law, lone, q, n, parameters=[xif]); e_minus = q
     fd = (e_plus(1) - e_minus(1)) / 2.0d-6
     if (rule == 2) then
        exact1 = (1.0_dp + h * xi) ** (-n)
        print '(a,a,es12.5,a,es12.5,a,es9.2,a,es9.2,a)', trim(names(rule)), ' & ', lambda(1), ' & ', fd, &
             & ' & ', abs(lambda(1) - exact1), ' & ', abs(lambda(1) - fd), ' \\'
     else
        print '(a,a,es12.5,a,es12.5,a,a,a,es9.2,a)', trim(names(rule)), ' & ', lambda(1), ' & ', fd, &
             & ' & ', '--', ' & ', abs(lambda(1) - fd), ' \\'
     end if
  end do

  ! ---- directional derivatives in xi to order 2, backward Euler, closed form
  print '(a)', 'TABLE directional'
  clock % rule = MARCH_BACKWARD
  n = 16; h = 1.0_dp / n; clock % step = h
  q = [1.0_dp]
  call clock % march(law, lone, q, n, trajectory=trajectory, parameters=[xif])
  xipath(1) % wrt = law % argument(2)
  allocate(xipath(1) % derivative(2))
  xipath(1) % derivative(1) % occupied = .true.
  xipath(1) % derivative(1) % direction = stored_field('path', cells, 1)
  call xipath(1) % derivative(1) % direction % set_real_vector([1.0_dp])
  xipath(1) % derivative(2) % occupied = .true.
  xipath(1) % derivative(2) % direction = stored_field('path', cells, 1)
  call xipath(1) % derivative(2) % direction % set_real_vector([0.0_dp])
  call clock % march_directional(law, lone, n, trajectory, 2, sens, parameters=[xif], paths=xipath)
  exact1 = -n * h * (1.0_dp + h * xi) ** (-n - 1)
  exact2 =  n * (n + 1) * h * h * (1.0_dp + h * xi) ** (-n - 2)
  print '(a,es14.7,a,es14.7,a,es9.2,a)', 'order 1 & ', sens(1, 1, n + 1), ' & ', exact1, ' & ', abs(sens(1, 1, n + 1) - exact1), ' \\'
  print '(a,es14.7,a,es14.7,a,es9.2,a)', 'order 2 & ', sens(1, 2, n + 1), ' & ', exact2, ' & ', abs(sens(1, 2, n + 1) - exact2), ' \\'

end program time_orders
