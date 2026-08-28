module gti_configuration
  use util_precision  , only : dp
  implicit none
  private
  public :: configuration, read_configuration, override, show
  public :: worded, lists, refuse_unknown, chosen_from
  ! the levels of this application's hierarchy, outermost first: what
  ! the tally files an amount under
  integer, parameter, public :: at_expansion = 1
  integer, parameter, public :: at_horizon   = 2
  integer, parameter, public :: at_block     = 3
  integer, parameter, public :: at_stage     = 4
  character(len=9), parameter, public :: hierarchy_levels(4) = &
       & ['expansion', 'horizon  ', 'block    ', 'stage    ']
  type :: configuration
     character(len=32) :: physics      = 'vanderpol'
     character(len=32) :: grid         = 'random'
     character(len=32) :: combinations = 'homogeneous'
     character(len=32) :: families     = 'bdf adams dirk'
     integer  :: state_degree             = 2
     integer  :: instants                 = 21
     integer  :: max_derivative_degree    = 3
     integer  :: max_discretization_order = 4
     integer  :: seed                     = 20260824
     integer  :: startup_refinement       = 4
     character(len=16) :: linear_solver   = 'direct'
     character(len=16) :: assembly        = 'matrix'
     character(len=16) :: storage         = 'dense'
     logical           :: multigrid       = .false.
     character(len=16) :: sweep           = 'space-time'
     character(len=64) :: designs         = 'physics'
     character(len=64) :: functionals     = 'energy'
     real(dp) :: time_duration = 7.0_dp
     real(dp) :: design        = 1.0_dp
     logical  :: automatic_order_conservation = .true.
     logical  :: mixed_orders                 = .false.
     character(len=128) :: initial_state = '1.0'
     character(len=16)  :: spatial_geometry = 'cartesian'
     character(len=64)  :: spatial_extent   = '1.0 1.0'
     character(len=64)  :: spatial_counts   = '0 0'
     character(len=16)  :: spatial_grid     = 'uniform'
     real(dp)           :: diffusion        = 0.0_dp
     integer            :: spatial_order    = 1
     character(len=16)  :: initial_field    = 'constant'
     character(len=16)  :: export           = 'none'
     character(len=128) :: export_path      = 'field'
     character(len=16)  :: check            = 'none'
     real(dp)          :: tolerance           = 1.0e-12_dp
     character(len=16) :: tolerance_criterion = 'relative'
     character(len=16) :: iteration_criterion = 'by_rate'
     integer           :: max_iterations      = 100
     integer           :: higher_order_jacobian_product = 1
     integer           :: krylov_restart        = 60
     integer           :: smoothing_sweeps      = 2
     integer           :: max_linear_iterations = 200
     logical  :: accounting                   = .false.
     character(len=256) :: measurements       = &
          & 'wall_time primal_loops tangent_loops adjoint_loops ' // &
          & 'newton_solves linear_solves factorisations'
  end type configuration
contains
  pure function levelled(text) result(name)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: name
    integer :: i, c
    name = trim(adjustl(text))
    do i = 1, len(name)
       if (name(i:i) == '-') name(i:i) = '_'
       c = iachar(name(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) name(i:i) = achar(c + 32)
    end do
  end function levelled
  pure function worded(text) result(list)
    character(len=*), intent(in) :: text
    character(len=32), allocatable :: list(:)
    character(len=32) :: held(32)
    integer :: i, first, n, last
    n    = 0
    i    = 1
    last = len_trim(text)
    do while (i <= last)
       if (text(i:i) == ' ') then
          i = i + 1
          cycle
       end if
       first = i
       do while (i <= last)
          if (text(i:i) == ' ') exit
          i = i + 1
       end do
       if (n == size(held)) exit
       n = n + 1
       held(n) = text(first:i-1)
    end do
    list = held(1:n)
  end function worded
  pure logical function lists(text, what) result(yes)
    character(len=*), intent(in) :: text, what
    character(len=32), allocatable :: list(:)
    integer :: i
    list = worded(text)
    yes  = .false.
    do i = 1, size(list)
       if (trim(list(i)) == what) yes = .true.
    end do
  end function lists
  subroutine refuse_unknown(text, every, subject)
    character(len=*), intent(in) :: text, every(:), subject
    character(len=32), allocatable :: list(:)
    integer :: i, j
    logical :: known
    list = worded(text)
    do i = 1, size(list)
       known = .false.
       do j = 1, size(every)
          if (trim(list(i)) == trim(every(j))) known = .true.
       end do
       if (.not. known) then
          write(*,'(a)') ' '
          write(*,'(a)') ' ' // subject // ' names ' // trim(list(i)) // &
               & ', which this program has nothing for.'
          error stop 'gti_configuration: a setting names something unknown'
       end if
    end do
  end subroutine refuse_unknown
  subroutine assign(cfg, name, value)
    type(configuration), intent(inout) :: cfg
    character(len=*)   , intent(in)    :: name, value
    select case (levelled(name))
    case ('physics')
       cfg % physics = value
    case ('grid')
       cfg % grid = value
    case ('combinations')
       cfg % combinations = value
    case ('families')
       cfg % families = value
    case ('state_degree')
       read(value, *) cfg % state_degree
    case ('instants')
       read(value, *) cfg % instants
    case ('max_derivative_degree')
       read(value, *) cfg % max_derivative_degree
    case ('max_discretization_order')
       read(value, *) cfg % max_discretization_order
    case ('seed')
       read(value, *) cfg % seed
    case ('startup_refinement')
       read(value, *) cfg % startup_refinement
    case ('linear_solver')
       cfg % linear_solver = value
    case ('sweep')
       cfg % sweep = value
    case ('designs')
       cfg % designs = value
    case ('functionals')
       cfg % functionals = value
    case ('assembly')
       cfg % assembly = value
    case ('storage')
       cfg % storage = value
    case ('multigrid')
       read(value, *) cfg % multigrid
    case ('time_duration')
       read(value, *) cfg % time_duration
    case ('design')
       read(value, *) cfg % design
    case ('automatic_order_conservation')
       read(value, *) cfg % automatic_order_conservation
    case ('mixed_orders')
       read(value, *) cfg % mixed_orders
    case ('initial_state')
       cfg % initial_state = value
    case ('spatial_geometry')
       cfg % spatial_geometry = value
    case ('spatial_extent')
       cfg % spatial_extent = value
    case ('spatial_counts')
       cfg % spatial_counts = value
    case ('spatial_grid')
       cfg % spatial_grid = value
    case ('diffusion')
       read(value, *) cfg % diffusion
    case ('spatial_order')
       read(value, *) cfg % spatial_order
    case ('initial_field')
       cfg % initial_field = value
    case ('export')
       cfg % export = value
    case ('export_path')
       cfg % export_path = value
    case ('check')
       cfg % check = value
    case ('tolerance')
       read(value, *) cfg % tolerance
    case ('tolerance_criterion')
       cfg % tolerance_criterion = value
    case ('iteration_criterion')
       cfg % iteration_criterion = value
    case ('krylov_restart')
       read(value, *) cfg % krylov_restart
    case ('smoothing_sweeps')
       read(value, *) cfg % smoothing_sweeps
    case ('max_linear_iterations')
       read(value, *) cfg % max_linear_iterations
    case ('max_iterations')
       read(value, *) cfg % max_iterations
    case ('higher_order_jacobian_product')
       read(value, *) cfg % higher_order_jacobian_product
    case ('accounting')
       read(value, *) cfg % accounting
    case ('measurements')
       cfg % measurements = value
    case default
       write(*,'(a)') ' this is not a setting: ' // levelled(name)
       error stop 'gti_configuration: every setting given is one that exists'
    end select
    call refuse_below(cfg % max_iterations, 1, &
         & 'max_iterations', 'an iteration budget is at least one')
    call refuse_below(cfg % higher_order_jacobian_product, 1, &
         & 'higher_order_jacobian_product', 'one is Newton, and there is no order below it')
    call refuse_below(cfg % spatial_order, 1, &
         & 'spatial_order', 'a form of degree below one fits no gradient')
    call refuse_below(cfg % max_derivative_degree, 0, &
         & 'max_derivative_degree', 'the value on its own is degree zero')
    call refuse_below(cfg % max_discretization_order, 1, &
         & 'max_discretization_order', 'no scheme is built below order one')
  end subroutine assign
  subroutine refuse_below(given, least, name, why)
    integer         , intent(in) :: given, least
    character(len=*), intent(in) :: name, why
    if (given >= least) return
    write(*,'(a)') ' '
    write(*,'(a,i0,a,i0,a)') ' ' // name // ' is ', given, ', under ', least, &
         & ': ' // why // '.'
    error stop 'gti_configuration: a setting is under the value it means anything at'
  end subroutine refuse_below
  subroutine take_line(cfg, line)
    type(configuration), intent(inout) :: cfg
    character(len=*)   , intent(in)    :: line
    character(len=:), allocatable :: text
    integer :: at
    text = trim(adjustl(line))
    at   = index(text, '#')
    if (at > 0) text = trim(text(:at - 1))
    if (len(text) == 0) return
    at = index(text, '=')
    if (at == 0) then
       write(*,'(a)') ' this line names no setting: ' // text
       error stop 'gti_configuration: a setting is named, then given its value'
    end if
    call assign(cfg, text(:at - 1), trim(adjustl(text(at + 1:))))
  end subroutine take_line
  subroutine read_configuration(name, cfg)
    character(len=*)   , intent(in)  :: name
    type(configuration), intent(out) :: cfg
    character(len=512) :: line
    integer :: unit, status
    open(newunit=unit, file='config/' // trim(name) // '.cfg', &
         & status='old', action='read', iostat=status)
    if (status /= 0) then
       write(*,'(a)') ' no such configuration: config/' // trim(name) // '.cfg'
       error stop 'gti_configuration: a configuration names a file that exists'
    end if
    do
       read(unit, '(a)', iostat=status) line
       if (status /= 0) exit
       call take_line(cfg, line)
    end do
    close(unit)
  end subroutine read_configuration
  subroutine override(cfg, argument)
    type(configuration), intent(inout) :: cfg
    character(len=*)   , intent(in)    :: argument
    character(len=:), allocatable :: text
    integer :: at
    text = trim(adjustl(argument))
    if (len(text) > 2) then
       if (text(1:2) == '--') text = text(3:)
    end if
    at = index(text, '=')
    if (at == 0) then
       write(*,'(a)') ' this argument names no setting: ' // text
       error stop 'gti_configuration: an argument is a setting and its value'
    end if
    call assign(cfg, text(:at - 1), trim(adjustl(text(at + 1:))))
  end subroutine override
  subroutine show(cfg)
    type(configuration), intent(in) :: cfg
    write(*,'(a)')         ' the run'
    write(*,'(a,a)')       '   physics                  ', trim(cfg % physics)
    write(*,'(a,i0)')      '   state degree             ', cfg % state_degree
    write(*,'(a,f0.4)')    '   time duration            ', cfg % time_duration
    write(*,'(a,i0)')      '   instants                 ', cfg % instants
    write(*,'(a,a)')       '   grid                     ', trim(cfg % grid)
    write(*,'(a,i0)')      '   seed                     ', cfg % seed
    write(*,'(a,i0)')      '   startup refinement       ', cfg % startup_refinement
    write(*,'(a,a)')       '   linear solver            ', trim(cfg % linear_solver)
    write(*,'(a,a)')       '   sweep                    ', trim(cfg % sweep)
    write(*,'(a,a)')       '   designs                  ', trim(cfg % designs)
    write(*,'(a,a)')       '   functionals              ', trim(cfg % functionals)
    write(*,'(a,a)')       '   assembly                 ', trim(cfg % assembly)
    write(*,'(a,a)')       '   storage                  ', trim(cfg % storage)
    write(*,'(a,l1)')      '   multigrid                ', cfg % multigrid
    write(*,'(a,f0.4)')    '   design                   ', cfg % design
    write(*,'(a,i0)')      '   max derivative degree    ', cfg % max_derivative_degree
    write(*,'(a,i0)')      '   max discretization order ', cfg % max_discretization_order
    write(*,'(a,a)')       '   families                 ', trim(cfg % families)
    write(*,'(a,a)')       '   combinations             ', trim(cfg % combinations)
    write(*,'(a,l1)')      '   automatic order conservation ', cfg % automatic_order_conservation
    write(*,'(a,l1)')      '   mixed orders             ', cfg % mixed_orders
    write(*,'(a,a)')       '   initial state, given     ', trim(cfg % initial_state)
    if (trim(cfg % spatial_counts) /= '0 0') then
       write(*,'(a,a)')    '   spatial geometry         ', trim(cfg % spatial_geometry)
       write(*,'(a,a)')    '   spatial extent           ', trim(cfg % spatial_extent)
       write(*,'(a,a)')    '   spatial counts           ', trim(cfg % spatial_counts)
       write(*,'(a,a)')    '   spatial grid             ', trim(cfg % spatial_grid)
       write(*,'(a,es9.2)')'   diffusion                ', cfg % diffusion
       write(*,'(a,i0)')   '   spatial order            ', cfg % spatial_order
       write(*,'(a,a)')    '   initial field            ', trim(cfg % initial_field)
       write(*,'(a,a)')    '   export                   ', trim(cfg % export)
       write(*,'(a,a)')    '   check                    ', trim(cfg % check)
    end if
    write(*,'(a,es9.2)')   '   tolerance                ', cfg % tolerance
    write(*,'(a,a)')       '   tolerance criterion      ', trim(cfg % tolerance_criterion)
    write(*,'(a,a)')       '   iteration criterion      ', trim(cfg % iteration_criterion)
    write(*,'(a,i0)')      '   max iterations           ', cfg % max_iterations
    write(*,'(a,i0)')      '   higher order jacobian product ', cfg % higher_order_jacobian_product
    write(*,'(a,i0)')      '   krylov restart           ', cfg % krylov_restart
    write(*,'(a,i0)')      '   smoothing sweeps         ', cfg % smoothing_sweeps
    write(*,'(a,i0)')      '   max linear iterations    ', cfg % max_linear_iterations
    write(*,'(a,l1)')      '   accounting               ', cfg % accounting
    if (cfg % accounting) then
       write(*,'(a,a)')    '   measurements             ', trim(cfg % measurements)
    end if
  end subroutine show
  integer function chosen_from(text, every, subject) result(which)
    character(len=*), intent(in) :: text, every(:), subject
    integer :: j
    call refuse_unknown(text, every, subject)
    which = 0
    do j = 1, size(every)
       if (trim(text) == trim(every(j))) which = j
    end do
  end function chosen_from
end module gti_configuration

!=====================================================================!
! The equations this application integrates, each stated once as an
! expression: the van der Pol residual at any degree, its energy, and
! the power its damping draws. The expression machinery knows none of
! them; they are this application's physics, so they live here.
!=====================================================================!

module gti_physics
  use util_precision      , only : dp
  use operation_expression, only : expression, unknown, design, derivative, stated, &
       & operator(+), operator(-), operator(*), operator(**)
  implicit none
  private
  public :: van_der_pol, van_der_pol_energy, van_der_pol_dissipation
contains
  function van_der_pol(degree) result(r)
    integer, intent(in) :: degree
    type(expression) :: r
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    r = stated(derivative(q, degree) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, degree - 1) &
         & + derivative(q, 0), degree, 'van der pol residual')
  end function van_der_pol

  function van_der_pol_energy(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    type(expression) :: q
    q = unknown()
    f = stated(0.5_dp * (derivative(q, 0)**2 + derivative(q, 1)**2), degree, 'van der pol energy')
  end function van_der_pol_energy

  function van_der_pol_dissipation(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    f = stated(nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, 1) * derivative(q, 1), &
         & degree, 'van der pol dissipation')
  end function van_der_pol_dissipation
end module gti_physics
module gti_sweeps
  use gti_configuration, only : refuse_unknown
  use util_precision  , only : dp, half_digits
  use graph_fractal   , only : graph
  use view_directed   , only : directed_graph
  use field_stored    , only : stored_field
  use operation_action, only : operation, applied, varied
  use operation_dense_direct, only : dense_direct
  use operation_multigrid   , only : multigrid
  use operation_gauss_seidel, only : gauss_seidel
  use operation_gmres       , only : gmres
  use operation_minimization, only : minimizer, relative, absolute, by_rate, by_count
  implicit none
  private
  public :: forward_route, reverse_route, route_of, route_substitutions, choose
  integer, parameter :: forward_route = 1
  integer, parameter :: reverse_route = 2
  public :: set_linear_solver, set_assembly, set_storage, set_multigrid
  public :: set_newton_order, newton_order
  public :: set_aggregates, set_coarse_nodes, coarse_nodes, assembly_present, multigrid_on
  public :: take_inner, keep_inner, forget_inner, set_linear_stopping, set_linear_budget
  public :: functional_of, functional_gradient
  character(len=16), save :: chosen_solver   = 'direct'
  character(len=16), save :: chosen_assembly = 'matrix'
  integer          , save :: chosen_newton_order = 1
  character(len=16), save :: chosen_storage  = 'dense'
  logical          , save :: chosen_multigrid = .false.
  integer, allocatable, save :: chosen_aggregates(:)
  real(dp), save :: linear_tolerance  = half_digits
  integer , save :: linear_criterion  = relative
  integer , save :: linear_budget     = by_rate
  integer , save :: linear_restart    = 60
  integer , save :: linear_sweeps     = 2
  integer , save :: linear_iterations = 200
  integer, allocatable, save :: chosen_coarse(:)
  class(minimizer), allocatable, save :: kept_inner
contains
  real(dp) function functional_of(integrand, instants, inputs, dt) result(f)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)

    real(dp), allocatable :: values(:)

    call applied(integrand, instants, inputs, values)
    f = sum(dt * values)

  end function functional_of

  subroutine functional_gradient(integrand, instants, inputs, dt, n, degrees, &
       & state_domain, g, along_state, along_design)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)
    integer              , intent(in) :: n, degrees
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: g(:)
    real(dp), intent(in), optional    :: along_state(:), along_design(:)

    real(dp), allocatable :: v(:), rate(:)
    integer :: d, k

    allocate(g(n * degrees), source=0.0_dp)
    allocate(v(n * degrees))
    do d = 0, degrees - 1
       v = 0.0_dp
       do k = 1, n
          v((k - 1) * degrees + d + 1) = 1.0_dp
       end do
       if (present(along_state)) then
          call varied(integrand, instants, inputs, 1, state_domain, v, rate, &
               & 1, state_domain, along_state)
       else if (present(along_design)) then
          call varied(integrand, instants, inputs, 1, state_domain, v, rate, &
               & 2, state_domain, along_design)
       else
          call varied(integrand, instants, inputs, 1, state_domain, v, rate)
       end if
       do k = 1, n
          g((k - 1) * degrees + d + 1) = dt(k) * rate(k)
       end do
    end do

  end subroutine functional_gradient
  subroutine set_linear_solver(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['direct   ', 'iterative'], 'linear_solver')
    chosen_solver = name
    call forget_inner()
  end subroutine set_linear_solver
  subroutine set_assembly(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['matrix', 'free  '], 'assembly')
    chosen_assembly = name
    call forget_inner()
  end subroutine set_assembly
  subroutine set_storage(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['dense ', 'sparse'], 'storage')
    chosen_storage = name
    call forget_inner()
  end subroutine set_storage
  subroutine set_multigrid(on)
    logical, intent(in) :: on
    chosen_multigrid = on
    call forget_inner()
  end subroutine set_multigrid
  subroutine set_newton_order(order)
    integer, intent(in) :: order
    if (order < 1) error stop 'gti_sweeps: one is Newton, and there is no order below it'
    chosen_newton_order = order
  end subroutine set_newton_order
  pure integer function newton_order() result(order)
    order = chosen_newton_order
  end function newton_order
  pure logical function assembly_present() result(yes)
    yes = trim(chosen_assembly) == 'matrix'
  end function assembly_present
  pure logical function multigrid_on() result(yes)
    yes = chosen_multigrid
  end function multigrid_on
  subroutine set_aggregates(aggregates)
    integer, intent(in), optional :: aggregates(:)
    if (allocated(chosen_aggregates)) deallocate(chosen_aggregates)
    if (present(aggregates)) chosen_aggregates = aggregates
  end subroutine set_aggregates
  subroutine set_coarse_nodes(cell)
    integer, intent(in), optional :: cell(:)
    if (allocated(chosen_coarse)) deallocate(chosen_coarse)
    if (present(cell)) chosen_coarse = cell
  end subroutine set_coarse_nodes
  function coarse_nodes(nodes) result(cell)
    integer, intent(in) :: nodes
    integer, allocatable :: cell(:)
    integer :: i
    if (allocated(chosen_coarse)) then
       if (size(chosen_coarse) < nodes) then
          error stop 'gti_sweeps: a coarse cell for every node'
       end if
       cell = chosen_coarse
    else
       cell = [(i, i = 1, nodes)]
    end if
  end function coarse_nodes
  subroutine set_linear_stopping(tolerance, criterion, budget)
    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, budget
    if (tolerance <= 0.0_dp) error stop 'gti_sweeps: a tolerance is positive'
    if (criterion /= relative .and. criterion /= absolute) then
       error stop 'gti_sweeps: a tolerance is measured relative or absolute'
    end if
    if (budget /= by_count .and. budget /= by_rate) then
       error stop 'gti_sweeps: a budget is counted or taken from the rate'
    end if
    linear_tolerance = tolerance
    linear_criterion = criterion
    linear_budget    = budget
  end subroutine set_linear_stopping
  subroutine set_linear_budget(restart, sweeps, iterations)
    integer, intent(in) :: restart, sweeps, iterations
    if (restart < 1 .or. sweeps < 1 .or. iterations < 1) then
       error stop 'gti_sweeps: a restart, a sweep count and a ceiling are positive'
    end if
    linear_restart    = restart
    linear_sweeps     = sweeps
    linear_iterations = iterations
  end subroutine set_linear_budget
  function inner_minimizer(count, width) result(inner)
    integer, intent(in) :: count, width
    class(minimizer), allocatable :: inner
    class(minimizer), allocatable :: named
    type(gmres)        :: krylov
    type(dense_direct) :: factorisation
    type(multigrid)    :: levels
    type(gauss_seidel) :: sweeps
    if (trim(chosen_assembly) == 'free' .and. trim(chosen_solver) == 'direct') then
       error stop 'gti_sweeps: a free assembly has no matrix to factorise; its solver iterates'
    end if
    if (trim(chosen_solver) == 'direct' .and. trim(chosen_storage) == 'sparse') then
       error stop 'gti_sweeps: a sparse direct solve is not built'
    end if
    if (trim(chosen_solver) == 'iterative' .and. trim(chosen_assembly) == 'matrix' &
         & .and. trim(chosen_storage) == 'dense') then
       error stop 'gti_sweeps: an iterative solve reads the sparse stencil; dense storage is for factorising'
    end if
    if (chosen_multigrid .and. trim(chosen_assembly) == 'free') then
       error stop 'gti_sweeps: multigrid coarsens a stencil, which a free assembly has not'
    end if
    select case (trim(chosen_solver))
    case ('direct')
       factorisation = dense_direct()
       factorisation % singular_reported = .true.
       allocate(named, source=factorisation)
    case ('iterative')
       krylov = gmres()
       krylov % restart        = min(count, linear_restart)
       krylov % tolerance      = linear_tolerance
       krylov % criterion      = linear_criterion
       krylov % budget         = linear_budget
       krylov % max_iterations = linear_iterations
       allocate(named, source=krylov)
    end select
    if (.not. chosen_multigrid) then
       call move_alloc(named, inner)
       return
    end if
    if (.not. allocated(chosen_aggregates)) then
       error stop 'gti_sweeps: multigrid coarsens by aggregates, and none were given'
    end if
    if (size(chosen_aggregates) /= count) then
       error stop 'gti_sweeps: one aggregate per unknown'
    end if
    sweeps % max_iterations = linear_sweeps
    sweeps % block_width    = width
    allocate(levels % smoother, source=sweeps)
    call move_alloc(named, levels % coarse)
    levels % block_width    = width
    levels % aggregates     = chosen_aggregates
    levels % tolerance      = linear_tolerance
    levels % criterion      = linear_criterion
    levels % budget         = linear_budget
    levels % max_iterations = linear_iterations
    allocate(inner, source=levels)
  end function inner_minimizer
  subroutine take_inner(inner, count, width)
    class(minimizer), allocatable, intent(out) :: inner
    integer                      , intent(in)  :: count, width
    if (allocated(kept_inner) .and. .not. chosen_multigrid) then
       call move_alloc(kept_inner, inner)
    else
       call forget_inner()
       allocate(inner, source=inner_minimizer(count, width))
    end if
  end subroutine take_inner
  subroutine keep_inner(inner)
    class(minimizer), allocatable, intent(inout) :: inner
    if (allocated(kept_inner)) deallocate(kept_inner)
    call move_alloc(inner, kept_inner)
  end subroutine keep_inner
  subroutine forget_inner()
    if (allocated(kept_inner)) deallocate(kept_inner)
  end subroutine forget_inner
  pure integer function route_of(num_designs, num_functionals, order) result(route)
    integer, intent(in) :: num_designs, num_functionals, order
    if (num_designs < 1 .or. num_functionals < 1 .or. order < 1) then
       error stop 'gti_sweeps: a route is chosen for at least one design, one functional and order one'
    end if
    if (num_designs <= order * num_functionals) then
       route = forward_route
    else
       route = reverse_route
    end if
  end function route_of
  pure integer function route_substitutions(route, num_designs, num_functionals, order) &
       & result(count)
    integer, intent(in) :: route, num_designs, num_functionals, order
    select case (route)
    case (forward_route)
       count = choose(num_designs + order - 1, order)
    case (reverse_route)
       count = (1 + num_functionals) * choose(num_designs + order - 2, order - 1)
    case default
       error stop 'gti_sweeps: a route is forward or reverse'
    end select
  end function route_substitutions
  pure integer function choose(n, k) result(c)
    integer, intent(in) :: n, k
    integer :: i
    c = 1
    do i = 1, k
       c = c * (n - k + i) / i
    end do
  end function choose
end module gti_sweeps
module gti_expansion
  use util_precision  , only : dp
  use graph_fractal         , only : graph, branch
  use view_level            , only : level_storage, level_consistent, &
       & level_is_leaf, level_members, level_couples, level_coupling
  use view_sequence         , only : sequence_empty, sequence_first, sequence_rest
  use view_relational       , only : relational_binding, relational_valid, &
       & num_relations, relation_at
  use relation_finitary     , only : relation
  use relation_binary       , only : csr_relation, binary_relation
  use map_value             , only : value_map, VALUE_UNKNOWN, VALUE_KNOWN
  use map_label             , only : label_map
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use operation_grid        , only : grid, partition
  use operation_coupling    , only : weights_of
  use operation_stencil     , only : stencil, combine_triples
  use operation_action      , only : variation
  use operation_weight      , only : scheme_weight
  use operation_expression     , only : expression
  implicit none
  private
  public :: expansion, family_holder
  public :: design_of_physics, design_of_steps
  public :: marches_by_stages
  public :: block_reach
  integer, parameter :: design_of_physics = 1
  integer, parameter :: design_of_steps   = 2
  type :: family_holder
     class(family), allocatable :: scheme
  end type family_holder
  type :: expansion
     type(level_storage)     , private :: nodes
     type(label_map)         , private :: labels
     type(value_map)         , private :: values
     type(set_map)           , private :: extents
     type(relational_binding), private :: bindings
     integer                 , private :: root_at = 0
     integer                 , private :: degrees = 0
     integer                 , private :: node_extent = 1
     integer                 , private :: spatial_coupling_at = 0
     integer, allocatable    , private :: design_at(:), design_kind(:)
     type(expression)        , private :: rule_kept
     class(grid), allocatable, private :: steps_kept
   contains
     procedure :: build
     procedure :: root
     procedure :: node
     procedure :: num_nodes
     procedure :: label_of
     procedure :: status_of
     procedure :: value_of
     procedure :: extent_of
     procedure :: consistent
     procedure :: tuples_of
     procedure :: rule
     procedure :: parameter
     procedure :: num_designs
     procedure :: design_kind_of
     procedure :: design_extent
     procedure :: design_value
     procedure :: step_partials
     procedure :: step_partial_along
     procedure, private :: weights_of_steps
     procedure, private :: refuse_assignment
     generic :: assignment(=) => refuse_assignment
  end type expansion
contains
  subroutine refuse_assignment(lhs, rhs)
    class(expansion), intent(out) :: lhs
    class(expansion), intent(in)  :: rhs
    associate (u1 => lhs, u2 => rhs); end associate
    error stop 'gti_expansion: an expansion is not assignable'
  end subroutine refuse_assignment
  pure integer function root(this)
    class(expansion), intent(in) :: this
    root = this % root_at
  end function root
  function node(this, at) result(g)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: at
    type(graph), pointer :: g
    g => this % nodes % node(at)
  end function node
  pure integer function num_nodes(this)
    class(expansion), intent(in) :: this
    num_nodes = this % nodes % num_nodes()
  end function num_nodes
  function label_of(this, g) result(text)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    character(len=:), allocatable :: text
    text = ''
    if (this % labels % labelled(g)) text = this % labels % label_of(g)
  end function label_of
  pure integer function status_of(this, g)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    status_of = this % values % status_of(g)
  end function status_of
  subroutine value_of(this, g, x)
    class(expansion)     , intent(in)  :: this
    type(graph)          , intent(in)  :: g
    real(dp), allocatable, intent(out) :: x(:)
    call this % values % value_of(g, x)
  end subroutine value_of
  integer function extent_of(this, g) result(n)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    n = 0
    if (this % extents % describes(g)) n = this % extents % num_members_of(g)
  end function extent_of
  subroutine tuples_of(this, coupling, table)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: coupling
    integer, allocatable, intent(out) :: table(:,:)
    class(relation), pointer :: r
    if (num_relations(coupling) /= 1) then
       error stop 'gti_expansion: a coupling holds one relation'
    end if
    r => relation_at(coupling, this % bindings, 1)
    select type (r)
    class is (binary_relation)
       call r % tuples(table)
    class default
       error stop 'gti_expansion: a coupling''s relation is binary'
    end select
  end subroutine tuples_of
  function in_relation_order(this, coupling, table, w) result(placed)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: coupling
    integer         , intent(in) :: table(:,:)
    real(dp)        , intent(in) :: w(:)
    real(dp), allocatable :: placed(:)
    integer, allocatable :: kept(:,:), at(:,:)
    integer :: e, n
    call this % tuples_of(coupling, kept)
    if (size(kept, 2) /= size(table, 2)) then
       error stop 'gti_expansion: a coupling names each tuple once'
    end if
    n = max(maxval(table(1, :)), maxval(table(2, :)))
    allocate(at(n, n), source=0)
    do e = 1, size(table, 2)
       at(table(1, e), table(2, e)) = e
    end do
    allocate(placed(size(kept, 2)))
    do e = 1, size(kept, 2)
       placed(e) = w(at(kept(1, e), kept(2, e)))
    end do
  end function in_relation_order
  subroutine build(this, physics, schemes, instants, steps, &
       & max_derivative_degree, parameter, nodes, spatial_discretization_stencil, weights, block_steps)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    class(grid)           , intent(in)    :: steps
    integer               , intent(in)    :: max_derivative_degree
    real(dp)              , intent(in)    :: parameter
    integer      , intent(in), optional   :: nodes
    type(stencil), intent(in), optional   :: spatial_discretization_stencil
    real(dp)     , intent(in), optional   :: weights(:), block_steps(:)
    real(dp), allocatable :: dt(:)
    integer , allocatable :: sweeps(:)
    integer :: s
    if (this % root_at /= 0) then
       error stop 'gti_expansion: an expansion is built once'
    end if
    if (size(schemes) /= size(instants)) then
       error stop 'gti_expansion: one family and one instant count per block'
    end if
    this % degrees = physics % equation_degree() + 1
    this % node_extent = 1
    if (present(nodes)) this % node_extent = nodes
    if (present(spatial_discretization_stencil)) this % spatial_coupling_at = spatial_discretization_coupling(this, spatial_discretization_stencil)
    this % rule_kept = physics
    allocate(this % steps_kept, source=steps)
    if (present(block_steps)) then
       if (size(block_steps) /= sum(instants) - 1) then
          error stop 'gti_expansion: one step per instant after the first'
       end if
       dt = [0.0_dp, block_steps]
    else if (present(weights)) then
       call partition(steps, sum(instants), weights, dt)
    else
       call partition(steps, sum(instants), [real(dp) ::], dt)
    end if
    allocate(sweeps(max_derivative_degree + 1))
    do s = 0, max_derivative_degree
       sweeps(s + 1) = one_sweep(this, physics, schemes, instants, dt, s)
    end do
    allocate(this % design_at(0), this % design_kind(0))
    call one_design(this, 'the physics'' parameter', [parameter], design_of_physics)
    if (present(weights)) call one_design(this, 'the weights of the steps', weights, design_of_steps)
    this % root_at = this % nodes % assemble([sweeps, this % nodes % assemble(this % design_at, 0)], 0)
    call this % labels % bind(this % node(this % root_at), &
         & 'expansion of ' // physics % name() // ' in the design')
  end subroutine build
  subroutine one_design(this, text, x, kind)
    class(expansion), intent(inout) :: this
    character(len=*), intent(in)    :: text
    real(dp)        , intent(in)    :: x(:)
    integer         , intent(in)    :: kind
    integer :: at
    at = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(at), text)
    call this % extents % bind(this % node(at), counted_set_representation(size(x)))
    call attach_known(this, at, x)
    this % design_at   = [this % design_at, at]
    this % design_kind = [this % design_kind, kind]
  end subroutine one_design
  pure integer function num_designs(this)
    class(expansion), intent(in) :: this
    num_designs = size(this % design_at)
  end function num_designs
  pure integer function design_kind_of(this, k)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: k
    design_kind_of = this % design_kind(k)
  end function design_kind_of
  integer function design_extent(this, k)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: k
    design_extent = this % extent_of(this % node(this % design_at(k)))
  end function design_extent
  subroutine design_value(this, k, x)
    class(expansion)     , intent(in)  :: this
    integer              , intent(in)  :: k
    real(dp), allocatable, intent(out) :: x(:)
    call this % value_of(this % node(this % design_at(k)), x)
  end subroutine design_value
  real(dp) function parameter(this)
    class(expansion), intent(in) :: this
    real(dp), allocatable :: x(:)
    call this % design_value(1, x)
    parameter = x(1)
  end function parameter
  function rule(this) result(r)
    class(expansion), intent(in) :: this
    type(expression) :: r
    r = this % rule_kept
  end function rule
  subroutine step_partial_along(this, weights_varied, u)
    class(expansion)     , intent(in)  :: this
    integer              , intent(in)  :: weights_varied(:)
    real(dp), allocatable, intent(out) :: u(:)
    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    type(stored_field), allocatable :: direction(:)
    type(variation)   , allocatable :: variations(:)
    class(field), allocatable :: out
    real(dp), allocatable :: weights(:), e(:)
    integer :: n, i
    call this % weights_of_steps(weights)
    if (any(weights_varied < 1) .or. any(weights_varied > size(weights))) then
       error stop 'gti_expansion: every weight varied is one of the grid''s'
    end if
    n = size(weights) + 1
    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), size(weights))
    call knobs % set_real_vector(weights)
    allocate(e(size(weights)), direction(size(weights_varied)), variations(size(weights_varied)))
    do i = 1, size(weights_varied)
       e = 0.0_dp
       e(weights_varied(i)) = 1.0_dp
       direction(i) = stored_field('direction', instants % vertex_set(), size(weights))
       call direction(i) % set_real_vector(e)
       variations(i) = variation(this % steps_kept % argument(1), direction(i))
    end do
    call this % steps_kept % partial_action(instants, [knobs], variations, out)
    call out % real_vector(u)
  end subroutine step_partial_along
  subroutine step_partials(this, v)
    class(expansion)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: v(:,:)
    real(dp), allocatable :: weights(:), column(:)
    integer :: j
    call this % weights_of_steps(weights)
    allocate(v(size(weights) + 1, size(weights)))
    do j = 1, size(weights)
       call this % step_partial_along([j], column)
       v(:, j) = column
    end do
  end subroutine step_partials
  subroutine weights_of_steps(this, weights)
    class(expansion)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: weights(:)
    integer :: k
    do k = 1, this % num_designs()
       if (this % design_kind(k) == design_of_steps) then
          call this % design_value(k, weights)
          return
       end if
    end do
    error stop 'gti_expansion: the steps of this tower are not designs'
  end subroutine weights_of_steps
  integer function spatial_discretization_coupling(this, spatial_discretization_stencil) result(at)
    class(expansion), intent(inout) :: this
    type(stencil)   , intent(in)    :: spatial_discretization_stencil
    integer , allocatable :: table(:,:), heads(:), tails(:), rows(:), cols(:)
    real(dp), allocatable :: given(:), w(:)
    integer :: read_nodes, entered_nodes, holder, e, ne
    if (spatial_discretization_stencil % pattern % num_vertices() /= this % node_extent) then
       error stop 'gti_expansion: the spatial discretization stencil is a stencil over the nodes'
    end if
    ne = spatial_discretization_stencil % pattern % num_edges()
    heads = [(spatial_discretization_stencil % pattern % edge_head(e), e = 1, ne)]
    tails = [(spatial_discretization_stencil % pattern % edge_tail(e), e = 1, ne)]
    call spatial_discretization_stencil % weights % real_vector(given)
    call combine_triples(this % node_extent, this % node_extent, heads, tails, given, &
         & rows, cols, w)
    allocate(table(2, size(rows)))
    table(1, :) = cols
    table(2, :) = rows
    read_nodes    = named_set(this, this % node_extent, 'the nodes the spatial discretization stencil reads')
    entered_nodes = named_set(this, this % node_extent, 'the nodes whose rows the spatial discretization stencil enters')
    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the spatial discretization stencil''s reach')
    at = this % nodes % couple([read_nodes, entered_nodes], [holder])
    call bind_carriers(this, [read_nodes, entered_nodes])
    call bind_reach(this, holder, read_nodes, entered_nodes, table)
    call this % labels % bind(this % node(at), 'the spatial discretization stencil')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))
  end function spatial_discretization_coupling
  subroutine attach_known(this, at, x)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: at
    real(dp)        , intent(in)    :: x(:)
    call this % values % attach_unknown(this % node(at))
    if (size(x) == 0) then
       call this % values % mark_known(this % node(at), [0.0_dp])
    else
       call this % values % mark_known(this % node(at), x)
    end if
  end subroutine attach_known
  integer function one_sweep(this, physics, schemes, instants, dt, sensitivity) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    real(dp)              , intent(in)    :: dt(:)
    integer               , intent(in)    :: sensitivity
    at = this % nodes % assemble([one_horizon(this, physics, schemes, instants, dt)], 0)
    if (sensitivity == 0) then
       call this % labels % bind(this % node(at), 'sweep 0, the functional itself')
    else
       call this % labels % bind(this % node(at), 'sweep ' // written(sensitivity) // &
            & ', derivative ' // written(sensitivity) // ' in the design')
    end if
    call this % values % attach_unknown(this % node(at))
  end function one_sweep
  integer function one_horizon(this, physics, schemes, instants, dt) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    real(dp)              , intent(in)    :: dt(:)
    integer, allocatable :: blocks(:)
    integer :: b, first
    allocate(blocks(size(instants)))
    first = 1
    do b = 1, size(instants)
       blocks(b) = one_block(this, physics, schemes(b) % scheme, &
            & first, first + instants(b) - 1, dt)
       first = first + instants(b)
    end do
    at = this % nodes % assemble(blocks, 0)
    call this % labels % bind(this % node(at), 'horizon of duration ' // &
         & written(sum(dt)))
  end function one_horizon
  integer function one_block(this, physics, scheme, first, last, dt) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: first, last
    real(dp)              , intent(in)    :: dt(:)
    integer, allocatable :: slices(:)
    integer :: k, coupling
    if (last - first + 1 <= scheme % history_depth(this % degrees - 1)) then
       error stop 'gti_expansion: a block holds more instants than its family reaches'
    end if
    allocate(slices(last - first + 1))
    do k = first, last
       slices(k - first + 1) = one_slice(this, physics, scheme, k, first, dt(k))
    end do
    if (marches_by_stages(scheme, this % degrees)) then
       coupling = carry_coupling(this, scheme, slices, first, last, dt)
    else
       coupling = block_coupling(this, physics, scheme, slices, first, last, dt)
    end if
    at       = this % nodes % assemble(slices, coupling)
    call this % labels % bind(this % node(at), scheme % name() // ' block')
    call attach_known(this, at, dt(first:last))
  end function one_block
  integer function one_slice(this, physics, scheme, instant, first, step) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: instant, first
    real(dp)              , intent(in)    :: step
    associate (u1 => physics); end associate
    if (marches_by_stages(scheme, this % degrees)) then
       at = staged_slice(this, scheme, instant, first, step)
    else
       at = plain_slice(this, scheme, instant, first)
    end if
    call this % labels % bind(this % node(at), 'slice at instant ' // written(instant))
  end function one_slice
  integer function plain_slice(this, scheme, instant, first) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first
    integer, allocatable :: components(:)
    integer :: d
    allocate(components(this % degrees))
    do d = 0, this % degrees - 1
       components(d + 1) = one_component(this, scheme, instant, first, d, .true.)
    end do
    at = this % nodes % assemble(components, 0)
  end function plain_slice
  logical function marches_by_stages(scheme, nd) result(staged)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd
    integer, allocatable :: offset(:), degrees_of(:)
    integer :: d, primary
    primary = scheme % primary_degree(nd - 1)
    staged  = .true.
    do d = 0, nd - 1
       if (d == primary) cycle
       call scheme % row_pattern(d, nd - 1, offset, degrees_of)
       if (size(offset) > 0) then
          staged = .false.
          return
       end if
    end do
  end function marches_by_stages
  integer function staged_slice(this, scheme, instant, first, step) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first
    real(dp)        , intent(in)    :: step
    integer, allocatable :: members(:)
    integer :: s, i
    s = scheme % num_stages()
    if (instant == first) then
       at = this % nodes % assemble([stage_node(this, scheme, instant, first, 0)], 0)
       return
    end if
    allocate(members(s + 1))
    do i = 1, s
       members(i) = stage_node(this, scheme, instant, first, i)
    end do
    members(s + 1) = stage_node(this, scheme, instant, first, 0)
    at = this % nodes % assemble(members, &
         & slice_coupling(this, scheme, members, s, step))
  end function staged_slice
  integer function stage_node(this, scheme, instant, first, index) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first, index
    integer, allocatable :: components(:)
    integer :: d
    allocate(components(this % degrees))
    do d = 0, this % degrees - 1
       components(d + 1) = one_component(this, scheme, instant, first, d, index > 0)
    end do
    at = this % nodes % assemble(components, 0)
    if (index == 0) then
       call this % labels % bind(this % node(at), 'the arriving instant')
    else
       call this % labels % bind(this % node(at), 'stage ' // written(index))
    end if
  end function stage_node
  integer function one_component(this, scheme, instant, first, degree, evaluated) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first, degree
    logical         , intent(in)    :: evaluated
    if (evaluated .and. degree == scheme % primary_degree(this % degrees - 1) &
         & .and. this % spatial_coupling_at > 0) then
       at = this % nodes % assemble([integer ::], this % spatial_coupling_at)
    else
       at = this % nodes % assemble([integer ::], 0)
    end if
    call this % labels % bind(this % node(at), 'component of degree ' // written(degree))
    call this % extents % bind(this % node(at), counted_set_representation(this % node_extent))
    if (instant - first < scheme % history_depth(this % degrees - 1)) then
       call attach_known(this, at, [0.0_dp])
    else
       call this % values % attach_unknown(this % node(at))
    end if
  end function one_component
  function written(n) result(text)
    class(*), intent(in) :: n
    character(len=:), allocatable :: text
    character(len=24) :: buffer
    select type (n)
    type is (integer)
       write(buffer,'(i0)') n
    type is (real(dp))
       write(buffer,'(f0.4)') n
    class default
       error stop 'gti_expansion: a label is written from a number'
    end select
    text = trim(buffer)
  end function written
  subroutine block_reach(scheme, nd, n, tails, heads, source_degree, determines)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: nd, n
    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)
    integer, allocatable :: offset(:), degrees_of(:)
    integer :: primary, kk, d, e, counted, at, pass
    primary = scheme % primary_degree(nd - 1)
    do pass = 1, 2
       counted = 0
       do kk = 1, n
          do d = 0, nd - 1
             if (d == primary) cycle
             call scheme % row_pattern(d, nd - 1, offset, degrees_of)
             if (size(offset) == 0) cycle
             if (kk - maxval(offset) < 1) cycle
             do e = 1, size(offset)
                counted = counted + 1
                if (pass == 2) then
                   at = counted
                   tails(at)         = kk - offset(e)
                   heads(at)         = kk
                   source_degree(at) = degrees_of(e)
                   determines(at)    = d
                end if
             end do
          end do
       end do
       if (pass == 1) allocate(tails(counted), heads(counted), &
            & source_degree(counted), determines(counted))
    end do
  end subroutine block_reach
  subroutine reach_weights(scheme, n, tails, heads, source_degree, determines, dt, w)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: n, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: dt(:)
    real(dp), allocatable, intent(out) :: w(:)
    call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)
  end subroutine reach_weights
  integer function block_coupling(this, physics, scheme, slices, first, last, dt) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: slices(:), first, last
    real(dp)              , intent(in)    :: dt(:)
    integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: w(:)
    integer :: n, nd, components, constraints, holder
    associate (u1 => physics); end associate
    n  = last - first + 1
    nd = this % degrees
    call block_reach(scheme, nd, n, tails, heads, source_degree, determines)
    call reach_weights(scheme, n, tails, heads, source_degree, determines, &
         & dt(first:last), w)
    components  = named_set(this, n * nd, 'the components of this block')
    constraints = named_set(this, n * nd, 'the constraint instances of this block')
    table       = tuples(nd, tails, heads, source_degree, determines)
    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the scheme reach')
    at = this % nodes % couple([slices, components, constraints], [holder])
    call bind_carriers(this, [slices, components, constraints])
    call bind_reach(this, holder, components, constraints, table)
    call this % labels % bind(this % node(at), scheme % name() // ' coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))
  end function block_coupling
  pure function tuples(nd, tails, heads, source_degree, determines) result(table)
    integer, intent(in) :: nd, tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    integer :: e
    allocate(table(2, size(tails)))
    table(1,:) = [((tails(e) - 1) * nd + source_degree(e) + 1, e = 1, size(tails))]
    table(2,:) = [((heads(e) - 1) * nd + determines(e) + 1, e = 1, size(heads))]
  end function tuples
  integer function named_set(this, n, text) result(at)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: n
    character(len=*), intent(in)    :: text
    at = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(at), text)
    call this % extents % bind(this % node(at), counted_set_representation(n))
  end function named_set
  subroutine bind_carriers(this, carriers)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: carriers(:)
    type(graph), pointer :: g
    integer :: i
    do i = 1, size(carriers)
       g => this % nodes % node(carriers(i))
       call this % bindings % bind_set(g, g)
    end do
  end subroutine bind_carriers
  subroutine bind_reach(this, holder, components, constraints, table)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: holder, components, constraints, table(:,:)
    type(graph), pointer :: g, from, into
    from => this % nodes % node(components)
    into => this % nodes % node(constraints)
    g    => this % nodes % node(holder)
    call this % bindings % bind_relation(g, &
         & csr_relation('scheme reach', from, into, table, this % extents))
  end subroutine bind_reach
  recursive logical function consistent(this, g) result(ok)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    type(graph), pointer :: coupling
    ok = level_consistent(g)
    if (.not. ok) return
    if (level_couples(g)) then
       coupling => level_coupling(g)
       ok = relational_valid(coupling, this % bindings)
       if (.not. ok) return
    end if
    if (level_is_leaf(g)) return
    ok = every_member(this, level_members(g))
  end function consistent
  recursive logical function every_member(this, members) result(ok)
    class(expansion), intent(in) :: this
    type(branch)    , intent(in) :: members
    type(graph), pointer :: first
    ok = .true.
    if (sequence_empty(members)) return
    first => sequence_first(members)
    ok = this % consistent(first)
    if (ok) ok = every_member(this, sequence_rest(members))
  end function every_member
  subroutine stage_reach(nd, s, tails, heads, source_degree, determines)
    integer, intent(in) :: nd, s
    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)
    integer :: d, i, j, at
    at = (nd - 1) * (s * (s + 1) / 2 + s) + s
    allocate(tails(at), heads(at), source_degree(at), determines(at))
    at = 0
    do d = 0, nd - 1
       do i = 1, s
          if (d == nd - 1) cycle
          do j = 1, i
             at = at + 1
             tails(at) = 1 + j
             heads(at) = 1 + i
             source_degree(at) = d + 1
             determines(at) = d
          end do
       end do
       do j = 1, s
          at = at + 1
          tails(at) = 1 + j
          heads(at) = 2 + s
          source_degree(at) = min(d + 1, nd - 1)
          determines(at) = d
       end do
    end do
  end subroutine stage_reach
  pure integer function stage_unknown(vertex, degree, s, nd) result(at)
    integer, intent(in) :: vertex, degree, s, nd
    integer :: member
    if (vertex == 2 + s) then
       member = s + 1
    else
       member = vertex - 1
    end if
    at = (member - 1) * nd + degree + 1
  end function stage_unknown
  subroutine stage_weights(scheme, s, tails, heads, source_degree, determines, step, w)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: s, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: step
    real(dp), allocatable, intent(out) :: w(:)
    call weights_of(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
         & source_degree, determines, w)
  end subroutine stage_weights
  integer function slice_coupling(this, scheme, members, s, step) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: members(:), s
    real(dp)        , intent(in)    :: step
    integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: w(:)
    integer :: nd, components, constraints, holder, e
    nd = this % degrees
    call stage_reach(nd, s, tails, heads, source_degree, determines)
    call stage_weights(scheme, s, tails, heads, source_degree, determines, step, w)
    components  = named_set(this, (s + 1) * nd, 'the components of this step')
    constraints = named_set(this, (s + 1) * nd, 'the constraint instances of this step')
    allocate(table(2, size(tails)))
    table(1,:) = [(stage_unknown(tails(e), source_degree(e), s, nd), e = 1, size(tails))]
    table(2,:) = [(stage_unknown(heads(e), determines(e), s, nd), e = 1, size(heads))]
    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the butcher reach')
    at = this % nodes % couple([members, components, constraints], [holder])
    call bind_carriers(this, [members, components, constraints])
    call bind_reach(this, holder, components, constraints, table)
    call this % labels % bind(this % node(at), scheme % name() // ' stage coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))
  end function slice_coupling
  pure integer function slice_base(kk, s, nd) result(at)
    integer, intent(in) :: kk, s, nd
    if (kk == 1) then
       at = 0
    else
       at = (1 + (kk - 2) * (s + 1)) * nd
    end if
  end function slice_base
  pure integer function closing_instant(kk, s, nd) result(at)
    integer, intent(in) :: kk, s, nd
    if (kk == 1) then
       at = slice_base(kk, s, nd)
    else
       at = slice_base(kk, s, nd) + s * nd
    end if
  end function closing_instant
  subroutine carry_reach(this, s, n, table, sources)
    class(expansion), intent(in) :: this
    integer         , intent(in) :: s, n
    integer, allocatable, intent(out) :: table(:,:)
    integer, allocatable, intent(out) :: sources(:)
    integer :: nd, kk, d, m, counted, pass, from, into
    nd = this % degrees
    do pass = 1, 2
       counted = 0
       do kk = 2, n
          from = closing_instant(kk - 1, s, nd)
          do d = 0, nd - 2
             do m = 1, s + 1
                counted = counted + 1
                into = slice_base(kk, s, nd) + (m - 1) * nd + d + 1
                if (pass == 2) then
                   table(1, counted) = from + d + 1
                   table(2, counted) = into
                   sources(counted)  = merge(2 + s, 1 + m, m == s + 1)
                end if
             end do
          end do
       end do
       if (pass == 1) then
          allocate(table(2, counted), sources(counted))
       end if
    end do
  end subroutine carry_reach
  integer function carry_coupling(this, scheme, slices, first, last, dt) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: slices(:), first, last
    real(dp)        , intent(in)    :: dt(:)
    integer, allocatable :: table(:,:), sources(:)
    real(dp), allocatable :: w(:)
    integer :: n, nd, s, unknowns, components, constraints, holder, e
    n  = last - first + 1
    nd = this % degrees
    s  = scheme % num_stages()
    call carry_reach(this, s, n, table, sources)
    unknowns = (1 + (n - 1) * (s + 1)) * nd
    call stage_weights(scheme, s, [(1, e = 1, size(sources))], sources, &
         & [((mod(table(1, e) - 1, nd)), e = 1, size(sources))], &
         & [((mod(table(2, e) - 1, nd)), e = 1, size(sources))], dt(first), w)
    components  = named_set(this, unknowns, 'the components of this block')
    constraints = named_set(this, unknowns, 'the constraint instances of this block')
    holder = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(holder), 'the carry between steps')
    at = this % nodes % couple([slices, components, constraints], [holder])
    call bind_carriers(this, [slices, components, constraints])
    call bind_reach(this, holder, components, constraints, table)
    call this % labels % bind(this % node(at), scheme % name() // ' carry coupling')
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))
  end function carry_coupling
end module gti_expansion
module gti_block
  use util_precision  , only : dp
  use operation_action     , only : operation, variation
  use view_directed        , only : directed_graph
  use view_directed_stored , only : stored_directed_graph
  use field_calculus       , only : field
  use field_stored         , only : stored_field
  use graph_fractal        , only : graph
  use gti_expansion        , only : expansion
  use view_level           , only : level_member, level_num_members, level_is_leaf
  use operation_stencil    , only : combine_triples, stencil
  use operation_family     , only : family
  use operation_weight     , only : scheme_weight
  use operation_coupling   , only : weights_terms
  use operation_expression    , only : expression, constant, stated
  use view_directed        , only : forward
  implicit none
  private
  public :: block_residual, coupling_reach
  type :: coupling_reach
     integer :: vertices = 0
     integer, allocatable :: step_of(:)
     integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
     integer, allocatable :: row(:), column(:)
  end type coupling_reach
  type, extends(operation) :: block_residual
     type(stencil)                       , private :: time_discretization_stencil
     type(expression)     , private :: physics
     type(stencil), allocatable, private :: spatial_discretization_stencil
     type(stored_directed_graph)         , private :: points
     integer , allocatable               , private :: at(:)
     integer , allocatable               , private :: carried(:)
     real(dp), allocatable               , private :: held(:)
     integer                             , private :: degrees  = 0
     integer                             , private :: unknowns = 0
     integer                             , private :: primary  = 0
     type(graph)    , pointer, private :: node  => null()
     type(expansion), pointer, private :: tower => null()
     integer, allocatable    , private :: kept(:)
     type(coupling_reach), allocatable, private :: reach(:)
   contains
     procedure :: name           => block_name
     procedure :: domain         => block_domain
     procedure :: apply          => block_apply
     procedure :: max_degree     => block_max_degree
     procedure :: partial_action => block_partial_action
     procedure :: compiled_tangent => block_compiled_tangent
     procedure :: restricted => block_restricted
     procedure :: placed_on
     procedure :: slice_of
     procedure :: node_of
     procedure :: moment_of
     procedure, private :: labels_of
     procedure :: num_nodes
     procedure :: spatial_discretization_laid
     procedure :: aggregates
     procedure :: with_reach
     procedure :: rows_terms
     procedure :: linear_block
     procedure :: member_order
     procedure :: num_unknowns
     procedure :: num_degrees
     procedure :: num_points
     procedure :: points_at
     procedure :: num_carried
     procedure :: carried_unknowns
     procedure :: held_values
     procedure :: first_held
  end type block_residual
  interface block_residual
     module procedure create
  end interface block_residual
contains
  function create(derived, physics, at, unknowns, degrees, primary, carried, held, &
       & spatial_discretization_stencil) result(this)
    type(stencil)         , intent(in) :: derived
    type(expression)      , intent(in) :: physics
    integer               , intent(in) :: at(:), unknowns, degrees, primary
    integer               , intent(in) :: carried(:)
    real(dp)              , intent(in) :: held(:)
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
    type(block_residual) :: this
    if (size(carried) /= size(held)) then
       error stop 'gti_block: one value per carried component'
    end if
    if (any(carried < 1) .or. any(carried > unknowns)) then
       error stop 'gti_block: every carried row names an unknown'
    end if
    if (any(at < 0) .or. any(at + degrees > unknowns)) then
       error stop 'gti_block: an evaluation point holds its degrees within the unknowns'
    end if
    this % time_discretization_stencil  = derived
    if (present(spatial_discretization_stencil)) this % spatial_discretization_stencil = spatial_discretization_stencil
    this % physics = physics
    this % at       = at
    this % unknowns = unknowns
    this % degrees  = degrees
    this % primary  = primary
    this % carried  = carried
    this % held     = held
    this % points = stored_directed_graph(size(at), tails=[integer ::], heads=[integer ::])
    call this % declare_arguments(2)
  end function create
  pure integer function num_unknowns(this)
    class(block_residual), intent(in) :: this
    num_unknowns = this % unknowns
  end function num_unknowns
  pure function carried_unknowns(this) result(c)
    class(block_residual), intent(in) :: this
    integer, allocatable :: c(:)
    c = this % carried
  end function carried_unknowns
  pure function held_values(this) result(h)
    class(block_residual), intent(in) :: this
    real(dp), allocatable :: h(:)
    h = this % held
  end function held_values
  pure integer function num_degrees(this)
    class(block_residual), intent(in) :: this
    num_degrees = this % degrees
  end function num_degrees
  pure integer function num_points(this)
    class(block_residual), intent(in) :: this
    num_points = size(this % at)
  end function num_points
  pure function points_at(this) result(at)
    class(block_residual), intent(in) :: this
    integer, allocatable :: at(:)
    at = this % at
  end function points_at
  pure function first_held(this) result(x)
    class(block_residual), intent(in) :: this
    real(dp), allocatable :: x(:)
    allocate(x(this % degrees), source=0.0_dp)
    if (size(this % held) >= this % degrees) x = this % held(1:this % degrees)
  end function first_held
  pure integer function num_carried(this)
    class(block_residual), intent(in) :: this
    num_carried = size(this % carried)
  end function num_carried
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
  pure integer function block_max_degree(this)
    class(block_residual), intent(in) :: this
    associate (u1 => this); end associate
    block_max_degree = 3
  end function block_max_degree
  pure function gathered(this, x) result(y)
    class(block_residual), intent(in) :: this
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: p
    allocate(y(size(this % at) * this % degrees))
    do p = 1, size(this % at)
       y((p - 1) * this % degrees + 1:p * this % degrees) = &
            & x(this % at(p) + 1:this % at(p) + this % degrees)
    end do
  end function gathered
  subroutine point_inputs(this, input_data, x, inputs)
    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    real(dp)             , intent(in) :: x(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)
    type(stored_field) :: state, design
    real(dp), allocatable :: knob(:)
    call input_data(2) % real_vector(knob)
    state = stored_field('state', this % points % vertex_set(), &
         & size(this % at) * this % degrees)
    call state % set_real_vector(gathered(this, x))
    design = stored_field('design', this % points % vertex_set(), size(knob))
    call design % set_real_vector(knob)
    inputs = [state, design]
  end subroutine point_inputs
  pure subroutine placed(this, governing, r)
    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: governing(:)
    real(dp)             , intent(inout) :: r(:)
    integer :: p
    do p = 1, size(this % at)
       r(this % at(p) + this % primary + 1) = &
            & r(this % at(p) + this % primary + 1) + governing(p)
    end do
  end subroutine placed
  pure subroutine carry(this, x, r)
    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: x(:)
    real(dp)             , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % carried)
       r(this % carried(i)) = x(this % carried(i)) - this % held(i)
    end do
  end subroutine carry
  pure subroutine carry_direction(this, v, r)
    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: v(:)
    real(dp)             , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % carried)
       r(this % carried(i)) = v(this % carried(i))
    end do
  end subroutine carry_direction
  pure subroutine carry_held(this, r)
    class(block_residual), intent(in)    :: this
    real(dp)             , intent(inout) :: r(:)
    integer :: i
    do i = 1, size(this % carried)
       r(this % carried(i)) = 0.0_dp
    end do
  end subroutine carry_held
  subroutine state_of(this, input_data, input_graph, x, state)
    class(block_residual), intent(in)  :: this
    class(field)         , intent(in)  :: input_data(:)
    class(directed_graph), intent(in)  :: input_graph
    real(dp), allocatable, intent(out) :: x(:)
    type(stored_field)   , intent(out) :: state
    if (size(input_data) < 2) then
       error stop 'gti_block: the state and the design are given'
    end if
    call input_data(1) % real_vector(x)
    if (size(x) /= this % num_unknowns()) then
       error stop 'gti_block: the state holds one component per degree per unknown point'
    end if
    state = stored_field('state', input_graph % vertex_set(), size(x))
    call state % set_real_vector(x)
  end subroutine state_of
  subroutine placed_output(this, input_graph, r, output)
    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    real(dp)             , intent(in) :: r(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: out
    out = stored_field(this % name(), input_graph % vertex_set(), size(r))
    call out % set_real_vector(r)
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)
  end subroutine placed_output
  subroutine block_apply(this, input_graph, input_data, output)
    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: state
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half
    real(dp), allocatable :: r(:), governing(:), x(:), coupled(:)
    if (.not. present(input_data)) then
       error stop 'gti_block: the state and the design are given'
    end if
    call state_of(this, input_data, input_graph, x, state)
    call point_inputs(this, input_data, x, inputs)
    call this % time_discretization_stencil % apply(input_graph, [state], half)
    call half % real_vector(r)
    if (allocated(this % spatial_discretization_stencil)) then
       call this % spatial_discretization_stencil % apply(input_graph, [state], half)
       call half % real_vector(coupled)
       r = r + coupled
    end if
    call this % physics % apply(this % points, inputs, half)
    call half % real_vector(governing)
    call placed(this, governing, r)
    call carry(this, x, r)
    call placed_output(this, input_graph, r, output)
  end subroutine block_apply
  subroutine placed_on(this, tower, node)
    class(block_residual), intent(inout)     :: this
    type(expansion)      , intent(in), target :: tower
    type(graph)          , intent(in), target :: node
    this % tower => tower
    this % node  => node
  end subroutine placed_on
  subroutine labels_of(this, slice, node, moment)
    class(block_residual), intent(in) :: this
    integer, allocatable , intent(out) :: slice(:), node(:), moment(:)
    type(graph), pointer :: one_slice, first
    integer, allocatable :: whole_slice(:), whole_node(:), whole_moment(:)
    integer :: n, k, j, members, moments, g, m, count, u, i, d
    if (.not. associated(this % node)) then
       error stop 'gti_block: the block has not been placed in the graph'
    end if
    n = level_num_members(this % node)
    moments = 0
    do k = 1, n
       moments = moments + members_of(level_member(this % node, k))
    end do
    first => level_member(this % node, 1)
    if (.not. level_is_leaf(level_member(first, 1))) first => level_member(first, 1)
    m     = this % tower % extent_of(level_member(first, 1))
    count = moments * m * this % degrees
    allocate(whole_slice(count), whole_node(count), whole_moment(count))
    g = 0
    do k = 1, n
       one_slice => level_member(this % node, k)
       members = members_of(one_slice)
       do j = 1, members
          g = g + 1
          do i = 1, m
             do d = 0, this % degrees - 1
                u = ((g - 1) * m + (i - 1)) * this % degrees + d + 1
                whole_slice(u)  = k
                whole_node(u)   = i
                whole_moment(u) = g
             end do
          end do
       end do
    end do
    if (allocated(this % kept)) then
       slice  = whole_slice(this % kept)
       node   = whole_node(this % kept)
       moment = whole_moment(this % kept)
    else
       call move_alloc(whole_slice , slice)
       call move_alloc(whole_node  , node)
       call move_alloc(whole_moment, moment)
    end if
  contains
    integer function members_of(one_slice)
      type(graph), intent(in) :: one_slice
      if (level_is_leaf(level_member(one_slice, 1))) then
         members_of = 1
      else
         members_of = level_num_members(one_slice)
      end if
    end function members_of
  end subroutine labels_of
  function slice_of(this) result(slice)
    class(block_residual), intent(in) :: this
    integer, allocatable :: slice(:)
    integer, allocatable :: node(:), moment(:)
    call this % labels_of(slice, node, moment)
  end function slice_of
  function node_of(this) result(node)
    class(block_residual), intent(in) :: this
    integer, allocatable :: node(:)
    integer, allocatable :: slice(:), moment(:)
    call this % labels_of(slice, node, moment)
  end function node_of
  function moment_of(this) result(moment)
    class(block_residual), intent(in) :: this
    integer, allocatable :: moment(:)
    integer, allocatable :: slice(:), node(:)
    call this % labels_of(slice, node, moment)
  end function moment_of
  integer function num_nodes(this)
    class(block_residual), intent(in) :: this
    integer, allocatable :: slice(:), node(:), moment(:)
    num_nodes = 1
    if (.not. associated(this % node)) return
    call this % labels_of(slice, node, moment)
    num_nodes = maxval(node)
  end function num_nodes
  subroutine spatial_discretization_laid(this, spatial_discretization_stencil)
    class(block_residual), intent(inout) :: this
    type(stencil)        , intent(in)    :: spatial_discretization_stencil
    integer , allocatable :: base(:,:), r(:), c(:), slice(:), node(:), moment(:)
    real(dp), allocatable :: lw(:), held(:), w(:)
    integer :: nodes, moments, p, u, e, g, ne, n, rc, cc
    call this % labels_of(slice, node, moment)
    nodes   = maxval(node)
    moments = maxval(moment)
    if (spatial_discretization_stencil % pattern % num_vertices() /= nodes) then
       error stop 'gti_block: the spatial discretization stencil is a stencil over the nodes'
    end if
    call spatial_discretization_stencil % constants % real_vector(held)
    if (any(abs(held) > 0.0_dp)) then
       error stop 'gti_block: the spatial discretization stencil carries no constant'
    end if
    call spatial_discretization_stencil % weights % real_vector(lw)
    allocate(base(nodes, moments), source=-1)
    do p = 1, size(this % at)
       u = this % at(p) + 1
       base(node(u), moment(u)) = this % at(p)
    end do
    ne = spatial_discretization_stencil % pattern % num_edges()
    allocate(r(ne * moments), c(ne * moments), w(ne * moments))
    n = 0
    do g = 1, moments
       do e = 1, ne
          rc = spatial_discretization_stencil % pattern % edge_head(e)
          cc = spatial_discretization_stencil % pattern % edge_tail(e)
          if (base(rc, g) < 0) cycle
          if (base(cc, g) < 0) then
             error stop 'gti_block: a moment holds every node or none'
          end if
          n    = n + 1
          r(n) = base(rc, g) + this % primary + 1
          c(n) = base(cc, g) + 1
          w(n) = lw(e)
       end do
    end do
    this % spatial_discretization_stencil = stencil(r(1:n), c(1:n), w(1:n), spread(0.0_dp, 1, this % unknowns), &
         & 'spatial discretization stencil')
  end subroutine spatial_discretization_laid
  subroutine with_reach(this, reach)
    class(block_residual), intent(inout) :: this
    type(coupling_reach) , intent(in)    :: reach(:)
    this % reach = reach
  end subroutine with_reach
  subroutine rows_terms(this, scheme, dt, seeds, r, c, w)
    class(block_residual), intent(in) :: this
    class(family)        , intent(in) :: scheme
    real(dp)             , intent(in) :: dt(:), seeds(:,:)
    integer , allocatable, intent(out) :: r(:), c(:)
    real(dp), allocatable, intent(out) :: w(:,:)
    real(dp), allocatable :: table(:,:)
    integer :: k, e, i, nodes, count, n
    if (.not. allocated(this % reach)) then
       error stop 'gti_block: the block was built without its reach'
    end if
    nodes = this % num_nodes()
    count = 0
    do k = 1, size(this % reach)
       count = count + size(this % reach(k) % tails) * nodes
    end do
    allocate(r(count), c(count), w(count, 0:size(seeds, 2)))
    n = 0
    do k = 1, size(this % reach)
       associate (reach => this % reach(k))
         call weights_terms(scheme_weight(scheme), reach % vertices, reach % tails, &
              & reach % heads, dt(reach % step_of), seeds(reach % step_of, :), &
              & reach % source_degree, reach % determines, table)
         do i = 1, nodes
            do e = 1, size(reach % tails)
               n       = n + 1
               r(n)    = reach % row(e)    + (i - 1) * this % degrees
               c(n)    = reach % column(e) + (i - 1) * this % degrees
               w(n, :) = -table(e, :)
            end do
         end do
       end associate
    end do
  end subroutine rows_terms
  function aggregates(this, cell) result(aggregate)
    class(block_residual), intent(in) :: this
    integer              , intent(in) :: cell(:)
    integer, allocatable :: aggregate(:)
    integer, allocatable :: numbered(:), slice(:), node(:), moment(:)
    integer :: u, coarse, key, count
    call this % labels_of(slice, node, moment)
    if (size(cell) < maxval(node)) then
       error stop 'gti_block: a coarse cell for every node'
    end if
    coarse = maxval(cell)
    allocate(aggregate(this % unknowns))
    allocate(numbered(maxval(moment) * coarse * this % degrees), source=0)
    count = 0
    do u = 1, this % unknowns
       key = ((moment(u) - 1) * coarse + cell(node(u)) - 1) * this % degrees &
            & + mod(u - 1, this % degrees) + 1
       if (numbered(key) == 0) then
          count         = count + 1
          numbered(key) = count
       end if
       aggregate(u) = numbered(key)
    end do
  end function aggregates
  function linear_block(this, input_graph, input_data, rhs, transposed, mark) result(lin)
    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field)         , intent(in) :: input_data(:)
    real(dp)             , intent(in) :: rhs(:)
    logical              , intent(in) :: transposed
    integer              , intent(in) :: mark
    type(block_residual) :: lin
    type(stencil) :: a
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical :: available
    if (size(rhs) /= this % unknowns) then
       error stop 'gti_block: one right side per unknown'
    end if
    call this % compiled_tangent(input_graph, input_data, 1, r, c, w, available)
    if (.not. available) then
       error stop 'gti_block: the tangent in the state compiles'
    end if
    a = stencil(r, c, w, spread(0.0_dp, 1, this % unknowns), 'frozen tangent')
    if (transposed) a = a % transpose()
    call a % constants % set_real_vector(-rhs)
    lin = block_residual(a, stated(constant(0.0_dp), this % degrees - 1, 'zero'), this % at, this % unknowns, &
         & this % degrees, this % primary, [integer ::], [real(dp) ::])
    call lin % stamped(mark, transposed=a % pattern % transposed())
    lin % tower => this % tower
    lin % node  => this % node
    if (allocated(this % kept)) lin % kept = this % kept
  end function linear_block
  function block_restricted(this, kept, values) result(sub)
    class(block_residual), intent(in) :: this
    integer              , intent(in) :: kept(:)
    real(dp)             , intent(in) :: values(:)
    type(block_residual) :: sub
    type(stencil) :: derived, spatial_discretization_stencil
    integer , allocatable :: sub_of(:), at(:), carried(:)
    real(dp), allocatable :: held(:)
    integer :: e, p, d, inside, npts, ncar
    allocate(sub_of(this % unknowns), source=0)
    do e = 1, size(kept)
       sub_of(kept(e)) = e
    end do
    npts = 0
    allocate(at(size(this % at)))
    do p = 1, size(this % at)
       inside = 0
       do d = 1, this % degrees
          if (sub_of(this % at(p) + d) > 0) inside = inside + 1
       end do
       if (inside == 0) cycle
       if (inside /= this % degrees) then
          error stop 'gti_block: a member holds whole points'
       end if
       npts     = npts + 1
       at(npts) = sub_of(this % at(p) + 1) - 1
    end do
    ncar = 0
    allocate(carried(size(this % carried)), held(size(this % carried)))
    do e = 1, size(this % carried)
       if (sub_of(this % carried(e)) == 0) cycle
       ncar          = ncar + 1
       carried(ncar) = sub_of(this % carried(e))
       held(ncar)    = this % held(e)
    end do
    derived = this % time_discretization_stencil % restricted(kept, values)
    if (allocated(this % spatial_discretization_stencil)) then
       spatial_discretization_stencil = this % spatial_discretization_stencil % restricted(kept, values)
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar), spatial_discretization_stencil=spatial_discretization_stencil)
    else
       sub = block_residual(derived, this % physics, at(1:npts), size(kept), &
            & this % degrees, this % primary, carried(1:ncar), held(1:ncar))
    end if
    sub % tower => this % tower
    sub % node  => this % node
    if (allocated(this % kept)) then
       sub % kept = this % kept(kept)
    else
       sub % kept = kept
    end if
  end function block_restricted
  subroutine block_compiled_tangent(this, input_graph, input_data, which, &
       & rows, columns, weights, available)
    class(block_residual), intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    class(field)         , intent(in)  :: input_data(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: available
    type(stored_field) :: state, direction
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: out
    real(dp), allocatable :: x(:), w(:), governing(:), v(:)
    integer , allocatable :: r(:), c(:)
    logical , allocatable :: is_carried(:)
    integer :: e, d, p, npts, n, kept, count
    available = which == 1
    if (.not. available) return
    n    = this % unknowns
    npts = size(this % at)
    allocate(is_carried(n), source=.false.)
    is_carried(this % carried) = .true.
    call state_of(this, input_data, input_graph, x, state)
    call point_inputs(this, input_data, x, inputs)
    count = this % time_discretization_stencil % pattern % num_edges() + npts * this % degrees + size(this % carried)
    if (allocated(this % spatial_discretization_stencil)) count = count + this % spatial_discretization_stencil % pattern % num_edges()
    allocate(r(count), c(count), w(count))
    kept = 0
    call stencil_triples(this % time_discretization_stencil, is_carried, r, c, w, kept)
    if (allocated(this % spatial_discretization_stencil)) call stencil_triples(this % spatial_discretization_stencil, is_carried, r, c, w, kept)
    allocate(v(npts * this % degrees))
    do d = 0, this % degrees - 1
       v = 0.0_dp
       do p = 1, npts
          v((p - 1) * this % degrees + d + 1) = 1.0_dp
       end do
       direction = stored_field('direction', this % points % vertex_set(), size(v))
       call direction % set_real_vector(v)
       call this % physics % partial_action(this % points, inputs, &
            & [variation(this % physics % argument(1), direction)], out)
       call out % real_vector(governing)
       do p = 1, npts
          if (is_carried(this % at(p) + this % primary + 1)) cycle
          kept    = kept + 1
          r(kept) = this % at(p) + this % primary + 1
          c(kept) = this % at(p) + d + 1
          w(kept) = governing(p)
       end do
    end do
    do e = 1, size(this % carried)
       kept    = kept + 1
       r(kept) = this % carried(e)
       c(kept) = this % carried(e)
       w(kept) = 1.0_dp
    end do
    call combine_triples(n, n, r(1:kept), c(1:kept), w(1:kept), rows, columns, weights)
  end subroutine block_compiled_tangent
  subroutine stencil_triples(op, is_carried, r, c, w, kept)
    type(stencil), intent(in)    :: op
    logical      , intent(in)    :: is_carried(:)
    integer      , intent(inout) :: r(:), c(:)
    real(dp)     , intent(inout) :: w(:)
    integer      , intent(inout) :: kept
    real(dp), allocatable :: weights(:)
    integer :: e, row
    call op % weights % real_vector(weights)
    do e = 1, op % pattern % num_edges()
       row = op % pattern % edge_head(e)
       if (is_carried(row)) cycle
       kept    = kept + 1
       r(kept) = row
       c(kept) = op % pattern % edge_tail(e)
       w(kept) = weights(e)
    end do
  end subroutine stencil_triples
  subroutine block_partial_action(this, input_graph, input_data, variations, output)
    class(block_residual), intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field)         , intent(in)        :: input_data(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output
    type(stored_field) :: state
    real(dp), allocatable :: r(:), governing(:), v(:), x(:)
    call this % require_owned(variations)
    if (size(variations) < 1 .or. size(variations) > this % max_degree()) then
       error stop 'gti_block: the requested order is within max_degree'
    end if
    call state_of(this, input_data, input_graph, x, state)
    if (size(variations) >= 2) then
       call second_tangent(this, input_data, variations, x, governing)
       allocate(r(this % num_unknowns()), source=0.0_dp)
       call placed(this, governing, r)
       call carry_held(this, r)
       call placed_output(this, input_graph, r, output)
       return
    end if
    call variations(1) % direction(v)
    if (variations(1) % argument_is(this % argument(1))) then
       call state_tangent(this, input_graph, input_data, variations, state, x, v, &
            & r, governing)
       call placed(this, governing, r)
       call carry_direction(this, v, r)
    else if (variations(1) % argument_is(this % argument(2))) then
       call design_tangent(this, input_data, variations, x, r, governing)
       call placed(this, governing, r)
       call carry_held(this, r)
    else
       error stop 'gti_block: a variation names the state or the design'
    end if
    call placed_output(this, input_graph, r, output)
  end subroutine block_partial_action
  subroutine state_tangent(this, input_graph, input_data, variations, state, x, v, &
       & r, governing)
    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    type(stored_field)   , intent(in) :: state
    real(dp)             , intent(in) :: x(:), v(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:)
    type(stored_field) :: direction
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half
    real(dp), allocatable :: coupled(:)
    call this % time_discretization_stencil % partial_action(input_graph, [state], &
         & [variations(1) % with_argument(this % time_discretization_stencil % argument(1))], half)
    call half % real_vector(r)
    if (allocated(this % spatial_discretization_stencil)) then
       call this % spatial_discretization_stencil % partial_action(input_graph, [state], &
            & [variations(1) % with_argument(this % spatial_discretization_stencil % argument(1))], half)
       call half % real_vector(coupled)
       r = r + coupled
    end if
    call point_inputs(this, input_data, x, inputs)
    direction = stored_field('direction', this % points % vertex_set(), &
         & size(this % at) * this % degrees)
    call direction % set_real_vector(gathered(this, v))
    call this % physics % partial_action(this % points, inputs, &
         & [variation(this % physics % argument(1), direction)], half)
    call half % real_vector(governing)
  end subroutine state_tangent
  subroutine second_tangent(this, input_data, variations, x, governing)
    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: governing(:)
    type(stored_field), allocatable :: inputs(:)
    type(variation), allocatable :: at_points(:)
    class(field), allocatable :: half
    integer :: i
    call point_inputs(this, input_data, x, inputs)
    allocate(at_points(size(variations)))
    do i = 1, size(variations)
       at_points(i) = physics_variation(this, variations(i))
    end do
    call this % physics % partial_action(this % points, inputs, at_points, half)
    call half % real_vector(governing)
  end subroutine second_tangent
  function physics_variation(this, given) result(at_points)
    class(block_residual), intent(in) :: this
    type(variation)      , intent(in) :: given
    type(variation) :: at_points
    type(stored_field) :: direction
    real(dp), allocatable :: v(:)
    call given % direction(v)
    if (given % argument_is(this % argument(1))) then
       direction = stored_field('direction', this % points % vertex_set(), &
            & size(this % at) * this % degrees)
       call direction % set_real_vector(gathered(this, v))
       at_points = variation(this % physics % argument(1), direction)
    else if (given % argument_is(this % argument(2))) then
       at_points = given % with_argument(this % physics % argument(2))
    else
       error stop 'gti_block: a variation names the state or the design'
    end if
  end function physics_variation
  subroutine design_tangent(this, input_data, variations, x, r, governing)
    class(block_residual), intent(in) :: this
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    real(dp)             , intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: r(:), governing(:)
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: half
    allocate(r(this % num_unknowns()), source=0.0_dp)
    call point_inputs(this, input_data, x, inputs)
    call this % physics % partial_action(this % points, inputs, &
         & [variations(1) % with_argument(this % physics % argument(2))], half)
    call half % real_vector(governing)
  end subroutine design_tangent
  subroutine member_order(this, by_instants, order)
    class(block_residual), intent(in)  :: this
    logical              , intent(in)  :: by_instants
    integer, allocatable , intent(out) :: order(:)
    integer, allocatable :: table(:,:), label(:), slice(:), node(:), moment(:)
    type(stored_directed_graph) :: coupling
    integer :: ne, e, n, t, h, k, members
    call this % labels_of(slice, node, moment)
    if (.not. by_instants) then
       order = [(k, k = 1, maxval(node))]
       return
    end if
    label   = slice
    members = maxval(label)
    ne      = this % time_discretization_stencil % pattern % num_edges()
    allocate(table(2, ne))
    n = 0
    do e = 1, ne
       t = label(this % time_discretization_stencil % pattern % edge_tail(e))
       h = label(this % time_discretization_stencil % pattern % edge_head(e))
       if (t == h) cycle
       n = n + 1
       table(:, n) = [t, h]
    end do
    coupling = stored_directed_graph(members, tails=table(1, 1:n), heads=table(2, 1:n))
    order    = coupling % loop(forward)
  end subroutine member_order
end module gti_block
module gti_march
  use util_precision  , only : dp, least_kind_for
  use iso_fortran_env , only : real128
  use operation_coupling      , only : weights_of
  use gti_configuration       , only : refuse_unknown
  use operation_weight        , only : scheme_weight
  use view_directed_stored    , only : stored_directed_graph
  use view_directed           , only : directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field
  use operation_action      , only : variation, jacobian_of
  use operation_stencil       , only : stencil
  use operation_newton        , only : newton
  use operation_minimization  , only : minimizer, relative, absolute, &
       & by_count, by_rate
  use operation_dense_direct  , only : dense_direct
  use operation_gmres         , only : gmres
  use operation_family        , only : family
  use operation_grid          , only : grid, partition, partitioned
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use operation_expression       , only : expression
  use gti_expansion           , only : family_holder, expansion, marches_by_stages
  use gti_block               , only : block_residual, coupling_reach
  use view_level              , only : level_member, level_num_members, level_coupling, &
       & level_couples
  use graph_fractal           , only : graph
  use map_value               , only : VALUE_KNOWN
  use gti_sweeps              , only : assembly_present, multigrid_on, newton_order, &
       & set_aggregates, coarse_nodes, take_inner, keep_inner, forget_inner, set_linear_stopping
  use util_tally              , only : tally_record, tangent_loops, adjoint_loops
  implicit none
  type :: imbalance
     logical  :: converged = .true.
     logical  :: diverging = .false.
     real(dp) :: norm      = 0.0_dp
     real(dp) :: began     = 0.0_dp
     real(dp), allocatable :: by_degree(:)
     integer  :: worst_slot = 0, worst_degree = 0
     integer  :: steepest_slot = 0, steepest_degree = 0
     real(dp) :: steepest = 0.0_dp
  end type imbalance
  character(len=16), save :: sweep_level = 'space-time'
  integer, save :: stamps_given = 0
  real(dp), save :: stopping_tolerance  = 1.0e-12_dp
  integer , save :: stopping_criterion  = relative
  integer , save :: stopping_budget     = by_rate
  integer , save :: stopping_iterations = 100
  private
  public :: solved, unknowns_graph
  public :: block_from
  public :: unknown, consistent_states, frozen_inputs
  public :: set_stopping
  public :: consistent_state
  public :: imbalance
  public :: swept, set_sweep, sweep_named
  public :: solved_linear, by_tangent, by_adjoint, fresh_stamp
  public :: weight_of, precision_needed
  public :: horizon_bounds
contains
  real(dp) function weight_of(scheme, degrees, step) result(w)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees
    real(dp)     , intent(in) :: step
    integer, allocatable :: offset(:), source_degree(:)
    integer :: d, reach, s, i, k
    logical :: any_pattern
    w = 1.0_dp
    any_pattern = .false.
    do d = 0, degrees - 1
       call scheme % row_pattern(d, degrees - 1, offset, source_degree)
       if (size(offset) == 0) cycle
       any_pattern = .true.
       reach = maxval(offset)
       w = max(w, 1.0_dp + row_weight(scheme, reach + 1, &
            & [(reach + 1 - offset(k), k = 1, size(offset))], reach + 1, &
            & source_degree, d, step))
    end do
    if (any_pattern) return
    s = scheme % num_stages()
    do d = 0, degrees - 2
       do i = 1, s
          w = max(w, 1.0_dp + row_weight(scheme, s + 2, &
               & [1, (1 + k, k = 1, i)], 1 + i, &
               & [d, (d + 1, k = 1, i)], d, step))
       end do
       w = max(w, 1.0_dp + row_weight(scheme, s + 2, &
            & [1, (1 + k, k = 1, s)], s + 2, &
            & [d, (d + 1, k = 1, s)], d, step))
    end do
  end function weight_of
  real(dp) function row_weight(scheme, num_vertices, tails, head, source_degree, &
       & determines, step) result(total)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: num_vertices, tails(:), head, source_degree(:), determines
    real(dp)     , intent(in) :: step
    real(dp), allocatable :: c(:)
    integer :: k
    call weights_of(scheme_weight(scheme), num_vertices, tails, [(head, k = 1, size(tails))], &
         & [(step, k = 1, num_vertices)], source_degree, [(determines, k = 1, size(tails))], c)
    total = sum(abs(c))
  end function row_weight
  subroutine precision_needed(weight, state_size, began, spacing_needed, least_kind)
    real(dp)        , intent(in)  :: weight, state_size, began
    real(real128)   , intent(out) :: spacing_needed
    character(len=:), allocatable, intent(out) :: least_kind
    real(dp) :: target
    select case (stopping_criterion)
    case (relative)
       target = stopping_tolerance * began
    case default
       target = stopping_tolerance
    end select
    spacing_needed = real(target, real128) / real(max(weight * state_size, tiny(1.0_dp)), real128)
    least_kind     = least_kind_for(spacing_needed)
  end subroutine precision_needed
  function consistent_state(physics, degrees, lower, design_value) result(q)
    type(expression)      , intent(in) :: physics
    integer               , intent(in) :: degrees
    real(dp)              , intent(in) :: lower(:), design_value
    real(dp), allocatable :: q(:)
    if (size(lower) /= degrees - 1) then
       error stop 'gti_march: the components below the highest are given, and no others'
    end if
    q = consistent_states(physics, degrees, reshape(lower, [degrees - 1, 1]), design_value)
  end function consistent_state
  function consistent_states(physics, degrees, lower, design_value, spatial_discretization_stencil) result(q)
    type(expression)      , intent(in)           :: physics
    integer               , intent(in)           :: degrees
    real(dp)              , intent(in)           :: lower(:,:), design_value
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
    real(dp), allocatable :: q(:)
    type(stored_directed_graph) :: points
    type(stored_field) :: state, knobs, direction
    class(field), allocatable :: out
    real(dp), allocatable :: below(:), r(:), slope(:), weights(:), e(:)
    real(dp) :: began, target
    integer  :: nodes, i, k, top, iteration
    nodes = size(lower, 2)
    top   = degrees - 1
    if (size(lower, 1) /= top) then
       error stop 'gti_march: the components below the highest are given at every node'
    end if
    allocate(below(nodes), source=0.0_dp)
    if (present(spatial_discretization_stencil)) then
       if (spatial_discretization_stencil % pattern % num_vertices() /= nodes) then
          error stop 'gti_march: the spatial discretization stencil is a stencil over the nodes'
       end if
       call spatial_discretization_stencil % weights % real_vector(weights)
       do k = 1, spatial_discretization_stencil % pattern % num_edges()
          below(spatial_discretization_stencil % pattern % edge_head(k)) = below(spatial_discretization_stencil % pattern % edge_head(k)) &
               & + weights(k) * lower(1, spatial_discretization_stencil % pattern % edge_tail(k))
       end do
    end if
    points = stored_directed_graph(nodes, tails=[integer ::], heads=[integer ::])
    allocate(q(nodes * degrees), source=0.0_dp)
    do i = 1, nodes
       q((i - 1) * degrees + 1:(i - 1) * degrees + top) = lower(:, i)
    end do
    allocate(e(nodes * degrees), source=0.0_dp)
    do i = 1, nodes
       e(i * degrees) = 1.0_dp
    end do
    knobs = stored_field('design', points % vertex_set(), nodes)
    call knobs % set_real_vector(spread(design_value, 1, nodes))
    direction = stored_field('direction', points % vertex_set(), nodes * degrees)
    call direction % set_real_vector(e)
    began = -1.0_dp
    do iteration = 1, stopping_iterations
       state = stored_field('state', points % vertex_set(), nodes * degrees)
       call state % set_real_vector(q)
       call physics % apply(points, [state, knobs], out)
       call out % real_vector(r)
       r = r + below
       if (began < 0.0_dp) began = norm2(r)
       if (stopping_criterion == relative) then
          target = stopping_tolerance * max(began, tiny(1.0_dp))
       else
          target = stopping_tolerance
       end if
       if (norm2(r) <= target) return
       call physics % partial_action(points, [state, knobs], &
            & [variation(physics % argument(1), direction)], out)
       call out % real_vector(slope)
       do i = 1, nodes
          q(i * degrees) = q(i * degrees) - r(i) / slope(i)
       end do
    end do
    write(*,'(a,es12.3)') ' the physics at the initial instant left a residual of ', norm2(r)
    error stop 'gti_march: the initial state is consistent with the physics'
  end function consistent_states
  subroutine frozen_inputs(q, design, num_points, unknowns, inputs)
    real(dp), intent(in) :: q(:), design
    integer , intent(in) :: num_points
    type(stored_directed_graph)    , intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)
    unknowns = stored_directed_graph(size(q), tails=[integer ::], heads=[integer ::])
    allocate(inputs(2))
    inputs(1) = stored_field('state' , unknowns % vertex_set(), size(q))
    inputs(2) = stored_field('design', unknowns % vertex_set(), num_points)
    call inputs(1) % set_real_vector(q)
    call inputs(2) % set_real_vector(spread(design, 1, num_points))
  end subroutine frozen_inputs
  subroutine set_stopping(tolerance, criterion, budget, iterations)
    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, budget, iterations
    if (tolerance <= 0.0_dp) then
       error stop 'gti_march: a tolerance is positive'
    end if
    if (criterion /= relative .and. criterion /= absolute) then
       error stop 'gti_march: a tolerance is measured relative or absolute'
    end if
    if (budget /= by_count .and. budget /= by_rate) then
       error stop 'gti_march: a budget is counted or taken from the rate'
    end if
    if (iterations < 1) then
       error stop 'gti_march: an iteration budget is positive'
    end if
    stopping_tolerance  = tolerance
    stopping_criterion  = criterion
    stopping_budget     = budget
    stopping_iterations = iterations
    call set_linear_stopping(tolerance, criterion, budget)
  end subroutine set_stopping
  pure integer function unknown(instant, degree, degrees, node, nodes) result(at)
    integer, intent(in)           :: instant, degree, degrees
    integer, intent(in), optional :: node, nodes
    integer :: i, m
    i = 1
    m = 1
    if (present(node))  i = node
    if (present(nodes)) m = nodes
    at = ((instant - 1) * m + (i - 1)) * degrees + degree + 1
  end function unknown
  function unknowns_graph(n, degrees) result(g)
    integer, intent(in) :: n, degrees
    type(stored_directed_graph) :: g
    g = stored_directed_graph(n * degrees, tails=[integer ::], heads=[integer ::])
  end function unknowns_graph
  subroutine block_from(tower, b, scheme, physics, held, rows, instants_at)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)  :: b
    class(family)         , intent(in)  :: scheme
    type(expression)      , intent(in)  :: physics
    real(dp)              , intent(in)  :: held(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)
    type(graph), pointer :: horizon, block, slice, moment_node, component, below
    type(coupling_reach), allocatable :: reach(:)
    integer , allocatable :: slice_of(:), member_of(:), members(:), at(:), carried(:)
    integer , allocatable :: r(:), c(:), table(:,:)
    real(dp), allocatable :: dt(:), w(:), dt_weights(:), spatial_weights(:)
    logical , allocatable :: point(:)
    integer :: m, nd, width, n, s, k, j, g, moments, i, d, u, count, e, npts, ncar
    logical :: staged
    nd = physics % equation_degree() + 1
    horizon => level_member(level_member(tower % node(tower % root()), 1), 1)
    block   => level_member(horizon, b)
    n       = level_num_members(block)
    call tower % value_of(block, dt)
    staged  = marches_by_stages(scheme, nd)
    s       = scheme % num_stages()
    m     = tower % extent_of(level_member(first_moment(block, staged), 1))
    width = nd * m
    allocate(members(n))
    do k = 1, n
       members(k) = merge(level_num_members(level_member(block, k)), 1, staged)
    end do
    moments = sum(members)
    allocate(slice_of(moments), member_of(moments), point(moments))
    g = 0
    do k = 1, n
       do j = 1, members(k)
          g = g + 1
          slice_of(g)  = k
          member_of(g) = j
          point(g)     = .not. staged .or. (k > 1 .and. j <= s)
       end do
    end do
    count = moments * width
    below => null()
    allocate(carried(count), at(moments * m))
    ncar = 0
    npts = 0
    do g = 1, moments
       slice => level_member(block, slice_of(g))
       if (staged) then
          moment_node => level_member(slice, member_of(g))
       else
          moment_node => slice
       end if
       do i = 1, m
          do d = 0, nd - 1
             component => level_member(moment_node, d + 1)
             if (tower % status_of(component) == VALUE_KNOWN) then
                ncar = ncar + 1
                carried(ncar) = (g - 1) * width + (i - 1) * nd + d + 1
             end if
          end do
       end do
       if (point(g)) then
          do i = 1, m
             npts = npts + 1
             at(npts) = (g - 1) * width + (i - 1) * nd
          end do
          component => level_member(moment_node, scheme % primary_degree(nd - 1) + 1)
          if (level_couples(component) .and. .not. associated(below)) then
             below => level_coupling(component)
          end if
       end if
    end do
    if (size(held) /= ncar) then
       error stop 'gti_march: one value per carried component'
    end if
    if (staged) then
       call stage_reach_of(tower, block, n, s, nd, width, slice_of, &
            & member_of, reach)
    else
       call block_reach_of(tower, block, n, nd, width, reach)
    end if
    count = 0
    do k = 1, size(reach)
       count = count + size(reach(k) % tails) * m
    end do
    allocate(r(count), c(count), w(count))
    e = 0
    do k = 1, size(reach)
       associate (one => reach(k))
         call weights_of(scheme_weight(scheme), one % vertices, one % tails, one % heads, &
              & dt(one % step_of), one % source_degree, one % determines, dt_weights)
         do i = 1, m
            do u = 1, size(one % tails)
               e    = e + 1
               r(e) = one % row(u)    + (i - 1) * nd
               c(e) = one % column(u) + (i - 1) * nd
               w(e) = dt_weights(u)
            end do
         end do
       end associate
    end do
    rows = block_residual(derived_constraints(r, c, w, moments * width, 'time discretization stencil'), &
         & physics, at(1:npts), moments * width, nd, scheme % primary_degree(nd - 1), &
         & carried(1:ncar), held)
    call rows % placed_on(tower, block)
    call rows % with_reach(reach)
    if (associated(below)) then
       call tower % tuples_of(below, table)
       call tower % value_of(below, spatial_weights)
       call rows % spatial_discretization_laid(stencil(table(2, :), table(1, :), spatial_weights, &
            & spread(0.0_dp, 1, m), 'spatial discretization stencil'))
    end if
    allocate(instants_at(n))
    g = 0
    do k = 1, n
       g = g + members(k)
       instants_at(k) = (g - 1) * width
    end do
  end subroutine block_from
  function first_moment(block, staged) result(moment)
    type(graph), intent(in) :: block
    logical    , intent(in) :: staged
    type(graph), pointer :: moment
    moment => level_member(block, 1)
    if (staged) moment => level_member(moment, 1)
  end function first_moment
  subroutine block_reach_of(tower, block, n, nd, width, reach)
    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: block
    integer        , intent(in) :: n, nd, width
    type(coupling_reach), allocatable, intent(out) :: reach(:)
    integer, allocatable :: table(:,:)
    integer :: e, ne
    call tower % tuples_of(level_coupling(block), table)
    ne = size(table, 2)
    allocate(reach(1))
    reach(1) % vertices = n
    reach(1) % step_of  = [(e, e = 1, n)]
    allocate(reach(1) % tails(ne), reach(1) % heads(ne), reach(1) % source_degree(ne), &
         &   reach(1) % determines(ne), reach(1) % row(ne), reach(1) % column(ne))
    do e = 1, ne
       reach(1) % tails(e)         = (table(1, e) - 1) / nd + 1
       reach(1) % source_degree(e) = mod(table(1, e) - 1, nd)
       reach(1) % heads(e)         = (table(2, e) - 1) / nd + 1
       reach(1) % determines(e)    = mod(table(2, e) - 1, nd)
       reach(1) % column(e) = (reach(1) % tails(e) - 1) * width + reach(1) % source_degree(e) + 1
       reach(1) % row(e)    = (reach(1) % heads(e) - 1) * width + reach(1) % determines(e) + 1
    end do
  end subroutine block_reach_of
  subroutine stage_reach_of(tower, block, n, s, nd, width, slice_of, &
       & member_of, reach)
    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: block
    integer        , intent(in) :: n, s, nd, width, slice_of(:), member_of(:)
    type(coupling_reach), allocatable, intent(out) :: reach(:)
    integer, allocatable :: table(:,:), carry(:,:), first_moment(:), counted(:), filled(:)
    integer :: kk, e, tail_moment, head_moment, vertex_tail, vertex_head
    allocate(reach(n - 1), first_moment(n), counted(n), filled(n))
    first_moment(1) = 1
    do kk = 2, n
       first_moment(kk) = first_moment(kk - 1) + merge(1, s + 1, kk - 1 == 1)
    end do
    counted = 0
    do kk = 2, n
       call tower % tuples_of(level_coupling(level_member(block, kk)), table)
       counted(kk) = size(table, 2)
    end do
    call tower % tuples_of(level_coupling(block), carry)
    do e = 1, size(carry, 2)
       head_moment = (carry(2, e) - 1) / nd + 1
       kk          = slice_of(head_moment)
       counted(kk) = counted(kk) + 1
    end do
    do kk = 2, n
       reach(kk - 1) % vertices = s + 2
       reach(kk - 1) % step_of  = spread(kk, 1, s + 2)
       allocate(reach(kk - 1) % tails(counted(kk)), reach(kk - 1) % heads(counted(kk)), &
            &   reach(kk - 1) % source_degree(counted(kk)), reach(kk - 1) % determines(counted(kk)), &
            &   reach(kk - 1) % row(counted(kk)), reach(kk - 1) % column(counted(kk)))
    end do
    filled = 0
    do kk = 2, n
       call tower % tuples_of(level_coupling(level_member(block, kk)), table)
       do e = 1, size(table, 2)
          filled(kk) = filled(kk) + 1
          vertex_tail = (table(1, e) - 1) / nd + 2
          vertex_head = (table(2, e) - 1) / nd + 2
          call put(reach(kk - 1), filled(kk), vertex_tail, vertex_head, &
               & mod(table(1, e) - 1, nd), mod(table(2, e) - 1, nd), &
               & (first_moment(kk) + vertex_tail - 2 - 1) * width, &
               & (first_moment(kk) + vertex_head - 2 - 1) * width)
       end do
    end do
    do e = 1, size(carry, 2)
       tail_moment = (carry(1, e) - 1) / nd + 1
       head_moment = (carry(2, e) - 1) / nd + 1
       kk          = slice_of(head_moment)
       filled(kk)  = filled(kk) + 1
       call put(reach(kk - 1), filled(kk), 1, member_of(head_moment) + 1, &
            & mod(carry(1, e) - 1, nd), mod(carry(2, e) - 1, nd), &
            & (tail_moment - 1) * width, (head_moment - 1) * width)
    end do
    if (any(filled /= counted)) then
       error stop 'gti_march: every edge of a step is placed once'
    end if
  contains
    subroutine put(one, e, tail, head, source_degree, determines, column_base, row_base)
      type(coupling_reach), intent(inout) :: one
      integer             , intent(in)    :: e, tail, head, source_degree, determines
      integer             , intent(in)    :: column_base, row_base
      one % tails(e)         = tail
      one % heads(e)         = head
      one % source_degree(e) = source_degree
      one % determines(e)    = determines
      one % column(e)        = column_base + source_degree + 1
      one % row(e)           = row_base + determines + 1
    end subroutine put
  end subroutine stage_reach_of
  subroutine solved(rows, design_value, q, achieved, left, seed)
    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: left
    real(dp)       , intent(in) , optional :: seed(:)
    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer :: count, width
    count    = rows % num_unknowns()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), rows % num_points())
    call design % set_real_vector(spread(design_value, 1, rows % num_points()))
    width = rows % num_degrees()
    if (multigrid_on()) call set_aggregates(rows % aggregates(coarse_nodes(rows % num_nodes())))
    call take_inner(solver % inner, count, width)
    call solver % attach(rows, unknowns, unknowns % vertex_set(), count, &
         & held_inputs = [design])
    if (present(seed)) then
       q = seed
    else
       q = at_first_instant(rows, count)
    end if
    solver % compiled       = assembly_present()
    solver % higher_order_jacobian_product = newton_order()
    solver % max_iterations = stopping_iterations
    solver % tolerance      = stopping_tolerance
    solver % criterion      = stopping_criterion
    solver % budget         = stopping_budget
    call solver % solve(spread(0.0_dp, 1, count), q, achieved)
    call keep_inner(solver % inner)
    if (present(left)) then
       left % converged = solver % converged(achieved)
       left % diverging = solver % diverging(achieved)
       left % norm      = achieved
       left % began     = solver % began()
       if (.not. left % converged) call by_aspect(rows, unknowns, q, design, left)
    end if
  end subroutine solved
  integer function fresh_stamp() result(mark)
    stamps_given = stamps_given + 1
    mark = stamps_given
  end function fresh_stamp
  subroutine solved_linear(rows, unknowns, inputs, rhs, transposed, mark, w)
    type(block_residual)       , intent(in)  :: rows
    class(directed_graph)      , intent(in)  :: unknowns
    type(stored_field)         , intent(in)  :: inputs(:)
    real(dp)                   , intent(in)  :: rhs(:)
    logical                    , intent(in)  :: transposed
    integer                    , intent(in)  :: mark
    real(dp), allocatable      , intent(out) :: w(:)
    type(block_residual) :: lin
    real(dp) :: achieved
    if (transposed) then
       call tally_record(adjoint_loops)
    else
       call tally_record(tangent_loops)
    end if
    lin = rows % linear_block(unknowns, inputs, rhs, transposed, mark)
    call swept(lin, 0.0_dp, w, achieved)
  end subroutine solved_linear
  real(dp) function by_tangent(rows, unknowns, inputs, g, design_rate, explicit, mark) &
       & result(df)
    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: mark
    real(dp), allocatable :: w(:)
    call solved_linear(rows, unknowns, inputs, -design_rate, .false., mark, w)
    df = explicit + dot_product(g, w)
  end function by_tangent
  real(dp) function by_adjoint(rows, unknowns, inputs, g, design_rate, explicit, mark) &
       & result(df)
    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: mark
    real(dp), allocatable :: lambda(:)
    call solved_linear(rows, unknowns, inputs, g, .true., mark, lambda)
    df = explicit - dot_product(lambda, design_rate)
  end function by_adjoint
  subroutine set_sweep(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['space-time', 'time      ', 'space     '], 'sweep')
    sweep_level = name
  end subroutine set_sweep
  pure function sweep_named() result(name)
    character(len=:), allocatable :: name
    name = trim(sweep_level)
  end function sweep_named
  subroutine swept(rows, design_value, q, achieved, left)
    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: left
    type(block_residual) :: sub
    type(newton) :: judge
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: design
    integer , allocatable :: member(:), order(:), label(:)
    real(dp), allocatable :: piece(:)
    logical , allocatable :: is_carried(:)
    real(dp) :: sub_achieved, before
    integer :: count, npts, members, m, mm, pass, k
    if (trim(sweep_level) == 'space-time') then
       call solved(rows, design_value, q, achieved, left)
       return
    end if
    count   = rows % num_unknowns()
    npts    = rows % num_points()
    if (trim(sweep_level) == 'time') then
       label = rows % slice_of()
    else
       label = rows % node_of()
    end if
    members = maxval(label)
    allocate(is_carried(count), source=.false.)
    is_carried(rows % carried_unknowns()) = .true.
    q = at_first_instant(rows, count)
    q(rows % carried_unknowns()) = rows % held_values()
    unknowns = stored_directed_graph(count, tails=[integer ::], heads=[integer ::])
    design   = stored_field('nu', unknowns % vertex_set(), npts)
    call design % set_real_vector(spread(design_value, 1, npts))
    judge % max_iterations = stopping_iterations
    judge % tolerance      = stopping_tolerance
    judge % criterion      = stopping_criterion
    judge % budget         = stopping_budget
    call judge % begin_imbalance()
    call rows % member_order(trim(sweep_level) == 'time', order)
    achieved = whole_residual(rows, unknowns, design, q)
    call judge % note_imbalance(achieved)
    do pass = 1, stopping_iterations
       before = achieved
       do mm = 1, members
          m = order(mm)
          member = pack([(k, k = 1, count)], label == m)
          if (all(is_carried(member))) cycle
          if (pass == 1 .and. trim(sweep_level) == 'time' .and. mm > 1) then
             call continued(q, member, pack([(k, k = 1, count)], label == order(mm - 1)))
          end if
          sub = rows % restricted(member, q)
          if (rows % stamp() /= 0) then
             call sub % stamped(abs(rows % stamp()) * members + m, rows % stamp_transposed())
          end if
          call solved(sub, design_value, piece, sub_achieved, seed=q(member))
          q(member) = piece
       end do
       achieved = whole_residual(rows, unknowns, design, q)
       call judge % note_imbalance(achieved)
       if (judge % converged(achieved)) exit
       if (judge % exhausted(pass)) exit
       if (achieved == before) exit
    end do
    if (present(left)) then
       left % converged = judge % converged(achieved)
       left % diverging = judge % diverging(achieved)
       left % norm      = achieved
       left % began     = judge % began()
       if (.not. left % converged) call by_aspect(rows, unknowns, q, design, left)
    end if
  end subroutine swept
  real(dp) function whole_residual(rows, unknowns, design, q) result(norm)
    type(block_residual)       , intent(in) :: rows
    type(stored_directed_graph), intent(in) :: unknowns
    type(stored_field)         , intent(in) :: design
    real(dp)                   , intent(in) :: q(:)
    type(stored_field) :: state
    class(field), allocatable :: out
    real(dp), allocatable :: r(:)
    state = stored_field('state', unknowns % vertex_set(), size(q))
    call state % set_real_vector(q)
    call rows % apply(unknowns, [state, design], out)
    call out % real_vector(r)
    norm = norm2(r)
  end function whole_residual
  subroutine continued(q, member, earlier)
    real(dp), intent(inout) :: q(:)
    integer , intent(in)    :: member(:), earlier(:)
    integer :: pieces, i, width
    if (size(member) == size(earlier)) then
       q(member) = q(earlier)
    else if (mod(size(member), size(earlier)) == 0) then
       width  = size(earlier)
       pieces = size(member) / width
       do i = 1, pieces
          q(member((i - 1) * width + 1:i * width)) = q(earlier)
       end do
    else if (mod(size(earlier), size(member)) == 0) then
       width = size(member)
       q(member) = q(earlier(size(earlier) - width + 1:))
    else
       error stop 'gti_march: a member is seeded from one of its extent, a multiple of it, or a divisor'
    end if
  end subroutine continued
  subroutine by_aspect(rows, unknowns, q, design, left)
    type(block_residual)       , intent(in)    :: rows
    type(stored_directed_graph), intent(in)    :: unknowns
    real(dp)                   , intent(in)    :: q(:)
    type(stored_field)         , intent(in)    :: design
    type(imbalance)            , intent(inout) :: left
    type(stored_field) :: state
    class(field), allocatable :: out
    real(dp), allocatable :: r(:), a(:,:), slope(:)
    integer :: i, d, nd, n
    n  = size(q)
    nd = rows % num_degrees()
    state = stored_field('state', unknowns % vertex_set(), n)
    call state % set_real_vector(q)
    call rows % apply(unknowns, [state, design], out)
    call out % real_vector(r)
    allocate(left % by_degree(0:nd - 1), source=0.0_dp)
    do i = 1, n
       d = mod(i - 1, nd)
       left % by_degree(d) = left % by_degree(d) + r(i) ** 2
    end do
    left % by_degree = sqrt(left % by_degree)
    i = maxloc(abs(r), dim=1)
    left % worst_slot   = (i - 1) / nd + 1
    left % worst_degree = mod(i - 1, nd)
    if (left % norm <= 0.0_dp) return
    call jacobian_of(rows, unknowns, [state, design], n, unknowns % vertex_set(), a)
    slope = matmul(r, a) / left % norm
    i = maxloc(abs(slope), dim=1)
    left % steepest_slot   = (i - 1) / nd + 1
    left % steepest_degree = mod(i - 1, nd)
    left % steepest        = slope(i)
  end subroutine by_aspect
  function at_first_instant(rows, count) result(q)
    type(block_residual), intent(in) :: rows
    integer             , intent(in) :: count
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: one(:)
    integer , allocatable :: at(:)
    integer :: p, nd
    one = rows % first_held()
    nd  = size(one)
    at  = rows % points_at()
    allocate(q(count), source=0.0_dp)
    do p = 1, size(at)
       q(at(p) + 1:at(p) + nd) = one
    end do
  end function at_first_instant
  subroutine horizon_bounds(schemes, added, equation_degree, first, last)
    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:), equation_degree
    integer, allocatable, intent(out) :: first(:), last(:)
    integer :: b
    if (size(schemes) /= size(added)) then
       error stop 'gti_march: one family and one instant count per block'
    end if
    allocate(first(size(added)), last(size(added)))
    do b = 1, size(added)
       if (b == 1) then
          first(b) = 1
          last(b)  = added(b)
       else
          first(b) = last(b - 1) - schemes(b) % scheme % history_depth(equation_degree) + 1
          last(b)  = last(b - 1) + added(b)
       end if
       if (added(b) <= schemes(b) % scheme % history_depth(equation_degree)) then
          error stop 'gti_march: a block adds more instants than its family reaches'
       end if
       if (first(b) < 1) then
          error stop 'gti_march: the horizon holds every instant its blocks reach back over'
       end if
    end do
  end subroutine horizon_bounds
end module gti_march
module gti_adaptive
  use util_precision   , only : dp
  use operation_family , only : family
  use operation_grid   , only : uniform_grid
  use operation_expression, only : expression
  use gti_block        , only : block_residual
  use gti_expansion    , only : expansion, family_holder
  use gti_march        , only : consistent_state, block_from, solved
  implicit none
  private
  public :: adaptive_partition
contains
  subroutine stepped(scheme, physics, degrees, state, h, n, design, arrived)
    class(family)   , intent(in)  :: scheme
    type(expression), intent(in)  :: physics
    integer         , intent(in)  :: degrees, n
    real(dp)        , intent(in)  :: state(:), h, design
    real(dp), allocatable, intent(out) :: arrived(:)
    type(expansion)      :: tower
    type(family_holder)  :: holder(1)
    type(block_residual) :: rows
    integer, allocatable :: at(:)
    real(dp), allocatable :: q(:)
    real(dp) :: achieved
    integer  :: last
    allocate(holder(1) % scheme, source=scheme)
    call tower % build(physics, holder, [n], uniform_grid(h), 0, design)
    call block_from(tower, 1, scheme, physics, state, rows, at)
    call solved(rows, design, q, achieved)
    last    = at(size(at))
    arrived = q(last + 1:last + degrees)
  end subroutine stepped
  pure real(dp) function estimate(coarse, fine, degrees, relative) result(e)
    real(dp), intent(in) :: coarse(:), fine(:)
    integer , intent(in) :: degrees
    logical , intent(in) :: relative
    e = norm2(coarse(1:degrees - 1) - fine(1:degrees - 1))
    if (relative) e = e / max(norm2(fine(1:degrees - 1)), tiny(1.0_dp))
  end function estimate
  function adaptive_partition(scheme, p, physics, degrees, duration, lower, design, &
       & tolerance, relative, rejects) result(dt)
    class(family)   , intent(in)  :: scheme
    integer         , intent(in)  :: p, degrees
    type(expression), intent(in)  :: physics
    real(dp)        , intent(in)  :: duration, lower(:), design, tolerance
    logical         , intent(in)  :: relative
    integer         , intent(out), optional :: rejects
    real(dp), allocatable :: dt(:)
    real(dp), parameter :: safety = 0.9_dp, grow = 5.0_dp, shrink = 0.2_dp
    real(dp), allocatable :: state(:), coarse(:), fine(:)
    real(dp) :: t, h, e, factor
    integer  :: attempt, rejected
    if (scheme % history_depth(degrees - 1) > 1) then
       error stop 'gti_adaptive: an adaptive march is a self-starting scheme'
    end if
    if (duration <= 0.0_dp .or. tolerance <= 0.0_dp) then
       error stop 'gti_adaptive: the duration and the tolerance are positive'
    end if
    state    = consistent_state(physics, degrees, lower, design)
    dt       = [real(dp) ::]
    t        = 0.0_dp
    h        = duration / 8.0_dp
    rejected = 0
    do while (t < duration * (1.0_dp - 1.0e-12_dp))
       h = min(h, duration - t)
       attempt = 0
       do
          attempt = attempt + 1
          call stepped(scheme, physics, degrees, state, h, 2, design, coarse)
          call stepped(scheme, physics, degrees, state, h, 3, design, fine)
          e      = estimate(coarse, fine, degrees, relative)
          factor = safety * (tolerance / max(e, tiny(1.0_dp))) ** (1.0_dp / real(p + 1, dp))
          factor = min(grow, max(shrink, factor))
          if (e <= tolerance .or. h <= duration * 1.0e-10_dp) exit
          rejected = rejected + 1
          h = h * factor
          if (attempt > 50) then
             error stop 'gti_adaptive: a step stays above the tolerance past fifty attempts'
          end if
       end do
       dt    = [dt, h]
       t     = t + h
       state = fine
       h     = h * factor
    end do
    if (present(rejects)) rejects = rejected
  end function adaptive_partition
end module gti_adaptive
module gti_space
  use util_precision            , only : dp
  use view_mesh                 , only : mesh
  use view_mesh_geometry        , only : mesh_from_incidence
  use field_stored              , only : stored_field
  use operation_stencil         , only : stencil
  use operation_diffusion       , only : diffusion_stencil
  use operation_conduction      , only : conduction
  use operation_robin_condition , only : robin_condition, neumann
  use field_forms               , only : polynomial_form
  use view_paraview_writer      , only : paraview_writer, polygon_cell
  use relation_binary           , only : ragged
  use util_string               , only : string
  use operation_grid            , only : uniform_grid, random_grid, partitioned
  use gti_configuration         , only : chosen_from
  implicit none
  private
  public :: room, spatial_mesh, spatial_operator, written_paraview, coarse_cells
  public :: cartesian, circular, elliptical, geometry_of
  integer, parameter :: cartesian  = 1
  integer, parameter :: circular   = 2
  integer, parameter :: elliptical = 3
  type :: room
     integer :: geometry  = cartesian
     integer :: num_cells = 0
     integer :: num_faces = 0
     type(mesh) :: m
     real(dp), allocatable :: corner(:,:)
     integer , allocatable :: first_corner(:)
     integer , allocatable :: cell_corner(:)
     real(dp), allocatable :: centre(:,:)
     real(dp), allocatable :: volume(:)
     integer , allocatable :: cell_ij(:,:)
     integer :: n1 = 0, n2 = 0
  end type room
  type :: face_record
     integer :: tail = 0, head = 0, corner_a = 0, corner_b = 0
  end type face_record
contains
  integer function geometry_of(name) result(geometry)
    character(len=*), intent(in) :: name
    geometry = chosen_from(name, ['cartesian ', 'circular  ', 'elliptical'], 'spatial_geometry')
  end function geometry_of
  pure function mapped(geometry, a, b, xi, eta) result(x)
    integer , intent(in) :: geometry
    real(dp), intent(in) :: a, b, xi, eta
    real(dp) :: x(2)
    real(dp) :: theta
    select case (geometry)
    case (cartesian)
       x = [a * xi, b * eta]
    case (circular)
       theta = 2.0_dp * acos(-1.0_dp) * eta
       x = [a * xi * cos(theta), a * xi * sin(theta)]
    case default
       theta = 2.0_dp * acos(-1.0_dp) * eta
       x = [a * xi * cos(theta), b * xi * sin(theta)]
    end select
  end function mapped
  function spatial_mesh(geometry, a, b, n1, n2, drawn, seed) result(this)
    integer , intent(in) :: geometry, n1, n2, seed
    real(dp), intent(in) :: a, b
    logical , intent(in) :: drawn
    type(room) :: this
    real(dp), allocatable :: xi(:), eta(:), dxi(:), deta(:)
    type(face_record), allocatable :: faces(:)
    integer :: i, j, c, f, polar, cells, ring
    if (n1 < 2 .or. n2 < 2) then
       error stop 'gti_space: at least two cells along each coordinate'
    end if
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
       error stop 'gti_space: an extent is positive'
    end if
    this % geometry = geometry
    polar = merge(1, 0, geometry /= cartesian)
    if (drawn) then
       call partitioned(random_grid(1.0_dp, seed),      n1 + 1, dxi,  xi)
       call partitioned(random_grid(1.0_dp, seed + n1), n2 + 1, deta, eta)
    else
       call partitioned(uniform_grid(1.0_dp), n1 + 1, dxi,  xi)
       call partitioned(uniform_grid(1.0_dp), n2 + 1, deta, eta)
    end if
    allocate(this % corner(2, (n1 + 1 - polar) * (n2 + 1)))
    do j = polar, n1
       do i = 0, n2
          this % corner(:, corner_index(i, j, n2, polar)) = mapped(geometry, a, b, xi(j + 1), eta(i + 1))
       end do
    end do
    if (polar == 1) then
       cells = 1 + (n1 - 1) * n2
    else
       cells = n1 * n2
    end if
    this % num_cells = cells
    allocate(this % first_corner(cells + 1))
    allocate(this % cell_corner(merge(n2 + 4 * (n1 - 1) * n2, 4 * n1 * n2, polar == 1)))
    allocate(this % centre(2, cells), this % volume(cells), this % cell_ij(2, cells))
    this % n1 = n1
    this % n2 = n2
    c = 0
    this % first_corner(1) = 1
    if (polar == 1) then
       c = 1
       do i = 0, n2 - 1
          this % cell_corner(i + 1) = corner_index(i, 1, n2, polar)
       end do
       this % first_corner(2) = n2 + 1
       this % cell_ij(:, 1) = [0, 1]
    end if
    do j = 1 + polar, n1
       do i = 1, n2
          c = c + 1
          call quad(this, c, i, j, n2)
          this % cell_ij(:, c) = [i, j]
       end do
    end do
    allocate(faces(2 * n1 * n2 + 2 * (n1 + n2) + n2))
    f = 0
    do j = 1 + polar, n1 - 1
       do i = 1, n2
          call face_between(faces, f, cell_index(i, j, n2, polar), &
               & cell_index(i, j + 1, n2, polar), corner_index(i - 1, j, n2, polar), &
               & corner_index(i, j, n2, polar))
       end do
    end do
    if (polar == 1) then
       do i = 1, n2
          call face_between(faces, f, 1, cell_index(i, 2, n2, polar), &
               & corner_index(i - 1, 1, n2, polar), corner_index(i, 1, n2, polar))
       end do
    end if
    do j = 1 + polar, n1
       do i = 1, n2 - 1 + polar
          ring = i + 1
          if (ring > n2) ring = 1
          call face_between(faces, f, cell_index(i, j, n2, polar), &
               & cell_index(ring, j, n2, polar), corner_index(i, j - 1, n2, polar), &
               & corner_index(i, j, n2, polar))
       end do
    end do
    do i = 1, n2
       call face_between(faces, f, cell_index(i, n1, n2, polar), 0, &
            & corner_index(i - 1, n1, n2, polar), corner_index(i, n1, n2, polar))
    end do
    if (polar == 0) then
       do i = 1, n2
          call face_between(faces, f, cell_index(i, 1, n2, polar), 0, &
               & corner_index(i - 1, 0, n2, polar), corner_index(i, 0, n2, polar))
       end do
       do j = 1, n1
          call face_between(faces, f, cell_index(1, j, n2, polar), 0, &
               & corner_index(0, j - 1, n2, polar), corner_index(0, j, n2, polar))
          call face_between(faces, f, cell_index(n2, j, n2, polar), 0, &
               & corner_index(n2, j - 1, n2, polar), corner_index(n2, j, n2, polar))
       end do
    end if
    this % num_faces = f
    call measured(this, faces(1:f))
  end function spatial_mesh
  pure integer function corner_index(i, j, n2, polar) result(c)
    integer, intent(in) :: i, j, n2, polar
    c = (j - polar) * (n2 + 1) + i + 1
  end function corner_index
  pure integer function cell_index(i, j, n2, polar) result(c)
    integer, intent(in) :: i, j, n2, polar
    if (polar == 1) then
       c = 1 + (j - 2) * n2 + i
    else
       c = (j - 1) * n2 + i
    end if
  end function cell_index
  subroutine quad(this, c, i, j, n2)
    type(room), intent(inout) :: this
    integer   , intent(in)    :: c, i, j, n2
    integer :: at, polar
    polar = merge(1, 0, this % geometry /= cartesian)
    at = this % first_corner(c)
    this % cell_corner(at)     = corner_index(i - 1, j - 1, n2, polar)
    this % cell_corner(at + 1) = corner_index(i,     j - 1, n2, polar)
    this % cell_corner(at + 2) = corner_index(i,     j,     n2, polar)
    this % cell_corner(at + 3) = corner_index(i - 1, j,     n2, polar)
    this % first_corner(c + 1) = at + 4
  end subroutine quad
  subroutine face_between(faces, f, tail, head, corner_a, corner_b)
    type(face_record), intent(inout) :: faces(:)
    integer          , intent(inout) :: f
    integer          , intent(in)    :: tail, head, corner_a, corner_b
    f = f + 1
    faces(f) = face_record(tail, head, corner_a, corner_b)
  end subroutine face_between
  subroutine measured(this, faces)
    type(room)       , intent(inout) :: this
    type(face_record), intent(in)    :: faces(:)
    integer , allocatable :: cell_vertices(:,:), num_cell_vertices(:)
    integer , allocatable :: face_vertices(:,:), num_face_vertices(:), face_cells(:,:), num_face_cells(:)
    character(len=4), allocatable :: tags(:)
    type(ragged) :: corners
    type(stored_field) :: measure
    real(dp), allocatable :: values(:)
    integer :: f, nf
    nf = size(faces)
    corners = ragged(this % first_corner, this % cell_corner)
    call corners % padded(cell_vertices, num_cell_vertices)
    allocate(face_vertices(2, nf), num_face_vertices(nf), face_cells(2, nf), num_face_cells(nf), tags(nf))
    do f = 1, nf
       face_vertices(:, f)  = [faces(f) % corner_a, faces(f) % corner_b]
       num_face_vertices(f) = 2
       face_cells(:, f)     = [faces(f) % tail, faces(f) % head]
       num_face_cells(f)    = merge(2, 1, faces(f) % head > 0)
       tags(f)              = merge('    ', 'wall', faces(f) % head > 0)
    end do
    this % m = mesh_from_incidence(2, this % corner, cell_vertices, num_cell_vertices, &
         & face_vertices, num_face_vertices, face_cells, num_face_cells, tags)
    measure = this % m % cell_centre()
    call measure % real_vector(values)
    this % centre = reshape(values, [2, this % num_cells])
    measure = this % m % cell_volume()
    call measure % real_vector(this % volume)
  end subroutine measured
  function spatial_operator(this, kappa, degree) result(op)
    type(room), intent(in) :: this
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op
    type(robin_condition) :: wall(1)
    if (degree < 1) then
       error stop 'gti_space: a form of degree below one fits no gradient'
    end if
    wall(1) = neumann('wall', 0.0_dp)
    op = diffusion_stencil(this % m, conduction(kappa), wall, polynomial_form(degree, this % m % dimension))
  end function spatial_operator
  function coarse_cells(this) result(aggregate)
    type(room), intent(in) :: this
    integer, allocatable :: aggregate(:)
    integer :: c, i, j, n2c, polar
    polar = merge(1, 0, this % geometry /= cartesian)
    n2c   = (this % n2 + 1) / 2
    allocate(aggregate(this % num_cells))
    do c = 1, this % num_cells
       i = this % cell_ij(1, c)
       j = this % cell_ij(2, c)
       if (polar == 1 .and. c == 1) then
          aggregate(c) = 1
       else
          aggregate(c) = polar + ((j - 1 - polar) / 2) * n2c + (i - 1) / 2 + 1
       end if
    end do
  end function coarse_cells
  subroutine written_paraview(this, path, names, values)
    type(room)      , intent(in) :: this
    character(len=*), intent(in) :: path, names(:)
    real(dp)        , intent(in) :: values(:,:)
    type(paraview_writer) :: writer
    if (size(values, 1) /= this % num_cells .or. size(values, 2) /= size(names)) then
       error stop 'gti_space: one value per cell per name'
    end if
    writer = paraview_writer(this % m, this % corner, &
         & ragged(this % first_corner, this % cell_corner), &
         & spread(polygon_cell, 1, this % num_cells))
    call writer % write(path, values, string(names))
  end subroutine written_paraview
end module gti_space
module gti_field
  use util_precision   , only : dp
  use operation_stencil, only : stencil
  use field_calculus   , only : field
  use field_stored     , only : stored_field
  use operation_expression, only : expression
  use gti_configuration, only : worded
  use gti_march        , only : consistent_states
  use gti_space        , only : room, spatial_operator, cartesian, written_paraview
  implicit none
  private
  public :: spatial_discretization_stencil_of, initial_field
  public :: against_the_laplacian, against_the_mode, export_instant
contains
  function spatial_discretization_stencil_of(space, kappa, degree) result(op)
    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op
    type(stencil) :: balance
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: lw(:), w(:), held(:)
    integer :: e, m
    balance = spatial_operator(space, kappa, degree)
    m       = balance % pattern % num_edges()
    call balance % weights % real_vector(lw)
    call balance % constants % real_vector(held)
    allocate(r(m), c(m), w(m))
    do e = 1, m
       r(e) = balance % pattern % edge_head(e)
       c(e) = balance % pattern % edge_tail(e)
       w(e) = -lw(e) / space % volume(r(e))
    end do
    op = stencil(r, c, w, held, 'spatial discretization stencil')
  end function spatial_discretization_stencil_of
  function initial_field(physics, degrees, kind, initial_state, design, spatial_discretization_stencil, space, a, b) &
       & result(q)
    type(expression)      , intent(in)           :: physics
    integer               , intent(in)           :: degrees
    character(len=*)      , intent(in)           :: kind, initial_state
    real(dp)              , intent(in)           :: design
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
    type(room)            , intent(in), optional :: space
    real(dp)              , intent(in), optional :: a, b
    real(dp), allocatable :: q(:)
    character(len=32), allocatable :: given(:)
    real(dp), allocatable :: lower(:,:)
    integer  :: nodes, i, d
    nodes = 1
    if (present(space)) nodes = space % num_cells
    allocate(lower(degrees - 1, nodes), source=0.0_dp)
    select case (trim(kind))
    case ('constant')
       given = worded(initial_state)
       if (size(given) > degrees - 1) then
          write(*,'(a,i0,a)') ' the initial state holds the ', degrees - 1, &
               & ' components below the highest, which the physics gives.'
          error stop 'gti_field: the initial state is given below the highest derivative'
       end if
       do d = 1, size(given)
          read(given(d), *) lower(d, 1)
       end do
       do i = 2, nodes
          lower(:, i) = lower(:, 1)
       end do
    case ('mode')
       if (.not. present(space)) error stop 'gti_field: the mode is a field over a mesh'
       if (space % geometry /= cartesian) error stop 'gti_field: the mode is the rectangle''s'
       lower(1, :) = mode_shape(space, a, b)
    case ('bump')
       if (.not. present(space)) error stop 'gti_field: the bump is a field over a mesh'
       lower(1, :) = 1.0_dp + 0.5_dp * mode_shape(space, a, b)
    case default
       error stop 'gti_field: an initial field is constant, the mode, or the bump'
    end select
    q = consistent_states(physics, degrees, lower, design, spatial_discretization_stencil)
  end function initial_field
  pure function mode_shape(space, a, b) result(shape)
    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b
    real(dp), allocatable :: shape(:)
    real(dp) :: pi
    integer  :: i
    pi = acos(-1.0_dp)
    shape = [(cos(pi * space % centre(1, i) / a) * cos(pi * space % centre(2, i) / b), &
         &    i = 1, space % num_cells)]
  end function mode_shape
  subroutine balance_of(space, kappa, degree, values, balanced)
    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa, values(:)
    integer   , intent(in) :: degree
    real(dp), allocatable, intent(out) :: balanced(:)
    type(stencil) :: op
    type(stored_field) :: given
    class(field), allocatable :: out
    op    = spatial_operator(space, kappa, degree)
    given = stored_field('values', op % pattern % vertex_set(), size(values))
    call given % set_real_vector(values)
    call op % apply(op % pattern, [given], out)
    call out % real_vector(balanced)
  end subroutine balance_of
  subroutine against_the_laplacian(space, a, b, kappa, degree)
    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b, kappa
    integer   , intent(in) :: degree
    real(dp), allocatable :: shape(:), balanced(:), exact(:)
    real(dp) :: pi, err(0:2), norm(0:2)
    integer  :: i, walls, count(0:2)
    if (space % geometry /= cartesian) then
       error stop 'gti_field: the laplacian check is the rectangle''s'
    end if
    pi    = acos(-1.0_dp)
    shape = mode_shape(space, a, b)
    exact = -kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2) * shape
    call balance_of(space, kappa, degree, shape, balanced)
    err   = 0.0_dp
    norm  = 0.0_dp
    count = 0
    do i = 1, space % num_cells
       walls = 0
       if (space % cell_ij(1, i) == 1 .or. space % cell_ij(1, i) == space % n2) walls = walls + 1
       if (space % cell_ij(2, i) == 1 .or. space % cell_ij(2, i) == space % n1) walls = walls + 1
       err(walls)   = err(walls)   + (balanced(i) / space % volume(i) - exact(i)) ** 2
       norm(walls)  = norm(walls)  + exact(i) ** 2
       count(walls) = count(walls) + 1
    end do
    write(*,'(a,i0,a,i0,a)') '   the operator against kappa laplacian of the mode, ', &
         & space % num_cells, ' cells, form degree ', degree, ':'
    write(*,'(a,3(a,es10.3))') '   relative rms error', &
         & '   interior ', sqrt(err(0) / max(norm(0), tiny(1.0_dp))), &
         & '   one wall ', sqrt(err(1) / max(norm(1), tiny(1.0_dp))), &
         & '   corner ',   sqrt(err(2) / max(norm(2), tiny(1.0_dp)))
  end subroutine against_the_laplacian
  subroutine against_the_mode(space, a, b, kappa, degree, design, t_last, x, degrees)
    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b, kappa, design, t_last, x(:)
    integer   , intent(in) :: degree, degrees
    real(dp) :: pi, omega, omega_h, exact, semi, e_exact, e_semi, area, mode
    real(dp), allocatable :: shape(:), balanced(:)
    integer  :: i
    if (space % geometry /= cartesian .or. design /= 0.0_dp) return
    pi    = acos(-1.0_dp)
    omega = sqrt(1.0_dp + kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2))
    shape = mode_shape(space, a, b)
    call balance_of(space, kappa, degree, shape, balanced)
    omega_h = sqrt(1.0_dp - dot_product(shape, balanced) / &
         & dot_product(shape, space % volume * shape))
    e_exact = 0.0_dp
    e_semi  = 0.0_dp
    area    = sum(space % volume)
    do i = 1, space % num_cells
       mode  = shape(i)
       exact = mode * cos(omega   * t_last)
       semi  = mode * cos(omega_h * t_last)
       e_exact = e_exact + space % volume(i) * (x((i - 1) * degrees + 1) - exact) ** 2
       e_semi  = e_semi  + space % volume(i) * (x((i - 1) * degrees + 1) - semi) ** 2
    end do
    write(*,'(a,es12.3,a,es12.3,a,f10.6,a,f10.6)') &
         & '      error at the last instant, against the mode ', sqrt(e_exact / area), &
         & '   semi-discrete ', sqrt(e_semi / area), '   omega ', omega, '   omega_h ', omega_h
  end subroutine against_the_mode
  subroutine export_instant(space, path, degrees, x)
    type(room)      , intent(in) :: space
    character(len=*), intent(in) :: path
    integer         , intent(in) :: degrees
    real(dp)        , intent(in) :: x(:)
    character(len=8), allocatable :: names(:)
    real(dp), allocatable :: values(:,:)
    integer :: i, d
    allocate(names(degrees), values(space % num_cells, degrees))
    do d = 0, degrees - 1
       write(names(d + 1),'(a,i0)') 'q', d
    end do
    do i = 1, space % num_cells
       do d = 0, degrees - 1
          values(i, d + 1) = x((i - 1) * degrees + d + 1)
       end do
    end do
    call written_paraview(space, path, names, values)
  end subroutine export_instant
end module gti_field
module gti_chain
  use util_precision  , only : dp
  use operation_family , only : family
  use operation_grid   , only : grid, partitioned
  use operation_expression, only : expression
  use gti_expansion    , only : family_holder, marches_by_stages, expansion, &
       & design_of_physics, design_of_steps
  use gti_block        , only : block_residual
  use gti_march        , only : imbalance, swept, solved_linear, fresh_stamp, horizon_bounds, &
       & frozen_inputs
  use gti_march        , only : block_from
  use gti_sweeps       , only : choose
  use util_derivative_terms, only : derivative_terms, coefficient, operator(*)
  use operation_stencil, only : stencil
  use operation_family_dirk, only : crouzeix_three_stage
  use view_directed_stored, only : stored_directed_graph
  use field_calculus   , only : field
  use field_stored     , only : stored_field
  use gti_sweeps       , only : route_of, route_substitutions, forward_route, reverse_route
  use util_tally            , only : tally_order, tally_enter, tally_leave
  use gti_configuration     , only : at_horizon, at_block, at_stage
  implicit none
  private
  public :: chain_block, march_chain, chain_expansion, instant_components
  public :: first_of, chain_derivative, asymmetry
  public :: multiset_count, multiset_rank, multiset_of, num_designs_of
  public :: chain_stamps
  public :: expansion_substitutions
  public :: sink_costates
  type :: chain_block
     type(block_residual)  :: rows
     integer , allocatable :: instants_at(:)
     real(dp), allocatable :: state(:)
     integer               :: first  = 0
     integer               :: last   = 0
     integer               :: stride = 1
     integer               :: given = 0
     integer               :: primary = 0
     integer :: width = 0
     integer :: nodes = 1
     logical :: staged = .false.
     integer , allocatable :: coarse_step(:)
     class(family), allocatable :: scheme
     real(dp)     , allocatable :: dt(:)
     real(dp)              :: fraction = 1.0_dp
     logical               :: counted = .true.
     real(dp) :: began = 0.0_dp
  end type chain_block
  type :: sink_costates
     integer , allocatable :: carried(:), last(:), interior(:)
     real(dp) :: departure = 0.0_dp
     real(dp) :: gradient  = 0.0_dp
     real(dp) :: costate   = 0.0_dp
     real(dp) :: unread    = 0.0_dp
  end type sink_costates
contains
  pure subroutine locate(chain, fine, held_by, local)
    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: held_by, local
    integer :: b
    held_by = 0
    local   = 0
    do b = size(chain), 1, -1
       if (fine < chain(b) % first .or. fine > chain(b) % last) cycle
       if (mod(fine - chain(b) % first, chain(b) % stride) /= 0) cycle
       held_by = b
       local   = (fine - chain(b) % first) / chain(b) % stride + 1
       return
    end do
  end subroutine locate
  pure function fine_components(chain, fine) result(x)
    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: fine
    real(dp), allocatable :: x(:)
    integer :: b, local, at
    call locate(chain, fine, b, local)
    if (b == 0) error stop 'gti_chain: that instant lies outside the chain'
    at = chain(b) % instants_at(local)
    x  = chain(b) % state(at + 1:at + chain(b) % width)
  end function fine_components
  pure function instant_components(chain, instant) result(x)
    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: instant
    real(dp), allocatable :: x(:)
    x = fine_components(chain, 1 + (instant - 1) * chain(size(chain)) % stride)
  end function instant_components
  pure function handed_over(earlier, first, stride, given) result(held)
    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, stride, given
    real(dp), allocatable :: held(:)
    integer :: i, width
    width = earlier(1) % width
    allocate(held(given * width))
    do i = 1, given
       held((i - 1) * width + 1:i * width) = fine_components(earlier, first + (i - 1) * stride)
    end do
  end function handed_over
  subroutine built(tower, b, scheme, physics, held, rows, instants_at)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)  :: b
    class(family)         , intent(in)  :: scheme
    type(expression)      , intent(in)  :: physics
    real(dp)              , intent(in)  :: held(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)
    call block_from(tower, b, scheme, physics, held, rows, instants_at)
  end subroutine built
  subroutine march_chain(schemes, added, physics, degrees, steps, &
       & design, initial, chain, tower, dt, t, achieved, grid_design, left, nodes, spatial_discretization_stencil, &
       & startup)
    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    type(expression)      , intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    type(expansion), allocatable, intent(inout), target :: tower
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)
    type(imbalance), intent(out), optional :: left
    integer        , intent(in) , optional :: nodes
    type(stencil)  , intent(in) , optional :: spatial_discretization_stencil
    integer        , intent(in) , optional :: startup
    type(family_holder), allocatable :: every(:)
    type(imbalance) :: one_left
    integer , allocatable :: first(:), last(:), spans(:)
    real(dp), allocatable :: fine(:), knobs(:)
    real(dp) :: one_achieved
    integer :: b, k, r, given, before
    logical :: with_startup
    if (size(added) < 1) then
       error stop 'gti_chain: a chain holds at least one block'
    end if
    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partitioned(steps, last(size(added)), dt, t, grid_design)
    given        = schemes(1) % scheme % history_depth(degrees - 1)
    with_startup = .false.
    r            = 1
    if (present(startup)) then
       if (given > 1) then
          with_startup = .true.
          r            = max(startup, 1)
       end if
    end if
    before = merge(1, 0, with_startup)
    allocate(chain(size(added) + before))
    knobs = [real(dp) ::]
    allocate(every(size(added) + before))
    if (with_startup) then
       fine  = [0.0_dp, (dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (given - 1) * r)]
       knobs = fine(2:)
       allocate(every(1) % scheme, source=crouzeix_three_stage())
    end if
    do b = 1, size(added)
       if (b == 1 .and. with_startup) then
          knobs = [knobs, fine(size(fine)), dt(2:last(1))]
       else
          knobs = [knobs, dt(first(b) + merge(1, 0, b == 1):last(b))]
       end if
       allocate(every(before + b) % scheme, source=schemes(b) % scheme)
    end do
    allocate(spans(size(every)))
    if (with_startup) spans(1) = (given - 1) * r + 1
    do b = 1, size(added)
       spans(before + b) = last(b) - first(b) + 1
    end do
    if (allocated(tower)) deallocate(tower)
    allocate(tower)
    call tower % build(physics, every, spans, steps, 0, design, nodes, spatial_discretization_stencil, &
         & weights=grid_design, block_steps=knobs)
    achieved = 0.0_dp
    call tally_enter(at_horizon)
    if (with_startup) then
       call one_block(chain, 1, tower, 1, every(1) % scheme, physics, degrees, 1, &
            & (given - 1) * r + 1, 1, fine, [0, (1 + (k - 1) / r + 1, k = 1, (given - 1) * r)], &
            & 1.0_dp / real(r, dp), .false., design, initial, one_achieved, one_left, nodes, &
            & spatial_discretization_stencil)
       achieved = one_achieved
       if (present(left)) left = one_left
    end if
    do b = 1, size(added)
       call one_block(chain, before + b, tower, before + b, schemes(b) % scheme, physics, degrees, &
            & 1 + (first(b) - 1) * r, 1 + (last(b) - 1) * r, r, dt(first(b):last(b)), &
            & [(k, k = first(b), last(b))], 1.0_dp, .true., design, initial, &
            & one_achieved, one_left, nodes, spatial_discretization_stencil)
       achieved = max(achieved, one_achieved)
       if (present(left)) then
          if (before + b == 1) left = one_left
          if (left % converged .and. .not. one_left % converged) left = one_left
       end if
    end do
    call tally_leave()
  end subroutine march_chain
  subroutine one_block(chain, b, tower, in_tower, scheme, physics, degrees, first, last, &
       & stride, dt, coarse_step, fraction, counted, design, initial, achieved, left, nodes, &
       & spatial_discretization_stencil)
    type(chain_block)     , intent(inout) :: chain(:)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)    :: b, in_tower, degrees, first, last, stride
    integer               , intent(in)    :: coarse_step(:)
    class(family)         , intent(in)    :: scheme
    type(expression)      , intent(in)    :: physics
    real(dp)              , intent(in)    :: dt(:), fraction, design, initial(:)
    logical               , intent(in)    :: counted
    real(dp)              , intent(out)   :: achieved
    type(imbalance)       , intent(out)   :: left
    integer      , intent(in), optional   :: nodes
    type(stencil), intent(in), optional   :: spatial_discretization_stencil
    real(dp), allocatable :: held(:)
    chain(b) % first    = first
    chain(b) % last     = last
    chain(b) % stride   = stride
    chain(b) % given    = scheme % history_depth(degrees - 1)
    chain(b) % primary  = scheme % primary_degree(degrees - 1)
    chain(b) % width    = degrees
    if (present(nodes)) then
       chain(b) % width = degrees * nodes
       chain(b) % nodes = nodes
    end if
    allocate(chain(b) % scheme, source=scheme)
    chain(b) % staged      = marches_by_stages(scheme, degrees)
    chain(b) % dt          = dt
    chain(b) % coarse_step = coarse_step
    chain(b) % fraction    = fraction
    chain(b) % counted     = counted
    if (b == 1) then
       if (size(initial) /= chain(b) % given * chain(b) % width) then
          error stop 'gti_chain: the initial state holds the first block''s given instants'
       end if
       held = initial
    else
       held = handed_over(chain(1:b - 1), first, stride, chain(b) % given)
    end if
    if (scheme % num_stages() > 1) then
       call tally_enter(at_stage)
    else
       call tally_enter(at_block)
    end if
    call built(tower, in_tower, scheme, physics, held, chain(b) % rows, &
         & chain(b) % instants_at)
    associate (u1 => nodes, u2 => spatial_discretization_stencil); end associate
    call swept(chain(b) % rows, design, chain(b) % state, achieved, left)
    chain(b) % began = left % began
    call tally_leave()
  end subroutine one_block
  pure subroutine owned(chain, b, from, to)
    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: b
    integer          , intent(out) :: from, to
    to = size(chain(b) % instants_at)
    if (.not. chain(b) % counted) then
       from = to + 1
       return
    end if
    from = 1
    if (b > 1) then
       if (chain(b - 1) % counted) from = 1 + chain(b) % given
    end if
  end subroutine owned
  subroutine chain_expansion(chain, tower, functionals, degrees, max_order, f, node_measure)
    type(chain_block)      , intent(in) :: chain(:)
    type(expansion)        , intent(in) :: tower
    type(expression)       , intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, max_order
    real(dp), allocatable  , intent(out) :: f(:,:)
    real(dp), intent(in), optional      :: node_measure(:)
    integer , allocatable :: marks(:)
    real(dp), allocatable :: by_order(:,:,:), table(:,:)
    integer :: m
    call chain_stamps(chain, tower, functionals, degrees, marks)
    call chain_derivative(chain, tower, marks, functionals, degrees, max_order, forward_route, &
         & table, node_measure, designs=1, by_order=by_order)
    allocate(f(0:max_order, size(functionals)))
    do m = 0, max_order
       f(m, :) = by_order(:, 1, m)
    end do
  end subroutine chain_expansion
  subroutine frozen_at(b, design, unknowns, inputs)
    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: design
    type(stored_directed_graph), intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)
    call frozen_inputs(b % state, design, b % rows % num_points(), unknowns, inputs)
  end subroutine frozen_at
  subroutine chain_stamps(chain, tower, functionals, degrees, marks, node_measure)
    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    type(expression) , intent(in) :: functionals(:)
    integer          , intent(in) :: degrees
    integer, allocatable, intent(out) :: marks(:)
    real(dp), intent(in), optional :: node_measure(:)
    integer :: b
    associate (u1 => tower, u2 => functionals, u3 => degrees, u4 => node_measure); end associate
    allocate(marks(size(chain)))
    do b = 1, size(chain)
       marks(b) = fresh_stamp()
    end do
  end subroutine chain_stamps
  integer function num_designs_of(tower)
    type(expansion), intent(in) :: tower
    real(dp), allocatable :: step_partials(:,:)
    real(dp) :: design
    call designs_of(tower, design, step_partials)
    num_designs_of = 1
    if (allocated(step_partials)) num_designs_of = 1 + size(step_partials, 2)
  end function num_designs_of
  subroutine designs_of(tower, design, step_partials)
    type(expansion), intent(in) :: tower
    real(dp)       , intent(out) :: design
    real(dp), allocatable, intent(out) :: step_partials(:,:)
    integer :: k
    if (tower % design_kind_of(1) /= design_of_physics) then
       error stop 'gti_chain: the physics'' parameter is the first design'
    end if
    design = tower % parameter()
    do k = 2, tower % num_designs()
       if (tower % design_kind_of(k) == design_of_steps) call tower % step_partials(step_partials)
    end do
  end subroutine designs_of
  pure real(dp) function first_of(table)
    real(dp), intent(in) :: table(:,:)
    first_of = table(1, 1)
  end function first_of
  pure function along_of(b, v) result(along)
    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: v(:)
    real(dp) :: along(size(b % dt))
    integer :: k
    along(1) = 0.0_dp
    do k = 2, size(b % dt)
       along(k) = v(b % coarse_step(k)) * b % fraction
    end do
  end function along_of
  pure subroutine holder_of(chain, fine, held_by, at)
    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: held_by, at
    integer :: local
    call locate(chain, fine, held_by, local)
    at = 0
    if (held_by > 0) at = chain(held_by) % instants_at(local)
  end subroutine holder_of
  subroutine chain_derivative(chain, tower, marks, functionals, degrees, order, route, &
       & table, node_measure, entries, designs, by_order, sinks)
    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    integer          , intent(in) :: marks(:)
    type(expression) , intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, order, route
    real(dp), allocatable  , intent(out) :: table(:,:)
    real(dp), intent(in), optional      :: node_measure(:)
    real(dp), allocatable, intent(out), optional :: entries(:,:,:)
    integer, intent(in), optional :: designs
    real(dp), allocatable, intent(out), optional :: by_order(:,:,:)
    type(sink_costates), intent(out), optional :: sinks
    real(dp), allocatable :: w(:,:,:,:), lambda(:,:,:,:,:), u(:,:,:)
    logical , allocatable :: is_sink(:,:)
    real(dp), allocatable :: diagonal(:,:)
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: step_partials(:,:), rhs(:,:), one(:), r(:), every(:,:,:)
    integer , allocatable :: s(:)
    type(expression) :: physics
    real(dp) :: design
    integer :: nf, nd, nb, top, widest, k, b, i, j, p, d, count, rank, held_by, at, instant
    if (order < 0) then
       error stop 'gti_chain: a derivative has an order of zero or more'
    end if
    if (route /= forward_route .and. route /= reverse_route) then
       error stop 'gti_chain: a route is forward or reverse'
    end if
    call designs_of(tower, design, step_partials)
    nd = 1
    if (allocated(step_partials)) nd = 1 + size(step_partials, 2)
    if (present(designs)) then
       if (designs < 1 .or. designs > nd) then
          error stop 'gti_chain: the designs run over are among the tower''s'
       end if
       nd = designs
    end if
    nf      = size(functionals)
    nb      = size(chain)
    physics = tower % rule()
    widest  = 0
    do b = 1, nb
       widest = max(widest, chain(b) % rows % num_unknowns())
    end do
    top = order
    if (route == reverse_route) top = max(order - 1, 0)
    call steps_along(tower, nd, max(order, 1), u)
    if (present(by_order)) then
       allocate(by_order(nf, multiset_count(nd, order), 0:order), source=0.0_dp)
    end if
    allocate(w(widest, nb, multiset_count(nd, max(top, 1)), max(top, 1)), source=0.0_dp)
    allocate(rhs(widest, nb))
    if (present(by_order)) call functional_tables(0)
    do k = 1, top
       call tally_order(k)
       call tally_enter(at_horizon)
       do rank = 1, multiset_count(nd, k)
          s = multiset_of(rank, k, nd)
          do b = 1, nb
             call tally_enter(at_block)
             count = chain(b) % rows % num_unknowns()
             call rows_along(chain, b, physics, degrees, design, s, w, u, nd, r)
             r = -r
             do i = 1, chain(b) % given * chain(b) % width
                instant = chain(b) % first + ((i - 1) / chain(b) % width) * chain(b) % stride
                d       = mod(i - 1, chain(b) % width)
                call holder_of(chain(1:b - 1), instant, held_by, at)
                if (held_by > 0) r(i) = w(at + d + 1, held_by, rank, k)
             end do
             call frozen_at(chain(b), design, unknowns, inputs)
             call solved_linear(chain(b) % rows, unknowns, inputs, r, .false., marks(b), one)
             w(1:count, b, rank, k) = one
             call tally_leave()
          end do
       end do
       call tally_leave()
       if (present(by_order) .and. k < order) call functional_tables(k)
    end do
    if (route == forward_route .or. order == 0) then
       call functional_tables(order)
       if (present(sinks)) then
          error stop 'gti_chain: the sinks are checked on the reverse route'
       end if
       if (present(by_order)) by_order(:, :, order) = table
       call tally_order(0)
       return
    end if
    if (present(sinks)) then
       allocate(sinks % carried(0:degrees - 1), source=0)
       allocate(sinks % last(0:degrees - 1), source=0)
       allocate(sinks % interior(0:degrees - 1), source=0)
       allocate(is_sink(widest, nb), source=.false.)
       allocate(diagonal(widest, nb), source=0.0_dp)
       do b = 1, nb
          call sinks_of(chain(b), degrees, design, is_sink(:, b), diagonal(:, b), sinks)
       end do
    end if
    allocate(lambda(widest, nb, nf, multiset_count(nd, max(top, 1)), 0:top), source=0.0_dp)
    do k = 0, top
       call tally_order(k + 1)
       call tally_enter(at_horizon)
       do rank = 1, multiset_count(nd, k)
          s = multiset_of(rank, k, nd)
          do i = 1, nf
             rhs = 0.0_dp
             do b = 1, nb
                count = chain(b) % rows % num_unknowns()
                call costate_rows(chain, b, physics, functionals(i), degrees, design, s, &
                     & w, lambda, u, nd, i, node_measure, r)
                rhs(1:count, b) = r
             end do
             do b = nb, 1, -1
                call tally_enter(at_block)
                count = chain(b) % rows % num_unknowns()
                call frozen_at(chain(b), design, unknowns, inputs)
                call solved_linear(chain(b) % rows, unknowns, inputs, rhs(1:count, b), .true., &
                     & marks(b), one)
                lambda(1:count, b, i, rank, k) = one
                if (present(sinks)) then
                   call sink_departure(is_sink(1:count, b), diagonal(1:count, b), &
                        & rhs(1:count, b), one, sinks)
                end if
                call tally_leave()
                do p = 1, chain(b) % given * chain(b) % width
                   instant = chain(b) % first + ((p - 1) / chain(b) % width) * chain(b) % stride
                   d       = mod(p - 1, chain(b) % width)
                   call holder_of(chain(1:b - 1), instant, held_by, at)
                   if (held_by > 0) rhs(at + d + 1, held_by) = rhs(at + d + 1, held_by) + one(p)
                end do
             end do
          end do
       end do
       call tally_leave()
    end do
    call tally_order(0)
    allocate(every(nf, nd, multiset_count(nd, top)), source=0.0_dp)
    do rank = 1, multiset_count(nd, top)
       s = multiset_of(rank, top, nd)
       do j = 1, nd
          do b = 1, nb
             do i = 1, nf
                every(i, j, rank) = every(i, j, rank) + entry_of(chain, b, physics, &
                     & functionals(i), degrees, design, s, j, w, lambda, u, nd, i, node_measure)
             end do
          end do
       end do
    end do
    allocate(table(nf, multiset_count(nd, order)))
    do rank = 1, multiset_count(nd, order)
       s = multiset_of(rank, order, nd)
       table(:, rank) = every(:, s(order), multiset_rank(s(1:order - 1), nd))
    end do
    if (present(entries)) entries = every
    if (present(by_order)) by_order(:, :, order) = table
  contains
    subroutine functional_tables(size_of)
      integer, intent(in) :: size_of
      real(dp), allocatable :: t(:,:)
      integer , allocatable :: s(:)
      integer :: rank, b, i
      allocate(t(nf, multiset_count(nd, size_of)), source=0.0_dp)
      do rank = 1, multiset_count(nd, size_of)
         s = multiset_of(rank, size_of, nd)
         do b = 1, nb
            do i = 1, nf
               t(i, rank) = t(i, rank) + functional_along(chain, b, functionals(i), &
                    & degrees, design, s, 0, w, u, nd, node_measure)
            end do
         end do
      end do
      if (size_of == order) then
         table = t
      else
         by_order(:, 1:size(t, 2), size_of) = t
      end if
    end subroutine functional_tables
  end subroutine chain_derivative
  subroutine steps_along(tower, nd, max_size, u)
    type(expansion), intent(in) :: tower
    integer        , intent(in) :: nd, max_size
    real(dp), allocatable, intent(out) :: u(:,:,:)
    real(dp), allocatable :: column(:)
    integer , allocatable :: s(:)
    integer :: k, rank, n
    n = 1
    if (nd > 1) then
       call tower % step_partial_along([1], column)
       n = size(column)
    end if
    allocate(u(n, multiset_count(nd, max_size), max_size), source=0.0_dp)
    do k = 1, max_size
       if (nd == 1) cycle
       do rank = 1, multiset_count(nd, k)
          s = multiset_of(rank, k, nd)
          if (any(s == 1)) cycle
          call tower % step_partial_along(s - 1, column)
          u(:, rank, k) = column
       end do
    end do
  end subroutine steps_along
  subroutine seeds_of(chain, b, s, open, with_full, w, u, nd, state_seed, step_seed, nu_seed)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, s(:), open, nd
    logical             , intent(in) :: with_full
    real(dp)            , intent(in) :: w(:,:,:,:), u(:,:,:)
    real(dp), allocatable, intent(out) :: state_seed(:,:), step_seed(:,:), nu_seed(:)
    integer, allocatable :: designs(:), t(:)
    integer :: n, full, mask, size_of, count, rank, i
    count = chain(b) % rows % num_unknowns()
    if (open > 0) then
       designs = [s, open]
    else
       designs = s
    end if
    n    = size(designs)
    full = 2**n - 1
    allocate(state_seed(count, 0:full), step_seed(size(chain(b) % dt), max(full, 1)), &
         & nu_seed(max(full, 1)))
    state_seed(:, 0) = chain(b) % state
    step_seed        = 0.0_dp
    nu_seed          = 0.0_dp
    do mask = 1, full
       size_of = popcnt(mask)
       t = sorted(pack(designs, [(btest(mask, i - 1), i = 1, n)]))
       if (size_of == 1 .and. t(1) == 1) nu_seed(mask) = 1.0_dp
       if (.not. any(t == 1)) then
          rank = multiset_rank(t, nd)
          step_seed(:, mask) = along_of(chain(b), u(:, rank, size_of))
       end if
       if (open > 0 .and. btest(mask, n - 1)) then
          state_seed(:, mask) = 0.0_dp
       else if (mask == full .and. .not. with_full) then
          state_seed(:, mask) = 0.0_dp
       else
          state_seed(:, mask) = w(1:count, b, multiset_rank(t, nd), size_of)
       end if
    end do
  end subroutine seeds_of
  pure function sorted(x) result(y)
    integer, intent(in) :: x(:)
    integer :: y(size(x))
    integer :: i, j, held
    y = x
    do i = 2, size(y)
       held = y(i)
       j    = i - 1
       do while (j >= 1)
          if (y(j) <= held) exit
          y(j + 1) = y(j)
          j        = j - 1
       end do
       y(j + 1) = held
    end do
  end function sorted
  function point_terms(rule, degrees, design, at, n, extra, state_seed, nu_seed) result(t)
    type(expression), intent(in) :: rule
    integer         , intent(in) :: degrees, at, n, extra
    real(dp)        , intent(in) :: design, state_seed(:, 0:), nu_seed(:)
    type(derivative_terms) :: t
    type(derivative_terms) :: q(0:degrees - 1), nu
    integer :: d, mask
    do d = 0, degrees - 1
       q(d) = derivative_terms(state_seed(at + d + 1, 0), n + extra)
       do mask = 1, 2**n - 1
          call q(d) % set_coefficient(mask, state_seed(at + d + 1, mask))
       end do
       if (extra > 0) call q(d) % set_direction(n + d + 1, 1.0_dp)
    end do
    nu = derivative_terms(design, n + extra)
    do mask = 1, 2**n - 1
       if (nu_seed(mask) /= 0.0_dp) call nu % set_coefficient(mask, nu_seed(mask))
    end do
    t = rule % at_instant(q, nu)
  end function point_terms
  subroutine quadrature_points(b, k, offset, weight)
    type(chain_block)    , intent(in)  :: b
    integer              , intent(in)  :: k
    integer , allocatable, intent(out) :: offset(:)
    real(dp), allocatable, intent(out) :: weight(:)
    integer :: s, i, width
    if (.not. b % staged) then
       offset = [b % instants_at(k)]
       weight = [1.0_dp]
       return
    end if
    if (k == 1) then
       allocate(offset(0), weight(0))
       return
    end if
    s     = b % scheme % num_stages()
    width = b % width
    allocate(offset(s), weight(s))
    do i = 1, s
       offset(i) = b % instants_at(k) - (s - i + 1) * width
       weight(i) = b % scheme % stage_weight(i)
    end do
  end subroutine quadrature_points
  function measure_terms(b, k, node, n, extra, step_seed, node_measure) result(t)
    type(chain_block), intent(in) :: b
    integer          , intent(in) :: k, node, n, extra
    real(dp)         , intent(in) :: step_seed(:,:)
    real(dp), intent(in), optional :: node_measure(:)
    type(derivative_terms) :: t
    real(dp) :: measure
    integer  :: mask
    measure = 1.0_dp
    if (present(node_measure)) then
       if (size(node_measure) /= b % nodes) then
          error stop 'gti_chain: one measure per node'
       end if
       measure = node_measure(node)
    end if
    t = derivative_terms(b % dt(k), n + extra)
    do mask = 1, 2**n - 1
       call t % set_coefficient(mask, step_seed(k, mask))
    end do
    t = measure * t
  end function measure_terms
  subroutine rows_along(chain, b, physics, degrees, design, s, w, u, nd, r)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b
    type(expression)    , intent(in) :: physics
    integer             , intent(in) :: degrees, s(:), nd
    real(dp)            , intent(in) :: design
    real(dp), intent(in) :: w(:,:,:,:), u(:,:,:)
    real(dp), allocatable, intent(out) :: r(:)
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: carried(:)
    integer :: n, full, e, mask, p, row
    n    = size(s)
    full = 2**n - 1
    call seeds_of(chain, b, s, 0, .false., w, u, nd, state_seed, step_seed, nu_seed)
    call chain(b) % rows % rows_terms(chain(b) % scheme, chain(b) % dt, step_seed, tr, tc, tw)
    carried = is_carried(chain(b))
    at      = chain(b) % rows % points_at()
    allocate(r(size(state_seed, 1)), source=0.0_dp)
    do e = 1, size(tr)
       if (carried(tr(e))) cycle
       do mask = 0, full - 1
          r(tr(e)) = r(tr(e)) + tw(e, ieor(full, mask)) * state_seed(tc(e), mask)
       end do
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (carried(row)) cycle
       r(row) = r(row) + coefficient(point_terms(physics, degrees, design, at(p), n, 0, &
            & state_seed, nu_seed), full)
    end do
  end subroutine rows_along
  pure function is_carried(b) result(carried)
    type(chain_block), intent(in) :: b
    logical, allocatable :: carried(:)
    allocate(carried(b % rows % num_unknowns()), source=.false.)
    carried(b % rows % carried_unknowns()) = .true.
  end function is_carried
  subroutine costates_at(chain, b, s, lambda, nd, i, lam)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, s(:), nd, i
    real(dp)            , intent(in) :: lambda(:,:,:,:,0:)
    real(dp), allocatable, intent(out) :: lam(:,:)
    integer, allocatable :: t(:)
    integer :: n, full, mask, count, k
    n     = size(s)
    full  = 2**n - 1
    count = chain(b) % rows % num_unknowns()
    allocate(lam(count, 0:full))
    do mask = 0, full
       t = pack(s, [(btest(mask, k - 1), k = 1, n)])
       lam(:, mask) = lambda(1:count, b, i, multiset_rank(t, nd), size(t))
    end do
  end subroutine costates_at
  subroutine sinks_of(b, degrees, design, is_sink, diagonal, sinks)
    type(chain_block)  , intent(in)    :: b
    integer            , intent(in)    :: degrees
    real(dp)           , intent(in)    :: design
    logical            , intent(out)   :: is_sink(:)
    real(dp)           , intent(out)   :: diagonal(:)
    type(sink_costates), intent(inout) :: sinks
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    integer , allocatable :: r(:), c(:), reads(:)
    real(dp), allocatable :: w(:)
    logical , allocatable :: has_diagonal(:), carried(:)
    logical :: available
    integer :: n, e, p, d
    call frozen_at(b, design, unknowns, inputs)
    call b % rows % compiled_tangent(unknowns, inputs, 1, r, c, w, available)
    if (.not. available) then
       error stop 'gti_chain: the block compiles its tangent in the state'
    end if
    n = b % rows % num_unknowns()
    allocate(reads(n), source=0)
    allocate(has_diagonal(n), source=.false.)
    is_sink  = .false.
    diagonal = 0.0_dp
    do e = 1, size(r)
       reads(c(e)) = reads(c(e)) + 1
       if (r(e) == c(e)) then
          has_diagonal(c(e)) = .true.
          diagonal(c(e))     = w(e)
       end if
    end do
    is_sink(1:n) = reads == 1 .and. has_diagonal
    carried = is_carried(b)
    do p = 1, n
       if (.not. is_sink(p)) cycle
       d = mod(p - 1, degrees)
       if (carried(p)) then
          sinks % carried(d) = sinks % carried(d) + 1
       else if (p > n - degrees) then
          sinks % last(d) = sinks % last(d) + 1
       else
          sinks % interior(d) = sinks % interior(d) + 1
       end if
    end do
  end subroutine sinks_of
  subroutine sink_departure(is_sink, diagonal, g, lambda, sinks)
    logical            , intent(in)    :: is_sink(:)
    real(dp)           , intent(in)    :: diagonal(:), g(:), lambda(:)
    type(sink_costates), intent(inout) :: sinks
    integer :: p
    sinks % gradient = max(sinks % gradient, maxval(abs(g)))
    sinks % costate  = max(sinks % costate , maxval(abs(lambda)))
    do p = 1, size(g)
       if (.not. is_sink(p)) cycle
       sinks % departure = max(sinks % departure, abs(diagonal(p) * lambda(p) - g(p)))
       if (g(p) == 0.0_dp) sinks % unread = max(sinks % unread, abs(lambda(p)))
    end do
  end subroutine sink_departure
  subroutine costate_rows(chain, b, physics, rule, degrees, design, s, w, lambda, u, nd, i, &
       & node_measure, g)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), nd, i
    type(expression)    , intent(in) :: physics, rule
    real(dp)            , intent(in) :: design
    real(dp), intent(in) :: w(:,:,:,:), lambda(:,:,:,:,0:), u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    real(dp), allocatable, intent(out) :: g(:)
    type(derivative_terms) :: t
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:), lam(:,:)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: carried(:)
    integer , allocatable :: offset(:)
    real(dp), allocatable :: beta(:)
    integer :: n, full, e, mask, p, d, row, k, node, from, to, point, count, pt
    n     = size(s)
    full  = 2**n - 1
    count = chain(b) % rows % num_unknowns()
    call seeds_of(chain, b, s, 0, .true., w, u, nd, state_seed, step_seed, nu_seed)
    allocate(g(count), source=0.0_dp)
    call owned(chain, b, from, to)
    do k = from, to
       call quadrature_points(chain(b), k, offset, beta)
       do pt = 1, size(offset)
          do node = 1, chain(b) % nodes
             point = offset(pt) + (node - 1) * degrees
             t = beta(pt) * measure_terms(chain(b), k, node, n, degrees, step_seed, node_measure) &
                  & * point_terms(rule, degrees, design, point, n, degrees, state_seed, nu_seed)
             do d = 0, degrees - 1
                g(point + d + 1) = g(point + d + 1) + coefficient(t, ior(full, shiftl(1, n + d)))
             end do
          end do
       end do
    end do
    if (n == 0) return
    call costates_at(chain, b, s, lambda, nd, i, lam)
    call chain(b) % rows % rows_terms(chain(b) % scheme, chain(b) % dt, step_seed, tr, tc, tw)
    carried = is_carried(chain(b))
    at      = chain(b) % rows % points_at()
    do e = 1, size(tr)
       if (carried(tr(e))) cycle
       do mask = 1, full
          g(tc(e)) = g(tc(e)) - tw(e, mask) * lam(tr(e), ieor(full, mask))
       end do
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (carried(row)) cycle
       t = point_terms(physics, degrees, design, at(p), n, degrees, state_seed, nu_seed)
       do mask = 1, full
          do d = 0, degrees - 1
             g(at(p) + d + 1) = g(at(p) + d + 1) &
                  & - coefficient(t, ior(mask, shiftl(1, n + d))) * lam(row, ieor(full, mask))
          end do
       end do
    end do
  end subroutine costate_rows
  real(dp) function entry_of(chain, b, physics, rule, degrees, design, s, j, w, lambda, u, nd, i, &
       & node_measure) result(part)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), j, nd, i
    type(expression)    , intent(in) :: physics, rule
    real(dp)            , intent(in) :: design
    real(dp), intent(in) :: w(:,:,:,:), lambda(:,:,:,:,0:), u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    type(derivative_terms) :: t
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:), lam(:,:)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: carried(:)
    integer , allocatable :: offset(:)
    real(dp), allocatable :: beta(:)
    integer :: n, fulln, jbit, full, e, mask, sub, p, row, k, node, from, to, point, pt
    n     = size(s)
    fulln = 2**n - 1
    jbit  = 2**n
    full  = 2**(n + 1) - 1
    call seeds_of(chain, b, s, j, .true., w, u, nd, state_seed, step_seed, nu_seed)
    part = 0.0_dp
    call owned(chain, b, from, to)
    do k = from, to
       call quadrature_points(chain(b), k, offset, beta)
       do pt = 1, size(offset)
          do node = 1, chain(b) % nodes
             point = offset(pt) + (node - 1) * degrees
             t = beta(pt) * measure_terms(chain(b), k, node, n + 1, 0, step_seed, node_measure) &
                  & * point_terms(rule, degrees, design, point, n + 1, 0, state_seed, nu_seed)
             part = part + coefficient(t, full)
          end do
       end do
    end do
    call costates_at(chain, b, s, lambda, nd, i, lam)
    call chain(b) % rows % rows_terms(chain(b) % scheme, chain(b) % dt, step_seed, tr, tc, tw)
    carried = is_carried(chain(b))
    at      = chain(b) % rows % points_at()
    do e = 1, size(tr)
       if (carried(tr(e))) cycle
       do mask = 0, fulln
          sub = mask
          do
             part = part - lam(tr(e), ieor(fulln, mask)) * tw(e, ior(ieor(mask, sub), jbit)) &
                  & * state_seed(tc(e), sub)
             if (sub == 0) exit
             sub = iand(sub - 1, mask)
          end do
       end do
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (carried(row)) cycle
       t = point_terms(physics, degrees, design, at(p), n + 1, 0, state_seed, nu_seed)
       do mask = 0, fulln
          part = part - lam(row, ieor(fulln, mask)) * coefficient(t, ior(mask, jbit))
       end do
    end do
  end function entry_of
  real(dp) function functional_along(chain, b, rule, degrees, design, s, open, w, u, nd, &
       & node_measure) result(part)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), open, nd
    type(expression)    , intent(in) :: rule
    real(dp)            , intent(in) :: design
    real(dp), intent(in) :: w(:,:,:,:), u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    type(derivative_terms) :: t
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:)
    integer , allocatable :: offset(:)
    real(dp), allocatable :: beta(:)
    integer :: n, full, k, node, from, to, point, pt
    n    = size(s) + merge(1, 0, open > 0)
    full = 2**n - 1
    call seeds_of(chain, b, s, open, .true., w, u, nd, state_seed, step_seed, nu_seed)
    part = 0.0_dp
    call owned(chain, b, from, to)
    do k = from, to
       call quadrature_points(chain(b), k, offset, beta)
       do pt = 1, size(offset)
          do node = 1, chain(b) % nodes
             point = offset(pt) + (node - 1) * degrees
             t = beta(pt) * measure_terms(chain(b), k, node, n, 0, step_seed, node_measure) &
                  & * point_terms(rule, degrees, design, point, n, 0, state_seed, nu_seed)
             part = part + coefficient(t, full)
          end do
       end do
    end do
  end function functional_along
  pure integer function multiset_count(designs, size_of)
    integer, intent(in) :: designs, size_of
    multiset_count = choose(designs + size_of - 1, size_of)
  end function multiset_count
  pure integer function multiset_rank(s, designs) result(rank)
    integer, intent(in) :: s(:), designs
    integer :: k, i, y, previous
    k = size(s)
    if (any(s < 1) .or. any(s > designs)) then
       error stop 'gti_chain: a multiset holds designs of the tower'
    end if
    do i = 2, k
       if (s(i) < s(i - 1)) error stop 'gti_chain: a multiset is nondecreasing'
    end do
    rank     = 1
    previous = 1
    do i = 1, k
       do y = previous, s(i) - 1
          rank = rank + choose(designs - y + k - i, k - i)
       end do
       previous = s(i)
    end do
  end function multiset_rank
  pure function multiset_of(rank, size_of, designs) result(s)
    integer, intent(in) :: rank, size_of, designs
    integer :: s(size_of)
    integer :: remaining, y, i, block
    if (rank < 1 .or. rank > multiset_count(designs, size_of)) then
       error stop 'gti_chain: a rank names one of the multisets'
    end if
    remaining = rank - 1
    y         = 1
    do i = 1, size_of
       do
          block = choose(designs - y + size_of - i, size_of - i)
          if (remaining < block) exit
          remaining = remaining - block
          y         = y + 1
       end do
       s(i) = y
    end do
  end function multiset_of
  pure real(dp) function asymmetry(entries, designs, order)
    real(dp), intent(in) :: entries(:,:,:)
    integer , intent(in) :: designs, order
    integer, allocatable :: s(:), rest(:)
    real(dp) :: lowest, highest, value
    integer  :: rank, i, position, k
    asymmetry = 0.0_dp
    do rank = 1, multiset_count(designs, order)
       s = multiset_of(rank, order, designs)
       do i = 1, size(entries, 1)
          lowest  =  huge(1.0_dp)
          highest = -huge(1.0_dp)
          do position = 1, order
             if (position > 1) then
                if (s(position) == s(position - 1)) cycle
             end if
             rest    = pack(s, [(k /= position, k = 1, order)])
             value   = entries(i, s(position), multiset_rank(rest, designs))
             lowest  = min(lowest, value)
             highest = max(highest, value)
          end do
          asymmetry = max(asymmetry, highest - lowest)
       end do
    end do
    asymmetry = asymmetry / max(tiny(1.0_dp), maxval(abs(entries)))
  end function asymmetry
  pure integer function expansion_substitutions(num_blocks, order) result(count)
    integer, intent(in) :: num_blocks, order
    count = num_blocks * route_substitutions(route_of(1, 1, order), 1, 1, order)
  end function expansion_substitutions
end module gti_chain
module gti_driver
  use iso_fortran_env  , only : int64
  use util_precision   , only : dp
  use operation_action , only : jacobian_of
  use operation_grid   , only : grid, uniform_grid, random_grid, partitioned
  use gti_configuration, only : configuration, read_configuration, override
  use view_directed_stored, only : stored_directed_graph
  use field_stored     , only : stored_field
  use operation_family , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage
  use operation_expression  , only : expression
  use gti_physics           , only : van_der_pol_energy, van_der_pol_dissipation
  use gti_chain             , only : chain_block
  implicit none
  private
  public :: settings, chosen_grid, steps_of, clock, cosine, dense_jacobian
  public :: family_named, functional_named
contains
  subroutine settings(default_name, cfg)
    character(len=*)   , intent(in)  :: default_name
    type(configuration), intent(out) :: cfg
    character(len=256) :: argument
    character(len=:), allocatable :: name
    integer :: i
    name = default_name
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(argument, '--config=') == 1) name = trim(argument(10:))
    end do
    call read_configuration(name, cfg)
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(argument, '--config=') == 1) cycle
       call override(cfg, argument)
    end do
  end subroutine settings
  function chosen_grid(cfg) result(steps)
    type(configuration), intent(in) :: cfg
    class(grid), allocatable :: steps
    select case (trim(cfg % grid))
    case ('uniform')
       allocate(steps, source=uniform_grid(cfg % time_duration))
    case ('random')
       allocate(steps, source=random_grid(cfg % time_duration, cfg % seed))
    case default
       error stop 'gti_driver: a grid is uniform or random'
    end select
  end function chosen_grid
  subroutine steps_of(cfg, dt, t)
    type(configuration), intent(in) :: cfg
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    call partitioned(chosen_grid(cfg), cfg % instants, dt, t)
  end subroutine steps_of
  real(dp) function clock() result(s)
    integer(int64) :: ticks, rate
    call system_clock(ticks, rate)
    s = real(ticks, dp) / real(rate, dp)
  end function clock
  pure real(dp) function cosine(d, t) result(q)
    integer , intent(in) :: d
    real(dp), intent(in) :: t
    select case (mod(d, 4))
    case (0)
       q =  cos(t)
    case (1)
       q = -sin(t)
    case (2)
       q = -cos(t)
    case default
       q =  sin(t)
    end select
  end function cosine
  subroutine dense_jacobian(chain, design, a)
    type(chain_block), intent(in) :: chain(:)
    real(dp)         , intent(in) :: design
    real(dp), allocatable, intent(out) :: a(:,:)
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, knobs
    integer :: n
    n        = chain(1) % rows % num_unknowns()
    unknowns = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), n)
    knobs    = stored_field('nu', unknowns % vertex_set(), chain(1) % rows % num_points())
    call state % set_real_vector(chain(1) % state)
    call knobs % set_real_vector(spread(design, 1, chain(1) % rows % num_points()))
    call jacobian_of(chain(1) % rows, unknowns, [state, knobs], n, unknowns % vertex_set(), a)
  end subroutine dense_jacobian
  subroutine family_named(name, order, scheme, ok)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: ok
    ok = .true.
    select case (name)
    case ('bdf')
       allocate(scheme, source=bdf_family(order))
    case ('adams')
       allocate(scheme, source=adams_family(order))
    case ('dirk')
       select case (order)
       case (2)
          allocate(scheme, source=implicit_midpoint())
       case (3)
          allocate(scheme, source=crouzeix_two_stage())
       case (4)
          allocate(scheme, source=crouzeix_three_stage())
       case default
          ok = .false.
       end select
    case default
       ok = .false.
    end select
  end subroutine family_named
  subroutine functional_named(name, degree, rule, ok)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: degree
    type(expression), intent(out) :: rule
    logical         , intent(out) :: ok
    ok = .true.
    select case (name)
    case ('energy')
       rule = van_der_pol_energy(degree)
    case ('dissipation')
       rule = van_der_pol_dissipation(degree)
    case default
       ok = .false.
    end select
  end subroutine functional_named
end module gti_driver
module gti_demos
  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use graph_fractal         , only : graph, branch, known_branch
  use view_sequence         , only : sequence_empty, sequence_first, sequence_rest
  use view_level            , only : level_storage, level_consistent, &
       & level_is_leaf, level_num_members, level_members, level_couples, level_coupling
  use view_relational       , only : relational_binding, num_member_sets, &
       & num_relations, relational_valid, relation_at
  use relation_finitary     , only : relation
  use relation_binary       , only : csr_relation
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation
  use map_set_store         , only : set_store
  use map_value             , only : value_map, VALUE_UNATTACHED, VALUE_UNKNOWN, &
       & VALUE_KNOWN
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_action      , only : variation, sweep_design_partial => design_partial
  use gti_sweeps            , only : functional_of, functional_gradient
  use operation_stencil     , only : stencil
  use operation_scheme_stencil, only : derived_constraints
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : dirk_family, implicit_midpoint, &
       & crouzeix_two_stage, crouzeix_three_stage
  use operation_grid        , only : grid, uniform_grid, random_grid, designed_grid, partition
  use operation_coupling    , only : weights_of, coupling_inputs
  use operation_weight      , only : scheme_weight
  use operation_expression  , only : expression
  use gti_physics           , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use operation_minimization, only : relative, by_rate
  use gti_expansion         , only : expansion, family_holder
  use gti_block             , only : block_residual
  use gti_march             , only : block_from, solved, unknowns_graph, &
       & horizon_bounds, set_stopping, consistent_state, imbalance, by_tangent, &
       & by_adjoint, fresh_stamp
  use gti_adaptive          , only : adaptive_partition
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & chain_stamps, chain_derivative, first_of, instant_components, &
       & asymmetry, multiset_count, multiset_rank, multiset_of
  use gti_sweeps            , only : route_of, forward_route, reverse_route
  use gti_driver            , only : clock, cosine, dense_jacobian
  implicit none
  private
  public :: demo_requested, run_demo
  character(len=24), parameter :: demo_names(21) = [character(len=24) :: &
       & 'adaptive_grid', 'assembled_tower', 'chained_horizon', &
       & 'constraint_rows', 'coupling_relation', 'expansion_check', &
       & 'family_coefficients', 'function_identities', 'grid_design_check', &
       & 'jacobian_shape', 'level_maps', 'level_shape', 'marched_block', &
       & 'marched_horizon', 'marched_stages', 'memory_shape', &
       & 'randomized_checks', 'scheme_weights', 'sensitivity', 'solve_cost', &
       & 'tolerance_form']
contains
  logical function demo_requested() result(yes)
    character(len=256) :: argument
    integer :: i
    yes = .false.
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(trim(argument), '--demo=') == 1 .or. trim(argument) == '--list-demos') yes = .true.
    end do
  end function demo_requested
  subroutine run_demo()
    character(len=:), allocatable :: name
    name = demo_name()
    select case (name)
    case ('list', '')
       call list_demos()
    case ('adaptive_grid')
       call demo_adaptive_grid()
    case ('assembled_tower')
       call demo_assembled_tower()
    case ('chained_horizon')
       call demo_chained_horizon()
    case ('constraint_rows')
       call demo_constraint_rows()
    case ('coupling_relation')
       call demo_coupling_relation()
    case ('expansion_check')
       call demo_expansion_check()
    case ('family_coefficients')
       call demo_family_coefficients()
    case ('function_identities')
       call demo_function_identities()
    case ('grid_design_check')
       call demo_grid_design_check()
    case ('jacobian_shape')
       call demo_jacobian_shape()
    case ('level_maps')
       call demo_level_maps()
    case ('level_shape')
       call demo_level_shape()
    case ('marched_block')
       call demo_marched_block()
    case ('marched_horizon')
       call demo_marched_horizon()
    case ('marched_stages')
       call demo_marched_stages()
    case ('memory_shape')
       call demo_memory_shape()
    case ('randomized_checks')
       call demo_randomized_checks()
    case ('scheme_weights')
       call demo_scheme_weights()
    case ('sensitivity')
       call demo_sensitivity()
    case ('solve_cost')
       call demo_solve_cost()
    case ('tolerance_form')
       call demo_tolerance_form()
    case default
       write(*,'(a)') ' unknown demo: ' // name
       call list_demos()
       error stop 'graph_time_integrator: unknown demo'
    end select
  end subroutine run_demo
  function demo_name() result(name)
    character(len=:), allocatable :: name
    character(len=256) :: argument
    integer :: i, j
    name = ''
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(trim(argument), '--demo=') == 1) name = trim(argument(8:))
       if (trim(argument) == '--list-demos') name = 'list'
    end do
    do j = 1, len(name)
       if (name(j:j) == '-') name(j:j) = '_'
    end do
  end function demo_name
  integer function demo_argument_count() result(count)
    character(len=256) :: argument
    integer :: i
    count = 0
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(trim(argument), '--demo=') == 1 .or. trim(argument) == '--list-demos') cycle
       count = count + 1
    end do
  end function demo_argument_count
  subroutine demo_argument(which, argument)
    integer, intent(in) :: which
    character(len=*), intent(out) :: argument
    character(len=256) :: held
    integer :: i, count
    argument = ''
    count = 0
    do i = 1, command_argument_count()
       call get_command_argument(i, held)
       if (index(trim(held), '--demo=') == 1 .or. trim(held) == '--list-demos') cycle
       count = count + 1
       if (count == which) then
          argument = held
          return
       end if
    end do
  end subroutine demo_argument
  subroutine list_demos()
    integer :: i
    write(*,'(a)') ' demos:'
    do i = 1, size(demo_names)
       write(*,'(a)') '   ' // trim(demo_names(i))
    end do
  end subroutine list_demos
  type(family_holder) function held_family(scheme) result(held)
    class(family), intent(in) :: scheme
    call set_family(held, scheme)
  end function held_family
  subroutine set_family(held, scheme)
    type(family_holder), intent(inout) :: held
    class(family)     , intent(in)    :: scheme
    if (allocated(held % scheme)) deallocate(held % scheme)
    allocate(held % scheme, source=scheme)
  end subroutine set_family
  function cosine_history(scheme, degrees, t) result(held)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees
    real(dp)     , intent(in) :: t(:)
    real(dp), allocatable :: held(:)
    integer :: k, d
    held = [((cosine(d, t(k)), d = 0, degrees - 1), &
         & k = 1, scheme % history_depth(degrees - 1))]
  end function cosine_history
  subroutine cosine_partition(scheme, degrees, duration, instants, held, dt, t)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, instants
    real(dp)     , intent(in) :: duration
    real(dp), allocatable, intent(out) :: held(:), dt(:), t(:)
    call partition(duration, instants, dt, t)
    held = cosine_history(scheme, degrees, t)
  end subroutine cosine_partition
  subroutine demo_adaptive_grid()
    implicit none
    integer , parameter :: degrees  = 3          ! van der Pol is degree two
    real(dp), parameter :: duration = 4.0_dp
    real(dp), parameter :: design   = 1.0_dp
    real(dp), parameter :: q0 = 2.0_dp, qd0 = 0.0_dp
    call report('implicit midpoint (order 2)', 2)
    call report('crouzeix two stage (order 3)', 3)
    call report('crouzeix three stage (order 4)', 4)
  contains
    function dirk_of(order) result(scheme)
      integer, intent(in) :: order
      class(family), allocatable :: scheme
      select case (order)
      case (2); allocate(scheme, source=implicit_midpoint())
      case (3); allocate(scheme, source=crouzeix_two_stage())
      case (4); allocate(scheme, source=crouzeix_three_stage())
      case default; error stop 'adaptive_grid: order two, three or four'
      end select
    end function dirk_of
    subroutine on_grid(scheme, dt, f, forward, reverse)
      class(family), intent(in)  :: scheme
      real(dp)     , intent(in)  :: dt(:)
      real(dp)     , intent(out) :: f, forward, reverse
      type(family_holder)     :: schemes(1)
      type(expression)       :: functionals(1)
      type(chain_block), allocatable :: chain(:)
      integer, allocatable :: marks(:)
      type(expansion), allocatable :: tower
      real(dp), allocatable :: grid_dt(:), t(:), fvals(:,:), df(:,:), other(:,:)
      real(dp) :: achieved
      integer  :: n
      n = size(dt) + 1
      call set_family(schemes(1), scheme)
      functionals(1) = van_der_pol_energy(degrees - 1)
      call march_chain(schemes, [n - 1], van_der_pol(degrees - 1), degrees, &
           & designed_grid(duration), design, &
           & consistent_state(van_der_pol(degrees - 1), degrees, [q0, qd0], design), &
           & chain, tower, grid_dt, t, achieved, grid_design=dt)
      call chain_expansion(chain, tower, functionals, degrees, 1, fvals)
      f = fvals(0, 1)
      call chain_stamps(chain, tower, functionals, degrees, marks)
      call chain_derivative(chain, tower, marks, functionals, degrees, 1, forward_route, df)
      call chain_derivative(chain, tower, marks, functionals, degrees, 1, reverse_route, other)
      forward = df(1, 1)
      reverse = other(1, 1)
    end subroutine on_grid
    subroutine report(title, order)
      character(len=*), intent(in) :: title
      integer         , intent(in) :: order
      class(family), allocatable :: scheme
      real(dp), allocatable :: dt(:)
      real(dp) :: tol, f, forward, reverse, span
      integer  :: rejects, level
      scheme = dirk_of(order)
      write(*,'(a)') ' '
      write(*,'(a)') ' ' // title
      write(*,'(a)') '   tolerance     steps   rejects        sum dt - T          functional     forward-reverse'
      do level = 1, 4
         tol = 10.0_dp ** (-3 - level)
         dt = adaptive_partition(scheme, order, van_der_pol(degrees - 1), degrees, duration, [q0, qd0], design, tol, .true., rejects)
         call on_grid(scheme, dt, f, forward, reverse)
         span = sum(dt) - duration
         write(*,'(a,es9.1,i9,i9,es18.2,f18.9,es18.2)') '   ', tol, size(dt), rejects, span, f, &
              & abs(forward - reverse) / max(1.0_dp, abs(forward))
      end do
    end subroutine report
  end subroutine demo_adaptive_grid
  subroutine demo_assembled_tower()
    implicit none
    real(dp), parameter :: duration = 7.0_dp
    integer , parameter :: seed = 20260824
    type(expansion) :: one, two
    type(family_holder) :: schemes(2)
    schemes = [held_family(bdf_family(2)), held_family(adams_family(3))]
    call one % build(van_der_pol(2), schemes, [5, 5], random_grid(duration, seed), 0, &
         & 0.0_dp)
    write(*,'(a)') ' a horizon of two blocks, marched by different families'
    call show(one, one % node(one % root()), 0)
    write(*,'(a)')      ' '
    write(*,'(a,i0)')   ' nodes owned                        ', one % num_nodes()
    write(*,'(a,l1)')   ' every level and coupling valid     ', &
         & one % consistent(one % node(one % root()))
    call two % build(van_der_pol(2), schemes, [5, 5], random_grid(duration, seed), 1, &
         & 0.0_dp)
    write(*,'(a,i0)')   ' nodes owned with one tangent sweep ', two % num_nodes()
    write(*,'(a,l1)')   ' every level and coupling valid     ', &
         & two % consistent(two % node(two % root()))
    call grid_partials()
    call gauss_exactness()
    call stage_block()
  contains
    recursive subroutine show(tower, g, depth)
      type(expansion), intent(in) :: tower
      type(graph)    , intent(in) :: g
      integer        , intent(in) :: depth
      type(graph), pointer :: coupling
      write(*,'(a,a,a,a)') repeat('   ', depth + 1), tower % label_of(g), &
           & status(tower, g), extent(tower, g)
      if (level_couples(g)) then
         coupling => level_coupling(g)
         write(*,'(a,a,a,a)') repeat('   ', depth + 2), tower % label_of(coupling), &
              & status(tower, coupling), extent(tower, coupling)
      end if
      if (level_is_leaf(g)) return
      call show_each(tower, level_members(g), depth + 1)
    end subroutine show
    recursive subroutine show_each(tower, members, depth)
      type(expansion), intent(in) :: tower
      type(branch)   , intent(in) :: members
      integer        , intent(in) :: depth
      type(graph), pointer :: first
      if (sequence_empty(members)) return
      first => sequence_first(members)
      call show(tower, first, depth)
      call show_each(tower, sequence_rest(members), depth)
    end subroutine show_each
    function status(tower, g) result(text)
      type(expansion), intent(in) :: tower
      type(graph)    , intent(in) :: g
      character(len=:), allocatable :: text
      real(dp), allocatable :: x(:)
      select case (tower % status_of(g))
      case (VALUE_KNOWN)
         call tower % value_of(g, x)
         text = '   holds ' // count_of(size(x))
      case (VALUE_UNKNOWN)
         text = '   not yet known'
      case (VALUE_UNATTACHED)
         text = ''
      case default
         error stop 'assembled_tower: a value status is one of the three'
      end select
    end function status
    function extent(tower, g) result(text)
      type(expansion), intent(in) :: tower
      type(graph)    , intent(in) :: g
      character(len=:), allocatable :: text
      text = ''
      if (tower % extent_of(g) > 0) text = '   extent ' // count_of(tower % extent_of(g))
    end function extent
    function count_of(n) result(text)
      integer, intent(in) :: n
      character(len=:), allocatable :: text
      character(len=12) :: buffer
      write(buffer,'(i0)') n
      text = trim(buffer)
    end function count_of
    subroutine grid_partials()
      integer , parameter :: num_instants = 6
      real(dp), parameter :: delta = 1.0e-6_dp
      type(grid) :: steps
      type(stored_directed_graph) :: instants
      type(stored_field) :: knobs, direction
      class(field), allocatable :: out
      real(dp) :: design(num_instants - 1), v(num_instants - 1)
      real(dp), allocatable :: dt(:), exact(:), plus(:), minus(:)
      design = [1.0_dp, 2.0_dp, 1.5_dp, 0.5_dp, 3.0_dp]
      v      = 0.0_dp
      v(2)   = 1.0_dp
      steps    = designed_grid(duration)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      knobs     = stored_field('design', instants % vertex_set(), size(design))
      direction = stored_field('v', instants % vertex_set(), size(design))
      call knobs     % set_real_vector(design)
      call direction % set_real_vector(v)
      call steps % apply(instants, [knobs], out)
      call out % real_vector(dt)
      call steps % partial_action(instants, [knobs], &
           & [variation(steps % argument(1), direction)], out)
      call out % real_vector(exact)
      call differenced(steps, instants, knobs, design, v, delta, plus, minus)
      call show_grid(dt, exact, (plus - minus) / (2.0_dp * delta))
    end subroutine grid_partials
    !-----------------------------------------------------------------!
    ! The quadrature kind of the same map: the n-point rule integrates
    ! t to every power below 2n exactly over [0, span]. The floor each
    ! difference is held to is the arithmetic's: n products, each
    ! within one epsilon of the integral's size.
    !-----------------------------------------------------------------!

    subroutine gauss_exactness()

      use operation_grid, only : gauss_grid, partitioned

      real(dp), parameter :: span = 2.0_dp
      real(dp), allocatable :: w(:), t(:)
      real(dp) :: quadrature, integral, floor, worst
      integer :: n, m, k

      write(*,'(a)') ' '
      write(*,'(a)') ' the gauss kind of the grid, exact to twice its points less one'
      do n = 3, 5, 2
         call partitioned(gauss_grid(span), n, w, t)
         worst = 0.0_dp
         do m = 0, 2 * n - 1
            quadrature = 0.0_dp
            do k = 1, n
               quadrature = quadrature + w(k) * t(k) ** m
            end do
            integral = span ** (m + 1) / real(m + 1, dp)
            floor    = real(n, dp) * epsilon(1.0_dp) * integral
            worst    = max(worst, abs(quadrature - integral) / floor)
         end do
         write(*,'(a,i0,a,i0,a,f8.3,a,es9.2)') '   ', n, ' points, powers 0 to ', &
              & 2 * n - 1, ':  sum of weights - span ', sum(w) - span, &
              & '   worst difference over its floor ', worst
      end do

    end subroutine gauss_exactness

    subroutine show_grid(dt, exact, differenced)
      real(dp), intent(in) :: dt(:), exact(:), differenced(:)
      write(*,'(a)')          ' '
      write(*,'(a)')          ' a designed grid of five steps'
      write(*,'(a,6f10.5)')   '   steps                       ', dt
      write(*,'(a,f10.5)')    '   their sum                   ', sum(dt)
      write(*,'(a,6f10.5)')   '   partial in design 2         ', exact
      write(*,'(a,6f10.5)')   '   central difference          ', differenced
      write(*,'(a,es10.2)')   '   the sum is unchanged along it', sum(exact)
    end subroutine show_grid
    subroutine differenced(steps, instants, knobs, design, v, delta, plus, minus)
      type(grid)                 , intent(in)    :: steps
      type(stored_directed_graph), intent(in)    :: instants
      type(stored_field)         , intent(inout) :: knobs
      real(dp)                   , intent(in)    :: design(:), v(:), delta
      real(dp), allocatable      , intent(out)   :: plus(:), minus(:)
      class(field), allocatable :: out
      call knobs % set_real_vector(design + delta * v)
      call steps % apply(instants, [knobs], out)
      call out % real_vector(plus)
      call knobs % set_real_vector(design - delta * v)
      call steps % apply(instants, [knobs], out)
      call out % real_vector(minus)
      call knobs % set_real_vector(design)
    end subroutine differenced
    subroutine stage_block()
      integer , parameter :: num_instants = 3
      real(dp), parameter :: gamma = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp
      type(expansion) :: staged
      type(family_holder) :: schemes(1)
      type(graph), pointer :: g, coupling
      real(dp), allocatable :: w(:)
      real(dp) :: step
      schemes = [held_family(crouzeix_two_stage())]
      call staged % build(van_der_pol(2), schemes, [num_instants], &
           & uniform_grid(duration), 0, 0.0_dp)
      write(*,'(a)') ' '
      write(*,'(a)') ' a block marched by a two-stage tableau'
      call show(staged, staged % node(staged % root()), 0)
      step = duration / real(num_instants - 1, dp)
      g => second_slice(staged)
      coupling => level_coupling(g)
      call staged % value_of(coupling, w)
      write(*,'(a)')        ' '
      write(*,'(a,l1)')     ' every level and coupling valid     ', &
           & staged % consistent(staged % node(staged % root()))
      write(*,'(a,i0)')     ' weights on one step                ', size(w)
      write(*,'(a,5f10.5)') '   value row, a_11 a_21 a_22 b_1 b_2', w(1:5)
      write(*,'(a,5f10.5)') '   step times the tableau           ', &
           & step * [gamma, 1.0_dp - 2.0_dp * gamma, gamma, 0.5_dp, 0.5_dp]
      write(*,'(a,2f10.5)') '   recovery of the highest degree   ', w(size(w) - 1:)
      write(*,'(a,2f10.5)') '   the tableau weights themselves   ', [0.5_dp, 0.5_dp]
    end subroutine stage_block
    function second_slice(tower) result(g)
      type(expansion), intent(in) :: tower
      type(graph), pointer :: g
      type(branch) :: members
      members = level_members(tower % node(tower % root()))   ! sweeps
      g => sequence_first(members)
      g => sequence_first(level_members(g))                   ! the horizon's blocks
      g => sequence_first(level_members(g))                   ! the block
      g => sequence_first(sequence_rest(level_members(g)))    ! its second slice
    end function second_slice
  end subroutine demo_assembled_tower
  subroutine demo_chained_horizon()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    integer , parameter :: max_order = 3
    real(dp), parameter :: duration = 2.0_dp
    real(dp), parameter :: delta = 1.0e-4_dp
    call splitting_changes_nothing()
    call across_families('bdf 2 then crouzeix two-stage', &
         & held_family(bdf_family(2)), held_family(crouzeix_two_stage()))
    call across_families('crouzeix two-stage then bdf 2', &
         & held_family(crouzeix_two_stage()), held_family(bdf_family(2)))
    call across_families('adams 3 then crouzeix two-stage', &
         & held_family(adams_family(3)), held_family(crouzeix_two_stage()))
  contains
    subroutine expanded(schemes, added, design, f)
      type(family_holder), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp)           , intent(in) :: design
      real(dp), allocatable, intent(out) :: f(:)
      real(dp), allocatable :: table(:,:)
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      real(dp), allocatable :: dt(:), t(:), held(:)
      real(dp) :: achieved
      call cosine_partition(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)
      call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(state_degree)], degrees, max_order, table)
      allocate(f(lbound(table, 1):ubound(table, 1)))
      f = table(:, 1)
    end subroutine expanded
    subroutine splitting_changes_nothing()
      type(family_holder) :: whole(1), split(2)
      real(dp), allocatable :: f_whole(:), f_split(:)
      whole = [held_family(bdf_family(2))]
      split = [held_family(bdf_family(2)), held_family(bdf_family(2))]
      call expanded(whole, [20], 1.0_dp, f_whole)
      call expanded(split, [10, 10], 1.0_dp, f_split)
      write(*,'(a)')        ' bdf 2 over twenty instants, in one block and in two'
      write(*,'(a,4f14.8)') '   whole                    ', f_whole
      write(*,'(a,4f14.8)') '   split                    ', f_split
      write(*,'(a,4es14.2)')'   apart                    ', abs(f_whole - f_split)
    end subroutine splitting_changes_nothing
    subroutine across_families(title, one, two)
      character(len=*)   , intent(in) :: title
      type(family_holder), intent(in) :: one, two
      type(family_holder) :: schemes(2)
      real(dp), allocatable :: f(:), plus(:), minus(:)
      real(dp) :: differenced(max_order)
      integer :: m
      schemes(1) = one
      schemes(2) = two
      call expanded(schemes, [10, 10], 1.0_dp, f)
      call expanded(schemes, [10, 10], 1.0_dp + delta, plus)
      call expanded(schemes, [10, 10], 1.0_dp - delta, minus)
      do m = 1, max_order
         differenced(m) = (plus(m - 1) - minus(m - 1)) / (2.0_dp * delta)
      end do
      write(*,'(a)')            ' '
      write(*,'(a)')            ' ' // title
      write(*,'(a,4i14)')       '   derivative order         ', [(m, m = 0, max_order)]
      write(*,'(a,4f14.8)')     '   from the expansion       ', f
      write(*,'(a,14x,3f14.8)') '   differenced              ', differenced
      write(*,'(a,14x,3es14.2)')'   apart                    ', &
           & abs(f(1:max_order) - differenced)
    end subroutine across_families
  end subroutine demo_chained_horizon
  subroutine demo_constraint_rows()
    implicit none
    integer , parameter :: order = 2
    integer , parameter :: num_instants = 5
    integer , parameter :: num_degrees = 3
    integer , parameter :: num_unknowns = num_instants * num_degrees
    real(dp), parameter :: nu = 1.0_dp
    call derived_rows()
    call physics_partials(2, [0.7_dp, -0.4_dp], [0.3_dp, 0.9_dp], [1.1_dp, -0.5_dp])
    call physics_partials(3, [0.7_dp, -0.4_dp], [0.3_dp, 0.9_dp], [1.1_dp, -0.5_dp])
  contains
    pure integer function unknown(instant, degree)
      integer, intent(in) :: instant, degree
      unknown = (instant - 1) * num_degrees + degree + 1
    end function unknown
    subroutine derived_rows()
      real(dp), parameter :: dt(num_instants) = [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp]
      type(stencil) :: rows
      type(stored_directed_graph) :: unknowns
      type(stored_field) :: state, direction
      class(field), allocatable :: out
      type(expression) :: physics
      real(dp), allocatable :: weight(:), residual(:), acted(:), governing(:)
      real(dp) :: q(num_unknowns), t(num_instants)
      integer , allocatable :: tails(:), heads(:), determines(:), source_degree(:)
      integer :: j, k
      call scheme_reach(tails, heads, source_degree, determines)
      call edge_weights(dt, tails, heads, source_degree, determines, weight)
      rows = derived_constraints( &
           & [(unknown(heads(j), determines(j)), j = 1, size(heads))], &
           & [(unknown(tails(j), source_degree(j)), j = 1, size(tails))], &
           & weight, num_unknowns, 'derived constraints')
      t(1) = 0.0_dp
      do k = 2, num_instants
         t(k) = t(k - 1) + dt(k)
      end do
      do k = 1, num_instants
         q(unknown(k, 0)) = t(k) * t(k)
         q(unknown(k, 1)) = 2.0_dp * t(k)
         q(unknown(k, 2)) = 2.0_dp
      end do
      unknowns = stored_directed_graph(num_unknowns, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', unknowns % vertex_set(), num_unknowns)
      call state % set_real_vector(q)
      call rows % apply(unknowns, [state], out)
      call out % real_vector(residual)
      direction = stored_field('v', unknowns % vertex_set(), num_unknowns)
      call direction % set_real_vector(q)
      call rows % partial_action(unknowns, [state], &
           & [variation(rows % argument(1), direction)], out)
      call out % real_vector(acted)
      physics = van_der_pol(2)
      call governing_rows(physics, q, governing)
      call show_block(governing, residual, maxval(abs(residual - acted)))
    end subroutine derived_rows
    subroutine scheme_reach(tails, heads, source_degree, determines)
      integer, allocatable, intent(out) :: tails(:), heads(:)
      integer, allocatable, intent(out) :: source_degree(:), determines(:)
      integer :: j, k
      tails = [((k - j, j = 0, order), k = order + 1, num_instants), &
           &   ((k - j, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
      heads = [((k, j = 0, order), k = order + 1, num_instants), &
           &   ((k, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
      determines = [((1, j = 0, order), k = order + 1, num_instants), &
           &        ((2, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
      allocate(source_degree(size(tails)), source=0)
    end subroutine scheme_reach
    subroutine show_block(governing, residual, jacobian_gap)
      real(dp), intent(in) :: governing(:), residual(:), jacobian_gap
      integer :: k
      write(*,'(a)') ' the block on a state sampled from t squared'
      write(*,'(a)') '   instant   governing    velocity   acceleration'
      do k = 1, num_instants
         write(*,'(i10,3es13.2)') k, governing(k), &
              & residual(unknown(k, 1)), residual(unknown(k, 2))
      end do
      write(*,'(a)')        ' '
      write(*,'(a)')        ' the derived rows are linear, so the stencil is its own jacobian:'
      write(*,'(a,es11.2)') '   largest difference between apply and partial action ', &
           & jacobian_gap
    end subroutine show_block
    subroutine governing_rows(physics, q, r)
      type(expression), intent(in) :: physics
      real(dp)         , intent(in) :: q(:)
      real(dp), allocatable, intent(out) :: r(:)
      type(stored_directed_graph) :: instants
      type(stored_field) :: state, design
      class(field), allocatable :: out
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', instants % vertex_set(), num_unknowns)
      design   = stored_field('nu', instants % vertex_set(), num_instants)
      call state  % set_real_vector(q)
      call design % set_real_vector(spread(nu, 1, num_instants))
      call physics % apply(instants, [state, design], out)
      call out % real_vector(r)
    end subroutine governing_rows
    subroutine physics_partials(degree, q0, q_top, design)
      integer , intent(in) :: degree
      real(dp), intent(in) :: q0(:), q_top(:), design(:)
      real(dp), parameter :: delta = 1.0e-6_dp
      integer , parameter :: instants = 2
      type(expression) :: physics
      type(stored_directed_graph) :: graph_of
      type(stored_field) :: state, nu_field, direction
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:)
      real(dp) :: q(instants * (degree + 1)), v(instants * (degree + 1))
      real(dp) :: closed(0:degree), taken(0:degree), differenced(0:degree)
      integer :: d, k, nd
      nd      = degree + 1
      physics = van_der_pol(degree)
      call sample(degree, q0, q_top, design, instants, q)
      graph_of = stored_directed_graph(instants, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', graph_of % vertex_set(), size(q))
      nu_field = stored_field('nu', graph_of % vertex_set(), instants)
      direction = stored_field('v', graph_of % vertex_set(), size(q))
      call state    % set_real_vector(q)
      call nu_field % set_real_vector(spread(nu, 1, instants))
      do d = 0, degree
         v = 0.0_dp
         do k = 1, instants
            v((k - 1) * nd + d + 1) = 1.0_dp
         end do
         call direction % set_real_vector(v)
         call physics % partial_action(graph_of, [state, nu_field], &
              & [variation(physics % argument(1), direction)], out)
         call out % real_vector(exact)
         taken(d) = exact(1)
         differenced(d) = state_difference(physics, graph_of, state, nu_field, q, v)
         closed(d)      = closed_form(degree, d, q(1:nd))
      end do
      call show_row(degree, taken, closed, differenced)
      call design_partial(physics, graph_of, state, nu_field, q, instants, nd)
    end subroutine physics_partials
    subroutine show_row(degree, taken, closed, differenced)
      integer , intent(in) :: degree
      real(dp), intent(in) :: taken(0:), closed(0:), differenced(0:)
      integer :: d
      write(*,'(a)')         ' '
      write(*,'(a,i0)')      ' van der pol at degree ', degree
      write(*,'(a,9i12)')    '   partial in degree      ', [(d, d = 0, degree)]
      write(*,'(a,9f12.6)')  '   from partial_action    ', taken
      write(*,'(a,9f12.6)')  '   closed form            ', closed
      write(*,'(a,9f12.6)')  '   central difference     ', differenced
    end subroutine show_row
    subroutine sample(degree, q0, q_top, design, instants, q)
      integer , intent(in)  :: degree, instants
      real(dp), intent(in)  :: q0(:), q_top(:), design(:)
      real(dp), intent(out) :: q(:)
      integer :: k, nd
      nd = degree + 1
      q  = 0.0_dp
      do k = 1, instants
         q((k - 1) * nd + 1)  = q0(k)
         q((k - 1) * nd + nd) = q_top(k)
         if (degree >= 2) q((k - 1) * nd + nd - 1) = design(k)
      end do
    end subroutine sample
    function state_difference(physics, graph_of, state, nu_field, q, v) result(d)
      type(expression)          , intent(in)    :: physics
      type(stored_directed_graph), intent(in)    :: graph_of
      type(stored_field)         , intent(inout) :: state
      type(stored_field)         , intent(in)    :: nu_field
      real(dp)                   , intent(in)    :: q(:), v(:)
      real(dp) :: d
      real(dp), parameter :: delta = 1.0e-6_dp
      class(field), allocatable :: out
      real(dp), allocatable :: plus(:), minus(:)
      call state % set_real_vector(q + delta * v)
      call physics % apply(graph_of, [state, nu_field], out)
      call out % real_vector(plus)
      call state % set_real_vector(q - delta * v)
      call physics % apply(graph_of, [state, nu_field], out)
      call out % real_vector(minus)
      call state % set_real_vector(q)
      d = (plus(1) - minus(1)) / (2.0_dp * delta)
    end function state_difference
    subroutine design_partial(physics, graph_of, state, nu_field, q, instants, nd)
      type(expression)          , intent(in)    :: physics
      type(stored_directed_graph), intent(in)    :: graph_of
      type(stored_field)         , intent(inout) :: state, nu_field
      real(dp)                   , intent(in)    :: q(:)
      integer                    , intent(in)    :: instants, nd
      real(dp), parameter :: delta = 1.0e-6_dp
      type(stored_field) :: direction
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:)
      real(dp) :: w(instants), q0, q_below
      w    = 1.0_dp
      q0      = q(1)
      q_below = q(nd - 1)
      direction = stored_field('w', graph_of % vertex_set(), instants)
      call direction % set_real_vector(w)
      call physics % partial_action(graph_of, [state, nu_field], &
           & [variation(physics % argument(2), direction)], out)
      call out % real_vector(exact)
      call nu_field % set_real_vector(spread(nu, 1, instants) + delta * w)
      call physics % apply(graph_of, [state, nu_field], out)
      call out % real_vector(plus)
      call nu_field % set_real_vector(spread(nu, 1, instants) - delta * w)
      call physics % apply(graph_of, [state, nu_field], out)
      call out % real_vector(minus)
      call nu_field % set_real_vector(spread(nu, 1, instants))
      write(*,'(a,f12.6)') '   partial in the design  ', exact(1)
      write(*,'(a,f12.6)') '   closed form            ', -(1.0_dp - q0 * q0) * q_below
      write(*,'(a,f12.6)') '   central difference     ', (plus(1) - minus(1)) / (2.0_dp * delta)
    end subroutine design_partial
    pure real(dp) function closed_form(degree, d, q) result(c)
      integer , intent(in) :: degree, d
      real(dp), intent(in) :: q(0:)
      if (d == degree) then
         c = 1.0_dp
      else if (d == degree - 1) then
         c = -nu * (1.0_dp - q(0) * q(0))
      else if (d == 0) then
         c = 2.0_dp * nu * q(0) * q(degree - 1) + 1.0_dp
      else
         c = 0.0_dp
      end if
    end function closed_form
    subroutine edge_weights(dt, tails, heads, source_degree, determines, w)
      real(dp), intent(in) :: dt(:)
      integer , intent(in) :: tails(:), heads(:), source_degree(:), determines(:)
      real(dp), allocatable, intent(out) :: w(:)
      call weights_of(scheme_weight(bdf_family(order)), num_instants, tails, heads, dt, &
           & source_degree, determines, w)
    end subroutine edge_weights
  end subroutine demo_constraint_rows
  subroutine demo_coupling_relation()
    implicit none
    integer , parameter :: order = 2
    integer , parameter :: num_instants = 5
    integer , parameter :: num_conditions = 2
    type(level_storage)      :: store
    type(relational_binding) :: binding
    type(set_map)            :: sets
    type(csr_relation)       :: reach
    integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: weight(:), dt(:)
    integer :: slices(num_instants), source_carrier, target_carrier
    integer :: relation_element, coupling, block
    integer :: j, k
    tails = [((k - j, j = 0, order), k = order + 1, num_instants), &
         &   ((k - j, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    heads = [((k, j = 0, order), k = order + 1, num_instants), &
         &   ((k, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    determines = [((1, j = 0, order), k = order + 1, num_instants), &
         &        ((2, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    allocate(source_degree(size(tails)), source=0)
    do k = 1, num_instants
       slices(k) = store % assemble([integer ::], 0)
    end do
    source_carrier   = store % assemble([integer ::], 0)
    target_carrier   = store % assemble([integer ::], 0)
    relation_element = store % assemble([integer ::], 0)
    coupling = store % couple([slices, source_carrier, target_carrier], &
         & [relation_element])
    block    = store % assemble(slices, coupling)
    call describe(source_carrier, num_instants)
    call describe(target_carrier, num_instants * num_conditions)
    allocate(table(2, size(tails)))
    table(1,:) = tails
    table(2,:) = (heads - 1) * num_conditions + determines
    reach = built_reach(table)
    do k = 1, num_instants
       call bind_carrier(slices(k))
    end do
    call bind_carrier(source_carrier)
    call bind_carrier(target_carrier)
    call bind_reach()
    write(*,'(a)')      ' the block and its coupling'
    write(*,'(a,i3)')   '   members of the block        ', level_num_members(store % node(block))
    write(*,'(a,l3)')   '   block is consistent         ', level_consistent(store % node(block))
    write(*,'(a,i3)')   '   carriers of the coupling    ', num_member_sets(store % node(coupling))
    write(*,'(a,i3)')   '   relations of the coupling   ', num_relations(store % node(coupling))
    write(*,'(a,l3)')   '   coupling is relationally valid', &
         & relational_valid(store % node(coupling), binding)
    call show_tuples()
  contains
    subroutine describe(at, n)
      integer, intent(in) :: at, n
      type(graph), pointer :: g
      g => store % node(at)
      call sets % bind(g, counted_set_representation(n))
    end subroutine describe
    function built_reach(tuples) result(r)
      integer, intent(in) :: tuples(:,:)
      type(csr_relation) :: r
      type(graph), pointer :: from, into
      from => store % node(source_carrier)
      into => store % node(target_carrier)
      r = csr_relation('scheme reach', from, into, tuples, sets)
    end function built_reach
    subroutine bind_carrier(at)
      integer, intent(in) :: at
      type(graph), pointer :: g
      g => store % node(at)
      call binding % bind_set(g, g)
    end subroutine bind_carrier
    subroutine bind_reach()
      type(graph), pointer :: g
      g => store % node(relation_element)
      call binding % bind_relation(g, reach)
    end subroutine bind_reach
    subroutine show_tuples()
      class(relation), pointer :: r
      integer, allocatable :: held(:,:)
      integer :: i, e, target_index
      dt = [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp]
      call edge_weights(dt, weight)
      r => relation_at(store % node(coupling), binding, 1)
      write(*,'(a)')    ' '
      write(*,'(a,i3)') ' tuples the relation holds     ', r % num_tuples()
      write(*,'(a)')    '   component   constraint   instant  determines      weight'
      call r % tuples(held)
      do i = 1, size(held, 2)
         do e = 1, size(tails)
            target_index = (heads(e) - 1) * num_conditions + determines(e)
            if (tails(e) == held(1, i) .and. target_index == held(2, i)) then
               write(*,'(i12,i13,i10,i12,f12.5)') held(1, i), held(2, i), &
                    & heads(e), determines(e), weight(e)
               exit
            end if
         end do
      end do
    end subroutine show_tuples
    subroutine edge_weights(steps, w)
      real(dp), intent(in) :: steps(:)
      real(dp), allocatable, intent(out) :: w(:)
      call weights_of(scheme_weight(bdf_family(order)), num_instants, tails, heads, steps, &
           & source_degree, determines, w)
    end subroutine edge_weights
  end subroutine demo_coupling_relation
  subroutine demo_expansion_check()
    implicit none
    character(len=32) :: argument
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    integer , parameter :: instants = 11
    integer , parameter :: max_order = 4
    real(dp), parameter :: duration = 2.0_dp
    real(dp), parameter :: design = 1.0_dp
    real(dp) :: delta, tau
    tau = 1.0e-12_dp
    call demo_argument(1, argument)
    if (len_trim(argument) > 0) read(argument, *) tau
    call set_stopping(tau, relative, by_rate, 100)
    delta = tau ** (1.0_dp / 3.0_dp)
    write(*,'(a,es9.2,a,es9.2,a,es9.2)') ' relative tolerance', tau, &
         & '   difference step', delta, '   expected agreement tau^(2/3)', tau ** (2.0_dp / 3.0_dp)
    call expansion_of('bdf 2', bdf_family(2))
    call expansion_of('adams-moulton 3', adams_family(3))
    call expansion_of('crouzeix two-stage', crouzeix_two_stage())
  contains
    subroutine expanded(scheme, design_value, f)
      class(family), intent(in) :: scheme
      real(dp)     , intent(in) :: design_value
      real(dp), allocatable, intent(out) :: f(:)
      type(family_holder) :: holder(1)
      type(chain_block), allocatable :: chain(:)
      type(expansion)  , allocatable, target :: tower
      real(dp), allocatable :: dt(:), t(:), held(:), table(:,:)
      real(dp) :: achieved
      call cosine_partition(scheme, degrees, duration, instants, held, dt, t)
      call set_family(holder(1), scheme)
      call march_chain(holder, [instants], van_der_pol(state_degree), degrees, uniform_grid(duration), &
           & design_value, held, chain, tower, dt, t, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(state_degree)], degrees, &
           & max_order, table)
      allocate(f(0:max_order))
      f(0:) = table(:, 1)
    end subroutine expanded
    subroutine expansion_of(title, scheme)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      real(dp), allocatable :: f(:), plus(:), minus(:)
      real(dp) :: differenced(max_order)
      integer :: m
      call expanded(scheme, design, f)
      call expanded(scheme, design + delta, plus)
      call expanded(scheme, design - delta, minus)
      do m = 1, max_order
         differenced(m) = (plus(m - 1) - minus(m - 1)) / (2.0_dp * delta)
      end do
      write(*,'(a)')          ' '
      write(*,'(a)')          ' ' // title // ', van der pol at a design of one'
      write(*,'(a,5i15)')     '   derivative order    ', [(m, m = 0, max_order)]
      write(*,'(a,5f15.8)')   '   from the expansion  ', f
      write(*,'(a,15x,4f15.8)') '   differenced       ', differenced
      write(*,'(a,15x,4es15.2)') '   apart             ', abs(f(1:max_order) - differenced)
    end subroutine expansion_of
  end subroutine demo_expansion_check
  subroutine demo_family_coefficients()
    implicit none
    integer, parameter :: order = 2
    integer :: k
    call bdf_on('uniform',     [(0.5_dp, k = 1, 2 * order + 1)])
    call bdf_on('non-uniform', [0.0_dp, 0.3_dp, 0.2_dp, 0.4_dp, 0.25_dp])
    call adams_on('uniform',     3, [(0.5_dp, k = 1, 3)])
    call adams_on('non-uniform', 3, [0.0_dp, 0.3_dp, 0.2_dp])
    call dirk_on(crouzeix_two_stage())
    call bdf_step_sensitivity([0.0_dp, 0.3_dp, 0.2_dp, 0.4_dp, 0.25_dp])
  contains
    subroutine coefficients(scheme, num_vertices, tails, head, source_degree, determines, dt, c)
      class(family), intent(in)  :: scheme
      integer      , intent(in)  :: num_vertices, tails(:), head, source_degree(:), determines(:)
      real(dp)     , intent(in)  :: dt(:)
      real(dp), allocatable, intent(out) :: c(:)
      call weights_of(scheme, num_vertices, tails, [(head, k = 1, size(tails))], dt, &
           & source_degree, determines, c)
    end subroutine coefficients
    subroutine bdf_on(label, dt)
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      integer, parameter :: last = 2 * order + 1
      real(dp), allocatable :: c(:)
      real(dp) :: h0, h1
      call coefficients(bdf_family(order), last, &
           & [(last - k, k = 0, order), (last - k, k = 0, 2 * order)], last, &
           & [(0, k = 0, order), (0, k = 0, 2 * order)], &
           & [(1, k = 0, order), (2, k = 0, 2 * order)], dt, c)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2 on a ' // label // ' grid'
      write(*,'(a,3f10.5)') '   velocity     alpha_0..2      ', c(1:order + 1)
      write(*,'(a,5f10.5)') '   acceleration beta_0..4       ', c(order + 2:)
      if (label == 'uniform') then
         write(*,'(a,3f10.5)') '   tabulated    alpha           ', [1.5_dp, -2.0_dp, 0.5_dp]
         write(*,'(a,5f10.5)') '   convolution  beta            ', [2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp]
      else
         h0 = dt(last)
         h1 = dt(last - 1)
         write(*,'(a,3f10.5)') '   set_bdf row  alpha           ', &
              & [(2.0_dp * h0 + h1) / (h0 + h1), -(h0 + h1) / h1, h0 * h0 / (h1 * (h0 + h1))]
      end if
    end subroutine bdf_on
    subroutine adams_on(label, p, dt)
      character(len=*), intent(in) :: label
      integer         , intent(in) :: p
      real(dp)        , intent(in) :: dt(:)
      real(dp), allocatable :: c(:)
      call coefficients(adams_family(p), p, [(p - k, k = 0, p - 1)], p, &
           & [(2, k = 0, p - 1)], [(1, k = 0, p - 1)], dt, c)
      write(*,'(a)') ' '
      write(*,'(a)') ' adams-moulton 3 on a ' // label // ' grid'
      write(*,'(a,3f10.5)') '   quadrature   alpha_0..2      ', c
      if (label == 'uniform') then
         write(*,'(a,3f10.5)') '   tabulated    alpha           ', [5.0_dp, 8.0_dp, -1.0_dp] / 12.0_dp
      end if
    end subroutine adams_on
    subroutine dirk_on(scheme)
      type(dirk_family), intent(in) :: scheme
      real(dp), allocatable :: c(:)
      integer :: s
      s = scheme % num_stages()
      call coefficients(scheme, 2 + s, [2, 2, 3], 3, [2, 2, 2], [1, 1, 1], &
           & [(0.5_dp, k = 1, 2 + s)], c)
      write(*,'(a)') ' '
      write(*,'(a)') ' crouzeix two-stage, stage 2 from stages 1, 1, 2'
      write(*,'(a,3f10.5)') '   a_21, a_21, a_22                ', c
      call coefficients(scheme, 2 + s, [2, 3], 2 + s, [2, 2], [2, 2], &
           & [(0.5_dp, k = 1, 2 + s)], c)
      write(*,'(a,2f10.5)') '   b_1, b_2 into the arriving instant', c
      write(*,'(a,3f10.5)') '   tableau gamma, 1 - 2 gamma, b   ', &
           & (3.0_dp + sqrt(3.0_dp)) / 6.0_dp, 1.0_dp - (3.0_dp + sqrt(3.0_dp)) / 3.0_dp, 0.5_dp
    end subroutine dirk_on
    subroutine bdf_step_sensitivity(dt)
      real(dp), intent(in) :: dt(:)
      integer , parameter :: last = 2 * order + 1
      real(dp), parameter :: delta = 1.0e-6_dp
      type(stored_directed_graph) :: coupling
      type(stored_field), allocatable :: inputs(:)
      type(stored_field) :: direction
      type(bdf_family) :: scheme
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)
      scheme = bdf_family(order)
      call coupling_inputs(last, [(last - k, k = 0, order)], [(last, k = 0, order)], dt, &
           & [(0, k = 0, order)], [(1, k = 0, order)], coupling, inputs)
      direction = stored_field('v', coupling % vertex_set(), last)
      allocate(v(last), source=0.0_dp)
      v(last) = 1.0_dp
      call direction % set_real_vector(v)
      call scheme % partial_action(coupling, inputs, [variation(scheme % argument(1), direction)], out)
      call out % real_vector(exact)
      call inputs(1) % set_real_vector(dt + delta * v)
      call scheme % apply(coupling, inputs, out)
      call out % real_vector(plus)
      call inputs(1) % set_real_vector(dt - delta * v)
      call scheme % apply(coupling, inputs, out)
      call out % real_vector(minus)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2, partial of alpha_0..2 in the last step, non-uniform grid'
      write(*,'(a,3f12.6)') '   partial_action, degree one    ', exact
      write(*,'(a,3f12.6)') '   central difference            ', (plus - minus) / (2.0_dp * delta)
      call bdf_second_partials(scheme, coupling, inputs(1), inputs(2), inputs(3), dt, v)
    end subroutine bdf_step_sensitivity
    subroutine bdf_second_partials(scheme, coupling, steps, degrees, conditions, dt, v)
      type(bdf_family)           , intent(in)    :: scheme
      type(stored_directed_graph), intent(in)    :: coupling
      type(stored_field)         , intent(inout) :: steps
      type(stored_field)         , intent(in)    :: degrees, conditions
      real(dp)                   , intent(in)    :: dt(:), v(:)
      real(dp), parameter :: delta = 1.0e-4_dp
      type(stored_field) :: along_v, along_w
      class(field), allocatable :: out
      real(dp), allocatable :: plus(:), at(:), minus(:), w(:)
      real(dp), allocatable :: second(:), mixed_partial(:)
      along_v = stored_field('v', coupling % vertex_set(), size(dt))
      along_w = stored_field('w', coupling % vertex_set(), size(dt))
      call along_v % set_real_vector(v)
      w = 0.0_dp * v
      w(size(dt) - 1) = 1.0_dp
      call along_w % set_real_vector(w)
      call steps % set_real_vector(dt)
      call scheme % partial_action(coupling, [steps, degrees, conditions], &
           & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_v)], out)
      call out % real_vector(second)
      call scheme % partial_action(coupling, [steps, degrees, conditions], &
           & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_w)], out)
      call out % real_vector(mixed_partial)
      call steps % set_real_vector(dt + delta * v)
      call scheme % apply(coupling, [steps, degrees, conditions], out)
      call out % real_vector(plus)
      call steps % set_real_vector(dt)
      call scheme % apply(coupling, [steps, degrees, conditions], out)
      call out % real_vector(at)
      call steps % set_real_vector(dt - delta * v)
      call scheme % apply(coupling, [steps, degrees, conditions], out)
      call out % real_vector(minus)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2, second partials of alpha_0..2, non-uniform grid'
      write(*,'(a,3f12.6)') '   partial_action, (last, last)  ', second
      write(*,'(a,3f12.6)') '   second central difference     ', (plus - 2.0_dp * at + minus) / delta**2
      call steps % set_real_vector(dt + delta * w)
      call scheme % partial_action(coupling, [steps, degrees, conditions], &
           & [variation(scheme % argument(1), along_v)], out)
      call out % real_vector(plus)
      call steps % set_real_vector(dt - delta * w)
      call scheme % partial_action(coupling, [steps, degrees, conditions], &
           & [variation(scheme % argument(1), along_v)], out)
      call out % real_vector(minus)
      write(*,'(a,3f12.6)') '   partial_action, (last, before)', mixed_partial
      write(*,'(a,3f12.6)') '   difference of degree one      ', (plus - minus) / (2.0_dp * delta)
    end subroutine bdf_second_partials
  end subroutine demo_family_coefficients
  subroutine demo_function_identities()
    use util_derivative_terms, only : derivative_terms, integer_power, mixed_partial, coefficient, &
         & operator(+), operator(-), operator(*), operator(/), operator(**), &
         & sin, cos, exp, log, sqrt
    implicit none
    integer :: n, failures
    failures = 0
    do n = 1, 5
       call identities(n, failures)
    end do
    if (failures > 0) then
       write(*,'(a,i0,a)') ' FAIL : ', failures, ' identities exceeded the floor'
       error stop
    end if
    write(*,'(a)') ' PASS : the elementary functions compose exactly to five directions'
  contains
    subroutine identities(n, failures)
      integer, intent(in)    :: n
      integer, intent(inout) :: failures
      type(derivative_terms) :: a, b, one
      real(dp) :: v(n), w(n), closed
      integer  :: i
      do i = 1, n
         v(i) = 0.3_dp + 0.1_dp * real(i, dp)
         w(i) = 0.7_dp - 0.1_dp * real(i, dp)
      end do
      a = seeded(1.2_dp, v)
      a = integer_power(a, n) + seeded(0.5_dp, w)
      b = seeded(0.8_dp, w)
      b = integer_power(b, n)
      one = derivative_terms(1.0_dp, n)
      call held(n, 'exp(a+b) = exp(a) exp(b)', exp(a + b), exp(a) * exp(b), exp(a + b), failures)
      call held(n, 'sin^2 + cos^2 = 1', sin(a) * sin(a) + cos(a) * cos(a), one, sin(a) * sin(a), failures)
      call held(n, 'sin(a+b) addition', sin(a + b), sin(a) * cos(b) + cos(a) * sin(b), sin(a) * cos(b), failures)
      call held(n, 'log(exp(a)) = a', log(exp(a)), a, exp(a), failures)
      call held(n, 'exp(log(a)) = a', exp(log(a)), a, a, failures)
      call held(n, 'sqrt(a) sqrt(a) = a', sqrt(a) * sqrt(a), a, a, failures)
      call held(n, 'a**0.5 = sqrt(a)', a ** 0.5_dp, sqrt(a), sqrt(a), failures)
      call held(n, 'a**3.0 = a a a', a ** 3.0_dp, a * a * a, a * a * a, failures)
      call held(n, 'a**(-1.0) = 1/a', a ** (-1.0_dp), one / a, one / a, failures)
      a = seeded(0.4_dp, v)
      closed = exp(0.4_dp) * product(v)
      call held_scalar(n, 'full partial of exp', mixed_partial(exp(a)), closed, failures)
      closed = sin(0.4_dp + real(n, dp) * acos(-1.0_dp) / 2.0_dp) * product(v)
      call held_scalar(n, 'full partial of sin', mixed_partial(sin(a)), closed, failures)
    end subroutine identities
    function seeded(x, v) result(a)
      real(dp), intent(in) :: x, v(:)
      type(derivative_terms) :: a
      integer :: i
      a = derivative_terms(x, size(v))
      do i = 1, size(v)
         call a % set_direction(i, v(i))
      end do
    end function seeded
    subroutine held(n, label, got, reference, operand, failures)
      integer               , intent(in)    :: n
      character(len=*)      , intent(in)    :: label
      type(derivative_terms), intent(in)    :: got, reference, operand
      integer               , intent(inout) :: failures
      real(dp) :: worst, scale, floor
      integer  :: m
      worst = 0.0_dp
      scale = 0.0_dp
      do m = 0, 2**n - 1
         worst = max(worst, abs(coefficient(got, m) - coefficient(reference, m)))
         scale = max(scale, abs(coefficient(operand, m)))
      end do
      floor = real(3**n, dp) * epsilon(1.0_dp) * scale
      call reported(n, label, worst, floor, failures)
    end subroutine held
    subroutine held_scalar(n, label, got, reference, failures)
      integer         , intent(in)    :: n
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: got, reference
      integer         , intent(inout) :: failures
      call reported(n, label, abs(got - reference), &
           & real(3**n, dp) * epsilon(1.0_dp) * abs(reference), failures)
    end subroutine held_scalar
    subroutine reported(n, label, worst, floor, failures)
      integer         , intent(in)    :: n
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: worst, floor
      integer         , intent(inout) :: failures
      write(*,'(a,i0,a,a28,a,es9.2,a,es9.2)') '   n = ', n, '  ', label, &
           & '  difference ', worst, '  floor ', floor
      if (worst > floor) failures = failures + 1
    end subroutine reported
  end subroutine demo_function_identities
  subroutine demo_grid_design_check()
    implicit none
    integer , parameter :: state_degree = 2, degrees = state_degree + 1
    integer , parameter :: instants = 21, checked(3) = [1, 7, 20]
    real(dp), parameter :: duration = 4.0_dp, design = 0.8_dp
    type(family_holder)     :: schemes(2)
    type(expression)       :: functionals(2)
    type(chain_block) , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    integer, allocatable :: marks(:)
    real(dp), allocatable :: p(:), dt(:), t(:), v(:,:), f(:,:), tangent(:,:), adjoint(:,:)
    real(dp), allocatable :: plus(:,:), minus(:,:), q0(:), table(:,:), entries(:,:,:)
    real(dp), allocatable :: below(:,:), above(:,:), by_class(:)
    real(dp) :: tau, delta, achieved, worst
    character(len=32) :: argument
    integer :: k, j, i, route, order, max_order, nd
    tau       = 1.0e-12_dp
    max_order = 3
    call demo_argument(1, argument)
    if (len_trim(argument) > 0) read(argument, *) tau
    call demo_argument(2, argument)
    if (len_trim(argument) > 0) read(argument, *) max_order
    call set_stopping(tau, relative, by_rate, 100)
    delta = tau ** (1.0_dp / 3.0_dp)
    schemes = [held_family(bdf_family(3)), held_family(adams_family(3))]
    functionals(1) = van_der_pol_energy(state_degree)
    functionals(2) = van_der_pol_dissipation(state_degree)
    p  = [(1.0_dp + 0.5_dp * sin(real(k, dp)), k = 1, instants - 1)]
    q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design)
    call marched(p, design, f)
    call tower % step_partials(v)
    call chain_stamps(chain, tower, functionals, degrees, marks)
    call chain_derivative(chain, tower, marks, functionals, degrees, 1, forward_route, tangent)
    call chain_derivative(chain, tower, marks, functionals, degrees, 1, reverse_route, adjoint)
    route   = route_of(size(tangent, 2), size(tangent, 1), 1)
    write(*,'(a,es9.2,a,es9.2,a,es9.2)') ' relative tolerance', tau, '   difference step', delta, &
         & '   expected agreement tau^(2/3)', tau ** (2.0_dp / 3.0_dp)
    write(*,'(a,i0,a,i0,a,a)') ' designs ', size(tangent, 2), '   functionals ', size(tangent, 1), &
         & '   the gate chooses the ', trim(merge('forward', 'reverse', route == forward_route))
    write(*,'(a,es10.2)') ' tangent against adjoint over the table, relative  ', &
         & maxval(abs(tangent - adjoint)) / maxval(abs(tangent))
    write(*,'(a,es10.2,a,es10.2)') ' physics column against the expansion, relative   ', &
         & abs(adjoint(1, 1) - f(1, 1)) / abs(f(1, 1)), '  ', abs(adjoint(2, 1) - f(1, 2)) / abs(f(1, 2))
    do i = 1, 2
       write(*,'(a,i0,a,es10.2,a,es10.2,a,es10.2)') ' homogeneity, functional ', i, &
            & ':  p . df/dp / |p||df/dp|  ', &
            & dot_product(p, adjoint(i, 2:)) / (norm2(p) * norm2(adjoint(i, 2:))), &
            & '   f ', f(0, i), '   |df/dp| ', norm2(adjoint(i, 2:))
    end do
    worst = 0.0_dp
    do k = 1, size(checked)
       j = checked(k)
       call marched(p + delta * unit(j), design, plus)
       call marched(p - delta * unit(j), design, minus)
       do i = 1, 2
          worst = max(worst, abs((plus(0, i) - minus(0, i)) / (2.0_dp * delta) - adjoint(i, 1 + j)) &
               & / max(1.0_dp, abs(adjoint(i, 1 + j))))
       end do
    end do
    write(*,'(a,es10.2)') ' differenced in three weights against the route, worst   ', worst
    call marched(p, design + delta, plus)
    call marched(p, design - delta, minus)
    write(*,'(a,es10.2)') ' differenced in the parameter against the route, worst   ', &
         & maxval(abs((plus(0, :) - minus(0, :)) / (2.0_dp * delta) - adjoint(:, 1)) / &
         &        max(1.0_dp, abs(adjoint(:, 1))))
    nd = size(tangent, 2)
    do order = 2, max_order
       call marched(p, design, f, order)
       call chain_stamps(chain, tower, functionals, degrees, marks)
       call chain_derivative(chain, tower, marks, functionals, degrees, order, reverse_route, &
            & table, entries=entries)
       write(*,'(a)') ' '
       write(*,'(a,i0,a,i0,a,i0,a,i0)') ' derivatives of order ', order, ' by the reverse route: ', &
            & size(table, 1), ' tables of ', size(table, 2), ' multisets over ', nd
       do i = 1, 2
          write(*,'(a,i0,a,es10.2,a,es10.2)') ' functional ', i, &
               & ':  departure among the entries of a multiset, relative ', &
               & asymmetry(entries(i:i, :, :), nd, order), &
               & '   parameter entry against the expansion ', &
               & abs(table(i, 1) - f(order, i)) / abs(f(order, i))
       end do
       allocate(by_class(0:order), source=0.0_dp)
       do k = 1, size(checked)
          j = checked(k)
          call differenced_table(p + delta * unit(j), design, order - 1, above)
          call differenced_table(p - delta * unit(j), design, order - 1, below)
          call classed((above - below) / (2.0_dp * delta), 1 + j)
       end do
       call differenced_table(p, design + delta, order - 1, above)
       call differenced_table(p, design - delta, order - 1, below)
       call classed((above - below) / (2.0_dp * delta), 1)
       write(*,'(a,i0,a,*(es10.2))') ' differenced tables of order ', order - 1, &
            & ' against the entries, worst by parameter count from ', &
            & by_class(order:0:-1) / max(1.0_dp, maxval(abs(table)))
       deallocate(by_class)
    end do
  contains
    subroutine marched(weights, nu, f, order)
      real(dp), intent(in) :: weights(:), nu
      real(dp), allocatable, intent(out) :: f(:,:)
      integer , intent(in), optional :: order
      type(imbalance) :: left
      integer :: m
      m = 1
      if (present(order)) m = order
      call march_chain(schemes, [11, 10], van_der_pol(state_degree), degrees, &
           & designed_grid(duration), nu, q0, chain, tower, dt, t, achieved, grid_design=weights, &
           & left=left, startup=4)
      if (.not. left % converged) error stop 'grid_design_check: the march converged'
      call chain_expansion(chain, tower, functionals, degrees, m, f)
    end subroutine marched
    subroutine differenced_table(weights, nu, order, t)
      real(dp), intent(in) :: weights(:), nu
      integer , intent(in) :: order
      real(dp), allocatable, intent(out) :: t(:,:)
      real(dp), allocatable :: f(:,:)
      call marched(weights, nu, f)
      call chain_stamps(chain, tower, functionals, degrees, marks)
      call chain_derivative(chain, tower, marks, functionals, degrees, order, reverse_route, t)
    end subroutine differenced_table
    subroutine classed(e, l)
      real(dp), intent(in) :: e(:,:)
      integer , intent(in) :: l
      integer, allocatable :: s(:), with(:)
      integer :: rank, c
      do rank = 1, size(e, 2)
         s    = multiset_of(rank, order - 1, nd)
         with = [s(1:count(s < l)), l, s(count(s < l) + 1:)]
         c    = count(with == 1)
         by_class(c) = max(by_class(c), maxval(abs(e(:, rank) - table(:, multiset_rank(with, nd)))))
      end do
    end subroutine classed
    pure function unit(j) result(e)
      integer, intent(in) :: j
      real(dp) :: e(instants - 1)
      e    = 0.0_dp
      e(j) = 1.0_dp
    end function unit
  end subroutine demo_grid_design_check
  subroutine demo_jacobian_shape()
    implicit none
    write(*,'(a)') ' '
    write(*,'(a)') '  scheme        unknowns    filled   below   above   per cent full' // &
         & '     largest    on diagonal    largest row'
    call shape_of('bdf 1',   bdf_family(1),        3, 21)
    call shape_of('bdf 2',   bdf_family(2),        3, 21)
    call shape_of('bdf 3',   bdf_family(3),        3, 21)
    call shape_of('adams 2', adams_family(2),      3, 21)
    call shape_of('adams 3', adams_family(3),      3, 21)
    call shape_of('dirk 2',  crouzeix_two_stage(), 3, 21)
    call shape_of('bdf 2',   bdf_family(2),        3, 61)
    call shape_of('bdf 2',   bdf_family(2),        4, 61)
  contains
    subroutine shape_of(label, scheme, degrees, instants)
      character(len=*), intent(in) :: label
      class(family)   , intent(in) :: scheme
      integer         , intent(in) :: degrees, instants
      type(family_holder), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      integer , allocatable :: added(:)
      real(dp), allocatable :: held(:), dt(:), t(:), a(:,:)
      real(dp) :: achieved, duration, design
      duration = 3.0_dp
      design   = 1.0_dp
      allocate(schemes(1))
      call set_family(schemes(1), scheme)
      added = [instants]
      call cosine_partition(scheme, degrees, duration, instants, held, dt, t)
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      call chain_stamps(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, marks)
      call dense_jacobian(chain, design, a)
      call reported(label, a)
    end subroutine shape_of
    subroutine reported(label, a)
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: a(:,:)
      integer  :: n, i, j, filled, below, above
      real(dp) :: biggest, least, on_diagonal, row_most
      n       = size(a, 1)
      biggest = maxval(abs(a))
      least   = 1.0e-12_dp * biggest
      filled = 0
      below  = 0
      above  = 0
      do j = 1, n
         do i = 1, n
            if (abs(a(i, j)) <= least) cycle
            filled = filled + 1
            below  = max(below, i - j)
            above  = max(above, j - i)
         end do
      end do
      on_diagonal = 0.0_dp
      do i = 1, n
         on_diagonal = max(on_diagonal, abs(a(i, i)))
      end do
      row_most = maxval(sum(abs(a), dim=2))
      write(*,'(a,a,i8,i10,i8,i8,f12.2,3es15.4)') '  ', label // repeat(' ', 12 - len(label)), &
           & n, filled, below, above, 100.0_dp * real(filled, dp) / real(n * n, dp), &
           & biggest, on_diagonal, row_most
    end subroutine reported
  end subroutine demo_jacobian_shape
  subroutine demo_level_maps()
    implicit none
    integer , parameter :: max_derivative_degree = 1     ! the primal and one tangent
    integer , parameter :: max_state_degree      = 2     ! q, q', q''
    integer , parameter :: num_instants          = 3     ! per block
    integer , parameter :: num_freedoms          = 1     ! an ordinary differential equation
    real(dp), parameter :: duration              = 1.15_dp
    type(level_storage) :: store
    type(value_map)     :: values
    type(set_store)     :: sets
    integer, allocatable :: sweeps(:)
    integer :: expansion, s
    sweeps = [(one_sweep(s), s = 0, max_derivative_degree)]
    expansion = store % assemble(sweeps, 0)
    call sets % name(store % node(expansion), 'expansion of the van der pol functional in nu')
    call sets % bind(store % node(expansion), counted_set_representation(1))
    call attach_known(store % node(expansion), [1.0_dp])
    write(*,'(a)') ' the tower, and what each level carries'
    call show(store % node(expansion), 0)
    write(*,'(a)')    ' '
    write(*,'(a,i0)') ' nodes owned by the storage      ', store % num_nodes()
    write(*,'(a,i0)') ' components not yet known        ', not_yet_known(store % node(expansion))
    call determine(store % node(expansion))
    write(*,'(a,i0)') ' after every block is solved     ', not_yet_known(store % node(expansion))
  contains
    subroutine attach_known(g, x)
      type(graph), intent(in) :: g
      real(dp)   , intent(in) :: x(:)
      call values % attach_unknown(g)
      call values % mark_known(g, x)
    end subroutine attach_known
    integer function one_component(degree, instant, block_index) result(at)
      integer, intent(in) :: degree, instant, block_index
      character(len=1) :: d
      at = store % assemble([integer ::], 0)
      write(d,'(i1)') degree
      call sets % name(store % node(at), 'component of degree ' // d)
      call sets % bind(store % node(at), counted_set_representation(num_freedoms))
      if (block_index == 1 .and. instant <= 2) then
         call attach_known(store % node(at), spread(0.0_dp, 1, num_freedoms))
      else
         call values % attach_unknown(store % node(at))
      end if
    end function one_component
    integer function one_slice(instant, block_index) result(at)
      integer, intent(in) :: instant, block_index
      character(len=2) :: k
      integer :: d
      at = store % assemble([(one_component(d, instant, block_index), &
           & d = 0, max_state_degree)], 0)
      write(k,'(i2)') instant
      call sets % name(store % node(at), 'slice at instant' // k)
    end function one_slice
    integer function one_block(block_index, family) result(at)
      integer         , intent(in) :: block_index
      character(len=*), intent(in) :: family
      real(dp) :: steps(num_instants)
      integer  :: k
      at = store % assemble([(one_slice(k, block_index), k = 1, num_instants)], 0)
      call sets % name(store % node(at), family)
      steps    = duration / real(2 * num_instants, dp)
      steps(1) = 0.0_dp
      call attach_known(store % node(at), steps)
    end function one_block
    integer function one_horizon() result(at)
      character(len=8) :: t
      at = store % assemble([one_block(1, 'bdf of order 2'), &
           &                 one_block(2, 'adams-moulton of order 3')], 0)
      write(t,'(f8.4)') duration
      call sets % name(store % node(at), 'horizon of duration' // t)
    end function one_horizon
    integer function one_sweep(sensitivity) result(at)
      integer, intent(in) :: sensitivity
      character(len=1) :: s
      at = store % assemble([one_horizon()], 0)
      write(s,'(i1)') sensitivity
      if (sensitivity == 0) then
         call sets % name(store % node(at), 'sweep 0, the functional itself')
      else
         call sets % name(store % node(at), 'sweep ' // s // ', derivative ' // s // ' in nu')
      end if
      call values % attach_unknown(store % node(at))
    end function one_sweep
    recursive subroutine show(g, depth)
      type(graph), intent(in) :: g
      integer    , intent(in) :: depth
      character(len=:), allocatable :: name
      name = ' '
      if (sets % labelled(g)) name = sets % label_of(g)
      write(*,'(a,a,a,a)') repeat('   ', depth + 1), name, &
           & '   [' // status_name(values % status_of(g)) // ']', extent_of(g)
      if (level_is_leaf(g)) return
      call show_each(level_members(g), depth + 1)
    end subroutine show
    recursive subroutine show_each(members, depth)
      type(branch), intent(in) :: members
      integer     , intent(in) :: depth
      type(graph), pointer :: first
      if (sequence_empty(members)) return
      first => sequence_first(members)
      call show(first, depth)
      call show_each(sequence_rest(members), depth)
    end subroutine show_each
    pure function status_name(status) result(name)
      integer, intent(in) :: status
      character(len=:), allocatable :: name
      select case (status)
      case (VALUE_KNOWN)
         name = 'known'
      case (VALUE_UNKNOWN)
         name = 'not yet known'
      case (VALUE_UNATTACHED)
         name = 'no value'
      case default
         error stop 'level_maps: a value status is one of the three'
      end select
    end function status_name
    function extent_of(g) result(text)
      type(graph), intent(in) :: g
      character(len=:), allocatable :: text
      character(len=3) :: n
      text = ''
      if (.not. sets % describes(g)) return
      write(n,'(i3)') sets % num_members_of(g)
      text = '   extent' // n
    end function extent_of
    recursive subroutine determine(g)
      type(graph), intent(in) :: g
      if (level_is_leaf(g)) then
         if (values % status_of(g) == VALUE_UNKNOWN) then
            call values % mark_known(g, spread(1.0_dp, 1, num_freedoms))
         end if
         return
      end if
      call determine_each(level_members(g))
    end subroutine determine
    recursive subroutine determine_each(members)
      type(branch), intent(in) :: members
      type(graph), pointer :: first
      if (sequence_empty(members)) return
      first => sequence_first(members)
      call determine(first)
      call determine_each(sequence_rest(members))
    end subroutine determine_each
    recursive integer function not_yet_known(g) result(n)
      type(graph), intent(in) :: g
      n = 0
      if (level_is_leaf(g)) then
         if (values % status_of(g) == VALUE_UNKNOWN) n = 1
         return
      end if
      n = counted(level_members(g))
    end function not_yet_known
    recursive integer function counted(members) result(n)
      type(branch), intent(in) :: members
      type(graph), pointer :: first
      n = 0
      if (sequence_empty(members)) return
      first => sequence_first(members)
      n = not_yet_known(first) + counted(sequence_rest(members))
    end function counted
  end subroutine demo_level_maps
  subroutine demo_level_shape()
    implicit none
    integer, parameter :: max_instants     = 2
    integer, parameter :: max_stages       = 2
    integer, parameter :: max_state_degree = 1
    type(level_storage) :: store
    type(graph), pointer :: root
    integer :: multistep, multistage
    multistep  = one_block(one_multistep_slice)
    multistage = one_block(one_multistage_slice)
    write(*,'(a)') ' a multistep block: slices hold components'
    root => store % node(multistep)
    call show(root, 1)
    write(*,'(a)') ' '
    write(*,'(a)') ' a multistage block: slices hold stages'
    root => store % node(multistage)
    call show(root, 1)
    write(*,'(a)')    ' '
    write(*,'(a,i0)') ' nodes owned by the storage   ', store % num_nodes()
    call the_stranger_is_refused()
  contains
    recursive function members_of(n, make) result(members)
      integer, intent(in) :: n
      interface
         integer function make()
         end function make
      end interface
      integer, allocatable :: members(:)
      integer :: first
      if (n == 0) then
         allocate(members(0))
         return
      end if
      first   = make()
      members = [first, members_of(n - 1, make)]
    end function members_of
    integer function one_leaf() result(leaf)
      leaf = store % assemble([integer ::], 0)
    end function one_leaf
    integer function coupled(members) result(level)
      integer, intent(in) :: members(:)
      level = store % assemble(members, store % assemble(members, 0))
    end function coupled
    function components() result(members)
      integer, allocatable :: members(:)
      members = members_of(max_state_degree + 1, one_leaf)
    end function components
    integer function one_multistep_slice() result(slice)
      slice = coupled(components())
    end function one_multistep_slice
    integer function one_stage() result(stage)
      stage = coupled(components())
    end function one_stage
    integer function one_multistage_slice() result(slice)
      slice = coupled(members_of(max_stages, one_stage))
    end function one_multistage_slice
    integer function one_block(make) result(block)
      interface
         integer function make()
         end function make
      end interface
      block = coupled(members_of(max_instants + 1, make))
    end function one_block
    recursive subroutine show(g, depth)
      type(graph), intent(in) :: g
      integer    , intent(in) :: depth
      if (level_is_leaf(g)) then
         write(*,'(a,a)') repeat('   ', depth), 'leaf'
         return
      end if
      write(*,'(a,a,i0,a,l1,a,l1)') repeat('   ', depth), 'members ', &
           & level_num_members(g), '   couples ', level_couples(g), &
           & '   consistent ', level_consistent(g)
      call show_each(level_members(g), depth + 1)
    end subroutine show
    recursive subroutine show_each(members, depth)
      type(branch), intent(in) :: members
      integer     , intent(in) :: depth
      type(graph), pointer :: first
      if (sequence_empty(members)) return
      first => sequence_first(members)
      call show(first, depth)
      call show_each(sequence_rest(members), depth)
    end subroutine show_each
    subroutine the_stranger_is_refused()
      type(graph), pointer :: subject_node, other
      integer, allocatable :: mine(:)
      integer :: stranger, alien, subject
      mine     = members_of(2, one_leaf)
      stranger = one_leaf()
      alien    = store % assemble([mine(1), stranger], 0)
      subject = coupled(mine)
      subject_node => store % node(subject)
      write(*,'(a)')    ' '
      write(*,'(a,l1)') ' coupling beginning with its own members is consistent ', &
           & level_consistent(subject_node)
      other => store % node(alien)
      subject_node % branch(2) = known_branch(other)
      write(*,'(a,i0)') ' the other coupling carries the same count             ', &
           & level_num_members(subject_node)
      write(*,'(a,l1)') ' and is refused by identity                            ', &
           & level_consistent(subject_node)
    end subroutine the_stranger_is_refused
  end subroutine demo_level_shape
  subroutine demo_marched_block()
    implicit none
    integer , parameter :: max_state_degree = 2
    integer , parameter :: degrees = max_state_degree + 1
    real(dp), parameter :: duration = 2.0_dp
    call order_of('bdf 2',           bdf_family(2),   2.0_dp)
    call order_of('bdf 3',           bdf_family(3),   3.0_dp)
    call order_of('adams-moulton 3', adams_family(3), 3.0_dp)
    call nonlinear()
  contains
    pure integer function unknown(instant, degree)
      integer, intent(in) :: instant, degree
      unknown = (instant - 1) * degrees + degree + 1
    end function unknown
    subroutine march(scheme, n, design_value, q, t, achieved)
      class(family), intent(in)  :: scheme
      integer      , intent(in)  :: n
      real(dp)     , intent(in)  :: design_value
      real(dp), allocatable, intent(out) :: q(:), t(:)
      real(dp)     , intent(out) :: achieved
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_holder) :: holder(1)
      integer, allocatable :: at(:)
      real(dp), allocatable :: dt(:), held(:)
      call cosine_partition(scheme, degrees, duration, n, held, dt, t)
      call set_family(holder(1), scheme)
      call tower % build(van_der_pol(max_state_degree), holder, [n], uniform_grid(duration), 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(max_state_degree), held, rows, at)
      call solved(rows, design_value, q, achieved)
    end subroutine march
    pure real(dp) function worst(q, t) result(e)
      real(dp), intent(in) :: q(:), t(:)
      integer :: k
      e = 0.0_dp
      do k = 1, size(t)
         e = max(e, abs(q(unknown(k, 0)) - cosine(0, t(k))))
      end do
    end function worst
    subroutine order_of(title, scheme, expected)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      real(dp)        , intent(in) :: expected
      real(dp), allocatable :: q(:), t(:)
      real(dp) :: e(4), achieved
      integer :: level, steps
      do level = 1, 4
         steps = 10 * 2 ** (level - 1)
         call march(scheme, steps + 1, 0.0_dp, q, t, achieved)
         e(level) = worst(q, t)
      end do
      write(*,'(a)')          ' '
      write(*,'(a)')          ' ' // title // ' on the harmonic oscillator'
      write(*,'(a,4i11)')     '   steps                      ', [(10 * 2 ** (level - 1), level = 1, 4)]
      write(*,'(a,4es11.3)')  '   worst error                ', e
      write(*,'(a,33x,3f11.3)') '   ratio                    ', e(1:3) / e(2:4)
      write(*,'(a,f11.3)')    '   two to the scheme order    ', 2.0_dp ** expected
      write(*,'(a,es11.3)')   '   residual newton achieved   ', achieved
    end subroutine order_of
    subroutine nonlinear()
      real(dp), allocatable :: q(:), t(:)
      real(dp) :: achieved
      call march(bdf_family(2), 41, 1.0_dp, q, t, achieved)
      write(*,'(a)')        ' '
      write(*,'(a)')        ' van der pol at a design of one, bdf 2, 40 steps'
      write(*,'(a,es11.3)') '   residual newton achieved   ', achieved
      write(*,'(a,3f11.5)') '   the last instant, q q'' q"  ', &
           & q(unknown(size(t), 0)), q(unknown(size(t), 1)), q(unknown(size(t), 2))
    end subroutine nonlinear
  end subroutine demo_marched_block
  subroutine demo_marched_horizon()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    real(dp), parameter :: duration = 2.0_dp
    call splitting_changes_nothing()
    call across_a_change_of_scheme()
    call sensitivity_across_the_junction()
  contains
    subroutine marched(schemes, added, q, achieved, design_value)
      type(family_holder), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp), allocatable, intent(out) :: q(:)
      real(dp)           , intent(out) :: achieved
      real(dp), intent(in), optional :: design_value
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: first(:), last(:)
      real(dp), allocatable :: dt(:), t(:), held(:)
      real(dp) :: design
      integer :: k, n
      design = 0.0_dp
      if (present(design_value)) design = design_value
      call horizon_bounds(schemes, added, degrees - 1, first, last)
      n = last(size(added))
      call cosine_partition(schemes(1) % scheme, degrees, duration, n, held, dt, t)
      call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
           & uniform_grid(duration), design, held, &
           & chain, tower, dt, t, achieved)
      allocate(q(n * degrees))
      do k = 1, n
         q((k - 1) * degrees + 1:k * degrees) = instant_components(chain, k)
      end do
    end subroutine marched
    subroutine splitting_changes_nothing()
      type(family_holder) :: whole(1), split(2)
      real(dp), allocatable :: q_whole(:), q_split(:)
      real(dp) :: achieved_whole, achieved_split
      whole = [held_family(bdf_family(2))]
      split = [held_family(bdf_family(2)), held_family(bdf_family(2))]
      call marched(whole, [40], q_whole, achieved_whole)
      call marched(split, [20, 20], q_split, achieved_split)
      write(*,'(a)')        ' bdf 2 over forty instants, in one block and in two'
      write(*,'(a,i0)')     '   unknowns, whole             ', size(q_whole)
      write(*,'(a,i0)')     '   unknowns, split             ', size(q_split)
      write(*,'(a,es11.2)') '   worst difference between them', maxval(abs(q_whole - q_split))
      write(*,'(a,es11.2)') '   residual, whole             ', achieved_whole
      write(*,'(a,es11.2)') '   residual, split             ', achieved_split
    end subroutine splitting_changes_nothing
    subroutine across_a_change_of_scheme()
      type(family_holder) :: schemes(2)
      real(dp), allocatable :: q(:), dt(:), t(:)
      real(dp) :: e(3), achieved
      integer :: level, added, k, n
      schemes = [held_family(bdf_family(2)), held_family(adams_family(3))]
      do level = 1, 3
         added = 10 * 2 ** (level - 1)
         call marched(schemes, [added, added], q, achieved)
         n = 2 * added
         call partition(duration, n, dt, t)
         e(level) = 0.0_dp
         do k = 1, n
            e(level) = max(e(level), abs(q((k - 1) * degrees + 1) - cosine(0, t(k))))
         end do
      end do
      write(*,'(a)')          ' '
      write(*,'(a)')          ' bdf 2 then adams-moulton 3, on the harmonic oscillator'
      write(*,'(a,3i11)')     '   instants                   ', [(20 * 2 ** (level - 1), level = 1, 3)]
      write(*,'(a,3es11.3)')  '   worst error                ', e
      write(*,'(a,22x,2f11.3)') '   ratio                    ', e(1:2) / e(2:3)
      write(*,'(a,es11.2)')   '   residual                   ', achieved
    end subroutine across_a_change_of_scheme
    real(dp) function energy_of(q, n, design_value) result(f)
      real(dp), intent(in) :: q(:)
      integer , intent(in) :: n
      real(dp), intent(in) :: design_value
      type(stored_directed_graph) :: unknowns, instants
      type(stored_field) :: state, knobs
      real(dp), allocatable :: dt(:), t(:)
      call partition(duration, n, dt, t)
      unknowns = unknowns_graph(n, degrees)
      instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', unknowns % vertex_set(), size(q))
      knobs    = stored_field('design', unknowns % vertex_set(), n)
      call state % set_real_vector(q)
      call knobs % set_real_vector(spread(design_value, 1, n))
      f = functional_of(van_der_pol_energy(state_degree), instants, [state, knobs], dt)
    end function energy_of
    subroutine sensitivity_across_the_junction()
      real(dp), parameter :: delta = 1.0e-6_dp
      real(dp), parameter :: design = 1.0_dp
      integer , parameter :: added(2) = [10, 10]
      type(family_holder) :: schemes(2)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      type(expression)       :: energy(1)
      real(dp), allocatable :: q(:), dt(:), table(:,:)
      real(dp) :: f, tangent, adjoint, differenced, achieved
      integer :: n
      schemes = [held_family(bdf_family(2)), held_family(adams_family(3))]
      n = sum(added)
      call marched(schemes, added, q, achieved, design)
      f = energy_of(q, n, design)
      call chained(schemes, added, design, chain, tower, dt)
      energy(1) = van_der_pol_energy(state_degree)
      call chain_stamps(chain, tower, energy, degrees, marks)
      call chain_derivative(chain, tower, marks, energy, degrees, 1, forward_route, table)
      tangent     = first_of(table)
      call chain_derivative(chain, tower, marks, energy, degrees, 1, reverse_route, table)
      adjoint     = first_of(table)
      differenced = differenced_energy(schemes, added, n, design, delta)
      write(*,'(a)')        ' '
      write(*,'(a)')        ' bdf 2 then adams-moulton 3, van der pol at a design of one'
      write(*,'(a,i0)')     '   blocks                      ', size(added)
      write(*,'(a,f16.10)') '   the functional             ', f
      write(*,'(a,f16.10)') '   sensitivity, tangent       ', tangent
      write(*,'(a,f16.10)') '   sensitivity, adjoint       ', adjoint
      write(*,'(a,f16.10)') '   sensitivity, differenced   ', differenced
      write(*,'(a,es16.2)') '   tangent against adjoint    ', abs(tangent - adjoint)
      write(*,'(a,es16.2)') '   tangent against difference ', abs(tangent - differenced)
    end subroutine sensitivity_across_the_junction
    subroutine chained(schemes, added, design, chain, tower, dt)
      type(family_holder), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp)           , intent(in) :: design
      type(chain_block), allocatable, intent(out) :: chain(:)
      type(expansion)  , allocatable, intent(inout), target :: tower
      real(dp)         , allocatable, intent(out) :: dt(:)
      integer , allocatable :: first(:), last(:)
      real(dp), allocatable :: t(:), held(:)
      real(dp) :: achieved
      call horizon_bounds(schemes, added, degrees - 1, first, last)
      call cosine_partition(schemes(1) % scheme, degrees, duration, last(size(added)), held, dt, t)
      call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
           & uniform_grid(duration), design, held, &
           & chain, tower, dt, t, achieved)
    end subroutine chained
    real(dp) function differenced_energy(schemes, added, n, design, delta) result(d)
      type(family_holder), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:), n
      real(dp)           , intent(in) :: design, delta
      real(dp), allocatable :: plus(:), minus(:)
      real(dp) :: achieved
      call marched(schemes, added, plus,  achieved, design + delta)
      call marched(schemes, added, minus, achieved, design - delta)
      d = (energy_of(plus, n, design + delta) - &
         & energy_of(minus, n, design - delta)) / (2.0_dp * delta)
    end function differenced_energy
  end subroutine demo_marched_horizon
  subroutine demo_marched_stages()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    real(dp), parameter :: duration = 2.0_dp
    call order_of('implicit midpoint',  implicit_midpoint())
    call order_of('crouzeix two-stage', crouzeix_two_stage())
    call order_of('crouzeix three-stage', crouzeix_three_stage())
  contains
    real(dp) function worst_error(scheme, n) result(e)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: n
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_holder) :: holder(1)
      real(dp), allocatable :: dt(:), t(:), q(:)
      integer , allocatable :: at(:)
      real(dp) :: achieved
      integer :: k, d
      call partition(duration, n, dt, t)
      call set_family(holder(1), scheme)
      call tower % build(van_der_pol(state_degree), holder, [n], uniform_grid(duration), 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), &
           & [(cosine(d, t(1)), d = 0, degrees - 1)], rows, at)
      call solved(rows, 0.0_dp, q, achieved)
      e = 0.0_dp
      do k = 1, n
         e = max(e, abs(q(at(k) + 1) - cosine(0, t(k))))
      end do
    end function worst_error
    subroutine order_of(title, scheme)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      real(dp) :: e(3)
      integer :: level
      do level = 1, 3
         e(level) = worst_error(scheme, 5 * 2 ** (level - 1) + 1)
      end do
      write(*,'(a)')            ' '
      write(*,'(a)')            ' ' // title // ' on the harmonic oscillator'
      write(*,'(a,i0,a)')       '   stages                     ', scheme % num_stages(), ''
      write(*,'(a,3i11)')       '   steps                      ', [(5 * 2 ** (level - 1), level = 1, 3)]
      write(*,'(a,3es11.3)')    '   worst error                ', e
      write(*,'(a,22x,2f11.3)') '   ratio                    ', e(1:2) / e(2:3)
    end subroutine order_of
  end subroutine demo_marched_stages
  subroutine demo_memory_shape()
    implicit none
    integer, parameter :: degrees = 3, order = 2
    type(stored_directed_graph) :: gr
    type(stored_field)          :: over
    type(block_residual)        :: rows
    type(expansion) :: tower
    type(family_holder) :: holder(1)
    integer, allocatable :: at(:)
    type(bdf_family)            :: scheme
    character(len=32) :: what, given
    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: dt(:), t(:), held(:)
    integer :: instants, n, m, h, band, i, j, e
    call demo_argument(1, what)
    call demo_argument(2, given)
    read(given,*) instants
    scheme = bdf_family(order)
    h      = scheme % history_depth(degrees - 1)
    n      = (instants - h) * degrees
    band   = h * degrees
    e = 0
    do i = 1, n
       do j = max(1, i - band), i
          e = e + 1
       end do
    end do
    m = e
    allocate(tails(m), heads(m))
    e = 0
    do i = 1, n
       do j = max(1, i - band), i
          e = e + 1
          tails(e) = j
          heads(e) = i
       end do
    end do
    select case (trim(what))
    case ('none')
       continue
    case ('vertices')
       gr = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    case ('edges')
       gr = stored_directed_graph(n, tails=tails, heads=heads)
    case ('field')
       gr   = stored_directed_graph(n, tails=tails, heads=heads)
       over = stored_field('x', gr % vertex_set(), n)
       call over % set_real_vector(spread(1.0_dp, 1, n))
    case ('block')
       call partition(3.0_dp, instants, dt, t)
       allocate(held(h * degrees), source=0.0_dp)
       call set_family(holder(1), scheme)
       call tower % build(van_der_pol(degrees - 1), holder, [instants], uniform_grid(3.0_dp), 0, 0.0_dp)
       call block_from(tower, 1, scheme, van_der_pol(degrees - 1), held, rows, at)
    case default
       error stop 'memory_shape: the part is none, vertices, edges, field or block'
    end select
    write(*,'(a,i8,i9,i9,f12.3)') trim(what), instants, n, m, peak()
  contains
    real(dp) function peak() result(mb)
      integer :: u, status, kb
      character(len=80) :: line
      mb = 0.0_dp
      open(newunit=u, file='/proc/self/status', action='read')
      do
         read(u,'(a)',iostat=status) line
         if (status /= 0) exit
         if (line(1:6) == 'VmHWM:') then
            read(line(7:),*) kb
            mb = real(kb, dp) / 1000.0_dp
            exit
         end if
      end do
      close(u)
    end function peak
  end subroutine demo_memory_shape
  subroutine demo_randomized_checks()
    implicit none
    integer , parameter :: max_order = 2
    integer :: seed, cases, i, failures, skipped
    character(len=32) :: argument
    call demo_argument(1, argument); read(argument,*) seed
    call demo_argument(2, argument); read(argument,*) cases
    failures = 0
    skipped  = 0
    write(*,'(a)') '  case  scheme          split          directions     routes'
    do i = 1, cases
       call one_case(seed + 7919 * i, i, failures, skipped)
    end do
    write(*,'(a)') ' '
    write(*,'(a)') '  case  chain                        orders'
    do i = 1, cases / 2
       call mixed_case(seed + 104729 * i, i, failures, skipped)
    end do
    write(*,'(a)')      ' '
    write(*,'(a,i0,a,i0,a,i0,a)') ' ', cases - skipped, ' cases checked, ', &
         & failures, ' failed, ', skipped, ' passed over'
    if (failures > 0) error stop 'randomized_checks: an invariant did not hold'
  contains
    logical function verbose()
      character(len=8) :: argument
      integer :: count
      count = demo_argument_count()
      verbose = .false.
      if (count >= 3) then
         call demo_argument(3, argument)
         verbose = trim(argument) == 'verbose'
      end if
    end function verbose
    integer function drawn(state, below) result(n)
      integer(int64), intent(inout) :: state
      integer       , intent(in)    :: below
      state = mod(1103515245_int64 * state + 12345_int64, 2147483648_int64)
      n = int(mod(state / 65536_int64, int(below, int64))) + 1
    end function drawn
    real(dp) function drawn_real(state, low, high) result(x)
      integer(int64), intent(inout) :: state
      real(dp)      , intent(in)    :: low, high
      state = mod(1103515245_int64 * state + 12345_int64, 2147483648_int64)
      x = low + (high - low) * real(state, dp) / 2147483648.0_dp
    end function drawn_real
    subroutine directions_of(schemes, added, degrees, duration, design, tangent, adjoint)
      type(family_holder), intent(in)  :: schemes(:)
      integer            , intent(in)  :: added(:), degrees
      real(dp)           , intent(in)  :: duration, design
      real(dp)           , intent(out) :: tangent, adjoint
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      type(expression)       :: energy(1)
      real(dp), allocatable :: held(:), dt(:), t(:), table(:,:)
      real(dp) :: achieved
      call cosine_partition(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      energy(1) = van_der_pol_energy(degrees - 1)
      call chain_stamps(chain, tower, energy, degrees, marks)
      call chain_derivative(chain, tower, marks, energy, degrees, 1, forward_route, table)
      tangent = first_of(table)
      call chain_derivative(chain, tower, marks, energy, degrees, 1, reverse_route, table)
      adjoint = first_of(table)
    end subroutine directions_of
    subroutine one_case(from, index, failures, skipped)
      integer, intent(in)    :: from, index
      integer, intent(inout) :: failures, skipped
      type(family_holder), allocatable :: whole(:), split(:)
      real(dp), allocatable :: f_whole(:), f_split(:)
      real(dp) :: duration, design, tangent, adjoint, achieved
      integer :: degrees, order, kind, instants, half
      character(len=16) :: label
      logical :: staged
      call draw(from, degrees, order, kind, instants, duration, design)
      staged = kind == 3
      label  = named(kind, order)
      allocate(whole(1), split(2))
      call fill(whole(1), kind, order)
      call fill(split(1), kind, order)
      call fill(split(2), kind, order)
      call halved(whole(1) % scheme, degrees, instants, half)
      if (verbose()) write(*,'(a,i0,a,i0,a,i0,a,f8.4,a,f8.4)') &
           & '        degrees ', degrees, '  instants ', instants, '  half ', half, &
           & '  duration ', duration, '  design ', design
      call both_ways(whole, split, instants, half, degrees, duration, design, &
           & f_whole, f_split, achieved)
      if (achieved > 1.0e-6_dp) then
         skipped = skipped + 1
         write(*,'(i6,2x,a16,a,es9.2)') index, label, &
              & '   a march did not converge, passed over: ', achieved
         return
      end if
      call directions_of(whole, [instants], degrees, duration, design, tangent, adjoint)
      if (verbose()) write(*,'(a,2es14.6)') '        f whole and split ', f_whole(0), f_split(0)
      call verdict(index, label, f_whole, f_split, tangent, adjoint, failures)
    end subroutine one_case
    subroutine both_ways(whole, split, instants, half, degrees, duration, design, &
         & f_whole, f_split, achieved)
      type(family_holder), intent(in)  :: whole(:), split(:)
      integer            , intent(in)  :: instants, half, degrees
      real(dp)           , intent(in)  :: duration, design
      real(dp), allocatable, intent(out) :: f_whole(:), f_split(:)
      real(dp)           , intent(out) :: achieved
      real(dp) :: one, two
      call expanded(whole, [instants], degrees, duration, design, f_whole, one)
      call expanded(split, [half, instants - half], degrees, duration, design, &
           & f_split, two)
      achieved = max(one, two)
    end subroutine both_ways
    subroutine mixed_case(from, index, failures, skipped)
      integer, intent(in)    :: from, index
      integer, intent(inout) :: failures, skipped
      real(dp), parameter :: delta = 1.0e-4_dp
      type(family_holder), allocatable :: schemes(:)
      integer , allocatable :: added(:)
      real(dp), allocatable :: f(:), plus(:), minus(:)
      real(dp) :: duration, design, achieved, gap, differenced
      integer :: degrees, blocks, m
      character(len=28) :: label
      blocks = 0
      call draw_chain(from, degrees, blocks, duration, design, schemes, added, label)
      call expanded(schemes, added, degrees, duration, design, f, achieved)
      if (achieved > 1.0e-6_dp) then
         skipped = skipped + 1
         write(*,'(i6,2x,a28,a)') index, label, '  march did not converge, passed over'
         return
      end if
      call expanded(schemes, added, degrees, duration, design + delta, plus, achieved)
      call expanded(schemes, added, degrees, duration, design - delta, minus, achieved)
      gap = 0.0_dp
      do m = 1, max_order
         differenced = (plus(m - 1) - minus(m - 1)) / (2.0_dp * delta)
         gap = max(gap, abs(f(m) - differenced) / max(1.0_dp, abs(f(m))))
      end do
      write(*,'(i6,2x,a28,es14.2)') index, label, gap
      if (gap > 1.0e-3_dp) failures = failures + 1
    end subroutine mixed_case
    subroutine draw_chain(from, degrees, blocks, duration, design, schemes, added, label)
      integer            , intent(in)  :: from
      integer            , intent(out) :: degrees, blocks
      real(dp)           , intent(out) :: duration, design
      type(family_holder), allocatable, intent(inout) :: schemes(:)
      integer            , allocatable, intent(inout) :: added(:)
      character(len=*)   , intent(out) :: label
      integer(int64) :: state
      integer :: b, kind, order, widest
      state    = int(from, int64)
      degrees  = drawn(state, 2) + 2
      blocks   = drawn(state, 2) + 1
      duration = drawn_real(state, 0.5_dp, 3.0_dp)
      design   = drawn_real(state, 0.0_dp, 1.5_dp)
      if (allocated(schemes)) deallocate(schemes)
      if (allocated(added))   deallocate(added)
      allocate(schemes(blocks), added(blocks))
      label = ''
      do b = 1, blocks
         kind  = drawn(state, 3)
         order = drawn(state, 3)
         call fill(schemes(b), kind, order)
         added(b) = schemes(b) % scheme % history_depth(degrees - 1) + 2 + drawn(state, 3)
         if (b > 1) label = trim(label) // '-'
         label = trim(label) // trim(named(kind, order))
      end do
      widest = 0
      do b = 2, blocks
         widest = max(widest, schemes(b) % scheme % history_depth(degrees - 1))
      end do
      added(1) = max(added(1), widest + 1)
    end subroutine draw_chain
    subroutine draw(from, degrees, order, kind, instants, duration, design)
      integer , intent(in)  :: from
      integer , intent(out) :: degrees, order, kind, instants
      real(dp), intent(out) :: duration, design
      integer(int64) :: state
      state    = int(from, int64)
      degrees  = drawn(state, 2) + 2
      order    = drawn(state, 3)
      kind     = drawn(state, 3)
      instants = 12 + 2 * drawn(state, 5)
      duration = drawn_real(state, 0.5_dp, 4.0_dp)
      design   = drawn_real(state, 0.0_dp, 1.5_dp)
    end subroutine draw
    subroutine halved(scheme, degrees, instants, half)
      class(family), intent(in)    :: scheme
      integer      , intent(in)    :: degrees
      integer      , intent(inout) :: instants
      integer      , intent(out)   :: half
      integer :: reach
      reach = scheme % history_depth(degrees - 1)
      half  = max(instants / 2, reach + 1)
      if (instants - half <= reach) then
         instants = 2 * (reach + 1)
         half     = instants / 2
      end if
    end subroutine halved
    subroutine fill(held, kind, order)
      type(family_holder), intent(out) :: held
      integer            , intent(in)  :: kind, order
      select case (kind)
      case (1)
         call set_family(held, bdf_family(order))
      case (2)
         call set_family(held, adams_family(order))
      case default
         call set_family(held, crouzeix_two_stage())
      end select
    end subroutine fill
    function named(kind, order) result(text)
      integer, intent(in) :: kind, order
      character(len=16) :: text
      character(len=1) :: digit
      write(digit,'(i1)') order
      select case (kind)
      case (1)
         text = 'bdf' // digit
      case (2)
         text = 'adams' // digit
      case default
         text = 'crouzeix2'
      end select
    end function named
    subroutine expanded(schemes, added, degrees, duration, design, f, achieved)
      type(family_holder), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:), degrees
      real(dp)           , intent(in) :: duration, design
      real(dp), allocatable, intent(out) :: f(:)
      real(dp), allocatable :: table(:,:)
      real(dp)           , intent(out) :: achieved
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      real(dp), allocatable :: held(:), dt(:), t(:)
      call cosine_partition(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, max_order, table)
      allocate(f(lbound(table, 1):ubound(table, 1)))
      f = table(:, 1)
    end subroutine expanded
    subroutine verdict(index, label, f_whole, f_split, tangent, adjoint, failures)
      integer         , intent(in)    :: index
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: f_whole(0:), f_split(0:), tangent, adjoint
      integer         , intent(inout) :: failures
      real(dp) :: split_gap, direction_gap, route_gap, scale
      scale         = max(1.0_dp, maxval(abs(f_whole)))
      split_gap     = maxval(abs(f_whole - f_split)) / scale
      direction_gap = abs(tangent - adjoint) / max(1.0_dp, abs(tangent))
      route_gap     = abs(tangent - f_whole(1)) / max(1.0_dp, abs(tangent))
      write(*,'(i6,2x,a16,3es15.2)') index, label, split_gap, direction_gap, route_gap
      if (split_gap > 1.0e-6_dp) failures = failures + 1
      if (direction_gap > 1.0e-8_dp) failures = failures + 1
      if (route_gap > 1.0e-6_dp) failures = failures + 1
    end subroutine verdict
  end subroutine demo_randomized_checks
  subroutine demo_scheme_weights()
    implicit none
    integer :: k
    call bdf_rows(2, 'uniform',     [0.0_dp, (0.5_dp, k = 2, 5)])
    call bdf_rows(2, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp])
    call adams_row(3, 'uniform',     [0.0_dp, 0.5_dp, 0.5_dp])
    call adams_row(3, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp])
    call weight_partials(2, [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp])
  contains
    pure function instants(dt) result(t)
      real(dp), intent(in) :: dt(:)
      real(dp) :: t(size(dt))
      integer :: i
      t(1) = 0.0_dp
      do i = 2, size(dt)
         t(i) = t(i - 1) + dt(i)
      end do
    end function instants
    pure real(dp) function power_derivative(m, d, t) result(q)
      integer , intent(in) :: m, d
      real(dp), intent(in) :: t
      integer :: i
      if (d > m) then
         q = 0.0_dp
         return
      end if
      q = 1.0_dp
      do i = 0, d - 1
         q = q * real(m - i, dp)
      end do
      q = q * t**(m - d)
    end function power_derivative
    subroutine row_fields(scheme, nv, tails, head, source_degree, determines, dt, &
         & tau, alpha, w)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: nv, tails(:), head, source_degree(:), determines(:)
      real(dp)     , intent(in) :: dt(:)
      real(dp), allocatable, intent(out) :: tau(:), alpha(:), w(:)
      type(stored_directed_graph) :: coupling
      type(stored_field), allocatable :: inputs(:)
      class(field), allocatable :: out
      type(scheme_weight) :: weights
      integer :: e
      call coupling_inputs(nv, tails, [(head, e = 1, size(tails))], dt, source_degree, determines, &
           & coupling, inputs)
      tau = [(step_power(dt(head), source_degree(e) - determines(e)), e = 1, size(tails))]
      call scheme % apply(coupling, inputs, out)
      call out % real_vector(alpha)
      weights = scheme_weight(scheme)
      call weights % apply(coupling, inputs, out)
      call out % real_vector(w)
    end subroutine row_fields
    pure real(dp) function step_power(h, n) result(p)
      real(dp), intent(in) :: h
      integer , intent(in) :: n
      integer :: i
      p = 1.0_dp
      do i = 1, abs(n)
         p = p * h
      end do
      if (n < 0) p = 1.0_dp / p
    end function step_power
    pure real(dp) function row_residual(w, tails, head, source_degree, determines, t, m) &
         & result(r)
      real(dp), intent(in) :: w(:), t(:)
      integer , intent(in) :: tails(:), head, source_degree(:), determines(:), m
      integer :: e
      r = -power_derivative(m, determines(1), t(head))
      do e = 1, size(w)
         r = r + w(e) * power_derivative(m, source_degree(e), t(tails(e)))
      end do
    end function row_residual
    subroutine one_row(title, scheme, nv, tails, head, source_degree, determines, dt, top)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      integer         , intent(in) :: nv, tails(:), head, source_degree(:), determines(:), top
      real(dp)        , intent(in) :: dt(:)
      real(dp), allocatable :: tau(:), alpha(:), w(:)
      real(dp) :: t(nv), residual(0:top)
      integer :: m
      call row_fields(scheme, nv, tails, head, source_degree, determines, dt, tau, alpha, w)
      t = instants(dt)
      write(*,'(a)') ' '
      write(*,'(a)')        ' ' // title
      write(*,'(a,9f11.5)') '   tau                        ', tau
      write(*,'(a,9f11.5)') '   alpha                      ', alpha
      write(*,'(a,9f11.5)') '   weight                     ', w
      do m = 0, top
         residual(m) = row_residual(w, tails, head, source_degree, determines, t, m)
      end do
      write(*,'(a,9i11)')     '   on t**m, m =              ', [(m, m = 0, top)]
      write(*,'(a,9es11.2)')  '   residual                  ', residual
    end subroutine one_row
    subroutine bdf_rows(p, label, dt)
      integer         , intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      integer :: last
      last = 2 * p + 1
      call one_row('bdf ' // digit(p) // ' velocity row, ' // label // ' grid', &
           & bdf_family(p), last, [(last - k, k = 0, p)], last, &
           & [(0, k = 0, p)], [(1, k = 0, p)], dt, p + 2)
      call one_row('bdf ' // digit(p) // ' acceleration row, ' // label // ' grid', &
           & bdf_family(p), last, [(last - k, k = 0, 2 * p)], last, &
           & [(0, k = 0, 2 * p)], [(2, k = 0, 2 * p)], dt, p + 2)
    end subroutine bdf_rows
    subroutine adams_row(p, label, dt)
      integer         , intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      call one_row('adams-moulton ' // digit(p) // ' velocity row, ' // label // ' grid', &
           & adams_family(p), p, [p - 1, (p - k, k = 0, p - 1)], p, &
           & [1, (2, k = 0, p - 1)], [(1, k = 0, p)], dt, p + 2)
    end subroutine adams_row
    subroutine weight_partials(p, dt)
      integer , intent(in) :: p
      real(dp), intent(in) :: dt(:)
      real(dp), parameter :: delta = 1.0e-6_dp
      type(stored_directed_graph) :: coupling
      type(stored_field), allocatable :: inputs(:)
      type(stored_field) :: direction
      type(scheme_weight) :: weights
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)
      integer , allocatable :: tails(:)
      integer :: last, e, j
      last  = 2 * p + 1
      tails = [(last - j, j = 0, p)]
      call coupling_inputs(last, tails, [(last, e = 1, size(tails))], dt, &
           & [(0, e = 1, size(tails))], [(1, e = 1, size(tails))], coupling, inputs)
      direction = stored_field('v', coupling % vertex_set(), last)
      allocate(v(last), source=0.0_dp)
      v(last) = 1.0_dp
      call direction % set_real_vector(v)
      weights = scheme_weight(bdf_family(p))
      call weights % partial_action(coupling, inputs, [variation(weights % argument(1), direction)], out)
      call out % real_vector(exact)
      call inputs(1) % set_real_vector(dt + delta * v)
      call weights % apply(coupling, inputs, out)
      call out % real_vector(plus)
      call inputs(1) % set_real_vector(dt - delta * v)
      call weights % apply(coupling, inputs, out)
      call out % real_vector(minus)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2 velocity weights, partial in the last step'
      write(*,'(a,3f13.5)') '   partial_action             ', exact
      write(*,'(a,3f13.5)') '   central difference         ', (plus - minus) / (2.0_dp * delta)
    end subroutine weight_partials
    pure function digit(n) result(c)
      integer, intent(in) :: n
      character(len=1) :: c
      write(c,'(i1)') n
    end function digit
  end subroutine demo_scheme_weights
  subroutine demo_sensitivity()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    integer , parameter :: num_instants = 21
    real(dp), parameter :: duration = 2.0_dp
    real(dp), parameter :: design = 1.0_dp
    call sensitivity_of('bdf 2', bdf_family(2))
    call sensitivity_of('adams-moulton 3', adams_family(3))
  contains
    real(dp) function marched(scheme, design_value, q) result(f)
      class(family), intent(in) :: scheme
      real(dp)     , intent(in) :: design_value
      real(dp), allocatable, intent(out) :: q(:)
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_holder) :: holder(1)
      integer, allocatable :: at(:)
      type(stored_directed_graph) :: unknowns, instants
      type(stored_field) :: state, knobs
      real(dp), allocatable :: dt(:), t(:), held(:)
      real(dp) :: achieved
      call cosine_partition(scheme, degrees, duration, num_instants, held, dt, t)
      call set_family(holder(1), scheme)
      call tower % build(van_der_pol(state_degree), holder, [num_instants], uniform_grid(duration), &
           & 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), held, rows, at)
      call solved(rows, design_value, q, achieved)
      unknowns = unknowns_graph(num_instants, degrees)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', unknowns % vertex_set(), size(q))
      knobs    = stored_field('nu', unknowns % vertex_set(), num_instants)
      call state % set_real_vector(q)
      call knobs % set_real_vector(spread(design_value, 1, num_instants))
      f = functional_of(van_der_pol_energy(state_degree), instants, [state, knobs], dt)
    end function marched
    subroutine sensitivity_of(title, scheme)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      real(dp), parameter :: delta = 1.0e-6_dp
      real(dp), allocatable :: q(:), plus(:), minus(:)
      real(dp) :: f, tangent, adjoint, differenced
      f = marched(scheme, design, q)
      call three_objects(scheme, q, tangent, adjoint)
      differenced = (marched(scheme, design + delta, plus) - &
           &         marched(scheme, design - delta, minus)) / (2.0_dp * delta)
      write(*,'(a)')        ' '
      write(*,'(a)')        ' ' // title // ', van der pol at a design of one'
      write(*,'(a,f16.10)') '   the functional             ', f
      write(*,'(a,f16.10)') '   sensitivity, tangent       ', tangent
      write(*,'(a,f16.10)') '   sensitivity, adjoint       ', adjoint
      write(*,'(a,f16.10)') '   sensitivity, differenced   ', differenced
      write(*,'(a,es16.2)') '   tangent against adjoint    ', abs(tangent - adjoint)
      write(*,'(a,es16.2)') '   tangent against difference ', abs(tangent - differenced)
    end subroutine sensitivity_of
    subroutine three_objects(scheme, q, tangent, adjoint)
      class(family), intent(in) :: scheme
      real(dp)     , intent(in) :: q(:)
      real(dp)     , intent(out) :: tangent, adjoint
      real(dp), allocatable :: g(:), rate(:)
      integer :: mark
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_holder) :: holder(1)
      integer, allocatable :: at(:)
      type(stored_directed_graph) :: unknowns, instants
      type(stored_field) :: state, knobs
      real(dp), allocatable :: dt(:), t(:), held(:)
      call cosine_partition(scheme, degrees, duration, num_instants, held, dt, t)
      call set_family(holder(1), scheme)
      call tower % build(van_der_pol(state_degree), holder, [num_instants], uniform_grid(duration), &
           & 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), held, rows, at)
      unknowns = unknowns_graph(num_instants, degrees)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      state    = stored_field('state', unknowns % vertex_set(), size(q))
      knobs    = stored_field('nu', unknowns % vertex_set(), num_instants)
      call state % set_real_vector(q)
      call knobs % set_real_vector(spread(design, 1, num_instants))
      call functional_gradient(van_der_pol_energy(state_degree), instants, &
           & [state, knobs], dt, num_instants, degrees, unknowns % vertex_set(), g)
      call sweep_design_partial(rows, unknowns, [state, knobs], num_instants, &
           & unknowns % vertex_set(), rate)
      mark    = fresh_stamp()
      tangent = by_tangent(rows, unknowns, [state, knobs], g, rate, 0.0_dp, mark)
      adjoint = by_adjoint(rows, unknowns, [state, knobs], g, rate, 0.0_dp, mark)
    end subroutine three_objects
  end subroutine demo_sensitivity
  subroutine demo_solve_cost()
    implicit none
    integer, parameter :: sizes(5) = [41, 61, 81, 101, 121]
    integer :: k
    write(*,'(a)') ' '
    write(*,'(a)') '  unknowns   march(s)   form(s)   solve(s)     achieved     newton tolerance'
    do k = 1, 5
       call cost_at(sizes(k))
    end do
  contains
    subroutine cost_at(instants)
      integer, intent(in) :: instants
      integer, parameter :: degrees = 3
      type(family_holder), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      type(expression)       :: energy(1)
      type(bdf_family) :: scheme
      integer , allocatable :: added(:)
      real(dp), allocatable :: held(:), dt(:), t(:), table(:,:)
      real(dp) :: achieved, duration, design, marched, formed, solved_in, tangent
      integer  :: n
      duration = 3.0_dp
      design   = 1.0_dp
      scheme   = bdf_family(2)
      allocate(schemes(1))
      call set_family(schemes(1), scheme)
      added = [instants]
      call cosine_partition(scheme, degrees, duration, instants, held, dt, t)
      marched = clock()
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      marched = clock() - marched
      energy(1) = van_der_pol_energy(degrees - 1)
      formed = clock()
      call chain_stamps(chain, tower, energy, degrees, marks)
      formed = clock() - formed
      n = chain(1) % rows % num_unknowns()
      solved_in = clock()
      call chain_derivative(chain, tower, marks, energy, degrees, 1, forward_route, table)
      tangent   = first_of(table)
      solved_in = clock() - solved_in
      write(*,'(i10,3f11.3,2es15.3)') n, marched, formed, solved_in, &
           & achieved, 1.0e-12_dp
    end subroutine cost_at
  end subroutine demo_solve_cost
  subroutine demo_tolerance_form()
    implicit none
    write(*,'(a)') ' '
    write(*,'(a)') '  the velocity row, and the d-th row it composes to'
    write(*,'(a)') '  scheme    d    sum|a|    sum|c(d)|   (sum|a|)^d   attained'
    call composed(1, 2)
    call composed(2, 2)
    call composed(3, 2)
    call composed(4, 2)
    call composed(2, 3)
    call composed(3, 3)
    write(*,'(a)') ' '
    write(*,'(a)') '  ||A||_inf, from the family against the assembled jacobian'
    write(*,'(a)') '  scheme    d   instants       dt      from family     assembled   agreement'
    call against(1, 3, 21)
    call against(2, 3, 21)
    call against(3, 3, 21)
    call against(2, 3, 61)
    call against(2, 4, 61)
    call against(3, 3, 41)
    write(*,'(a)') ' '
    write(*,'(a)') '  conditioning, and what the same form removes from it'
    write(*,'(a)') '  the row determining degree d is divided by dt^-d, which is the'
    write(*,'(a)') '  weight the family says it carries'
    write(*,'(a)') ' '
    write(*,'(a)') '  scheme    d   instants       dt      kappa(A)   kappa(DA)     ratio'
    call conditioned(2, 3, 21)
    call conditioned(2, 3, 41)
    call conditioned(2, 3, 61)
    call conditioned(2, 3, 81)
    call conditioned(3, 3, 41)
    call conditioned(2, 4, 41)
  contains
    real(dp) function row_sum(scheme, order, determines) result(total)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: order, determines
      real(dp), allocatable :: c(:)
      integer :: reach, last, k
      reach = determines * order
      last  = reach + 1
      call weights_of(scheme, last, [(last - k, k = 0, reach)], [(last, k = 0, reach)], &
           & [(1.0_dp, k = 1, last)], [(0, k = 0, reach)], [(determines, k = 0, reach)], c)
      total = sum(abs(c))
    end function row_sum
    subroutine composed(order, determines)
      integer, intent(in) :: order, determines
      real(dp) :: velocity, derived, powered
      character(len=8) :: named
      velocity = row_sum(bdf_family(order), order, 1)
      derived  = row_sum(bdf_family(order), order, determines)
      powered  = velocity ** determines
      write(named,'(a,i0)') 'bdf ', order
      write(*,'(a,a,i5,3f12.4,a)') '  ', named, determines, velocity, derived, powered, &
           & merge('   yes', '    no', abs(derived - powered) <= 1.0e-10_dp * powered)
    end subroutine composed
    subroutine against(order, degrees, instants)
      integer, intent(in) :: order, degrees, instants
      type(family_holder), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      type(bdf_family) :: scheme
      integer , allocatable :: added(:)
      real(dp), allocatable :: held(:), dt(:), t(:), a(:,:)
      real(dp) :: achieved, duration, design, predicted, assembled
      character(len=8) :: named
      integer :: k, d, top
      duration = 3.0_dp
      design   = 1.0_dp
      scheme   = bdf_family(order)
      top      = degrees - 1
      allocate(schemes(1))
      call set_family(schemes(1), scheme)
      added = [instants]
      call partition(duration, instants, dt, t)
      held = [((0.0_dp, d = 0, degrees - 1), k = 1, &
           &   scheme % history_depth(degrees - 1))]
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      call chain_stamps(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, marks)
      predicted = 1.0_dp + row_sum(scheme, order, top) / dt(size(dt)) ** top
      call dense_jacobian(chain, design, a)
      assembled = largest_row(a)
      write(named,'(a,i0)') 'bdf ', order
      write(*,'(a,a,i5,i10,f10.5,2es15.5,f11.3,a)') '  ', named, top, instants, &
           & dt(size(dt)), predicted, assembled, &
           & 100.0_dp * (1.0_dp - abs(predicted - assembled) / assembled), ' %'
    end subroutine against
    subroutine conditioned(order, degrees, instants)
      integer, intent(in) :: order, degrees, instants
      type(family_holder), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: marks(:)
      type(bdf_family) :: scheme
      integer , allocatable :: added(:)
      real(dp), allocatable :: held(:), dt(:), t(:), b(:,:), a(:,:)
      real(dp) :: achieved, duration, design, bare, scaled, step
      character(len=8) :: named
      integer :: k, d, i, n
      duration = 3.0_dp
      design   = 1.0_dp
      scheme   = bdf_family(order)
      allocate(schemes(1))
      call set_family(schemes(1), scheme)
      added = [instants]
      call partition(duration, instants, dt, t)
      held = [((0.0_dp, d = 0, degrees - 1), k = 1, &
           &   scheme % history_depth(degrees - 1))]
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
      call chain_stamps(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, marks)
      step = dt(size(dt))
      call dense_jacobian(chain, design, a)
      bare = kappa(a)
      n = size(a, 1)
      allocate(b(n, n))
      do i = 1, n
         d = mod(i - 1, degrees)
         b(i, :) = a(i, :) * step ** d
      end do
      scaled = kappa(b)
      write(named,'(a,i0)') 'bdf ', order
      write(*,'(a,a,i5,i10,f10.5,2es13.4,f10.1)') '  ', named, degrees - 1, &
           & instants, step, bare, scaled, bare / scaled
    end subroutine conditioned
    real(dp) function kappa(a) result(k)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: w(:,:), inverse(:,:), row(:)
      real(dp) :: pivot, factor
      integer :: n, i, j, p
      n = size(a, 1)
      allocate(w(n, n), inverse(n, n), row(n))
      w       = a
      inverse = 0.0_dp
      do i = 1, n
         inverse(i, i) = 1.0_dp
      end do
      do j = 1, n
         p = j - 1 + maxloc(abs(w(j:n, j)), dim=1)
         if (p /= j) then
            row          = w(j, :)
            w(j, :)      = w(p, :)
            w(p, :)      = row
            row          = inverse(j, :)
            inverse(j, :) = inverse(p, :)
            inverse(p, :) = row
         end if
         pivot = w(j, j)
         if (abs(pivot) <= tiny(1.0_dp)) then
            k = huge(1.0_dp)
            return
         end if
         w(j, :)       = w(j, :) / pivot
         inverse(j, :) = inverse(j, :) / pivot
         do i = 1, n
            if (i == j) cycle
            factor = w(i, j)
            w(i, :)       = w(i, :) - factor * w(j, :)
            inverse(i, :) = inverse(i, :) - factor * inverse(j, :)
         end do
      end do
      k = largest_row(a) * largest_row(inverse)
    end function kappa
    real(dp) function largest_row(a) result(most)
      real(dp), intent(in) :: a(:,:)
      most = maxval(sum(abs(a), dim=2))
    end function largest_row
  end subroutine demo_tolerance_form
end module gti_demos
program graph_time_integrator
  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : uniform_grid, random_grid, designed_grid, fixed_grid
  use operation_expression  , only : expression
  use gti_physics           , only : van_der_pol, van_der_pol_energy
  use operation_grid        , only : grid
  use gti_march             , only : set_stopping, imbalance, set_sweep, weight_of, precision_needed
  use gti_adaptive          , only : adaptive_partition
  use operation_family_dirk , only : crouzeix_three_stage
  use operation_stencil     , only : stencil
  use gti_space             , only : room, spatial_mesh, geometry_of, coarse_cells
  use gti_field             , only : spatial_discretization_stencil_of, initial_field, against_the_laplacian, &
       & against_the_mode, export_instant
  use util_precision        , only : precision_named
  use iso_fortran_env       , only : real128
  use gti_expansion         , only : family_holder, expansion
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & expansion_substitutions, chain_stamps, num_designs_of, &
       & instant_components, chain_derivative, asymmetry, sink_costates
  use gti_sweeps            , only : set_linear_solver, set_assembly, set_storage, set_multigrid, &
       & set_coarse_nodes, set_linear_budget, set_newton_order
  use gti_sweeps            , only : route_of, forward_route, reverse_route
  use operation_minimization, only : relative, absolute, by_count, by_rate
  use gti_driver            , only : settings, chosen_grid, steps_of, family_named, clock, &
       & functional_named
  use gti_configuration     , only : configuration, read_configuration, override, show, &
       & lists, refuse_unknown, worded
  use util_tally            , only : tally_open, tally_close, tally_order, &
       & tally_enter, tally_leave, tally_amount, tally_event_of, &
       & tally_num_levels, tally_level_name, tally_event_name, wall_time
  use gti_configuration     , only : at_expansion, at_horizon, hierarchy_levels
  use gti_demos           , only : demo_requested, run_demo
  implicit none
  type(configuration) :: cfg
  type(room)   , allocatable :: space
  type(stencil), allocatable :: spatial_discretization_stencil
  real(dp)     , allocatable :: volume(:), q0(:)
  real(dp) :: extent_a = 0.0_dp, extent_b = 0.0_dp
  integer  :: nodes = 1
  logical  :: over_field = .false.
  type(expression)      , allocatable :: functionals(:)
  logical :: grid_designed = .false.
  logical :: grid_adaptive = .false.
  real(dp), allocatable :: adaptive_weights(:)
  if (demo_requested()) then
     call run_demo()
     stop
  end if
  call settings('homogeneous', cfg)
  call show(cfg)
  call set_linear_solver(cfg % linear_solver)
  call set_newton_order(cfg % higher_order_jacobian_product)
  call set_assembly(cfg % assembly)
  call set_storage(cfg % storage)
  call set_multigrid(cfg % multigrid)
  call set_sweep(cfg % sweep)
  call field_context(cfg)
  call chosen_functionals(cfg)
  call table(cfg)
contains
  integer function widest_reach(cfg) result(widest)
    type(configuration), intent(in) :: cfg
    character(len=8) :: every(3)
    class(family), allocatable :: scheme
    logical :: staged, ok
    integer :: i, order, reach
    every  = ['bdf     ', 'adams   ', 'dirk    ']
    widest = 0
    do i = 1, 3
       if (index(cfg % families, trim(every(i))) == 0) cycle
       do order = 1, cfg % max_discretization_order
          call chosen(trim(every(i)), order, scheme, staged, ok)
          if (.not. ok) cycle
          reach = scheme % history_depth(cfg % state_degree)
          if (reach < cfg % instants) widest = max(widest, reach)
       end do
    end do
  end function widest_reach
  subroutine chosen(name, order, scheme, staged, ok)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: staged, ok
    call family_named(name, order, scheme, ok)
    staged = .false.
    if (ok) staged = scheme % num_stages() > 1
  end subroutine chosen
  function labelled(names, orders) result(text)
    character(len=*), intent(in) :: names(:)
    integer         , intent(in) :: orders(:)
    character(len=:), allocatable :: text
    character(len=2) :: digit
    integer :: b
    text = ''
    do b = 1, size(names)
       write(digit,'(i0)') orders(b)
       if (b > 1) text = text // '-'
       text = text // trim(names(b)) // trim(digit)
    end do
  end function labelled
  subroutine shown_initial(cfg)
    type(configuration), intent(in) :: cfg
    real(dp), allocatable :: q(:)
    character(len=:), allocatable :: line
    character(len=18) :: cell
    integer :: i
    q = q0(1:cfg % state_degree + 1)
    line = '   initial state, consistent'
    do i = 1, size(q)
       write(cell,'(es18.10)') q(i)
       line = line // cell
    end do
    if (over_field) write(cell,'(a,i0)') '   at node 1 of ', nodes
    if (over_field) line = line // trim(cell)
    write(*,'(a)') line
  end subroutine shown_initial
  subroutine heading(cfg)
    type(configuration), intent(in) :: cfg
    character(len=:), allocatable :: line, name
    integer :: m
    line = '  scheme' // repeat(' ', 14) // 'solved' // repeat(' ', 10)
    do m = 0, cfg % max_derivative_degree
       name = order_named(m)
       line = line // repeat(' ', 20 - len(name)) // name // ' '
    end do
    write(*,'(a)') ' '
    write(*,'(a)') line
  end subroutine heading
  subroutine show_row(label, solved, f, left, columns)
    character(len=*), intent(in) :: label
    integer         , intent(in) :: solved
    real(dp)        , intent(in) :: f(0:)
    type(imbalance) , intent(in) :: left
    integer         , intent(in) :: columns
    character(len=21) :: cell
    character(len=6)  :: counted
    character(len=:), allocatable :: line
    integer :: m
    line = '  ' // label // repeat(' ', max(2, 20 - len(label)))
    write(counted,'(i6)') solved
    line = line // counted // repeat(' ', 10)
    do m = 0, columns
       if (m <= ubound(f, 1)) then
          write(cell,'(es20.11)') f(m)
       else
          write(cell,'(a20)') '-'
       end if
       line = line // cell
    end do
    if (.not. left % converged) then
       if (left % diverging) then
          line = line // '   diverging'
       else
          line = line // '   unconverged'
       end if
    end if
    write(*,'(a)') line
    if (.not. left % converged) call shown_aspect(left)
  end subroutine show_row
  subroutine shown_precision(nd, chain, cfg)
    integer            , intent(in) :: nd
    type(chain_block)  , intent(in) :: chain(:)
    type(configuration), intent(in) :: cfg
    real(dp) :: weight, state_size
    real(real128) :: needed
    character(len=:), allocatable :: least
    logical :: shown
    integer :: b
    shown = cfg % accounting
    do b = 1, size(chain)
       call precision_needed(weight_of(chain(b) % scheme, nd, minval(chain(b) % dt(2:))), &
            & maxval(abs(chain(b) % state)), chain(b) % began, needed, least)
       if (least /= precision_named() .and. least /= 'single') shown = .true.
    end do
    if (.not. shown) return
    do b = 1, size(chain)
       weight     = weight_of(chain(b) % scheme, nd, minval(chain(b) % dt(2:)))
       state_size = maxval(abs(chain(b) % state))
       call precision_needed(weight, state_size, chain(b) % began, needed, least)
       write(*,'(a,i0,a,es9.2,a,es9.2,a,es9.2,a,a,a,a)') '      precision, block ', b, &
            & '  ||A|| ', weight, '  ||q|| ', state_size, '  spacing needed ', real(needed, dp), &
            & '  least kind ', least, '  this build ', precision_named()
    end do
  end subroutine shown_precision
  subroutine shown_aspect(left)
    type(imbalance), intent(in) :: left
    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: d
    write(*,'(a,es10.3,a,es10.3,a)') '      imbalance ', left % norm, &
         & ' against ', left % began, ' where the march began'
    line = '      by degree '
    do d = 0, ubound(left % by_degree, 1)
       write(cell,'(es14.3)') left % by_degree(d)
       line = line // cell
    end do
    write(*,'(a)') line
    write(*,'(a,i0,a,i0)') '      largest entry at slot ', left % worst_slot, &
         & ' degree ', left % worst_degree
    write(*,'(a,i0,a,i0,a,es10.3)') '      steepest in the state at slot ', &
         & left % steepest_slot, ' degree ', left % steepest_degree, &
         & ', d||r||/dq = ', left % steepest
  end subroutine shown_aspect
  subroutine one_row(cfg, names, orders, printed)
    type(configuration), intent(in)    :: cfg
    character(len=*)   , intent(in)    :: names(:)
    integer            , intent(in)    :: orders(:)
    integer            , intent(inout) :: printed
    type(family_holder), allocatable :: schemes(:)
    type(chain_block)  , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    integer , allocatable :: added(:)
    real(dp), allocatable :: dt(:), t(:), f(:,:), weights(:)
    type(imbalance) :: left
    real(dp) :: achieved
    integer :: nd, width, given, reported, i, m
    logical :: ok
    character(len=20) :: cell
    character(len=:), allocatable :: line
    nd    = cfg % state_degree + 1
    width = nd * nodes
    allocate(schemes(size(names)), added(size(names)))
    call assembled(cfg, names, orders, schemes, added, ok)
    if (.not. ok) return
    call grid_partition(cfg, dt, t)
    given = schemes(1) % scheme % history_depth(nd - 1)
    call tally_enter(at_expansion)
    call tally_order(0)
    if (grid_designed) then
       weights = dt(2:cfg % instants)
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & designed_grid(cfg % time_duration), cfg % design, q0, chain, tower, dt, t, &
            & achieved, grid_design=weights, left=left, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else if (grid_adaptive) then
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & fixed_grid(adaptive_weights), cfg % design, q0, chain, tower, dt, t, achieved, &
            & left=left, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & chosen_grid(cfg), cfg % design, q0, chain, tower, dt, t, achieved, left=left, &
            & nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, startup=cfg % startup_refinement)
    end if
    if (.not. left % converged) then
       reported = 0
    else
       reported = cfg % max_derivative_degree
    end if
    call chain_expansion(chain, tower, functionals, nd, reported, f, node_measure=volume)
    call tally_leave()
    call show_row(labelled(names, orders), cfg % instants - given, f(:, 1), left, &
         & cfg % max_derivative_degree)
    do i = 2, size(functionals)
       line = '      ' // functionals(i) % name()
       line = line // repeat(' ', max(1, 36 - len(line)))
       do m = 0, ubound(f, 1)
          write(cell,'(es20.11)') f(m, i)
          line = line // cell
       end do
       write(*,'(a)') line
    end do
    call shown_precision(nd, chain, cfg)
    if (reported >= 1 .and. (grid_designed .or. lists(cfg % check, 'routes') &
         & .or. lists(cfg % check, 'sinks'))) then
       call first_derivatives(cfg, chain, tower, nd, dt, f)
    end if
    if (over_field) then
       if (lists(cfg % check, 'ode')) call against_the_ode(cfg, schemes, added, f(:, 1))
       if (lists(cfg % check, 'mode')) then
          call against_the_mode(space, extent_a, extent_b, cfg % diffusion, cfg % spatial_order, &
               & cfg % design, t(cfg % instants), instant_components(chain, cfg % instants), nd)
       end if
       if (trim(cfg % export) == 'paraview') call exported(cfg, chain, labelled(names, orders), nd)
    end if
    printed = printed + 1
  end subroutine one_row
  subroutine first_derivatives(cfg, chain, tower, nd, dt, f)
    type(configuration), intent(in) :: cfg
    type(chain_block)  , intent(in) :: chain(:)
    type(expansion)    , intent(in) :: tower
    integer            , intent(in) :: nd
    real(dp)           , intent(in) :: dt(:), f(0:, :)
    integer, allocatable :: marks(:)
    real(dp), allocatable :: p(:), df(:,:), other(:,:), table(:,:), entries(:,:,:)
    type(sink_costates) :: sinks
    real(dp) :: euler
    integer  :: num_designs, num_functionals, route, i, order
    num_functionals = size(functionals)
    call chain_stamps(chain, tower, functionals, nd, marks, node_measure=volume)
    num_designs = num_designs_of(tower)
    route = route_of(num_designs, num_functionals, 1)
    call chain_derivative(chain, tower, marks, functionals, nd, 1, route, df, node_measure=volume)
    write(*,'(a,a,a,i0,a,i0,a,es10.2)') '      first derivatives by the ', &
         & trim(merge('forward', 'reverse', route == forward_route)), ' route, designs ', &
         & num_designs, ' functionals ', num_functionals, &
         & ':  physics column against the expansion ', &
         & maxval(abs(df(:, 1) - f(1, :)) / max(1.0_dp, abs(f(1, :))))
    if (grid_designed) then
       p = dt(2:cfg % instants)
       do i = 1, num_functionals
          euler = dot_product(p, df(i, 2:)) / max(tiny(1.0_dp), norm2(p) * norm2(df(i, 2:)))
          write(*,'(a,i0,a,es12.4,a,es10.2)') '      grid design, functional ', i, &
               & ':  |df/dp| ', norm2(df(i, 2:)), '   p . df/dp / |p||df/dp| (theory 0) ', euler
       end do
    end if
    if (lists(cfg % check, 'routes')) then
       call chain_derivative(chain, tower, marks, functionals, nd, 1, &
            & merge(reverse_route, forward_route, route == forward_route), other, node_measure=volume)
       write(*,'(a,es10.2)') '      tangent against adjoint over the table, relative ', &
            & maxval(abs(df - other)) / max(1.0_dp, maxval(abs(df)))
    end if
    if (lists(cfg % check, 'sinks')) then
       call chain_derivative(chain, tower, marks, functionals, nd, 1, reverse_route, other, &
            & node_measure=volume, sinks=sinks)
       call shown_sinks(sinks, nd)
    end if
    if (grid_designed) then
       do order = 2, ubound(f, 1)
          route = route_of(num_designs, num_functionals, order)
          call chain_derivative(chain, tower, marks, functionals, nd, order, route, table, &
               & node_measure=volume, entries=entries)
          do i = 1, num_functionals
             if (route == reverse_route) then
                write(*,'(a,i0,a,a,i0,a,es12.4,a,es10.2,a,es10.2)') '      derivatives of order ', &
                     & order, ' by the reverse route, functional ', '', i, ':  |T| ', &
                     & maxval(abs(table(i, :))), '   departure among the entries of a multiset ', &
                     & asymmetry(entries, num_designs, order), &
                     & '   parameter entry against the expansion ', &
                     & abs(table(i, 1) - f(order, i)) / max(1.0_dp, abs(f(order, i)))
             else
                write(*,'(a,i0,a,i0,a,es12.4,a,es10.2)') '      derivatives of order ', order, &
                     & ' by the forward route, functional ', i, ':  |T| ', maxval(abs(table(i, :))), &
                     & '   parameter entry against the expansion ', &
                     & abs(table(i, 1) - f(order, i)) / max(1.0_dp, abs(f(order, i)))
             end if
          end do
       end do
    end if
  end subroutine first_derivatives
  subroutine shown_sinks(sinks, nd)
    type(sink_costates), intent(in) :: sinks
    integer            , intent(in) :: nd
    write(*,'(a,a,a,a,a,a)') '      sink costates: unknowns no row of their block reads, by degree 0..', &
         & trim(counted(nd, sinks % interior)), ' interior', trim(counted(nd, sinks % last)), &
         & ' at the last point', trim(counted(nd, sinks % carried)), ' carried'
    write(*,'(a,es10.2,a,es10.2)') '         J_ii lambda_i - g_i on the sinks, relative ', &
         & sinks % departure / max(1.0_dp, sinks % gradient), &
         & '   lambda where the functional reads nothing, relative ', &
         & sinks % unread / max(1.0_dp, sinks % costate)
  end subroutine shown_sinks
  function counted(nd, per_degree) result(line)
    integer, intent(in) :: nd, per_degree(0:)
    character(len=:), allocatable :: line
    character(len=16) :: word
    integer :: d
    line = ''
    do d = 0, nd - 1
       write(word, '(i0)') per_degree(d)
       line = line // merge(' ', '/', d == 0) // trim(word)
    end do
  end function counted
  subroutine chosen_functionals(cfg)
    type(configuration), intent(in) :: cfg
    character(len=32), allocatable :: names(:)
    logical :: ok
    integer :: i
    call refuse_unknown(cfg % designs, ['physics', 'grid   '], 'designs')
    call refuse_unknown(cfg % functionals, ['energy     ', 'dissipation'], 'functionals')
    if (.not. lists(cfg % designs, 'physics')) then
       error stop 'graph_time_integrator: the physics'' parameter is the first design'
    end if
    grid_designed = lists(cfg % designs, 'grid')
    names = worded(cfg % functionals)
    allocate(functionals(size(names)))
    do i = 1, size(names)
       call functional_named(trim(names(i)), cfg % state_degree, functionals(i), ok)
    end do
  end subroutine chosen_functionals
  subroutine against_the_ode(cfg, schemes, added, f_field)
    type(configuration), intent(in) :: cfg
    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:)
    real(dp)           , intent(in) :: f_field(0:)
    type(chain_block), allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    real(dp), allocatable :: f(:,:), dt(:), t(:)
    real(dp) :: achieved, area
    integer  :: nd, d
    character(len=:), allocatable :: line
    character(len=20) :: cell
    nd   = cfg % state_degree + 1
    area = sum(volume)
    call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, chosen_grid(cfg), &
         & cfg % design, q0(1:nd), chain, tower, dt, t, achieved, startup=cfg % startup_refinement)
    call chain_expansion(chain, tower, functionals, nd, ubound(f_field, 1), f)
    line = '      field / area over the node, less one:'
    do d = lbound(f, 1), ubound(f, 1)
       write(cell,'(es14.2)') f_field(d) / area / f(d, 1) - 1.0_dp
       line = line // cell
    end do
    write(*,'(a)') line
  end subroutine against_the_ode
  subroutine exported(cfg, chain, label, nd)
    type(configuration), intent(in) :: cfg
    type(chain_block)  , intent(in) :: chain(:)
    character(len=*)   , intent(in) :: label
    integer            , intent(in) :: nd
    character(len=len(label)) :: name
    character(len=256) :: path
    integer :: k, i
    name = label
    do i = 1, len(name)
       if (name(i:i) == ' ') name(i:i) = '_'
    end do
    do k = 1, cfg % instants
       write(path,'(a,a,a,a,i4.4,a)') trim(cfg % export_path), '_', trim(name), '_', k, '.vtu'
       call export_instant(space, trim(path), nd, instant_components(chain, k))
    end do
    write(*,'(a,i0,a,a,a)') '      written ', cfg % instants, ' files ', &
         & trim(cfg % export_path) // '_' // trim(name), '_*.vtu'
  end subroutine exported
  subroutine field_context(cfg)
    type(configuration), intent(in) :: cfg
    real(dp) :: x, y, began
    integer  :: n1, n2
    call refuse_unknown(cfg % initial_field, ['constant', 'mode    ', 'bump    '], 'initial_field')
    call refuse_unknown(cfg % export, ['none    ', 'paraview'], 'export')
    call refuse_unknown(cfg % check, ['none    ', 'ode     ', 'mode    ', 'operator', 'routes  ', &
         & 'sinks   '], &
         & 'check')
    call pair_of(cfg % spatial_counts, x, y, 'counts')
    n1 = nint(x)
    n2 = nint(y)
    if (real(n1, dp) /= x .or. real(n2, dp) /= y) then
       error stop 'graph_time_integrator: a count of cells is whole'
    end if
    over_field = n1 > 0 .or. n2 > 0
    if (over_field .and. (n1 <= 0 .or. n2 <= 0)) then
       error stop 'graph_time_integrator: a mesh has cells along both coordinates'
    end if
    if (over_field) then
       call refuse_unknown(cfg % spatial_grid, ['uniform', 'random '], 'spatial_grid')
       call pair_of(cfg % spatial_extent, extent_a, extent_b, 'extents')
       began = clock()
       allocate(space)
       space = spatial_mesh(geometry_of(cfg % spatial_geometry), extent_a, extent_b, n1, n2, &
            & trim(cfg % spatial_grid) == 'random', cfg % seed)
       write(*,'(a,i0,a,i0,a,f12.6,a,i0,a,f9.3,a)') '   spatial mesh: cells ', &
            & space % num_cells, '   faces ', space % num_faces, '   area ', sum(space % volume), &
            & '   form degree ', cfg % spatial_order, '   built in ', clock() - began, ' s'
       spatial_discretization_stencil = spatial_discretization_stencil_of(space, cfg % diffusion, cfg % spatial_order)
       call set_coarse_nodes(coarse_cells(space))
       nodes  = space % num_cells
       volume = space % volume
       if (lists(cfg % check, 'operator')) then
          call against_the_laplacian(space, extent_a, extent_b, cfg % diffusion, cfg % spatial_order)
       end if
    else
       nodes  = 1
       volume = [1.0_dp]
    end if
    q0 = initial_field(van_der_pol(cfg % state_degree), cfg % state_degree + 1, &
         & cfg % initial_field, cfg % initial_state, cfg % design, &
         & spatial_discretization_stencil=spatial_discretization_stencil, space=space, a=extent_a, b=extent_b)
  end subroutine field_context
  subroutine pair_of(text, x, y, subject)
    character(len=*), intent(in)  :: text, subject
    real(dp)        , intent(out) :: x, y
    character(len=32), allocatable :: w(:)
    w = worded(text)
    if (size(w) /= 2) error stop 'graph_time_integrator: two ' // subject // ', one per coordinate'
    read(w(1), *) x
    read(w(2), *) y
  end subroutine pair_of
  subroutine assembled(cfg, names, orders, schemes, added, ok)
    type(configuration), intent(in)  :: cfg
    character(len=*)   , intent(in)  :: names(:)
    integer            , intent(in)  :: orders(:)
    type(family_holder), intent(inout) :: schemes(:)
    integer            , intent(inout) :: added(:)
    logical            , intent(out)   :: ok
    class(family), allocatable :: scheme
    logical :: staged, exists
    integer :: b, blocks, share
    blocks = size(names)
    ok = .true.
    share = cfg % instants / blocks
    added = share
    added(1) = cfg % instants - share * (blocks - 1)
    do b = 1, blocks
       call chosen(names(b), orders(b), scheme, staged, exists)
       if (.not. exists) then
          ok = .false.
          cycle
       end if
       allocate(schemes(b) % scheme, source=scheme)
       deallocate(scheme)
       if (added(b) <= schemes(b) % scheme % history_depth(cfg % state_degree)) ok = .false.
    end do
  end subroutine assembled
  subroutine adaptive_context(cfg)
    type(configuration), intent(inout) :: cfg
    integer :: nd, rejects
    if (trim(cfg % grid) /= 'adaptive') return
    if (over_field) then
       error stop 'graph_time_integrator: an adaptive grid is over time alone'
    end if
    nd = cfg % state_degree + 1
    adaptive_weights = adaptive_partition(crouzeix_three_stage(), 4, &
         & van_der_pol(cfg % state_degree), nd, cfg % time_duration, &
         & q0(1:cfg % state_degree), cfg % design, cfg % tolerance, &
         & trim(cfg % tolerance_criterion) == 'relative', rejects)
    cfg % instants = size(adaptive_weights) + 1
    grid_adaptive  = .true.
    write(*,'(a,i0,a,es9.2,a,i0,a)') '   adaptive grid: ', size(adaptive_weights), &
         & ' steps to tolerance ', cfg % tolerance, ' (', rejects, ' rejected)'
  end subroutine adaptive_context
  subroutine grid_partition(cfg, dt, t)
    type(configuration), intent(in) :: cfg
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    integer :: k
    if (grid_adaptive) then
       allocate(dt(cfg % instants), t(cfg % instants))
       dt(1)  = 0.0_dp
       dt(2:) = adaptive_weights
       t(1)   = 0.0_dp
       do k = 2, cfg % instants
          t(k) = t(k - 1) + dt(k)
       end do
    else
       call steps_of(cfg, dt, t)
    end if
  end subroutine grid_partition
  subroutine table(cfg)
    type(configuration), intent(inout) :: cfg
    real(dp), allocatable :: dt(:), t(:)
    integer :: widest, printed
    call refuse_unknown(cfg % physics, ['vanderpol'], 'physics')
    call refuse_unknown(cfg % tolerance_criterion, ['relative', 'absolute'], &
         & 'tolerance_criterion')
    call refuse_unknown(cfg % iteration_criterion, ['by_rate ', 'by_count'], &
         & 'iteration_criterion')
    call set_stopping(cfg % tolerance, &
         & merge(relative, absolute, trim(cfg % tolerance_criterion) == 'relative'), &
         & merge(by_rate, by_count, trim(cfg % iteration_criterion) == 'by_rate'), &
         & cfg % max_iterations)
    call set_linear_budget(cfg % krylov_restart, cfg % smoothing_sweeps, &
         & cfg % max_linear_iterations)
    call adaptive_context(cfg)
    if (cfg % accounting) then
       call refuse_unknown(cfg % measurements, &
            & ['wall_time     ', 'primal_loops  ', 'tangent_loops ', &
            &  'adjoint_loops ', 'newton_solves ', 'linear_solves ', &
            &  'factorisations'], 'measurements')
    end if
    call refuse_unknown(cfg % families, ['bdf     ', 'adams   ', 'dirk    '], 'families')
    call refuse_unknown(cfg % combinations, &
         & ['homogeneous', 'pairs      ', 'triples    '], 'combinations')
    widest = widest_reach(cfg)
    if (.not. cfg % automatic_order_conservation) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' the widest row here looks back over instants: ', widest
       write(*,'(a)')    ' filling them by any other means leaves the rows solving different'
       write(*,'(a)')    ' problems from different starting states, and no table read across'
       write(*,'(a)')    ' such rows means anything.'
       error stop 'graph_time_integrator: order conservation is the only startup built'
    end if
    if (widest == 0) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' every family and order asked for looks further back than the'
       write(*,'(a,i0)') ' horizon holds, which is instants: ', cfg % instants
       error stop 'graph_time_integrator: no row fits in this horizon'
    end if
    call grid_partition(cfg, dt, t)
    call shown_initial(cfg)
    write(*,'(a,a)') '   precision of this build  ', precision_named()
    call heading(cfg)
    if (cfg % accounting) call tally_open(cfg % max_derivative_degree, hierarchy_levels)
    printed = 0
    if (asked(cfg, 'homogeneous')) call tuple_rows(cfg, 1, printed)
    if (asked(cfg, 'pairs'))       call tuple_rows(cfg, 2, printed)
    if (asked(cfg, 'triples'))     call tuple_rows(cfg, 3, printed)
    if (cfg % accounting) then
       call tally_close()
       call accounted(cfg)
    end if
    if (printed == 0) then
       write(*,'(a)') ' '
       write(*,'(a)') ' no row was built. A family has no scheme at every order - a stage'
       write(*,'(a)') ' family has none below order two - and a row whose blocks would add'
       write(*,'(a)') ' no more instants than they look back over is not built either.'
    end if
  end subroutine table
  pure logical function asked(cfg, what) result(yes)
    type(configuration), intent(in) :: cfg
    character(len=*)   , intent(in) :: what
    yes = lists(cfg % combinations, what)
  end function asked
  function listed(cfg) result(list)
    type(configuration), intent(in) :: cfg
    character(len=8), allocatable :: list(:)
    character(len=8) :: every(3)
    integer :: i, n
    every = ['bdf     ', 'adams   ', 'dirk    ']
    n = 0
    do i = 1, 3
       if (lists(cfg % families, trim(every(i)))) n = n + 1
    end do
    allocate(list(n))
    n = 0
    do i = 1, 3
       if (lists(cfg % families, trim(every(i)))) then
          n = n + 1
          list(n) = every(i)
       end if
    end do
  end function listed
  subroutine tuple_rows(cfg, arity, printed)
    type(configuration), intent(in)    :: cfg
    integer            , intent(in)    :: arity
    integer            , intent(inout) :: printed
    character(len=8), allocatable :: names(:)
    integer :: which(arity), orders(arity)
    integer :: m, code, k, r, order
    names = listed(cfg)
    m     = size(names)
    do code = 0, m ** arity - 1
       r = code
       do k = arity, 1, -1
          which(k) = mod(r, m) + 1
          r        = r / m
       end do
       if (any([(any(which(1:k-1) == which(k)), k = 2, arity)])) cycle
       if (cfg % mixed_orders .and. arity >= 2) then
          code_of_orders: block
            integer :: oc
            do oc = 0, cfg % max_discretization_order ** arity - 1
               r = oc
               do k = arity, 1, -1
                  orders(k) = mod(r, cfg % max_discretization_order) + 1
                  r         = r / cfg % max_discretization_order
               end do
               call one_row(cfg, names(which), orders, printed)
            end do
          end block code_of_orders
       else
          do order = 1, cfg % max_discretization_order
             orders = order
             call one_row(cfg, names(which), orders, printed)
          end do
       end if
    end do
  end subroutine tuple_rows
  subroutine accounted(cfg)
    type(configuration), intent(in) :: cfg
    character(len=32), allocatable :: wanted(:)
    integer :: i, event
    wanted = worded(cfg % measurements)
    do i = 1, size(wanted)
       event = tally_event_of(trim(wanted(i)))
       call one_measurement(cfg, event)
    end do
    call route_note(cfg)
    call cliff_note(cfg)
  end subroutine accounted
  subroutine route_note(cfg)
    type(configuration), intent(in) :: cfg
    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: m
    real(dp) :: counted, rows
    rows = over_levels(0, 5)
    if (rows <= 0.0_dp) return
    write(*,'(a)') ' '
    write(*,'(a)') '   tangent substitutions per row, the model against the count'
    line = '   model            '
    do m = 1, cfg % max_derivative_degree
       write(cell,'(i14)') expansion_substitutions(1, m)
       line = line // cell
    end do
    write(*,'(a)') line
    line = '   counted          '
    do m = 1, cfg % max_derivative_degree
       counted = over_levels(m, 3)
       write(cell,'(f14.2)') counted / rows
       line = line // cell
    end do
    write(*,'(a)') line
  end subroutine route_note
  subroutine one_measurement(cfg, event)
    type(configuration), intent(in) :: cfg
    integer            , intent(in) :: event
    real(dp), allocatable :: whole(:)
    character(len=:), allocatable :: line
    integer :: level, m, top
    top = cfg % max_derivative_degree
    allocate(whole(0:top), source=0.0_dp)
    do m = 0, top
       if (event == wall_time) then
          whole(m) = tally_amount(at_horizon, m, event)
       else
          whole(m) = over_levels(m, event)
       end if
    end do
    write(*,'(a)') ' '
    write(*,'(a)') ' accounting: ' // tally_event_name(event)
    if (event == wall_time) then
       write(*,'(a)') '   seconds. A level holds the levels opened inside it, and the'
       write(*,'(a)') '   expansion spans every order, so the ratios are the horizon.'
    end if
    write(*,'(a)') ' '
    line = '   at each level    '
    do m = 0, top
       line = line // right(order_named(m))
    end do
    write(*,'(a)') line
    do level = 1, tally_num_levels()
       line = '   ' // tally_level_name(level) // &
            & repeat(' ', max(1, 18 - len(tally_level_name(level))))
       do m = 0, top
          line = line // right(amount_text(tally_amount(level, m, event), event))
       end do
       write(*,'(a)') line
    end do
    if (event == wall_time) then
       line = '   per order        '
    else
       line = '   whole run        '
    end if
    do m = 0, top
       line = line // right(amount_text(whole(m), event))
    end do
    write(*,'(a)') line
    call ratio_matrix(whole, top)
  end subroutine one_measurement
  subroutine ratio_matrix(whole, top)
    real(dp), intent(in) :: whole(0:)
    integer , intent(in) :: top
    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: i, j
    write(*,'(a)') ' '
    line = '   row over column  '
    do j = 0, top
       line = line // right(order_named(j))
    end do
    write(*,'(a)') line
    do i = 0, top
       line = '   ' // order_named(i) // repeat(' ', max(1, 18 - len(order_named(i))))
       do j = 0, top
          if (whole(j) > 0.0_dp) then
             write(cell,'(f14.2)') whole(i) / whole(j)
          else
             write(cell,'(a14)') '-'
          end if
          line = line // cell
       end do
       write(*,'(a)') line
    end do
  end subroutine ratio_matrix
  subroutine cliff_note(cfg)
    type(configuration), intent(in) :: cfg
    real(dp) :: loops, solves
    integer  :: m
    loops  = 0.0_dp
    solves = 0.0_dp
    do m = 0, cfg % max_derivative_degree
       loops  = loops  + over_levels(m, 2)
       solves = solves + over_levels(m, 5)
    end do
    if (solves <= 0.0_dp) return
    write(*,'(a)') ' '
    write(*,'(a,f8.1)') '   primal loops per newton solve      ', loops / solves
    if (loops / solves >= 39.0_dp) then
       write(*,'(a)') '   at the iteration budget: the march is not converging, so'
       write(*,'(a)') '   these ratios report the budget and not the derivative order.'
    end if
  end subroutine cliff_note
  function amount_text(spent, event) result(text)
    real(dp), intent(in) :: spent
    integer , intent(in) :: event
    character(len=:), allocatable :: text
    character(len=14) :: cell
    if (event == wall_time) then
       write(cell,'(f14.4)') spent
    else
       write(cell,'(i14)') nint(spent)
    end if
    text = trim(adjustl(cell))
  end function amount_text
  real(dp) function over_levels(m, event) result(total)
    integer, intent(in) :: m, event
    integer :: level
    total = 0.0_dp
    do level = 1, tally_num_levels()
       total = total + tally_amount(level, m, event)
    end do
  end function over_levels
  function order_named(m) result(named)
    integer, intent(in) :: m
    character(len=:), allocatable :: named
    character(len=2) :: digit
    write(digit,'(i0)') m
    if (m == 0) then
       named = 'f'
    else if (m == 1) then
       named = 'dfdx'
    else
       named = 'd' // trim(digit) // 'fdx' // trim(digit)
    end if
  end function order_named
  function right(text) result(cell)
    character(len=*), intent(in) :: text
    character(len=14) :: cell
    write(cell,'(a14)') text
  end function right
end program graph_time_integrator
