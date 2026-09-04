module gti_configuration
  use util_precision  , only : dp
  implicit none
  private
  public :: configuration, read_configuration, override, show
  public :: words_of, lists, refuse_unknown, chosen_from, argument_values
  ! WHICH ARGUMENTS SELECT WHAT TO RUN RATHER THAN HOW. --config=name
  ! names the configuration file; --demo=name and --list-demos select
  ! a demonstration, the latter valued list; every other argument is
  ! a setting, or a demonstration's own, and is its own value.
  integer, parameter, public :: names_setting = 1
  integer, parameter, public :: names_config  = 2
  integer, parameter, public :: names_demo    = 3
  ! the levels of this application's hierarchy, outermost first: the
  ! level an accounting count is recorded under
  integer, parameter, public :: at_expansion = 1
  integer, parameter, public :: at_horizon   = 2
  integer, parameter, public :: at_block     = 3
  integer, parameter, public :: at_stage     = 4
  character(len=9), parameter, public :: hierarchy_levels(4) = &
       & ['expansion', 'horizon  ', 'block    ', 'stage    ']
  type :: configuration
     character(len=32) :: physics      = 'vanderpol'
     character(len=32) :: grid         = 'random'
     character(len=32) :: combinations = '1'

     ! HOW AN ORDER OF ACCURACY IS ESTIMATED FROM REFINED GRIDS. The coarsest
     ! grid has coarsest_instants per window; each grid after it
     ! has refinement_ratio times the previous count, rounded to a whole
     ! number of instants; and all of them are compared against a
     ! reference grid refined reference_refinement times beyond the finest.
     !
     ! THE RATIO NEED NOT BE TWO. The same order estimated at two
     ! different ratios is evidence that the estimate is asymptotic
     ! rather than dependent on the particular grid spacings.
     real(dp) :: refinement_ratio     = 2.0_dp
     integer  :: refinement_grids     = 4
     integer  :: coarsest_instants    = 10
     real(dp) :: reference_refinement = 4.0_dp

     ! HOW FAR AN ESTIMATE MAY DEVIATE AND STILL MATCH THE ORDER, and how far
     ! the pairwise estimates may differ and still count as one power of h.
     ! Grids coarse enough to compute at low cost do not determine a high order to two
     ! decimal places, and both tolerances state how much deviation is
     ! accepted.
     real(dp) :: order_tolerance      = 0.5_dp
     real(dp) :: spread_tolerance     = 1.0_dp

     ! ONE CHAIN, STATED EXPLICITLY. A window per word, each a family and
     ! the order required of it, so `bdf:2 dirk:3 adams:3 bdf:2` is a
     ! chain of four windows. Any length is accepted, and a family may
     ! occupy more than one window - neither of which a survey over
     ! window counts can express. Empty means no such chain is
     ! requested, which is the default.
     character(len=256) :: chain = ''
     character(len=32) :: families     = 'bdf adams dirk'
     integer  :: state_degree             = 2
     integer  :: instants                 = 21
     integer  :: max_derivative_degree    = 3
     integer  :: max_discretization_order = 4
     integer  :: seed                     = 20260824
     integer  :: startup_refinement       = 4
     character(len=16) :: linear_solver   = 'direct'
     character(len=16) :: jacobian        = 'matrix'

     ! THE ROWS THE LINEAR SOLVER READS. Each kind named here is assembled
     ! as rows of its own and its values are unknowns of the
     ! solve. A kind left out has been eliminated, and elimination
     ! states how. The states are always assembled.
     character(len=64) :: rows            = 'states state-time-derivatives'

     ! How a kind left out of rows was removed. symbolic substitutes the
     ! equation out as the system is formed; numerical assembles its rows
     ! and eliminates them before the solve. Both solve the same
     ! equations.
     character(len=16) :: elimination     = 'symbolic'
     character(len=16) :: storage         = 'dense'
     logical           :: multigrid       = .false.
     character(len=16) :: space           = 'coupled'
     character(len=16) :: time            = 'sequential'
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
     character(len=16) :: adaptive_check       = 'step_doubling'
     character(len=16) :: iteration_criterion = 'by_rate'
     integer           :: max_iterations      = 100
     integer           :: higher_order_jacobian_product = 1
     integer           :: krylov_restart        = 60
     integer           :: smoothing_sweeps      = 2
     integer           :: max_linear_iterations = 200
     logical  :: accounting                   = .false.
     character(len=256) :: measurements       = &
          & 'elapsed_time primal_loops tangent_loops adjoint_loops ' // &
          & 'newton_solves linear_solves factorisations'
  end type configuration
contains
  pure function levelled(given) result(name)
    character(len=*), intent(in) :: given
    character(len=:), allocatable :: name
    integer :: i, c
    name = trim(adjustl(given))
    do i = 1, len(name)
       if (name(i:i) == '-') name(i:i) = '_'
       c = iachar(name(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) name(i:i) = achar(c + 32)
    end do
  end function levelled
  function argument_values(kinds) result(values)
    integer, intent(in) :: kinds(:)
    character(len=256), allocatable :: values(:)
    character(len=256) :: argument, found(command_argument_count())
    integer :: i, n, kind
    n = 0
    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(argument, '--config=') == 1) then
          kind     = names_config
          argument = argument(10:)
       else if (index(argument, '--demo=') == 1) then
          kind     = names_demo
          argument = argument(8:)
       else if (trim(argument) == '--list-demos') then
          kind     = names_demo
          argument = 'list'
       else
          kind = names_setting
       end if
       if (.not. any(kinds == kind)) cycle
       n        = n + 1
       found(n) = argument
    end do
    values = found(1:n)
  end function argument_values
  pure function words_of(phrase) result(list)
    character(len=*), intent(in) :: phrase
    character(len=32), allocatable :: list(:)
    character(len=32) :: words(32)
    integer :: i, first, n, last
    n    = 0
    i    = 1
    last = len_trim(phrase)
    do while (i <= last)
       if (phrase(i:i) == ' ') then
          i = i + 1
          cycle
       end if
       first = i
       do while (i <= last)
          if (phrase(i:i) == ' ') exit
          i = i + 1
       end do
       if (n == size(words)) exit
       n = n + 1
       words(n) = phrase(first:i-1)
    end do
    list = words(1:n)
  end function words_of
  pure logical function lists(phrase, what) result(yes)
    character(len=*), intent(in) :: phrase, what
    character(len=32), allocatable :: list(:)
    integer :: i
    list = words_of(phrase)
    yes  = .false.
    do i = 1, size(list)
       if (trim(list(i)) == what) yes = .true.
    end do
  end function lists
  subroutine refuse_unknown(given, every, subject, which)
    character(len=*), intent(in) :: given, every(:), subject
    integer, intent(out), optional :: which
    character(len=32), allocatable :: list(:)
    integer :: i, j, at
    list = words_of(given)
    at   = 0
    do i = 1, size(list)
       at = 0
       do j = 1, size(every)
          if (trim(list(i)) == trim(every(j))) at = j
       end do
       if (at == 0) then
          write(*,'(a)') ' '
          write(*,'(a)') ' ' // subject // ' names ' // trim(list(i)) // &
               & ', which this program does not define.'
          error stop 'gti_configuration: a setting names something unknown'
       end if
    end do
    if (present(which)) which = at
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
    case ('chain')
       cfg % chain = value
    case ('refinement_ratio')
       read(value, *) cfg % refinement_ratio
    case ('refinement_grids')
       read(value, *) cfg % refinement_grids
    case ('coarsest_instants')
       read(value, *) cfg % coarsest_instants
    case ('reference_refinement')
       read(value, *) cfg % reference_refinement
    case ('order_tolerance')
       read(value, *) cfg % order_tolerance
    case ('spread_tolerance')
       read(value, *) cfg % spread_tolerance
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
    case ('space')
       cfg % space = value
    case ('time')
       cfg % time = value
    case ('designs')
       cfg % designs = value
    case ('functionals')
       cfg % functionals = value
    case ('jacobian')
       cfg % jacobian = value
    case ('rows')
       cfg % rows = value
    case ('elimination')
       cfg % elimination = value
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
    case ('adaptive_check')
       cfg % adaptive_check = value
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
         & 'max_iterations', 'an iteration limit is at least one')
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
    error stop 'gti_configuration: a setting is below its least meaningful value'
  end subroutine refuse_below
  subroutine read_line(cfg, line)
    type(configuration), intent(inout) :: cfg
    character(len=*)   , intent(in)    :: line
    character(len=:), allocatable :: setting
    integer :: at
    setting = trim(adjustl(line))
    at   = index(setting, '#')
    if (at > 0) setting = trim(setting(:at - 1))
    if (len(setting) == 0) return
    at = index(setting, '=')
    if (at == 0) then
       write(*,'(a)') ' this line names no setting: ' // setting
       error stop 'gti_configuration: a setting is named, then given its value'
    end if
    call assign(cfg, setting(:at - 1), trim(adjustl(setting(at + 1:))))
  end subroutine read_line
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
       call read_line(cfg, line)
    end do
    close(unit)
  end subroutine read_configuration
  subroutine override(cfg, argument)
    type(configuration), intent(inout) :: cfg
    character(len=*)   , intent(in)    :: argument
    character(len=:), allocatable :: setting
    integer :: at
    setting = trim(adjustl(argument))
    if (len(setting) > 2) then
       if (setting(1:2) == '--') setting = setting(3:)
    end if
    at = index(setting, '=')
    if (at == 0) then
       write(*,'(a)') ' this argument names no setting: ' // setting
       error stop 'gti_configuration: an argument is a setting and its value'
    end if
    call assign(cfg, setting(:at - 1), trim(adjustl(setting(at + 1:))))
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
    write(*,'(a,a)')       '   space                    ', trim(cfg % space)
    write(*,'(a,a)')       '   time                     ', trim(cfg % time)
    write(*,'(a,a)')       '   designs                  ', trim(cfg % designs)
    write(*,'(a,a)')       '   functionals              ', trim(cfg % functionals)
    write(*,'(a,a)')       '   jacobian                 ', trim(cfg % jacobian)
    write(*,'(a,a)')       '   rows                     ', trim(cfg % rows)
    write(*,'(a,a)')       '   elimination              ', trim(cfg % elimination)
    write(*,'(a,a)')       '   storage                  ', trim(cfg % storage)
    write(*,'(a,l1)')      '   multigrid                ', cfg % multigrid
    write(*,'(a,f0.4)')    '   design                   ', cfg % design
    write(*,'(a,i0)')      '   max derivative degree    ', cfg % max_derivative_degree
    write(*,'(a,i0)')      '   max discretization order ', cfg % max_discretization_order
    write(*,'(a,a)')       '   families                 ', trim(cfg % families)
    write(*,'(a,a)')       '   combinations             ', trim(cfg % combinations)
    if (len_trim(cfg % chain) > 0) then
       write(*,'(a,a)')    '   chain                    ', trim(cfg % chain)
    end if
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
    if (trim(cfg % grid) == 'adaptive') then
       write(*,'(a,a)')    '   adaptive check           ', trim(cfg % adaptive_check)
       if (trim(cfg % adaptive_check) == 'goal_oriented') then
          write(*,'(a)')   '     targeting van der Pol energy'
       end if
    end if
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
  integer function chosen_from(name, every, subject) result(which)
    character(len=*), intent(in) :: name, every(:), subject
    call refuse_unknown(name, every, subject, which)
  end function chosen_from
end module gti_configuration

!=====================================================================!
! The equations this application integrates, each stated once as an
! expression: the van der Pol residual at any degree, its energy, and
! the power its damping dissipates. The expression module defines none of
! them; they are this application's physics, so they are defined here.
!=====================================================================!

module gti_physics
  use util_precision      , only : dp
  use operation_expression, only : expression, unknown, design, derivative, stated, stated_over, &
       & euler_lagrange, at_zero, &
       & FIRST_COORDINATE, &
       & operator(+), operator(-), operator(*), operator(**)
  implicit none
  private
  public :: van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  ! the multiplier field of the van der pol Lagrangian: the costate
  integer, parameter :: COSTATE = 2
contains
  !===================================================================!
  ! The van der pol residual in the state alone, unstated.
  !===================================================================!
  function residual_rule(degree) result(r)
    integer, intent(in) :: degree
    type(expression) :: r
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    r = derivative(q, degree) - nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, degree - 1) &
         & + derivative(q, 0)
  end function residual_rule
  !===================================================================!
  ! The Lagrangian L = F + lambda R over the state and the costate,
  ! for a functional F in the state. The residual is its stationarity
  ! in the costate and F is the Lagrangian at zero costate.
  !===================================================================!
  function lagrangian(functional, degree, label) result(l)
    type(expression), intent(in) :: functional
    integer         , intent(in) :: degree
    character(len=*), intent(in) :: label
    type(expression) :: l
    l = stated(functional + unknown(COSTATE) * residual_rule(degree), degree, label)
  end function lagrangian
  function energy_rule() result(f)
    type(expression) :: f
    type(expression) :: q
    q = unknown()
    f = 0.5_dp * (derivative(q, 0)**2 + derivative(q, 1)**2)
  end function energy_rule
  function dissipation_rule() result(f)
    type(expression) :: f
    type(expression) :: q, nu
    q  = unknown()
    nu = design()
    f = nu * (1.0_dp - derivative(q, 0)**2) * derivative(q, 1) * derivative(q, 1)
  end function dissipation_rule
  function van_der_pol(degree) result(r)
    integer, intent(in) :: degree
    type(expression) :: r
    r = euler_lagrange(lagrangian(energy_rule(), degree, 'van der pol lagrangian'), COSTATE, &
         & 'van der pol residual')
  end function van_der_pol
  function van_der_pol_energy(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    f = at_zero(lagrangian(energy_rule(), degree, 'van der pol lagrangian'), COSTATE, 'van der pol energy')
  end function van_der_pol_energy
  function van_der_pol_dissipation(degree) result(f)
    integer, intent(in) :: degree
    type(expression) :: f
    f = at_zero(lagrangian(dissipation_rule(), degree, 'van der pol lagrangian'), COSTATE, &
         & 'van der pol dissipation')
  end function van_der_pol_dissipation
end module gti_physics
module gti_sweeps
  use gti_configuration, only : refuse_unknown, lists
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
  public :: forward_pass, reverse_pass, pass_of, pass_substitutions, choose
  integer, parameter :: forward_pass = 1
  integer, parameter :: reverse_pass = 2
  public :: set_linear_solver, set_jacobian, set_storage, set_multigrid
  public :: set_rows, set_elimination
  public :: set_newton_order, newton_order
  public :: set_aggregates, set_coarse_nodes, coarse_nodes, jacobian_present, multigrid_on
  public :: read_inner, store_inner, clear_inner, set_linear_stopping, set_linear_budget
  public :: stopping_applied
  public :: functional_of, functional_gradient
  character(len=16), save :: chosen_solver   = 'direct'
  character(len=16), save :: chosen_jacobian = 'matrix'
  character(len=64), save :: chosen_rows        = 'states state-time-derivatives'
  character(len=16), save :: chosen_elimination = 'symbolic'
  integer          , save :: chosen_newton_order = 1
  character(len=16), save :: chosen_storage  = 'dense'
  logical          , save :: chosen_multigrid = .false.
  integer, allocatable, save :: chosen_aggregates(:)
  real(dp), save :: linear_tolerance  = half_digits
  integer , save :: linear_criterion  = relative
  integer , save :: linear_limit_kind     = by_rate
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
    call clear_inner()
  end subroutine set_linear_solver
  subroutine set_rows(name)
    character(len=*), intent(in) :: name
    logical :: states, in_time, in_space
    call refuse_unknown(name, ['states                   ', &
         &                      'state-time-derivatives   ', &
         &                      'state-spatial-derivatives'], 'rows')
    states   = lists(name, 'states')
    in_time  = lists(name, 'state-time-derivatives')
    in_space = lists(name, 'state-spatial-derivatives')
    if (.not. states) then
       error stop 'gti_sweeps: the states are always assembled; rows names them'
    end if
    ! The scheme states one row per derived degree, so the time
    ! derivatives are assembled in the current implementation and cannot yet be left out. The
    ! spatial law is substituted into the residual, so its derivatives have
    ! no rows to assemble yet.
    if (.not. in_time) then
       error stop 'gti_sweeps: leaving the time derivatives out of rows is not implemented; &
            &the scheme assembles one row for each of them'
    end if
    if (in_space) then
       error stop 'gti_sweeps: assembling rows for the spatial derivatives is not implemented; &
            &the spatial law is stored in the state row'
    end if
    chosen_rows = name
    call clear_inner()
  end subroutine set_rows
  subroutine set_elimination(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['symbolic ', 'numerical'], 'elimination')
    ! A numerical elimination assembles the rows it then removes, and
    ! no kind left out of rows has rows to assemble yet.
    if (trim(name) == 'numerical') then
       error stop 'gti_sweeps: a numerical elimination assembles the rows it removes, &
            &and no kind final_imbalance out of rows is assembled yet'
    end if
    chosen_elimination = name
    call clear_inner()
  end subroutine set_elimination
  pure function rows_named() result(name)
    character(len=:), allocatable :: name
    name = trim(chosen_rows) // ', eliminated ' // trim(chosen_elimination)
  end function rows_named
  subroutine set_jacobian(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['matrix', 'free  '], 'jacobian')
    chosen_jacobian = name
    call clear_inner()
  end subroutine set_jacobian
  subroutine set_storage(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['dense ', 'sparse'], 'storage')
    chosen_storage = name
    call clear_inner()
  end subroutine set_storage
  subroutine set_multigrid(on)
    logical, intent(in) :: on
    chosen_multigrid = on
    call clear_inner()
  end subroutine set_multigrid
  subroutine set_newton_order(order)
    integer, intent(in) :: order
    if (order < 1) error stop 'gti_sweeps: one is Newton, and there is no order below it'
    chosen_newton_order = order
  end subroutine set_newton_order
  pure integer function newton_order() result(order)
    order = chosen_newton_order
  end function newton_order
  pure logical function jacobian_present() result(yes)
    yes = trim(chosen_jacobian) == 'matrix'
  end function jacobian_present
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
  subroutine set_linear_stopping(tolerance, criterion, limit_kind)
    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, limit_kind
    if (tolerance <= 0.0_dp) error stop 'gti_sweeps: a tolerance is positive'
    if (criterion /= relative .and. criterion /= absolute) then
       error stop 'gti_sweeps: a tolerance is measured relative or absolute'
    end if
    if (limit_kind /= by_count .and. limit_kind /= by_rate) then
       error stop 'gti_sweeps: an iteration limit is counted or taken from the rate'
    end if
    linear_tolerance = tolerance
    linear_criterion = criterion
    linear_limit_kind    = limit_kind
  end subroutine set_linear_stopping
  subroutine set_linear_budget(restart, sweeps, iterations)
    integer, intent(in) :: restart, sweeps, iterations
    if (restart < 1 .or. sweeps < 1 .or. iterations < 1) then
       error stop 'gti_sweeps: a restart, a sweep count and an iteration limit are positive'
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
    if (trim(chosen_jacobian) == 'free' .and. trim(chosen_solver) == 'direct') then
       error stop 'gti_sweeps: a matrix-free jacobian has no matrix to factorise; its solver iterates'
    end if
    if (trim(chosen_solver) == 'direct' .and. trim(chosen_storage) == 'sparse') then
       error stop 'gti_sweeps: a sparse direct solve is not implemented'
    end if
    if (trim(chosen_solver) == 'iterative' .and. trim(chosen_jacobian) == 'matrix' &
         & .and. trim(chosen_storage) == 'dense') then
       error stop 'gti_sweeps: an iterative solve reads the sparse stencil; dense storage is for factorising'
    end if
    if (chosen_multigrid .and. trim(chosen_jacobian) == 'free') then
       error stop 'gti_sweeps: multigrid coarsens a stencil, which a matrix-free jacobian does not store'
    end if
    select case (trim(chosen_solver))
    case ('direct')
       factorisation = dense_direct()
       factorisation % singular_reported = .true.
       allocate(named, source=factorisation)
    case ('iterative')
       krylov = gmres()
       krylov % restart = min(count, linear_restart)
       call stopping_applied(krylov, linear_tolerance, linear_criterion, linear_limit_kind, &
            & linear_iterations)
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
    call stopping_applied(levels, linear_tolerance, linear_criterion, linear_limit_kind, &
         & linear_iterations)
    allocate(inner, source=levels)
  end function inner_minimizer
  subroutine stopping_applied(m, tolerance, criterion, limit_kind, iterations)
    class(minimizer), intent(inout) :: m
    real(dp)        , intent(in)    :: tolerance
    integer         , intent(in)    :: criterion, limit_kind, iterations
    m % tolerance      = tolerance
    m % criterion      = criterion
    m % limit_kind     = limit_kind
    m % max_iterations = iterations
  end subroutine stopping_applied
  subroutine read_inner(inner, count, width)
    class(minimizer), allocatable, intent(out) :: inner
    integer                      , intent(in)  :: count, width
    if (allocated(kept_inner) .and. .not. chosen_multigrid) then
       call move_alloc(kept_inner, inner)
    else
       call clear_inner()
       allocate(inner, source=inner_minimizer(count, width))
    end if
  end subroutine read_inner
  subroutine store_inner(inner)
    class(minimizer), allocatable, intent(inout) :: inner
    if (allocated(kept_inner)) deallocate(kept_inner)
    call move_alloc(inner, kept_inner)
  end subroutine store_inner
  subroutine clear_inner()
    if (allocated(kept_inner)) deallocate(kept_inner)
  end subroutine clear_inner
  pure integer function pass_of(num_designs, num_functionals, order) result(pass_kind)
    integer, intent(in) :: num_designs, num_functionals, order
    if (num_designs < 1 .or. num_functionals < 1 .or. order < 1) then
       error stop 'gti_sweeps: a differentiation pass is chosen for at least one design, one functional and order one'
    end if
    if (num_designs <= order * num_functionals) then
       pass_kind = forward_pass
    else
       pass_kind = reverse_pass
    end if
  end function pass_of
  pure integer function pass_substitutions(pass_kind, num_designs, num_functionals, order) &
       & result(count)
    integer, intent(in) :: pass_kind, num_designs, num_functionals, order
    select case (pass_kind)
    case (forward_pass)
       count = choose(num_designs + order - 1, order)
    case (reverse_pass)
       count = (1 + num_functionals) * choose(num_designs + order - 2, order - 1)
    case default
       error stop 'gti_sweeps: a pass is forward or reverse'
    end select
  end function pass_substitutions
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
  use view_sequence         , only : sequence_num_elements, sequence_element
  use view_relational       , only : relational_binding, relational_valid, &
       & num_relations, relation_at
  use relation_finitary     , only : relation
  use relation_binary       , only : csr_relation, binary_relation
  use map_value             , only : value_map, VALUE_UNKNOWN, VALUE_KNOWN
  use map_label             , only : label_map
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation
  use view_directed_stored  , only : stored_directed_graph
  use view_directed_connectivity, only : connectivity_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field, typed_field_domain
  use operation_family      , only : family
  use operation_grid        , only : grid, partition
  use operation_coupling    , only : weights_of
  use operation_stencil     , only : stencil, combine_triples
  use operation_action      , only : variation
  use operation_weight      , only : scheme_weight
  use operation_expression     , only : expression
  use operation_domain      , only : continuous_domain
  implicit none
  private
  public :: expansion, family_container
  public :: design_of_physics, design_of_steps
  public :: marches_by_stages
  integer, parameter :: design_of_physics = 1
  integer, parameter :: design_of_steps   = 2
  type :: family_container
     class(family), allocatable :: scheme
  end type family_container
  type :: expansion
     type(level_storage)     , private :: nodes
     type(label_map)         , private :: labels
     type(value_map)         , private :: values
     type(set_map)           , private :: extents
     type(relational_binding), private :: bindings
     integer                 , private :: root_at = 0
     ! degrees is the EQUATION'S degree count, which is what every
     ! scheme query reads. spatial_degrees counts the components the
     ! spatial law determines. The layout stride is their sum, and is
     ! requested by name so the two are never confused.
     integer                 , private :: degrees = 0
     integer                 , private :: spatial_degrees = 0
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
     procedure :: stride
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
  function label_of(this, g) result(label)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    character(len=:), allocatable :: label
    label = ''
    if (this % labels % labelled(g)) label = this % labels % label_of(g)
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
       error stop 'gti_expansion: a coupling contains one relation'
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
    type(family_container)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    class(grid)           , intent(in)    :: steps
    integer               , intent(in)    :: max_derivative_degree
    real(dp)              , intent(in)    :: parameter
    integer      , intent(in), optional   :: nodes
    type(stencil), intent(in), optional   :: spatial_discretization_stencil
    real(dp)     , intent(in), optional   :: weights(:), block_steps(:)
    real(dp), allocatable :: dt(:)
    type(continuous_domain) :: continuous
    integer , allocatable :: sweeps(:)
    integer :: s
    if (this % root_at /= 0) then
       error stop 'gti_expansion: an expansion is built once'
    end if
    if (size(schemes) /= size(instants)) then
       error stop 'gti_expansion: one family and one instant count per block'
    end if
    continuous = continuous_domain(physics)
    ! degrees is the marching coordinate's own count, which is what
    ! every scheme query reads; the rest of the point's components
    ! belong to the other coordinates
    this % degrees         = continuous % equation_degree() + 1
    this % spatial_degrees = continuous % num_components() - this % degrees
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
  subroutine one_design(this, label, x, kind)
    class(expansion), intent(inout) :: this
    character(len=*), intent(in)    :: label
    real(dp)        , intent(in)    :: x(:)
    integer         , intent(in)    :: kind
    integer :: at
    at = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(at), label)
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
    type(stored_field) :: designs
    type(stored_field), allocatable :: direction(:)
    type(typed_field_domain) :: instant_scalars
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
    instant_scalars = typed_field_domain(instants % vertex_set(), size(weights))
    designs  = instant_scalars % design(weights)
    allocate(e(size(weights)), direction(size(weights_varied)), variations(size(weights_varied)))
    do i = 1, size(weights_varied)
       e = 0.0_dp
       e(weights_varied(i)) = 1.0_dp
       direction(i) = instant_scalars % direction(e)
       variations(i) = variation(this % steps_kept % argument(1), direction(i))
    end do
    call this % steps_kept % partial_action(instants, this % steps_kept % bind([designs]), variations, out)
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
    error stop 'gti_expansion: the steps of this expansion are not designs'
  end subroutine weights_of_steps
  integer function spatial_discretization_coupling(this, spatial_discretization_stencil) result(at)
    class(expansion), intent(inout) :: this
    type(stencil)   , intent(in)    :: spatial_discretization_stencil
    integer , allocatable :: table(:,:), heads(:), tails(:), rows(:), cols(:)
    real(dp), allocatable :: given(:), w(:)
    integer :: e, ne
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
    at = coupled(this, [integer ::], this % node_extent, &
         & 'the nodes the spatial discretization stencil reads', &
         & 'the nodes whose rows the spatial discretization stencil enters', &
         & 'the spatial discretization stencil''s connectivity', 'the spatial discretization stencil', table, w)
  end function spatial_discretization_coupling
  !===================================================================!
  ! A COUPLING OVER A CONNECTIVITY: the set the connectivity reads and
  ! the set whose rows it enters, each of extent n and labelled; the
  ! connectivity, labelled connectivity_label, bound as a csr relation
  ! over the table; the coupling of the members with both sets,
  ! labelled coupling_label, storing the weights in the relation's tuple order.
  !===================================================================!
  integer function coupled(this, members, n, read_label, entered_label, connectivity_label, coupling_label, table, w) &
       & result(at)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: members(:), n, table(:,:)
    character(len=*), intent(in)    :: read_label, entered_label, connectivity_label, coupling_label
    real(dp)        , intent(in)    :: w(:)
    integer :: from, into, owner
    from  = named_set(this, n, read_label)
    into  = named_set(this, n, entered_label)
    owner = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(owner), connectivity_label)
    at = this % nodes % couple([members, from, into], [owner])
    call bind_carriers(this, [members, from, into])
    call bind_coupling(this, owner, from, into, table)
    call this % labels % bind(this % node(at), coupling_label)
    call attach_known(this, at, in_relation_order(this, this % node(at), table, w))
  end function coupled
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
    type(family_container)   , intent(in)    :: schemes(:)
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
    type(family_container)   , intent(in)    :: schemes(:)
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
       error stop 'gti_expansion: a block contains more instants than its family''s history depth'
    end if
    allocate(slices(last - first + 1))
    do k = first, last
       slices(k - first + 1) = one_slice(this, physics, scheme, k, first, dt(k))
    end do
    if (marches_by_stages(scheme, this % degrees)) then
       coupling = add_coupling(this, scheme, slices, first, last, dt)
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
       at = unstaged_slice(this, scheme, instant, first)
    end if
    call this % labels % bind(this % node(at), 'slice at instant ' // written(instant))
  end function one_slice
  integer function unstaged_slice(this, scheme, instant, first) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first
    integer, allocatable :: components(:)
    integer :: d
    allocate(components(this % stride()))
    do d = 0, this % stride() - 1
       components(d + 1) = one_component(this, scheme, instant, first, d, .true.)
    end do
    at = this % nodes % assemble(components, 0)
  end function unstaged_slice
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
    allocate(components(this % stride()))
    do d = 0, this % stride() - 1
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
  function written(n) result(name)
    class(*), intent(in) :: n
    character(len=:), allocatable :: name
    character(len=24) :: buffer
    select type (n)
    type is (integer)
       write(buffer,'(i0)') n
    type is (real(dp))
       write(buffer,'(f0.4)') n
    class default
       error stop 'gti_expansion: a label is written from a number'
    end select
    name = trim(buffer)
  end function written
  integer function block_coupling(this, physics, scheme, slices, first, last, dt) result(at)
    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    class(family)         , intent(in)    :: scheme
    integer               , intent(in)    :: slices(:), first, last
    real(dp)              , intent(in)    :: dt(:)
    type(connectivity_graph) :: connectivity
    real(dp), allocatable :: w(:)
    integer :: n, nd
    associate (u1 => physics); end associate
    n  = last - first + 1
    nd = this % degrees
    connectivity = scheme % block_connectivity(nd, n)
    call weights_of(scheme_weight(scheme), connectivity, dt(first:last), w)
    at = coupled(this, slices, n * nd, 'the components of this block', &
         & 'the constraint instances of this block', 'the scheme connectivity', &
         & scheme % name() // ' coupling', tuples(nd, connectivity), w)
  end function block_coupling
  pure function tuples(nd, connectivity) result(table)
    integer                  , intent(in) :: nd
    type(connectivity_graph) , intent(in) :: connectivity
    integer, allocatable :: table(:,:)
    integer :: e, ne
    ne = connectivity % num_edges()
    allocate(table(2, ne))
    table(1,:) = [((connectivity % edge_tail(e) - 1) * nd + connectivity % tail_degree(e) + 1, e = 1, ne)]
    table(2,:) = [((connectivity % edge_head(e) - 1) * nd + connectivity % head_degree(e) + 1, e = 1, ne)]
  end function tuples
  integer function named_set(this, n, label) result(at)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: n
    character(len=*), intent(in)    :: label
    at = this % nodes % assemble([integer ::], 0)
    call this % labels % bind(this % node(at), label)
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
  subroutine bind_coupling(this, owner, components, constraints, table)
    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: owner, components, constraints, table(:,:)
    type(graph), pointer :: g, from, into
    from => this % nodes % node(components)
    into => this % nodes % node(constraints)
    g    => this % nodes % node(owner)
    call this % bindings % bind_relation(g, &
         & csr_relation('scheme coupling', from, into, table, this % extents))
  end subroutine bind_coupling
  recursive logical function consistent(this, g) result(passes_check)
    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g
    type(graph), pointer :: coupling
    type(branch) :: members
    integer :: k
    passes_check = level_consistent(g)
    if (.not. passes_check) return
    if (level_couples(g)) then
       coupling => level_coupling(g)
       passes_check = relational_valid(coupling, this % bindings)
       if (.not. passes_check) return
    end if
    if (level_is_leaf(g)) return
    members = level_members(g)
    do k = 1, sequence_num_elements(members)
       passes_check = this % consistent(sequence_element(members, k))
       if (.not. passes_check) return
    end do
  end function consistent
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
  integer function slice_coupling(this, scheme, members, s, step) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: members(:), s
    real(dp)        , intent(in)    :: step
    type(connectivity_graph) :: connectivity
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: w(:)
    integer :: nd, e, ne
    nd = this % degrees
    connectivity = scheme % stage_connectivity(nd)
    ne = connectivity % num_edges()
    call weights_of(scheme_weight(scheme), connectivity, spread(step, 1, s + 2), w)
    allocate(table(2, ne))
    table(1,:) = [(stage_unknown(connectivity % edge_tail(e), connectivity % tail_degree(e), s, nd), e = 1, ne)]
    table(2,:) = [(stage_unknown(connectivity % edge_head(e), connectivity % head_degree(e), s, nd), e = 1, ne)]
    at = coupled(this, members, (s + 1) * nd, 'the components of this step', &
         & 'the constraint instances of this step', 'the butcher connectivity', &
         & scheme % name() // ' stage coupling', table, w)
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
  subroutine transfer_table(this, s, n, table, sources)
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
  end subroutine transfer_table
  !===================================================================!
  ! THE LAYOUT STRIDE: how many components one node stores at one
  ! moment. The equation's degrees come first and the spatial law's
  ! come after, so a scheme query reads degrees and a layout
  ! query reads this.
  !===================================================================!
  pure integer function stride(this)
    class(expansion), intent(in) :: this
    stride = this % degrees + this % spatial_degrees
  end function stride
  integer function add_coupling(this, scheme, slices, first, last, dt) result(at)
    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: slices(:), first, last
    real(dp)        , intent(in)    :: dt(:)
    integer, allocatable :: table(:,:), sources(:)
    type(connectivity_graph) :: connectivity
    real(dp), allocatable :: w(:)
    integer :: n, nd, s, e
    n  = last - first + 1
    nd = this % degrees
    s  = scheme % num_stages()
    call transfer_table(this, s, n, table, sources)
    connectivity = connectivity_graph(s + 2, [(1, e = 1, size(sources))], sources, &
         & [((mod(table(1, e) - 1, nd)), e = 1, size(sources))], &
         & [((mod(table(2, e) - 1, nd)), e = 1, size(sources))])
    call weights_of(scheme_weight(scheme), connectivity, spread(dt(first), 1, s + 2), w)
    at = coupled(this, slices, (1 + (n - 1) * (s + 1)) * nd, 'the components of this block', &
         & 'the constraint instances of this block', 'the transfer between steps', &
         & scheme % name() // ' transfer coupling', table, w)
  end function add_coupling
end module gti_expansion
module gti_block
  use util_precision  , only : dp
  use operation_action     , only : operation, variation, contract
  use operation_action     , only : binding
  use view_directed        , only : directed_graph
  use view_directed_stored , only : stored_directed_graph
  use field_calculus       , only : field, FIELD_REAL
  use field_stored         , only : stored_field
  use graph_fractal        , only : graph
  use gti_expansion        , only : expansion
  use view_level           , only : level_member, level_num_members, level_is_leaf
  use operation_stencil    , only : combine_triples, stencil
  use operation_residual   , only : residual_operator
  use operation_family     , only : family
  use operation_coupling   , only : matrix_scheme_connectivity, connectivity_terms
  use operation_expression    , only : expression
  use view_directed        , only : forward
  implicit none
  private
  public :: block_residual
  ! the slice, node and moment of every unknown of a block, read from
  ! the tower once when the block is placed
  type :: block_layout
     integer, allocatable :: slice(:), node(:), moment(:)
  end type block_layout
  type, extends(residual_operator) :: block_residual
     type(block_layout)      , private :: layout
     integer, allocatable    , private :: free(:)
     type(matrix_scheme_connectivity), allocatable, private :: connectivity(:)
   contains
     procedure :: constrained_block => block_constrained
     procedure :: linear_block => block_linear_block
     procedure :: placed_on
     procedure :: slice_of
     procedure :: node_of
     procedure :: moment_of
     procedure, private :: labels_of
     procedure :: num_nodes
     procedure :: spatial_discretization_laid
     procedure :: aggregates
     procedure :: with_connectivity
     procedure :: rows_terms
     procedure :: member_order
     procedure :: sweep_labels
  end type block_residual
  interface block_residual
     module procedure create
  end interface block_residual
contains
  function create(derived, physics, at, unknowns, degrees, primary, fixed_rows, fixed, &
       & spatial_discretization_stencil) result(this)
    type(stencil)         , intent(in) :: derived
    type(expression)      , intent(in) :: physics
    integer               , intent(in) :: at(:), unknowns, degrees, primary
    integer               , intent(in) :: fixed_rows(:)
    real(dp)              , intent(in) :: fixed(:)
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
    type(block_residual) :: this
    this % residual_operator = residual_operator(derived, physics, at, unknowns, degrees, primary, &
         & fixed_rows, fixed, spatial_discretization_stencil)
  end function create
  subroutine placed_on(this, tower, node)
    class(block_residual), intent(inout)     :: this
    type(expansion)      , intent(in)        :: tower
    type(graph)          , intent(in), target :: node
    type(graph), pointer :: one_slice, first
    integer :: n, k, j, members, moments, g, m, count, u, i, d, degrees
    degrees = this % num_degrees()
    n = level_num_members(node)
    moments = 0
    do k = 1, n
       moments = moments + members_of(level_member(node, k))
    end do
    first => level_member(node, 1)
    if (.not. level_is_leaf(level_member(first, 1))) first => level_member(first, 1)
    m     = tower % extent_of(level_member(first, 1))
    count = moments * m * degrees
    allocate(this % layout % slice(count), this % layout % node(count), this % layout % moment(count))
    g = 0
    do k = 1, n
       one_slice => level_member(node, k)
       members = members_of(one_slice)
       do j = 1, members
          g = g + 1
          do i = 1, m
             do d = 0, degrees - 1
                u = ((g - 1) * m + (i - 1)) * degrees + d + 1
                this % layout % slice(u)  = k
                this % layout % node(u)   = i
                this % layout % moment(u) = g
             end do
          end do
       end do
    end do
  contains
    integer function members_of(one_slice)
      type(graph), intent(in) :: one_slice
      if (level_is_leaf(level_member(one_slice, 1))) then
         members_of = 1
      else
         members_of = level_num_members(one_slice)
      end if
    end function members_of
  end subroutine placed_on
  subroutine labels_of(this, slice, node, moment)
    class(block_residual), intent(in) :: this
    integer, allocatable , intent(out) :: slice(:), node(:), moment(:)
    if (.not. allocated(this % layout % moment)) then
       error stop 'gti_block: the block has not been placed in the graph'
    end if
    if (allocated(this % free)) then
       slice  = this % layout % slice(this % free)
       node   = this % layout % node(this % free)
       moment = this % layout % moment(this % free)
    else
       slice  = this % layout % slice
       node   = this % layout % node
       moment = this % layout % moment
    end if
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
    if (.not. allocated(this % layout % node)) return
    call this % labels_of(slice, node, moment)
    num_nodes = maxval(node)
  end function num_nodes
  subroutine spatial_discretization_laid(this, spatial_discretization_stencil)
    class(block_residual), intent(inout) :: this
    type(stencil)        , intent(in)    :: spatial_discretization_stencil
    integer , allocatable :: base(:,:), r(:), c(:), slice(:), node(:), moment(:), at(:)
    real(dp), allocatable :: lw(:), fixed(:), w(:)
    integer :: nodes, moments, p, u, e, g, ne, n, rc, cc
    call this % labels_of(slice, node, moment)
    at = this % points_at()
    nodes   = maxval(node)
    moments = maxval(moment)
    if (spatial_discretization_stencil % pattern % num_vertices() /= nodes) then
       error stop 'gti_block: the spatial discretization stencil is a stencil over the nodes'
    end if
    call spatial_discretization_stencil % constants % real_vector(fixed)
    if (any(abs(fixed) > 0.0_dp)) then
       error stop 'gti_block: the spatial discretization stencil has no constant term'
    end if
    call spatial_discretization_stencil % weights % real_vector(lw)
    allocate(base(nodes, moments), source=-1)
    do p = 1, size(at)
       u = at(p) + 1
       base(node(u), moment(u)) = at(p)
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
             error stop 'gti_block: a moment contains every node or none'
          end if
          n    = n + 1
          r(n) = base(rc, g) + this % primary_degree() + 1
          c(n) = base(cc, g) + 1
          w(n) = lw(e)
       end do
    end do
    call this % attach_connected_stencil(stencil(r(1:n), c(1:n), w(1:n), spread(0.0_dp, 1, this % num_unknowns()), &
         & 'spatial discretization stencil'))
  end subroutine spatial_discretization_laid
  subroutine with_connectivity(this, connectivity)
    class(block_residual), intent(inout) :: this
    type(matrix_scheme_connectivity) , intent(in)    :: connectivity(:)
    this % connectivity = connectivity
  end subroutine with_connectivity
  subroutine rows_terms(this, scheme, dt, seeds, r, c, w)
    class(block_residual), intent(in) :: this
    class(family)        , intent(in) :: scheme
    real(dp)             , intent(in) :: dt(:), seeds(:,:)
    integer , allocatable, intent(out) :: r(:), c(:)
    real(dp), allocatable, intent(out) :: w(:,:)
    if (.not. allocated(this % connectivity)) then
       error stop 'gti_block: the block was built without its connectivity'
    end if
    call connectivity_terms(scheme, this % connectivity, this % num_nodes(), this % num_degrees(), dt, seeds, r, c, w)
  end subroutine rows_terms
  function aggregates(this, cell) result(aggregate)
    class(block_residual), intent(in) :: this
    integer              , intent(in) :: cell(:)
    integer, allocatable :: aggregate(:)
    integer, allocatable :: numbered(:), slice(:), node(:), moment(:)
    integer :: u, coarse, key, count, unknowns, degrees
    unknowns = this % num_unknowns()
    degrees  = this % num_degrees()
    call this % labels_of(slice, node, moment)
    if (size(cell) < maxval(node)) then
       error stop 'gti_block: a coarse cell for every node'
    end if
    coarse = maxval(cell)
    allocate(aggregate(unknowns))
    allocate(numbered(maxval(moment) * coarse * degrees), source=0)
    count = 0
    do u = 1, unknowns
       key = ((moment(u) - 1) * coarse + cell(node(u)) - 1) * degrees &
            & + mod(u - 1, degrees) + 1
       if (numbered(key) == 0) then
          count         = count + 1
          numbered(key) = count
       end if
       aggregate(u) = numbered(key)
    end do
  end function aggregates
  !===================================================================!
  ! THE EXPLICIT TANGENT, FROZEN INTO A BLOCK RESIDUAL. this %
  ! residual_operator % linearize (src/operation_residual.f90) is
  ! Element's own operation (Ch. 4.6.3 of the dissertation); a block
  ! residual's own layout and free pass across unchanged, since
  ! freezing the tangent changes no point and drops no unknown.
  !===================================================================!
  function block_linear_block(this, input_graph, inputs, rhs, transposed, mark) result(lin)
    class(block_residual), intent(in) :: this
    class(directed_graph), intent(in) :: input_graph
    type(binding)         , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: rhs(:)
    logical              , intent(in) :: transposed
    integer              , intent(in) :: mark
    type(block_residual) :: lin
    lin % residual_operator = this % residual_operator % linearize(input_graph, inputs, rhs, transposed, mark)
    lin % layout = this % layout
    if (allocated(this % free)) lin % free = this % free
  end function block_linear_block
  !===================================================================!
  ! THE BLOCK CONSTRAINED TO A SUBSET OF ITS OWN UNKNOWNS. this %
  ! residual_operator % constrain (src/operation_residual.f90) is
  ! Element's own operation (Ch. 4.6.3 of the dissertation) - it
  ! eliminates unknowns from stencils, points and fixed rows alone.
  ! A block residual's own operation additionally passes its layout
  ! across unchanged and composes its own free with the constraint
  ! just applied, so the two are named separately (constrained_block
  ! here, constrain on residual_operator) rather than sharing one
  ! name for two different operations.
  !===================================================================!
  function block_constrained(this, free, values) result(sub)
    class(block_residual), intent(in) :: this
    integer              , intent(in) :: free(:)
    real(dp)             , intent(in) :: values(:)
    type(block_residual) :: sub
    sub % residual_operator = this % residual_operator % constrain(free, values)
    sub % layout = this % layout
    if (allocated(this % free)) then
       sub % free = this % free(free)
    else
       sub % free = free
    end if
  end function block_constrained
  !===================================================================!
  ! THE LABELLING A SWEEP ITERATES OVER, and the order of traversal.
  ! Space and time are chosen separately: a dimension left coupled places
  ! every one of its members in the same partition member, and a dimension made
  ! sequential separates them. Time is traversed in the order its own
  ! discretisation couples the moments, so a member is solved after
  ! every member it reads.
  !===================================================================!

  subroutine sweep_labels(this, sequential_space, sequential_time, label, order)
    class(block_residual), intent(in)  :: this
    logical              , intent(in)  :: sequential_space, sequential_time
    integer, allocatable , intent(out) :: label(:), order(:)
    integer, allocatable :: slice(:), node(:), moment(:), in_time(:)
    integer :: nodes, k, t, m

    call this % labels_of(slice, node, moment)
    nodes = maxval(node)

    if (sequential_time) then
       call this % member_order(moment, order=in_time)
    else
       in_time = [1]
    end if

    if (sequential_time .and. sequential_space) then
       label = (moment - 1) * nodes + node
       allocate(order(size(in_time) * nodes))
       m = 0
       do t = 1, size(in_time)
          do k = 1, nodes
             m        = m + 1
             order(m) = (in_time(t) - 1) * nodes + k
          end do
       end do
    else if (sequential_time) then
       label = moment
       order = in_time
    else
       label = node
       order = [(k, k = 1, nodes)]
    end if

  end subroutine sweep_labels

  subroutine member_order(this, label, order)
    class(block_residual), intent(in)  :: this
    integer              , intent(in)  :: label(:)
    integer, allocatable , intent(out) :: order(:)
    integer, allocatable :: table(:,:)
    type(stencil) :: primary_law
    type(stored_directed_graph) :: coupling
    integer :: ne, e, n, t, h, members
    primary_law = this % primary_stencil()
    members = maxval(label)
    ne      = primary_law % pattern % num_edges()
    allocate(table(2, ne))
    n = 0
    do e = 1, ne
       t = label(primary_law % pattern % edge_tail(e))
       h = label(primary_law % pattern % edge_head(e))
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
  use operation_coupling      , only : weights_of, matrix_scheme_connectivity, connectivity_terms
  use gti_configuration       , only : refuse_unknown
  use operation_weight        , only : scheme_weight
  use view_directed_stored    , only : stored_directed_graph
  use view_directed_connectivity, only : connectivity_graph
  use view_directed           , only : directed_graph
  use field_calculus          , only : field
  use field_stored            , only : stored_field, typed_field_domain
  use operation_action      , only : variation, jacobian_of
  use operation_stencil       , only : stencil
  use operation_newton        , only : newton
  use operation_minimization  , only : minimizer, relative, by_rate
  use operation_temporal_minimization, only : temporal_minimizer
  use operation_dense_direct  , only : dense_direct
  use operation_gmres         , only : gmres
  use operation_family        , only : family
  use operation_grid          , only : grid, partition, partitioned
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use operation_expression       , only : expression
  use operation_domain        , only : continuous_domain, discrete_domain
  use gti_expansion           , only : family_container, expansion, marches_by_stages
  use gti_block               , only : block_residual
  use view_level              , only : level_member, level_num_members, level_coupling, &
       & level_couples
  use graph_fractal           , only : graph
  use map_value               , only : VALUE_KNOWN
  use gti_sweeps              , only : jacobian_present, multigrid_on, newton_order, &
       & set_aggregates, coarse_nodes, read_inner, store_inner, clear_inner, set_linear_stopping, &
       & stopping_applied
  use util_tally              , only : tally_record, tangent_loops, adjoint_loops
  implicit none
  type :: imbalance
     logical  :: converged = .true.
     logical  :: diverging = .false.
     real(dp) :: norm      = 0.0_dp
     real(dp) :: began     = 0.0_dp
     real(dp), allocatable :: by_degree(:)
     integer  :: largest_slot = 0, largest_degree = 0
     integer  :: steepest_slot = 0, steepest_degree = 0
     real(dp) :: steepest = 0.0_dp
  end type imbalance
  ! One step's edges, filled incrementally by stage_connectivity
  ! before the connectivity_graph they describe can be built - a
  ! graph is built whole, not edge by edge, so the raw lists are
  ! stored here until every edge of the step is placed.
  type :: raw_connectivity
     integer, allocatable :: tails(:), heads(:), tail_degree(:), head_degree(:)
   contains
     procedure :: place => raw_connectivity_place
  end type raw_connectivity
  type :: embedded_edge
     integer :: tail = 0, head = 0
     integer :: tail_degree = 0, head_degree = 0
     integer :: column_base = 0, row_base = 0
  end type embedded_edge
  type :: block_embedding
     type(expansion), pointer :: tower => null()
     type(graph)    , pointer :: block => null()
     integer :: num_slices = 0, num_stages = 0, num_degrees = 0
     integer :: moment_width = 0
     integer, allocatable :: slice_of(:), member_of(:)
   contains
     procedure :: block_connectivity => block_embedding_connectivity
     procedure :: stage_connectivity => stage_embedding_connectivity
  end type block_embedding
  ! Space and time are configured independently. A coupled dimension
  ! places all of its members in one system; a sequential one solves them one after another.
  character(len=16), save :: space_coupling = 'coupled'
  character(len=16), save :: time_coupling  = 'sequential'
  integer, save :: versions_given = 0
  real(dp), save :: stopping_tolerance  = 1.0e-12_dp
  integer , save :: stopping_criterion  = relative
  integer , save :: stopping_budget     = by_rate
  integer , save :: stopping_iterations = 100
  private
  public :: solved, unknowns_graph
  public :: block_from, instants_at_of
  public :: unknown, consistent_states, frozen_inputs
  public :: set_stopping
  public :: consistent_state
  public :: imbalance
  public :: swept, set_space_coupling, set_time_coupling, coupling_named
  public :: solve_linear, by_tangent, by_adjoint, next_version
  public :: weight_of, precision_needed
  public :: horizon_bounds
contains
  real(dp) function weight_of(scheme, degrees, step) result(w)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees
    real(dp)     , intent(in) :: step
    integer, allocatable :: offset(:), tail_degree(:)
    type(connectivity_graph) :: edges
    real(dp), allocatable :: c(:)
    integer :: d, depth, s, i, k
    logical :: any_pattern
    w = 1.0_dp
    any_pattern = .false.
    do d = 0, degrees - 1
       call scheme % row_pattern(d, degrees - 1, offset, tail_degree)
       if (size(offset) == 0) cycle
       any_pattern = .true.
       depth = maxval(offset)
       edges = connectivity_graph(depth + 1, &
            & [(depth + 1 - offset(k), k = 1, size(offset))], [(depth + 1, k = 1, size(offset))], &
            & tail_degree, [(d, k = 1, size(offset))])
       call weights_of(scheme_weight(scheme), edges, [(step, k = 1, depth + 1)], c)
       w = max(w, 1.0_dp + sum(abs(c)))
    end do
    if (any_pattern) return
    s = scheme % num_stages()
    do d = 0, degrees - 2
       do i = 1, s
          edges = connectivity_graph(s + 2, [1, (1 + k, k = 1, i)], &
               & [(1 + i, k = 0, i)], [d, (d + 1, k = 1, i)], [(d, k = 0, i)])
          call weights_of(scheme_weight(scheme), edges, [(step, k = 1, s + 2)], c)
          w = max(w, 1.0_dp + sum(abs(c)))
       end do
       edges = connectivity_graph(s + 2, [1, (1 + k, k = 1, s)], &
            & [(s + 2, k = 0, s)], [d, (d + 1, k = 1, s)], [(d, k = 0, s)])
       call weights_of(scheme_weight(scheme), edges, [(step, k = 1, s + 2)], c)
       w = max(w, 1.0_dp + sum(abs(c)))
    end do
  end function weight_of
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
    type(stored_field) :: state, design_field, direction
    type(typed_field_domain) :: point_scalars, point_states
    type(continuous_domain) :: continuous
    type(discrete_domain) :: point_domain
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
    continuous    = continuous_domain(physics)
    point_domain  = continuous % discrete(points)
    point_scalars = point_domain % design_fields()
    point_states  = point_domain % state_fields()
    design_field  = point_scalars % design(spread(design_value, 1, nodes))
    direction     = point_states % direction(e)
    began = -1.0_dp
    do iteration = 1, stopping_iterations
       state = point_states % state(q)
       call physics % apply(points, physics % bind([state, design_field]), out)
       call out % real_vector(r)
       r = r + below
       if (began < 0.0_dp) began = norm2(r)
       if (stopping_criterion == relative) then
          target = stopping_tolerance * max(began, tiny(1.0_dp))
       else
          target = stopping_tolerance
       end if
       if (norm2(r) <= target) return
       call physics % partial_action(points, physics % bind([state, design_field]), &
            & [variation(physics % argument(1), direction)], out)
       call out % real_vector(slope)
       do i = 1, nodes
          q(i * degrees) = q(i * degrees) - r(i) / slope(i)
       end do
    end do
    write(*,'(a,es12.3)') ' the residual of the physics at the initial instant is ', norm2(r)
    error stop 'gti_march: the initial state is consistent with the physics'
  end function consistent_states
  subroutine frozen_inputs(q, design, num_points, unknowns, inputs)
    real(dp), intent(in) :: q(:), design
    integer , intent(in) :: num_points
    type(stored_directed_graph)    , intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)
    type(typed_field_domain) :: states, designs
    unknowns = stored_directed_graph(size(q), tails=[integer ::], heads=[integer ::])
    allocate(inputs(2))
    states    = typed_field_domain(unknowns % vertex_set(), size(q))
    designs   = typed_field_domain(unknowns % vertex_set(), num_points)
    inputs(1) = states % state(q)
    inputs(2) = designs % design(spread(design, 1, num_points))
  end subroutine frozen_inputs
  subroutine set_stopping(tolerance, criterion, limit_kind, iterations)
    real(dp), intent(in) :: tolerance
    integer , intent(in) :: criterion, limit_kind, iterations
    call set_linear_stopping(tolerance, criterion, limit_kind)
    if (iterations < 1) then
       error stop 'gti_march: an iteration limit is positive'
    end if
    stopping_tolerance  = tolerance
    stopping_criterion  = criterion
    stopping_budget     = limit_kind
    stopping_iterations = iterations
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
  !===================================================================!
  ! THE OFFSET OF EACH INSTANT OF A BLOCK in the block's own state.
  ! The offset is determined by the tower's structure alone - the
  ! number of moments in each slice, and the width of a moment - so
  ! it is computed before the block is solved, and before the block
  ! is built. An assembler that must state which values one block
  ! passes to another calls this, and requires no prior evaluation.
  !===================================================================!

  function instants_at_of(tower, b, scheme, physics) result(instants_at)

    type(expansion) , intent(in), target :: tower
    integer         , intent(in) :: b
    class(family)   , intent(in) :: scheme
    type(expression), intent(in) :: physics
    integer, allocatable :: instants_at(:)

    type(graph), pointer :: horizon, block
    type(continuous_domain) :: continuous
    integer :: n, k, g, m, width, nd, stride
    logical :: staged

    continuous = continuous_domain(physics)
    nd     = continuous % equation_degree() + 1
    stride = continuous % num_components()
    horizon => level_member(level_member(tower % node(tower % root()), 1), 1)
    block   => level_member(horizon, b)
    n       = level_num_members(block)
    staged  = marches_by_stages(scheme, nd)
    m       = tower % extent_of(level_member(first_moment(block, staged), 1))
    width   = stride * m

    allocate(instants_at(n))
    g = 0
    do k = 1, n
       g = g + merge(level_num_members(level_member(block, k)), 1, staged)
       instants_at(k) = (g - 1) * width
    end do

  end function instants_at_of

  subroutine block_from(tower, b, scheme, physics, fixed, rows, instants_at)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)  :: b
    class(family)         , intent(in)  :: scheme
    type(expression)      , intent(in)  :: physics
    real(dp)              , intent(in)  :: fixed(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)
    type(graph), pointer :: horizon, block, slice, moment_node, component, below
    type(matrix_scheme_connectivity), allocatable :: connectivity(:)
    type(block_embedding) :: embedding
    type(continuous_domain) :: continuous
    integer , allocatable :: slice_of(:), member_of(:), members(:), at(:), fixed_rows(:)
    integer , allocatable :: r(:), c(:), table(:,:)
    real(dp), allocatable :: dt(:), w(:,:), seeds(:,:), spatial_weights(:)
    logical , allocatable :: point(:)
    integer :: m, nd, stride, width, n, s, k, j, g, moments, i, d, count, npts, ncar
    logical :: staged
    ! nd is the marching coordinate's degree count, which the scheme
    ! reads; stride is the point's whole component count, which the
    ! layout reads. The two differ once a rule names a second
    ! coordinate, so they are kept as distinct names.
    continuous = continuous_domain(physics)
    nd     = continuous % equation_degree() + 1
    stride = continuous % num_components()
    horizon => level_member(level_member(tower % node(tower % root()), 1), 1)
    block   => level_member(horizon, b)
    n       = level_num_members(block)
    call tower % value_of(block, dt)
    staged  = marches_by_stages(scheme, nd)
    s       = scheme % num_stages()
    m     = tower % extent_of(level_member(first_moment(block, staged), 1))
    width = stride * m
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
    allocate(fixed_rows(count), at(moments * m))
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
          do d = 0, stride - 1
             component => level_member(moment_node, d + 1)
             if (tower % status_of(component) == VALUE_KNOWN) then
                ncar = ncar + 1
                fixed_rows(ncar) = (g - 1) * width + (i - 1) * stride + d + 1
             end if
          end do
       end do
       if (point(g)) then
          do i = 1, m
             npts = npts + 1
             at(npts) = (g - 1) * width + (i - 1) * stride
          end do
          component => level_member(moment_node, scheme % primary_degree(nd - 1) + 1)
          if (level_couples(component) .and. .not. associated(below)) then
             below => level_coupling(component)
          end if
       end if
    end do
    if (size(fixed) /= ncar) then
       error stop 'gti_march: one fixed value per known component'
    end if
    embedding % tower => tower
    embedding % block => block
    embedding % num_slices = n
    embedding % num_stages = s
    embedding % num_degrees = nd
    embedding % moment_width = width
    embedding % slice_of = slice_of
    embedding % member_of = member_of
    if (staged) then
       call embedding % stage_connectivity(connectivity)
    else
       call embedding % block_connectivity(connectivity)
    end if
    allocate(seeds(size(dt), 0))
    call connectivity_terms(scheme, connectivity, m, nd, dt, seeds, r, c, w)
    rows = block_residual(derived_constraints(r, c, -w(:, 0), moments * width, 'time discretization stencil'), &
         & physics, at(1:npts), moments * width, nd, scheme % primary_degree(nd - 1), &
         & fixed_rows(1:ncar), fixed)
    call rows % placed_on(tower, block)
    call rows % with_connectivity(connectivity)
    if (associated(below)) then
       call tower % tuples_of(below, table)
       call tower % value_of(below, spatial_weights)
       call rows % spatial_discretization_laid(stencil(table(2, :), table(1, :), spatial_weights, &
            & spread(0.0_dp, 1, m), 'spatial discretization stencil'))
    end if
    instants_at = instants_at_of(tower, b, scheme, physics)
  end subroutine block_from
  function first_moment(block, staged) result(moment)
    type(graph), intent(in) :: block
    logical    , intent(in) :: staged
    type(graph), pointer :: moment
    moment => level_member(block, 1)
    if (staged) moment => level_member(moment, 1)
  end function first_moment
  subroutine block_embedding_connectivity(this, connectivity)
    class(block_embedding), intent(in) :: this
    type(matrix_scheme_connectivity), allocatable, intent(out) :: connectivity(:)
    integer, allocatable :: table(:,:)
    integer, allocatable :: tails(:), heads(:), tail_degree(:), head_degree(:)
    integer :: e, ne
    call this % tower % tuples_of(level_coupling(this % block), table)
    ne = size(table, 2)
    allocate(connectivity(1))
    connectivity(1) % step_of  = [(e, e = 1, this % num_slices)]
    allocate(tails(ne), heads(ne), tail_degree(ne), head_degree(ne))
    allocate(connectivity(1) % row(ne), connectivity(1) % column(ne))
    do e = 1, ne
       tails(e)       = (table(1, e) - 1) / this % num_degrees + 1
       tail_degree(e) = mod(table(1, e) - 1, this % num_degrees)
       heads(e)       = (table(2, e) - 1) / this % num_degrees + 1
       head_degree(e) = mod(table(2, e) - 1, this % num_degrees)
       connectivity(1) % column(e) = (tails(e) - 1) * this % moment_width + tail_degree(e) + 1
       connectivity(1) % row(e)    = (heads(e) - 1) * this % moment_width + head_degree(e) + 1
    end do
    connectivity(1) % graph = connectivity_graph(this % num_slices, tails, heads, tail_degree, head_degree)
  end subroutine block_embedding_connectivity
  subroutine stage_embedding_connectivity(this, connectivity)
    class(block_embedding), intent(in) :: this
    type(matrix_scheme_connectivity), allocatable, intent(out) :: connectivity(:)
    type(raw_connectivity), allocatable :: raw(:)
    integer, allocatable :: table(:,:), accumulate_state(:,:), first_moment(:), counted(:), filled(:)
    type(embedded_edge) :: edge
    integer :: kk, e, tail_moment, head_moment, vertex_tail, vertex_head
    allocate(connectivity(this % num_slices - 1), raw(this % num_slices - 1), &
         & first_moment(this % num_slices), counted(this % num_slices), filled(this % num_slices))
    first_moment(1) = 1
    do kk = 2, this % num_slices
       first_moment(kk) = first_moment(kk - 1) + merge(1, this % num_stages + 1, kk - 1 == 1)
    end do
    counted = 0
    do kk = 2, this % num_slices
       call this % tower % tuples_of(level_coupling(level_member(this % block, kk)), table)
       counted(kk) = size(table, 2)
    end do
    call this % tower % tuples_of(level_coupling(this % block), accumulate_state)
    do e = 1, size(accumulate_state, 2)
       head_moment = (accumulate_state(2, e) - 1) / this % num_degrees + 1
       kk          = this % slice_of(head_moment)
       counted(kk) = counted(kk) + 1
    end do
    do kk = 2, this % num_slices
       connectivity(kk - 1) % step_of  = spread(kk, 1, this % num_stages + 2)
       allocate(raw(kk - 1) % tails(counted(kk)), raw(kk - 1) % heads(counted(kk)), &
            &   raw(kk - 1) % tail_degree(counted(kk)), raw(kk - 1) % head_degree(counted(kk)), &
            &   connectivity(kk - 1) % row(counted(kk)), connectivity(kk - 1) % column(counted(kk)))
    end do
    filled = 0
    do kk = 2, this % num_slices
       call this % tower % tuples_of(level_coupling(level_member(this % block, kk)), table)
       do e = 1, size(table, 2)
          filled(kk) = filled(kk) + 1
          vertex_tail = (table(1, e) - 1) / this % num_degrees + 2
          vertex_head = (table(2, e) - 1) / this % num_degrees + 2
          edge % tail        = vertex_tail
          edge % head        = vertex_head
          edge % tail_degree = mod(table(1, e) - 1, this % num_degrees)
          edge % head_degree = mod(table(2, e) - 1, this % num_degrees)
          edge % column_base = (first_moment(kk) + vertex_tail - 3) * this % moment_width
          edge % row_base    = (first_moment(kk) + vertex_head - 3) * this % moment_width
          call raw(kk - 1) % place(connectivity(kk - 1), filled(kk), edge)
       end do
    end do
    do e = 1, size(accumulate_state, 2)
       tail_moment = (accumulate_state(1, e) - 1) / this % num_degrees + 1
       head_moment = (accumulate_state(2, e) - 1) / this % num_degrees + 1
       kk          = this % slice_of(head_moment)
       filled(kk)  = filled(kk) + 1
       edge % tail        = 1
       edge % head        = this % member_of(head_moment) + 1
       edge % tail_degree = mod(accumulate_state(1, e) - 1, this % num_degrees)
       edge % head_degree = mod(accumulate_state(2, e) - 1, this % num_degrees)
       edge % column_base = (tail_moment - 1) * this % moment_width
       edge % row_base    = (head_moment - 1) * this % moment_width
       call raw(kk - 1) % place(connectivity(kk - 1), filled(kk), edge)
    end do
    if (any(filled /= counted)) then
       error stop 'gti_march: every edge of a step is placed once'
    end if
    do kk = 2, this % num_slices
       connectivity(kk - 1) % graph = connectivity_graph(this % num_stages + 2, raw(kk - 1) % tails, raw(kk - 1) % heads, &
            & raw(kk - 1) % tail_degree, raw(kk - 1) % head_degree)
    end do
  end subroutine stage_embedding_connectivity
  subroutine raw_connectivity_place(this, placed, e, edge)
    class(raw_connectivity), intent(inout) :: this
    type(matrix_scheme_connectivity), intent(inout) :: placed
    integer            , intent(in) :: e
    type(embedded_edge), intent(in) :: edge
    this % tails(e)       = edge % tail
    this % heads(e)       = edge % head
    this % tail_degree(e) = edge % tail_degree
    this % head_degree(e) = edge % head_degree
    placed % column(e)    = edge % column_base + edge % tail_degree + 1
    placed % row(e)       = edge % row_base    + edge % head_degree + 1
  end subroutine raw_connectivity_place
  subroutine solved(rows, design_value, q, achieved, final_imbalance, seed)
    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: final_imbalance
    real(dp)       , intent(in) , optional :: seed(:)
    type(newton) :: solver
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    integer :: count, width
    count = rows % num_unknowns()
    if (present(seed)) then
       q = seed
    else
       q = at_first_instant(rows, count)
    end if
    call frozen_inputs(q, design_value, rows % num_points(), unknowns, inputs)
    width = rows % num_degrees()
    if (multigrid_on()) call set_aggregates(rows % aggregates(coarse_nodes(rows % num_nodes())))
    call read_inner(solver % inner, count, width)
    call solver % state(rows, unknowns, unknowns % vertex_set(), count, &
         & stored_inputs = [inputs(2)])
    solver % explicit       = jacobian_present()
    solver % higher_order_jacobian_product = newton_order()
    call stopping_applied(solver, stopping_tolerance, stopping_criterion, stopping_budget, &
         & stopping_iterations)
    call solver % solve(spread(0.0_dp, 1, count), q, achieved)
    call store_inner(solver % inner)
    if (present(final_imbalance)) call imbalance_of(solver, achieved, rows, unknowns, q, inputs(2), &
         & final_imbalance)
  end subroutine solved
  !===================================================================!
  ! THE IMBALANCE A SOLVE ENDED AT, and where it lies when the solve
  ! did not converge.
  !===================================================================!
  subroutine imbalance_of(solver, achieved, rows, unknowns, q, design, final_imbalance)
    class(minimizer)           , intent(in)  :: solver
    real(dp)                   , intent(in)  :: achieved, q(:)
    type(block_residual)       , intent(in)  :: rows
    type(stored_directed_graph), intent(in)  :: unknowns
    type(stored_field)         , intent(in)  :: design
    type(imbalance)            , intent(out) :: final_imbalance
    final_imbalance % converged = solver % converged(achieved)
    final_imbalance % diverging = solver % diverging(achieved)
    final_imbalance % norm      = achieved
    final_imbalance % began     = solver % began()
    if (.not. final_imbalance % converged) call by_aspect(rows, unknowns, q, design, final_imbalance)
  end subroutine imbalance_of
  integer function next_version() result(version)
    versions_given = versions_given + 1
    version = versions_given
  end function next_version
  subroutine solve_linear(rows, unknowns, inputs, rhs, transposed, version, w)
    type(block_residual)       , intent(in)  :: rows
    class(directed_graph)      , intent(in)  :: unknowns
    type(stored_field)         , intent(in)  :: inputs(:)
    real(dp)                   , intent(in)  :: rhs(:)
    logical                    , intent(in)  :: transposed
    integer                    , intent(in)  :: version
    real(dp), allocatable      , intent(out) :: w(:)
    type(block_residual) :: lin
    real(dp) :: achieved
    if (transposed) then
       call tally_record(adjoint_loops)
    else
       call tally_record(tangent_loops)
    end if
    lin = rows % linear_block(unknowns, rows % bind(inputs), rhs, transposed, version)
    call swept(lin, 0.0_dp, w, achieved)
  end subroutine solve_linear
  real(dp) function by_tangent(rows, unknowns, inputs, g, design_rate, explicit, version) &
       & result(df)
    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: version
    real(dp), allocatable :: w(:)
    type(typed_field_domain) :: unknown_fields
    type(stored_field) :: gradient, tangent
    call solve_linear(rows, unknowns, inputs, -design_rate, .false., version, w)
    unknown_fields = typed_field_domain(unknowns % vertex_set(), size(g))
    gradient       = unknown_fields % real_field('functional gradient', g)
    tangent        = unknown_fields % tangent(w)
    df = explicit + gradient % inner_product(tangent)
  end function by_tangent
  real(dp) function by_adjoint(rows, unknowns, inputs, g, design_rate, explicit, version) &
       & result(df)
    type(block_residual) , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: g(:), design_rate(:), explicit
    integer              , intent(in) :: version
    real(dp), allocatable :: lambda(:)
    type(typed_field_domain) :: unknown_fields
    type(stored_field) :: costate, forcing
    call solve_linear(rows, unknowns, inputs, g, .true., version, lambda)
    unknown_fields = typed_field_domain(unknowns % vertex_set(), size(g))
    costate        = unknown_fields % costate(lambda)
    forcing        = unknown_fields % forcing(design_rate)
    df = explicit - costate % inner_product(forcing)
  end function by_adjoint
  subroutine set_space_coupling(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['coupled   ', 'sequential'], 'space')
    space_coupling = name
  end subroutine set_space_coupling
  subroutine set_time_coupling(name)
    character(len=*), intent(in) :: name
    call refuse_unknown(name, ['coupled   ', 'sequential'], 'time')
    time_coupling = name
  end subroutine set_time_coupling
  pure function coupling_named() result(name)
    character(len=:), allocatable :: name
    name = 'space ' // trim(space_coupling) // ', time ' // trim(time_coupling)
  end function coupling_named
  subroutine swept(rows, design_value, q, achieved, final_imbalance)
    type(block_residual), intent(in)  :: rows
    real(dp)            , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)            , intent(out) :: achieved
    type(imbalance), intent(out), optional :: final_imbalance
    type(temporal_minimizer) :: solver
    type(newton) :: newton_solver
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    integer , allocatable :: order(:), label(:)
    real(dp), allocatable :: rhs(:)
    integer :: count, width
    logical :: sequential_space, sequential_time
    sequential_space = trim(space_coupling) == 'sequential'
    sequential_time  = trim(time_coupling)  == 'sequential'
    if (.not. sequential_space .and. .not. sequential_time) then
       call solved(rows, design_value, q, achieved, final_imbalance)
       return
    end if
    count = rows % num_unknowns()
    call rows % sweep_labels(sequential_space, sequential_time, label, order)
    q = at_first_instant(rows, count)
    q(rows % fixed_unknowns()) = rows % fixed_values()
    call frozen_inputs(q, design_value, rows % num_points(), unknowns, inputs)
    width = rows % num_degrees()
    if (multigrid_on()) call set_aggregates(rows % aggregates(coarse_nodes(rows % num_nodes())))
    call read_inner(newton_solver % inner, count, width)
    call newton_solver % state(rows, unknowns, unknowns % vertex_set(), count, &
         & stored_inputs = [inputs(2)])
    newton_solver % explicit       = jacobian_present()
    newton_solver % higher_order_jacobian_product = newton_order()
    call stopping_applied(newton_solver, stopping_tolerance, stopping_criterion, stopping_budget, &
         & stopping_iterations)
    allocate(solver % inner, source=newton_solver)
    call solver % state(rows, unknowns, unknowns % vertex_set(), count, &
         & stored_inputs = [inputs(2)])
    call stopping_applied(solver, stopping_tolerance, stopping_criterion, stopping_budget, &
         & stopping_iterations)
    call solver % partition(label, order, seed_from_previous=sequential_time)
    allocate(rhs(count), source=0.0_dp)
    call solver % solve(rhs, q, achieved)
    call clear_inner()
    if (present(final_imbalance)) call imbalance_of(solver, achieved, rows, unknowns, q, &
         & inputs(2), final_imbalance)
  end subroutine swept
  subroutine by_aspect(rows, unknowns, q, design, final_imbalance)
    type(block_residual)       , intent(in)    :: rows
    type(stored_directed_graph), intent(in)    :: unknowns
    real(dp)                   , intent(in)    :: q(:)
    type(stored_field)         , intent(in)    :: design
    type(imbalance)            , intent(inout) :: final_imbalance
    type(stored_field) :: state
    type(typed_field_domain) :: states
    class(field), allocatable :: out
    real(dp), allocatable :: r(:), a(:,:), slope(:)
    integer :: i, d, nd, n
    n  = size(q)
    nd = rows % num_degrees()
    states = typed_field_domain(unknowns % vertex_set(), n)
    state  = states % state(q)
    call rows % apply(unknowns, rows % bind([state, design]), out)
    call out % real_vector(r)
    allocate(final_imbalance % by_degree(0:nd - 1), source=0.0_dp)
    do i = 1, n
       d = mod(i - 1, nd)
       final_imbalance % by_degree(d) = final_imbalance % by_degree(d) + r(i) ** 2
    end do
    final_imbalance % by_degree = sqrt(final_imbalance % by_degree)
    i = maxloc(abs(r), dim=1)
    final_imbalance % largest_slot   = (i - 1) / nd + 1
    final_imbalance % largest_degree = mod(i - 1, nd)
    if (final_imbalance % norm <= 0.0_dp) return
    call jacobian_of(rows, unknowns, [state, design], n, unknowns % vertex_set(), a)
    slope = matmul(r, a) / final_imbalance % norm
    i = maxloc(abs(slope), dim=1)
    final_imbalance % steepest_slot   = (i - 1) / nd + 1
    final_imbalance % steepest_degree = mod(i - 1, nd)
    final_imbalance % steepest        = slope(i)
  end subroutine by_aspect
  function at_first_instant(rows, count) result(q)
    type(block_residual), intent(in) :: rows
    integer             , intent(in) :: count
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: one(:)
    integer , allocatable :: at(:)
    integer :: p, nd
    one = rows % first_fixed()
    nd  = size(one)
    at  = rows % points_at()
    allocate(q(count), source=0.0_dp)
    do p = 1, size(at)
       q(at(p) + 1:at(p) + nd) = one
    end do
  end function at_first_instant
  subroutine horizon_bounds(schemes, added, equation_degree, first, last)
    type(family_container), intent(in) :: schemes(:)
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
       ! A single-instant block is accepted where the family reaches
       ! back one instant: a first-order-history or multistage family.
       ! The horizon then advances one step per block at the family's
       ! own order. A deeper family still adds more instants than it
       ! reaches back over, so its given instants are in one block.
       if (added(b) < 1) then
          error stop 'gti_march: a block adds an instant at least'
       end if
       if ((b == 1 .or. schemes(b) % scheme % history_depth(equation_degree) > 1) .and. &
            & added(b) <= schemes(b) % scheme % history_depth(equation_degree)) then
          error stop 'gti_march: a block adds more instants than its family reaches'
       end if
       if (first(b) < 1) then
          error stop 'gti_march: the horizon contains every instant its blocks reach back over'
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
  use gti_expansion    , only : expansion, family_container
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
    type(family_container)  :: owner(1)
    type(block_residual) :: rows
    integer, allocatable :: at(:)
    real(dp), allocatable :: q(:)
    real(dp) :: achieved
    integer  :: last
    allocate(owner(1) % scheme, source=scheme)
    call tower % build(physics, owner, [n], uniform_grid(h), 0, design)
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
    real(dp), parameter :: safety = 0.9_dp, growth = 5.0_dp, shrinkage = 0.2_dp
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
          factor = min(growth, max(shrinkage, factor))
          if (e <= tolerance .or. h <= duration * 1.0e-10_dp) exit
          rejected = rejected + 1
          h = h * factor
          if (attempt > 50) then
             error stop 'gti_adaptive: a step remains above the tolerance after fifty attempts'
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
  public :: spatial_domain, spatial_mesh, spatial_operator, written_paraview, coarse_cells
  public :: cartesian, circular, elliptical, geometry_of
  integer, parameter :: cartesian  = 1
  integer, parameter :: circular   = 2
  integer, parameter :: elliptical = 3
  type :: spatial_domain
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
  end type spatial_domain
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
    type(spatial_domain) :: this
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
    type(spatial_domain), intent(inout) :: this
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
    type(spatial_domain)       , intent(inout) :: this
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
       tags(f)              = merge('    ', 'edge', faces(f) % head > 0)
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
    type(spatial_domain), intent(in) :: this
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op
    type(robin_condition) :: boundary_condition(1)
    if (degree < 1) then
       error stop 'gti_space: a form of degree below one fits no gradient'
    end if
    boundary_condition(1) = neumann('edge', 0.0_dp)
    op = diffusion_stencil(this % m, conduction(kappa), boundary_condition, polynomial_form(degree, this % m % dimension))
  end function spatial_operator
  function coarse_cells(this) result(aggregate)
    type(spatial_domain), intent(in) :: this
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
    type(spatial_domain)      , intent(in) :: this
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
  use gti_configuration, only : words_of
  use gti_march        , only : consistent_states
  use gti_space        , only : spatial_domain, spatial_operator, cartesian, written_paraview
  implicit none
  private
  public :: spatial_discretization_stencil_of, initial_field
  public :: against_the_laplacian, against_the_mode, export_instant
contains
  function spatial_discretization_stencil_of(space, kappa, degree) result(op)
    type(spatial_domain), intent(in) :: space
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op
    type(stencil) :: balance
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: lw(:), w(:), fixed(:)
    integer :: e, m
    balance = spatial_operator(space, kappa, degree)
    m       = balance % pattern % num_edges()
    call balance % weights % real_vector(lw)
    call balance % constants % real_vector(fixed)
    allocate(r(m), c(m), w(m))
    do e = 1, m
       r(e) = balance % pattern % edge_head(e)
       c(e) = balance % pattern % edge_tail(e)
       w(e) = -lw(e) / space % volume(r(e))
    end do
    op = stencil(r, c, w, fixed, 'spatial discretization stencil')
  end function spatial_discretization_stencil_of
  function initial_field(physics, degrees, kind, initial_state, design, spatial_discretization_stencil, space, a, b) &
       & result(q)
    type(expression)      , intent(in)           :: physics
    integer               , intent(in)           :: degrees
    character(len=*)      , intent(in)           :: kind, initial_state
    real(dp)              , intent(in)           :: design
    type(stencil)         , intent(in), optional :: spatial_discretization_stencil
    type(spatial_domain)            , intent(in), optional :: space
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
       given = words_of(initial_state)
       if (size(given) > degrees - 1) then
          write(*,'(a,i0,a)') ' the initial state contains the ', degrees - 1, &
               & ' components below the highest; the physics determines the highest.'
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
       if (space % geometry /= cartesian) error stop 'gti_field: the mode is defined on the rectangle'
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
    type(spatial_domain), intent(in) :: space
    real(dp)  , intent(in) :: a, b
    real(dp), allocatable :: shape(:)
    real(dp) :: pi
    integer  :: i
    pi = acos(-1.0_dp)
    shape = [(cos(pi * space % centre(1, i) / a) * cos(pi * space % centre(2, i) / b), &
         &    i = 1, space % num_cells)]
  end function mode_shape
  subroutine balance_of(space, kappa, degree, values, balanced)
    type(spatial_domain), intent(in) :: space
    real(dp)  , intent(in) :: kappa, values(:)
    integer   , intent(in) :: degree
    real(dp), allocatable, intent(out) :: balanced(:)
    type(stencil) :: op
    type(stored_field) :: given
    class(field), allocatable :: out
    op    = spatial_operator(space, kappa, degree)
    given = stored_field('values', op % pattern % vertex_set(), size(values))
    call given % set_real_vector(values)
    call op % apply(op % pattern, op % bind([given]), out)
    call out % real_vector(balanced)
  end subroutine balance_of
  subroutine against_the_laplacian(space, a, b, kappa, degree)
    type(spatial_domain), intent(in) :: space
    real(dp)  , intent(in) :: a, b, kappa
    integer   , intent(in) :: degree
    real(dp), allocatable :: shape(:), balanced(:), exact(:)
    real(dp) :: pi, err(0:2), norm(0:2)
    integer  :: i, boundary_count, count(0:2)
    if (space % geometry /= cartesian) then
       error stop 'gti_field: the laplacian check is defined on the rectangle'
    end if
    pi    = acos(-1.0_dp)
    shape = mode_shape(space, a, b)
    exact = -kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2) * shape
    call balance_of(space, kappa, degree, shape, balanced)
    err   = 0.0_dp
    norm  = 0.0_dp
    count = 0
    do i = 1, space % num_cells
       boundary_count = 0
       if (space % cell_ij(1, i) == 1 .or. space % cell_ij(1, i) == space % n2) boundary_count = boundary_count + 1
       if (space % cell_ij(2, i) == 1 .or. space % cell_ij(2, i) == space % n1) boundary_count = boundary_count + 1
       err(boundary_count)   = err(boundary_count)   + (balanced(i) / space % volume(i) - exact(i)) ** 2
       norm(boundary_count)  = norm(boundary_count)  + exact(i) ** 2
       count(boundary_count) = count(boundary_count) + 1
    end do
    write(*,'(a,i0,a,i0,a)') '   the operator compared with kappa times the laplacian of the mode, ', &
         & space % num_cells, ' cells, form degree ', degree, ':'
    write(*,'(a,3(a,es10.3))') '   relative rms error', &
         & '   interior ', sqrt(err(0) / max(norm(0), tiny(1.0_dp))), &
         & '   one boundary ', sqrt(err(1) / max(norm(1), tiny(1.0_dp))), &
         & '   corner ',   sqrt(err(2) / max(norm(2), tiny(1.0_dp)))
  end subroutine against_the_laplacian
  subroutine against_the_mode(space, a, b, kappa, degree, design, t_last, x, degrees)
    type(spatial_domain), intent(in) :: space
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
         & '      error at the last instant, compared with the mode ', sqrt(e_exact / area), &
         & '   semi-discrete ', sqrt(e_semi / area), '   omega ', omega, '   omega_h ', omega_h
  end subroutine against_the_mode
  !===================================================================!
  ! The name of a derivative of the state: the state itself is q,
  ! and each order appends the letter of the coordinate it is taken
  ! along, so a time derivative is named qt, qtt and a spatial one qx,
  ! qxx. The order is read from the name rather than counted from it.
  !===================================================================!
  pure function derivative_named(along, order) result(name)
    character(len=*), intent(in) :: along
    integer         , intent(in) :: order
    character(len=:), allocatable :: name
    name = 'q'
    if (order > 0) name = name // repeat(along, order)
  end function derivative_named
  subroutine export_instant(space, path, degrees, x)
    type(spatial_domain)      , intent(in) :: space
    character(len=*), intent(in) :: path
    integer         , intent(in) :: degrees
    real(dp)        , intent(in) :: x(:)
    character(len=8), allocatable :: names(:)
    real(dp), allocatable :: values(:,:)
    integer :: i, d
    allocate(names(degrees), values(space % num_cells, degrees))
    do d = 0, degrees - 1
       names(d + 1) = derivative_named('t', d)
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
  use operation_grid   , only : grid, partitioned, designed_grid
  use operation_expression, only : expression
  use gti_expansion    , only : family_container, marches_by_stages, expansion, &
       & design_of_physics, design_of_steps
  use gti_block        , only : block_residual
  use gti_march        , only : imbalance, swept, solve_linear, next_version, horizon_bounds, &
       & frozen_inputs
  use gti_march        , only : block_from, consistent_state
  use gti_sweeps       , only : choose
  use util_derivative_terms, only : derivative_terms, coefficient, mixed_partial, leibniz_parts, &
       & inner_product, operator(+), operator(-), operator(*)
  use operation_stencil, only : stencil
  use operation_family , only : crouzeix_three_stage
  use view_directed_stored, only : stored_directed_graph
  use field_calculus   , only : field, FIELD_REAL
  use operation_action , only : operation, emit, contract
  use operation_action , only : binding, is_bound, bound_value
  use operation_driver , only : driver, rule_graph, data_graph, pairing
  use operation_temporal_minimization, only : temporal_minimizer
  use view_read_write  , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use view_directed    , only : forward
  use view_directed    , only : directed_graph
  use field_stored     , only : stored_field, typed_field_domain
  use gti_sweeps       , only : pass_of, pass_substitutions, forward_pass, reverse_pass
  use util_tally            , only : tally_order, tally_enter, tally_leave
  use gti_configuration     , only : at_horizon, at_block, at_stage
  implicit none
  private
  public :: chain_block, march_chain, chain_expansion, instant_components
  public :: first_of, chain_derivative, asymmetry
  public :: multiset_count, multiset_rank, multiset_of, num_designs_of
  public :: chain_versions
  public :: expansion_substitutions
  public :: sink_costates
  public :: goal_oriented_partition
  public :: chain_incidence
  !===================================================================!
  ! THE LAYOUT OF ONE BLOCK'S STATE: the instants it covers, as the
  ! chain numbers them, from first to last by stride; the given
  ! instants its scheme reaches back over; the width of one instant,
  ! its degrees at every node; and the offset at which each covered
  ! instant begins inside the values.
  !
  !     first        first+stride      first+2*stride     instants
  !     |            |                 |
  !     [ ...... ] [ ...... ] [ ...... ]                 values
  !     ^          ^          ^
  !     instants_at(1)        instants_at(3)
  !===================================================================!
  type :: block_layout
     integer :: first = 0, last = 0, stride = 1, given = 0, width = 0, nodes = 1
     integer, allocatable :: instants_at(:)
  end type block_layout
  type, extends(block_layout) :: chain_block
     type(block_residual)  :: rows
     real(dp), allocatable :: state(:)
     integer               :: primary = 0
     logical :: staged = .false.
     ! the block storing each transferred unknown, zero for none, and
     ! the unknown's position in that block's state
     integer , allocatable :: source_block(:), source_at(:)
     integer , allocatable :: coarse_step(:)
     class(family), allocatable :: scheme
     real(dp)     , allocatable :: dt(:)
     real(dp)              :: fraction = 1.0_dp
     logical               :: counted = .true.
     real(dp) :: began = 0.0_dp
  end type chain_block
  type :: sink_costates
     integer , allocatable :: fixed_rows(:), last(:), interior(:)
     real(dp) :: departure = 0.0_dp
     real(dp) :: gradient  = 0.0_dp
     real(dp) :: costate   = 0.0_dp
     real(dp) :: unread    = 0.0_dp
  end type sink_costates
  !===================================================================!
  ! THE TANGENT TOWER OF ONE BLOCK: the derivatives of its state along
  ! every multiset of designs, of every order, w(unknown, multiset,
  ! order). Stored per block, so that a block's tower can be deallocated
  ! once the last block reading it has been solved.
  !===================================================================!
  type :: tangent_tower
     real(dp), allocatable :: w(:,:,:)
  end type tangent_tower
  !===================================================================!
  ! WHAT THE PIPELINED MARCH SHARES ACROSS ITS BLOCKS: the physics and
  ! the functionals, the derivative order, the towers with their last
  ! readers, the accumulated tables, and the storage counters. One
  ! design, the physics' parameter; a designed grid is rejected.
  !===================================================================!
  type :: taylor_context
     integer  :: order = 0, nf = 0, degrees = 0
     real(dp) :: design = 0.0_dp
     type(expression)              :: physics
     type(expression), allocatable :: functionals(:)
     real(dp)        , allocatable :: u(:,:,:)
     type(tangent_tower), allocatable :: w(:)
     integer , allocatable :: last_reader(:), versions(:)
     real(dp), allocatable :: table(:,:), by_order(:,:,:)
     integer :: tower_live = 0, tower_high = 0, tower_total = 0
     integer :: state_live = 0, state_high = 0, state_total = 0
  end type taylor_context
  !===================================================================!
  ! ONE BLOCK OF THE CHAIN, AS A RULE. A vertex of the operation
  ! graph stores an operation, and solving a block was a subroutine
  ! with twenty arguments, so here that subroutine becomes a rule.
  !
  ! EVERY COMPONENT OF THIS TYPE IS A DEPENDENCY, MADE EXPLICIT. A
  ! rule without data would store the scheme and the physics and
  ! nothing else; this one stores the chain it writes into, the tower
  ! it reads and the steps it marches over, because at present those
  ! are stored inside the rule rather than beside it. The separation
  ! of data from rule passes through this type, and the type is
  ! written this way so that the separation has a defined start.
  !
  !      (block 1) --> (block 2) --> (block 3)      the arcs
  !         |             |             |
  !       [state]       [state]       [state]       the output of each
  !===================================================================!

  ! THE PIPELINED MARCH'S CONTEXT, stored by the module for the one
  ! march in progress: a rule is copied through the driver, so the rule
  ! stores a flag and never a pointer.
  type(taylor_context), allocatable :: pipelined

  ! the two parts of the bipartite digraph, named for this caller
  integer, parameter :: BLOCKS = FIRST_PART
  integer, parameter :: STATES = SECOND_PART

  !===================================================================!
  ! ONE BLOCK'S STATE, AS A DATUM THAT STORES ITS OWN LAYOUT. The
  ! values a block solved for, and the layout needed to read an instant
  ! from them. A reader requests an instant, not an offset. So no
  ! offset is fixed before the march, and a block passes its state to
  ! the next reader without a layout agreed between the two beforehand.
  !===================================================================!

  type, extends(stored_field) :: block_state

     ! which block solved these values, so that where two blocks
     ! cover one instant the later of them is the one read
     integer :: at = 0

     type(block_layout) :: layout

   contains

     procedure :: covers                    ! whether an instant is covered
     procedure :: values_at                ! the values at one instant
     procedure :: assign_in => block_state_assign_in

  end type block_state

  type, extends(operation) :: block_rule

     type(chain_block), pointer :: chain(:) => null()
     type(expansion)  , pointer :: tower    => null()

     class(family)   , allocatable :: scheme
     type(expression), allocatable :: physics
     real(dp)        , allocatable :: dt(:), initial(:)
     integer         , allocatable :: coarse_step(:)

     ! the instants, stride and node count of the block; the given
     ! count and the offsets are the block's own once it is built
     type(block_layout) :: layout

     integer  :: at = 0, in_tower = 0, degrees = 0
     real(dp) :: fraction = 1.0_dp, design = 0.0_dp
     logical  :: counted = .true.

     ! the pipelined derivative, when one is requested of the march
     logical :: taylor = .false.

     ! the result the solve recorded, read after the graph is evaluated
     real(dp)        :: achieved = 0.0_dp
     type(imbalance) :: final_imbalance

   contains

     procedure :: name  => block_rule_name
     procedure :: apply => block_rule_apply

  end type block_rule

contains
  pure subroutine locate(chain, fine, owner_block, local)
    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: owner_block, local
    integer :: b
    owner_block = 0
    local   = 0
    do b = size(chain), 1, -1
       if (fine < chain(b) % first .or. fine > chain(b) % last) cycle
       if (mod(fine - chain(b) % first, chain(b) % stride) /= 0) cycle
       owner_block = b
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
    if (b == 0) error stop 'gti_chain: the instant lies outside the chain'
    at = chain(b) % instants_at(local)
    x  = chain(b) % state(at + 1:at + chain(b) % width)
  end function fine_components
  pure function instant_components(chain, instant) result(x)
    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: instant
    real(dp), allocatable :: x(:)
    x = fine_components(chain, 1 + (instant - 1) * chain(size(chain)) % stride)
  end function instant_components
  pure function transferred(earlier, first, stride, given) result(fixed)
    type(chain_block), intent(in) :: earlier(:)
    integer          , intent(in) :: first, stride, given
    real(dp), allocatable :: fixed(:)
    integer :: i, width
    width = earlier(1) % width
    allocate(fixed(given * width))
    do i = 1, given
       fixed((i - 1) * width + 1:i * width) = fine_components(earlier, first + (i - 1) * stride)
    end do
  end function transferred
  subroutine built(tower, b, scheme, physics, fixed, rows, instants_at)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)  :: b
    class(family)         , intent(in)  :: scheme
    type(expression)      , intent(in)  :: physics
    real(dp)              , intent(in)  :: fixed(:)
    type(block_residual)  , intent(out) :: rows
    integer, allocatable  , intent(out) :: instants_at(:)
    call block_from(tower, b, scheme, physics, fixed, rows, instants_at)
  end subroutine built
  subroutine march_chain(schemes, added, physics, degrees, steps, &
       & design, initial, chain, tower, dt, t, achieved, grid_design, final_imbalance, nodes, spatial_discretization_stencil, &
       & startup, functionals, derivative_order, f, tower_storage, state_storage)
    type(family_container)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    type(expression)      , intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    type(expansion), allocatable, intent(inout), target :: tower
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)
    type(imbalance), intent(out), optional :: final_imbalance
    integer        , intent(in) , optional :: nodes
    type(stencil)  , intent(in) , optional :: spatial_discretization_stencil
    integer        , intent(in) , optional :: startup
    ! THE PIPELINED DERIVATIVE. With functionals and an order given,
    ! every block's tangent tower is solved immediately after the block,
    ! the tables accumulate block by block, and states and towers are
    ! deallocated once their last reader has been solved. f is returned
    ! as (0:order, functional); each storage pair is [maximum live,
    ! total], towers and states separately.
    type(expression), intent(in), optional :: functionals(:)
    integer         , intent(in), optional :: derivative_order
    real(dp), allocatable, intent(out), optional :: f(:,:)
    integer, intent(out), optional :: tower_storage(2), state_storage(2)
    logical :: fused
    type(family_container), allocatable :: every(:)
    type(imbalance) :: one_imbalance
    integer , allocatable :: first(:), last(:), spans(:)
    real(dp), allocatable :: fine(:), design_field(:)
    real(dp) :: one_achieved
    integer :: b, k, r, given, before, m
    logical :: with_startup
    if (size(added) < 1) then
       error stop 'gti_chain: a chain contains at least one block'
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
    design_field = [real(dp) ::]
    allocate(every(size(added) + before))
    if (with_startup) then
       fine  = [0.0_dp, (dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (given - 1) * r)]
       design_field = fine(2:)
       allocate(every(1) % scheme, source=crouzeix_three_stage())
    end if
    do b = 1, size(added)
       if (b == 1 .and. with_startup) then
          design_field = [design_field, fine(size(fine)), dt(2:last(1))]
       else
          design_field = [design_field, dt(first(b) + merge(1, 0, b == 1):last(b))]
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
         & weights=grid_design, block_steps=design_field)
    fused = present(functionals) .and. present(derivative_order)
    if (fused) then
       if (allocated(pipelined)) deallocate(pipelined)
       allocate(pipelined)
       call taylor_prepare(pipelined, tower, functionals, derivative_order, degrees, &
            & schemes, added, r, first, last, before)
    end if
    achieved = 0.0_dp
    m = 1
    if (present(nodes)) m = nodes
    call tally_enter(at_horizon)
    if (with_startup) then
       call one_block(chain, 1, tower, 1, every(1) % scheme, physics, degrees, &
            & block_layout(1, (given - 1) * r + 1, 1, nodes=m), fine, &
            & [0, (1 + (k - 1) / r + 1, k = 1, (given - 1) * r)], &
            & 1.0_dp / real(r, dp), .false., design, initial, one_achieved, one_imbalance)
       achieved = one_achieved
       if (present(final_imbalance)) final_imbalance = one_imbalance
       if (fused) call taylor_block(pipelined, chain, 1)
    end if
    ! THE BLOCKS ARE A GRAPH, AND THE DRIVER EVALUATES IT. Block b
    ! reads block b - 1 in its one slot, so the arcs are the chain's
    ! own order and no loop here specifies which block is next.
    call marched_by_driver(chain, tower, schemes, added, physics, degrees, r, &
         & first, last, dt, design, initial, before, achieved, final_imbalance, nodes, taylor=fused)
    call tally_leave()
    if (fused) then
       pipelined % by_order(:, :, derivative_order) = pipelined % table
       if (present(f)) then
          allocate(f(0:derivative_order, size(functionals)))
          do k = 0, derivative_order
             f(k, :) = pipelined % by_order(:, 1, k)
          end do
       end if
       if (present(tower_storage)) tower_storage = [pipelined % tower_high, pipelined % tower_total]
       if (present(state_storage)) state_storage = [pipelined % state_high, pipelined % state_total]
       deallocate(pipelined)
    end if
  end subroutine march_chain
  !===================================================================!
  ! PREPARE THE PIPELINED MARCH: the shared context, and the last block
  ! reading each block's data, determined from the horizon's layout
  ! before any block is solved. A tower with more designs than the
  ! physics' parameter is rejected.
  !===================================================================!
  subroutine taylor_prepare(context, tower, functionals, order, degrees, &
       & schemes, added, r, first, last, before)
    type(taylor_context), intent(inout) :: context
    type(expansion)     , intent(in)    :: tower
    type(expression)    , intent(in)    :: functionals(:)
    integer             , intent(in)    :: order, degrees, r, before
    type(family_container) , intent(in)    :: schemes(:)
    integer             , intent(in)    :: added(:), first(:), last(:)
    real(dp), allocatable :: step_partials(:,:)
    integer , allocatable :: bfirst(:), blast(:), bstride(:), bgiven(:)
    integer :: nbb, bb, b, i, instant, owner, e
    if (order < 0) then
       error stop 'gti_chain: a derivative has an order of zero or more'
    end if
    if (num_designs_of(tower) /= 1) then
       error stop 'gti_chain: the pipelined march runs over the physics'' parameter alone'
    end if
    context % order   = order
    context % nf      = size(functionals)
    context % degrees = degrees
    call designs_of(tower, context % design, step_partials)
    context % physics     = tower % rule()
    context % functionals = functionals
    call steps_along(tower, 1, max(order, 1), context % u)
    nbb = size(added) + before
    allocate(context % w(nbb), context % versions(nbb), context % last_reader(nbb))
    allocate(context % table(context % nf, 1), context % by_order(context % nf, 1, 0:order), source=0.0_dp)
    allocate(bfirst(nbb), blast(nbb), bstride(nbb), bgiven(nbb))
    if (before == 1) then
       bfirst(1)  = 1
       blast(1)   = 1 + (schemes(1) % scheme % history_depth(degrees - 1) - 1) * r
       bstride(1) = 1
       bgiven(1)  = 0
    end if
    do b = 1, size(added)
       bfirst(before + b)  = 1 + (first(b) - 1) * r
       blast(before + b)   = 1 + (last(b) - 1) * r
       bstride(before + b) = r
       bgiven(before + b)  = schemes(b) % scheme % history_depth(degrees - 1)
    end do
    do bb = 1, nbb
       context % last_reader(bb) = bb
    end do
    do bb = 2, nbb
       do i = 1, bgiven(bb)
          instant = bfirst(bb) + (i - 1) * bstride(bb)
          owner  = 0
          do e = bb - 1, 1, -1
             if (instant < bfirst(e) .or. instant > blast(e)) cycle
             if (mod(instant - bfirst(e), bstride(e)) /= 0) cycle
             owner = e
             exit
          end do
          if (owner > 0) context % last_reader(owner) = max(context % last_reader(owner), bb)
       end do
    end do
  end subroutine taylor_prepare
  !===================================================================!
  ! THE PIPELINED TAYLOR STATE MARCH, ONE BLOCK'S CONTRIBUTION: solve
  ! the block's tangent tower with its own factorisation, add the
  ! block's contribution to every table, and deallocate each state and
  ! tower whose last reader has now been solved. The storage counters
  ! record the live size.
  !===================================================================!
  subroutine taylor_block(context, chain, at)
    type(taylor_context), intent(inout) :: context
    type(chain_block)   , intent(inout) :: chain(:)
    integer             , intent(in)    :: at
    integer :: h
    context % versions(at)   = next_version()
    context % state_live  = context % state_live + size(chain(at) % state)
    context % state_total = context % state_total + size(chain(at) % state)
    context % state_high  = max(context % state_high, context % state_live)
    call forward_block(chain, at, context % physics, context % functionals, context % degrees, &
         & context % design, 1, context % order, context % order, 0, context % order, context % versions, &
         & context % u, context % w, context % tower_live, context % tower_total, context % tower_high, &
         & context % table, context % by_order, last_reader=context % last_reader)
    call tally_order(0)
    do h = 1, at
       if (context % last_reader(h) /= at .or. .not. allocated(chain(h) % state)) cycle
       context % state_live = context % state_live - size(chain(h) % state)
       deallocate(chain(h) % state)
    end do
  end subroutine taylor_block
  !===================================================================!
  ! ONE BLOCK OF THE TAYLOR STATE MARCH. The block's tower: for every
  ! order to top and every multiset of the nd designs, one linear
  ! solve against the block's factorisation, the right side the rows'
  ! derivative along the multiset with the lower orders and the towers
  ! of the blocks storing the given instants substituted. Then the
  ! block's share of every functional at the multiset sizes from_size
  ! to to_size, the size order into table and every other into
  ! by_order. Then, when last_reader is given, every tower whose last
  ! reader is this block is released. live, total and peak_storage
  ! count the tower numbers stored.
  !===================================================================!
  subroutine forward_block(chain, b, physics, functionals, degrees, design, nd, top, order, &
       & from_size, to_size, versions, u, w, live, total, peak_storage, table, by_order, node_measure, &
       & last_reader)
    type(chain_block)  , intent(in)    :: chain(:)
    integer            , intent(in)    :: b, degrees, nd, top, order, from_size, to_size, versions(:)
    type(expression)   , intent(in)    :: physics, functionals(:)
    real(dp)           , intent(in)    :: design, u(:,:,:)
    type(tangent_tower), intent(inout) :: w(:)
    integer            , intent(inout) :: live, total, peak_storage
    real(dp), intent(inout), optional  :: table(:,:), by_order(:,:,0:)
    real(dp), intent(in)   , optional  :: node_measure(:)
    integer , intent(in)   , optional  :: last_reader(:)
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: r(:), one(:)
    integer , allocatable :: s(:)
    real(dp) :: share
    integer :: count, k, rank, i, h, size_of
    count = chain(b) % rows % num_unknowns()
    allocate(w(b) % w(count, multiset_count(nd, max(top, 1)), max(top, 1)), source=0.0_dp)
    live         = live + size(w(b) % w)
    total        = total + size(w(b) % w)
    peak_storage = max(peak_storage, live)
    do k = 1, top
       call tally_order(k)
       call tally_enter(at_horizon)
       do rank = 1, multiset_count(nd, k)
          s = multiset_of(rank, k, nd)
          call tally_enter(at_block)
          call forcing_of(chain, b, physics, degrees, design, s, w, u, nd, r)
          r = -r
          do i = 1, size(chain(b) % source_at)
             if (chain(b) % source_block(i) > 0) then
                r(i) = w(chain(b) % source_block(i)) % w(chain(b) % source_at(i), rank, k)
             end if
          end do
          call frozen_at(chain(b), design, unknowns, inputs)
          call solve_linear(chain(b) % rows, unknowns, inputs, r, .false., versions(b), one)
          w(b) % w(1:count, rank, k) = one
          call tally_leave()
       end do
       call tally_leave()
    end do
    do size_of = from_size, to_size
       do rank = 1, multiset_count(nd, size_of)
          s = multiset_of(rank, size_of, nd)
          do i = 1, size(functionals)
             share = functional_along(chain, b, functionals(i), degrees, design, s, 0, w, u, nd, &
                  & node_measure)
             if (size_of == order) then
                table(i, rank) = table(i, rank) + share
             else
                by_order(i, rank, size_of) = by_order(i, rank, size_of) + share
             end if
          end do
       end do
    end do
    if (.not. present(last_reader)) return
    do h = 1, b
       if (last_reader(h) == b .and. allocated(w(h) % w)) then
          live = live - size(w(h) % w)
          deallocate(w(h) % w)
       end if
    end do
  end subroutine forward_block
  !===================================================================!
  ! THE CHAIN AS A BIPARTITE DIGRAPH, AND NOTHING ELSE. Which block
  ! writes which state, which block reads which, and the transpose the
  ! derivative is taken over. The incidence computes nothing and
  ! marches nothing, so the evaluation order and every datum's lifetime
  ! can be read from it before any block is solved.
  !===================================================================!

  function chain_incidence(schemes, added, degrees, r, first, last) result(incidence)

    type(family_container), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:), degrees, r, first(:), last(:)
    type(bipartite_digraph) :: incidence

    integer, allocatable :: from_part(:), from_vertex(:), to_part(:), to_vertex(:)
    integer :: nb, b, a, n_arc, forward_arcs, given, i, e, instant, owner

    nb = size(added)

    ! THE ARCS ARE READS AND WRITES, and nothing here states which
    ! block follows which. Block b writes datum b; block b reads
    ! datum b - 1. The order over the blocks is the projection of
    ! those arcs onto the blocks, and the driver derives it.
    !
    !      (b1) --> [d1] --> (b2) --> [d2] --> (b3)
    !       write     read    write    read
    !
    ! at most one write per block, and at most one read per instant
    ! its scheme reaches back over
    n_arc = nb
    do b = 2, nb
       n_arc = n_arc + schemes(b) % scheme % history_depth(degrees - 1)
    end do
    ! capacity for the forward arcs, their transpose, and the state each
    ! block of the transpose reads
    allocate(from_part(2 * n_arc + nb), from_vertex(2 * n_arc + nb), &
         & to_part(2 * n_arc + nb), to_vertex(2 * n_arc + nb))
    a = 0
    do b = 1, nb
       a = a + 1                                  ! block b writes datum b
       from_part(a)   = BLOCKS
       from_vertex(a) = b
       to_part(a)     = STATES
       to_vertex(a)   = b
    end do
    do b = 2, nb
       ! block b is passed the instants its scheme reaches back over,
       ! and each of them is stored by the latest earlier block that
       ! covers it - the block locate returns. That block, and no
       ! other, is read.
       given = schemes(b) % scheme % history_depth(degrees - 1)
       do i = 1, given
          instant = 1 + (first(b) - 1) * r + (i - 1) * r
          owner  = 0
          do e = b - 1, 1, -1
             if (instant < 1 + (first(e) - 1) * r) cycle
             if (instant > 1 + (last(e)  - 1) * r) cycle
             if (mod(instant - (1 + (first(e) - 1) * r), r) /= 0) cycle
             owner = e
             exit
          end do
          if (owner == 0) cycle
          if (any(from_vertex(1:a) == owner .and. to_vertex(1:a) == b &
               & .and. from_part(1:a) == STATES)) cycle
          a = a + 1
          from_part(a)   = STATES
          from_vertex(a) = owner
          to_part(a)     = BLOCKS
          to_vertex(a)   = b
       end do
    end do
    forward_arcs = a

    ! THE TRANSPOSE, AND WHY IT IS HERE. The derivative is taken by a
    ! reverse pass over these same blocks, and that pass reads every
    ! block's state. Until those reads are arcs the graph does not
    ! contain them, and a lifetime read from it would deallocate a
    ! state the reverse pass still reads. So block b of the transpose
    ! is placed at nb + b, and three families of arc define it:
    !
    !   forward     (b) --> [d_b] --> (c)        c reads what b wrote
    !
    !   the state   [d_b] --> (nb+b)             the reverse pass at
    !               read again, so d_b is        b reads b's own state
    !               allocated until step nb+b
    !
    !   transposed  (nb+c) --> [e_c] --> (nb+b)  every forward arc,
    !               with its ends exchanged      reversed
    !
    ! The transposed arcs place nb+nb ... nb+1 after every forward
    ! block and in decreasing order, which is the order the reverse
    ! pass runs in. So last_reader_of(d_b) is step nb+b, and the
    ! states are deallocated one at a time as the reverse pass
    ! completes them rather than being retained to the end.
    !
    ! No rule is placed at a transposed vertex yet: chain_derivative
    ! still computes the costates. What is placed here is the
    ! lifetime, which is the quantity the driver requires.
    do e = 1, forward_arcs
       a = a + 1
       if (from_part(e) == BLOCKS) then
          from_part(a)   = BLOCKS
          from_vertex(a) = nb + from_vertex(e)
          to_part(a)     = STATES
          to_vertex(a)   = nb + to_vertex(e)
       else
          from_part(a)   = STATES
          from_vertex(a) = nb + to_vertex(e)
          to_part(a)     = BLOCKS
          to_vertex(a)   = nb + from_vertex(e)
       end if
    end do
    do b = 1, nb
       a = a + 1
       from_part(a)   = STATES
       from_vertex(a) = b
       to_part(a)     = BLOCKS
       to_vertex(a)   = nb + b
    end do
    n_arc = a
    incidence = bipartite_digraph(2 * nb, 2 * nb, from_part(1:n_arc), from_vertex(1:n_arc), &
         & to_part(1:n_arc), to_vertex(1:n_arc))

  end function chain_incidence

  !===================================================================!
  ! ASSEMBLE THE CHAIN AS A GRAPH AND PASS IT TO THE DRIVER. The
  ! incidence is a value built separately, and all that is added here
  ! is a rule at each forward vertex: the transposed vertices store
  ! none, because the reverse pass they represent is still
  ! chain_derivative's. What the driver does with any of it - the order,
  ! the lifetimes, and in future which of them run concurrently - is not stated here.
  !===================================================================!

  subroutine marched_by_driver(chain, tower, schemes, added, physics, degrees, r, &
       & first, last, dt, design, initial, before, achieved, final_imbalance, nodes, &
       & data_stored, taylor)

    type(chain_block)  , intent(inout), target :: chain(:)
    type(expansion)    , intent(in)   , target :: tower
    type(family_container), intent(in)   :: schemes(:)
    integer            , intent(in)   :: added(:), degrees, r, before
    type(expression)   , intent(in)   :: physics
    integer            , intent(in)   :: first(:), last(:)
    real(dp)           , intent(in)   :: dt(:), design, initial(:)
    real(dp)           , intent(inout):: achieved
    type(imbalance), intent(inout), optional :: final_imbalance
    integer        , intent(in)   , optional :: nodes
    ! WHERE THE STATES ARE STORED. The driver places every block's
    ! state at its data vertex, and a caller that requires them
    ! requests the data graph rather than the chain.
    type(data_graph), intent(out) , optional :: data_stored
    logical, intent(in), optional :: taylor

    type(data_graph)            :: values
    type(pairing)               :: pairs
    type(rule_graph)       :: rules
    type(bipartite_digraph)     :: incidence
    type(driver)                :: schedule
    type(temporal_minimizer)    :: executor
    type(block_rule)            :: one
    type(contract), allocatable :: contracts(:)
    type(stored_directed_graph) :: bare
    integer, allocatable :: reads(:)
    real(dp), allocatable :: no_rhs(:), no_solution(:)
    real(dp) :: driver_achieved
    integer :: nb, b, k

    nb = size(added)

    incidence = chain_incidence(schemes, added, degrees, r, first, last)

    ! one rule per forward vertex, storing what solving that block
    ! requires; the transposed vertices store none and compute nothing
    allocate(rules % at(2 * nb), values % at(2 * nb))
    do b = 1, nb
       one % chain    => chain
       one % tower    => tower
       one % at       = before + b
       one % in_tower = before + b
       allocate(one % scheme , source=schemes(b) % scheme)
       allocate(one % physics, source=physics)
       one % degrees     = degrees
       one % layout      = block_layout(1 + (first(b) - 1) * r, 1 + (last(b) - 1) * r, r)
       if (present(nodes)) one % layout % nodes = nodes
       one % dt          = dt(first(b):last(b))
       one % coarse_step = [(k, k = first(b), last(b))]
       one % fraction    = 1.0_dp
       one % counted     = .true.
       if (present(taylor)) one % taylor = taylor
       one % design      = design
       one % initial     = initial
       call incidence % in_neighbourhood(FIRST_PART, b, reads)
       allocate(contracts(size(reads)), source=contract(FIELD_REAL, 1))
       call one % declare_arguments(size(reads), contracts)
       allocate(rules % at(b) % rule, source=one)
       deallocate(contracts)
       deallocate(one % scheme, one % physics)
    end do

    bare = stored_directed_graph(nb, tails=[integer ::], heads=[integer ::])
    schedule = driver(rules % at(1) % rule, incidence, forward)
    call executor % state(schedule, bare, bare % vertex_set(), 0)
    call executor % pair_with(rules % pair(values))
    allocate(no_rhs(0), no_solution(0))
    call executor % solve(no_rhs, no_solution, driver_achieved)

    ! WHERE THE STATES ARE STORED. The driver placed each block's
    ! state at its data vertex; a caller that requires them requests
    ! the data branch rather than reading the chain.
    pairs = executor % pairing_of()
    if (present(data_stored)) data_stored = pairs % stored_data()


    ! the result each block recorded, read from the rules
    do b = 1, nb
       select type (solved_rule => rules % at(b) % rule)
       type is (block_rule)
          achieved = max(achieved, solved_rule % achieved)
          if (present(final_imbalance)) then
             if (before + b == 1) final_imbalance = solved_rule % final_imbalance
             if (final_imbalance % converged .and. .not. solved_rule % final_imbalance % converged) final_imbalance = solved_rule % final_imbalance
          end if
       end select
    end do

  end subroutine marched_by_driver

  !===================================================================!
  ! WHETHER THIS STATE COVERS AN INSTANT. The block's instants run
  ! from first to last by stride, and no others are stored.
  !===================================================================!

  pure logical function covers(this, instant)
    class(block_state), intent(in) :: this
    integer           , intent(in) :: instant
    covers = .false.
    if (instant < this % layout % first) return
    if (instant > this % layout % last)  return
    if (this % layout % stride < 1) return
    if (mod(instant - this % layout % first, this % layout % stride) /= 0) return
    covers = .true.
  end function covers

  !===================================================================!
  ! THE VALUES THIS STATE STORES AT ONE INSTANT. The instant's position
  ! among the covered ones gives its offset, and the caller's width
  ! specifies how many values begin there.
  !===================================================================!

  subroutine values_at(this, instant, width, values)
    class(block_state)   , intent(in)  :: this
    integer              , intent(in)  :: instant, width
    real(dp), allocatable, intent(out) :: values(:)
    real(dp), allocatable :: whole(:)
    integer :: local, offset
    if (.not. this % covers(instant)) then
       error stop 'gti_chain: a state covers the instant requested of it'
    end if
    local = (instant - this % layout % first) / this % layout % stride + 1
    if (.not. allocated(this % layout % instants_at)) then
       error stop 'gti_chain: a state stores the offsets of the instants it covers'
    end if
    if (local < 1 .or. local > size(this % layout % instants_at)) then
       error stop 'gti_chain: a state stores the offsets of the instants it covers'
    end if
    offset = this % layout % instants_at(local)
    call this % real_vector(whole)
    if (offset < 0 .or. offset + width > size(whole)) then
       error stop 'gti_chain: an instant lies inside the values a state stores'
    end if
    values = whole(offset + 1:offset + width)
  end subroutine values_at

  !===================================================================!
  ! Place this state at a location that is a block state. A location
  ! of any other type would lose the layout and is an error.
  !===================================================================!

  subroutine block_state_assign_in(this, location)
    class(block_state), intent(in)    :: this
    class(field)      , intent(inout) :: location
    select type (location)
    type is (block_state)
       location = this
    class default
       error stop 'gti_chain: a block state is placed at a block state'
    end select
  end subroutine block_state_assign_in

  pure function block_rule_name(this) result(name)
    class(block_rule), intent(in) :: this
    character(len=:), allocatable :: name
    character(len=12) :: digits
    write(digits,'(i0)') this % at
    name = 'block ' // trim(digits) // ' of the chain'
  end function block_rule_name

  !===================================================================!
  ! Solving the block this rule is placed at, and storing its state as
  ! the datum of the vertex. The inputs a driver gathers are not
  ! read yet: what this block requires from the one before it is still
  ! stored inside the chain the rule points at.
  !===================================================================!

  subroutine block_rule_apply(this, input_graph, inputs, output)
    class(block_rule)        , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)             , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output
    class(field), allocatable :: value
    type(block_state) :: state
    type(stored_directed_graph) :: state_domain
    type(typed_field_domain) :: states
    real(dp), allocatable :: transferred_values(:), one_datum(:)
    real(dp) :: achieved
    type(imbalance) :: final_imbalance
    integer :: i, e, given, width, instant, owner_block, taken
    if (.not. associated(this % chain) .or. .not. associated(this % tower)) then
       error stop 'gti_chain: a block rule is placed at a block of a chain'
    end if
    ! THE TRANSFER, READ FROM THE DATA THEMSELVES. Each instant the
    ! scheme reaches back over is requested of every datum passed in,
    ! and a block's state returns it from the layout stored inside it.
    ! Where two states cover one instant the later block is read, which
    ! is the chain's own order. Nothing is read from the chain, and no
    ! offset is fixed before the march: this is the separation of data
    ! from rule, and it is complete here.
    given = 0
    if (present(inputs)) then
       if (size(inputs) > 0) given = this % scheme % history_depth(this % degrees - 1)
    end if
    if (given > 0) then
       width = this % degrees * this % layout % nodes
       allocate(transferred_values(given * width))
       do i = 1, given
          instant = this % layout % first + (i - 1) * this % layout % stride
          owner_block = 0
          taken   = 0
          do e = 1, this % num_arguments()
             if (.not. is_bound(inputs, this % argument(e))) cycle
             call bound_value(inputs, this % argument(e), value)
             select type (datum => value)
             type is (block_state)
                if (.not. datum % covers(instant)) cycle
                if (datum % at < owner_block) cycle
                owner_block = datum % at
                taken   = e
             end select
          end do
          if (taken < 1) then
             error stop 'gti_chain: a block is passed the data its scheme reaches back over'
          end if
          call bound_value(inputs, this % argument(taken), value)
          select type (datum => value)
          type is (block_state)
             call datum % values_at(instant, width, one_datum)
          end select
          transferred_values((i - 1) * width + 1:i * width) = one_datum
       end do
    end if

    ! an unallocated transfer is an absent argument
    call one_block(this % chain, this % at, this % tower, this % in_tower, this % scheme, &
         & this % physics, this % degrees, this % layout, this % dt, this % coarse_step, &
         & this % fraction, this % counted, this % design, this % initial, achieved, &
         & final_imbalance, transferred_values)
    ! THE DATUM'S DOMAIN IS THE BLOCK'S, NOT THE SCHEDULE'S. The graph
    ! a driver evaluates over specifies which rule runs when; it
    ! specifies nothing about how many points a state stores, and the
    ! two counts are unrelated. Substituting one for the other would
    ! give a field a domain it does not have.
    state_domain = stored_directed_graph(this % chain(this % at) % rows % num_points(), &
         & tails=[integer ::], heads=[integer ::])
    states = typed_field_domain(state_domain % vertex_set(), &
         & size(this % chain(this % at) % state))
    state % stored_field = states % state(this % chain(this % at) % state)
    state % at     = this % at
    state % layout = this % chain(this % at) % block_layout
    call emit(state, output)
    ! THE PIPELINED DERIVATIVE. With a context attached, the block's
    ! tower is solved immediately after the block, and every datum
    ! with no later reader is deallocated.
    if (this % taylor .and. allocated(pipelined)) then
       call taylor_block(pipelined, this % chain, this % at)
    end if
  end subroutine block_rule_apply

  subroutine one_block(chain, b, tower, in_tower, scheme, physics, degrees, layout, &
       & dt, coarse_step, fraction, counted, design, initial, achieved, final_imbalance, &
       & transferred_values)
    type(chain_block)     , intent(inout) :: chain(:)
    type(expansion)       , intent(in), target :: tower
    integer               , intent(in)    :: b, in_tower, degrees
    type(block_layout)    , intent(in)    :: layout
    integer               , intent(in)    :: coarse_step(:)
    class(family)         , intent(in)    :: scheme
    type(expression)      , intent(in)    :: physics
    real(dp)              , intent(in)    :: dt(:), fraction, design, initial(:)
    logical               , intent(in)    :: counted
    real(dp)              , intent(out)   :: achieved
    type(imbalance)       , intent(out)   :: final_imbalance
    real(dp)     , intent(in), optional   :: transferred_values(:)
    real(dp), allocatable :: fixed(:)
    chain(b) % block_layout = layout
    chain(b) % given        = scheme % history_depth(degrees - 1)
    chain(b) % primary      = scheme % primary_degree(degrees - 1)
    chain(b) % width        = degrees * layout % nodes
    allocate(chain(b) % scheme, source=scheme)
    chain(b) % staged      = marches_by_stages(scheme, degrees)
    chain(b) % dt          = dt
    chain(b) % coarse_step = coarse_step
    chain(b) % fraction    = fraction
    chain(b) % counted     = counted
    call transfer_layout(chain, b)
    ! WHAT IS PASSED TO THIS BLOCK. The first block is passed the
    ! initial state; every other is passed the instants its scheme
    ! reaches back over, read from the earlier block that stores
    ! each. A caller that already has the transfer passes it, and then
    ! nothing here reads the chain - that argument is the separation
    ! between the data and the rule, and on this side it is incomplete.
    if (present(transferred_values)) then
       fixed = transferred_values
    else if (b == 1) then
       if (size(initial) /= chain(b) % given * chain(b) % width) then
          error stop 'gti_chain: the initial state contains the first block''s given instants'
       end if
       fixed = initial
    else
       fixed = transferred(chain(1:b - 1), layout % first, layout % stride, chain(b) % given)
    end if
    if (scheme % num_stages() > 1) then
       call tally_enter(at_stage)
    else
       call tally_enter(at_block)
    end if
    call built(tower, in_tower, scheme, physics, fixed, chain(b) % rows, &
         & chain(b) % instants_at)
    call swept(chain(b) % rows, design, chain(b) % state, achieved, final_imbalance)
    chain(b) % began = final_imbalance % began
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
    integer , allocatable :: versions(:)
    real(dp), allocatable :: by_order(:,:,:), table(:,:)
    integer :: m
    call chain_versions(chain, tower, functionals, degrees, versions)
    call chain_derivative(chain, tower, versions, functionals, degrees, max_order, forward_pass, &
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
  subroutine chain_versions(chain, tower, functionals, degrees, versions, node_measure)
    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    type(expression) , intent(in) :: functionals(:)
    integer          , intent(in) :: degrees
    integer, allocatable, intent(out) :: versions(:)
    real(dp), intent(in), optional :: node_measure(:)
    integer :: b
    associate (u1 => tower, u2 => functionals, u3 => degrees, u4 => node_measure); end associate
    allocate(versions(size(chain)))
    do b = 1, size(chain)
       versions(b) = next_version()
    end do
  end subroutine chain_versions
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
  !===================================================================!
  ! WHERE EACH TRANSFERRED UNKNOWN OF BLOCK b IS STORED: for the i-th
  ! component of its given instants, the latest earlier block covering
  ! that instant - zero for none - and the component's position in
  ! that block's state. Read from the earlier blocks once, when block
  ! b is built.
  !===================================================================!
  pure subroutine transfer_layout(chain, b)
    type(chain_block), intent(inout) :: chain(:)
    integer          , intent(in)    :: b
    integer :: i, n, local, owner_block
    n = chain(b) % given * chain(b) % width
    allocate(chain(b) % source_block(n), chain(b) % source_at(n), source=0)
    do i = 1, n
       call locate(chain(1:b - 1), chain(b) % first + ((i - 1) / chain(b) % width) * chain(b) % stride, &
            & owner_block, local)
       chain(b) % source_block(i) = owner_block
       if (owner_block > 0) then
          chain(b) % source_at(i) = chain(owner_block) % instants_at(local) + mod(i - 1, chain(b) % width) + 1
       end if
    end do
  end subroutine transfer_layout
  subroutine chain_derivative(chain, tower, versions, functionals, degrees, order, pass_kind, &
       & table, node_measure, entries, designs, by_order, sinks, leibniz, tower_storage)
    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    integer          , intent(in) :: versions(:)
    type(expression) , intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, order, pass_kind
    real(dp), allocatable  , intent(out) :: table(:,:)
    real(dp), intent(in), optional      :: node_measure(:)
    real(dp), allocatable, intent(out), optional :: entries(:,:,:)
    integer, intent(in), optional :: designs
    real(dp), allocatable, intent(out), optional :: by_order(:,:,:)
    type(sink_costates), intent(out), optional :: sinks
    real(dp), allocatable, intent(out), optional :: leibniz(:,:,:,:)
    ! the maximum live and the total tower sizes: the storage the pass used
    integer, intent(out), optional :: tower_storage(2)
    type(tangent_tower), allocatable :: w(:)
    real(dp), allocatable :: lambda(:,:,:,:,:), u(:,:,:)
    integer , allocatable :: last_reader(:)
    integer :: live, peak_storage, total
    type(derivative_terms) :: l
    type(derivative_terms), allocatable :: products(:,:,:)
    real(dp), allocatable :: split(:)
    integer :: m, rank_s, mask, from_size, to_size
    logical , allocatable :: is_sink(:,:)
    real(dp), allocatable :: diagonal(:,:)
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: step_partials(:,:), rhs(:,:), one(:), r(:), every(:,:,:)
    integer , allocatable :: s(:)
    type(expression) :: physics
    real(dp) :: design
    logical  :: forward
    integer :: nf, nd, nb, top, widest, k, b, i, j, p, count, rank, owner_block
    if (order < 0) then
       error stop 'gti_chain: a derivative has an order of zero or more'
    end if
    if (pass_kind /= forward_pass .and. pass_kind /= reverse_pass) then
       error stop 'gti_chain: a pass is forward or reverse'
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
    if (pass_kind == reverse_pass) top = max(order - 1, 0)
    forward = pass_kind == forward_pass .or. order == 0
    call steps_along(tower, nd, max(order, 1), u)
    if (present(by_order)) then
       allocate(by_order(nf, multiset_count(nd, order), 0:order), source=0.0_dp)
    end if
    allocate(w(nb))
    allocate(rhs(widest, nb))
    if (forward) then
       allocate(table(nf, multiset_count(nd, order)), source=0.0_dp)
       ! THE LAST BLOCK THAT READS EACH BLOCK'S TOWER: the latest block
       ! whose given instants that block stores, and the block itself.
       last_reader = [(b, b = 1, nb)]
       do b = 2, nb
          do i = 1, size(chain(b) % source_block)
             owner_block = chain(b) % source_block(i)
             if (owner_block > 0) last_reader(owner_block) = max(last_reader(owner_block), b)
          end do
       end do
    end if
    ! THE TAYLOR STATE MARCH: block outer, order inner. A block's tower
    ! of order k reads the towers of the blocks storing its given
    ! instants and its own lower orders, all solved already. So the
    ! forward pass accumulates each block's contribution to every
    ! table block by block and deallocates a tower once its last reader
    ! has been solved; the live storage is the reach in blocks, for any
    ! horizon length. The reverse pass reads every tower again in its
    ! reverse pass, so it retains them all, and accumulates the
    ! functionals themselves alone.
    from_size = merge(0, order, present(by_order))
    to_size   = order
    if (.not. forward) to_size = merge(0, -1, present(by_order))
    live         = 0
    peak_storage = 0
    total        = 0
    do b = 1, nb
       call forward_block(chain, b, physics, functionals, degrees, design, nd, top, order, &
            & from_size, to_size, versions, u, w, live, total, peak_storage, table, by_order, &
            & node_measure, last_reader)
    end do
    call tally_order(0)
    if (present(tower_storage)) tower_storage = [peak_storage, total]
    if (forward) then
       if (present(sinks)) then
          error stop 'gti_chain: the sinks are checked on the reverse pass'
       end if
       if (present(by_order)) by_order(:, :, order) = table
       return
    end if
    if (present(sinks)) then
       allocate(sinks % fixed_rows(0:degrees - 1), source=0)
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
                call solve_linear(chain(b) % rows, unknowns, inputs, rhs(1:count, b), .true., &
                     & versions(b), one)
                lambda(1:count, b, i, rank, k) = one
                if (present(sinks)) then
                   call sink_residual(is_sink(1:count, b), diagonal(1:count, b), &
                        & rhs(1:count, b), one, sinks)
                end if
                call tally_leave()
                do p = 1, size(chain(b) % source_at)
                   owner_block = chain(b) % source_block(p)
                   if (owner_block > 0) then
                      rhs(chain(b) % source_at(p), owner_block) = rhs(chain(b) % source_at(p), owner_block) + one(p)
                   end if
                end do
             end do
          end do
       end do
       call tally_leave()
    end do
    call tally_order(0)
    allocate(every(nf, nd, multiset_count(nd, top)), source=0.0_dp)
    if (present(leibniz))  allocate(leibniz(nf, nd, multiset_count(nd, top), 0:top + 1), source=0.0_dp)
    if (present(by_order)) allocate(products(nf, nd, multiset_count(nd, top)))
    allocate(split(0:top + 1))
    do rank = 1, multiset_count(nd, top)
       s = multiset_of(rank, top, nd)
       do j = 1, nd
          do i = 1, nf
             l     = derivative_terms(0.0_dp, top + 1)
             split = 0.0_dp
             do b = 1, nb
                l = l + lagrangian_term(chain, b, physics, functionals(i), degrees, design, s, j, &
                     & w, lambda, u, nd, i, node_measure, split)
             end do
             every(i, j, rank) = mixed_partial(l)
             if (present(leibniz))  leibniz(i, j, rank, :) = split
             if (present(by_order)) products(i, j, rank)   = l
          end do
       end do
    end do
    allocate(table(nf, multiset_count(nd, order)))
    do rank = 1, multiset_count(nd, order)
       s = multiset_of(rank, order, nd)
       table(:, rank) = every(:, s(order), multiset_rank(s(1:order - 1), nd))
    end do
    if (present(entries)) entries = every
    if (present(by_order)) then
       by_order(:, :, order) = table
       ! the lower orders, read from the same products: a multiset of
       ! size m - 1 with the explicit design is a subset of one of size
       ! top containing it, and the coefficient there is its entry
       do m = 1, order - 1
          do rank = 1, multiset_count(nd, m)
             s = multiset_of(rank, m, nd)
             call embedding_of(s(1:m - 1), top, nd, rank_s, mask)
             do i = 1, nf
                by_order(i, rank, m) = coefficient(products(i, s(m), rank_s), ior(mask, 2**top))
             end do
          end do
       end do
    end if
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
    type(tangent_tower) , intent(in) :: w(:)
    real(dp)            , intent(in) :: u(:,:,:)
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
          state_seed(:, mask) = w(b) % w(1:count, multiset_rank(t, nd), size_of)
       end if
    end do
  end subroutine seeds_of
  pure function sorted(x) result(y)
    integer, intent(in) :: x(:)
    integer :: y(size(x))
    integer :: i, j, fixed
    y = x
    do i = 2, size(y)
       fixed = y(i)
       j    = i - 1
       do while (j >= 1)
          if (y(j) <= fixed) exit
          y(j + 1) = y(j)
          j        = j - 1
       end do
       y(j + 1) = fixed
    end do
  end function sorted
  function point_terms(rule, degrees, design, at, n, additional, state_seed, nu_seed) result(t)
    type(expression), intent(in) :: rule
    integer         , intent(in) :: degrees, at, n, additional
    real(dp)        , intent(in) :: design, state_seed(:, 0:), nu_seed(:)
    type(derivative_terms) :: t
    type(derivative_terms) :: q(0:degrees - 1), nu
    integer :: d, mask
    do d = 0, degrees - 1
       q(d) = derivative_terms(state_seed(at + d + 1, 0), n + additional)
       do mask = 1, 2**n - 1
          call q(d) % set_coefficient(mask, state_seed(at + d + 1, mask))
       end do
       if (additional > 0) call q(d) % set_direction(n + d + 1, 1.0_dp)
    end do
    nu = derivative_terms(design, n + additional)
    do mask = 1, 2**n - 1
       if (nu_seed(mask) /= 0.0_dp) call nu % set_coefficient(mask, nu_seed(mask))
    end do
    t = rule % at_instant(q, nu)
  end function point_terms
  !===================================================================!
  ! THE POINTS ONE STEP IS INTEGRATED OVER, AND THEIR WEIGHTS.
  !
  ! A STAGE FAMILY places its points inside the step. The stages are
  ! already computed and the tableau weights them, which is a rule of
  ! the tableau's order, so a multistage family requires no history.
  !
  !        |---- step k ----|
  !        x    o    o    o          o the stages, weighted by beta
  !
  ! A MULTISTEP FAMILY has no points inside the step, so the rule is
  ! placed on the instants the step reaches back over - which its own
  ! stencil already contains. The family returns their weights, and p
  ! instants give order p.
  !
  !   o----o----o----|---- step k ----|
  !   k-3  k-2  k-1                   k       the history, weighted by
  !                                           the family's own rule
  !===================================================================!

  subroutine quadrature_points(b, k, steps, offset, weight)
    type(chain_block)     , intent(in)  :: b
    integer               , intent(in)  :: k
    type(derivative_terms), intent(in)  :: steps(:)
    integer , allocatable , intent(out) :: offset(:)
    type(derivative_terms), allocatable, intent(out) :: weight(:)
    integer :: s, i, width, nodes
    if (.not. b % staged) then
       call b % scheme % step_quadrature(steps, k, weight)
       nodes = size(weight)
       allocate(offset(nodes))
       do i = 1, nodes
          if (k - i + 1 < 1) then
             error stop 'gti_chain: a quadrature reads instants the block stores'
          end if
          offset(i) = b % instants_at(k - i + 1)
       end do
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
       weight(i) = derivative_terms(b % scheme % stage_weight(i), steps(k))
    end do
  end subroutine quadrature_points

  !===================================================================!
  ! THE BLOCK'S STEPS AS DERIVATIVE TERMS, each storing the design
  ! partials of its own width. A quadrature weight built on a designed
  ! grid varies with the design, so the weights are terms and not
  ! numbers, and the width matches the quantity the weights multiply.
  !===================================================================!

  function stepped_terms(b, step_seed, n, additional) result(steps)
    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: step_seed(:,:)
    integer          , intent(in) :: n, additional
    type(derivative_terms), allocatable :: steps(:)
    integer :: k, mask
    allocate(steps(size(b % dt)))
    do k = 1, size(b % dt)
       steps(k) = derivative_terms(b % dt(k), n + additional)
       do mask = 1, 2**n - 1
          call steps(k) % set_coefficient(mask, step_seed(k, mask))
       end do
    end do
  end function stepped_terms
  function measure_terms(b, k, node, n, additional, step_seed, node_measure) result(t)
    type(chain_block), intent(in) :: b
    integer          , intent(in) :: k, node, n, additional
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
    t = derivative_terms(b % dt(k), n + additional)
    do mask = 1, 2**n - 1
       call t % set_coefficient(mask, step_seed(k, mask))
    end do
    t = measure * t
  end function measure_terms
  subroutine forcing_of(chain, b, physics, degrees, design, s, w, u, nd, r)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b
    type(expression)    , intent(in) :: physics
    integer             , intent(in) :: degrees, s(:), nd
    real(dp)            , intent(in) :: design
    type(tangent_tower), intent(in) :: w(:)
    real(dp), intent(in) :: u(:,:,:)
    real(dp), allocatable, intent(out) :: r(:)
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: fixed_rows(:)
    integer :: n, full, e, mask, p, row
    n    = size(s)
    full = 2**n - 1
    call seeds_of(chain, b, s, 0, .false., w, u, nd, state_seed, step_seed, nu_seed)
    call chain(b) % rows % rows_terms(chain(b) % scheme, chain(b) % dt, step_seed, tr, tc, tw)
    fixed_rows = chain(b) % rows % fixed_mask()
    at      = chain(b) % rows % points_at()
    allocate(r(size(state_seed, 1)), source=0.0_dp)
    do e = 1, size(tr)
       if (fixed_rows(tr(e))) cycle
       do mask = 0, full - 1
          r(tr(e)) = r(tr(e)) + tw(e, ieor(full, mask)) * state_seed(tc(e), mask)
       end do
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (fixed_rows(row)) cycle
       r(row) = r(row) + coefficient(point_terms(physics, degrees, design, at(p), n, 0, &
            & state_seed, nu_seed), full)
    end do
  end subroutine forcing_of
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
    logical , allocatable :: has_diagonal(:), fixed_rows(:)
    logical :: available
    integer :: n, e, p, d
    call frozen_at(b, design, unknowns, inputs)
    call b % rows % explicit_tangent(unknowns, b % rows % bind(inputs), 1, r, c, w, available)
    if (.not. available) then
       error stop 'gti_chain: the block tangent in the state is explicit'
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
    fixed_rows = b % rows % fixed_mask()
    do p = 1, n
       if (.not. is_sink(p)) cycle
       d = mod(p - 1, degrees)
       if (fixed_rows(p)) then
          sinks % fixed_rows(d) = sinks % fixed_rows(d) + 1
       else if (p > n - degrees) then
          sinks % last(d) = sinks % last(d) + 1
       else
          sinks % interior(d) = sinks % interior(d) + 1
       end if
    end do
  end subroutine sinks_of
  subroutine sink_residual(is_sink, diagonal, g, lambda, sinks)
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
  end subroutine sink_residual
  subroutine costate_rows(chain, b, physics, rule, degrees, design, s, w, lambda, u, nd, i, &
       & node_measure, g)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), nd, i
    type(expression)    , intent(in) :: physics, rule
    real(dp)            , intent(in) :: design
    type(tangent_tower), intent(in) :: w(:)
    real(dp), intent(in) :: lambda(:,:,:,:,0:), u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    real(dp), allocatable, intent(out) :: g(:)
    type(derivative_terms) :: t
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:), lam(:,:)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: fixed_rows(:)
    integer , allocatable :: offset(:)
    type(derivative_terms), allocatable :: beta(:), steps(:)
    integer :: n, full, e, mask, p, d, row, k, node, from, to, point, count, pt
    n     = size(s)
    full  = 2**n - 1
    count = chain(b) % rows % num_unknowns()
    call seeds_of(chain, b, s, 0, .true., w, u, nd, state_seed, step_seed, nu_seed)
    allocate(g(count), source=0.0_dp)
    call owned(chain, b, from, to)
    steps = stepped_terms(chain(b), step_seed, n, degrees)
    do k = from, to
       call quadrature_points(chain(b), k, steps, offset, beta)
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
    fixed_rows = chain(b) % rows % fixed_mask()
    at      = chain(b) % rows % points_at()
    do e = 1, size(tr)
       if (fixed_rows(tr(e))) cycle
       do mask = 1, full
          g(tc(e)) = g(tc(e)) - tw(e, mask) * lam(tr(e), ieor(full, mask))
       end do
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (fixed_rows(row)) cycle
       t = point_terms(physics, degrees, design, at(p), n, degrees, state_seed, nu_seed)
       do mask = 1, full
          do d = 0, degrees - 1
             g(at(p) + d + 1) = g(at(p) + d + 1) &
                  & - coefficient(t, ior(mask, shiftl(1, n + d))) * lam(row, ieor(full, mask))
          end do
       end do
    end do
  end subroutine costate_rows
  !===================================================================!
  ! THE LAGRANGIAN'S TERMS AT ONE BLOCK, for the multiset s of designs
  ! the state is differentiated along and the design j it is not:
  ! over n + 1 directions, the first n storing the state's
  ! derivatives along s and the last storing j alone,
  !
  !      L  =  F  -  sum over the rows  lambda_row * R_row
  !
  ! with lambda_row the costates along the subsets of s and constant
  ! along j. The coefficient of the full subset is the block's contribution
  ! to the table entry of order n + 1 along s with j; the coefficient
  ! of a subset containing j and m of the first n directions is its contribution
  ! of the entry of order m + 1 along that sub-multiset. One product
  ! per row yields every order. parts accumulates the same product
  ! split by the order of the costate's factor - the terms of Leibniz,
  ! a subset of size k representing its binomial count once - and the
  ! functional's own term at n + 1. A parts of any other length stops
  ! the program.
  !===================================================================!
  function lagrangian_term(chain, b, physics, rule, degrees, design, s, j, w, lambda, u, nd, i, &
       & node_measure, parts) result(l)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), j, nd, i
    type(expression)    , intent(in) :: physics, rule
    real(dp)            , intent(in) :: design
    type(tangent_tower), intent(in) :: w(:)
    real(dp), intent(in) :: lambda(:,:,:,:,0:), u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    real(dp), intent(inout)          :: parts(0:)
    type(derivative_terms) :: l
    type(derivative_terms) :: f
    type(derivative_terms), allocatable :: residual(:), costate(:), beta(:), steps(:)
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:), tw(:,:), lam(:,:)
    real(dp) :: along(0:2**(size(s) + 1) - 1), split(0:size(s) + 1)
    integer , allocatable :: tr(:), tc(:), at(:)
    logical , allocatable :: fixed_rows(:)
    integer , allocatable :: offset(:)
    integer :: n, fulln, e, p, row, k, node, from, to, point, pt, count
    n     = size(s)
    fulln = 2**n - 1
    count = chain(b) % rows % num_unknowns()
    if (size(parts) /= n + 2) then
       error stop 'gti_chain: the parts number the costate orders and the functional''s own term'
    end if
    call seeds_of(chain, b, s, j, .true., w, u, nd, state_seed, step_seed, nu_seed)
    f = derivative_terms(0.0_dp, n + 1)
    call owned(chain, b, from, to)
    steps = stepped_terms(chain(b), step_seed, n + 1, 0)
    do k = from, to
       call quadrature_points(chain(b), k, steps, offset, beta)
       do pt = 1, size(offset)
          do node = 1, chain(b) % nodes
             point = offset(pt) + (node - 1) * degrees
             f = f + beta(pt) * measure_terms(chain(b), k, node, n + 1, 0, step_seed, node_measure) &
                  & * point_terms(rule, degrees, design, point, n + 1, 0, state_seed, nu_seed)
          end do
       end do
    end do
    call chain(b) % rows % rows_terms(chain(b) % scheme, chain(b) % dt, step_seed, tr, tc, tw)
    fixed_rows = chain(b) % rows % fixed_mask()
    at      = chain(b) % rows % points_at()
    allocate(residual(count))
    residual = derivative_terms(0.0_dp, n + 1)
    do e = 1, size(tr)
       if (fixed_rows(tr(e))) cycle
       residual(tr(e)) = residual(tr(e)) + derivative_terms(tw(e, :)) * derivative_terms(state_seed(tc(e), :))
    end do
    do p = 1, size(at)
       row = at(p) + chain(b) % primary + 1
       if (fixed_rows(row)) cycle
       residual(row) = residual(row) + point_terms(physics, degrees, design, at(p), n + 1, 0, state_seed, nu_seed)
    end do
    call costates_at(chain, b, s, lambda, nd, i, lam)
    allocate(costate(count))
    do row = 1, count
       along          = 0.0_dp
       along(0:fulln) = lam(row, :)
       costate(row)   = derivative_terms(along)
    end do
    l            = f - inner_product(costate, residual, active=.not. fixed_rows)
    parts(n + 1) = parts(n + 1) + mixed_partial(f)
    do row = 1, count
       if (fixed_rows(row)) cycle
       split      = leibniz_parts(costate(row), residual(row))
       parts(0:n) = parts(0:n) - split(0:n)
    end do
  end function lagrangian_term
  real(dp) function functional_along(chain, b, rule, degrees, design, s, open, w, u, nd, &
       & node_measure) result(part)
    type(chain_block)   , intent(in) :: chain(:)
    integer             , intent(in) :: b, degrees, s(:), open, nd
    type(expression)    , intent(in) :: rule
    real(dp)            , intent(in) :: design
    type(tangent_tower), intent(in) :: w(:)
    real(dp), intent(in) :: u(:,:,:)
    real(dp), intent(in), optional   :: node_measure(:)
    type(derivative_terms) :: t
    real(dp), allocatable :: state_seed(:,:), step_seed(:,:), nu_seed(:)
    integer , allocatable :: offset(:)
    type(derivative_terms), allocatable :: beta(:), steps(:)
    integer :: n, full, k, node, from, to, point, pt
    n    = size(s) + merge(1, 0, open > 0)
    full = 2**n - 1
    call seeds_of(chain, b, s, open, .true., w, u, nd, state_seed, step_seed, nu_seed)
    part = 0.0_dp
    call owned(chain, b, from, to)
    steps = stepped_terms(chain(b), step_seed, n, 0)
    do k = from, to
       call quadrature_points(chain(b), k, steps, offset, beta)
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
       error stop 'gti_chain: a multiset contains designs of the tower'
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
  !===================================================================!
  ! A multiset of the given size containing the sub-multiset t - t with
  ! copies of the first design added - as its rank, and the positions
  ! t occupies in it as a mask, so that a coefficient along t is read
  ! from the terms formed along the larger multiset. A t larger than
  ! the size stops the program.
  !===================================================================!
  pure subroutine embedding_of(t, size_of, designs, rank, mask)
    integer, intent(in)  :: t(:), size_of, designs
    integer, intent(out) :: rank, mask
    integer :: padded(size_of)
    integer :: k, p
    if (size(t) > size_of) then
       error stop 'gti_chain: a sub-multiset is no larger than the multiset containing it'
    end if
    padded = sorted([t, (1, k = 1, size_of - size(t))])
    rank   = multiset_rank(padded, designs)
    mask   = 0
    p      = 1
    do k = 1, size(t)
       do while (padded(p) /= t(k))
          p = p + 1
       end do
       mask = ibset(mask, p - 1)
       p    = p + 1
    end do
  end subroutine embedding_of
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
    count = num_blocks * pass_substitutions(pass_of(1, 1, order), 1, 1, order)
  end function expansion_substitutions

  !===================================================================!
  ! A grid built for one functional, not the state alone.
  !
  ! Once the steps are given as a design - gti_expansion's
  ! design_of_steps - the reverse pass already returns dF/d(weight_i)
  ! for every step i, one adjoint solve over the whole trajectory. It
  ! is an adjoint-weighted sensitivity to the grid's position, not the
  ! residual-weighted local defect a Becker-Rannacher estimator forms
  ! from a comparison-order scheme - DIRK here has no embedded pair
  ! to form one from. The two are related, not the same object.
  !
  ! A weight is a fraction of a fixed duration, so F, as a function of
  ! the weights, is unchanged by scaling all of them by one factor: F
  ! is homogeneous of degree zero in w, and Euler's identity makes
  ! sum(w_i dF/dw_i) vanish identically for every grid - a
  ! reparametrisation identity, true of an accurate grid and an
  ! inaccurate one alike, and no acceptance check can be read from it.
  ! What remains once that direction is projected away is the part of
  ! dF/dw that does contain local information: the first-order change
  ! in F from assigning a step more of the duration and every other step less.
  !
  ! The grid is accepted once the spread of that projected gradient,
  ! scaled by one average step's duration, is within tolerance
  ! relative to F (relative) or as an absolute value (absolute) - the
  ! change in F from moving one such fraction from the least to the
  ! most sensitive step. On rejection, every step whose projected
  ! gradient is at least the mean is halved - the steps with more than
  ! the mean sensitivity - and the whole trajectory, adjoint included,
  ! is computed again.
  !===================================================================!

  function goal_oriented_partition(scheme, physics, functional, degrees, duration, &
       & lower, design, tolerance, relative, rejects) result(dt)

    class(family)   , intent(in)  :: scheme
    type(expression), intent(in)  :: physics, functional
    integer         , intent(in)  :: degrees
    real(dp)        , intent(in)  :: duration, lower(:), design, tolerance
    logical         , intent(in)  :: relative
    integer         , intent(out), optional :: rejects
    real(dp), allocatable :: dt(:)

    integer, parameter :: seed_instants = 9

    type(family_container)  :: schemes(1)
    type(expression)     :: functionals(1)
    type(chain_block), allocatable :: chain(:)
    type(expansion)  , allocatable :: tower
    integer , allocatable :: versions(:)
    real(dp), allocatable :: state(:), resolved(:), t(:), fvals(:,:), table(:,:), eta(:)
    logical , allocatable :: split(:)
    real(dp) :: achieved, f, e, threshold
    integer  :: attempt, rejected

    allocate(schemes(1) % scheme, source=scheme)
    functionals(1) = functional
    state = consistent_state(physics, degrees, lower, design)

    dt = spread(duration / real(seed_instants - 1, dp), 1, seed_instants - 1)

    rejected = 0
    attempt  = 0
    do
       attempt = attempt + 1

       call march_chain(schemes, [size(dt)], physics, degrees, designed_grid(duration), &
            & design, state, chain, tower, resolved, t, achieved, grid_design=dt)

       call chain_expansion(chain, tower, functionals, degrees, 0, fvals)
       f = fvals(0, 1)

       call chain_versions(chain, tower, functionals, degrees, versions)
       call chain_derivative(chain, tower, versions, functionals, degrees, 1, reverse_pass, table)

       ! the scale direction, projected away: see the banner of this function.
       eta = table(1, 2:size(table, 2))
       eta = eta - sum(eta * dt) / sum(dt * dt) * dt
       e   = (maxval(eta) - minval(eta)) * (duration / real(size(eta), dp))
       if (relative) e = e / max(abs(f), tiny(1.0_dp))

       if (e <= tolerance) exit

       rejected  = rejected + 1
       threshold = sum(eta) / real(size(eta), dp)
       split     = eta >= threshold
       dt        = halved(dt, split)

       if (attempt > 50) then
          error stop 'gti_chain: a goal-oriented grid remains above the tolerance after fifty attempts'
       end if
    end do

    if (present(rejects)) rejects = rejected

  end function goal_oriented_partition

  pure function halved(dt, split) result(refined)
    real(dp), intent(in) :: dt(:)
    logical , intent(in) :: split(:)
    real(dp), allocatable :: refined(:)
    integer :: i
    refined = [real(dp) ::]
    do i = 1, size(dt)
       if (split(i)) then
          refined = [refined, dt(i) / 2.0_dp, dt(i) / 2.0_dp]
       else
          refined = [refined, dt(i)]
       end if
    end do
  end function halved

end module gti_chain
module gti_driver
  use iso_fortran_env  , only : int64
  use util_precision   , only : dp
  use operation_action , only : jacobian_of
  use operation_grid   , only : grid, uniform_grid, random_grid, partitioned
  use gti_configuration, only : configuration, read_configuration, override, argument_values, &
       & names_config, names_setting
  use view_directed_stored, only : stored_directed_graph
  use field_stored     , only : stored_field
  use gti_march        , only : frozen_inputs
  use operation_family , only : family
  use operation_family      , only : bdf_family
  use operation_family      , only : adams_family
  use operation_family      , only : implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage, &
       & newmark_family, taylor_newmark_family
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
    character(len=256), allocatable :: given(:)
    character(len=:), allocatable :: name
    integer :: i
    name  = default_name
    given = argument_values([names_config])
    if (size(given) > 0) name = trim(given(size(given)))
    call read_configuration(name, cfg)
    given = argument_values([names_setting])
    do i = 1, size(given)
       call override(cfg, given(i))
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
    type(stored_field), allocatable :: inputs(:)
    call frozen_inputs(chain(1) % state, design, chain(1) % rows % num_points(), unknowns, inputs)
    call jacobian_of(chain(1) % rows, unknowns, inputs, chain(1) % rows % num_unknowns(), &
         & unknowns % vertex_set(), a)
  end subroutine dense_jacobian
  subroutine family_named(name, order, scheme, passes_check)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: passes_check
    passes_check = .true.
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
          passes_check = .false.
       end select
    case ('newmark')
       select case (order)
       case (1)
          allocate(scheme, source=newmark_family(0.0_dp, 0.0_dp, 1))
       case (2)
          allocate(scheme, source=newmark_family(0.25_dp, 0.5_dp, 2))
       case (3)
          allocate(scheme, source=newmark_family(1.0_dp / 12.0_dp, 0.5_dp, 3))
       case default
          passes_check = .false.
       end select
    case ('taylor-newmark')
       if (order == 1) then
          allocate(scheme, source=taylor_newmark_family())
       else
          passes_check = .false.
       end if
    case default
       passes_check = .false.
    end select
  end subroutine family_named
  subroutine functional_named(name, degree, rule, passes_check)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: degree
    type(expression), intent(out) :: rule
    logical         , intent(out) :: passes_check
    passes_check = .true.
    select case (name)
    case ('energy')
       rule = van_der_pol_energy(degree)
    case ('dissipation')
       rule = van_der_pol_dissipation(degree)
    case default
       passes_check = .false.
    end select
  end subroutine functional_named
end module gti_driver
module gti_demos
  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use graph_fractal         , only : graph, branch, known_branch
  use view_sequence         , only : sequence_first, sequence_rest, sequence_num_elements, &
       & sequence_element
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
  use view_directed_connectivity, only : connectivity_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field, typed_field_domain
  use operation_action      , only : variation, sweep_design_partial => design_partial
  use gti_sweeps            , only : functional_of, functional_gradient
  use operation_stencil     , only : stencil
  use operation_scheme_stencil, only : derived_constraints
  use operation_family      , only : family
  use operation_family      , only : bdf_family
  use operation_family      , only : adams_family
  use operation_family      , only : dirk_family, implicit_midpoint, &
       & crouzeix_two_stage, crouzeix_three_stage, hairer_wanner_five_stage, &
       & newmark_family, taylor_newmark_family
  use operation_grid        , only : grid, uniform_grid, random_grid, designed_grid, partition
  use operation_coupling    , only : weights_of, coupling_inputs
  use operation_weight      , only : scheme_weight
  use operation_expression  , only : expression
  use operation_domain      , only : continuous_domain, discrete_domain
  use operation_temporal_minimization, only : temporal_minimizer
  use gti_physics           , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use operation_minimization, only : relative, by_rate
  use gti_expansion         , only : expansion, family_container
  use gti_block             , only : block_residual
  use gti_march             , only : block_from, solved, unknowns_graph, &
       & horizon_bounds, set_stopping, consistent_state, imbalance, by_tangent, &
       & by_adjoint, next_version, instants_at_of
  use gti_adaptive          , only : adaptive_partition
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & chain_versions, chain_derivative, first_of, instant_components, &
       & asymmetry, multiset_count, multiset_rank, multiset_of
  use gti_chain             , only : chain_incidence
  use operation_driver      , only : rule_graph, data_graph, pairing
  use gti_sweeps            , only : pass_of, forward_pass, reverse_pass, choose
  use gti_driver            , only : clock, cosine, dense_jacobian, family_named, settings
  use gti_configuration     , only : configuration, argument_values, names_config, names_setting, &
       & names_demo
  use view_read_write       , only : bipartite_digraph, FIRST_PART, SECOND_PART
  use operation_driver      , only : driver
  use view_directed         , only : forward
  implicit none
  private
  public :: demo_requested, run_demo
  ! which part of the bipartite digraph this demonstration reads as which
  integer, parameter :: BLOCKS_PART = FIRST_PART
  integer, parameter :: DATA_PART   = SECOND_PART
  character(len=5) , parameter :: family_names(3) = ['bdf  ', 'adams', 'dirk ']
  character(len=24), parameter :: demo_names(27) = [character(len=24) :: &
       & 'adaptive_grid', 'assembled_tower', 'chained_horizon', &
       & 'constraint_rows', 'coupling_relation', 'expansion_check', &
       & 'family_coefficients', 'function_identities', 'grid_design_check', &
       & 'jacobian_shape', 'lagrangian_expansion', 'level_maps', 'level_shape', 'marched_block', &
       & 'marched_horizon', 'marched_stages', 'memory_shape', &
       & 'transfer_offsets', 'randomized_checks', 'read_write_graph', &
       & 'order_of_accuracy', 'scheme_weights', 'sensitivity', 'solve_cost', 'taylor_state', &
       & 'tolerance_form', 'transposed_dependencies']
contains
  logical function demo_requested() result(yes)
    yes = size(argument_values([names_demo])) > 0
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
    case ('transfer_offsets')
       call demo_transfer_offsets()
    case ('jacobian_shape')
       call demo_jacobian_shape()
    case ('lagrangian_expansion')
       call demo_lagrangian_expansion()
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
    case ('order_of_accuracy')
       call demo_order_of_accuracy()
    case ('read_write_graph')
       call demo_read_write_graph()
    case ('transposed_dependencies')
       call demo_transposed_dependencies()
    case ('scheme_weights')
       call demo_scheme_weights()
    case ('sensitivity')
       call demo_sensitivity()
    case ('solve_cost')
       call demo_solve_cost()
    case ('taylor_state')
       call demo_taylor_state()
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
    character(len=256), allocatable :: given(:)
    integer :: j
    name  = ''
    given = argument_values([names_demo])
    if (size(given) > 0) name = trim(given(size(given)))
    do j = 1, len(name)
       if (name(j:j) == '-') name(j:j) = '_'
    end do
  end function demo_name
  ! the demonstration's own arguments: every argument that selects no
  ! demonstration, in the order given; blank beyond the last
  subroutine demo_argument(which, argument)
    integer, intent(in) :: which
    character(len=*), intent(out) :: argument
    character(len=256), allocatable :: given(:)
    given    = argument_values([names_config, names_setting])
    argument = ''
    if (which <= size(given)) argument = given(which)
  end subroutine demo_argument
  ! the k-th argument read as a number, or the default where none is given
  real(dp) function demo_real(k, default) result(x)
    integer , intent(in) :: k
    real(dp), intent(in) :: default
    character(len=32) :: argument
    x = default
    call demo_argument(k, argument)
    if (len_trim(argument) > 0) read(argument, *) x
  end function demo_real
  subroutine list_demos()
    integer :: i
    write(*,'(a)') ' demos:'
    do i = 1, size(demo_names)
       write(*,'(a)') '   ' // trim(demo_names(i))
    end do
  end subroutine list_demos
  type(family_container) function stored_family(scheme) result(fixed)
    class(family), intent(in) :: scheme
    call set_family(fixed, scheme)
  end function stored_family
  subroutine set_family(fixed, scheme)
    type(family_container), intent(inout) :: fixed
    class(family)     , intent(in)    :: scheme
    if (allocated(fixed % scheme)) deallocate(fixed % scheme)
    allocate(fixed % scheme, source=scheme)
  end subroutine set_family
  function cosine_history(scheme, degrees, t) result(fixed)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees
    real(dp)     , intent(in) :: t(:)
    real(dp), allocatable :: fixed(:)
    integer :: k, d
    fixed = [((cosine(d, t(k)), d = 0, degrees - 1), &
         & k = 1, scheme % history_depth(degrees - 1))]
  end function cosine_history
  subroutine cosine_partition(scheme, degrees, duration, instants, fixed, dt, t)
    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, instants
    real(dp)     , intent(in) :: duration
    real(dp), allocatable, intent(out) :: fixed(:), dt(:), t(:)
    call partition(duration, instants, dt, t)
    fixed = cosine_history(scheme, degrees, t)
  end subroutine cosine_partition
  !===================================================================!
  ! VAN DER POL MARCHED FROM THE COSINE HISTORY on a uniform grid of
  ! the duration: the first family's history depth of instants is
  ! sampled from the cosine, and the chain adds the rest.
  !===================================================================!
  subroutine marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
    type(family_container), intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    real(dp)              , intent(in) :: duration, design
    type(chain_block), allocatable, intent(out)   :: chain(:)
    type(expansion)  , allocatable, intent(inout), target :: tower
    real(dp)              , intent(out) :: achieved
    real(dp), allocatable :: fixed(:), dt(:), t(:)
    call cosine_partition(schemes(1) % scheme, degrees, duration, sum(added), fixed, dt, t)
    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, uniform_grid(duration), &
         & design, fixed, chain, tower, dt, t, achieved)
  end subroutine marched_cosine
  ! the first derivatives of the functionals in the designs, by the
  ! forward pass and by the reverse pass
  subroutine directions(chain, tower, functionals, degrees, tangent, adjoint)
    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    type(expression) , intent(in) :: functionals(:)
    integer          , intent(in) :: degrees
    real(dp), allocatable, intent(out) :: tangent(:,:), adjoint(:,:)
    integer, allocatable :: versions(:)
    call chain_versions(chain, tower, functionals, degrees, versions)
    call chain_derivative(chain, tower, versions, functionals, degrees, 1, forward_pass, tangent)
    call chain_derivative(chain, tower, versions, functionals, degrees, 1, reverse_pass, adjoint)
  end subroutine directions
  ! the named family at the given order; a family with no scheme at
  ! that order stops the program
  function container_named(family_of, order) result(h)
    character(len=*), intent(in) :: family_of
    integer         , intent(in) :: order
    type(family_container) :: h
    logical :: passes_check
    call family_named(family_of, order, h % scheme, passes_check)
    if (.not. passes_check) error stop 'gti_demos: the named family has a scheme at that order'
  end function container_named
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
    subroutine on_grid(scheme, dt, f, forward, reverse)
      class(family), intent(in)  :: scheme
      real(dp)     , intent(in)  :: dt(:)
      real(dp)     , intent(out) :: f, forward, reverse
      type(family_container)     :: schemes(1)
      type(expression)       :: functionals(1)
      type(chain_block), allocatable :: chain(:)
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
      call directions(chain, tower, functionals, degrees, df, other)
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
      logical  :: passes_check
      call family_named('dirk', order, scheme, passes_check)
      if (.not. passes_check) error stop 'adaptive_grid: order two, three or four'
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
    type(family_container) :: schemes(2)
    schemes = [stored_family(bdf_family(2)), stored_family(adams_family(3))]
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
      type(branch) :: members
      integer :: k
      write(*,'(a,a,a,a)') repeat('   ', depth + 1), tower % label_of(g), &
           & status(tower, g), extent(tower, g)
      if (level_couples(g)) then
         coupling => level_coupling(g)
         write(*,'(a,a,a,a)') repeat('   ', depth + 2), tower % label_of(coupling), &
              & status(tower, coupling), extent(tower, coupling)
      end if
      if (level_is_leaf(g)) return
      members = level_members(g)
      do k = 1, sequence_num_elements(members)
         call show(tower, sequence_element(members, k), depth + 1)
      end do
    end subroutine show
    function status(tower, g) result(label)
      type(expansion), intent(in) :: tower
      type(graph)    , intent(in) :: g
      character(len=:), allocatable :: label
      real(dp), allocatable :: x(:)
      select case (tower % status_of(g))
      case (VALUE_KNOWN)
         call tower % value_of(g, x)
         label = '   stores ' // count_of(size(x))
      case (VALUE_UNKNOWN)
         label = '   not yet known'
      case (VALUE_UNATTACHED)
         label = ''
      case default
         error stop 'assembled_tower: a value status is one of the three'
      end select
    end function status
    function extent(tower, g) result(label)
      type(expansion), intent(in) :: tower
      type(graph)    , intent(in) :: g
      character(len=:), allocatable :: label
      label = ''
      if (tower % extent_of(g) > 0) label = '   extent ' // count_of(tower % extent_of(g))
    end function extent
    function count_of(n) result(name)
      integer, intent(in) :: n
      character(len=:), allocatable :: name
      character(len=12) :: buffer
      write(buffer,'(i0)') n
      name = trim(buffer)
    end function count_of
    subroutine grid_partials()
      integer , parameter :: num_instants = 6
      real(dp), parameter :: delta = 1.0e-6_dp
      type(grid) :: steps
      type(stored_directed_graph) :: instants
      type(stored_field) :: design_field, direction
      type(typed_field_domain) :: instant_scalars
      class(field), allocatable :: out
      real(dp) :: design(num_instants - 1), v(num_instants - 1)
      real(dp), allocatable :: dt(:), exact(:), plus(:), minus(:)
      design = [1.0_dp, 2.0_dp, 1.5_dp, 0.5_dp, 3.0_dp]
      v      = 0.0_dp
      v(2)   = 1.0_dp
      steps    = designed_grid(duration)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      instant_scalars = typed_field_domain(instants % vertex_set(), size(design))
      design_field    = instant_scalars % design(design)
      direction       = instant_scalars % direction(v)
      call steps % apply(instants, steps % bind([design_field]), out)
      call out % real_vector(dt)
      call steps % partial_action(instants, steps % bind([design_field]), &
           & [variation(steps % argument(1), direction)], out)
      call out % real_vector(exact)
      call differenced(steps, instants, design_field, design, v, delta, plus, minus)
      call show_grid(dt, exact, (plus - minus) / (2.0_dp * delta))
    end subroutine grid_partials
    !-----------------------------------------------------------------!
    ! The quadrature kind of the same map: the n-point rule integrates
    ! t to every power below 2n exactly over [0, span]. The bound each
    ! difference is bound to is the arithmetic's: n products, each
    ! within one epsilon of the integral's size.
    !-----------------------------------------------------------------!

    subroutine gauss_exactness()

      use operation_grid, only : gauss_grid, partitioned

      real(dp), parameter :: span = 2.0_dp
      real(dp), allocatable :: w(:), t(:)
      real(dp) :: quadrature, integral, rounding_bound, maximum_departure
      integer :: n, m, k

      write(*,'(a)') ' '
      write(*,'(a)') ' the gauss kind of the grid, exact to twice its points less one'
      do n = 3, 5, 2
         call partitioned(gauss_grid(span), n, w, t)
         maximum_departure = 0.0_dp
         do m = 0, 2 * n - 1
            quadrature = 0.0_dp
            do k = 1, n
               quadrature = quadrature + w(k) * t(k) ** m
            end do
            integral = span ** (m + 1) / real(m + 1, dp)
            rounding_bound    = real(n, dp) * epsilon(1.0_dp) * integral
            maximum_departure    = max(maximum_departure, abs(quadrature - integral) / rounding_bound)
         end do
         write(*,'(a,i0,a,i0,a,f8.3,a,es9.2)') '   ', n, ' points, powers 0 to ', &
              & 2 * n - 1, ':  sum of weights - span ', sum(w) - span, &
              & '   largest difference over its bound ', maximum_departure
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
      write(*,'(a,es10.2)')   '   sum of the partial, theory 0 ', sum(exact)
    end subroutine show_grid
    subroutine differenced(steps, instants, design_field, design, v, delta, plus, minus)
      type(grid)                 , intent(in)    :: steps
      type(stored_directed_graph), intent(in)    :: instants
      type(stored_field)         , intent(inout) :: design_field
      real(dp)                   , intent(in)    :: design(:), v(:), delta
      real(dp), allocatable      , intent(out)   :: plus(:), minus(:)
      class(field), allocatable :: out
      call design_field % set_real_vector(design + delta * v)
      call steps % apply(instants, steps % bind([design_field]), out)
      call out % real_vector(plus)
      call design_field % set_real_vector(design - delta * v)
      call steps % apply(instants, steps % bind([design_field]), out)
      call out % real_vector(minus)
      call design_field % set_real_vector(design)
    end subroutine differenced
    subroutine stage_block()
      integer , parameter :: num_instants = 3
      real(dp), parameter :: gamma = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp
      type(expansion) :: staged
      type(family_container) :: schemes(1)
      type(graph), pointer :: g, coupling
      real(dp), allocatable :: w(:)
      real(dp) :: step
      schemes = [stored_family(crouzeix_two_stage())]
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
      write(*,'(a,2f10.5)') '   the tableau weights b_1 b_2      ', [0.5_dp, 0.5_dp]
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
    call splitting_is_invariant()
    call across_families('bdf 2 then crouzeix two-stage', &
         & stored_family(bdf_family(2)), stored_family(crouzeix_two_stage()))
    call across_families('crouzeix two-stage then bdf 2', &
         & stored_family(crouzeix_two_stage()), stored_family(bdf_family(2)))
    call across_families('adams 3 then crouzeix two-stage', &
         & stored_family(adams_family(3)), stored_family(crouzeix_two_stage()))
  contains
    subroutine expanded(schemes, added, design, f)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp)           , intent(in) :: design
      real(dp), allocatable, intent(out) :: f(:)
      real(dp), allocatable :: table(:,:)
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      real(dp) :: achieved
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(state_degree)], degrees, max_order, table)
      allocate(f(lbound(table, 1):ubound(table, 1)))
      f = table(:, 1)
    end subroutine expanded
    subroutine splitting_is_invariant()
      type(family_container) :: whole(1), split(2)
      real(dp), allocatable :: f_whole(:), f_split(:)
      whole = [stored_family(bdf_family(2))]
      split = [stored_family(bdf_family(2)), stored_family(bdf_family(2))]
      call expanded(whole, [20], 1.0_dp, f_whole)
      call expanded(split, [10, 10], 1.0_dp, f_split)
      write(*,'(a)')        ' bdf 2 over twenty instants, in one block and in two'
      write(*,'(a,4f14.8)') '   whole                    ', f_whole
      write(*,'(a,4f14.8)') '   split                    ', f_split
      write(*,'(a,4es14.2)')'   difference               ', abs(f_whole - f_split)
    end subroutine splitting_is_invariant
    subroutine across_families(title, one, two)
      character(len=*)   , intent(in) :: title
      type(family_container), intent(in) :: one, two
      type(family_container) :: schemes(2)
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
      write(*,'(a,14x,3es14.2)')'   difference               ', &
           & abs(f(1:max_order) - differenced)
    end subroutine across_families
  end subroutine demo_chained_horizon

  !===================================================================!
  ! DOES A CHAIN OF DIFFERENT SCHEMES PRESERVE THE ORDER OF ACCURACY,
  ! AND IS THE ORDER PRESERVED ON EVERY DERIVATIVE OF THE TIME FUNCTIONAL?
  !
  ! Halving the step should shrink the error by a factor 2**p, where p
  ! is the order of the scheme. So marching the same problem on two
  ! grids and comparing the errors yields p:
  !
  !          observed order = log2( error(h) / error(h/2) )
  !
  ! There is no exact solution to measure against, so a reference grid
  ! far finer than any of the three measured grids is the reference.
  !
  ! Each window retains the same share of the horizon on every grid, so
  ! refining never moves a window boundary. Without that the measured
  ! order would depend on the boundary positions and not on the schemes.
  !
  ! A chain has order no higher than its lowest-order window, so the
  ! formal order of a chain is the minimum among the schemes in it.
  !
  ! WHY THE CHAINS HERE ARE ALL RUNGE-KUTTA. The time functional is
  ! integrated by one point per step for a multistep family, which is
  ! a rectangle rule and first order whatever the scheme is. Any chain
  ! containing such a window measures order one and gives no information
  ! about the chaining. The multistep families are still shown, at the
  ! end and on their own, so that limit is displayed.
  !===================================================================!

  subroutine demo_order_of_accuracy()

    implicit none
    type(configuration) :: cfg
    integer  :: state_degree, degrees, max_order, per_window, grids
    real(dp) :: duration, ratio, finer

    real(dp) :: allowed, settled

    ! THE GRIDS ARE THE CONFIGURATION'S, so the refinement ratio can be
    ! changed and the same table recomputed. An order that changes
    ! with the ratio is not asymptotic.
    call settings('order_of_accuracy', cfg)
    state_degree = cfg % state_degree
    degrees      = state_degree + 1
    max_order    = cfg % max_derivative_degree
    duration     = cfg % time_duration
    per_window   = cfg % coarsest_instants
    grids        = cfg % refinement_grids
    ratio        = cfg % refinement_ratio
    finer        = cfg % reference_refinement
    allowed      = cfg % order_tolerance
    settled      = cfg % spread_tolerance
    if (ratio <= 1.0_dp) then
       error stop 'gti_demos: a refinement ratio is above one'
    end if
    if (grids < 3) then
       error stop 'gti_demos: three grids at least, or no spread can be computed'
    end if
    if (finer <= 1.0_dp) then
       error stop 'gti_demos: a reference grid is finer than the finest measured'
    end if

    write(*,'(a)') ' '
    write(*,'(a)') ' DOES CHAINING DIFFERENT SCHEMES PRESERVE THE ORDER OF ACCURACY?'
    write(*,'(a)') ' '
    write(*,'(a)') ' THREE DISTINCT NUMBERS APPEAR BELOW, AND NO TWO OF THEM COINCIDE.'
    write(*,'(a)') ' '
    write(*,'(a)') '   p, THE ORDER OF ACCURACY OF THE SCHEME. The rate at which the scheme''s'
    write(*,'(a)') '      own error decreases as the step h decreases: the error is O(h**p).'
    write(*,'(a)') '      BDF-2 has p = 2 and SDIRK-4 has p = 4. A chain has order no higher'
    write(*,'(a)') '      than its lowest-order window, so a chain''s p is the minimum among them.'
    write(*,'(a)') '      Written below as "scheme order p". p is independent of r.'
    write(*,'(a)') ' '
    write(*,'(a)') '   r, THE DEGREE OF THE DESIGN DERIVATIVE. Which derivative of the time'
    write(*,'(a)') '      functional F with respect to the design nu is being measured:'
    write(*,'(a)') '      r = 0 is F itself, r = 1 is dF/dnu, r = 2 is d2F/dnu2, and so on'
    write(*,'(a)') '      to r = 6. Nothing about r refers to h or to any discretisation.'
    write(*,'(a)') '      It is the column heading below.'
    write(*,'(a)') ' '
    write(*,'(a)') '   THE OBSERVED ORDER, which is p measured rather than formal, and'
    write(*,'(a)') '      measured separately for the derivative of each degree r. That is'
    write(*,'(a)') '      the number in the table. The check here is whether it equals p'
    write(*,'(a)') '      at every r, or only at the low ones.'
    write(*,'(a)') ' '
    write(*,'(a)') '   A FOURTH NUMBER is fixed throughout and appears nowhere below. Van'
    write(*,'(a,i0,a,i0,a)') '      der Pol is of order ', state_degree, &
         & ' in time, so every instant stores ', degrees, ' state components'
    write(*,'(a)') '      (q, q-dot, q-double-dot). That is neither an order of accuracy'
    write(*,'(a)') '      nor a derivative degree, and it is fixed here.'
    write(*,'(a)') ' '
    write(*,'(a)') ' HOW THE OBSERVED ORDER IS MEASURED'
    write(*,'(a,i0,a,i0,a)') '   Each chain is marched on ', grids, &
         & ' grids, the coarsest having ', per_window, ' instants per window,'
    write(*,'(a,f0.3,a)') '   each grid after it ', ratio, ' times the previous,'
    write(*,'(a,f0.2,a)') '   and all of them against a reference grid ', finer, &
         & ' times finer than the finest.'
    write(*,'(a)') '   An error that is one power of h obeys e = C h**q, so log e is a'
    write(*,'(a)') '   straight line in log N, of slope -q. The order printed'
    write(*,'(a)') '   is that slope fitted by least squares over every grid at once,'
    write(*,'(a)') '   separately for each derivative degree r.'
    write(*,'(a)') '   Beside it is the spread of the slopes taken pair by pair. Where'
    write(*,'(a)') '   the error is a single power they all agree and the spread is near'
    write(*,'(a)') '   zero. Where a leading coefficient changes sign inside the range of'
    write(*,'(a)') '   grids the error passes through zero, the points leave the line, and'
    write(*,'(a)') '   the spread records this instead of accepting a meaningless slope.'
    write(*,'(a)') '   Every window retains its share of the horizon as the grid refines, so'
    write(*,'(a)') '   no window boundary moves and the measured order depends on the schemes alone.'
    write(*,'(a,f4.2,a)') '   An order within ', allowed, ' of p is counted as reaching p, and a'
    write(*,'(a,f4.2,a)') '   spread above ', settled, ' is counted as no single power, printed "?".'
    write(*,'(a)') '   A dash means an error reached round-off, and no order can be measured.'
    write(*,'(a)') ' '
    write(*,'(a)') '   THE RATIO IS THE CONFIGURATION''S. Measuring the same order at a second'
    write(*,'(a)') '   ratio shows the measurement is the scheme''s order and not an'
    write(*,'(a)') '   artefact of the grid positions:'
    write(*,'(a)') '     --demo=order_of_accuracy --refinement_ratio=1.618 --refinement_grids=5'
    write(*,'(a)') ' '
    write(*,'(a)') ' WHY THE FAMILIES ARE GROUPED'
    write(*,'(a)') '   The time functional is integrated over each step by the points the'
    write(*,'(a)') '   scheme places there. A Runge-Kutta step contributes its stages and the'
    write(*,'(a)') '   tableau weights; a multistep step contributes the instants of its own'
    write(*,'(a)') '   stencil and the interpolatory weights on them. Both rules have the'
    write(*,'(a)') '   order p of the scheme they belong to.'

    write(*,'(a)') ' '
    write(*,'(a)') ' ==== DIRK ALONE ===='
    call order_row('DIRK-2   implicit midpoint', [2], [stored_family(implicit_midpoint())])
    call order_row('DIRK-3   Crouzeix two-stage', [3], [stored_family(crouzeix_two_stage())])
    call order_row('DIRK-4   Crouzeix three-stage', [4], [stored_family(crouzeix_three_stage())])
    call order_row('SDIRK-4  Hairer-Wanner five-stage', [4], &
         & [stored_family(hairer_wanner_five_stage())])

    write(*,'(a)') ' '
    write(*,'(a)') ' ==== BDF ALONE ===='
    call order_row('BDF-1', [1], [stored_family(bdf_family(1))])
    call order_row('BDF-2', [2], [stored_family(bdf_family(2))])
    call order_row('BDF-3', [3], [stored_family(bdf_family(3))])
    call order_row('BDF-4', [4], [stored_family(bdf_family(4))])

    write(*,'(a)') ' '
    write(*,'(a)') ' ==== ABM ALONE ===='
    call order_row('ABM-2', [2], [stored_family(adams_family(2))])
    call order_row('ABM-3', [3], [stored_family(adams_family(3))])
    call order_row('ABM-4', [4], [stored_family(adams_family(4))])

    write(*,'(a)') ' '
    write(*,'(a)') ' ==== CHAINS ACROSS THE THREE FAMILIES ===='
    call order_row('SDIRK-4 -> BDF-4 -> ABM-4   (the composite of the paper)', [4, 4, 4], &
         & [stored_family(hairer_wanner_five_stage()), stored_family(bdf_family(4)), &
         &  stored_family(adams_family(4))])
    call order_row('DIRK-3 -> BDF-2', [3, 2], &
         & [stored_family(crouzeix_two_stage()), stored_family(bdf_family(2))])
    call order_row('BDF-2 -> DIRK-3', [2, 3], &
         & [stored_family(bdf_family(2)), stored_family(crouzeix_two_stage())])
    call order_row('ABM-3 -> DIRK-4', [3, 4], &
         & [stored_family(adams_family(3)), stored_family(crouzeix_three_stage())])
    call order_row('BDF-2 -> ABM-3 -> DIRK-4 -> SDIRK-4', [2, 3, 4, 4], &
         & [stored_family(bdf_family(2)), stored_family(adams_family(3)), &
         &  stored_family(crouzeix_three_stage()), stored_family(hairer_wanner_five_stage())])

    write(*,'(a)') ' '
    write(*,'(a)') ' ==== CHAINS WITHIN DIRK, WHERE THE QUADRATURE IS NOT THE LIMIT ===='
    call order_row('DIRK-3 -> DIRK-2', [3, 2], &
         & [stored_family(crouzeix_two_stage()), stored_family(implicit_midpoint())])
    call order_row('DIRK-3 -> DIRK-4', [3, 4], &
         & [stored_family(crouzeix_two_stage()), stored_family(crouzeix_three_stage())])
    call order_row('DIRK-4 -> SDIRK-4', [4, 4], &
         & [stored_family(crouzeix_three_stage()), stored_family(hairer_wanner_five_stage())])
    call order_row('DIRK-4 -> DIRK-3 -> DIRK-2', [4, 3, 2], &
         & [stored_family(crouzeix_three_stage()), stored_family(crouzeix_two_stage()), &
         &  stored_family(implicit_midpoint())])
    call order_row('DIRK-2 -> DIRK-3 -> DIRK-4 -> SDIRK-4', [2, 3, 4, 4], &
         & [stored_family(implicit_midpoint()), stored_family(crouzeix_two_stage()), &
         &  stored_family(crouzeix_three_stage()), stored_family(hairer_wanner_five_stage())])
    call order_row('SDIRK-4 -> DIRK-4 -> DIRK-3 -> DIRK-2 -> DIRK-4', [4, 4, 3, 2, 4], &
         & [stored_family(hairer_wanner_five_stage()), stored_family(crouzeix_three_stage()), &
         &  stored_family(crouzeix_two_stage()), stored_family(implicit_midpoint()), &
         &  stored_family(crouzeix_three_stage())])

  contains

    !----------------------------------------------------------------!
    ! One chain on one grid: the time functional and its derivatives
    ! in the design, with the windows in a fixed share of the whole.
    !----------------------------------------------------------------!

    subroutine on_grid(schemes, instants, f)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: instants
      real(dp), allocatable, intent(out) :: f(:)
      real(dp), allocatable :: table(:,:)
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer , allocatable :: added(:)
      real(dp) :: achieved
      integer :: windows, share
      windows = size(schemes)
      share   = instants / windows
      allocate(added(windows))
      added    = share
      added(1) = instants - share * (windows - 1)
      call marched_cosine(schemes, added, degrees, duration, 1.0_dp, chain, tower, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(state_degree)], degrees, &
           & max_order, table)
      allocate(f(lbound(table, 1):ubound(table, 1)))
      f = table(:, 1)
    end subroutine on_grid

    subroutine order_row(title, formal, schemes)
      character(len=*)   , intent(in) :: title
      integer            , intent(in) :: formal(:)
      type(family_container), intent(in) :: schemes(:)
      real(dp), allocatable :: f(:), reference(:)
      real(dp), allocatable :: discrepancy(:,:), counts(:)
      real(dp) :: fitted(0:max_order), spread(0:max_order)
      logical  :: readable(0:max_order)
      character(len=:), allocatable :: line, orders, spreads, outcomes
      character(len=16) :: cell
      integer :: g, m, windows, instants, expected, kept, counted, reaches

      windows  = size(schemes)
      expected = minval(formal)

      ! A WINDOW TOO SHORT FOR THE FAMILY ASSIGNED TO IT cannot be
      ! marched at all, and the coarsest grid determines this.
      ! A configuration that lowers coarsest_instants is reported the
      ! chains that are excluded rather than stopped.
      reaches = 0
      do g = 1, windows
         reaches = max(reaches, schemes(g) % scheme % history_depth(degrees - 1))
      end do
      if (per_window <= reaches) then
         write(*,'(a)') ' '
         write(*,'(a,a,i0,a,i0,a)') '   ' // title, '   -   not measured: the coarsest grid gives ', &
              & per_window, ' instants per window and a family here has history depth ', &
              & reaches, '.'
         write(*,'(a,i0,a)') '     raise coarsest_instants above ', reaches, ' to measure this chain.'
         return
      end if

      allocate(discrepancy(grids, 0:max_order), counts(grids))

      ! the grids, and the reference beyond the finest of them
      do g = 1, grids
         counts(g) = real(windows * per_window, dp) * ratio ** (g - 1)
      end do
      call on_grid(schemes, nint(counts(grids) * finer), reference)
      do g = 1, grids
         instants = nint(counts(g))
         counts(g) = real(instants, dp)
         call on_grid(schemes, instants, f)
         do m = 0, max_order
            discrepancy(g, m) = abs(f(m) - reference(m))
         end do
         deallocate(f)
      end do

      do m = 0, max_order
         readable(m) = all(discrepancy(:, m) > 0.0_dp)
         fitted(m)   = 0.0_dp
         spread(m)   = 0.0_dp
         if (readable(m)) call order_of_grids(counts, discrepancy(:, m), fitted(m), spread(m))
      end do

      kept    = count(readable .and. spread <= settled &
           &          .and. fitted >= real(expected, dp) - allowed)
      counted = count(readable .and. spread <= settled)

      write(*,'(a)') ' '
      write(*,'(a,a,i0,a,i0)') '   ' // title, &
           & '   -   scheme order p = ', expected, ', windows ', windows
      write(cell,'(i0)') expected
      line     = '     r, derivative degree of F in nu  '
      orders   = '     order fitted over the grids     '
      spreads  = '     spread of pairwise slopes       '
      outcomes = '     order reaches p = ' // trim(cell) // ' ?'
      outcomes = outcomes // repeat(' ', max(1, 38 - len(outcomes)))
      do m = 0, max_order
         write(cell,'(i7)') m
         line = line // cell(1:7)
         if (.not. readable(m)) then
            orders   = orders   // '      -'
            spreads  = spreads  // '      -'
            outcomes = outcomes // '      -'
         else
            write(cell,'(f7.2)') fitted(m)
            orders = orders // cell(1:7)
            write(cell,'(f7.2)') spread(m)
            spreads = spreads // cell(1:7)
            if (spread(m) > settled) then
               outcomes = outcomes // '      ?'
            else if (fitted(m) >= real(expected, dp) - allowed) then
               outcomes = outcomes // '    yes'
            else
               outcomes = outcomes // '     no'
            end if
         end if
      end do
      write(*,'(a)') line
      write(*,'(a)') orders
      write(*,'(a)') spreads
      write(*,'(a)') outcomes
      write(cell,'(i0)') kept
      line = '     reaches p at ' // trim(cell) // ' of '
      write(cell,'(i0)') counted
      line = line // trim(cell) // ' derivative degrees where one power applies'
      write(*,'(a)') line
    end subroutine order_row

    !----------------------------------------------------------------!
    ! THE ORDER A SET OF GRIDS IMPLIES, AND WHETHER ONE POWER OF h
    ! EXPLAINS THEM.
    !
    ! An error that is a single power of the step satisfies
    !
    !        e = C h**q,   so   log e = log C - q log N
    !
    ! since h is proportional to 1/N. So a straight line through the
    ! points (log N, log e) has slope -q, and the least squares slope
    ! over every grid is the order. Fitting over all the grids at once
    ! rather than one pair prevents a single outlying point from
    ! determining the result.
    !
    ! WHETHER THE FIT IS VALID is a separate check, and the spread of
    ! the pairwise slopes decides it. Where the error is a single
    ! power, every consecutive pair gives the same order and the
    ! spread is near zero. Where a coefficient changes sign inside the
    ! range of grids the error passes through zero, one point falls
    ! far below the line, and the pairs disagree by a large amount - so
    ! a large spread marks the fit as invalid rather than accepting it.
    !
    !     log e                    log e
    !       |  .                     |  .
    !       |     .                  |     .
    !       |        .               |            .        <- crossing
    !       |           .            |        .
    !       +-------------- log N    +-------------- log N
    !        one power, small spread   not one power, large spread
    !----------------------------------------------------------------!

    pure subroutine order_of_grids(counts, discrepancy, fitted, spread)
      real(dp), intent(in)  :: counts(:), discrepancy(:)
      real(dp), intent(out) :: fitted, spread
      real(dp) :: x(size(counts)), y(size(counts)), pair(size(counts) - 1)
      real(dp) :: mean_x, mean_y, top, bottom
      integer  :: i, n
      n = size(counts)
      x = log(counts)
      y = log(discrepancy)
      mean_x = sum(x) / real(n, dp)
      mean_y = sum(y) / real(n, dp)
      top    = sum((x - mean_x) * (y - mean_y))
      bottom = sum((x - mean_x) ** 2)
      fitted = 0.0_dp
      if (bottom > 0.0_dp) fitted = -top / bottom
      do i = 1, n - 1
         pair(i) = -(y(i + 1) - y(i)) / (x(i + 1) - x(i))
      end do
      spread = maxval(pair) - minval(pair)
    end subroutine order_of_grids

  end subroutine demo_order_of_accuracy
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
      type(typed_field_domain) :: unknown_fields
      class(field), allocatable :: out
      type(expression) :: physics
      real(dp), allocatable :: weight(:), residual(:), acted(:), governing(:)
      real(dp) :: q(num_unknowns), t(num_instants)
      integer , allocatable :: tails(:), heads(:), head_degree(:), tail_degree(:)
      type(connectivity_graph) :: edges
      integer :: j, k
      call scheme_edges(tails, heads, tail_degree, head_degree)
      edges = connectivity_graph(num_instants, tails, heads, tail_degree, head_degree)
      call weights_of(scheme_weight(bdf_family(order)), edges, dt, weight)
      rows = derived_constraints( &
           & [(unknown(heads(j), head_degree(j)), j = 1, size(heads))], &
           & [(unknown(tails(j), tail_degree(j)), j = 1, size(tails))], &
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
      unknown_fields = typed_field_domain(unknowns % vertex_set(), num_unknowns)
      state          = unknown_fields % state(q)
      call rows % apply(unknowns, rows % bind([state]), out)
      call out % real_vector(residual)
      direction = unknown_fields % direction(q)
      call rows % partial_action(unknowns, rows % bind([state]), &
           & [variation(rows % argument(1), direction)], out)
      call out % real_vector(acted)
      physics = van_der_pol(2)
      call governing_rows(physics, q, governing)
      call show_block(governing, residual, maxval(abs(residual - acted)))
    end subroutine derived_rows
    !--------------------------------------------------------------!
    ! The rows the family itself states, rather than a pattern written
    ! out here beside it: which degree each row reads, and its history
    ! depth, are declared by the family and not by this demonstration.
    !--------------------------------------------------------------!

    subroutine scheme_edges(tails, heads, tail_degree, head_degree)
      integer, allocatable, intent(out) :: tails(:), heads(:)
      integer, allocatable, intent(out) :: tail_degree(:), head_degree(:)
      type(family) :: scheme
      integer, allocatable :: offset(:), degrees_of(:)
      integer :: d, k, e, counted, pass
      scheme = bdf_family(order)
      do pass = 1, 2
         counted = 0
         do d = 1, num_degrees - 1
            call scheme % row_pattern(d, num_degrees - 1, offset, degrees_of)
            if (size(offset) == 0) cycle
            do k = maxval(offset) + 1, num_instants
               do e = 1, size(offset)
                  counted = counted + 1
                  if (pass == 2) then
                     tails(counted)         = k - offset(e)
                     heads(counted)         = k
                     tail_degree(counted) = degrees_of(e)
                     head_degree(counted)    = d
                  end if
               end do
            end do
         end do
         if (pass == 1) allocate(tails(counted), heads(counted), &
              & tail_degree(counted), head_degree(counted))
      end do
    end subroutine scheme_edges
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
      type(typed_field_domain) :: instant_states, instant_scalars
      type(continuous_domain) :: continuous
      type(discrete_domain) :: domain
      class(field), allocatable :: out
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      continuous      = continuous_domain(physics)
      domain          = continuous % discrete(instants)
      instant_states  = domain % state_fields()
      instant_scalars = domain % design_fields()
      state           = instant_states % state(q)
      design          = instant_scalars % design(spread(nu, 1, num_instants))
      call physics % apply(instants, physics % bind([state, design]), out)
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
      type(typed_field_domain) :: point_states, point_scalars
      type(continuous_domain) :: continuous
      type(discrete_domain) :: domain
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:)
      real(dp) :: q(instants * (degree + 1)), v(instants * (degree + 1))
      real(dp) :: closed(0:degree), taken(0:degree), differenced(0:degree)
      integer :: d, k, nd
      nd      = degree + 1
      physics = van_der_pol(degree)
      call sample(degree, q0, q_top, design, instants, q)
      graph_of = stored_directed_graph(instants, tails=[integer ::], heads=[integer ::])
      continuous    = continuous_domain(physics)
      domain        = continuous % discrete(graph_of)
      point_states  = domain % state_fields()
      point_scalars = domain % design_fields()
      state         = point_states % state(q)
      nu_field      = point_scalars % design(spread(nu, 1, instants))
      do d = 0, degree
         v = 0.0_dp
         do k = 1, instants
            v((k - 1) * nd + d + 1) = 1.0_dp
         end do
         direction = point_states % direction(v)
         call physics % partial_action(graph_of, physics % bind([state, nu_field]), &
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
      call physics % apply(graph_of, physics % bind([state, nu_field]), out)
      call out % real_vector(plus)
      call state % set_real_vector(q - delta * v)
      call physics % apply(graph_of, physics % bind([state, nu_field]), out)
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
      type(typed_field_domain) :: design_fields
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:)
      real(dp) :: w(instants), q0, q_below
      w    = 1.0_dp
      q0      = q(1)
      q_below = q(nd - 1)
      design_fields = typed_field_domain(graph_of % vertex_set(), instants)
      direction     = design_fields % direction(w)
      call physics % partial_action(graph_of, physics % bind([state, nu_field]), &
           & [variation(physics % argument(2), direction)], out)
      call out % real_vector(exact)
      call nu_field % set_real_vector(spread(nu, 1, instants) + delta * w)
      call physics % apply(graph_of, physics % bind([state, nu_field]), out)
      call out % real_vector(plus)
      call nu_field % set_real_vector(spread(nu, 1, instants) - delta * w)
      call physics % apply(graph_of, physics % bind([state, nu_field]), out)
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
  end subroutine demo_constraint_rows
  subroutine demo_coupling_relation()
    implicit none
    integer , parameter :: order = 2
    integer , parameter :: num_instants = 5
    integer , parameter :: num_conditions = 2
    type(level_storage)      :: store
    type(relational_binding) :: binding
    type(set_map)            :: sets
    type(csr_relation)       :: coupling_relation
    integer, allocatable :: tails(:), heads(:), tail_degree(:), head_degree(:)
    integer, allocatable :: table(:,:)
    real(dp), allocatable :: weight(:), dt(:)
    integer :: slices(num_instants), source_carrier, target_carrier
    integer :: relation_element, coupling, block
    integer :: j, k
    ! both derived rows reach over order instants, the velocity's on
    ! the value and the acceleration's on the velocity
    tails = [((k - j, j = 0, order), k = order + 1, num_instants), &
         &   ((k - j, j = 0, order), k = order + 1, num_instants)]
    heads = [((k, j = 0, order), k = order + 1, num_instants), &
         &   ((k, j = 0, order), k = order + 1, num_instants)]
    head_degree = [((1, j = 0, order), k = order + 1, num_instants), &
         &        ((2, j = 0, order), k = order + 1, num_instants)]
    tail_degree = [((0, j = 0, order), k = order + 1, num_instants), &
         &           ((1, j = 0, order), k = order + 1, num_instants)]
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
    table(2,:) = (heads - 1) * num_conditions + head_degree
    coupling_relation = built_coupling_relation(table)
    do k = 1, num_instants
       call bind_carrier(slices(k))
    end do
    call bind_carrier(source_carrier)
    call bind_carrier(target_carrier)
    call bind_coupling_relation()
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
    function built_coupling_relation(tuples) result(r)
      integer, intent(in) :: tuples(:,:)
      type(csr_relation) :: r
      type(graph), pointer :: from, into
      from => store % node(source_carrier)
      into => store % node(target_carrier)
      r = csr_relation('scheme coupling', from, into, tuples, sets)
    end function built_coupling_relation
    subroutine bind_carrier(at)
      integer, intent(in) :: at
      type(graph), pointer :: g
      g => store % node(at)
      call binding % bind_set(g, g)
    end subroutine bind_carrier
    subroutine bind_coupling_relation()
      type(graph), pointer :: g
      g => store % node(relation_element)
      call binding % bind_relation(g, coupling_relation)
    end subroutine bind_coupling_relation
    subroutine show_tuples()
      class(relation), pointer :: r
      type(connectivity_graph) :: edges
      integer, allocatable :: fixed(:,:)
      integer :: i, e, target_index
      dt = [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp]
      edges = connectivity_graph(num_instants, tails, heads, tail_degree, head_degree)
      call weights_of(scheme_weight(bdf_family(order)), edges, dt, weight)
      r => relation_at(store % node(coupling), binding, 1)
      write(*,'(a)')    ' '
      write(*,'(a,i3)') ' tuples the relation contains  ', r % num_tuples()
      write(*,'(a)')    '   component   constraint   instant  head_degree      weight'
      call r % tuples(fixed)
      do i = 1, size(fixed, 2)
         do e = 1, size(tails)
            target_index = (heads(e) - 1) * num_conditions + head_degree(e)
            if (tails(e) == fixed(1, i) .and. target_index == fixed(2, i)) then
               write(*,'(i12,i13,i10,i12,f12.5)') fixed(1, i), fixed(2, i), &
                    & heads(e), head_degree(e), weight(e)
               exit
            end if
         end do
      end do
    end subroutine show_tuples
  end subroutine demo_coupling_relation
  subroutine demo_expansion_check()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    integer , parameter :: instants = 11
    integer , parameter :: max_order = 4
    real(dp), parameter :: duration = 2.0_dp
    real(dp), parameter :: design = 1.0_dp
    real(dp) :: delta, tau
    tau = demo_real(1, 1.0e-12_dp)
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
      type(family_container) :: owner(1)
      type(chain_block), allocatable :: chain(:)
      type(expansion)  , allocatable, target :: tower
      real(dp), allocatable :: table(:,:)
      real(dp) :: achieved
      call set_family(owner(1), scheme)
      call marched_cosine(owner, [instants], degrees, duration, design_value, chain, tower, achieved)
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
      write(*,'(a,15x,4es15.2)') '   difference        ', abs(f(1:max_order) - differenced)
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
    call newmark_on('newmark average acceleration', newmark_family(0.25_dp, 0.5_dp, 2))
    call newmark_on('taylor-newmark', taylor_newmark_family())
    call bdf_step_sensitivity([0.0_dp, 0.3_dp, 0.2_dp, 0.4_dp, 0.25_dp])
  contains
    subroutine bdf_on(label, dt)
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      integer, parameter :: last = 2 * order + 1
      type(connectivity_graph) :: edges
      real(dp), allocatable :: c(:)
      real(dp) :: h0, h1
      ! both rows are the same difference operator, the velocity's on
      ! the value and the acceleration's on the velocity, so both
      ! reach over order instants and read the degree below their own
      edges = connectivity_graph(last, &
           & [(last - k, k = 0, order), (last - k, k = 0, order)], [(last, k = 0, 2 * order + 1)], &
           & [(0, k = 0, order), (1, k = 0, order)], &
           & [(1, k = 0, order), (2, k = 0, order)])
      call weights_of(bdf_family(order), edges, dt, c)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2 on a ' // label // ' grid'
      write(*,'(a,3f10.5)') '   velocity     alpha_0..2      ', c(1:order + 1)
      write(*,'(a,5f10.5)') '   acceleration beta_0..2       ', c(order + 2:)
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
      type(connectivity_graph) :: edges
      real(dp), allocatable :: c(:)
      edges = connectivity_graph(p, [(p - k, k = 0, p - 1)], [(p, k = 0, p - 1)], &
           & [(2, k = 0, p - 1)], [(1, k = 0, p - 1)])
      call weights_of(adams_family(p), edges, dt, c)
      write(*,'(a)') ' '
      write(*,'(a)') ' adams-moulton 3 on a ' // label // ' grid'
      write(*,'(a,3f10.5)') '   quadrature   alpha_0..2      ', c
      if (label == 'uniform') then
         write(*,'(a,3f10.5)') '   tabulated    alpha           ', [5.0_dp, 8.0_dp, -1.0_dp] / 12.0_dp
      end if
    end subroutine adams_on
    subroutine dirk_on(scheme)
      type(family), intent(in) :: scheme
      type(connectivity_graph) :: edges
      real(dp), allocatable :: c(:)
      integer :: s
      s = scheme % num_stages()
      edges = connectivity_graph(2 + s, [2, 2, 3], [3, 3, 3], [2, 2, 2], [1, 1, 1])
      call weights_of(scheme, edges, [(0.5_dp, k = 1, 2 + s)], c)
      write(*,'(a)') ' '
      write(*,'(a)') ' crouzeix two-stage, stage 2 from stages 1, 1, 2'
      write(*,'(a,3f10.5)') '   a_21, a_21, a_22                ', c
      edges = connectivity_graph(2 + s, [2, 3], [2 + s, 2 + s], [2, 2], [2, 2])
      call weights_of(scheme, edges, [(0.5_dp, k = 1, 2 + s)], c)
      write(*,'(a,2f10.5)') '   b_1, b_2 into the end instant   ', c
      write(*,'(a,3f10.5)') '   tableau gamma, 1 - 2 gamma, b   ', &
           & (3.0_dp + sqrt(3.0_dp)) / 6.0_dp, 1.0_dp - (3.0_dp + sqrt(3.0_dp)) / 3.0_dp, 0.5_dp
    end subroutine dirk_on
    subroutine newmark_on(label, scheme)
      character(len=*), intent(in) :: label
      type(family)    , intent(in) :: scheme
      type(connectivity_graph) :: edges
      real(dp), allocatable :: c(:)
      edges = connectivity_graph(3, [2, 2, 2, 3], [3, 3, 3, 3], &
           & [0, 1, 2, 2], [0, 0, 0, 0])
      call weights_of(scheme, edges, [0.0_dp, 0.5_dp, 0.5_dp], c)
      write(*,'(a)') ' '
      write(*,'(a)') ' ' // label // ', value row'
      write(*,'(a,4f10.5)') '   q, qdot, qddot behind/ahead  ', c
      edges = connectivity_graph(3, [2, 2, 3], [3, 3, 3], [1, 2, 2], [1, 1, 1])
      call weights_of(scheme, edges, [0.0_dp, 0.5_dp, 0.5_dp], c)
      write(*,'(a)') ' ' // label // ', velocity row'
      write(*,'(a,3f10.5)') '   qdot, qddot behind/ahead     ', c
    end subroutine newmark_on
    subroutine bdf_step_sensitivity(dt)
      real(dp), intent(in) :: dt(:)
      integer , parameter :: last = 2 * order + 1
      real(dp), parameter :: delta = 1.0e-6_dp
      type(stored_directed_graph) :: coupling
      type(connectivity_graph) :: edges
      type(stored_field), allocatable :: inputs(:)
      type(stored_field) :: direction
      type(typed_field_domain) :: coupling_fields
      type(family) :: scheme
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)
      scheme = bdf_family(order)
      edges = connectivity_graph(last, [(last - k, k = 0, order)], [(last, k = 0, order)], &
           & [(0, k = 0, order)], [(1, k = 0, order)])
      call coupling_inputs(edges, dt, inputs)
      coupling = edges % stored_directed_graph
      allocate(v(last), source=0.0_dp)
      v(last) = 1.0_dp
      coupling_fields = typed_field_domain(coupling % vertex_set(), last)
      direction       = coupling_fields % direction(v)
      call scheme % partial_action(coupling, scheme % bind(inputs), [variation(scheme % argument(1), direction)], out)
      call out % real_vector(exact)
      call inputs(1) % set_real_vector(dt + delta * v)
      call scheme % apply(coupling, scheme % bind(inputs), out)
      call out % real_vector(plus)
      call inputs(1) % set_real_vector(dt - delta * v)
      call scheme % apply(coupling, scheme % bind(inputs), out)
      call out % real_vector(minus)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2, partial of alpha_0..2 in the last step, non-uniform grid'
      write(*,'(a,3f12.6)') '   partial_action, degree one    ', exact
      write(*,'(a,3f12.6)') '   central difference            ', (plus - minus) / (2.0_dp * delta)
      call bdf_second_partials(scheme, coupling, inputs(1), inputs(2), inputs(3), dt, v)
    end subroutine bdf_step_sensitivity
    subroutine bdf_second_partials(scheme, coupling, steps, degrees, conditions, dt, v)
      type(family)           , intent(in)    :: scheme
      type(stored_directed_graph), intent(in)    :: coupling
      type(stored_field)         , intent(inout) :: steps
      type(stored_field)         , intent(in)    :: degrees, conditions
      real(dp)                   , intent(in)    :: dt(:), v(:)
      real(dp), parameter :: delta = 1.0e-4_dp
      type(stored_field) :: along_v, along_w
      type(typed_field_domain) :: coupling_fields
      class(field), allocatable :: out
      real(dp), allocatable :: plus(:), at(:), minus(:), w(:)
      real(dp), allocatable :: second(:), mixed_partial(:)
      w = 0.0_dp * v
      w(size(dt) - 1) = 1.0_dp
      coupling_fields = typed_field_domain(coupling % vertex_set(), size(dt))
      along_v = coupling_fields % direction(v)
      along_w = coupling_fields % direction(w)
      call steps % set_real_vector(dt)
      call scheme % partial_action(coupling, scheme % bind([steps, degrees, conditions]), &
           & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_v)], out)
      call out % real_vector(second)
      call scheme % partial_action(coupling, scheme % bind([steps, degrees, conditions]), &
           & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_w)], out)
      call out % real_vector(mixed_partial)
      call steps % set_real_vector(dt + delta * v)
      call scheme % apply(coupling, scheme % bind([steps, degrees, conditions]), out)
      call out % real_vector(plus)
      call steps % set_real_vector(dt)
      call scheme % apply(coupling, scheme % bind([steps, degrees, conditions]), out)
      call out % real_vector(at)
      call steps % set_real_vector(dt - delta * v)
      call scheme % apply(coupling, scheme % bind([steps, degrees, conditions]), out)
      call out % real_vector(minus)
      write(*,'(a)') ' '
      write(*,'(a)') ' bdf 2, second partials of alpha_0..2, non-uniform grid'
      write(*,'(a,3f12.6)') '   partial_action, (last, last)  ', second
      write(*,'(a,3f12.6)') '   second central difference     ', (plus - 2.0_dp * at + minus) / delta**2
      call steps % set_real_vector(dt + delta * w)
      call scheme % partial_action(coupling, scheme % bind([steps, degrees, conditions]), &
           & [variation(scheme % argument(1), along_v)], out)
      call out % real_vector(plus)
      call steps % set_real_vector(dt - delta * w)
      call scheme % partial_action(coupling, scheme % bind([steps, degrees, conditions]), &
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
       write(*,'(a,i0,a)') ' FAIL : ', failures, ' identities exceeded the rounding bound'
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
      call agrees(n, 'exp(a+b) = exp(a) exp(b)', exp(a + b), exp(a) * exp(b), exp(a + b), failures)
      call agrees(n, 'sin^2 + cos^2 = 1', sin(a) * sin(a) + cos(a) * cos(a), one, sin(a) * sin(a), failures)
      call agrees(n, 'sin(a+b) addition', sin(a + b), sin(a) * cos(b) + cos(a) * sin(b), sin(a) * cos(b), failures)
      call agrees(n, 'log(exp(a)) = a', log(exp(a)), a, exp(a), failures)
      call agrees(n, 'exp(log(a)) = a', exp(log(a)), a, a, failures)
      call agrees(n, 'sqrt(a) sqrt(a) = a', sqrt(a) * sqrt(a), a, a, failures)
      call agrees(n, 'a**0.5 = sqrt(a)', a ** 0.5_dp, sqrt(a), sqrt(a), failures)
      call agrees(n, 'a**3.0 = a a a', a ** 3.0_dp, a * a * a, a * a * a, failures)
      call agrees(n, 'a**(-1.0) = 1/a', a ** (-1.0_dp), one / a, one / a, failures)
      a = seeded(0.4_dp, v)
      closed = exp(0.4_dp) * product(v)
      call agrees_scalar(n, 'full partial of exp', mixed_partial(exp(a)), closed, failures)
      closed = sin(0.4_dp + real(n, dp) * acos(-1.0_dp) / 2.0_dp) * product(v)
      call agrees_scalar(n, 'full partial of sin', mixed_partial(sin(a)), closed, failures)
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
    subroutine agrees(n, label, obtained, reference, operand, failures)
      integer               , intent(in)    :: n
      character(len=*)      , intent(in)    :: label
      type(derivative_terms), intent(in)    :: obtained, reference, operand
      integer               , intent(inout) :: failures
      real(dp) :: maximum_departure, scale, rounding_bound
      integer  :: m
      maximum_departure = 0.0_dp
      scale = 0.0_dp
      do m = 0, 2**n - 1
         maximum_departure = max(maximum_departure, abs(coefficient(obtained, m) - coefficient(reference, m)))
         scale = max(scale, abs(coefficient(operand, m)))
      end do
      rounding_bound = real(3**n, dp) * epsilon(1.0_dp) * scale
      call reported(n, label, maximum_departure, rounding_bound, failures)
    end subroutine agrees
    subroutine agrees_scalar(n, label, obtained, reference, failures)
      integer         , intent(in)    :: n
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: obtained, reference
      integer         , intent(inout) :: failures
      call reported(n, label, abs(obtained - reference), &
           & real(3**n, dp) * epsilon(1.0_dp) * abs(reference), failures)
    end subroutine agrees_scalar
    subroutine reported(n, label, maximum_departure, rounding_bound, failures)
      integer         , intent(in)    :: n
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: maximum_departure, rounding_bound
      integer         , intent(inout) :: failures
      write(*,'(a,i0,a,a28,a,es9.2,a,es9.2)') '   n = ', n, '  ', label, &
           & '  difference ', maximum_departure, '  bound ', rounding_bound
      if (maximum_departure > rounding_bound) failures = failures + 1
    end subroutine reported
  end subroutine demo_function_identities
  subroutine demo_grid_design_check()
    implicit none
    integer , parameter :: state_degree = 2, degrees = state_degree + 1
    integer , parameter :: instants = 21, checked(3) = [1, 7, 20]
    real(dp), parameter :: duration = 4.0_dp, design = 0.8_dp
    type(family_container)     :: schemes(2)
    type(expression)       :: functionals(2)
    type(chain_block) , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    integer, allocatable :: versions(:)
    real(dp), allocatable :: p(:), dt(:), t(:), v(:,:), f(:,:), tangent(:,:), adjoint(:,:)
    real(dp), allocatable :: plus(:,:), minus(:,:), q0(:), table(:,:), entries(:,:,:)
    real(dp), allocatable :: below(:,:), above(:,:), by_class(:)
    real(dp) :: tau, delta, achieved, maximum_departure
    integer :: k, j, i, pass_kind, order, max_order, nd
    tau       = demo_real(1, 1.0e-12_dp)
    max_order = nint(demo_real(2, 3.0_dp))
    call set_stopping(tau, relative, by_rate, 100)
    delta = tau ** (1.0_dp / 3.0_dp)
    schemes = [stored_family(bdf_family(3)), stored_family(adams_family(3))]
    functionals(1) = van_der_pol_energy(state_degree)
    functionals(2) = van_der_pol_dissipation(state_degree)
    p  = [(1.0_dp + 0.5_dp * sin(real(k, dp)), k = 1, instants - 1)]
    q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design)
    call marched(p, design, f)
    call tower % step_partials(v)
    call directions(chain, tower, functionals, degrees, tangent, adjoint)
    pass_kind   = pass_of(size(tangent, 2), size(tangent, 1), 1)
    write(*,'(a,es9.2,a,es9.2,a,es9.2)') ' relative tolerance', tau, '   difference step', delta, &
         & '   expected agreement tau^(2/3)', tau ** (2.0_dp / 3.0_dp)
    write(*,'(a,i0,a,i0,a,a)') ' designs ', size(tangent, 2), '   functionals ', size(tangent, 1), &
         & '   pass_of selects the ', trim(merge('forward', 'reverse', pass_kind == forward_pass))
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
    maximum_departure = 0.0_dp
    do k = 1, size(checked)
       j = checked(k)
       call marched(p + delta * unit(j), design, plus)
       call marched(p - delta * unit(j), design, minus)
       do i = 1, 2
          maximum_departure = max(maximum_departure, abs((plus(0, i) - minus(0, i)) / (2.0_dp * delta) - adjoint(i, 1 + j)) &
               & / max(1.0_dp, abs(adjoint(i, 1 + j))))
       end do
    end do
    write(*,'(a,es10.2)') ' differenced in three weights against the pass, largest ', maximum_departure
    call marched(p, design + delta, plus)
    call marched(p, design - delta, minus)
    write(*,'(a,es10.2)') ' differenced in the parameter against the pass, largest ', &
         & maxval(abs((plus(0, :) - minus(0, :)) / (2.0_dp * delta) - adjoint(:, 1)) / &
         &        max(1.0_dp, abs(adjoint(:, 1))))
    nd = size(tangent, 2)
    do order = 2, max_order
       call marched(p, design, f, order)
       call chain_versions(chain, tower, functionals, degrees, versions)
       call chain_derivative(chain, tower, versions, functionals, degrees, order, reverse_pass, &
            & table, entries=entries)
       write(*,'(a)') ' '
       write(*,'(a,i0,a,i0,a,i0,a,i0)') ' derivatives of order ', order, ' by the reverse pass: ', &
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
            & ' against the entries, largest by parameter count from ', &
            & by_class(order:0:-1) / max(1.0_dp, maxval(abs(table)))
       deallocate(by_class)
    end do
  contains
    subroutine marched(weights, nu, f, order)
      real(dp), intent(in) :: weights(:), nu
      real(dp), allocatable, intent(out) :: f(:,:)
      integer , intent(in), optional :: order
      type(imbalance) :: final_imbalance
      integer :: m
      m = 1
      if (present(order)) m = order
      call march_chain(schemes, [11, 10], van_der_pol(state_degree), degrees, &
           & designed_grid(duration), nu, q0, chain, tower, dt, t, achieved, grid_design=weights, &
           & final_imbalance=final_imbalance, startup=4)
      if (.not. final_imbalance % converged) error stop 'grid_design_check: the march converged'
      call chain_expansion(chain, tower, functionals, degrees, m, f)
    end subroutine marched
    subroutine differenced_table(weights, nu, order, t)
      real(dp), intent(in) :: weights(:), nu
      integer , intent(in) :: order
      real(dp), allocatable, intent(out) :: t(:,:)
      real(dp), allocatable :: f(:,:)
      call marched(weights, nu, f)
      call chain_versions(chain, tower, functionals, degrees, versions)
      call chain_derivative(chain, tower, versions, functionals, degrees, order, reverse_pass, t)
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
  !===================================================================!
  ! THE LEIBNIZ EXPANSION OF THE LAGRANGIAN. With L = F - lambda . R
  ! and R(Q(nu); nu) = 0 along the march, the derivative of order n of
  ! the functional in the design is the explicit derivative of the
  ! Lagrangian of order n - 1,
  !
  !    d^n F     d  [  (n-1)    n-1  (n-1)     (k)    (n-1-k) ]
  !    ----- = ---- [ F      -  sum  (   ) lambda  . R        ]
  !    dnu^n   dnu  [           k=0  ( k )                    ]
  !
  ! the outer derivative holding every derivative of the state fixed,
  ! and the costate of order k solving the transposed block with the
  ! lower costates on its right-hand side. The reverse pass forms
  ! this as one product per row over n directions, n - 1 for the
  ! state's derivatives and one for the design alone: a subset of
  ! size k of the n - 1 represents the binomial count once, and the
  ! coefficient of a smaller subset containing the explicit direction is
  ! the entry of a lower order, so one reverse pass at order N yields
  ! every order to N. The costates come from the transposed solves and
  ! the state's derivatives from the forward ones as before: the one
  ! product is the assembly, not the solves.
  !
  ! Bounds, each stated: the terms add to the table to the rounding of
  ! their count; the table of every order from one reverse pass agrees
  ! with the forward pass to what two solves at the march's relative
  ! tolerance tau admit, tau^(2/3).
  !===================================================================!
  subroutine demo_lagrangian_expansion()
    integer , parameter :: state_degree = 2, degrees = state_degree + 1, instants = 21
    real(dp), parameter :: duration = 4.0_dp, design = 0.8_dp
    type(family_container) :: schemes(1)
    type(expression)    :: functionals(2)
    type(chain_block), allocatable :: chain(:)
    type(expansion)  , allocatable, target :: tower
    type(imbalance)  :: final_imbalance
    integer , allocatable :: versions(:)
    real(dp), allocatable :: q0(:), dt(:), t(:), f(:,:), table(:,:), by_order(:,:,:), terms(:,:,:,:)
    real(dp) :: tau, achieved, agreement, departure, scale, rounding_bound
    integer :: max_order, order, n, k, i, b, rows, failures
    failures  = 0
    tau       = demo_real(1, 1.0e-12_dp)
    max_order = nint(demo_real(2, 4.0_dp))
    call set_stopping(tau, relative, by_rate, 100)
    agreement      = tau ** (2.0_dp / 3.0_dp)
    schemes(1)     = stored_family(bdf_family(3))
    functionals(1) = van_der_pol_energy(state_degree)
    functionals(2) = van_der_pol_dissipation(state_degree)
    q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design)
    call march_chain(schemes, [instants], van_der_pol(state_degree), degrees, uniform_grid(duration), &
         & design, q0, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, startup=4)
    if (.not. final_imbalance % converged) error stop 'lagrangian_expansion: the march converged'
    rows = 0
    do b = 1, size(chain)
       rows = rows + chain(b) % rows % num_unknowns()
    end do
    call chain_expansion(chain, tower, functionals, degrees, max_order, f)
    call chain_versions(chain, tower, functionals, degrees, versions)
    call chain_derivative(chain, tower, versions, functionals, degrees, max_order, reverse_pass, table, &
         & by_order=by_order)
    write(*,'(a)') ' bdf 3 over a crouzeix start, van der pol, the physics'' parameter the one design'
    write(*,'(a,es9.2,a,es9.2)') ' relative tolerance', tau, &
         & '   agreement of two solves at it, tau^(2/3)', agreement
    write(*,'(a)') ' '
    write(*,'(a,i0,a)') ' every order from one reverse pass at order ', max_order, ', against the forward pass'
    write(*,'(a)') '   order  functional       forward pass         reverse pass    departure'
    do n = 0, max_order
       do i = 1, size(functionals)
          departure = abs(by_order(i, 1, n) - f(n, i)) / max(1.0_dp, abs(f(n, i)))
          write(*,'(i8,i12,2es20.10,es13.2)') n, i, f(n, i), by_order(i, 1, n), departure
          if (departure > agreement) failures = failures + 1
       end do
    end do
    do order = 1, max_order
       call chain_derivative(chain, tower, versions, functionals, degrees, order, reverse_pass, table, &
            & leibniz=terms)
       n = order - 1
       write(*,'(a)') ' '
       write(*,'(a,i0,a,i0,a)') ' order ', order, ': the explicit derivative of the lagrangian of order ', &
            & n, ', term by term'
       do i = 1, size(functionals)
          scale = sum(abs(terms(i, 1, 1, :)))
          rounding_bound = real(rows, dp) * real(2**order, dp) * epsilon(1.0_dp) * scale
          write(*,'(a,i0)') '   functional ', i
          do k = 0, n
             write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,es20.10)') '     C(', n, ',', k, ') = ', choose(n, k), &
                  & '   lambda^(', k, ') . d/dnu R^(', n - k, ')   ', terms(i, 1, 1, k)
          end do
          write(*,'(a,i0,a,es20.10)') '     d/dnu F^(', n, ')                            ', terms(i, 1, 1, n + 1)
          departure = abs(sum(terms(i, 1, 1, :)) - table(i, 1))
          write(*,'(a,es20.10,a,es9.2,a,es9.2)') '     the terms summed, against the table  ', table(i, 1), &
               & '   difference ', departure, '   bound ', rounding_bound
          if (departure > rounding_bound) failures = failures + 1
          departure = abs(table(i, 1) - f(order, i)) / max(1.0_dp, abs(f(order, i)))
          write(*,'(a,es20.10,a,es9.2)') '     the forward pass                       ', f(order, i), &
               & '   departure ', departure
          if (departure > agreement) failures = failures + 1
       end do
    end do
    write(*,'(a)') ' '
    if (failures > 0) then
       write(*,'(a,i0,a)') ' FAIL : ', failures, ' checks exceeded their bound'
       error stop
    end if
    write(*,'(a)') ' PASS : the reverse pass is the leibniz expansion of the lagrangian, every order from one pass'
  end subroutine demo_lagrangian_expansion
  !===================================================================!
  ! THE TAYLOR STATE MARCH. The forward pass solves, block by block,
  ! the tower of the state's derivatives in the design to the requested
  ! order: the zeroth order by Newton in the march, every higher order
  ! by one linear solve against the block's own factorisation. A
  ! block's tower is read by the blocks whose given instants lie in it
  ! and by no later block, so it is released once the last such block
  ! has been solved. The live storage is then the history depth in
  ! blocks, independent of the horizon; the reverse pass reads every
  ! tower during the reverse traversal and retains them all.
  !
  ! One horizon is marched in one block, two, four and eight. The
  ! tables agree to what two solves at the tolerance admit, tau^(2/3);
  ! the maximum live tower storage equals the largest sum of two
  ! consecutive towers, an exact count, so its bound is zero.
  !===================================================================!
  subroutine demo_taylor_state()
    integer , parameter :: state_degree = 2, degrees = state_degree + 1, solved = 40
    integer , parameter :: splits(4) = [1, 2, 4, 8]
    real(dp), parameter :: duration = 8.0_dp, design = 0.8_dp
    type(family_container), allocatable :: schemes(:)
    type(expression)    :: functionals(2)
    type(chain_block), allocatable :: chain(:)
    type(expansion)  , allocatable, target :: tower
    type(imbalance)  :: final_imbalance
    integer , allocatable :: versions(:), added(:), sizes(:)
    real(dp), allocatable :: q0(:), dt(:), t(:), table(:,:), whole(:,:), by_order(:,:,:)
    real(dp) :: tau, achieved, agreement, departure
    integer :: max_order, storage(2), expected, nb, b, k, failures
    failures  = 0
    tau       = demo_real(1, 1.0e-12_dp)
    max_order = nint(demo_real(2, 4.0_dp))
    call set_stopping(tau, relative, by_rate, 100)
    agreement      = tau ** (2.0_dp / 3.0_dp)
    functionals(1) = van_der_pol_energy(state_degree)
    functionals(2) = van_der_pol_dissipation(state_degree)
    q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design)
    write(*,'(a,i0,a,i0,a)') ' bdf 3 over a crouzeix start, van der pol, ', solved, &
         & ' instants, derivatives to order ', max_order, ' in the physics'' parameter'
    write(*,'(a,es9.2,a,es9.2)') ' relative tolerance', tau, &
         & '   agreement of two solves at it, tau^(2/3)', agreement
    write(*,'(a)') ' '
    write(*,'(a)') '   blocks   tower numbers live at most   of all tower numbers   expected' // &
         & '   departure of the tables from one block'
    do k = 1, size(splits)
       nb = splits(k)
       allocate(schemes(nb), added(nb))
       do b = 1, nb
          schemes(b) = stored_family(bdf_family(3))
          added(b)   = solved / nb
       end do
       call march_chain(schemes, added, van_der_pol(state_degree), degrees, uniform_grid(duration), &
            & design, q0, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, startup=4)
       if (.not. final_imbalance % converged) error stop 'taylor_state: the march converged'
       call chain_versions(chain, tower, functionals, degrees, versions)
       call chain_derivative(chain, tower, versions, functionals, degrees, max_order, forward_pass, &
            & table, designs=1, by_order=by_order, tower_storage=storage)
       allocate(sizes(size(chain)))
       do b = 1, size(chain)
          sizes(b) = chain(b) % rows % num_unknowns() * max_order
       end do
       expected = sizes(1)
       do b = 2, size(chain)
          expected = max(expected, sizes(b - 1) + sizes(b))
       end do
       if (k == 1) whole = by_order(:, 1, :)
       departure = maxval(abs(by_order(:, 1, :) - whole) / max(1.0_dp, abs(whole)))
       write(*,'(i9,i22,i22,i11,es30.2)') nb, storage(1), storage(2), expected, departure
       if (storage(1) /= expected) failures = failures + 1
       if (departure > agreement) failures = failures + 1
       deallocate(schemes, added, sizes)
    end do
    call chain_derivative(chain, tower, versions, functionals, degrees, max_order, reverse_pass, &
         & table, designs=1, tower_storage=storage)
    write(*,'(a)') ' '
    write(*,'(a,i0,a,i0)') ' the reverse pass over the last chain retains every tower: live at most ', &
         & storage(1), ' of ', storage(2)
    if (storage(1) /= storage(2)) failures = failures + 1
    write(*,'(a)') ' '
    call pipelined()
    write(*,'(a)') ' '
    if (failures > 0) then
       write(*,'(a,i0,a)') ' FAIL : ', failures, ' checks exceeded their bound'
       error stop
    end if
    write(*,'(a)') ' PASS : the taylor state march stores one block''s history depth of towers, and the tables are unchanged'
  contains
    !----------------------------------------------------------------!
    ! THE PIPELINED MARCH: dirk 3, one instant per block, so the
    ! horizon advances one instant per block at the family's order. Each block's tower is
    ! solved with the block and released past its last reader, states
    ! included, inside the march itself. Doubling the horizon must
    ! leave the live storage unchanged - an integer, bound zero -
    ! while the totals double; the tables agree with the
    ! whole-horizon march of the same grid to tau^(2/3).
    !----------------------------------------------------------------!
    subroutine pipelined()
      integer, parameter :: horizons(3) = [20, 40, 80]
      type(family_container), allocatable :: steps(:)
      real(dp), allocatable :: fp(:,:), fr(:,:)
      integer :: ts(2), ss(2), previous_ts, previous_ss, n, j
      previous_ts = 0
      previous_ss = 0
      write(*,'(a,i0)') ' dirk 3, one instant per block: the pipelined march, derivatives to order ', max_order
      write(*,'(a)') '   instants    towers live   of all    states live   of all   departure from the whole march'
      do j = 1, size(horizons)
         n = horizons(j)
         allocate(steps(n - 1))
         do b = 1, n - 1
            steps(b) = stored_family(crouzeix_three_stage())
         end do
         call march_chain(steps, [2, (1, b = 2, n - 1)], van_der_pol(state_degree), degrees, &
              & uniform_grid(duration), design, q0, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, &
              & functionals=functionals, derivative_order=max_order, f=fp, &
              & tower_storage=ts, state_storage=ss)
         if (.not. final_imbalance % converged) error stop 'taylor_state: the pipelined march converged'
         deallocate(steps)
         allocate(steps(1))
         steps(1) = stored_family(crouzeix_three_stage())
         call march_chain(steps, [n], van_der_pol(state_degree), degrees, uniform_grid(duration), &
              & design, q0, chain, tower, dt, t, achieved, final_imbalance=final_imbalance)
         if (.not. final_imbalance % converged) error stop 'taylor_state: the whole march converged'
         call chain_expansion(chain, tower, functionals, degrees, max_order, fr)
         departure = maxval(abs(fp - fr) / max(1.0_dp, abs(fr)))
         write(*,'(i11,i15,i9,i15,i9,es22.2)') n, ts(1), ts(2), ss(1), ss(2), departure
         if (departure > agreement) failures = failures + 1
         if (j > 1 .and. (ts(1) /= previous_ts .or. ss(1) /= previous_ss)) failures = failures + 1
         previous_ts = ts(1)
         previous_ss = ss(1)
         deallocate(steps)
      end do
    end subroutine pipelined
  end subroutine demo_taylor_state
  !===================================================================!
  ! WHERE ONE BLOCK'S VALUES ARE STORED INSIDE ANOTHER'S STATE. The
  ! assembler computes every transfer offset before a block is built,
  ! from the tower alone; the block that is then built computes its
  ! instant offsets independently. The two must agree at every instant
  ! of every block, or a block would receive the wrong values and solve
  ! a different problem without detection. The departure counts offsets
  ! that differ, so the bound is zero.
  !===================================================================!

  subroutine demo_transfer_offsets()
    implicit none
    integer , parameter :: state_degree = 2
    integer , parameter :: degrees = state_degree + 1
    real(dp), parameter :: duration = 2.0_dp
    ! a block adds more instants than its family's history depth,
    ! so the counts exceed every history depth used
    call checked('bdf 2 alone      ', [container_named('bdf', 2)], [8])
    call checked('bdf 2 then bdf 1 ', [container_named('bdf', 2), container_named('bdf', 1)], [8, 8])
    call checked('adams 2 then bdf2', [container_named('adams', 2), container_named('bdf', 2)], [8, 8])
    call checked('dirk, staged     ', [container_named('dirk', 3)], [8])
  contains
    subroutine checked(title, schemes, added)
      character(len=*)   , intent(in) :: title
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: first(:), last(:), settled(:)
      real(dp), allocatable :: dt(:), t(:), fixed(:)
      real(dp) :: achieved
      integer :: b, k, n, discrepancy, counted
      call horizon_bounds(schemes, added, degrees - 1, first, last)
      n = last(size(added))
      call cosine_partition(schemes(1) % scheme, degrees, duration, n, fixed, dt, t)
      call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
           & uniform_grid(duration), 0.0_dp, fixed, chain, tower, dt, t, achieved)
      discrepancy   = 0
      counted = 0
      do b = 1, size(chain)
         settled = instants_at_of(tower, b, chain(b) % scheme, van_der_pol(state_degree))
         if (size(settled) /= size(chain(b) % instants_at)) then
            discrepancy   = discrepancy + 1
            counted = counted + 1
            cycle
         end if
         do k = 1, size(settled)
            counted = counted + 1
            if (settled(k) /= chain(b) % instants_at(k)) discrepancy = discrepancy + 1
         end do
      end do
      write(*,'(a,a,i4,a,i0,a,i0)') '   ', title, counted, &
           & ' offsets computed before building, differing ', discrepancy, ', bound ', 0
    end subroutine checked
  end subroutine demo_transfer_offsets

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
      type(family_container), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: versions(:)
      real(dp), allocatable :: a(:,:)
      real(dp) :: achieved, duration, design
      duration = 3.0_dp
      design   = 1.0_dp
      allocate(schemes(1))
      call set_family(schemes(1), scheme)
      call marched_cosine(schemes, [instants], degrees, duration, design, chain, tower, achieved)
      call chain_versions(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, versions)
      call dense_jacobian(chain, design, a)
      call reported(label, a)
    end subroutine shape_of
    subroutine reported(label, a)
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: a(:,:)
      integer  :: n, i, j, filled, below, above
      real(dp) :: largest_entry, least, on_diagonal, row_most
      n       = size(a, 1)
      largest_entry = maxval(abs(a))
      least   = 1.0e-12_dp * largest_entry
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
           & largest_entry, on_diagonal, row_most
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
    write(*,'(a)') ' the tower, and what each level stores'
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
      type(branch) :: members
      integer :: k
      name = ' '
      if (sets % labelled(g)) name = sets % label_of(g)
      write(*,'(a,a,a,a)') repeat('   ', depth + 1), name, &
           & '   [' // status_name(values % status_of(g)) // ']', extent_of(g)
      if (level_is_leaf(g)) return
      members = level_members(g)
      do k = 1, sequence_num_elements(members)
         call show(sequence_element(members, k), depth + 1)
      end do
    end subroutine show
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
    function extent_of(g) result(label)
      type(graph), intent(in) :: g
      character(len=:), allocatable :: label
      character(len=3) :: n
      label = ''
      if (.not. sets % describes(g)) return
      write(n,'(i3)') sets % num_members_of(g)
      label = '   extent' // n
    end function extent_of
    recursive subroutine determine(g)
      type(graph), intent(in) :: g
      type(branch) :: members
      integer :: k
      if (level_is_leaf(g)) then
         if (values % status_of(g) == VALUE_UNKNOWN) then
            call values % mark_known(g, spread(1.0_dp, 1, num_freedoms))
         end if
         return
      end if
      members = level_members(g)
      do k = 1, sequence_num_elements(members)
         call determine(sequence_element(members, k))
      end do
    end subroutine determine
    recursive integer function not_yet_known(g) result(n)
      type(graph), intent(in) :: g
      type(branch) :: members
      integer :: k
      n = 0
      if (level_is_leaf(g)) then
         if (values % status_of(g) == VALUE_UNKNOWN) n = 1
         return
      end if
      members = level_members(g)
      do k = 1, sequence_num_elements(members)
         n = n + not_yet_known(sequence_element(members, k))
      end do
    end function not_yet_known
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
    write(*,'(a)') ' a multistep block: slices contain components'
    root => store % node(multistep)
    call show(root, 1)
    write(*,'(a)') ' '
    write(*,'(a)') ' a multistage block: slices contain stages'
    root => store % node(multistage)
    call show(root, 1)
    write(*,'(a)')    ' '
    write(*,'(a,i0)') ' nodes owned by the storage   ', store % num_nodes()
    call foreign_member_is_refused()
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
      type(branch) :: members
      integer :: k
      if (level_is_leaf(g)) then
         write(*,'(a,a)') repeat('   ', depth), 'leaf'
         return
      end if
      write(*,'(a,a,i0,a,l1,a,l1)') repeat('   ', depth), 'members ', &
           & level_num_members(g), '   couples ', level_couples(g), &
           & '   consistent ', level_consistent(g)
      members = level_members(g)
      do k = 1, sequence_num_elements(members)
         call show(sequence_element(members, k), depth + 1)
      end do
    end subroutine show
    subroutine foreign_member_is_refused()
      type(graph), pointer :: subject_node, other
      integer, allocatable :: own(:)
      integer :: non_member, foreign, subject
      own        = members_of(2, one_leaf)
      non_member = one_leaf()
      foreign    = store % assemble([own(1), non_member], 0)
      subject = coupled(own)
      subject_node => store % node(subject)
      write(*,'(a)')    ' '
      write(*,'(a,l1)') ' coupling beginning with its own members is consistent ', &
           & level_consistent(subject_node)
      other => store % node(foreign)
      subject_node % branch(2) = known_branch(other)
      write(*,'(a,i0)') ' the other coupling has the same count                 ', &
           & level_num_members(subject_node)
      write(*,'(a,l1)') ' and is refused by identity                            ', &
           & level_consistent(subject_node)
    end subroutine foreign_member_is_refused
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
      type(family_container) :: owner(1)
      integer, allocatable :: at(:)
      real(dp), allocatable :: dt(:), fixed(:)
      call cosine_partition(scheme, degrees, duration, n, fixed, dt, t)
      call set_family(owner(1), scheme)
      call tower % build(van_der_pol(max_state_degree), owner, [n], uniform_grid(duration), 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(max_state_degree), fixed, rows, at)
      call solved(rows, design_value, q, achieved)
    end subroutine march
    pure real(dp) function maximum_departure(q, t) result(e)
      real(dp), intent(in) :: q(:), t(:)
      integer :: k
      e = 0.0_dp
      do k = 1, size(t)
         e = max(e, abs(q(unknown(k, 0)) - cosine(0, t(k))))
      end do
    end function maximum_departure
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
         e(level) = maximum_departure(q, t)
      end do
      write(*,'(a)')          ' '
      write(*,'(a)')          ' ' // title // ' on the harmonic oscillator'
      write(*,'(a,4i11)')     '   steps                      ', [(10 * 2 ** (level - 1), level = 1, 4)]
      write(*,'(a,4es11.3)')  '   maximum error              ', e
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
    call splitting_is_invariant()
    call across_a_change_of_scheme()
    call sensitivity_across_the_block_boundary()
  contains
    subroutine marched(schemes, added, q, achieved, design_value)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp), allocatable, intent(out) :: q(:)
      real(dp)           , intent(out) :: achieved
      real(dp), intent(in), optional :: design_value
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      real(dp) :: design
      integer :: k, n
      design = 0.0_dp
      if (present(design_value)) design = design_value
      n = sum(added)
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      allocate(q(n * degrees))
      do k = 1, n
         q((k - 1) * degrees + 1:k * degrees) = instant_components(chain, k)
      end do
    end subroutine marched
    subroutine splitting_is_invariant()
      type(family_container) :: whole(1), split(2)
      real(dp), allocatable :: q_whole(:), q_split(:)
      real(dp) :: achieved_whole, achieved_split
      whole = [stored_family(bdf_family(2))]
      split = [stored_family(bdf_family(2)), stored_family(bdf_family(2))]
      call marched(whole, [40], q_whole, achieved_whole)
      call marched(split, [20, 20], q_split, achieved_split)
      write(*,'(a)')        ' bdf 2 over forty instants, in one block and in two'
      write(*,'(a,i0)')     '   unknowns, whole             ', size(q_whole)
      write(*,'(a,i0)')     '   unknowns, split             ', size(q_split)
      write(*,'(a,es11.2)') '   maximum difference           ', maxval(abs(q_whole - q_split))
      write(*,'(a,es11.2)') '   residual, whole             ', achieved_whole
      write(*,'(a,es11.2)') '   residual, split             ', achieved_split
    end subroutine splitting_is_invariant
    subroutine across_a_change_of_scheme()
      type(family_container) :: schemes(2)
      real(dp), allocatable :: q(:), dt(:), t(:)
      real(dp) :: e(3), achieved
      integer :: level, added, k, n
      schemes = [stored_family(bdf_family(2)), stored_family(adams_family(3))]
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
      write(*,'(a,3es11.3)')  '   maximum error              ', e
      write(*,'(a,22x,2f11.3)') '   ratio                    ', e(1:2) / e(2:3)
      write(*,'(a,es11.2)')   '   residual                   ', achieved
    end subroutine across_a_change_of_scheme
    subroutine sensitivity_across_the_block_boundary()
      real(dp), parameter :: delta = 1.0e-6_dp
      real(dp), parameter :: design = 1.0_dp
      integer , parameter :: added(2) = [10, 10]
      type(family_container) :: schemes(2)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      type(expression)       :: energy(1)
      real(dp), allocatable :: q(:), table(:,:), other(:,:)
      real(dp) :: f, tangent, adjoint, differenced, achieved
      integer :: n
      schemes = [stored_family(bdf_family(2)), stored_family(adams_family(3))]
      n = sum(added)
      call marched(schemes, added, q, achieved, design)
      f = chained_energy(schemes, added, design)
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      energy(1) = van_der_pol_energy(state_degree)
      call directions(chain, tower, energy, degrees, table, other)
      tangent     = first_of(table)
      adjoint     = first_of(other)
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
    end subroutine sensitivity_across_the_block_boundary
    !--------------------------------------------------------------!
    ! THE FUNCTIONAL THE CHAIN ITSELF INTEGRATES. A multistep family
    ! weights the instants within its own stencil's history, so a sum
    ! of the instants at their step sizes is a different rule and a
    ! different number - one order lower. A tangent is the derivative
    ! of the functional it was taken from, so the finite difference
    ! below is taken of this functional and of no other.
    !--------------------------------------------------------------!

    real(dp) function chained_energy(schemes, added, design) result(f)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      real(dp)           , intent(in) :: design
      type(chain_block), allocatable :: chain(:)
      type(expansion)  , allocatable, target :: tower
      type(expression) :: energy(1)
      real(dp), allocatable :: table(:,:)
      real(dp) :: achieved
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      energy(1) = van_der_pol_energy(state_degree)
      call chain_expansion(chain, tower, energy, degrees, 0, table)
      f = table(0, 1)
    end function chained_energy

    real(dp) function differenced_energy(schemes, added, n, design, delta) result(d)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:), n
      real(dp)           , intent(in) :: design, delta
      associate (u1 => n); end associate
      d = (chained_energy(schemes, added, design + delta) - &
         & chained_energy(schemes, added, design - delta)) / (2.0_dp * delta)
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
    real(dp) function maximum_error(scheme, n) result(e)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: n
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_container) :: owner(1)
      real(dp), allocatable :: dt(:), t(:), q(:)
      integer , allocatable :: at(:)
      real(dp) :: achieved
      integer :: k, d
      call partition(duration, n, dt, t)
      call set_family(owner(1), scheme)
      call tower % build(van_der_pol(state_degree), owner, [n], uniform_grid(duration), 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), &
           & [(cosine(d, t(1)), d = 0, degrees - 1)], rows, at)
      call solved(rows, 0.0_dp, q, achieved)
      e = 0.0_dp
      do k = 1, n
         e = max(e, abs(q(at(k) + 1) - cosine(0, t(k))))
      end do
    end function maximum_error
    subroutine order_of(title, scheme)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      real(dp) :: e(3)
      integer :: level
      do level = 1, 3
         e(level) = maximum_error(scheme, 5 * 2 ** (level - 1) + 1)
      end do
      write(*,'(a)')            ' '
      write(*,'(a)')            ' ' // title // ' on the harmonic oscillator'
      write(*,'(a,i0,a)')       '   stages                     ', scheme % num_stages(), ''
      write(*,'(a,3i11)')       '   steps                      ', [(5 * 2 ** (level - 1), level = 1, 3)]
      write(*,'(a,3es11.3)')    '   maximum error              ', e
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
    type(family_container) :: owner(1)
    integer, allocatable :: at(:)
    type(family)            :: scheme
    character(len=32) :: what, given
    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: dt(:), t(:), fixed(:)
    integer :: instants, n, m, h, band, i, j, e
    what     = 'block'
    call demo_argument(1, given)
    if (len_trim(given) > 0) what = given
    instants = nint(demo_real(2, 41.0_dp))
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
       allocate(fixed(h * degrees), source=0.0_dp)
       call set_family(owner(1), scheme)
       call tower % build(van_der_pol(degrees - 1), owner, [instants], uniform_grid(3.0_dp), 0, 0.0_dp)
       call block_from(tower, 1, scheme, van_der_pol(degrees - 1), fixed, rows, at)
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
    seed  = nint(demo_real(1, 7.0_dp))
    cases = nint(demo_real(2, 2.0_dp))
    failures = 0
    skipped  = 0
    write(*,'(a)') '  case  scheme          split          directions     passes'
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
         & failures, ' failed, ', skipped, ' skipped'
    if (failures > 0) error stop 'randomized_checks: an invariant did not hold'
  contains
    logical function verbose()
      character(len=8) :: argument
      call demo_argument(3, argument)
      verbose = trim(argument) == 'verbose'
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
      type(family_container), intent(in)  :: schemes(:)
      integer            , intent(in)  :: added(:), degrees
      real(dp)           , intent(in)  :: duration, design
      real(dp)           , intent(out) :: tangent, adjoint
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      real(dp), allocatable :: table(:,:), other(:,:)
      real(dp) :: achieved
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      call directions(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, table, other)
      tangent = first_of(table)
      adjoint = first_of(other)
    end subroutine directions_of
    subroutine one_case(from, index, failures, skipped)
      integer, intent(in)    :: from, index
      integer, intent(inout) :: failures, skipped
      type(family_container), allocatable :: whole(:), split(:)
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
      call whole_and_split(whole, split, instants, half, degrees, duration, design, &
           & f_whole, f_split, achieved)
      if (achieved > 1.0e-6_dp) then
         skipped = skipped + 1
         write(*,'(i6,2x,a16,a,es9.2)') index, label, &
              & '   a march did not converge, skipped: ', achieved
         return
      end if
      call directions_of(whole, [instants], degrees, duration, design, tangent, adjoint)
      if (verbose()) write(*,'(a,2es14.6)') '        f whole and split ', f_whole(0), f_split(0)
      call check_agreement(index, label, f_whole, f_split, tangent, adjoint, failures)
    end subroutine one_case
    subroutine whole_and_split(whole, split, instants, half, degrees, duration, design, &
         & f_whole, f_split, achieved)
      type(family_container), intent(in)  :: whole(:), split(:)
      integer            , intent(in)  :: instants, half, degrees
      real(dp)           , intent(in)  :: duration, design
      real(dp), allocatable, intent(out) :: f_whole(:), f_split(:)
      real(dp)           , intent(out) :: achieved
      real(dp) :: one, two
      call expanded(whole, [instants], degrees, duration, design, f_whole, one)
      call expanded(split, [half, instants - half], degrees, duration, design, &
           & f_split, two)
      achieved = max(one, two)
    end subroutine whole_and_split
    subroutine mixed_case(from, index, failures, skipped)
      integer, intent(in)    :: from, index
      integer, intent(inout) :: failures, skipped
      real(dp), parameter :: delta = 1.0e-4_dp
      type(family_container), allocatable :: schemes(:)
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
         write(*,'(i6,2x,a28,a)') index, label, '  march did not converge, skipped'
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
      type(family_container), allocatable, intent(inout) :: schemes(:)
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
      integer :: depth
      depth = scheme % history_depth(degrees - 1)
      half  = max(instants / 2, depth + 1)
      if (instants - half <= depth) then
         instants = 2 * (depth + 1)
         half     = instants / 2
      end if
    end subroutine halved
    ! kind 3 is the crouzeix two-stage tableau, the dirk of order three
    subroutine fill(fixed, kind, order)
      type(family_container), intent(out) :: fixed
      integer            , intent(in)  :: kind, order
      fixed = container_named(trim(family_names(kind)), merge(order, 3, kind < 3))
    end subroutine fill
    function named(kind, order) result(name)
      integer, intent(in) :: kind, order
      character(len=16) :: name
      character(len=1) :: digit
      write(digit,'(i1)') order
      name = trim(family_names(kind)) // digit
      if (kind == 3) name = 'crouzeix2'
    end function named
    subroutine expanded(schemes, added, degrees, duration, design, f, achieved)
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:), degrees
      real(dp)           , intent(in) :: duration, design
      real(dp), allocatable, intent(out) :: f(:)
      real(dp), allocatable :: table(:,:)
      real(dp)           , intent(out) :: achieved
      type(chain_block), allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      call marched_cosine(schemes, added, degrees, duration, design, chain, tower, achieved)
      call chain_expansion(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, max_order, table)
      allocate(f(lbound(table, 1):ubound(table, 1)))
      f = table(:, 1)
    end subroutine expanded
    subroutine check_agreement(index, label, f_whole, f_split, tangent, adjoint, failures)
      integer         , intent(in)    :: index
      character(len=*), intent(in)    :: label
      real(dp)        , intent(in)    :: f_whole(0:), f_split(0:), tangent, adjoint
      integer         , intent(inout) :: failures
      real(dp) :: split_gap, direction_gap, pass_gap, scale
      scale         = max(1.0_dp, maxval(abs(f_whole)))
      split_gap     = maxval(abs(f_whole - f_split)) / scale
      direction_gap = abs(tangent - adjoint) / max(1.0_dp, abs(tangent))
      pass_gap     = abs(tangent - f_whole(1)) / max(1.0_dp, abs(tangent))
      write(*,'(i6,2x,a16,3es15.2)') index, label, split_gap, direction_gap, pass_gap
      if (split_gap > 1.0e-6_dp) failures = failures + 1
      if (direction_gap > 1.0e-8_dp) failures = failures + 1
      if (pass_gap > 1.0e-6_dp) failures = failures + 1
    end subroutine check_agreement
  end subroutine demo_randomized_checks
  !===================================================================!
  ! THE ARCS DETERMINE THE ORDER, NOT THE NUMBERING. A bipartite digraph
  ! is built whose blocks are chained 2 -> 4 -> 1 -> 3, so the order
  ! the arcs imply differs from the order of the labels. The
  ! projection onto the blocks must recover that chain, previous and
  ! next must be transposes of one another, and every arc must be
  ! found from both of its ends. Each departure counts an exact
  ! disagreement, so the bound is zero and no tolerance is chosen.
  !===================================================================!

  subroutine demo_read_write_graph()

    integer, parameter :: blocks = 4
    integer, parameter :: data_of = 3

    ! block 2 writes datum 1; block 4 reads it and writes datum 2;
    ! block 1 reads that and writes datum 3; block 3 reads datum 3
    integer, parameter :: from_part(7)   = [BLOCKS_PART, DATA_PART, BLOCKS_PART, &
         & DATA_PART, BLOCKS_PART, DATA_PART, BLOCKS_PART]
    integer, parameter :: from_vertex(7) = [2, 1, 4, 2, 1, 3, 3]
    integer, parameter :: to_part(7)     = [DATA_PART, BLOCKS_PART, DATA_PART, &
         & BLOCKS_PART, DATA_PART, BLOCKS_PART, DATA_PART]
    integer, parameter :: to_vertex(7)   = [1, 4, 2, 1, 3, 3, 3]

    integer, parameter :: chained(4) = [2, 4, 1, 3]

    type(bipartite_digraph) :: b
    type(stored_directed_graph) :: among
    integer, allocatable :: order(:), forth(:), back(:), reads(:)
    integer :: u, w, k, e, discrepancy_order, discrepancy_transpose, discrepancy_ends

    b = bipartite_digraph(blocks, data_of, from_part, from_vertex, to_part, to_vertex)

    write(*,'(a)')      ' '
    write(*,'(a)')      ' a bipartite digraph over blocks and the data between them'
    write(*,'(a,i0,a,i0,a,i0)') '   blocks ', b % order_of_part(BLOCKS_PART), &
         & '   data ', b % order_of_part(DATA_PART), &
         & '   arcs ', b % size_of_digraph()

    ! the projection recovers the chain the arcs imply
    among = b % projection(BLOCKS_PART)
    order = among % loop(forward)
    discrepancy_order = count(order /= chained)
    write(*,'(a)')      ' '
    write(*,'(a,4i4)')  '   the arcs chain the blocks   ', chained
    write(*,'(a,4i4)')  '   the projection orders them  ', order
    write(*,'(a,i0,a,i0)') '   positions differing ', discrepancy_order, ', bound ', 0

    ! previous and next are one relation and its transpose
    discrepancy_transpose = 0
    do u = 1, blocks
       call b % next(BLOCKS_PART, u, forth)
       do k = 1, size(forth)
          call b % previous(BLOCKS_PART, forth(k), back)
          if (.not. any(back == u)) discrepancy_transpose = discrepancy_transpose + 1
       end do
    end do
    write(*,'(a,i0,a,i0)') '   next without a previous ', discrepancy_transpose, ', bound ', 0

    ! every arc is found from its head
    discrepancy_ends = 0
    do e = 1, b % size_of_digraph()
       w = to_vertex(e)
       if (to_part(e) == BLOCKS_PART) then
          call b % in_neighbourhood(BLOCKS_PART, w, reads)
       else
          call b % in_neighbourhood(DATA_PART, w, reads)
       end if
       if (.not. any(reads == from_vertex(e))) discrepancy_ends = discrepancy_ends + 1
    end do
    write(*,'(a,i0,a,i0)') '   arcs not found from the head ', discrepancy_ends, ', bound ', 0

    ! a block and the one two steps along it share no datum
    write(*,'(a)')      ' '
    write(*,'(a,l1,a)') '   blocks 2 and 4 share a datum   ', &
         & b % share_a_neighbour(BLOCKS_PART, 2, 4), '   (they are chained)'
    write(*,'(a,l1,a)') '   blocks 2 and 3 share a datum   ', &
         & b % share_a_neighbour(BLOCKS_PART, 2, 3), '   (no arc joins them)'

  end subroutine demo_read_write_graph

  !===================================================================!
  ! STATE DEPENDENCE IN THE TRANSPOSE. A chain of three blocks, and
  ! the transpose the derivative is taken over.
  !
  !   forward   (1) --> [1] --> (2) --> [2] --> (3) --> [3]
  !                      |               |               |
  !                      v               v               v
  !                     (4) <-- [5] <-- (5) <-- [6] <-- (6)
  !   the transpose, releasing the states in decreasing b
  !
  ! The three downward arcs state that the reverse pass depends on
  ! each block's state. Without them the sweep would end at block 3:
  ! state 1 would have its final dependence at step 2 and state 3
  ! would have no dependent vertex at all, so a driver would release
  ! both while the derivative still required them.
  !
  ! Three departures are counted, and each has bound zero. A STATE
  ! WITHOUT A DEPENDENT is one the graph would let a driver release
  ! immediately after it was written. A STATE RELEASED ELSEWHERE is
  ! one whose final dependent is not its own transpose block, or which
  ! is never released at all: block nb + b is the final dependent of
  ! state b exactly because every forward dependence on state b sends
  ! a transposed arc into it. A STEP OUT OF REVERSE is one where the
  ! transpose does not retrace the forward sweep backwards - which is
  ! the purpose of exchanging the ends of every forward arc, and which
  ! a transpose that retained its ends would not satisfy while still
  ! depending on every state.
  !
  ! The fourth count is the driver's result. No rule is bound at a
  ! transposed vertex, and a step that computes nothing still
  ! releases whatever has final dependence there - so a traversal that
  ! placed a value at every state must leave none of them stored. STATES
  ! STILL STORED counts the ones a traversal did not release.
  !===================================================================!

  subroutine demo_transposed_dependencies()

    implicit none
    integer, parameter :: state_degree = 2
    integer, parameter :: degrees = state_degree + 1

    write(*,'(a)') ' '
    write(*,'(a)') ' the state dependencies of a reverse pass, and the release step of each'
    write(*,'(a)') ' '
    call checked('bdf 1 alone       ', [container_named('bdf', 1)], [8])
    call checked('bdf 2 then bdf 1  ', [container_named('bdf', 2), container_named('bdf', 1)], [8, 8])
    call checked('adams 2 then bdf 2', [container_named('adams', 2), container_named('bdf', 2)], [8, 8])
    call checked('three of bdf 2    ', [container_named('bdf', 2), container_named('bdf', 2), &
         & container_named('bdf', 2)], [8, 8, 8])

  contains

    subroutine checked(title, schemes, added)
      character(len=*)   , intent(in) :: title
      type(family_container), intent(in) :: schemes(:)
      integer            , intent(in) :: added(:)
      type(bipartite_digraph) :: incidence
      type(driver)     :: schedule
      type(temporal_minimizer) :: executor
      type(expression) :: immaterial
      type(rule_graph) :: rules
      type(data_graph) :: values, remaining
      type(pairing)    :: pairs
      type(stored_directed_graph) :: one_point, bare
      type(stored_field) :: datum
      type(typed_field_domain) :: point_fields
      integer, allocatable :: first(:), last(:), order(:), releasable(:)
      real(dp), allocatable :: no_rhs(:), no_solution(:)
      real(dp) :: driver_achieved
      integer :: nb, b, k, no_dependent, elsewhere, released, out_of_reverse, still_stored

      nb = size(added)
      call horizon_bounds(schemes, added, degrees - 1, first, last)
      incidence = chain_incidence(schemes, added, degrees, 1, first, last)

      ! THE LIFETIMES ARE THE GRAPH'S AND THE ORDER'S, and no rule
      ! enters either result, so the rule passed here is never
      ! applied and nothing is marched.
      bare = stored_directed_graph(nb, tails=[integer ::], heads=[integer ::])
      schedule = driver(immaterial, incidence, forward)
      call executor % state(schedule, bare, bare % vertex_set(), 0)
      order    = executor % visits()

      no_dependent = 0
      elsewhere = 0
      do b = 1, nb
         if (executor % last_dependent_of(b) < 1) then
            no_dependent = no_dependent + 1
            cycle
         end if
         if (order(executor % last_dependent_of(b)) /= nb + b) elsewhere = elsewhere + 1
      end do

      ! the transpose retraces the forward sweep backwards, step for
      ! step, so the last nb steps are the first nb read in reverse
      out_of_reverse = 0
      do k = 1, nb
         if (order(size(order) - k + 1) /= nb + order(k)) out_of_reverse = out_of_reverse + 1
      end do

      ! every state is released by the end of the transpose
      releasable = executor % released_after(size(order))
      released   = count(releasable <= nb)

      ! A TRAVERSAL RELEASES WHAT IT VISITS. Every state is given a
      ! value and no rule is bound anywhere, so nothing is computed and
      ! the traversal's only effect is to apply the lifetimes.
      allocate(rules % at(2 * nb), values % at(2 * nb))
      one_point = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
      point_fields = typed_field_domain(one_point % vertex_set(), 1)
      datum        = point_fields % state([0.0_dp])
      do b = 1, nb
         allocate(values % at(b) % datum, source=datum)
      end do
      call executor % pair_with(rules % pair(values))
      allocate(no_rhs(0), no_solution(0))
      call executor % solve(no_rhs, no_solution, driver_achieved)
      pairs     = executor % pairing_of()
      remaining  = pairs % stored_data()
      still_stored = count([(remaining % at(b) % written(), b = 1, nb)])

      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '   ', title, &
           & '  blocks ', nb, '  steps ', size(order), &
           & '   no dependent ', no_dependent, ', bound ', 0, &
           & ';  released elsewhere ', elsewhere + (nb - released), ', bound ', 0, &
           & ';  steps out of reverse ', out_of_reverse, ', bound ', 0, &
           & ';  still stored ', still_stored, ', bound ', 0
    end subroutine checked

  end subroutine demo_transposed_dependencies

  subroutine demo_scheme_weights()
    implicit none
    integer :: k
    call bdf_rows(2, 'uniform',     [0.0_dp, (0.5_dp, k = 2, 5)])
    call bdf_rows(2, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp])
    call adams_row(3, 'uniform',     [0.0_dp, 0.5_dp, 0.5_dp])
    call adams_row(3, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp])
    call newmark_rows('newmark average acceleration', newmark_family(0.25_dp, 0.5_dp, 2), 2)
    call newmark_rows('newmark Fox-Goodwin', newmark_family(1.0_dp / 12.0_dp, 0.5_dp, 3), 2)
    call newmark_rows('taylor-newmark', taylor_newmark_family(), 2)
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
    subroutine row_fields(scheme, nv, tails, head, tail_degree, head_degree, dt, &
         & tau, alpha, w)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: nv, tails(:), head, tail_degree(:), head_degree(:)
      real(dp)     , intent(in) :: dt(:)
      real(dp), allocatable, intent(out) :: tau(:), alpha(:), w(:)
      type(stored_directed_graph) :: coupling
      type(connectivity_graph) :: edges
      type(stored_field), allocatable :: inputs(:)
      class(field), allocatable :: out
      type(scheme_weight) :: weights
      integer :: e
      edges = connectivity_graph(nv, tails, [(head, e = 1, size(tails))], tail_degree, head_degree)
      call coupling_inputs(edges, dt, inputs)
      coupling = edges % stored_directed_graph
      tau = [(step_power(dt(head), tail_degree(e) - head_degree(e)), e = 1, size(tails))]
      call scheme % apply(coupling, scheme % bind(inputs), out)
      call out % real_vector(alpha)
      weights = scheme_weight(scheme)
      call weights % apply(coupling, weights % bind(inputs), out)
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
    pure real(dp) function row_residual(w, tails, head, tail_degree, head_degree, t, m) &
         & result(r)
      real(dp), intent(in) :: w(:), t(:)
      integer , intent(in) :: tails(:), head, tail_degree(:), head_degree(:), m
      integer :: e
      r = -power_derivative(m, head_degree(1), t(head))
      do e = 1, size(w)
         r = r + w(e) * power_derivative(m, tail_degree(e), t(tails(e)))
      end do
    end function row_residual
    subroutine one_row(title, scheme, nv, tails, head, tail_degree, head_degree, dt, top, exact_top)
      character(len=*), intent(in) :: title
      class(family)   , intent(in) :: scheme
      integer         , intent(in) :: nv, tails(:), head, tail_degree(:), head_degree(:), top
      real(dp)        , intent(in) :: dt(:)
      integer         , intent(in), optional :: exact_top
      real(dp), allocatable :: tau(:), alpha(:), w(:)
      real(dp) :: t(nv), residual(0:top)
      integer :: m
      call row_fields(scheme, nv, tails, head, tail_degree, head_degree, dt, tau, alpha, w)
      t = instants(dt)
      write(*,'(a)') ' '
      write(*,'(a)')        ' ' // title
      write(*,'(a,9f11.5)') '   tau                        ', tau
      write(*,'(a,9f11.5)') '   alpha                      ', alpha
      write(*,'(a,9f11.5)') '   weight                     ', w
      do m = 0, top
         residual(m) = row_residual(w, tails, head, tail_degree, head_degree, t, m)
      end do
      write(*,'(a,9i11)')     '   on t**m, m =              ', [(m, m = 0, top)]
      write(*,'(a,9es11.2)')  '   residual                  ', residual
      if (present(exact_top)) call require_exact(title, residual, exact_top)
    end subroutine one_row
    subroutine require_exact(title, residual, exact_top)
      character(len=*), intent(in) :: title
      real(dp)        , intent(in) :: residual(0:)
      integer         , intent(in) :: exact_top
      real(dp) :: floor
      if (exact_top > ubound(residual, 1)) then
         error stop 'gti_demos: an exactness check lies inside the reported degree range'
      end if
      floor = real(128 * max(1, exact_top + 1), dp) * epsilon(1.0_dp)
      if (maxval(abs(residual(0:exact_top))) > floor) then
         write(*,'(a)') ' '
         write(*,'(a)') ' row failed: ' // title
         write(*,'(a,es12.4)') '   floor       ', floor
         write(*,'(a,9es12.4)') '   residual    ', residual(0:exact_top)
         error stop 'gti_demos: a scheme row reproduces its polynomial class'
      end if
    end subroutine require_exact
    subroutine bdf_rows(p, label, dt)
      integer         , intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      integer :: last
      last = 2 * p + 1
      call one_row('bdf ' // digit(p) // ' velocity row, ' // label // ' grid', &
           & bdf_family(p), last, [(last - k, k = 0, p)], last, &
           & [(0, k = 0, p)], [(1, k = 0, p)], dt, p + 2)
      ! the acceleration row differences the velocity, so it reads
      ! degree one over the same p instants the velocity row reads
      ! degree zero over
      call one_row('bdf ' // digit(p) // ' acceleration row, ' // label // ' grid', &
           & bdf_family(p), last, [(last - k, k = 0, p)], last, &
           & [(1, k = 0, p)], [(2, k = 0, p)], dt, p + 2)
    end subroutine bdf_rows
    subroutine adams_row(p, label, dt)
      integer         , intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp)        , intent(in) :: dt(:)
      call one_row('adams-moulton ' // digit(p) // ' velocity row, ' // label // ' grid', &
           & adams_family(p), p, [p - 1, (p - k, k = 0, p - 1)], p, &
           & [1, (2, k = 0, p - 1)], [(1, k = 0, p)], dt, p + 2)
    end subroutine adams_row
    subroutine newmark_rows(label, scheme, exact_top)
      character(len=*), intent(in) :: label
      type(family)    , intent(in) :: scheme
      integer         , intent(in) :: exact_top
      real(dp), parameter :: dt(3) = [0.0_dp, 0.5_dp, 0.5_dp]
      call one_row(label // ' value row, uniform grid', &
           & scheme, 3, [2, 2, 2, 3], 3, [0, 1, 2, 2], &
           & [0, 0, 0, 0], dt, exact_top, exact_top=exact_top)
      call one_row(label // ' velocity row, uniform grid', &
           & scheme, 3, [2, 2, 3], 3, [1, 2, 2], &
           & [1, 1, 1], dt, exact_top, exact_top=exact_top)
    end subroutine newmark_rows
    subroutine weight_partials(p, dt)
      integer , intent(in) :: p
      real(dp), intent(in) :: dt(:)
      real(dp), parameter :: delta = 1.0e-6_dp
      type(stored_directed_graph) :: coupling
      type(connectivity_graph) :: edges
      type(stored_field), allocatable :: inputs(:)
      type(stored_field) :: direction
      type(typed_field_domain) :: coupling_fields
      type(scheme_weight) :: weights
      class(field), allocatable :: out
      real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)
      integer , allocatable :: tails(:)
      integer :: last, e, j
      last  = 2 * p + 1
      tails = [(last - j, j = 0, p)]
      edges = connectivity_graph(last, tails, [(last, e = 1, size(tails))], &
           & [(0, e = 1, size(tails))], [(1, e = 1, size(tails))])
      call coupling_inputs(edges, dt, inputs)
      coupling = edges % stored_directed_graph
      allocate(v(last), source=0.0_dp)
      v(last) = 1.0_dp
      coupling_fields = typed_field_domain(coupling % vertex_set(), last)
      direction       = coupling_fields % direction(v)
      weights = scheme_weight(bdf_family(p))
      call weights % partial_action(coupling, weights % bind(inputs), [variation(weights % argument(1), direction)], out)
      call out % real_vector(exact)
      call inputs(1) % set_real_vector(dt + delta * v)
      call weights % apply(coupling, weights % bind(inputs), out)
      call out % real_vector(plus)
      call inputs(1) % set_real_vector(dt - delta * v)
      call weights % apply(coupling, weights % bind(inputs), out)
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
      type(family_container) :: owner(1)
      integer, allocatable :: at(:)
      type(stored_directed_graph) :: unknowns, instants
      type(stored_field) :: state, design_field
      type(typed_field_domain) :: energy_states, instant_scalars
      type(continuous_domain) :: energy_domain
      type(discrete_domain) :: energy_points
      real(dp), allocatable :: dt(:), t(:), fixed(:)
      real(dp) :: achieved
      type(expression) :: energy
      call cosine_partition(scheme, degrees, duration, num_instants, fixed, dt, t)
      call set_family(owner(1), scheme)
      call tower % build(van_der_pol(state_degree), owner, [num_instants], uniform_grid(duration), &
           & 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), fixed, rows, at)
      call solved(rows, design_value, q, achieved)
      unknowns = unknowns_graph(num_instants, degrees)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      energy   = van_der_pol_energy(state_degree)
      energy_domain  = continuous_domain(energy)
      energy_points  = energy_domain % discrete(instants)
      energy_states  = energy_points % functional_state_fields()
      instant_scalars = energy_points % design_fields()
      state           = energy_states % functional_state(q)
      design_field    = instant_scalars % design(spread(design_value, 1, num_instants))
      f = functional_of(energy, instants, [state, design_field], dt)
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
      integer :: version
      type(block_residual) :: rows
      type(expansion) :: tower
      type(family_container) :: owner(1)
      integer, allocatable :: at(:)
      type(stored_directed_graph) :: unknowns, instants
      type(stored_field) :: state, design_field, energy_state, energy_design
      type(typed_field_domain) :: unknown_fields, unknown_designs, energy_states, instant_scalars
      type(continuous_domain) :: energy_domain
      type(discrete_domain) :: energy_points
      type(expression) :: energy
      real(dp), allocatable :: dt(:), t(:), fixed(:)
      call cosine_partition(scheme, degrees, duration, num_instants, fixed, dt, t)
      call set_family(owner(1), scheme)
      call tower % build(van_der_pol(state_degree), owner, [num_instants], uniform_grid(duration), &
           & 0, 0.0_dp)
      call block_from(tower, 1, scheme, van_der_pol(state_degree), fixed, rows, at)
      unknowns = unknowns_graph(num_instants, degrees)
      instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
      unknown_fields = typed_field_domain(unknowns % vertex_set(), size(q))
      unknown_designs = typed_field_domain(unknowns % vertex_set(), num_instants)
      state        = unknown_fields % state(q)
      design_field = unknown_designs % design(spread(design, 1, num_instants))
      energy       = van_der_pol_energy(state_degree)
      energy_domain  = continuous_domain(energy)
      energy_points  = energy_domain % discrete(instants)
      energy_states  = energy_points % functional_state_fields()
      instant_scalars = energy_points % design_fields()
      energy_state    = energy_states % functional_state(q)
      energy_design   = instant_scalars % design(spread(design, 1, num_instants))
      call functional_gradient(energy, instants, &
           & [energy_state, energy_design], dt, num_instants, degrees, unknowns % vertex_set(), g)
      call sweep_design_partial(rows, unknowns, [state, design_field], num_instants, &
           & unknowns % vertex_set(), rate)
      version    = next_version()
      tangent = by_tangent(rows, unknowns, [state, design_field], g, rate, 0.0_dp, version)
      adjoint = by_adjoint(rows, unknowns, [state, design_field], g, rate, 0.0_dp, version)
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
      type(family_container), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: versions(:)
      type(expression)       :: energy(1)
      real(dp), allocatable :: table(:,:)
      real(dp) :: achieved, duration, design, marched, formed, solved_in, tangent
      integer  :: n
      duration = 3.0_dp
      design   = 1.0_dp
      allocate(schemes(1))
      call set_family(schemes(1), bdf_family(2))
      marched = clock()
      call marched_cosine(schemes, [instants], degrees, duration, design, chain, tower, achieved)
      marched = clock() - marched
      energy(1) = van_der_pol_energy(degrees - 1)
      formed = clock()
      call chain_versions(chain, tower, energy, degrees, versions)
      formed = clock() - formed
      n = chain(1) % rows % num_unknowns()
      solved_in = clock()
      call chain_derivative(chain, tower, versions, energy, degrees, 1, forward_pass, table)
      tangent   = first_of(table)
      solved_in = clock() - solved_in
      write(*,'(i10,3f11.3,2es15.3)') n, marched, formed, solved_in, &
           & achieved, 1.0e-12_dp
    end subroutine cost_at
  end subroutine demo_solve_cost
  subroutine demo_tolerance_form()
    implicit none
    write(*,'(a)') ' '
    write(*,'(a)') '  the velocity row, and the d-th row that repeats it'
    write(*,'(a)') '  scheme    d    sum|a|    sum|c(d)|      sum|a|    the same'
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
    write(*,'(a)') '  weight the family declares for it'
    write(*,'(a)') ' '
    write(*,'(a)') '  scheme    d   instants       dt      kappa(A)   kappa(DA)     ratio'
    call conditioned(2, 3, 21)
    call conditioned(2, 3, 41)
    call conditioned(2, 3, 61)
    call conditioned(2, 3, 81)
    call conditioned(3, 3, 41)
    call conditioned(2, 4, 41)
  contains
    real(dp) function row_sum(scheme, order, head_degree) result(total)
      class(family), intent(in) :: scheme
      integer      , intent(in) :: order, head_degree
      type(connectivity_graph) :: edges
      real(dp), allocatable :: c(:)
      integer :: depth, last, k
      depth = order
      last  = depth + 1
      edges = connectivity_graph(last, [(last - k, k = 0, depth)], [(last, k = 0, depth)], &
           & [(head_degree - 1, k = 0, depth)], [(head_degree, k = 0, depth)])
      call weights_of(scheme, edges, [(1.0_dp, k = 1, last)], c)
      total = sum(abs(c))
    end function row_sum
    subroutine composed(order, head_degree)
      integer, intent(in) :: order, head_degree
      real(dp) :: velocity, derived, powered
      character(len=8) :: named
      ! EVERY DERIVED ROW IS THE SAME OPERATOR, the velocity's on the
      ! value and each higher one on the degree below it, so every row
      ! sums its coefficients to the same total. Composing the
      ! operator on the value instead would raise that total to the
      ! power of the degree, and double the stencil's history depth.
      velocity = row_sum(bdf_family(order), order, 1)
      derived  = row_sum(bdf_family(order), order, head_degree)
      powered  = velocity
      write(named,'(a,i0)') 'bdf ', order
      write(*,'(a,a,i5,3f12.4,a)') '  ', named, head_degree, velocity, derived, powered, &
           & merge('   yes', '    no', abs(derived - powered) <= 1.0e-10_dp * powered)
    end subroutine composed
    subroutine against(order, degrees, instants)
      integer, intent(in) :: order, degrees, instants
      type(family_container), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: versions(:)
      type(family) :: scheme
      integer , allocatable :: added(:)
      real(dp), allocatable :: fixed(:), dt(:), t(:), a(:,:)
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
      fixed = [((0.0_dp, d = 0, degrees - 1), k = 1, &
           &   scheme % history_depth(degrees - 1))]
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, fixed, chain, tower, dt, t, achieved)
      call chain_versions(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, versions)
      predicted = 1.0_dp + row_sum(scheme, order, top) / dt(size(dt))
      call dense_jacobian(chain, design, a)
      assembled = largest_row(a)
      write(named,'(a,i0)') 'bdf ', order
      write(*,'(a,a,i5,i10,f10.5,2es15.5,f11.3,a)') '  ', named, top, instants, &
           & dt(size(dt)), predicted, assembled, &
           & 100.0_dp * (1.0_dp - abs(predicted - assembled) / assembled), ' %'
    end subroutine against
    subroutine conditioned(order, degrees, instants)
      integer, intent(in) :: order, degrees, instants
      type(family_container), allocatable :: schemes(:)
      type(chain_block) , allocatable :: chain(:)
      type(expansion), allocatable, target :: tower
      integer, allocatable :: versions(:)
      type(family) :: scheme
      integer , allocatable :: added(:)
      real(dp), allocatable :: fixed(:), dt(:), t(:), b(:,:), a(:,:)
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
      fixed = [((0.0_dp, d = 0, degrees - 1), k = 1, &
           &   scheme % history_depth(degrees - 1))]
      call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
           & uniform_grid(duration), design, fixed, chain, tower, dt, t, achieved)
      call chain_versions(chain, tower, [van_der_pol_energy(degrees - 1)], degrees, versions)
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
  use operation_family      , only : bdf_family
  use operation_family      , only : adams_family
  use operation_grid        , only : uniform_grid, random_grid, designed_grid, fixed_grid
  use operation_expression  , only : expression, stated_over
  use gti_physics           , only : van_der_pol, van_der_pol_energy
  use operation_grid        , only : grid
  use gti_march             , only : set_stopping, imbalance, set_space_coupling, set_time_coupling, weight_of, precision_needed
  use gti_adaptive          , only : adaptive_partition
  use operation_family      , only : crouzeix_three_stage
  use operation_stencil     , only : stencil
  use operation_domain      , only : continuous_domain
  use gti_space             , only : spatial_domain, spatial_mesh, geometry_of, coarse_cells
  use gti_field             , only : spatial_discretization_stencil_of, initial_field, against_the_laplacian, &
       & against_the_mode, export_instant
  use util_precision        , only : precision_named
  use iso_fortran_env       , only : real128
  use gti_expansion         , only : family_container, expansion
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & expansion_substitutions, chain_versions, num_designs_of, &
       & instant_components, chain_derivative, asymmetry, sink_costates, &
       & goal_oriented_partition
  use gti_sweeps            , only : set_linear_solver, set_jacobian, set_storage, set_multigrid, &
       & set_rows, set_elimination, &
       & set_coarse_nodes, set_linear_budget, set_newton_order
  use gti_sweeps            , only : pass_of, forward_pass, reverse_pass
  use operation_minimization, only : relative, absolute, by_count, by_rate
  use gti_driver            , only : settings, chosen_grid, steps_of, family_named, clock, &
       & functional_named
  use gti_configuration     , only : configuration, read_configuration, override, show, &
       & lists, refuse_unknown, words_of
  use util_tally            , only : tally_open, tally_close, tally_order, &
       & tally_enter, tally_leave, tally_amount, tally_event_of, &
       & tally_num_levels, tally_level_name, tally_event_name, elapsed_time
  use gti_configuration     , only : at_expansion, at_horizon, hierarchy_levels
  use gti_demos           , only : demo_requested, run_demo
  implicit none
  type(configuration) :: cfg
  type(spatial_domain)   , allocatable :: space
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
  call set_jacobian(cfg % jacobian)
  call set_rows(cfg % rows)
  call set_elimination(cfg % elimination)
  call set_storage(cfg % storage)
  call set_multigrid(cfg % multigrid)
  call set_space_coupling(cfg % space)
  call set_time_coupling(cfg % time)
  call field_context(cfg)
  call chosen_functionals(cfg)
  call table(cfg)
contains
  integer function widest_depth(cfg) result(widest)
    type(configuration), intent(in) :: cfg
    character(len=16) :: every(5)
    class(family), allocatable :: scheme
    logical :: staged, passes_check
    integer :: i, order, depth
    every  = [character(len=16) :: 'bdf', 'adams', 'dirk', 'newmark', 'taylor-newmark']
    widest = 0
    do i = 1, size(every)
       if (.not. lists(cfg % families, trim(every(i)))) cycle
       do order = 1, cfg % max_discretization_order
          call chosen(trim(every(i)), order, scheme, staged, passes_check)
          if (.not. passes_check) cycle
          depth = scheme % history_depth(cfg % state_degree)
          if (depth < cfg % instants) widest = max(widest, depth)
       end do
    end do
  end function widest_depth
  subroutine chosen(name, order, scheme, staged, passes_check)
    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: staged, passes_check
    call family_named(name, order, scheme, passes_check)
    staged = .false.
    if (passes_check) staged = scheme % num_stages() > 1
  end subroutine chosen
  function labelled(names, orders) result(label)
    character(len=*), intent(in) :: names(:)
    integer         , intent(in) :: orders(:)
    character(len=:), allocatable :: label
    character(len=2) :: digit
    integer :: b
    label = ''
    do b = 1, size(names)
       write(digit,'(i0)') orders(b)
       if (b > 1) label = label // '-'
       label = label // trim(names(b)) // trim(digit)
    end do
  end function labelled
  subroutine shown_initial(cfg)
    type(configuration), intent(in) :: cfg
    real(dp), allocatable :: q(:)
    character(len=:), allocatable :: line
    character(len=18) :: cell
    character(len=12) :: counted
    integer :: i
    q = q0(1:cfg % state_degree + 1)
    line = '   initial state, consistent'
    do i = 1, size(q)
       write(cell,'(es18.10)') q(i)
       line = line // cell
    end do
    if (over_field) then
       write(counted,'(i0)') nodes
       line = line // '   at node 1 of ' // trim(counted)
    end if
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
  subroutine show_row(label, solved, f, final_imbalance, columns)
    character(len=*), intent(in) :: label
    integer         , intent(in) :: solved
    real(dp)        , intent(in) :: f(0:)
    type(imbalance) , intent(in) :: final_imbalance
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
    if (.not. final_imbalance % converged) then
       if (final_imbalance % diverging) then
          line = line // '   diverging'
       else
          line = line // '   unconverged'
       end if
    end if
    write(*,'(a)') line
    if (.not. final_imbalance % converged) call shown_aspect(final_imbalance)
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
  subroutine shown_aspect(final_imbalance)
    type(imbalance), intent(in) :: final_imbalance
    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: d
    write(*,'(a,es10.3,a,es10.3,a)') '      imbalance ', final_imbalance % norm, &
         & ' against ', final_imbalance % began, ' at the start of the march'
    line = '      by degree '
    do d = 0, ubound(final_imbalance % by_degree, 1)
       write(cell,'(es14.3)') final_imbalance % by_degree(d)
       line = line // cell
    end do
    write(*,'(a)') line
    write(*,'(a,i0,a,i0)') '      largest entry at slot ', final_imbalance % largest_slot, &
         & ' degree ', final_imbalance % largest_degree
    write(*,'(a,i0,a,i0,a,es10.3)') '      steepest in the state at slot ', &
         & final_imbalance % steepest_slot, ' degree ', final_imbalance % steepest_degree, &
         & ', d||r||/dq = ', final_imbalance % steepest
  end subroutine shown_aspect
  !===================================================================!
  ! The governing law the run marches. One procedure defines it, so
  ! every march, tangent and adjoint reads the same expression.
  !===================================================================!
  function physics_of(cfg) result(r)
    type(configuration), intent(in) :: cfg
    type(expression) :: r
    r = van_der_pol(cfg % state_degree)
  end function physics_of
  !===================================================================!
  ! The coordinates the run's state is declared over, read from the
  ! law it marches. Every other rule the run evaluates - a functional,
  ! a measure - is declared over the same, since all of them read the
  ! same state.
  !===================================================================!
  function state_degrees_of(cfg) result(degrees)
    type(configuration), intent(in) :: cfg
    integer, allocatable :: degrees(:)
    type(expression) :: law
    type(continuous_domain) :: continuous
    integer :: beside
    law        = physics_of(cfg)
    continuous = continuous_domain(law)
    beside     = continuous % num_components() - (continuous % equation_degree() + 1)
    if (beside > 0) then
       degrees = [continuous % equation_degree(), beside]
    else
       degrees = [continuous % equation_degree()]
    end if
  end function state_degrees_of
  subroutine one_row(cfg, names, orders, printed)
    type(configuration), intent(in)    :: cfg
    character(len=*)   , intent(in)    :: names(:)
    integer            , intent(in)    :: orders(:)
    integer            , intent(inout) :: printed
    type(family_container), allocatable :: schemes(:)
    type(chain_block)  , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    integer , allocatable :: added(:)
    real(dp), allocatable :: dt(:), t(:), f(:,:), weights(:)
    type(imbalance) :: final_imbalance
    real(dp) :: achieved
    integer :: nd, width, given, reported, i, m
    logical :: passes_check
    character(len=20) :: cell
    character(len=:), allocatable :: line
    nd    = cfg % state_degree + 1
    width = nd * nodes
    allocate(schemes(size(names)), added(size(names)))
    call assembled(cfg, names, orders, schemes, added, passes_check)
    if (.not. passes_check) return
    call grid_partition(cfg, dt, t)
    given = schemes(1) % scheme % history_depth(nd - 1)
    call tally_enter(at_expansion)
    call tally_order(0)
    if (grid_designed) then
       weights = dt(2:cfg % instants)
       call march_chain(schemes, added, physics_of(cfg), nd, &
            & designed_grid(cfg % time_duration), cfg % design, q0, chain, tower, dt, t, &
            & achieved, grid_design=weights, final_imbalance=final_imbalance, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else if (grid_adaptive) then
       call march_chain(schemes, added, physics_of(cfg), nd, &
            & fixed_grid(adaptive_weights), cfg % design, q0, chain, tower, dt, t, achieved, &
            & final_imbalance=final_imbalance, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else
       call march_chain(schemes, added, physics_of(cfg), nd, &
            & chosen_grid(cfg), cfg % design, q0, chain, tower, dt, t, achieved, final_imbalance=final_imbalance, &
            & nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, startup=cfg % startup_refinement)
    end if
    if (.not. final_imbalance % converged) then
       reported = 0
    else
       reported = cfg % max_derivative_degree
    end if
    call chain_expansion(chain, tower, functionals, nd, reported, f, node_measure=volume)
    call tally_leave()
    call show_row(labelled(names, orders), cfg % instants - given, f(:, 1), final_imbalance, &
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
    if (reported >= 1 .and. (grid_designed .or. lists(cfg % check, 'passes') &
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
    integer, allocatable :: versions(:)
    real(dp), allocatable :: p(:), df(:,:), other(:,:), table(:,:), entries(:,:,:)
    type(sink_costates) :: sinks
    real(dp) :: euler
    integer  :: num_designs, num_functionals, pass_kind, i, order
    num_functionals = size(functionals)
    call chain_versions(chain, tower, functionals, nd, versions, node_measure=volume)
    num_designs = num_designs_of(tower)
    pass_kind = pass_of(num_designs, num_functionals, 1)
    call chain_derivative(chain, tower, versions, functionals, nd, 1, pass_kind, df, node_measure=volume)
    write(*,'(a,a,a,i0,a,i0,a,es10.2)') '      first derivatives by the ', &
         & trim(merge('forward', 'reverse', pass_kind == forward_pass)), ' pass, designs ', &
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
    if (lists(cfg % check, 'passes')) then
       call chain_derivative(chain, tower, versions, functionals, nd, 1, &
            & merge(reverse_pass, forward_pass, pass_kind == forward_pass), other, node_measure=volume)
       write(*,'(a,es10.2)') '      tangent against adjoint over the table, relative ', &
            & maxval(abs(df - other)) / max(1.0_dp, maxval(abs(df)))
    end if
    if (lists(cfg % check, 'sinks')) then
       call chain_derivative(chain, tower, versions, functionals, nd, 1, reverse_pass, other, &
            & node_measure=volume, sinks=sinks)
       call shown_sinks(sinks, nd)
    end if
    if (grid_designed) then
       do order = 2, ubound(f, 1)
          pass_kind = pass_of(num_designs, num_functionals, order)
          call chain_derivative(chain, tower, versions, functionals, nd, order, pass_kind, table, &
               & node_measure=volume, entries=entries)
          do i = 1, num_functionals
             if (pass_kind == reverse_pass) then
                write(*,'(a,i0,a,a,i0,a,es12.4,a,es10.2,a,es10.2)') '      derivatives of order ', &
                     & order, ' by the reverse pass, functional ', '', i, ':  |T| ', &
                     & maxval(abs(table(i, :))), '   departure among the entries of a multiset ', &
                     & asymmetry(entries, num_designs, order), &
                     & '   parameter entry against the expansion ', &
                     & abs(table(i, 1) - f(order, i)) / max(1.0_dp, abs(f(order, i)))
             else
                write(*,'(a,i0,a,i0,a,es12.4,a,es10.2)') '      derivatives of order ', order, &
                     & ' by the forward pass, functional ', i, ':  |T| ', maxval(abs(table(i, :))), &
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
         & ' at the last point', trim(counted(nd, sinks % fixed_rows)), ' in fixed rows'
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
    logical :: passes_check
    integer :: i
    call refuse_unknown(cfg % designs, ['physics', 'grid   '], 'designs')
    call refuse_unknown(cfg % functionals, ['energy     ', 'dissipation'], 'functionals')
    if (.not. lists(cfg % designs, 'physics')) then
       error stop 'graph_time_integrator: the physics'' parameter is the first design'
    end if
    grid_designed = lists(cfg % designs, 'grid')
    names = words_of(cfg % functionals)
    allocate(functionals(size(names)))
    do i = 1, size(names)
       call functional_named(trim(names(i)), cfg % state_degree, functionals(i), passes_check)
       if (passes_check) functionals(i) = stated_over(functionals(i), state_degrees_of(cfg), functionals(i) % name())
    end do
  end subroutine chosen_functionals
  subroutine against_the_ode(cfg, schemes, added, f_field)
    type(configuration), intent(in) :: cfg
    type(family_container), intent(in) :: schemes(:)
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
    call march_chain(schemes, added, physics_of(cfg), nd, chosen_grid(cfg), &
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
    type(expression) :: law
    type(continuous_domain) :: continuous
    real(dp) :: x, y, began
    integer  :: n1, n2
    call refuse_unknown(cfg % initial_field, ['constant', 'mode    ', 'bump    '], 'initial_field')
    call refuse_unknown(cfg % export, ['none    ', 'paraview'], 'export')
    call refuse_unknown(cfg % check, ['none    ', 'ode     ', 'mode    ', 'operator', 'passes  ', &
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
    ! the initial state stores one value per component the rule reads,
    ! a count the rule itself declares
    law = physics_of(cfg)
    continuous = continuous_domain(law)
    q0 = initial_field(law, continuous % num_components(), &
         & cfg % initial_field, cfg % initial_state, cfg % design, &
         & spatial_discretization_stencil=spatial_discretization_stencil, space=space, a=extent_a, b=extent_b)
  end subroutine field_context
  subroutine pair_of(pair, x, y, subject)
    character(len=*), intent(in)  :: pair, subject
    real(dp)        , intent(out) :: x, y
    character(len=32), allocatable :: w(:)
    w = words_of(pair)
    if (size(w) /= 2) error stop 'graph_time_integrator: two ' // subject // ', one per coordinate'
    read(w(1), *) x
    read(w(2), *) y
  end subroutine pair_of
  subroutine assembled(cfg, names, orders, schemes, added, passes_check)
    type(configuration), intent(in)  :: cfg
    character(len=*)   , intent(in)  :: names(:)
    integer            , intent(in)  :: orders(:)
    type(family_container), intent(inout) :: schemes(:)
    integer            , intent(inout) :: added(:)
    logical            , intent(out)   :: passes_check
    class(family), allocatable :: scheme
    logical :: staged, exists
    integer :: b, blocks, share
    blocks = size(names)
    passes_check = .true.
    share = cfg % instants / blocks
    added = share
    added(1) = cfg % instants - share * (blocks - 1)
    do b = 1, blocks
       call chosen(names(b), orders(b), scheme, staged, exists)
       if (.not. exists) then
          passes_check = .false.
          cycle
       end if
       allocate(schemes(b) % scheme, source=scheme)
       deallocate(scheme)
       if (added(b) <= schemes(b) % scheme % history_depth(cfg % state_degree)) passes_check = .false.
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
    if (trim(cfg % adaptive_check) == 'goal_oriented') then
       adaptive_weights = goal_oriented_partition(crouzeix_three_stage(), &
            & physics_of(cfg), van_der_pol_energy(cfg % state_degree), nd, &
            & cfg % time_duration, q0(1:cfg % state_degree), cfg % design, cfg % tolerance, &
            & trim(cfg % tolerance_criterion) == 'relative', rejects)
    else
       adaptive_weights = adaptive_partition(crouzeix_three_stage(), 4, &
            & physics_of(cfg), nd, cfg % time_duration, &
            & q0(1:cfg % state_degree), cfg % design, cfg % tolerance, &
            & trim(cfg % tolerance_criterion) == 'relative', rejects)
    end if
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
    call refuse_unknown(cfg % adaptive_check, ['step_doubling', 'goal_oriented'], &
         & 'adaptive_check')
    call set_stopping(cfg % tolerance, &
         & merge(relative, absolute, trim(cfg % tolerance_criterion) == 'relative'), &
         & merge(by_rate, by_count, trim(cfg % iteration_criterion) == 'by_rate'), &
         & cfg % max_iterations)
    call set_linear_budget(cfg % krylov_restart, cfg % smoothing_sweeps, &
         & cfg % max_linear_iterations)
    call adaptive_context(cfg)
    if (cfg % accounting) then
       call refuse_unknown(cfg % measurements, &
            & ['elapsed_time  ', 'primal_loops  ', 'tangent_loops ', &
            &  'adjoint_loops ', 'newton_solves ', 'linear_solves ', &
            &  'factorisations'], 'measurements')
    end if
    call refuse_unknown(cfg % families, &
         & [character(len=16) :: 'bdf', 'adams', 'dirk', 'newmark', 'taylor-newmark'], 'families')
    call refuse_unwindowed(cfg % combinations)
    widest = widest_depth(cfg)
    if (.not. cfg % automatic_order_conservation) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' the widest row here has history depth in instants: ', widest
       write(*,'(a)')    ' filling them by any other means leaves the rows solving different'
       write(*,'(a)')    ' problems from different initial states, and no table compared across'
       write(*,'(a)')    ' such rows is meaningful.'
       error stop 'graph_time_integrator: order conservation is the only startup built'
    end if
    if (widest == 0) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' every family and order requested has a history depth beyond the'
       write(*,'(a,i0)') ' horizon, whose instant count is: ', cfg % instants
       error stop 'graph_time_integrator: no row fits in this horizon'
    end if
    call grid_partition(cfg, dt, t)
    call shown_initial(cfg)
    write(*,'(a,a)') '   precision of this build  ', precision_named()
    call heading(cfg)
    if (cfg % accounting) call tally_open(cfg % max_derivative_degree, hierarchy_levels)
    printed = 0
    call every_window_count(cfg, printed)
    call the_named_chain(cfg, printed)
    if (cfg % accounting) then
       call tally_close()
       call accounted(cfg)
    end if
    if (printed == 0) then
       write(*,'(a)') ' '
       write(*,'(a)') ' no row was built. A family has no scheme at every order - a stage'
       write(*,'(a)') ' family has none below order two - and a row whose blocks would add'
       write(*,'(a)') ' no more instants than their history depth is not built either.'
    end if
  end subroutine table
  !===================================================================!
  ! HOW MANY WINDOWS EACH SURVEYED CHAIN HAS. The setting lists
  ! the counts, so a survey over chains of one, two and three windows
  ! reads `combinations = 1 2 3`, and any count may be requested. A
  ! count above the number of families named yields no chain, because
  ! a surveyed chain gives each window a family of its own; a chain
  ! that repeats a family is specified explicitly through `chain`.
  !===================================================================!

  subroutine every_window_count(cfg, printed)
    type(configuration), intent(in)    :: cfg
    integer            , intent(inout) :: printed
    character(len=32), allocatable :: counts(:)
    integer :: i, windows
    counts = words_of(cfg % combinations)
    do i = 1, size(counts)
       read(counts(i), *) windows
       call tuple_rows(cfg, windows, printed)
    end do
  end subroutine every_window_count

  !===================================================================!
  ! THE CHAIN SPECIFIED EXPLICITLY, IF ONE IS. One window per word, each
  ! reading family:order. The instants divide evenly among however
  ! many windows are named, and a window too short for the family
  ! assigned to it is reported rather than skipped - a chain requested
  ! by name is not silently omitted the way a surveyed one is.
  !===================================================================!

  subroutine the_named_chain(cfg, printed)
    type(configuration), intent(in)    :: cfg
    integer            , intent(inout) :: printed
    character(len=16), allocatable :: names(:)
    integer         , allocatable :: orders(:)
    integer :: windows, before
    if (len_trim(cfg % chain) < 1) return
    call windows_of(cfg % chain, names, orders)
    windows = size(names)
    if (cfg % instants / windows <= 0) then
       write(*,'(a)') ' '
       write(*,'(a,i0,a,i0,a)') ' the chain names ', windows, ' windows and the horizon contains ', &
            & cfg % instants, ' instants, so a window would contain none.'
       error stop 'graph_time_integrator: a window of a named chain contains instants'
    end if
    before = printed
    call one_row(cfg, names, orders, printed)
    if (printed == before) then
       write(*,'(a)') ' '
       write(*,'(a)') ' the chain ' // trim(cfg % chain) // ' built no row. A family has no'
       write(*,'(a)') ' scheme at the order requested of it, or a window adds no more instants'
       write(*,'(a)') ' than the history depth of the family assigned to it.'
       error stop 'graph_time_integrator: a named chain builds its row'
    end if
  end subroutine the_named_chain

  !===================================================================!
  ! THE WINDOWS A CHAIN NAMES, read as family and order. A word
  ! without a colon, or an order that is not a whole number, is
  ! refused where it is written.
  !===================================================================!

  subroutine windows_of(specification, names, orders)
    character(len=*), intent(in) :: specification
    character(len=16), allocatable, intent(out) :: names(:)
    integer         , allocatable, intent(out) :: orders(:)
    character(len=32), allocatable :: words(:)
    integer :: i, mark, failed
    words = words_of(specification)
    if (size(words) < 1) then
       error stop 'gti_configuration: a chain names a window at least'
    end if
    allocate(names(size(words)), orders(size(words)))
    do i = 1, size(words)
       mark = index(words(i), ':')
       if (mark < 2 .or. mark >= len_trim(words(i))) then
          write(*,'(a)') ' '
          write(*,'(a)') ' the chain names ' // trim(words(i)) // &
               & ', which is not a family and an order.'
          error stop 'gti_configuration: a setting names something unknown'
       end if
       names(i) = words(i)(1:mark - 1)
       read(words(i)(mark + 1:), *, iostat=failed) orders(i)
       if (failed /= 0 .or. orders(i) < 1) then
          write(*,'(a)') ' '
          write(*,'(a)') ' the chain requests ' // trim(words(i)) // &
               & ' for an order that is not a whole number of one or more.'
          error stop 'gti_configuration: a setting names something unknown'
       end if
    end do
    call refuse_unknown(names_phrase(names), &
         & [character(len=16) :: 'bdf', 'adams', 'dirk', 'newmark', 'taylor-newmark'], 'chain')
  end subroutine windows_of

  pure function names_phrase(names) result(phrase)
    character(len=*), intent(in) :: names(:)
    character(len=:), allocatable :: phrase
    integer :: i
    phrase = ''
    do i = 1, size(names)
       phrase = phrase // ' ' // trim(names(i))
    end do
  end function names_phrase

  !===================================================================!
  ! A window count is a whole number of one or more. Anything else is
  ! refused where it is written rather than where it would be used.
  !===================================================================!

  subroutine refuse_unwindowed(combinations)
    character(len=*), intent(in) :: combinations
    character(len=32), allocatable :: counts(:)
    integer :: i, windows, failed
    counts = words_of(combinations)
    if (size(counts) < 1) then
       error stop 'gti_configuration: the combinations name a window count'
    end if
    do i = 1, size(counts)
       read(counts(i), *, iostat=failed) windows
       if (failed /= 0 .or. windows < 1) then
          write(*,'(a)') ' '
          write(*,'(a)') ' combinations names ' // trim(counts(i)) // &
               & ', which is not a count of windows.'
          error stop 'gti_configuration: a setting names something unknown'
       end if
    end do
  end subroutine refuse_unwindowed
  function listed(cfg) result(list)
    type(configuration), intent(in) :: cfg
    character(len=16), allocatable :: list(:)
    character(len=16) :: every(5)
    integer :: i, n
    every = [character(len=16) :: 'bdf', 'adams', 'dirk', 'newmark', 'taylor-newmark']
    n = 0
    do i = 1, size(every)
       if (lists(cfg % families, trim(every(i)))) n = n + 1
    end do
    allocate(list(n))
    n = 0
    do i = 1, size(every)
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
    character(len=16), allocatable :: names(:)
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
    character(len=32), allocatable :: requested(:)
    integer :: i, event
    requested = words_of(cfg % measurements)
    do i = 1, size(requested)
       event = tally_event_of(trim(requested(i)))
       call one_measurement(cfg, event)
    end do
    call substitution_note(cfg)
    call iteration_limit_note(cfg)
  end subroutine accounted
  subroutine substitution_note(cfg)
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
  end subroutine substitution_note
  subroutine one_measurement(cfg, event)
    type(configuration), intent(in) :: cfg
    integer            , intent(in) :: event
    real(dp), allocatable :: whole(:)
    character(len=:), allocatable :: line
    integer :: level, m, top
    top = cfg % max_derivative_degree
    allocate(whole(0:top), source=0.0_dp)
    do m = 0, top
       if (event == elapsed_time) then
          whole(m) = tally_amount(at_horizon, m, event)
       else
          whole(m) = over_levels(m, event)
       end if
    end do
    write(*,'(a)') ' '
    write(*,'(a)') ' accounting: ' // tally_event_name(event)
    if (event == elapsed_time) then
       write(*,'(a)') '   seconds. A level contains the levels opened inside it, and the'
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
          line = line // right(amount_cell(tally_amount(level, m, event), event))
       end do
       write(*,'(a)') line
    end do
    if (event == elapsed_time) then
       line = '   per order        '
    else
       line = '   whole run        '
    end if
    do m = 0, top
       line = line // right(amount_cell(whole(m), event))
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
  subroutine iteration_limit_note(cfg)
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
       write(*,'(a)') '   at the iteration limit: the march is not converging, so'
       write(*,'(a)') '   these ratios report the limit and not the derivative order.'
    end if
  end subroutine iteration_limit_note
  function amount_cell(elapsed, event) result(entry)
    real(dp), intent(in) :: elapsed
    integer , intent(in) :: event
    character(len=:), allocatable :: entry
    character(len=14) :: cell
    if (event == elapsed_time) then
       write(cell,'(f14.4)') elapsed
    else
       write(cell,'(i14)') nint(elapsed)
    end if
    entry = trim(adjustl(cell))
  end function amount_cell
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
  function right(entry) result(cell)
    character(len=*), intent(in) :: entry
    character(len=14) :: cell
    write(cell,'(a14)') entry
  end function right
end program graph_time_integrator
