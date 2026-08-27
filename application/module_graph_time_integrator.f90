!=====================================================================!
! Packed graph-time-integrator application.
!
! The finite-volume library remains in ../src. This source holds the
! graph-time-integrator application modules followed by the one main
! program.
!=====================================================================!

!=====================================================================!
! packed from physics/physics_vanderpol.f90
!=====================================================================!
!=====================================================================!
! The van der Pol oscillator as a governing constraint, and two
! functional integrands beside it, each stated once at one instant.
!
! The equation, at degree N,
!
!      R  =  q^(N)  -  nu (1 - q^2) q^(N-1)  +  q  =  0
!
! is the ordinary oscillator at N = 2 and its higher-degree
! continuation above that: the damping acts on the derivative one
! below the highest, and the restoring term on the value. The energy
! is
!
!      F  =  ( q^2 + (q')^2 ) / 2
!
! and the dissipation, the power the damping term draws,
!
!      F  =  nu (1 - q^2) (q')^2 ,
!
! a functional that reads the design itself, so its own partial in
! the design is not zero.
!
! Each is an expression over the unknown and the design, and its
! partials in either are taken by evaluating it; none is written out
! here. A functional reads the velocity, so a degree-zero problem has
! none: stating it at degree zero stops the program.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_vanderpol

  use util_precision    , only : dp
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

end module physics_vanderpol

!=====================================================================!
! packed from application/gti_configuration.f90
!=====================================================================!
!=====================================================================!
! What a run was asked for.
!
! One record of settings, read from a named file and then overridden
! by whatever the command line repeats. A setting names itself, so a
! file reads
!
!      time_duration = 7.0
!      max_derivative_degree = 4
!
! and the same setting on the command line reads
!
!      --time-duration=7.0
!
! the two spellings differing only in the separator, so that neither
! has to be learnt twice.
!
!             WHAT IS REFUSED
!
! A setting that is not one of those below, whether it comes from a
! file or from the command line: a run that silently ignores what it
! was told is worse than one that stops. A file that cannot be
! opened. A value that is not of the setting's kind.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_configuration

  use util_precision  , only : dp

  implicit none

  private
  public :: configuration, read_configuration, override, show
  public :: worded, lists, refuse_unknown, chosen_from

  type :: configuration

     character(len=32) :: physics      = 'vanderpol'
     ! grid is uniform, random, or adaptive - the last discovered by an
     ! error-controlled march to the tolerance, the steps then frozen
     character(len=32) :: grid         = 'random'
     character(len=32) :: combinations = 'homogeneous'
     character(len=32) :: families     = 'bdf adams dirk'

     integer  :: state_degree             = 2
     integer  :: instants                 = 21
     integer  :: max_derivative_degree    = 3
     integer  :: max_discretization_order = 4
     integer  :: seed                     = 20260824
     integer  :: startup_refinement       = 4
     ! HOW A LINEAR SYSTEM IS SOLVED, four specifications each its
     ! own: linear_solver direct or iterative; assembly matrix or
     ! free; storage dense or sparse; multigrid on or off, the solver
     ! named then smoothing a two-grid over aggregates.
     character(len=16) :: linear_solver   = 'direct'
     character(len=16) :: assembly        = 'matrix'
     character(len=16) :: storage         = 'dense'
     logical           :: multigrid       = .false.

     ! WHICH LEVEL IS SWEPT to solve a block: space-time solves the
     ! whole block at once, time sweeps the instants, space sweeps the
     ! nodes.
     character(len=16) :: sweep           = 'space-time'

     ! THE DESIGNS AND THE FUNCTIONALS. designs lists what the
     ! functionals are differentiated in: physics, the equation's
     ! parameter, and grid, the weights of the steps on the simplex.
     ! functionals lists what is integrated over the horizon.
     character(len=64) :: designs         = 'physics'
     character(len=64) :: functionals     = 'energy'

     real(dp) :: time_duration = 7.0_dp
     real(dp) :: design        = 1.0_dp

     logical  :: automatic_order_conservation = .true.
     logical  :: mixed_orders                 = .false.

     ! THE STATE AT THE FIRST INSTANT, below the highest derivative.
     ! The highest is not given: it is what the physics says it is
     ! there, solved for. Fewer values than degrees are taken as the
     ! rest being zero.
     character(len=128) :: initial_state = '1.0'

     ! THE LEVEL BELOW. A spatial mesh under every instant: its
     ! geometry, its extents along the two coordinates, the cell
     ! counts along them, and its spacing, uniform or drawn from the
     ! seed as the time grid is. Counts of zero mean no spatial level.
     ! diffusion is the kappa on the laplacian. initial_field is
     ! constant, from initial_state at every node, or the rectangle's
     ! mode. export writes every instant for paraview; check names a
     ! comparison the run makes against something known: ode, mode,
     ! operator, routes, or sinks - the costates of the unknowns no row
     ! reads, J_ii lambda_i = g_i, at any degree.
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

     ! HOW A MARCH STOPS. The tolerance is a ratio where the
     ! criterion is relative, which is the question asked whenever the
     ! target is a reduction in imbalance, and a number in its own
     ! right where the criterion is absolute, which is the question
     ! asked only where that number means something on its own. The
     ! iteration criterion says whether the budget is a count or is
     ! taken from the rate the march itself shows.
     real(dp)          :: tolerance           = 1.0e-12_dp
     character(len=16) :: tolerance_criterion = 'relative'
     character(len=16) :: iteration_criterion = 'by_rate'
     integer           :: max_iterations      = 100

     ! THE INNER SOLVES honour the same tolerance, criterion and budget
     ! kind. What they add is declared here: the width of a Krylov
     ! space before it restarts, the smoothing sweeps a multigrid
     ! cycle takes, and the ceiling on cycles or restarts.
     integer           :: krylov_restart        = 60
     integer           :: smoothing_sweeps      = 2
     integer           :: max_linear_iterations = 200

     ! Whether the run counts what it spends, and which of the counts
     ! it prints. Counting is off unless it is asked for.
     logical  :: accounting                   = .false.
     character(len=256) :: measurements       = &
          & 'wall_time primal_loops tangent_loops adjoint_loops ' // &
          & 'newton_solves linear_solves factorisations'

  end type configuration

contains

  !===================================================================!
  ! A setting named the way a file names it: the separators a command
  ! line uses become the ones a file uses, and the case is levelled,
  ! so one name is read from either.
  !===================================================================!

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

  !===================================================================!
  ! One setting given a value. A name that is not a setting stops the
  ! program.
  !===================================================================!

  !===================================================================!
  ! The blank-separated words of a list, and whether one of them is a
  ! given word. Containment is not the test: a word is what lies
  ! between blanks, so a list naming bdfx does not thereby name bdf.
  !===================================================================!

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

  !===================================================================!
  ! A word this program has nothing for stops it. Left to run, the
  ! setting would take effect nowhere and the run would report a
  ! reason that is not the one.
  !===================================================================!

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
    case ('accounting')
       read(value, *) cfg % accounting
    case ('measurements')
       cfg % measurements = value
    case default
       write(*,'(a)') ' this is not a setting: ' // levelled(name)
       error stop 'gti_configuration: every setting given is one that exists'
    end select

    ! A count below the least value it means anything at is refused
    ! here, where it is given, rather than where it is first indexed
    ! by or counted over.
    call refuse_below(cfg % max_iterations, 1, &
         & 'max_iterations', 'an iteration budget is at least one')
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

  !===================================================================!
  ! One line of a file: a comment, a blank, or a setting and its
  ! value either side of an equals sign.
  !===================================================================!

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

  !===================================================================!
  ! Every setting a named file gives. The file is looked for under
  ! config, named for the configuration and ending in cfg.
  !===================================================================!

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

  !===================================================================!
  ! One command-line argument, which is a setting spelt with dashes
  ! and given its value after an equals sign.
  !===================================================================!

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
    write(*,'(a,i0)')      '   krylov restart           ', cfg % krylov_restart
    write(*,'(a,i0)')      '   smoothing sweeps         ', cfg % smoothing_sweeps
    write(*,'(a,i0)')      '   max linear iterations    ', cfg % max_linear_iterations
    write(*,'(a,l1)')      '   accounting               ', cfg % accounting
    if (cfg % accounting) then
       write(*,'(a,a)')    '   measurements             ', trim(cfg % measurements)
    end if

  end subroutine show

  !===================================================================!
  ! Which of the words listed a setting names, as its place in the
  ! list; a word the list has not stops the program the way every
  ! unknown word does.
  !===================================================================!

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
! packed from application/gti_sweeps.f90
!=====================================================================!
!=====================================================================!
! The sensitivity of a functional to a design, by the tangent and
! by the adjoint.
!
! A block's statement R(q, x) = 0 determines q from x, and the
! functional is a sum over the instants,
!
!      f  =  sum over k of  dt_k F(q_k, x) .
!
! Differentiating the statement gives J dq/dx = -dR/dx, so
!
!      tangent    solve J w = -dR/dx  once, then  df/dx = f_x + g.w
!      adjoint    solve J^T l = g     once, then  df/dx = f_x - l.dR/dx
!
! with J the jacobian in the state, g the gradient of f in the state
! and f_x its own partial in the design. The two read the same
! three objects and must agree to round-off; that they do is what
! this module exists to make checkable.
!
! Each of the three comes from a partial action that is exact: the
! scheme's rows are linear and are their own jacobian, and the
! physics and the integrand differentiate their own rules. Nothing
! here is differenced.
!
!             HOW THE JACOBIAN IS FORMED
!
! For the adjoint, column by column: one partial action per unknown,
! and then solved densely. That does not scale, and what would is
! the stencil's compiled transpose rather than a dense one.
!
! For the tangent no matrix is formed at all past a large block. The
! statement's partial action is already a matvec, so freezing it at
! the trajectory gives a linear operation a krylov solver can be
! driven with directly. Below the same threshold gti_march uses, a
! dense factorisation is the faster of the two and is taken.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_sweeps

  use gti_configuration, only : refuse_unknown
  use util_precision  , only : dp, half_digits
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use graph_fractal         , only : graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_dense_direct, only : dense_direct
  use operation_multigrid   , only : multigrid
  use operation_gauss_seidel, only : gauss_seidel
  use operation_gmres       , only : gmres
  use operation_minimization, only : minimizer, relative, absolute, by_rate, by_count

  implicit none

  private
  public :: functional_of, functional_gradient
  public :: design_partial, jacobian_of
  public :: forward_route, reverse_route, route_of, route_substitutions, choose

  integer, parameter :: forward_route = 1
  integer, parameter :: reverse_route = 2
  public :: set_linear_solver, set_assembly, set_storage, set_multigrid
  public :: set_aggregates, set_coarse_nodes, coarse_nodes, assembly_present, multigrid_on
  public :: take_inner, keep_inner, forget_inner, set_linear_stopping, set_linear_budget

  !===================================================================!
  ! HOW A LINEAR SYSTEM IS SOLVED: four specifications, each its own.
  !
  !      linear_solver   direct | iterative    factorise, or iterate
  !      assembly        matrix | free         a matrix is formed, or
  !                                            only a matvec is attached
  !      storage         dense | sparse        the matrix formed is a
  !                                            square, or a stencil
  !      multigrid       yes | no              a two-grid over aggregates:
  !                                            block gauss-seidel over a
  !                                            point's components smooths,
  !                                            the solver named serves
  !                                            the coarse level
  !
  ! The corners that mean nothing are refused by name: direct on a
  ! free assembly, sparse direct (not built), dense iterative, and
  ! multigrid on a free assembly, its coarse operator being read
  ! through the aggregates from a stencil.
  !===================================================================!

  character(len=16), save :: chosen_solver   = 'direct'
  character(len=16), save :: chosen_assembly = 'matrix'
  character(len=16), save :: chosen_storage  = 'dense'
  logical          , save :: chosen_multigrid = .false.
  integer, allocatable, save :: chosen_aggregates(:)

  ! WHAT THE INNER SOLVES STOP AT: the tolerance, its criterion and
  ! the budget kind are the march's own, given once; the Krylov
  ! restart, the smoothing sweeps and the ceiling are declared.
  real(dp), save :: linear_tolerance  = half_digits
  integer , save :: linear_criterion  = relative
  integer , save :: linear_budget     = by_rate
  integer , save :: linear_restart    = 60
  integer , save :: linear_sweeps     = 2
  integer , save :: linear_iterations = 200
  ! the coarse cell of every node, which a block's aggregates are read from
  integer, allocatable, save :: chosen_coarse(:)

  ! The inner minimizer kept between solves, so that a direct one
  ! keeps its factors across the statements stamped alike.
  class(minimizer), allocatable, save :: kept_inner

contains

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

  pure logical function assembly_present() result(yes)

    yes = trim(chosen_assembly) == 'matrix'

  end function assembly_present

  pure logical function multigrid_on() result(yes)

    yes = chosen_multigrid

  end function multigrid_on

  !===================================================================!
  ! The aggregates multigrid coarsens by: one block per unknown. A
  ! caller with no field clears them.
  !===================================================================!

  subroutine set_aggregates(aggregates)

    integer, intent(in), optional :: aggregates(:)

    if (allocated(chosen_aggregates)) deallocate(chosen_aggregates)
    if (present(aggregates)) chosen_aggregates = aggregates

  end subroutine set_aggregates

  !===================================================================!
  ! The coarse cell of every node, from which a block reads the
  ! aggregates it is coarsened by. None given, every node is its own
  ! coarse cell and the coarse level is the fine one.
  !===================================================================!

  subroutine set_coarse_nodes(cell)

    integer, intent(in), optional :: cell(:)

    if (allocated(chosen_coarse)) deallocate(chosen_coarse)
    if (present(cell)) chosen_coarse = cell

  end subroutine set_coarse_nodes

  function coarse_nodes(nodes) result(cell)

    ! nodes: the largest node label the map must reach
    integer, intent(in) :: nodes
    integer, allocatable :: cell(:)

    integer :: i

    ! a member of a block keeps the block's node labels, so the map
    ! given must reach the largest of them
    if (allocated(chosen_coarse)) then
       if (size(chosen_coarse) < nodes) then
          error stop 'gti_sweeps: a coarse cell for every node'
       end if
       cell = chosen_coarse
    else
       cell = [(i, i = 1, nodes)]
    end if

  end function coarse_nodes

  !===================================================================!
  ! What the inner solves stop at: the march's own tolerance,
  ! criterion and budget kind, handed over by the march's stopping;
  ! and what is declared, the restart, the sweeps and the ceiling.
  ! Invalid input: a tolerance that is not positive, a criterion or
  ! budget kind that is neither, a restart, sweep count or ceiling
  ! below one.
  !===================================================================!

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

  !===================================================================!
  ! The minimizer the specifications name, built for a system of the
  ! given count. A singular pivot in a direct one is reported, the
  ! matrix being a tangent at an intermediate iterate.
  !===================================================================!

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

    ! gauss-seidel smooths, a point's components at a time - the
    ! coupling within a point being what no point smoother damps -
    ! and the solver named serves the coarse level
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

  !===================================================================!
  ! The inner minimizer taken for a solve and kept after it, so that
  ! what it holds - a direct solver's factors - outlives one solve.
  !===================================================================!

  subroutine take_inner(inner, count, width)

    class(minimizer), allocatable, intent(out) :: inner
    integer                      , intent(in)  :: count, width

    ! multigrid is built afresh for every statement, its aggregates
    ! being the statement's; a kept one would carry another's
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

  subroutine applied(action, on, inputs, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp), allocatable, intent(out) :: y(:)

    class(field), allocatable :: out

    call action % apply(on, inputs, out)
    call out % real_vector(y)

  end subroutine applied





  !===================================================================!
  ! The partial action of a statement along one direction in one of
  ! its arguments.
  !===================================================================!

  subroutine varied(action, on, inputs, which, domain, v, y, which2, domain2, v2)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: which
    type(graph)          , intent(in) :: domain
    real(dp)             , intent(in) :: v(:)
    real(dp), allocatable, intent(out) :: y(:)
    ! a second variation: the mixed second partial along both
    integer    , intent(in), optional :: which2
    type(graph), intent(in), optional :: domain2
    real(dp)   , intent(in), optional :: v2(:)

    type(stored_field) :: direction, second
    class(field), allocatable :: out

    direction = stored_field('direction', domain, size(v))
    call direction % set_real_vector(v)
    if (present(which2)) then
       second = stored_field('direction', domain2, size(v2))
       call second % set_real_vector(v2)
       call action % partial_action(on, inputs, &
            & [variation(action % argument(which), direction), &
            &  variation(action % argument(which2), second)], out)
    else
       call action % partial_action(on, inputs, &
            & [variation(action % argument(which), direction)], out)
    end if
    call out % real_vector(y)

  end subroutine varied

  !===================================================================!
  ! The functional: the integrand at every instant, weighted by the
  ! step that ends there. The first instant carries no step and so
  ! contributes nothing.
  !===================================================================!

  real(dp) function functional_of(integrand, instants, inputs, dt) result(f)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)

    real(dp), allocatable :: values(:)

    call applied(integrand, instants, inputs, values)
    f = sum(dt * values)

  end function functional_of

  !===================================================================!
  ! Its gradient in the state. The integrand at one instant reads
  ! only that instant's components, so one partial action per degree
  ! gives the whole gradient rather than one per unknown.
  !===================================================================!

  subroutine functional_gradient(integrand, instants, inputs, dt, n, degrees, &
       & state_domain, g, along_state, along_design)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)
    integer              , intent(in) :: n, degrees
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: g(:)
    ! given, the gradient's own partial along a state direction over
    ! the points, or along the design at every point
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

  !===================================================================!
  ! The statement's partial in the design, along the direction that
  ! varies every instant's design value together, which is a single
  ! design number for the whole block.
  !===================================================================!

  subroutine design_partial(rows, unknowns, inputs, n, design_domain, d)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: n
    type(graph)          , intent(in) :: design_domain
    real(dp), allocatable, intent(out) :: d(:)

    call varied(rows, unknowns, inputs, 2, design_domain, spread(1.0_dp, 1, n), d)

  end subroutine design_partial

  !===================================================================!
  ! The jacobian in the state, one column per unknown.
  !===================================================================!

  subroutine jacobian_of(rows, unknowns, inputs, num_unknowns, state_domain, a)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: num_unknowns
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: a(:,:)

    real(dp), allocatable :: v(:), column(:), w(:)
    integer , allocatable :: r(:), c(:)
    logical :: available
    integer :: j, e

    allocate(a(num_unknowns, num_unknowns), source=0.0_dp)

    ! A statement that compiles its tangent hands over its triples,
    ! and the square is filled from them; any other is probed one
    ! column at a time, a partial action each.
    call rows % compiled_tangent(unknowns, inputs, 1, r, c, w, available)
    if (available) then
       do e = 1, size(r)
          a(r(e), c(e)) = a(r(e), c(e)) + w(e)
       end do
       return
    end if

    allocate(v(num_unknowns), source=0.0_dp)
    do j = 1, num_unknowns
       v    = 0.0_dp
       v(j) = 1.0_dp
       call varied(rows, unknowns, inputs, 1, state_domain, v, column)
       a(:, j) = column
    end do

  end subroutine jacobian_of

  !===================================================================!
  ! THE GATE. Which route computes the m-th derivative of n_f
  ! functionals in n_d designs, by the count of substitutions each
  ! costs against one kept factorisation:
  !
  !      forward, m times            C(n_d + m - 1, m)
  !      forward m-1 times over      (1 + n_f) C(n_d + m - 2, m - 1)
  !      one reverse
  !
  ! The ratio of the first to the second is (n_d + m - 1) / (m (1 + n_f)),
  ! so the forward route is the cheaper exactly where n_d <= m n_f. At
  ! m = 1 that is the rule for a gradient, n_d <= n_f. At one design
  ! and one functional the forward route costs m substitutions and the
  ! other 2 m, at every order.
  !
  ! A count below one, or an order below one, stops the program.
  !===================================================================!

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

  !===================================================================!
  ! The substitutions a route costs at one order, per block. What the
  ! accounting layer counts is compared against this.
  !===================================================================!

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

!=====================================================================!
! packed from application/gti_expansion.f90
!=====================================================================!
!=====================================================================!
! The expansion: the graph a time integrator is, with what hangs on
! it.
!
! One object, one identity. Its branches are the six levels -
! expansion, sweep, horizon, block, slice, component - and everything
! that is not structure is kept in a map keyed on a node's identity:
! what it is called, what it holds, and how many members the set it
! denotes has. That is the arrangement view_mesh uses, one step
! further: a mesh attaches its measurements as components because a
! mesh has a fixed set of them, and the levels here do not, so they
! are attached by identity instead.
!
!             THE SLOTS
!
! Handed a physics, one family per block, the instants each block
! covers, a grid and a design, this builds a hierarchy in which every
! level is consistent and every coupling is relationally valid. The
! slots are the only things it is told; nothing else about the
! problem is written here.
!
!             WHAT IS NOT ASSIGNABLE
!
! The storage lends pointers into its own nodes, so an expansion is
! refused assignment: a copy would share them and either release
! would strand the other. That is why it is built into a variable
! rather than returned from one - a constructor's result would be
! assigned, and the assignment is what is refused.
!
!             THE ROWS A BLOCK CARRIES
!
! A block holds only the rows that fit inside it. A row on the d-th
! derivative at instant k reads a fixed number of instants back, and
! at the first instants of a block there are not that many, so those
! rows are absent and their components are carried in instead: they
! are marked as holding a value from the start, and every component
! after them waits on a march. Joining one block's last instants to
! the next block's first is the junction constraint, which is not
! built here.
!
!             WHAT IS REFUSED
!
! A block with fewer instants than its family reaches; a design the
! grid cannot read; an expansion built twice. Every level is checked
! as it is assembled, so a coupling that names another level's
! members stops the build where it is made.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

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
  use operation_grid        , only : grid
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

  !===================================================================!
  ! One family per block. Families of different kinds cannot share an
  ! array, so each is held in its own allocatable slot.
  !===================================================================!

  ! what a design leaf is read by: the physics, or the grid
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
     ! THE NODES a component holds - one for an equation at a point,
     ! the cells of a mesh for a field - and the spatial discretization stencil, one
     ! coupling over the nodes, laid on every component the physics
     ! sits on; zero where there is none
     integer                 , private :: node_extent = 1
     integer                 , private :: spatial_coupling_at = 0
     ! THE DESIGNS: leaves of their own level under the root, one per
     ! design - the physics' parameter, and the weights of the steps
     ! when they are designs - each holding its value and its extent;
     ! and what reads them, the physics and the grid, kept here so
     ! that a partial in a design is asked of the tower
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

  !===================================================================!
  ! An expansion lends pointers into its own storage, so a copy would
  ! share them.
  !===================================================================!

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

  !===================================================================!
  ! How many members the set a node denotes has, or zero where no
  ! extent was recorded.
  !===================================================================!

  integer function extent_of(this, g) result(n)

    class(expansion), intent(in) :: this
    type(graph)     , intent(in) :: g

    n = 0
    if (this % extents % describes(g)) n = this % extents % num_members_of(g)

  end function extent_of

  !===================================================================!
  ! The tuples of a coupling's relation, in the relation's own order,
  ! which is the order the coupling's value holds its weights in.
  ! Invalid input: a node that is not a coupling of one relation.
  !===================================================================!

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

  !===================================================================!
  ! The weights of a coupling in the order its relation holds the
  ! tuples: the relation groups them by source and keeps each once,
  ! so a weight computed per tuple as given is placed where the
  ! relation put its tuple. Invalid input: a tuple given twice, which
  ! would leave one weight with no place.
  !===================================================================!

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

  !===================================================================!
  ! THE BUILD.
  !===================================================================!

  subroutine build(this, physics, schemes, instants, steps, &
       & max_derivative_degree, parameter, nodes, spatial_discretization_stencil, weights, block_steps)

    class(expansion)      , intent(inout) :: this
    type(expression)      , intent(in)    :: physics
    type(family_holder)   , intent(in)    :: schemes(:)
    integer               , intent(in)    :: instants(:)
    class(grid)           , intent(in)    :: steps
    integer               , intent(in)    :: max_derivative_degree
    ! the physics' parameter, the first design
    real(dp)              , intent(in)    :: parameter
    ! given, every component holds one freedom per node, and the
    ! level below - a stencil over the nodes - is laid on every
    ! component the physics sits on
    integer      , intent(in), optional   :: nodes
    type(stencil), intent(in), optional   :: spatial_discretization_stencil
    ! given, the weights of the steps are designs too, read by the
    ! grid; and given, the steps of the horizon's instants are these
    ! rather than the grid's partition
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

    ! the designs, leaves of their own level, the last member of the root
    allocate(this % design_at(0), this % design_kind(0))
    call one_design(this, 'the physics'' parameter', [parameter], design_of_physics)
    if (present(weights)) call one_design(this, 'the weights of the steps', weights, design_of_steps)
    this % root_at = this % nodes % assemble([sweeps, this % nodes % assemble(this % design_at, 0)], 0)
    call this % labels % bind(this % node(this % root_at), &
         & 'expansion of ' // physics % name() // ' in the design')

  end subroutine build

  !===================================================================!
  ! One design as a leaf: its value, its extent, and what reads it.
  !===================================================================!

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

  !===================================================================!
  ! The designs read back: how many, what reads each, its extent and
  ! its value; the physics' parameter by itself; the rule the tower
  ! was built for.
  !===================================================================!

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

  !===================================================================!
  ! The steps' partials in the weights of the designed grid: the
  ! total derivative of every step along the weights listed, one
  ! variation per entry so that a repeated weight is a repeated
  ! derivative, exact from the grid's own derivative terms; and the
  ! first partials as a matrix, one column per weight.
  !===================================================================!

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

  !===================================================================!
  ! THE SPATIAL DISCRETIZATION STENCIL: one coupling over the nodes, shared by every
  ! component the physics sits on. Its carriers are the nodes read
  ! and the nodes whose rows are entered; its relation is the
  ! stencil's pattern, a node read into a node's row; its value the
  ! stencil's weights in the relation's order. Invalid input: a
  ! stencil over other than the nodes.
  !===================================================================!

  integer function spatial_discretization_coupling(this, spatial_discretization_stencil) result(at)

    class(expansion), intent(inout) :: this
    type(stencil)   , intent(in)    :: spatial_discretization_stencil

    integer , allocatable :: table(:,:), heads(:), tails(:), rows(:), cols(:)
    real(dp), allocatable :: given(:), w(:)
    integer :: read_nodes, entered_nodes, holder, e, ne

    if (spatial_discretization_stencil % pattern % num_vertices() /= this % node_extent) then
       error stop 'gti_expansion: the spatial discretization stencil is a stencil over the nodes'
    end if
    ! a stencil may name a pair of nodes more than once, its entries
    ! adding; a relation names a pair once, so the entries are combined
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

  !===================================================================!
  ! The steps, from the grid, over every instant of the horizon.
  !===================================================================!

  subroutine partition(steps, num_instants, design, dt)

    class(grid), intent(in) :: steps
    integer    , intent(in) :: num_instants
    real(dp)   , intent(in) :: design(:)
    real(dp), allocatable, intent(out) :: dt(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    class(field), allocatable :: out

    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    knobs    = stored_field('design', instants % vertex_set(), max(size(design), 1))
    call knobs % set_real_vector(padded(design))

    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)

  end subroutine partition

  pure function padded(design) result(x)

    real(dp), intent(in) :: design(:)
    real(dp), allocatable :: x(:)

    if (size(design) == 0) then
       allocate(x(1), source=0.0_dp)
    else
       x = design
    end if

  end function padded

  !===================================================================!
  ! A row attached and marked in one step.
  !===================================================================!

  subroutine attach_known(this, at, x)

    class(expansion), intent(inout) :: this
    integer         , intent(in)    :: at
    real(dp)        , intent(in)    :: x(:)

    call this % values % attach_unknown(this % node(at))
    call this % values % mark_known(this % node(at), padded(x))

  end subroutine attach_known

  !===================================================================!
  ! One sweep, over a horizon of its own. Sweep zero is the primal
  ! and the sweeps above it are the tangents; each owns its horizon,
  ! because the value rows are keyed on identity and one sweep's
  ! state must not stand for another's.
  !===================================================================!

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

  !===================================================================!
  ! One horizon: the blocks that partition the instants, in order.
  !===================================================================!

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

  !===================================================================!
  ! One block: its slices, the family that marches them, its steps,
  ! and the coupling that carries the scheme reach.
  !===================================================================!

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

  !===================================================================!
  ! One slice: the components of every degree at one instant.
  !===================================================================!

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

  !===================================================================!
  ! A slice of a family whose rows run between instants: the
  ! components of every degree, and no stage level at all.
  !===================================================================!

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

  !===================================================================!
  ! Whether a family's time discretization stencil rows run between the stages of one
  ! step rather than between instants. A stage family gives an empty
  ! pattern at every degree, which is how it says its rows are not
  ! offsets back through the instants; that question is asked here,
  ! so no family declares its kind twice.
  !===================================================================!

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

  !===================================================================!
  ! A slice of a stage family: the stages of the step arriving at
  ! this instant, then the instant itself, which is the step's
  ! closing evaluation and holds components like any stage. The
  ! block's first instant has no step arriving at it and so holds no
  ! stages; its components are the initial data.
  !===================================================================!

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

  !===================================================================!
  ! One stage, or the arriving instant when the index is zero: the
  ! components of every degree held at one evaluation point.
  !===================================================================!

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

  !===================================================================!
  ! One component: a leaf holding its freedoms, one per node, and on
  ! the physics' own degree of an evaluated moment the spatial discretization stencil.
  ! The instants a block reaches back over carry their values from
  ! the start; every component after them waits on a march.
  !===================================================================!

  integer function one_component(this, scheme, instant, first, degree, evaluated) result(at)

    class(expansion), intent(inout) :: this
    class(family)   , intent(in)    :: scheme
    integer         , intent(in)    :: instant, first, degree
    ! whether the physics is evaluated at this component's moment,
    ! where the spatial discretization stencil is laid on the physics' own degree
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

  !===================================================================!
  ! An integer and a real as text, for the labels.
  !===================================================================!

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

  !===================================================================!
  ! THE COUPLING OF A BLOCK.
  !
  ! How many edges the rows that fit inside this block make, and, on
  ! the second pass, what they are. A row on degree d at local
  ! instant kk fits when every source it reads lies at or after the
  ! block's first instant.
  !===================================================================!

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

  !===================================================================!
  ! The weights the reach carries, from the family and the steps of
  ! this block alone: a row that fits reads no instant before the
  ! block's first, so the block's own steps are all it needs.
  !===================================================================!

  subroutine reach_weights(scheme, n, tails, heads, source_degree, determines, dt, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: n, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: dt(:)
    real(dp), allocatable, intent(out) :: w(:)

    call weights_of(scheme_weight(scheme), n, tails, heads, dt, source_degree, determines, w)

  end subroutine reach_weights

  !===================================================================!
  ! The coupling itself: this block's slices as carriers, then the
  ! two sets the relation runs between, then the relation. The
  ! weights are the coupling's own value, so the sparsity and the
  ! numbers hang on one identity.
  !===================================================================!

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

  !===================================================================!
  ! Each edge as a tuple: the unknown its source is, and the unknown
  ! its constraint determines. Rows and columns share one index
  ! space, so the block is square.
  !===================================================================!

  pure function tuples(nd, tails, heads, source_degree, determines) result(table)

    integer, intent(in) :: nd, tails(:), heads(:), source_degree(:), determines(:)
    integer, allocatable :: table(:,:)

    integer :: e

    allocate(table(2, size(tails)))

    table(1,:) = [((tails(e) - 1) * nd + source_degree(e) + 1, e = 1, size(tails))]
    table(2,:) = [((heads(e) - 1) * nd + determines(e) + 1, e = 1, size(heads))]

  end function tuples

  !===================================================================!
  ! A carrier denoting a set of the given extent.
  !===================================================================!

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

  !===================================================================!
  ! Every level consistent and every coupling relationally valid.
  ! Each level is already checked as it is assembled; this checks the
  ! built hierarchy once more, from the outside.
  !===================================================================!

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

  !===================================================================!
  ! THE ROWS OF ONE STEP.
  !
  ! In the numbering a stage family uses, vertex one is the instant
  ! the step leaves from, vertices two to one plus s are the stages,
  ! and the last is the instant it arrives at. The rows made here are
  ! the ones inside the step, so no edge leaves vertex one: those
  ! carry the previous instant across and belong to the block.
  !
  !      degree d below the highest
  !          stage i reads stage j at degree d+1, for j at or before i
  !          the arriving instant reads every stage at degree d+1
  !      the highest degree
  !          the arriving instant reads every stage at that degree,
  !          which is the recovery
  !===================================================================!

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

  !===================================================================!
  ! A vertex of the family's numbering as an unknown of this step:
  ! stage i is member i, the arriving instant is member s+1.
  !===================================================================!

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

  !===================================================================!
  ! The weights of one step. The step is read at every vertex of the
  ! family's numbering, the stages included, because that is where
  ! the scaling reads it.
  !===================================================================!

  subroutine stage_weights(scheme, s, tails, heads, source_degree, determines, step, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: s, tails(:), heads(:), source_degree(:), determines(:)
    real(dp)     , intent(in) :: step
    real(dp), allocatable, intent(out) :: w(:)

    call weights_of(scheme_weight(scheme), s + 2, tails, heads, spread(step, 1, s + 2), &
         & source_degree, determines, w)

  end subroutine stage_weights

  !===================================================================!
  ! The coupling of one step: its stages and its arriving instant as
  ! carriers, then the two sets the tableau runs between, then the
  ! relation. The weights are the coupling's own value.
  !===================================================================!

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

  !===================================================================!
  ! THE CARRY BETWEEN STEPS.
  !
  ! Where a block's rows run inside its steps, what joins one step to
  ! the next is the instant they share: every row below the highest
  ! degree reads the previous instant's own degree, carried across
  ! unchanged. Those are the edges this coupling holds, and the
  ! family gives their weight rather than this module assuming it.
  !===================================================================!

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

!=====================================================================!
! packed from application/gti_block.f90
!=====================================================================!
!=====================================================================!
! The residual of one block, as one operation.
!
! A block's rows are of three kinds and they are added here into a
! single statement, because a minimizer drives one statement to zero:
!
!      derived     the scheme's own rows, already assembled as a
!                  stencil; linear in the state, so that stencil is
!                  also their jacobian
!      governing   the physics, at each evaluation point's primary
!                  degree - the one degree no time discretization stencil row determines
!      carried     the instants a block reaches back over, whose
!                  components are known before it starts; their rows
!                  are the identity less what they hold, so the block
!                  is square and nonsingular
!
! Its two arguments are the state and the design, in that order,
! which is what a minimizer supplies when the design is handed to it
! as a held input.
!
!             WHERE THE PHYSICS IS EVALUATED
!
! At the points given, and nowhere else. A multistep block evaluates
! at its instants, and its points are the instants in order, so the
! components it hands the physics are the state unchanged. A stage
! block evaluates at its stages and recovers its instants from them,
! so its points are the stages and the components are gathered out
! from between them. The physics is nodal either way and never learns
! which it is being asked about.
!
!             THE JACOBIAN
!
! Both halves carry exact partials - the stencil by being linear, the
! physics by differentiating its own rule - so the tangent is exact
! and nothing is differenced. A variation arrives named for this
! statement's argument and is renamed for each half before it is
! passed on, since each half checks the variation against its own,
! and a variation in the state is gathered to the points along with
! the state itself.
!
! A variation in the design is answered too, and it is a different
! statement: the scheme's rows are frozen at the steps they were
! built from and the carried rows hold given numbers, so neither
! varies with the design and only the governing rows do. That partial
! is what a sensitivity reads, by either the tangent or the adjoint.
!
!             WHAT IS REFUSED
!
! A state that is not one component per degree per unknown point; a
! missing argument; a carried row outside the unknowns; an evaluation
! point whose degrees run past them.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

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

  !-------------------------------------------------------------------!
  ! THE REACH a time discretization stencil row was built from, in the family's own
  ! numbering of vertices, kept so that the rows can be weighted again
  ! along a direction in the steps: which of the block's steps each
  ! vertex reads, each edge's tail and head vertex and degrees, and
  ! the block unknown each edge determines and reads at the first
  ! node - every other node lies a degrees' width further on.
  !-------------------------------------------------------------------!
  type :: coupling_reach
     integer :: vertices = 0
     integer, allocatable :: step_of(:)
     integer, allocatable :: tails(:), heads(:), source_degree(:), determines(:)
     integer, allocatable :: row(:), column(:)
  end type coupling_reach

  type, extends(operation) :: block_residual

     type(stencil)                       , private :: time_discretization_stencil
     type(expression)     , private :: physics

     ! THE SPATIAL DISCRETIZATION STENCIL. A stencil over the same unknowns coupling the
     ! components of one moment across the nodes of a spatial mesh:
     ! the spatial operator, linear in the state and independent of
     ! the design, laid on the block by spatial_discretization_laid. It adds to the
     ! time discretization stencil rows in the apply and in the tangent, and nowhere else,
     ! having no design partial and no partial above the first.
     ! Absent, the block is one node's.
     type(stencil), allocatable, private :: spatial_discretization_stencil
     type(stored_directed_graph)         , private :: points
     integer , allocatable               , private :: at(:)
     integer , allocatable               , private :: carried(:)
     real(dp), allocatable               , private :: held(:)
     integer                             , private :: degrees  = 0
     integer                             , private :: unknowns = 0
     integer                             , private :: primary  = 0

     ! WHERE THE BLOCK LIES in the graph: its node of the expansion,
     ! and the expansion that node belongs to. Where each unknown lies
     ! - the member of the time level, an instant or a step; its node
     ! of the space level; its moment, the instant or stage whose
     ! values it is among - is read from the graph whenever asked and
     ! held nowhere else: a sweep reads its members and their coupling
     ! from it, the spatial discretization stencil is laid on the moments, the aggregates
     ! a multigrid coarsens by are read off it. A member of a block
     ! keeps the block's node and the unknowns it was restricted to.
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

  !===================================================================!
  ! What the first instant a block was given holds. A solver starting
  ! from it begins near the trajectory rather than at nothing, which
  ! for a state of any size is much the same thing as starting at the
  ! wrong end of it.
  !===================================================================!

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

  !===================================================================!
  ! The stencil supplies one order and the physics many, so the
  ! statement supplies one.
  !===================================================================!

  pure integer function block_max_degree(this)

    class(block_residual), intent(in) :: this

    associate (u1 => this); end associate
    block_max_degree = 3

  end function block_max_degree

  !===================================================================!
  ! The components held at the evaluation points, taken out from
  ! among the unknowns so that a nodal rule reads them one point at a
  ! time.
  !===================================================================!

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

  !===================================================================!
  ! What the physics is handed: the gathered components and the
  ! design, both over the points.
  !===================================================================!

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

  !===================================================================!
  ! The governing value at each point, on the row its primary degree
  ! holds.
  !===================================================================!

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

  pure subroutine carry_direction(this, v, r)

    class(block_residual), intent(in)    :: this
    real(dp)             , intent(in)    :: v(:)
    real(dp)             , intent(inout) :: r(:)

    integer :: i

    do i = 1, size(this % carried)
       r(this % carried(i)) = v(this % carried(i))
    end do

  end subroutine carry_direction

  !===================================================================!
  ! A carried row holds a given number, which varies with nothing.
  !===================================================================!

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

  !===================================================================!
  ! The tangent: each half differentiated in its own argument, the
  ! variation renamed for it and, in the state, gathered to the
  ! points the physics reads.
  !===================================================================!

  !===================================================================!
  ! Where the unknowns lie: one time member and one space member per
  ! unknown. Labels of the wrong extent, or one below one, stop the
  ! program.
  !===================================================================!

  subroutine placed_on(this, tower, node)

    class(block_residual), intent(inout)     :: this
    type(expansion)      , intent(in), target :: tower
    type(graph)          , intent(in), target :: node

    this % tower => tower
    this % node  => node

  end subroutine placed_on

  !-------------------------------------------------------------------!
  ! Where every unknown lies, read from the block's node: the slices
  ! are its members; a slice whose first member is a leaf is one
  ! moment, any other slice's members are its moments; a moment holds
  ! one freedom per node of its first component's extent at every
  ! degree, node by node, degrees within a node; the moments lie one
  ! after another. A member of a block reads the block's and keeps
  ! the unknowns it was restricted to. Invalid input: a block placed
  ! nowhere.
  !-------------------------------------------------------------------!

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

  ! the largest node label: a member of a block keeps the block's
  ! numbering, so this is the extent a map over the nodes must reach
  integer function num_nodes(this)

    class(block_residual), intent(in) :: this

    integer, allocatable :: slice(:), node(:), moment(:)

    num_nodes = 1
    if (.not. associated(this % node)) return
    call this % labels_of(slice, node, moment)
    num_nodes = maxval(node)

  end function num_nodes

  !-------------------------------------------------------------------!
  ! THE SPATIAL DISCRETIZATION STENCIL, laid on this block: a stencil over the nodes is
  ! placed at every moment the block evaluates its physics at, on the
  ! row the physics sits on, and reads the values of that moment. A
  ! moment with no evaluation point - an instant a stage block
  ! recovers - takes no spatial rows, since no physics is stated
  ! there. Invalid input: a stencil over other than the nodes; a
  ! stencil carrying a constant, which would be a source the block
  ! has no place for; a moment holding some nodes and not others.
  !-------------------------------------------------------------------!

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

    ! where each node's components lie at each moment with a point
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

  !-------------------------------------------------------------------!
  ! The reach the time discretization stencil rows were built from, given to the block.
  !-------------------------------------------------------------------!

  subroutine with_reach(this, reach)

    class(block_residual), intent(inout) :: this
    type(coupling_reach) , intent(in)    :: reach(:)

    this % reach = reach

  end subroutine with_reach

  !-------------------------------------------------------------------!
  ! The time discretization stencil rows' weights and every total
  ! derivative of the weights along subsets of n directions in the block's
  ! steps, the steps seeded subset by subset (seeds(k, m) the total
  ! derivative of step k along the subset with mask m): one triple
  ! per node per edge of every reach, the row and column among the
  ! unknowns, and w(:, m) the weight's total derivative along mask m,
  ! the weight at m = 0, signed as the rows are. The family's
  ! weight action carries the partials in the steps, and the
  ! determined component, entering with one, takes no part. Invalid
  ! input: a block built without its reach.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The aggregates a multigrid coarsens this block by: the coarse
  ! cell of each unknown's node, at its own moment and degree, and
  ! numbered from one in the order met. Invalid input: a map that
  ! does not reach every node label.
  !-------------------------------------------------------------------!

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

  !===================================================================!
  ! THE LINEAR BLOCK: the tangent in the state at the inputs given,
  ! frozen, as a block of its own, so that a linear statement A w = b
  ! - or A^T w = b - goes through the same solve as the block it came
  ! from and converges in one newton step. Its derived stencil is A
  ! with -b as its constant, its physics is zero, its points and
  ! degrees are this block's, and it carries no rows: the identities
  ! on the carried components are already in A. A negative stamp is
  ! given to the transpose, which a direct solver reads as the same
  ! factors the other way round.
  !===================================================================!

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

  !===================================================================!
  ! THE BLOCK RESTRICTED to a member of a level - some of its unknowns
  ! - with the rest held at the values given. The derived and the
  ! spatial rows restrict as stencils do, the outside taken into
  ! their constants; the points whose components all lie inside stay
  ! points; the carried rows inside stay carried. A point half inside
  ! stops the program, a member being whole points or nothing.
  !===================================================================!

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
    ! the member lies where the block lies, on the unknowns kept
    sub % tower => this % tower
    sub % node  => this % node
    if (allocated(this % kept)) then
       sub % kept = this % kept(kept)
    else
       sub % kept = kept
    end if

  end function block_restricted

  !===================================================================!
  ! THE COMPILED TANGENT in the state. The block knows its own
  ! structure: the time discretization stencil rows and the spatial rows are stencils
  ! already, and the physics is nodal, so its tangent at every point
  ! comes from one partial action per degree - a direction of one on
  ! that degree at every point at once, the points being independent.
  ! A carried row is an identity. Triples landing on one entry are
  ! combined. Only the state's tangent is compiled; any other
  ! argument is not available.
  !===================================================================!

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
    integer :: e, d, p, ne, npts, n, kept, count

    available = which == 1
    if (.not. available) return

    n    = this % unknowns
    npts = size(this % at)
    allocate(is_carried(n), source=.false.)
    is_carried(this % carried) = .true.

    call state_of(this, input_data, input_graph, x, state)
    call point_inputs(this, input_data, x, inputs)

    ! room for the derived and spatial triples, the physics's degrees
    ! per point, and the carried identities
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

    associate (u1 => ne); end associate

  end subroutine block_compiled_tangent

  !-------------------------------------------------------------------!
  ! A stencil's triples, less those on carried rows, appended.
  !-------------------------------------------------------------------!

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

    ! THE SECOND AND THIRD PARTIALS. The time discretization stencil
    ! rows, the spatial discretization stencil and the carried rows
    ! are linear in the state and read no design, so only the physics
    ! has one: the partial at every point along every direction given,
    ! on the row the equation is imposed on.
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

  !===================================================================!
  ! The physics' second partial along two variations, each on the
  ! argument its variation names: a state direction is gathered over
  ! the points, a design direction is read as given.
  !===================================================================!

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

  !===================================================================!
  ! The design half: the scheme's rows are frozen and the carried
  ! rows hold given numbers, so only the governing rows vary.
  !===================================================================!

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

  !===================================================================!
  ! THE ORDER A LEVEL IS SWEPT IN, derived and not declared. Swept by
  ! instants, the members are coupled by the time discretization stencil rows - a row at
  ! one instant reading a point at another - and that coupling is
  ! acyclic for any march, since every scheme reads backward; its
  ! loop is the sweep. The transposed block's pattern is the same
  ! graph read the other way, so it sweeps from the last instant by
  ! the same rule and no one says so. Swept by nodes, the
  ! coupling is the spatial discretization stencil's and symmetric, so no node is before
  ! another and they are swept as they lie.
  !===================================================================!

  subroutine member_order(this, by_instants, order)

    class(block_residual), intent(in)  :: this
    logical              , intent(in)  :: by_instants
    integer, allocatable , intent(out) :: order(:)

    integer, allocatable :: table(:,:), label(:), slice(:), node(:), moment(:)
    type(stored_directed_graph) :: coupling
    integer :: ne, e, n, t, h, k, members

    call this % labels_of(slice, node, moment)

    ! the space level's members couple both ways through the mesh,
    ! which has no loop, and are swept as numbered
    if (.not. by_instants) then
       order = [(k, k = 1, maxval(node))]
       return
    end if

    ! the time level's coupling: a time discretization stencil row at one member that
    ! reads an unknown at another, which for every family looks one way
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

!=====================================================================!
! packed from application/gti_march.f90
!=====================================================================!
!=====================================================================!
! Building one block and solving it.
!
! The pieces are the same ones the assembly uses - the rows a family
! reaches over, the weights on them, the stencil they make, and the
! statement that adds the governing and carried rows to it - gathered
! here so that a caller marching a block and a caller differentiating
! one write them once.
!
!             WHERE THE BLOCKS OF A HORIZON SIT
!
! horizon_bounds says which instants each block of a chain spans. A
! block reaches back over instants that begin before it does, so
! every block after the first overlaps what came before it by
! exactly what its family reaches. What is done with that overlap -
! the junction, and the layouts either side of it - belongs to
! gti_chain, which marches them.
!
! A block must add more instants than its family reaches back over,
! or it would consist of nothing but what it was given.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

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
  use operation_action      , only : variation
  use operation_stencil       , only : stencil
  use operation_newton        , only : newton
  use operation_minimization  , only : minimizer, relative, absolute, &
       & by_count, by_rate
  use operation_dense_direct  , only : dense_direct
  use operation_gmres         , only : gmres
  use operation_family        , only : family
  use operation_grid          , only : grid, uniform_grid
  use operation_weight        , only : scheme_weight
  use operation_scheme_stencil, only : derived_constraints
  use operation_expression       , only : expression
  use gti_expansion           , only : family_holder, expansion, marches_by_stages
  use gti_block               , only : block_residual, coupling_reach
  use view_level              , only : level_member, level_num_members, level_coupling, &
       & level_couples
  use graph_fractal           , only : graph
  use map_value               , only : VALUE_KNOWN
  use gti_sweeps              , only : jacobian_of, assembly_present, multigrid_on, &
       & set_aggregates, coarse_nodes, take_inner, keep_inner, forget_inner, set_linear_stopping
  use util_tally              , only : tally_record, tangent_loops, adjoint_loops

  implicit none

  !-------------------------------------------------------------------!
  ! WHAT AN UNCONVERGED MARCH LEFT, by aspect. The imbalance is a
  ! vector with one entry per unknown, and the unknowns lie in slots
  ! of one instant's - or one stage's - components each, so its norm
  ! splits exactly, ||r||^2 = sum over slots and degrees of r^2, and
  ! its steepest direction in the state is the gradient of the norm,
  ! d||r||/dq = A^T r / ||r||, one transposed matvec. The largest entry
  ! of each names where the imbalance sits and which state drives it.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! WHICH LEVEL IS SWEPT. The block is one nonlinear statement over
  ! every instant, node and component; solving it is a choice of
  ! which level's members are solved exactly inside and which level
  ! is swept over them with the rest held:
  !
  !      space-time    no level: the whole block at once
  !      time          the instants, in order: each instant's nodes
  !                    and components solved with the instants before
  !                    it held - the classical step, exact in one pass
  !                    since every scheme looks backward
  !      space         the nodes: each node's whole history solved
  !                    with its neighbours' histories held, the sweep
  !                    repeated until the coupling agrees
  !
  ! The same fixed point in all three. Nothing below this loop knows
  ! which was chosen.
  !-------------------------------------------------------------------!

  character(len=16), save :: sweep_level = 'space-time'

  ! Stamps handed out to statements, so that a direct solver can tell
  ! a statement it has factorised from a new one.
  integer, save :: stamps_given = 0

  !-------------------------------------------------------------------!
  ! HOW A MARCH STOPS. A caller that sets nothing gets a tolerance
  ! measured against the imbalance the march began at, and a budget
  ! taken from the rate the march itself shows. The count is a
  ! backstop and not the operative limit.
  !-------------------------------------------------------------------!

  real(dp), save :: stopping_tolerance  = 1.0e-12_dp
  integer , save :: stopping_criterion  = relative
  integer , save :: stopping_budget     = by_rate
  integer , save :: stopping_iterations = 100

  private
  public :: partition, partitioned, solved, unknowns_graph
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

  !===================================================================!
  ! THE WEIGHT A BLOCK CARRIES: ||A||_inf from the family and the step,
  ! without the matrix. A row determining degree d reads the sources
  ! its pattern names, each weighted alpha dt^(sigma - d) by the same
  ! scheme_weight that builds the block, so the row's absolute sum is
  ! one apply on a coupling of that pattern at the step given, plus
  ! the one the row carries on the column it determines. The largest
  ! over the degrees is the norm. A stage family has no pattern in
  ! instants and its rows lie within one step: the incoming instant
  ! and the stages at or before, read the same way.
  !===================================================================!

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

  !===================================================================!
  ! THE PRECISION A TARGET NEEDS. The floor a march reaches is
  ! eps ||A|| ||q||, so a target is reachable at a kind whose spacing
  ! is under target / (||A|| ||q||). The target is the tolerance times
  ! the starting imbalance where the criterion is relative, and the
  ! tolerance itself where it is absolute.
  !===================================================================!

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

  !===================================================================!
  ! THE CONSISTENT INITIAL STATE. Given the components below the
  ! highest at one instant, the highest is what the physics says it
  ! is there: q^(N) with R(q, q', ..., q^(N)) = 0, solved at that one
  ! instant with everything below it held.
  !
  ! This is the smallest block there is - one evaluation point, no
  ! scheme rows, the lower components carried and the highest the one
  ! unknown - and it is solved by the same newton as every other
  ! block. Nothing about the physics is assumed: whatever R is, its
  ! zero at the instant is what comes back. A lower vector of the
  ! wrong extent stops the program.
  !===================================================================!

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

  !===================================================================!
  ! The state at the first instant, consistent with the physics: the
  ! components below the highest are given at every node, and the
  ! highest is what the physics then requires, with the spatial discretization stencil
  ! - laid on the given values - entering its row. No block is laid
  ! for it: the physics is a rule at one point, so the highest
  ! component solves node by node, and the physics' partial in it,
  ! which the rule carries, is the slope. Invalid input: components
  ! for other than every degree below the highest; a spatial discretization stencil over
  ! other than the nodes.
  !===================================================================!

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
    integer  :: nodes, i, d, k, top, iteration

    nodes = size(lower, 2)
    top   = degrees - 1
    if (size(lower, 1) /= top) then
       error stop 'gti_march: the components below the highest are given at every node'
    end if

    ! the spatial discretization stencil on the given values: what it adds to each
    ! node's row of the highest degree
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

    ! the state over the nodes as points, the highest component from
    ! the rule's own linear part
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

    ! newton on the highest component, node by node at once: the
    ! residual and its partial in that component at every node
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
    associate (u1 => d); end associate

  end function consistent_states

  !===================================================================!
  ! A state and a design as the two inputs a block's rows read: the
  ! state on its unknowns, one design value per point.
  !===================================================================!

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

  !===================================================================!
  ! How every march that follows stops. A criterion or a budget that
  ! is neither of its two stops the program.
  !===================================================================!

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
    ! the inner solves honour the same tolerance, criterion and budget
    call set_linear_stopping(tolerance, criterion, budget)

  end subroutine set_stopping

  !===================================================================!
  ! A uniform partition of the duration, and the instants it makes.
  !===================================================================!

  subroutine partition(duration, n, dt, t)

    real(dp), intent(in) :: duration
    integer , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)

    call partitioned(uniform_grid(duration), n, dt, t)

  end subroutine partition

  !===================================================================!
  ! The instants a grid makes over the duration it was given.
  !===================================================================!

  subroutine partitioned(steps, n, dt, t, design)

    class(grid), intent(in) :: steps
    integer    , intent(in) :: n
    real(dp), allocatable, intent(out) :: dt(:), t(:)
    real(dp), intent(in), optional :: design(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs
    class(field), allocatable :: out
    integer :: k

    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])

    if (present(design)) then
       knobs = stored_field('design', instants % vertex_set(), size(design))
       call knobs % set_real_vector(design)
    else
       knobs = stored_field('design', instants % vertex_set(), 1)
       call knobs % set_real_vector([0.0_dp])
    end if

    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)

    allocate(t(n))
    t(1) = 0.0_dp
    do k = 2, n
       t(k) = t(k - 1) + dt(k)
    end do

  end subroutine partitioned



  !===================================================================!
  ! Where a component lies: instants follow one another, nodes lie
  ! within an instant, and the components of one point stay together,
  !
  !      ((instant - 1) nodes + (node - 1)) degrees + degree + 1
  !
  ! which at one node is (instant - 1) degrees + degree + 1, the
  ! ordinary block. A field over a mesh is this at nodes > 1 and
  ! nothing else.
  !===================================================================!

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
  ! THE BLOCK FROM ITS GRAPH NODE. Block b of the horizon of an
  ! expansion, read under the level view: the slices are the block's
  ! members; each slice's components - or, for a stage family, its
  ! stages and arriving instant, each with components - are the
  ! moments, laid one after another, each a width of nodes times
  ! degrees, node by node within a moment; the couplings' relations
  ! give the time discretization stencil rows, their weights recomputed from the family in
  ! the relations' own tuple order; a component the graph holds as
  ! known is carried; the physics is evaluated at every moment of a
  ! difference family and at the stages of a stage family. The reach
  ! is kept on the block, the spatial discretization stencil laid on the moments. The
  ! block's steps are the graph node's own value. Invalid input: a
  ! held value for other than every carried component.
  !===================================================================!

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

    ! the nodes a component holds, read from the first component
    m     = tower % extent_of(level_member(first_moment(block, staged), 1))
    width = nd * m

    ! the moments in order: which slice each lies in, and which member
    ! of it; a difference family's slice is one moment, a stage
    ! family's the stages then the arriving instant, the first slice
    ! the instant alone
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

    ! the carried components: known in the graph, at every node; and
    ! the spatial discretization stencil, one coupling over the nodes shared by every
    ! evaluated moment's physics component, read once
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
       ! carried at every node, node by node, degrees within a node
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

    ! the reach: one coupling over the instants for a difference
    ! family; for a stage family one per step, its stages' tableau
    ! and the carry from the instant before
    if (staged) then
       call stage_reach_of(tower, block, scheme, n, s, nd, width, moments, slice_of, &
            & member_of, reach)
    else
       call block_reach_of(tower, block, n, nd, width, reach)
    end if

    ! the time discretization stencil rows: every coupling's edges weighted by the family
    ! at the block's steps, replicated per node
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

    ! the block lies at its node of the graph, which says where every
    ! unknown lies
    call rows % placed_on(tower, block)
    call rows % with_reach(reach)

    ! the spatial discretization stencil, laid on every moment the physics sits at: the
    ! coupling's relation is the stencil's pattern, a node read into
    ! a node's row, its value the weights in that order
    if (associated(below)) then
       call tower % tuples_of(below, table)
       call tower % value_of(below, spatial_weights)
       call rows % spatial_discretization_laid(stencil(table(2, :), table(1, :), spatial_weights, &
            & spread(0.0_dp, 1, m), 'spatial discretization stencil'))
    end if

    ! where each instant lies: a slice's last moment
    allocate(instants_at(n))
    g = 0
    do k = 1, n
       g = g + members(k)
       instants_at(k) = (g - 1) * width
    end do

  end subroutine block_from

  !===================================================================!
  ! The first moment of a block: its first slice for a difference
  ! family, that slice's one member for a stage family.
  !===================================================================!

  function first_moment(block, staged) result(moment)

    type(graph), intent(in) :: block
    logical    , intent(in) :: staged
    type(graph), pointer :: moment

    moment => level_member(block, 1)
    if (staged) moment => level_member(moment, 1)

  end function first_moment

  !===================================================================!
  ! The reach of a difference family's block: its coupling's tuples,
  ! each a source component and the component it determines in the
  ! slice-major numbering with one node, read back into instants and
  ! degrees; every vertex reads its own step.
  !===================================================================!

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

  !===================================================================!
  ! The reach of a stage family's block: one coupling per step, its
  ! vertices numbered as the family numbers them - the instant the
  ! step leaves from, its stages, the instant it arrives at - every
  ! vertex taking the step's own size. The step's own coupling on
  ! its slice gives the tableau's edges in the step's numbering of
  ! members; the block's coupling gives the carry from the instant
  ! before, in the block's numbering with one node.
  !===================================================================!

  subroutine stage_reach_of(tower, block, scheme, n, s, nd, width, moments, slice_of, &
       & member_of, reach)

    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: block
    class(family)  , intent(in) :: scheme
    integer        , intent(in) :: n, s, nd, width, moments, slice_of(:), member_of(:)
    type(coupling_reach), allocatable, intent(out) :: reach(:)

    integer, allocatable :: table(:,:), carry(:,:), first_moment(:), counted(:), filled(:)
    integer :: kk, e, g, tail_moment, head_moment, vertex_tail, vertex_head

    associate (u1 => scheme); end associate
    allocate(reach(n - 1), first_moment(n), counted(n), filled(n))
    first_moment(1) = 1
    do kk = 2, n
       first_moment(kk) = first_moment(kk - 1) + merge(1, s + 1, kk - 1 == 1)
    end do

    ! count the edges of each step: the tableau's plus the carry's
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

    ! the tableau's edges: member j of the step is vertex j + 1
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

    ! the carry's edges: the instant before is vertex one
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
    associate (u2 => moments); end associate

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

  !===================================================================!
  ! Newton over the whole block. The design is held while the state
  ! varies, which is what a minimizer supplies as an extra input.
  !
  !             WHAT COUNTS AS SOLVED
  !
  ! A scheme's rows carry a power of the step, so a difference on the
  ! second derivative weighs its sources by the inverse square of it.
  ! Refining the grid therefore raises the size of a residual for the
  ! same trajectory, and the smallest one reachable in the arithmetic
  ! rises with it: at a hundredth of a unit it is near ten to the
  ! minus thirteen, and finer than that it passes any fixed target.
  !
  ! Asked for a fixed one, newton reaches the trajectory in two steps
  ! and then spends its whole budget failing to better it. Measured on a
  ! degree-two problem over three units: a hundred and twenty instants
  ! took a hundred and sixty seconds to produce what forty iterations
  ! produce in a sixth of one, to the same six digits.
  !
  ! So the target is set against the residual the first guess gives,
  ! which is the only scale in the problem that is known before it is
  ! solved, and the budget is a backstop rather than a cost.
  !===================================================================!

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

    ! every unknown lies in a point of degrees consecutive components,
    ! a stage's as much as an instant's, and a point is smoothed whole
    width = rows % num_degrees()
    ! multigrid coarsens by aggregates read off the block: the coarse
    ! cell of each unknown's node, at its own moment and degree
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

  !===================================================================!
  ! A stamp no statement has had before.
  !===================================================================!

  integer function fresh_stamp() result(mark)

    stamps_given = stamps_given + 1
    mark = stamps_given

  end function fresh_stamp

  !===================================================================!
  ! A LINEAR SYSTEM IN THE TANGENT, A w = rhs or A^T w = rhs, solved
  ! as the block it came from is solved: as a linear block through the
  ! sweep, where newton stops after one step. The stamp given is the
  ! tangent's; every right side against the same tangent gives the
  ! same stamp, and a direct solver then factorises once.
  !===================================================================!

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

  !===================================================================!
  ! The gradient in the design by the tangent - one solve in the
  ! state, the gradient read along it - and by the adjoint - one solve
  ! against the transpose, the design partial read along it. Both
  ! through the sweep, against one tangent, stamped once.
  !===================================================================!

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

  !===================================================================!
  ! Which level the marches sweep. A name that is none of the three
  ! stops the program.
  !===================================================================!

  subroutine set_sweep(name)

    character(len=*), intent(in) :: name

    call refuse_unknown(name, ['space-time', 'time      ', 'space     '], 'sweep')
    sweep_level = name

  end subroutine set_sweep

  pure function sweep_named() result(name)

    character(len=:), allocatable :: name

    name = trim(sweep_level)

  end function sweep_named

  !===================================================================!
  ! THE SWEEP. The block's points lie instant by instant, nodes within
  ! an instant, so a member of the time level is one instant's points
  ! and a member of the space level is one node's points across the
  ! instants. Each member is solved as a block of its own, restricted
  ! from the whole with the rest held at the current state, and its
  ! solution written back; a member with nothing to solve - every
  ! component carried - is passed over. A pass is judged on the whole
  ! block's residual by the same criteria as any iteration, so the
  ! time sweep, exact after one pass, stops on the second, and the
  ! space sweep stops where the coupling has settled.
  !===================================================================!

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

    ! the level's members, read off the block's own labels: instants
    ! or steps for the time level, nodes for the space level
    if (trim(sweep_level) == 'time') then
       label = rows % slice_of()
    else
       label = rows % node_of()
    end if
    members = maxval(label)



    allocate(is_carried(count), source=.false.)
    is_carried(rows % carried_unknowns()) = .true.

    ! the seed, and the carried components at what they are held at:
    ! a member that is all carried is then already solved
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


    ! the order the members are swept in is the coupling's own: a
    ! transposed statement, upper triangular in time, sweeps from the
    ! last instant because its pattern says so
    call rows % member_order(trim(sweep_level) == 'time', order)


    ! The residual where the sweep begins is what a relative target
    ! is measured against, as the first residual is for any march.
    achieved = whole_residual(rows, unknowns, design, q)
    call judge % note_imbalance(achieved)

    do pass = 1, stopping_iterations

       before = achieved

       do mm = 1, members
          m = order(mm)
          member = pack([(k, k = 1, count)], label == m)
          if (all(is_carried(member))) cycle

          ! a member not yet solved is seeded from the one before it in
          ! the order swept, which is continuation, the seed every step
          ! of a march has: a member of the same extent is copied, and
          ! a step's stages and arriving instant each take the instant
          ! before them
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

       ! A pass that left the residual exactly where it was has
       ! reached the sweep's fixed point; another would do the same.
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

  !-------------------------------------------------------------------!
  ! The unknowns of the m-th member: the points of instant m, or the
  ! points of node m across the instants, each point's components.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The seed of a member from the member solved before it. Of the
  ! same extent, the values are copied; otherwise the last point of
  ! the earlier member - the instant a step arrives at - is laid on
  ! every point of the later one, which is where a step's stages and
  ! its own arriving instant begin.
  !-------------------------------------------------------------------!

  subroutine continued(q, member, earlier)

    real(dp), intent(inout) :: q(:)
    integer , intent(in)    :: member(:), earlier(:)

    integer :: pieces, i, width

    ! a member of the same extent is copied; a step's stages and
    ! arriving instant each take the instant before them; an instant
    ! after a step takes the step's last piece
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

  !===================================================================!
  ! The aspects of what was left: the norm split by degree, the
  ! largest entry, and the largest entry of A^T r / ||r||. A is formed
  ! here in full, which is O(n^2) and is paid only on a march that
  ! did not converge.
  !===================================================================!

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

  !===================================================================!
  ! How large a residual the first guess gives, which is the scale
  ! the target is set against.
  !===================================================================!


  !===================================================================!
  ! A first guess: every point of the block holding what its first
  ! instant was given. It costs nothing to form and it starts newton
  ! near the trajectory rather than at zero, which for a state of any
  ! size is far away and is where a jacobian is most likely to be
  ! singular.
  !===================================================================!

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

  !===================================================================!
  ! The linear solver inside newton. A dense factorisation forms the
  ! jacobian column by column - one application of the statement per
  ! unknown - and then costs the cube of the count to factor, so it
  ! wins while the count is small and loses badly once it is not. The
  ! statement supplies a matvec through its partial action, so a
  ! krylov solver forms no matrix at all.
  !
  ! Where the crossing sits, and why it is settable, is stated in
  ! gti_sweeps, which owns it.
  !===================================================================!


  !===================================================================!
  ! Where each block begins and ends. A block adds the instants given
  ! for it and reaches back over its predecessor's last, so the
  ! blocks overlap by exactly what each family reaches.
  !===================================================================!

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

!=====================================================================!
! packed from application/gti_adaptive.f90
!=====================================================================!
!=====================================================================!
! Adaptive time stepping, phase one: an error-controlled forward march
! that discovers a grid. The grid it returns is an ordinary partition
! of the duration; the chain runs on it unchanged, so its sensitivities
! are formed exactly as on a fixed grid. The march itself takes no
! derivatives.
!
!             THE ERROR ESTIMATE
!
! Step doubling: one step of size h against two of h/2, both from the
! same state with the same scheme. Their difference on the solution
! components below the highest estimates the local error of the h step,
! which for an order-p scheme is O(h^(p+1)). A step is accepted when the
! estimate, measured relative to the state or absolute as a solve's
! tolerance is, is at or below the tolerance; the next step is
! h (tolerance / estimate)^(1/(p+1)), bounded so one step neither grows
! nor shrinks without limit, and clamped so the last lands on the
! duration exactly. The state advanced is the two-half-step one.
!
! The families marched are the diagonally implicit ones, self-starting
! and one step at a time, so a step's coefficients do not depend on the
! steps around it. A variable-step multistep family, whose coefficients
! do, is a separate construction and is refused here by its history
! reaching past one instant.
!
!             WHAT IS REFUSED
!
! A scheme that reaches back over more than one instant; a nonpositive
! duration, tolerance or first step; a step that stays above the
! tolerance past fifty attempts, where the estimate is not falling.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

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

  !-------------------------------------------------------------------!
  ! One step of size h from a state by n - 1 uniform steps: n = 2 is
  ! the single step, n = 3 the two half steps. A block over n instants
  ! of a uniform grid of extent h has n - 1 steps of h / (n - 1); the
  ! arriving instant's components come back.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The estimate over the solution components below the highest, the
  ! highest being algebraically determined: relative to the state, or
  ! absolute.
  !-------------------------------------------------------------------!

  pure real(dp) function estimate(coarse, fine, degrees, relative) result(e)

    real(dp), intent(in) :: coarse(:), fine(:)
    integer , intent(in) :: degrees
    logical , intent(in) :: relative

    e = norm2(coarse(1:degrees - 1) - fine(1:degrees - 1))
    if (relative) e = e / max(norm2(fine(1:degrees - 1)), tiny(1.0_dp))

  end function estimate

  !-------------------------------------------------------------------!
  ! The accepted steps of a scheme of order p over [0, duration] to a
  ! tolerance, from the state consistent with the lower components. The
  ! first step is a fraction of the duration; the controller is bounded
  ! and the last step clamped to the duration.
  !-------------------------------------------------------------------!

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

!=====================================================================!
! packed from application/gti_space.f90
!=====================================================================!
! THE SPATIAL LEVEL: the framework's mesh over a two-dimensional
! domain, and the framework's operator on it.
!
! The mesh is structured in two parametric coordinates on the unit
! square and mapped to the domain by its geometry:
!
!      cartesian     x = a xi,          y = b eta
!      circular      x = a xi cos 2 pi eta, y = a xi sin 2 pi eta
!      elliptical    x = a xi cos 2 pi eta, y = b xi sin 2 pi eta
!
! so one indexing serves every shape and the shape is a mapping. The
! spacing along each coordinate IS the time grid's, uniform or drawn
! from a seed through the same partition, so a seed means the same
! thing in space as in time - and every cell is the
! polygon of its mapped corners: its area, centroid and face geometry
! are read off those corners, whatever the shape. A polar mapping
! collapses the inner corners onto the origin, so the innermost ring
! is one cell, a polygon of the first ring's corners.
!
! What is built is the mesh view_mesh already defines - cells as
! vertices, faces as edges, a boundary face an edge without a head
! and tagged - and the operator is operation_diffusion's: the fitted
! balance, a polynomial form of the order asked for fitted on each
! face's neighbourhood and aimed along the face normal, so a skewed
! face is measured in the right direction and the order is the
! form's. The outer boundary holds no flux. Nothing spatial is
! discretised here; this module maps a shape and hands the mesh on.
!
! The operator is the flux balance, integrated over each cell; a
! caller wanting the laplacian divides each row by its cell's area.
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
  use operation_grid            , only : uniform_grid, random_grid
  use gti_configuration         , only : chosen_from
  use gti_march                 , only : partitioned

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

     ! the framework's mesh
     type(mesh) :: m

     ! corners, and each cell's corners in order, ragged - kept for
     ! writing the cells out as polygons
     real(dp), allocatable :: corner(:,:)
     integer , allocatable :: first_corner(:)
     integer , allocatable :: cell_corner(:)

     real(dp), allocatable :: centre(:,:)
     real(dp), allocatable :: volume(:)

     ! each cell's place on the parametric grid, for coarsening
     integer , allocatable :: cell_ij(:,:)
     integer :: n1 = 0, n2 = 0

  end type room

  ! a face while the mesh is being built: its cells, its corners
  type :: face_record
     integer :: tail = 0, head = 0, corner_a = 0, corner_b = 0
  end type face_record

contains

  !-------------------------------------------------------------------!
  ! The geometry a name denotes; an unknown name stops the program.
  !-------------------------------------------------------------------!

  integer function geometry_of(name) result(geometry)

    character(len=*), intent(in) :: name

    ! the constants are the places in this list
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

  !-------------------------------------------------------------------!
  ! The mesh: counts along the two coordinates, the extents a and b,
  ! the spacing, and the geometry. A count below two along either
  ! coordinate, or an extent that is not positive, stops the program.
  !-------------------------------------------------------------------!

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

    ! the spacing along each coordinate is the time grid's own draw,
    ! scaled to the unit interval; the second coordinate continues
    ! the draw past the first, as a seed offset by the first's count
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

    !----------------------------------------------------------------!
    ! Faces. Interior ones along the first coordinate between rings
    ! or columns and along the second between neighbours, periodic
    ! in the second for a polar mapping; boundary ones on the outer
    ! edge of the domain, headless.
    !----------------------------------------------------------------!

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

    ! a polar mesh's first row of corners is its innermost ring's;
    ! the row below it, every point the origin, is never made
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

  !-------------------------------------------------------------------!
  ! The mesh from the corners, the cells and the faces enumerated
  ! above, through the framework's one ending of every mesh pipeline:
  ! areas, centroids, normals, deltas and weights are its, computed
  ! as for a mesh read from a file. The outer boundary is the wall.
  ! The cell centres and areas the level reads are then the mesh's.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The spatial operator: the diffusion statement on the mesh, the
  ! conductivity kappa through every face, no flux at the wall, and
  ! the polynomial form of the degree asked for, each fit taking as
  ! many rings as its members need. What comes back is the flux
  ! balance per cell, integrated.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The cells coarsened by pairs along each parametric coordinate:
  ! one aggregate per block of two by two, the innermost polar cell
  ! its own. The aggregate of every cell, numbered from one.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! One instant written for paraview by the framework's writer: the
  ! cells as polygons in the order their corners lie, one scalar per
  ! cell for each name. A numbered series of these is read as steps
  ! in time.
  !-------------------------------------------------------------------!

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

!=====================================================================!
! packed from application/gti_field.f90
!=====================================================================!
! THE FIELD: what a field adds to a march over a spatial mesh, and
! what is checked about it.
!
! A block over a field holds every node's components at every moment,
! node by node within a moment, laid out by the constructors in
! gti_march and gti_stage exactly as a single node's block is with
! nodes = 1. What the field adds is the spatial discretization stencil: a stencil over
! the nodes carrying minus the framework's diffusion operator, each
! row divided by its cell's area so that the flux balance becomes
! kappa times the laplacian. A block lays it on every moment it
! evaluates its physics at, where it adds to the physics' row. It is
! linear in the state and knows nothing of the design, so it enters
! the apply and the tangent and nothing else. Its order is the form's,
! and the form's degree is given.
!
! The functional over a field is the integral over the domain and the
! duration, so the measure of a node is its cell's area, and the
! chain weights each point by the step times that area.
!
! Three things are checked, each against something known:
!
!   the operator  the balance of the rectangle's mode against kappa
!                 times its laplacian, cell by cell, walls apart
!   the mode      at nu = 0 on a rectangle the field q = cos(pi x / a)
!                 cos(pi y / b) cos(omega t) is exact, with omega^2 =
!                 1 + kappa pi^2 (1/a^2 + 1/b^2), so the last instant
!                 is measured against it, and against the semi-discrete
!                 solution that isolates the time error
!   the ode       at kappa = 0 with a constant field every node is one
!                 node's equation - checked by the program, which
!                 marches the node
!
! and any instant may be written as a vtu file for paraview.
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

  !-------------------------------------------------------------------!
  ! The spatial discretization stencil as a stencil over the nodes: -kappa times the
  ! laplacian, one row per cell. A wall holding a value would enter
  ! as a source, which a block has no place for; the wall here holds
  ! no flux, and the block refuses a constant if one arrives.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The state at the first instant: the components below the highest
  ! at every node - constant from the words given, or the rectangle's
  ! mode, or one plus half of it - and the highest solved from the
  ! physics with the spatial discretization stencil laid on, so that the state is
  ! consistent with the equation rather than merely plausible. One
  ! node with no mesh is one node's equation. Invalid input: more
  ! components than lie below the highest; a mode with no rectangle.
  !-------------------------------------------------------------------!

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
       ! one plus half the rectangle's mode, on any geometry: a field
       ! that is not uniform, so the spatial discretization stencil has something to do
       if (.not. present(space)) error stop 'gti_field: the bump is a field over a mesh'
       lower(1, :) = 1.0_dp + 0.5_dp * mode_shape(space, a, b)
    case default
       error stop 'gti_field: an initial field is constant, the mode, or the bump'
    end select

    q = consistent_states(physics, degrees, lower, design, spatial_discretization_stencil)

  end function initial_field

  !-------------------------------------------------------------------!
  ! The rectangle's mode at every cell centre, cos(pi x / a) cos(pi y
  ! / b): the shape every check on the rectangle reads.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The operator applied to a field over the cells: the integrated
  ! flux balance of that field.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The operator alone against the laplacian of the mode, cell by
  ! cell: the balance over the area against kappa times minus pi^2
  ! (1/a^2 + 1/b^2) times the mode, which has no normal derivative at
  ! any wall. The error is reported over the cells that touch no
  ! wall, those that touch one, and those that touch two, so a wall's
  ! treatment is told apart from the interior's.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! At nu = 0 on a rectangle, the last instant against the exact mode
  ! and against the semi-discrete mode, whose frequency carries the
  ! operator's eigenvalue as built, by the rayleigh quotient of the
  ! mode. The first error holds space and time, the second time
  ! alone. The instant's components arrive node by node, degrees
  ! within a node. Nothing is said on any other geometry or design.
  !-------------------------------------------------------------------!

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

    ! minus the mode against its own balance, over the mode against
    ! itself by area
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

  !-------------------------------------------------------------------!
  ! One instant as one vtu file: every degree of every node, the
  ! components arriving node by node, degrees within a node.
  !-------------------------------------------------------------------!

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

!=====================================================================!
! packed from application/gti_chain.f90
!=====================================================================!
!=====================================================================!
! A horizon of blocks whose layouts need not agree.
!
! A multistep block holds one set of components per instant. A stage
! block holds, for each step, its stages and then the instant it
! arrives at. So a horizon that changes from one family to the other
! cannot keep its state in one array indexed the same way throughout,
! and the junction between two blocks stops being a contiguous copy.
!
! What it becomes is an index map, and the only thing it needs is a
! question each block already resolves for itself: where among its
! unknowns its k-th instant sits. A block hands its successor
! components, never rows, so that question is the whole of the
! interface between them.
!
! A block may reach back further than the block before it is long - a
! stage family spans one instant and a backward difference of order
! three on a degree-three equation looks back over nine - so what a
! block is given is gathered from whichever earlier block computed
! each instant, not from the one immediately before it.
!
!             WHOSE INSTANT IS IT
!
! An instant shared by two blocks is computed by the earlier and
! carried by the later, so the functional counts it once, under the
! block that computed it. The first block additionally owns the
! instants it was given, whose values are initial conditions: they
! contribute to the functional and not to any derivative of it,
! because they do not move.
!
!             THE EXPANSION ALONG A CHAIN
!
! Every order travels the junction the way the trajectory does. At
! order m each block solves against its own jacobian for a right side
! its own physics determines, with its carried rows set to what its
! predecessor found at those instants at that same order. Order zero
! is the march itself.
!
!             THE ORDER A COEFFICIENT SITS AT
!
! A series is indexed from zero, the order being the index, and every
! array that holds one is allocated that way on purpose. A section of
! such an array is indexed from one, so copying one into a fresh
! array and then striking out an order by its number strikes out the
! order below it. Where that happens the trajectory itself is wiped
! and every derivative comes out exactly zero, which is what it did.
!
!             WHAT IS REFUSED
!
! A chain of no blocks; a block that adds no more instants than its
! family reaches back over; and initial conditions that are not one
! value per degree over the instants the first block was given.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_chain

  use util_precision  , only : dp
  use operation_family , only : family
  use operation_grid   , only : grid
  use operation_expression, only : expression
  use gti_expansion    , only : family_holder, marches_by_stages, expansion, &
       & design_of_physics, design_of_steps
  use gti_block        , only : block_residual
  use gti_march        , only : imbalance, swept, solved_linear, fresh_stamp, partitioned, horizon_bounds, &
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
  use util_tally            , only : tally_order, tally_enter, tally_leave, at_horizon, at_block, at_stage

  implicit none

  private
  public :: chain_block, march_chain, chain_expansion, instant_components
  public :: first_of, chain_derivative, asymmetry
  public :: multiset_count, multiset_rank, multiset_of, num_designs_of
  public :: chain_stamps
  public :: expansion_substitutions
  public :: sink_costates


  !===================================================================!
  ! One block of a chain: its statement, where its instants sit among
  ! its unknowns, which global instants it spans, how many it was
  ! given, and what it computed.
  !===================================================================!

  type :: chain_block
     type(block_residual)  :: rows
     integer , allocatable :: instants_at(:)
     real(dp), allocatable :: state(:)
     ! WHERE THE BLOCK LIES on the horizon: its first and last instant
     ! as fine indices, and its stride, the fine indices between its
     ! own instants - one for a startup block over refined steps, the
     ! refinement for every block over the march's own steps, and one
     ! for both when there is no startup
     integer               :: first  = 0
     integer               :: last   = 0
     integer               :: stride = 1
     integer               :: given = 0
     integer               :: primary = 0
     ! the components one instant holds: the degrees at every node
     integer :: width = 0
     integer :: nodes = 1
     ! whether the block marches by stages: its functional is the
     ! stage quadrature, dt times the sum over stages of the tableau
     ! weight times the integrand at the stage; a multistep block's is
     ! the integrand at the instant
     logical :: staged = .false.
     ! the family and the steps the block was built from, so that its
     ! rows can be differentiated along a direction in the steps; each
     ! step as a fraction of the march's step it lies in, and which
     ! one, so that a direction in the march's steps reads on its own
     integer , allocatable :: coarse_step(:)
     class(family), allocatable :: scheme
     real(dp)     , allocatable :: dt(:)
     real(dp)              :: fraction = 1.0_dp
     ! whether its instants count in the functional: a startup's do
     ! not, its instants being what the first block is given
     logical               :: counted = .true.
     ! the imbalance its solve began at, against which its relative
     ! tolerance was measured
     real(dp) :: began = 0.0_dp
  end type chain_block

  !===================================================================!
  ! THE COSTATE OF A SINK. An unknown read by no row but its own is a
  ! sink of the block's reads graph: its column of the jacobian holds
  ! the diagonal alone, so the costate equation J^T lambda = g gives
  ! J_ii lambda_i = g_i on it exactly, and lambda_i = 0 wherever the
  ! functional does not read the unknown. Which unknowns are sinks is
  ! read off the compiled pattern, not declared, at any degree: in a
  ! stage block the arriving instant's highest degree is one, since
  ! the next step's stages read the lower degrees and the governing
  ! rows sit at the stages; in a multistep block none is, the
  ! governing row at the same instant reading every degree. The
  ! departure checks the transposition, the solve and the placing of
  ! the functional's gradient together, at the cost of reading the
  ! pattern once per block.
  !===================================================================!

  type :: sink_costates
     ! sinks at each degree over the chain, 0 to degrees - 1, of three
     ! kinds: carried unknowns, whose rows are identities; the last
     ! point of a block, read by the next block through the junction,
     ! which enters the right side and not the pattern; and the rest,
     ! the interior, where theory allows sinks in a stage block at the
     ! highest degree alone and in a multistep block none
     integer , allocatable :: carried(:), last(:), interior(:)
     ! max |J_ii lambda_i - g_i| over the sinks, every functional and multiset
     real(dp) :: departure = 0.0_dp
     ! max |g_i| and max |lambda_i| over the rows of the same solves
     real(dp) :: gradient  = 0.0_dp
     real(dp) :: costate   = 0.0_dp
     ! max |lambda_i| over the sinks the functional does not read
     real(dp) :: unread    = 0.0_dp
  end type sink_costates

contains

  !===================================================================!
  ! The block holding a fine instant, the latest that does - a block
  ! recomputes the instants it was given, so the latest is what the
  ! next block reads - and where the instant lies in it. None holds
  ! it, and the block is zero.
  !===================================================================!

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

  !===================================================================!
  ! The components a chain holds at one fine instant, and at one of
  ! the march's own instants, which is the fine one at the march's
  ! stride.
  !===================================================================!

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


  !===================================================================!
  ! What a block is given at the instants it shares with the ones
  ! before it, laid out the way its own carried rows expect: instant
  ! by instant at its own stride, the components within an instant.
  !===================================================================!

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

  !===================================================================!
  ! One block's statement and where its instants lie, read from its
  ! node of the expansion graph.
  !===================================================================!

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

  !===================================================================!
  ! The whole chain built and marched, block after block, each given
  ! what its predecessor computed at the instants they share.
  !===================================================================!

  subroutine march_chain(schemes, added, physics, degrees, steps, &
       & design, initial, chain, tower, dt, t, achieved, grid_design, left, nodes, spatial_discretization_stencil, &
       & startup)

    type(family_holder)   , intent(in) :: schemes(:)
    integer               , intent(in) :: added(:), degrees
    type(expression)      , intent(in) :: physics
    real(dp)              , intent(in) :: design, initial(:)
    class(grid)           , intent(in) :: steps
    type(chain_block), allocatable, intent(out) :: chain(:)
    ! THE GRAPH the chain is read from, the caller's, built here and
    ! outliving the march: every block lies at its node of it
    type(expansion), allocatable, intent(inout), target :: tower
    real(dp)         , allocatable, intent(out) :: dt(:), t(:)
    real(dp)              , intent(out) :: achieved
    real(dp), intent(in), optional     :: grid_design(:)
    type(imbalance), intent(out), optional :: left
    integer        , intent(in) , optional :: nodes
    type(stencil)  , intent(in) , optional :: spatial_discretization_stencil
    ! given, the first block's given instants are marched first by a
    ! stage family of order four, every step split this many ways,
    ! as block one of the chain: a startup that is part of the chain
    ! and so of every derivative, and reads the initial state at the
    ! first instant alone
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

    ! THE GRAPH the blocks are read from: one node per block, slice
    ! and component, with the couplings' relations. The expansion lays
    ! its blocks end to end, where the chain's blocks share the
    ! instants one hands the next, so the tower is built over each
    ! block's own span with each block's own steps laid end to end -
    ! a block reads only its own steps, and the sharing is the
    ! chain's junction. The startup, over its refined steps, is the
    ! first block of the same tower.
    knobs = [real(dp) ::]
    allocate(every(size(added) + before))
    if (with_startup) then
       fine  = [0.0_dp, (dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (given - 1) * r)]
       knobs = fine(2:)
       allocate(every(1) % scheme, source=crouzeix_three_stage())
    end if
    do b = 1, size(added)
       ! the step ending at a block's first instant: none at the
       ! horizon's first, and after a startup the two share an instant,
       ! so a positive placeholder no row reads takes the place of the zero
       ! a partition refuses
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
       ! The report kept is the first block's that did not converge:
       ! every block after it reads a state it never reached.
       if (present(left)) then
          if (before + b == 1) left = one_left
          if (left % converged .and. .not. one_left % converged) left = one_left
       end if
    end do
    call tally_leave()

  end subroutine march_chain

  !===================================================================!
  ! One block of a chain: given what its predecessor computed at the
  ! instants they share, or the initial conditions if it is first,
  ! then built and solved.
  !===================================================================!

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

    ! A block whose scheme keeps stages within a step is filed under
    ! the stage level, every other under the block level.
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

  !===================================================================!
  ! The instants one block owns, as its own local indices: the ones
  ! it computed, and the ones it was given as well when no counted
  ! block before it holds them. A startup block owns none.
  !===================================================================!

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

  !===================================================================!
  ! The functional and every derivative of the functional in the physics' design
  ! alone, to the order given: the recursion over one design, by the
  ! forward route, which at one design is what the gate chooses at
  ! every order. Order zero is the functional.
  !===================================================================!

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

  !===================================================================!
  ! What one block's tangent is frozen at: its own trajectory and the
  ! design, over its own unknowns.
  !===================================================================!

  subroutine frozen_at(b, design, unknowns, inputs)

    type(chain_block), intent(in) :: b
    real(dp)         , intent(in) :: design
    type(stored_directed_graph), intent(out) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)

    call frozen_inputs(b % state, design, b % rows % num_points(), unknowns, inputs)

  end subroutine frozen_at




  !===================================================================!
  ! Every block's tangent stamped at the trajectory already marched.
  !===================================================================!

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

  !===================================================================!
  ! The designs a chain's derivatives run over: the physics' parameter
  ! and, when the steps are designs, every weight of the grid.
  !===================================================================!

  integer function num_designs_of(tower)

    type(expansion), intent(in) :: tower

    real(dp), allocatable :: step_partials(:,:)
    real(dp) :: design

    call designs_of(tower, design, step_partials)
    num_designs_of = 1
    if (allocated(step_partials)) num_designs_of = 1 + size(step_partials, 2)

  end function num_designs_of

  !===================================================================!
  ! What the tower says the designs are: the physics' parameter, and
  ! the steps' partials in the weights when the weights are designs.
  ! Invalid input: a tower whose first design is not the parameter.
  !===================================================================!

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


  !===================================================================!
  ! The first entry of a table of derivatives, for a caller with one
  ! functional and one design.
  !===================================================================!

  pure real(dp) function first_of(table)

    real(dp), intent(in) :: table(:,:)

    first_of = table(1, 1)

  end function first_of

  !===================================================================!
  ! A direction in the march's steps read on a block's own: each of
  ! its steps takes the direction of the march's step it lies in,
  ! times its fraction of that step.
  !===================================================================!

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
  ! Which block holds one fine instant, the latest that does, and
  ! where in it.
  !===================================================================!

  pure subroutine holder_of(chain, fine, held_by, at)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: fine
    integer          , intent(out) :: held_by, at

    integer :: local

    call locate(chain, fine, held_by, local)
    at = 0
    if (held_by > 0) at = chain(held_by) % instants_at(local)

  end subroutine holder_of

  !===================================================================!
  ! THE DERIVATIVES OF EVERY ORDER BY ONE RECURSION.
  !
  ! With R(q, p) = 0 and F = f(q, p) over designs p, write x for what
  ! R reads - the state q, the parameter nu, the steps dt - and, for a
  ! multiset S of designs, x_S for the total derivative of x along S:
  ! the tangent w_S of the state, the grid's u_S of the steps, one for
  ! the parameter when S is the parameter alone. The total derivative
  ! of R along S is the sum over the set partitions of S of the
  ! partial of R along one x_B per block B, and the sum is zero. The
  ! partition with one block is A w_S with A = R_q, so
  !
  !    A w_S = -(the sum over the partitions with two or more blocks)
  !
  ! and that sum is one coefficient: R evaluated over derivative terms
  ! whose subsets are seeded with the x_T, T within S, the full subset
  ! seeded with zero - the product rule on subsets lists the
  ! partitions. The costate of F for S solves, by the Leibniz rule on
  ! A^T lambda = f_q,
  !
  !    A^T lambda_S = f_q along S - sum over T within S, T not S, of
  !                   (A along S less T)^T lambda_T
  !
  ! where A along a subset U is the state gradient of R along U: the
  ! coefficient of U with one more direction on each state component.
  ! The entry of the table for design j and multiset S is, by the
  ! reverse route,
  !
  !    T_jS = f_pj along S - sum over T within S of lambda_T^T (R_pj
  !           along S less T)
  !
  ! every term the coefficient of S with j as one more direction and
  ! no state seed on a subset holding j; by the forward route the
  ! entry for S is the coefficient of the full subset of f seeded with
  ! w_S. The sums run over the subsets of the positions of S, a
  ! repeated design being two positions, which counts the multinomial
  ! factors of a repeated derivative. The costate of order one is S
  ! empty, the hessian is order two, and no order is written out by
  ! hand.
  !
  !             THE COST
  !
  ! By the reverse route, one tangent per multiset of size below the
  ! order and one costate per functional and multiset of the same,
  ! then one contraction per design and multiset; by the forward
  ! route one tangent per multiset up to the order and one contraction
  ! each. The gate chooses by the top size, C(D + m - 1, m) against
  ! (1 + F) C(D + m - 2, m - 1). Along a chain the tangents are handed
  ! forward and the costates back exactly as at order one. The rows
  ! are linear in the state, and the block's carried rows hold given
  ! numbers, so neither varies with a design; the physics' partials
  ! are read from the rule at the points, the weights' from the
  ! family's action, the steps' from the grid, all exact.
  !
  !             WHAT IS REFUSED
  !
  ! An order below one, a route that is neither, and a tower whose
  ! first design is not the parameter.
  !===================================================================!

  subroutine chain_derivative(chain, tower, marks, functionals, degrees, order, route, &
       & table, node_measure, entries, designs, by_order, sinks)

    type(chain_block), intent(in) :: chain(:)
    type(expansion)  , intent(in) :: tower
    integer          , intent(in) :: marks(:)
    type(expression) , intent(in) :: functionals(:)
    integer                , intent(in) :: degrees, order, route
    ! one column per multiset of designs of the order's size, in the
    ! lexicographic order multiset_of names
    real(dp), allocatable  , intent(out) :: table(:,:)
    real(dp), intent(in), optional      :: node_measure(:)
    ! by the reverse route, given: every entry T_jS before one is
    ! chosen for the table, one per design j and multiset S of the
    ! size below; the departure among the entries of one multiset is
    ! the check on the route
    real(dp), allocatable, intent(out), optional :: entries(:,:,:)
    ! given, the count of designs run over from the first: one for the
    ! physics' parameter alone
    integer, intent(in), optional :: designs
    ! given, every order's table up to the order: the ones below by
    ! the forward route from the tangents in hand, the order asked for
    ! by the route given
    real(dp), allocatable, intent(out), optional :: by_order(:,:,:)
    ! given, the costates of the sinks checked over every solve; the
    ! forward route makes no costate and is refused
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

    ! THE TANGENTS of every multiset up to the top size, in increasing
    ! size, each block given what an earlier block found at the
    ! instants the block carries
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

    ! THE COSTATES of every functional and multiset up to the size
    ! below the order, the empty multiset's the costate of order one,
    ! each handed back along the chain
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

    ! THE ENTRIES, every design against every multiset of the size
    ! below; the table takes, for each multiset of the order's size,
    ! the entry whose design is the multiset's largest
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

    ! the table of one size by the forward route, from the tangents
    ! in hand: into the table at the order asked for, into by_order
    ! below
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

  !===================================================================!
  ! The steps' total derivatives along every multiset of designs up
  ! to one size, from the grid: zero along a multiset holding the
  ! parameter, and zero throughout when the steps are no design.
  !===================================================================!

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

  !===================================================================!
  ! THE SEEDS of one block for the positions of a multiset s, and one
  ! more position for an open design when one is given: for every
  ! nonempty subset of the positions, the state's total derivative
  ! along the designs at those positions - the tangent of that
  ! multiset, zero on a subset holding the open position, and zero on
  ! the full subset unless asked for - the steps' from the grid, and
  ! the parameter's, one on the subset of a single position holding
  ! the parameter. Column zero of the state's is the state.
  !===================================================================!

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

  !===================================================================!
  ! A nodal rule at one point over the seeded terms, with one more
  ! direction on each state component when asked: the derivative
  ! terms of the rule's value there, every coefficient a total
  ! derivative along the subset the mask names.
  !===================================================================!

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

  !===================================================================!
  ! The quadrature points of one owned instant of a block, and the
  ! weight each carries. A multistep block has one, the instant
  ! itself, weight one. A stage block has its s stages, laid before
  ! the arriving instant at offsets instants_at(k) - (s - i + 1) width,
  ! each weighted by the tableau weight of its stage; the arriving
  ! instant is not a quadrature point.
  !===================================================================!

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

    ! a stage block's first slice is the initial instant, which has no
    ! stages and no step ending at it, so it is not a quadrature point
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

  !===================================================================!
  ! The step ending at one of a block's instants over the seeded
  ! terms, times the measure of a node: the measure of an owned point
  ! with every total derivative of the measure.
  !===================================================================!

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

  !===================================================================!
  ! The partitions with two or more blocks of the total derivative
  ! of one block's rows along a multiset: the coefficient of the full
  ! subset with the state's full-subset seed zero. The time
  ! discretization rows are linear in the state, so theirs is the
  ! weights' derivative along each complement applied to the state's
  ! derivative along the rest; the physics' is read at the points;
  ! the carried rows hold given numbers and take no part.
  !===================================================================!

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

  !===================================================================!
  ! The costates at the subsets: for every subset of the positions of
  ! s, the costate of functional i for the designs at those positions,
  ! over one block's unknowns; the empty subset's is the costate of
  ! order one and the full subset's the costate of s. A
  ! contraction reads the costate of a subset's complement from here.
  !===================================================================!

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

  !===================================================================!
  ! One block's sinks from its compiled pattern: the columns holding
  ! their diagonal and nothing else. Zero weights are structural
  ! entries and count as reads, so a partial that happens to vanish
  ! at the frozen state makes no sink. A block that does not compile
  ! its tangent stops the program.
  !===================================================================!

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

  !===================================================================!
  ! The identity J_ii lambda_i = g_i on one block's sinks after one
  ! costate solve with right side g, the departure accumulated.
  !===================================================================!

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

  !===================================================================!
  ! The right side of one block for the costate of functional i and
  ! multiset s: the functional's gradient along s over the owned
  ! points, the measure carried as terms, less the rows' derivatives
  ! along every nonempty subset of the positions transposed against
  ! the costate of the complement.
  !===================================================================!

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

  !===================================================================!
  ! One block's part of the entry for design j and multiset s by the
  ! reverse route: the functional along s with j as one more
  ! direction over the owned points, less every costate of a
  ! complement against the rows along the rest with j.
  !===================================================================!

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

  !===================================================================!
  ! One block's part of the functional's total derivative along a
  ! multiset by the forward route: the coefficient of the full subset
  ! over the owned points, the tangent of the multiset seeded too.
  !===================================================================!

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

  !===================================================================!
  ! MULTISETS OF DESIGNS: the count of those of one size, the rank of
  ! one in lexicographic order, and the one at a rank. A multiset is
  ! a nondecreasing list of design indices. Invalid input: a list that
  ! is not nondecreasing within one to the count of designs.
  !===================================================================!

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

  !===================================================================!
  ! The largest departure among the entries of one multiset - T_jS
  ! over the distinct designs j of a multiset of the order's size,
  ! with S the rest - relative to the largest entry. The entries agree
  ! in theory and are not made to, which makes the departure the check
  ! on the reverse route.
  !===================================================================!

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

  !===================================================================!
  ! What the model says an expansion to the given order costs in
  ! substitutions, for the accounting layer to set beside what it
  ! counted: per order, one per block by the forward route at one
  ! design and one functional.
  !===================================================================!

  pure integer function expansion_substitutions(num_blocks, order) result(count)

    integer, intent(in) :: num_blocks, order

    count = num_blocks * route_substitutions(route_of(1, 1, order), 1, 1, order)

  end function expansion_substitutions




end module gti_chain

!=====================================================================!
! packed from application/gti_driver.f90
!=====================================================================!
!=====================================================================!
! What every driver in this directory needs and none should state
! for itself: the configuration named on the command line with the
! arguments after it applied over it, the time grid a configuration
! asks for, a clock, and the derivatives of a cosine.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_driver

  use iso_fortran_env  , only : int64
  use util_precision   , only : dp
  use operation_grid   , only : grid, uniform_grid, random_grid
  use gti_configuration, only : configuration, read_configuration, override
  use gti_march        , only : partitioned
  use gti_sweeps       , only : jacobian_of
  use view_directed_stored, only : stored_directed_graph
  use field_stored     , only : stored_field
  use operation_family , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : implicit_midpoint, crouzeix_two_stage, crouzeix_three_stage
  use operation_expression  , only : expression
  use physics_vanderpol     , only : van_der_pol_energy, van_der_pol_dissipation
  use gti_chain             , only : chain_block

  implicit none

  private
  public :: settings, chosen_grid, steps_of, clock, cosine, dense_jacobian
  public :: family_named, functional_named

contains

  !-------------------------------------------------------------------!
  ! The configuration named on the command line, --config=<name>, or
  ! the default given, then every other argument applied over it.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The time grid a configuration names, and the instants it makes.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The d-th derivative of the cosine at t.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The dense jacobian of the first chain block at its solved state,
  ! formed from the compiled tangent when the block provides one.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! One family, by name and order. A name or an order no family is
  ! built for is reported rather than refused, so a table may pass
  ! it over.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! One functional, by name, over an equation of the given degree. A
  ! name nothing is built for is reported rather than refused.
  !-------------------------------------------------------------------!

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

!=====================================================================!
! packed from application/graph_time_integrator.f90
!=====================================================================!
! The graph time integrator.
!
! It assembles what src holds and prints one table: a row for every
! scheme a configuration asks for, and a column for the functional
! and each of its derivatives in the design.
!
! Every row is one march and one expansion. The march solves the
! block whole, over a partition of the duration that is not uniform -
! the steps are drawn from a seed and scaled so that they sum to the
! duration exactly - and the expansion then reads the functional and
! its derivatives off the same jacobian, one solve per order.
!
!      ./graph_time_integrator --config=homogeneous
!      ./graph_time_integrator --config=homogeneous --max-derivative-degree=5
!
! A setting on the command line overrides the one in the file, and a
! setting that is not a setting stops the run rather than being
! passed over.
!
!             WHERE EVERY ROW STARTS
!
! A scheme cannot take its first step until it has instants behind it
! to look back at, and how many differs by family: a backward
! difference of order four wants eight, an Adams quadrature of order
! one wants a single one, and a stage family wants a single one
! whatever its order. Whatever fills those instants is not solved by
! the scheme; it is handed to it.
!
! If it were filled from a formula the rows would not be comparable.
! The widest scheme would hold a third of the horizon at numbers that
! are not a trajectory, would only integrate what remained, and would
! begin that from a state the equation would never have produced. Its
! functional would be mostly the formula and the narrowest scheme's
! mostly a solution, and the two would have no reason to agree.
!
! So the instants are integrated rather than invented. A stage family
! needs one instant and therefore no filler at all, so one is marched
! first over a refined grid across the startup, and every row takes
! its own reach from what that produced. That march is itself a chain
! of short blocks rather than one long one, because a block is solved
! whole and a long one costs far more than the several it could have
! been - which is the same junction the table's own rows use. Every row then begins from
! the same trajectory, holding one instant of it or eight is equally
! sound, and what separates the rows is how well each integrates,
! which is what the table is for.
!
! How many instants a row is handed and how many it works out is
! printed beside it, because they differ sharply: a backward
! difference of order four on an equation of degree four looks back
! over sixteen, so on a horizon of twenty-one it integrates five. Its
! functional is then mostly what it was given, and a reader who did
! not know that would take it for a peer of a row that integrated
! twenty.
!
! That is what automatic_order_conservation asks for. Turned off, the
! startup would have to be filled some other way, and there is no
! other way here that keeps the rows comparable, so the run says so
! and stops rather than printing a table that cannot be read across.
program graph_time_integrator

  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : uniform_grid, random_grid, designed_grid, fixed_grid
  use operation_expression  , only : expression
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
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
       & set_coarse_nodes, set_linear_budget
  use gti_sweeps            , only : route_of, forward_route, reverse_route
  use operation_minimization, only : relative, absolute, by_count, by_rate
  use gti_driver            , only : settings, chosen_grid, steps_of, family_named, clock, &
       & functional_named
  use gti_configuration     , only : configuration, read_configuration, override, show, &
       & lists, refuse_unknown, worded
  use util_tally            , only : tally_open, tally_close, tally_order, &
       & tally_enter, tally_leave, tally_amount, tally_event_of, &
       & tally_num_levels, tally_level_name, tally_event_name, &
       & at_expansion, at_horizon, wall_time

  implicit none


  type(configuration) :: cfg

  ! THE FIELD, when the configuration names a mesh: its room, the
  ! level below as a stencil over the nodes, the measure of each node,
  ! and the state at the first instant over every node. With no mesh
  ! there is one node, no spatial discretization stencil, and a measure of one: one
  ! node's equation, marched by the same chain.
  type(room)   , allocatable :: space
  type(stencil), allocatable :: spatial_discretization_stencil
  real(dp)     , allocatable :: volume(:), q0(:)
  real(dp) :: extent_a = 0.0_dp, extent_b = 0.0_dp
  integer  :: nodes = 1
  logical  :: over_field = .false.

  ! THE FUNCTIONALS the configuration names, and whether the grid's
  ! step weights are designs beside the physics' parameter
  type(expression)      , allocatable :: functionals(:)
  logical :: grid_designed = .false.
  ! the adaptive grid, discovered once and frozen: not a design
  logical :: grid_adaptive = .false.
  real(dp), allocatable :: adaptive_weights(:)

  call settings('homogeneous', cfg)
  call show(cfg)
  call set_linear_solver(cfg % linear_solver)
  call set_assembly(cfg % assembly)
  call set_storage(cfg % storage)
  call set_multigrid(cfg % multigrid)
  call set_sweep(cfg % sweep)
  call field_context(cfg)
  call chosen_functionals(cfg)
  call table(cfg)

contains

  !-------------------------------------------------------------------!
  ! The one instant a stage family needs, and it is consistent with
  ! the equation rather than merely plausible: the value and every
  ! derivative below the highest are chosen, and the highest is what
  ! the governing constraint then requires.
  !-------------------------------------------------------------------!


  !-------------------------------------------------------------------!
  ! How many instants the widest row that fits looks back over. A row
  ! that reaches past the horizon is not built, so it does not decide
  ! how long a startup the others need; zero means none of them fit.
  !-------------------------------------------------------------------!

  integer function widest_reach(cfg) result(widest)

    type(configuration), intent(in) :: cfg

    character(len=8) :: every(3)
    class(family), allocatable :: scheme
    logical :: staged, ok
    integer :: i, order, reach

    every  = ['bdf     ', 'adams   ', 'dirk    ']
    widest = 0

    ! how far back a family looks is the family's own answer
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

  !-------------------------------------------------------------------!
  ! The instants every row starts from, integrated rather than
  ! invented: a stage family over the startup, on a grid refined
  ! within each of its steps, sampled back at the coarse instants.

  !-------------------------------------------------------------------!
  ! One family, by name and order. A stage family says that it is
  ! one, since its block is laid out differently.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The heading, then one line per row.
  !-------------------------------------------------------------------!

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

    ! Each name sits over its own column, right against the digits.
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

    ! A column the row holds no expansion for is left empty rather
    ! than filled, there being no number to state under it.
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

  !-------------------------------------------------------------------!
  ! What the march left, by aspect, beneath the row that did not
  ! converge: how the norm splits by degree, where the largest entry
  ! sits, and which state the norm is steepest in.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The precision each block of a row needs: from its own family and
  ! smallest step, the norm of its tangent in closed form; from its
  ! own state, the size of what is subtracted; from the imbalance its
  ! solve began at, the target. The spacing those ask for names the
  ! least kind, block by block, since a chain may need more precision
  ! in one block than in another. Nothing is said when every block's
  ! least kind is this build's or below it, unless accounting is on.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! One row: a chain of blocks, marched and then expanded. A chain of
  ! one is a homogeneous row and takes the same path.
  !-------------------------------------------------------------------!

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

    ! the chain from the state at the first instant: a startup block
    ! over the first given instants where the family reaches back
    ! over more than one, then the row's own blocks
    call tally_enter(at_expansion)
    call tally_order(0)
    if (grid_designed) then
       ! the steps as the weights of a designed grid, which give the
       ! same steps back, so that the weights are designs of the tower
       weights = dt(2:cfg % instants)
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & designed_grid(cfg % time_duration), cfg % design, q0, chain, tower, dt, t, &
            & achieved, grid_design=weights, left=left, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else if (grid_adaptive) then
       ! the discovered steps as a given partition, frozen: no grid
       ! design, so the tower carries the physics' parameter alone
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & fixed_grid(adaptive_weights), cfg % design, q0, chain, tower, dt, t, achieved, &
            & left=left, nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, &
            & startup=cfg % startup_refinement)
    else
       call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
            & chosen_grid(cfg), cfg % design, q0, chain, tower, dt, t, achieved, left=left, &
            & nodes=nodes, spatial_discretization_stencil=spatial_discretization_stencil, startup=cfg % startup_refinement)
    end if
    ! Every derivative is taken at the state the march reached, so a
    ! row that did not converge has none to take and only its value is
    ! expanded.
    if (.not. left % converged) then
       reported = 0
    else
       reported = cfg % max_derivative_degree
    end if
    call chain_expansion(chain, tower, functionals, nd, reported, f, node_measure=volume)
    call tally_leave()

    call show_row(labelled(names, orders), cfg % instants - given, f(:, 1), left, &
         & cfg % max_derivative_degree)
    ! every functional after the first, under the row, expanded from
    ! the same state series
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

    ! the first derivatives by the routes: where the grid is a design
    ! they are the only account of it, and where the routes are checked
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

  !-------------------------------------------------------------------!
  ! The derivatives of every functional in every design by the routes.
  ! At first order the gate chooses from the counts: forward, one
  ! solve per design, or reverse, one per functional; at second order
  ! the reverse route, when the gate chooses it. With the grid's weights
  ! among the designs the steps are homogeneous of degree zero in
  ! them, so the weights against the gradient sum to zero - a check
  ! of the whole chain rule through the grid - and the physics' column
  ! is the expansion's first order. Asked for, the other route is run
  ! too and the two are compared over the whole table. The first
  ! instants a family reaches back over are held as given, so their
  ! own dependence on the steps is not carried.
  !-------------------------------------------------------------------!

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

    ! the designs are the tower's: the parameter, and the weights of
    ! the steps when the grid was designed
    num_functionals = size(functionals)
    call chain_stamps(chain, tower, functionals, nd, marks, node_measure=volume)
    num_designs = num_designs_of(tower)
    if (grid_designed) p = dt(2:cfg % instants)

    route = route_of(num_designs, num_functionals, 1)
    call chain_derivative(chain, tower, marks, functionals, nd, 1, route, df, node_measure=volume)

    write(*,'(a,a,a,i0,a,i0,a,es10.2)') '      first derivatives by the ', &
         & trim(merge('forward', 'reverse', route == forward_route)), ' route, designs ', &
         & num_designs, ' functionals ', num_functionals, &
         & ':  physics column against the expansion ', &
         & maxval(abs(df(:, 1) - f(1, :)) / max(1.0_dp, abs(f(1, :))))
    if (grid_designed) then
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

    ! the costates of the sinks: J_ii lambda_i = g_i on every unknown
    ! no row reads, and lambda_i = 0 where the functional does not
    ! read it either - in theory the highest degree at the arriving
    ! instants of a stage block, and no unknown of a multistep block
    if (lists(cfg % check, 'sinks')) then
       call chain_derivative(chain, tower, marks, functionals, nd, 1, reverse_route, other, &
            & node_measure=volume, sinks=sinks)
       call shown_sinks(sinks, nd)
    end if

    ! the derivatives of every order above one, when the grid is
    ! designed, by the route the gate chooses: one table per order,
    ! one column per multiset of designs; by the reverse route the
    ! entries of one multiset agree in theory and are not made to; the
    ! entry of the parameter alone is the expansion's coefficient
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

  !-------------------------------------------------------------------!
  ! The sinks by degree and the two departures, each relative to the
  ! largest entry of the solves it was read from.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The functionals the configuration names, in its order, and the
  ! designs: the physics' parameter always, the grid's weights when
  ! named.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! At kappa = 0 with a constant field every node is one node's
  ! equation: the field's functional over the area is the node's,
  ! order by order. The node's march is the same chain, startup
  ! included, from the first node's own first instant.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! Every instant as one vtu file, numbered, so paraview reads the
  ! series as time.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The field the configuration names, or one node when it names no
  ! mesh: the room, the spatial discretization stencil, the coarse cells a multigrid
  ! coarsens the nodes by, the measure of each node, and the state at
  ! the first instant. The operator alone is checked here when asked,
  ! before any march.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! Two numbers from a setting, one per coordinate.
  !-------------------------------------------------------------------!

  subroutine pair_of(text, x, y, subject)

    character(len=*), intent(in)  :: text, subject
    real(dp)        , intent(out) :: x, y

    character(len=32), allocatable :: w(:)

    w = worded(text)
    if (size(w) /= 2) error stop 'graph_time_integrator: two ' // subject // ', one per coordinate'
    read(w(1), *) x
    read(w(2), *) y

  end subroutine pair_of

  !-------------------------------------------------------------------!
  ! The families a row names, and the instants split among them. A
  ! row whose blocks would add no more instants than their families
  ! reach back over is not built.
  !-------------------------------------------------------------------!

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


  !-------------------------------------------------------------------!
  ! Every row the configuration asks for.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The adaptive grid: when the grid is adaptive, an order-four
  ! diagonally implicit march to the tolerance discovers the steps,
  ! and they are frozen as the run''s grid - its instant count set from
  ! them. The grid is over time alone; a spatial field is refused. The
  ! steps are not designs: the frozen grid is an ordinary grid to the
  ! expansion, so no grid sensitivity is taken.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The steps and their times: the discovered partition when the grid
  ! is adaptive, the chosen grid's otherwise.
  !-------------------------------------------------------------------!

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

    ! Before anything is measured against the horizon, since a word
    ! this program has nothing for reaches back over nothing and
    ! would be reported as a horizon too narrow to hold it.
    !
    ! physics is refused rather than dispatched on because one
    ! integrand is built. Were it neither, the run would state a
    ! physics in its heading and integrate a different one.
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

    if (cfg % accounting) call tally_open(cfg % max_derivative_degree)

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


  !-------------------------------------------------------------------!
  ! The names a configuration lists, in the order it lists them.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! Every ordered tuple of this many distinct families, at every
  ! order - one order for the whole tuple, or, when mixed orders are
  ! asked for, every tuple of orders. One arity serves the
  ! homogeneous rows, the pairs and the triples alike.
  !-------------------------------------------------------------------!

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
       ! the tuple of families, the last position varying fastest
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

  !-------------------------------------------------------------------!
  ! What the run spent, one table per measurement asked for: the
  ! amount at each level of the hierarchy against the derivative order
  ! it was spent on, and then the same amounts as ratios of one order
  ! to another.
  !
  ! The ratio is what a higher order costs against a lower one, so the
  ! entry at row i and column j is the amount at order i over the
  ! amount at order j. A column whose order spent nothing leaves its
  ! ratio empty rather than dividing by it.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! What the route's cost model said each order would cost in
  ! substitutions, beside what was counted. One block per row is what
  ! the homogeneous table builds, so the model is read at one block.
  !-------------------------------------------------------------------!

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

    ! Levels nest, so a level's time already holds the time of the
    ! levels opened inside it and a sum over levels would count the
    ! same seconds again. The expansion is opened once for a whole
    ! row and closed after every order has been taken, so its time
    ! belongs to no single order and is filed where the row began.
    ! The horizon is opened once per order, which is what a time
    ! against an order means, so it is the one the ratios are taken
    ! from. A count is filed at one level only and does sum.
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

  !-------------------------------------------------------------------!
  ! Row over column, the whole run. An order that spent nothing is no
  ! denominator, and its column is left empty.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! Which side of the iteration cap the run sits on. Past it every
  ! order spends the whole budget instead of converging, and a ratio
  ! measured there reports the cap and not the order.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! An amount summed over every level of the hierarchy, for one order
  ! and one event.
  !-------------------------------------------------------------------!

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
