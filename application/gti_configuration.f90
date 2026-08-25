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
     ! comparison the run makes against something known.
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
