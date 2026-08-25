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

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: configuration, read_configuration, override, show

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
     integer  :: krylov_above             = huge(1)

     real(dp) :: time_duration = 7.0_dp
     real(dp) :: design        = 1.0_dp

     logical  :: automatic_order_conservation = .true.
     logical  :: mixed_orders                 = .false.

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
    case ('krylov_above')
       read(value, *) cfg % krylov_above
    case ('time_duration')
       read(value, *) cfg % time_duration
    case ('design')
       read(value, *) cfg % design
    case ('automatic_order_conservation')
       read(value, *) cfg % automatic_order_conservation
    case ('mixed_orders')
       read(value, *) cfg % mixed_orders
    case default
       write(*,'(a)') ' this is not a setting: ' // levelled(name)
       error stop 'gti_configuration: every setting given is one that exists'
    end select

    ! A count below the least value it means anything at is refused
    ! here, where it is given, rather than where it is first indexed
    ! by or counted over.
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
    write(*,'(a,i0)')      '   krylov above             ', cfg % krylov_above
    write(*,'(a,f0.4)')    '   design                   ', cfg % design
    write(*,'(a,i0)')      '   max derivative degree    ', cfg % max_derivative_degree
    write(*,'(a,i0)')      '   max discretization order ', cfg % max_discretization_order
    write(*,'(a,a)')       '   families                 ', trim(cfg % families)
    write(*,'(a,a)')       '   combinations             ', trim(cfg % combinations)
    write(*,'(a,l1)')      '   automatic order conservation ', cfg % automatic_order_conservation
    write(*,'(a,l1)')      '   mixed orders             ', cfg % mixed_orders

  end subroutine show

end module gti_configuration
