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

  implicit none

  private
  public :: settings, chosen_grid, steps_of, clock, cosine

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

end module gti_driver
