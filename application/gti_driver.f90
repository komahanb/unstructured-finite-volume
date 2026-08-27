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
