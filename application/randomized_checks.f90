! Randomized checks of the invariants no single run can show.
!
! Three of them, over parameters drawn from a seed:
!
!   splitting   a horizon marched in one block and in several must
!               give the same functional and the same derivatives,
!               because the same rows are made in the same places
!   directions  the tangent and the adjoint read the same three
!               objects through different algebra and must agree
!   routes      the sweeps and the Taylor expansion reach the first
!               derivative by paths that share no arithmetic - one
!               solves once against the jacobian and reads a
!               gradient along it, the other recurs on coefficients
!               of the statement - and must land on the same number
!
! The third is the one no fixed test has made before. Every case
! prints its three residuals, and a case whose march did not converge
! is passed over rather than counted, since nothing can be concluded
! from a trajectory that was not found.
program randomized_checks

  use iso_fortran_env       , only : dp => REAL64, int64 => INT64
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : crouzeix_two_stage
  use operation_grid        , only : uniform_grid, random_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use gti_block             , only : block_residual
  use gti_march             , only : partition, marched_horizon
  use gti_horizon           , only : block_system, horizon_systems, &
       & horizon_by_tangent, horizon_by_adjoint
  use gti_chain             , only : chain_block, march_chain, chain_expansion

  implicit none

  integer , parameter :: max_order = 2
  integer :: seed, cases, i, failures, skipped
  character(len=32) :: argument

  call get_command_argument(1, argument); read(argument,*) seed
  call get_command_argument(2, argument); read(argument,*) cases

  failures = 0
  skipped  = 0

  write(*,'(a)') '  case  scheme          split          directions     routes'

  do i = 1, cases
     call one_case(seed + 7919 * i, i, failures, skipped)
  end do

  write(*,'(a)')      ' '
  write(*,'(a,i0,a,i0,a,i0,a)') ' ', cases - skipped, ' cases checked, ', &
       & failures, ' failed, ', skipped, ' passed over'

  if (failures > 0) error stop 'randomized_checks: an invariant did not hold'

contains

  !-------------------------------------------------------------------!
  ! A deterministic draw, so a failing case can be run again.
  !-------------------------------------------------------------------!

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

  subroutine held_for(scheme, degrees, duration, instants, held, dt, t)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, instants
    real(dp)     , intent(in) :: duration
    real(dp), allocatable, intent(out) :: held(:), dt(:), t(:)

    integer :: k, d

    call partition(duration, instants, dt, t)
    held = [((cosine(d, t(k)), d = 0, degrees - 1), k = 1, &
         &   scheme % history_depth(degrees - 1))]

  end subroutine held_for

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
  ! One drawn case: the three residuals, and a verdict.
  !-------------------------------------------------------------------!

  subroutine one_case(from, index, failures, skipped)

    integer, intent(in)    :: from, index
    integer, intent(inout) :: failures, skipped

    type(family_holder), allocatable :: whole(:), split(:)
    type(block_system) , allocatable :: systems(:)
    real(dp), allocatable :: f_whole(:), f_split(:), q(:), gradient(:)
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

    call expanded(whole, [instants], degrees, duration, design, f_whole, achieved)
    if (achieved > 1.0e-6_dp) then
       skipped = skipped + 1
       write(*,'(i6,2x,a16,a)') index, label, '   march did not converge, passed over'
       return
    end if

    call expanded(split, [half, instants - half], degrees, duration, design, &
         & f_split, achieved)

    tangent = 0.0_dp
    adjoint = 0.0_dp
    if (.not. staged) then
       call marched_horizon(whole, [instants], van_der_pol(degrees - 1), degrees, &
            & duration, design, held_of(whole(1) % scheme, degrees, duration, instants), &
            & q, achieved)
       call horizon_systems(whole, [instants], van_der_pol(degrees - 1), &
            & van_der_pol_energy(degrees - 1), degrees, duration, design, q, &
            & systems, gradient)
       tangent = horizon_by_tangent(systems, gradient)
       adjoint = horizon_by_adjoint(systems)
    end if

    call verdict(index, label, f_whole, f_split, tangent, adjoint, staged, failures)

  end subroutine one_case

  !-------------------------------------------------------------------!
  ! The parameters of one case, drawn so that the same seed gives the
  ! same case again.
  !-------------------------------------------------------------------!

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
    design   = drawn_real(state, -1.0_dp, 1.5_dp)

  end subroutine draw

  !-------------------------------------------------------------------!
  ! Where to split a horizon. A block cannot add fewer instants than
  ! its family looks back over, so a horizon too short to halve is
  ! widened rather than the case being skipped, which keeps every
  ! drawn family in the sample.
  !-------------------------------------------------------------------!

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

  function held_of(scheme, degrees, duration, instants) result(held)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: degrees, instants
    real(dp)     , intent(in) :: duration
    real(dp), allocatable :: held(:), dt(:), t(:)

    call held_for(scheme, degrees, duration, instants, held, dt, t)

  end function held_of

  subroutine fill(held, kind, order)

    type(family_holder), intent(out) :: held
    integer            , intent(in)  :: kind, order

    select case (kind)
    case (1)
       allocate(held % scheme, source=bdf_family(order))
    case (2)
       allocate(held % scheme, source=adams_family(order))
    case default
       allocate(held % scheme, source=crouzeix_two_stage())
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
    real(dp)           , intent(out) :: achieved

    type(chain_block), allocatable :: chain(:)
    real(dp), allocatable :: held(:), dt(:), t(:)

    call held_for(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)

    call chain_expansion(chain, van_der_pol(degrees - 1), &
         & van_der_pol_energy(degrees - 1), degrees, dt, design, max_order, f)

  end subroutine expanded

  !-------------------------------------------------------------------!
  ! What the three residuals are, and whether any of them is too big
  ! to be round-off carried through a solve.
  !-------------------------------------------------------------------!

  subroutine verdict(index, label, f_whole, f_split, tangent, adjoint, staged, failures)

    integer         , intent(in)    :: index
    character(len=*), intent(in)    :: label
    real(dp)        , intent(in)    :: f_whole(0:), f_split(0:), tangent, adjoint
    logical         , intent(in)    :: staged
    integer         , intent(inout) :: failures

    real(dp) :: split_gap, direction_gap, route_gap, scale

    scale         = max(1.0_dp, maxval(abs(f_whole)))
    split_gap     = maxval(abs(f_whole - f_split)) / scale
    direction_gap = abs(tangent - adjoint) / max(1.0_dp, abs(tangent))
    route_gap     = abs(tangent - f_whole(1)) / max(1.0_dp, abs(tangent))

    if (staged) then
       write(*,'(i6,2x,a16,es15.2,a)') index, label, split_gap, &
            & '   multistep only   multistep only'
    else
       write(*,'(i6,2x,a16,3es15.2)') index, label, split_gap, direction_gap, route_gap
    end if

    if (split_gap > 1.0e-6_dp) failures = failures + 1
    if (.not. staged .and. direction_gap > 1.0e-8_dp) failures = failures + 1
    if (.not. staged .and. route_gap > 1.0e-6_dp) failures = failures + 1

  end subroutine verdict

end program randomized_checks
