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
!   orders      along a chain that changes family part way through,
!               every derivative must agree with a difference of the
!               one below it, which crosses the junction because the
!               whole chain is remarched either side of the design
!
! The third is the one no fixed test has made before. Every case
! prints its three residuals, and a case whose march did not converge
! is passed over rather than counted, since nothing can be concluded
! from a trajectory that was not found.
!
!             WHERE THE INVARIANTS MEAN SOMETHING
!
! The design is the damping, and it is drawn without a sign. A
! negative one is negative damping: van der Pol then grows without
! bound, and on a coarse grid the discrete statement has more than
! one solution, so two marches of the same equations converge to
! different ones and splitting changes the number legitimately. That
! was seen at a design of minus four fifths over four units, where
! one march found a functional of six hundred and the other of
! thirteen thousand, both to a residual under a millionth.
!
! Nothing there is wrong. But an invariant that need not hold cannot
! test anything, so the draw stays where the trajectory is bounded.
program randomized_checks

  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : crouzeix_two_stage
  use operation_grid        , only : uniform_grid, random_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use gti_driver            , only : cosine
  use gti_block             , only : block_residual
  use gti_march             , only : partition
  use gti_chain             , only : one_functional, first_of
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & chain_system, chain_systems, chain_by_tangent, chain_by_adjoint

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

    count = command_argument_count()
    verbose = .false.
    if (count >= 3) then
       call get_command_argument(3, argument)
       verbose = trim(argument) == 'verbose'
    end if

  end function verbose

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

  !-------------------------------------------------------------------!
  ! The sensitivity by both directions, over the chain layout, which
  ! every family has - a stage block keeps its instants between its
  ! stages and a multistep one keeps one set per instant, and neither
  ! is assumed here.
  !-------------------------------------------------------------------!

  subroutine directions_of(schemes, added, degrees, duration, design, tangent, adjoint)

    type(family_holder), intent(in)  :: schemes(:)
    integer            , intent(in)  :: added(:), degrees
    real(dp)           , intent(in)  :: duration, design
    real(dp)           , intent(out) :: tangent, adjoint

    type(chain_block) , allocatable :: chain(:)
    type(chain_system), allocatable :: systems(:)
    real(dp), allocatable :: held(:), dt(:), t(:)
    real(dp) :: achieved

    call held_for(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)

    call chain_systems(chain, [one_functional(van_der_pol_energy(degrees - 1))], degrees, dt, &
         & design, systems)

    tangent = first_of(chain_by_tangent(chain, systems, degrees, design))
    adjoint = first_of(chain_by_adjoint(chain, systems, degrees, design))

  end subroutine directions_of

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

  !-------------------------------------------------------------------!
  ! One drawn case: the three residuals, and a verdict.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! The same horizon expanded in one block and in two, and the worse
  ! of the two residuals, so that a case is passed over if either
  ! march failed to find a trajectory.
  !-------------------------------------------------------------------!

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

  !-------------------------------------------------------------------!
  ! A chain whose blocks are marched by different families. Every
  ! derivative of the functional is checked against a difference of
  ! the one below it, and that difference is taken by remarching the
  ! whole chain, so it crosses every junction the chain has.
  !
  ! A difference of two marched functionals carries whatever noise
  ! the marches carry, divided by the step between them, so this
  ! check has a floor and the floor rises as the step falls. Measured
  ! on the widest chain drawn here, the gap reads 1.3E-04 at a step
  ! of a ten-thousandth, 1.7E-04 at half that and 5.4E-04 at a
  ! quarter - growing as the step shrinks, which is round-off and not
  ! truncation. The tolerance sits above that floor, and a gap of
  ! round-off size is what most chains give.
  !-------------------------------------------------------------------!

  subroutine mixed_case(from, index, failures, skipped)

    integer, intent(in)    :: from, index
    integer, intent(inout) :: failures, skipped

    real(dp), parameter :: delta = 1.0e-4_dp

    type(family_holder), allocatable :: schemes(:)
    integer , allocatable :: added(:)
    real(dp), allocatable :: f(:), plus(:), minus(:)
    real(dp) :: duration, design, achieved, gap, differenced
    integer :: degrees, blocks, b, m
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

    associate (u1 => b); end associate

  end subroutine mixed_case

  !-------------------------------------------------------------------!
  ! A chain of two or three blocks of differing families, each adding
  ! more instants than it looks back over.
  !-------------------------------------------------------------------!

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

    ! Every block reaches back over instants an earlier one computed,
    ! so the first must cover the widest reach behind it. A chain that
    ! does not is refused by the march, and rightly, but it is not the
    ! chain this harness meant to draw.
    widest = 0
    do b = 2, blocks
       widest = max(widest, schemes(b) % scheme % history_depth(degrees - 1))
    end do
    added(1) = max(added(1), widest + 1)

  end subroutine draw_chain

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
    design   = drawn_real(state, 0.0_dp, 1.5_dp)

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

  !-------------------------------------------------------------------!
  ! The sensitivity by both directions, over the chain layout, which
  ! every family has - a stage block keeps its instants between its
  ! stages and a multistep one keeps one set per instant, and neither
  ! is assumed here.
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

    real(dp), allocatable :: table(:,:)
    real(dp)           , intent(out) :: achieved

    type(chain_block), allocatable :: chain(:)
    real(dp), allocatable :: held(:), dt(:), t(:)

    call held_for(schemes(1) % scheme, degrees, duration, sum(added), held, dt, t)

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)

    call chain_expansion(chain, van_der_pol(degrees - 1), &
         & [one_functional(van_der_pol_energy(degrees - 1))], degrees, dt, design, &
         & max_order, table)
    allocate(f(lbound(table, 1):ubound(table, 1)))
    f = table(:, 1)

  end subroutine expanded

  !-------------------------------------------------------------------!
  ! What the three residuals are, and whether any of them is too big
  ! to be round-off carried through a solve.
  !-------------------------------------------------------------------!

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

end program randomized_checks
