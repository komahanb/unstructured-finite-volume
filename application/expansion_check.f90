! The functional and its derivatives in the design, to fourth order,
! each checked against a difference of the one below it.
!
! The expansion is exact given the orders beneath it, so the m-th
! derivative and a central difference of the (m-1)-th must agree to
! the accuracy of the difference and no better. Running the whole
! expansion at three designs gives every check at once, and each one
! is independent: the m-th derivative is computed from the statement,
! never from the (m-1)-th.
program expansion_check

  use util_precision  , only : dp
  use operation_family     , only : family
  use operation_family_bdf , only : bdf_family
  use operation_family_adams, only : adams_family
  use physics_vanderpol    , only : van_der_pol, van_der_pol_energy
  use gti_march            , only : partition, set_stopping
  use operation_minimization, only : relative, by_rate
  use gti_expansion        , only : expansion, family_holder
  use operation_grid       , only : uniform_grid
  use operation_family_dirk, only : crouzeix_two_stage
  use gti_chain            , only : chain_block, march_chain, chain_expansion, one_functional

  implicit none

  character(len=32) :: argument

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  integer , parameter :: instants = 11
  integer , parameter :: max_order = 4
  real(dp), parameter :: duration = 2.0_dp
  real(dp), parameter :: design = 1.0_dp
  ! The difference step is not chosen. A central difference of a
  ! quantity known to a relative error tau has round-off tau f / delta
  ! and truncation delta^2 f''' / 6, balanced at delta = tau^(1/3), where
  ! the error is tau^(2/3). tau is the march's own relative tolerance,
  ! which is what the expanded functional is known to.
  real(dp) :: delta, tau

  ! An optional relative tolerance; the default otherwise.
  tau = 1.0e-12_dp
  call get_command_argument(1, argument)
  if (len_trim(argument) > 0) read(argument, *) tau
  call set_stopping(tau, relative, by_rate, 100)
  delta = tau ** (1.0_dp / 3.0_dp)
  write(*,'(a,es9.2,a,es9.2,a,es9.2)') ' relative tolerance', tau, &
       & '   difference step', delta, '   expected agreement tau^(2/3)', tau ** (2.0_dp / 3.0_dp)

  call expansion_of('bdf 2', bdf_family(2), .false.)
  call expansion_of('adams-moulton 3', adams_family(3), .false.)
  call expansion_of('crouzeix two-stage', crouzeix_two_stage(), .true.)

contains

  pure real(dp) function exact(d, t) result(q)

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

  end function exact

  subroutine expanded(scheme, staged, design_value, f)

    class(family), intent(in) :: scheme
    logical      , intent(in) :: staged
    real(dp)     , intent(in) :: design_value
    real(dp), allocatable, intent(out) :: f(:)

    type(family_holder) :: holder(1)
    type(chain_block), allocatable :: chain(:)
    type(expansion)  , allocatable, target :: tower
    real(dp), allocatable :: dt(:), t(:), held(:), table(:,:)
    real(dp) :: achieved
    integer :: k, d

    call partition(duration, instants, dt, t)
    held = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, scheme % history_depth(degrees - 1))]
    allocate(holder(1) % scheme, source=scheme)
    associate (u1 => staged); end associate

    ! one block, marched as a chain of one, then expanded by the
    ! recursion over the one design
    call march_chain(holder, [instants], van_der_pol(state_degree), degrees, uniform_grid(duration), &
         & design_value, held, chain, tower, dt, t, achieved)
    call chain_expansion(chain, tower, [one_functional(van_der_pol_energy(state_degree))], degrees, &
         & max_order, table)
    allocate(f(0:max_order))
    f(0:) = table(:, 1)

  end subroutine expanded

  subroutine expansion_of(title, scheme, staged)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme
    logical         , intent(in) :: staged

    real(dp), allocatable :: f(:), plus(:), minus(:)
    real(dp) :: differenced(max_order)
    integer :: m

    call expanded(scheme, staged, design, f)
    call expanded(scheme, staged, design + delta, plus)
    call expanded(scheme, staged, design - delta, minus)

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

end program expansion_check
