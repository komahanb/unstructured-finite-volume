! What one formation and one solve cost, apart from the iteration
! that calls them.
!
! chain_systems forms the jacobian, the design rate and the gradient.
! The factorisation is kept with the system, so a solve here is one
! substitution against it. Timing them apart from the march says which
! of the two carries the growth, and whether the march's own growth is
! more than the two of them together.
program solve_cost

  use iso_fortran_env       , only : dp => REAL64, int64 => INT64
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_grid        , only : uniform_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use util_factorisation    , only : dense_factorisation
  use gti_march             , only : partition
  use gti_chain             , only : chain_block, march_chain, &
       & chain_system, chain_systems

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
    type(chain_system), allocatable :: systems(:)
    type(bdf_family) :: scheme
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:), rhs(:), x(:)
    real(dp) :: achieved, duration, design, marched, formed, solved_in
    integer  :: n, d, j

    duration = 3.0_dp
    design   = 1.0_dp
    scheme   = bdf_family(2)

    allocate(schemes(1))
    allocate(schemes(1) % scheme, source=scheme)
    added = [instants]

    call partition(duration, instants, dt, t)
    held = [((cosine(d, t(j)), d = 0, degrees - 1), j = 1, &
         &   scheme % history_depth(degrees - 1))]

    marched = clock()
    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)
    marched = clock() - marched

    formed = clock()
    call chain_systems(chain, van_der_pol_energy(degrees - 1), degrees, dt, &
         & design, systems)
    formed = clock() - formed

    n = size(systems(1) % a, 1)
    allocate(rhs(n), source=1.0_dp)

    solved_in = clock()
    call systems(1) % factor % substitute(rhs, x, transposed=.false.)
    solved_in = clock() - solved_in

    write(*,'(i10,3f11.3,2es15.3)') n, marched, formed, solved_in, &
         & achieved, 1.0e-12_dp

  end subroutine cost_at

  real(dp) function clock() result(s)

    integer(int64) :: count, rate

    call system_clock(count, rate)
    s = real(count, dp) / real(rate, dp)

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

end program solve_cost
