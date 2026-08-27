! What one formation and one solve cost, apart from the iteration
! that calls them.
!
! chain_stamps stamps the jacobian; the tangent route forms the rest.
! The factorisation is kept with the system, so a solve here is one
! substitution against it. Timing them apart from the march says which
! of the two carries the growth, and whether the march's own growth is
! more than the two of them together.
program solve_cost

  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_grid        , only : uniform_grid
  use operation_expression  , only : expression
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder, expansion
  use gti_driver            , only : clock, cosine
  use gti_march             , only : partition
  use gti_chain             , only : first_of
  use gti_sweeps            , only : forward_route
  use gti_chain             , only : chain_block, march_chain, &
       & chain_stamps, chain_derivative

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
    type(expansion), allocatable, target :: tower
    integer, allocatable :: marks(:)
    type(expression)       :: energy(1)
    type(bdf_family) :: scheme
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:), table(:,:)
    real(dp) :: achieved, duration, design, marched, formed, solved_in, tangent
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
         & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
    marched = clock() - marched

    energy(1) = van_der_pol_energy(degrees - 1)
    formed = clock()
    call chain_stamps(chain, tower, energy, degrees, marks)
    formed = clock() - formed

    n = chain(1) % rows % num_unknowns()

    solved_in = clock()
    call chain_derivative(chain, tower, marks, energy, degrees, 1, forward_route, table)
    tangent   = first_of(table)
    solved_in = clock() - solved_in

    write(*,'(i10,3f11.3,2es15.3)') n, marched, formed, solved_in, &
         & achieved, 1.0e-12_dp

  end subroutine cost_at

end program solve_cost
