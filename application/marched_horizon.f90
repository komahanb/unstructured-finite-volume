! A horizon marched block after block, across a change of scheme.
!
! A block reaches back over instants that begin before it does, so
! every block after the first overlaps its predecessor by exactly
! what its family reaches, and carries there what the predecessor
! computed. That hand-over is the junction, and it is what lets a
! horizon change scheme part way through.
!
! Two checks. Splitting a horizon into two blocks of the same family
! must change nothing at all: the same rows are made, in the same
! places, from the same steps, so the two states must agree to the
! tolerance newton was asked for. Then a horizon whose first block is
! a backward difference and whose second is an Adams quadrature is
! marched against the harmonic oscillator, where the solution is the
! cosine, and the error is watched as the steps are halved.
program marched_horizon

  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use physics_vanderpol     , only : van_der_pol
  use gti_expansion         , only : family_holder, expansion
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use physics_vanderpol     , only : van_der_pol_energy
  use gti_march             , only : horizon_bounds, partition, unknowns_graph
  use gti_sweeps            , only : functional_of
  use gti_chain             , only : one_functional, first_of
  use gti_chain             , only : chain_block, march_chain, instant_components, &
       & chain_system, chain_systems, chain_by_tangent, chain_by_adjoint
  use operation_grid        , only : uniform_grid

  implicit none

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  real(dp), parameter :: duration = 2.0_dp

  call splitting_changes_nothing()
  call across_a_change_of_scheme()
  call sensitivity_across_the_junction()

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

  !-------------------------------------------------------------------!
  ! The initial conditions a first block carries: the cosine and its
  ! derivatives over the instants that block reaches back across.
  !-------------------------------------------------------------------!

  function initial_for(scheme, instants) result(held)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: instants
    real(dp), allocatable :: held(:)

    real(dp), allocatable :: dt(:), t(:)
    integer :: k, d

    call partition(duration, instants, dt, t)
    held = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, scheme % history_depth(degrees - 1))]

  end function initial_for

  !-------------------------------------------------------------------!
  ! One horizon marched. The design of zero makes the physics the
  ! harmonic oscillator, whose solution is the cosine.
  !-------------------------------------------------------------------!

  subroutine marched(schemes, added, q, achieved, design_value)

    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:)
    real(dp), allocatable, intent(out) :: q(:)
    real(dp)           , intent(out) :: achieved
    real(dp), intent(in), optional :: design_value

    type(chain_block), allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    integer, allocatable :: first(:), last(:)
    real(dp), allocatable :: dt(:), t(:)
    real(dp) :: design
    integer :: k, n

    design = 0.0_dp
    if (present(design_value)) design = design_value

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    n = last(size(added))
    call partition(duration, n, dt, t)

    call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
         & uniform_grid(duration), design, initial_for(schemes(1) % scheme, n), &
         & chain, tower, dt, t, achieved)

    ! The chain keeps its state one block at a time, and a stage
    ! block keeps its instants between its stages, so the instants
    ! are gathered here into the one layout the checks below read.
    allocate(q(n * degrees))
    do k = 1, n
       q((k - 1) * degrees + 1:k * degrees) = instant_components(chain, k)
    end do

  end subroutine marched

  function bdf_of(order) result(held)

    integer, intent(in) :: order
    type(family_holder) :: held

    allocate(held % scheme, source=bdf_family(order))

  end function bdf_of

  function adams_of(order) result(held)

    integer, intent(in) :: order
    type(family_holder) :: held

    allocate(held % scheme, source=adams_family(order))

  end function adams_of

  !-------------------------------------------------------------------!
  ! Splitting a horizon must change nothing: the same rows are made,
  ! in the same places, from the same steps.
  !-------------------------------------------------------------------!

  subroutine splitting_changes_nothing()

    type(family_holder) :: whole(1), split(2)
    real(dp), allocatable :: q_whole(:), q_split(:)
    real(dp) :: achieved_whole, achieved_split

    whole(1) = bdf_of(2)
    split(1) = bdf_of(2)
    split(2) = bdf_of(2)

    call marched(whole, [40], q_whole, achieved_whole)
    call marched(split, [20, 20], q_split, achieved_split)

    write(*,'(a)')        ' bdf 2 over forty instants, in one block and in two'
    write(*,'(a,i0)')     '   unknowns, whole             ', size(q_whole)
    write(*,'(a,i0)')     '   unknowns, split             ', size(q_split)
    write(*,'(a,es11.2)') '   worst difference between them', maxval(abs(q_whole - q_split))
    write(*,'(a,es11.2)') '   residual, whole             ', achieved_whole
    write(*,'(a,es11.2)') '   residual, split             ', achieved_split

  end subroutine splitting_changes_nothing

  !-------------------------------------------------------------------!
  ! A backward difference handing over to an Adams quadrature.
  !-------------------------------------------------------------------!

  subroutine across_a_change_of_scheme()

    type(family_holder) :: schemes(2)
    real(dp), allocatable :: q(:), dt(:), t(:)
    real(dp) :: e(3), achieved
    integer :: level, added, k, n

    schemes(1) = bdf_of(2)
    schemes(2) = adams_of(3)

    do level = 1, 3
       added = 10 * 2 ** (level - 1)
       call marched(schemes, [added, added], q, achieved)

       n = 2 * added
       call partition(duration, n, dt, t)

       e(level) = 0.0_dp
       do k = 1, n
          e(level) = max(e(level), abs(q((k - 1) * degrees + 1) - exact(0, t(k))))
       end do
    end do

    write(*,'(a)')          ' '
    write(*,'(a)')          ' bdf 2 then adams-moulton 3, on the harmonic oscillator'
    write(*,'(a,3i11)')     '   instants                   ', [(20 * 2 ** (level - 1), level = 1, 3)]
    write(*,'(a,3es11.3)')  '   worst error                ', e
    write(*,'(a,22x,2f11.3)') '   ratio                    ', e(1:2) / e(2:3)
    write(*,'(a,es11.2)')   '   residual                   ', achieved

  end subroutine across_a_change_of_scheme

  !-------------------------------------------------------------------!
  ! The functional over a whole marched horizon.
  !-------------------------------------------------------------------!

  real(dp) function energy_of(q, n, design_value) result(f)

    real(dp), intent(in) :: q(:)
    integer , intent(in) :: n
    real(dp), intent(in) :: design_value

    type(stored_directed_graph) :: unknowns, instants
    type(stored_field) :: state, knobs
    real(dp), allocatable :: dt(:), t(:)

    call partition(duration, n, dt, t)

    unknowns = unknowns_graph(n, degrees)
    instants = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), size(q))
    knobs    = stored_field('design', unknowns % vertex_set(), n)
    call state % set_real_vector(q)
    call knobs % set_real_vector(spread(design_value, 1, n))

    f = functional_of(van_der_pol_energy(state_degree), instants, [state, knobs], dt)

  end function energy_of

  !-------------------------------------------------------------------!
  ! The sensitivity of that functional to the design, across the
  ! junction, by the tangent, by the adjoint, and by marching the
  ! whole horizon either side of the design.
  !-------------------------------------------------------------------!

  subroutine sensitivity_across_the_junction()

    real(dp), parameter :: delta = 1.0e-6_dp
    real(dp), parameter :: design = 1.0_dp
    integer , parameter :: added(2) = [10, 10]

    type(family_holder) :: schemes(2)
    type(chain_block) , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    type(chain_system), allocatable :: systems(:)
    real(dp), allocatable :: q(:), dt(:)
    real(dp) :: f, tangent, adjoint, differenced, achieved
    integer :: n

    schemes(1) = bdf_of(2)
    schemes(2) = adams_of(3)
    n = sum(added)

    call marched(schemes, added, q, achieved, design)
    f = energy_of(q, n, design)

    call chained(schemes, added, design, chain, tower, dt)
    call chain_systems(chain, tower, [one_functional(van_der_pol_energy(state_degree))], degrees, systems)

    tangent     = first_of(chain_by_tangent(chain, systems, degrees, design))
    adjoint     = first_of(chain_by_adjoint(chain, systems, degrees, design))
    differenced = differenced_energy(schemes, added, n, design, delta)

    write(*,'(a)')        ' '
    write(*,'(a)')        ' bdf 2 then adams-moulton 3, van der pol at a design of one'
    write(*,'(a,i0)')     '   blocks                      ', size(added)
    write(*,'(a,f16.10)') '   the functional             ', f
    write(*,'(a,f16.10)') '   sensitivity, tangent       ', tangent
    write(*,'(a,f16.10)') '   sensitivity, adjoint       ', adjoint
    write(*,'(a,f16.10)') '   sensitivity, differenced   ', differenced
    write(*,'(a,es16.2)') '   tangent against adjoint    ', abs(tangent - adjoint)
    write(*,'(a,es16.2)') '   tangent against difference ', abs(tangent - differenced)

  end subroutine sensitivity_across_the_junction

  !-------------------------------------------------------------------!
  ! The chain itself, and the steps it was built over.
  !-------------------------------------------------------------------!

  subroutine chained(schemes, added, design, chain, tower, dt)

    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:)
    real(dp)           , intent(in) :: design
    type(chain_block), allocatable, intent(out) :: chain(:)
    type(expansion)  , allocatable, intent(inout), target :: tower
    real(dp)         , allocatable, intent(out) :: dt(:)

    integer , allocatable :: first(:), last(:)
    real(dp), allocatable :: t(:)
    real(dp) :: achieved

    call horizon_bounds(schemes, added, degrees - 1, first, last)
    call partition(duration, last(size(added)), dt, t)

    call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
         & uniform_grid(duration), design, &
         & initial_for(schemes(1) % scheme, last(size(added))), &
         & chain, tower, dt, t, achieved)

  end subroutine chained

  !-------------------------------------------------------------------!
  ! The whole horizon marched either side of the design and the
  ! functional differenced, which shares no arithmetic with either
  ! sweep.
  !-------------------------------------------------------------------!

  real(dp) function differenced_energy(schemes, added, n, design, delta) result(d)

    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:), n
    real(dp)           , intent(in) :: design, delta

    real(dp), allocatable :: plus(:), minus(:)
    real(dp) :: achieved

    call marched(schemes, added, plus,  achieved, design + delta)
    call marched(schemes, added, minus, achieved, design - delta)

    d = (energy_of(plus, n, design + delta) - &
       & energy_of(minus, n, design - delta)) / (2.0_dp * delta)

  end function differenced_energy

end program marched_horizon
