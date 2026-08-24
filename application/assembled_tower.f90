! A time integrator assembled from its slots.
!
! The slots are a physics, one family per block, the instants each
! block covers, a grid and a design. Nothing else about the problem
! is given, and the two blocks here are marched by different families
! to show that a horizon does not have to be homogeneous.
!
! The tower is printed by one traversal that never asks which level
! it is at, and then checked from the outside: every level consistent,
! every coupling relationally valid.
!
! The last part is the grid's own partial. The steps are normalised
! so that they sum to the duration, so a design that changes one step
! changes them all, and the partial has to carry that. It is compared
! against a central difference, and the sum of the steps is printed
! before and after to show the constraint holding.
program assembled_tower

  use iso_fortran_env       , only : dp => REAL64
  use graph_fractal         , only : graph, branch
  use view_sequence         , only : sequence_empty, sequence_first, sequence_rest
  use view_level            , only : level_is_leaf, level_members, level_couples, &
       & level_coupling
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_action      , only : variation
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : random_grid, designed_grid, uniform_grid
  use operation_family_dirk , only : crouzeix_two_stage
  use physics_vanderpol     , only : van_der_pol
  use map_value             , only : VALUE_KNOWN, VALUE_UNKNOWN, VALUE_UNATTACHED
  use gti_expansion         , only : expansion, family_holder

  implicit none

  real(dp), parameter :: duration = 7.0_dp
  integer , parameter :: seed = 20260824

  type(expansion) :: one, two
  type(family_holder) :: schemes(2)

  allocate(schemes(1) % scheme, source=bdf_family(2))
  allocate(schemes(2) % scheme, source=adams_family(3))

  call one % build(van_der_pol(2), schemes, [5, 5], random_grid(duration, seed), 0, &
       & [real(dp) ::])

  write(*,'(a)') ' a horizon of two blocks, marched by different families'
  call show(one, one % node(one % root()), 0)

  write(*,'(a)')      ' '
  write(*,'(a,i0)')   ' nodes owned                        ', one % num_nodes()
  write(*,'(a,l1)')   ' every level and coupling valid     ', &
       & one % consistent(one % node(one % root()))

  call two % build(van_der_pol(2), schemes, [5, 5], random_grid(duration, seed), 1, &
       & [real(dp) ::])
  write(*,'(a,i0)')   ' nodes owned with one tangent sweep ', two % num_nodes()
  write(*,'(a,l1)')   ' every level and coupling valid     ', &
       & two % consistent(two % node(two % root()))

  call grid_partials()
  call stage_block()

contains

  !-------------------------------------------------------------------!
  ! One level, then each of its members. A coupling is printed under
  ! the level that holds it, and the numbers it carries are its own
  ! value.
  !-------------------------------------------------------------------!

  recursive subroutine show(tower, g, depth)

    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: g
    integer        , intent(in) :: depth

    type(graph), pointer :: coupling

    write(*,'(a,a,a,a)') repeat('   ', depth + 1), tower % label_of(g), &
         & status(tower, g), extent(tower, g)

    if (level_couples(g)) then
       coupling => level_coupling(g)
       write(*,'(a,a,a,a)') repeat('   ', depth + 2), tower % label_of(coupling), &
            & status(tower, coupling), extent(tower, coupling)
    end if

    if (level_is_leaf(g)) return
    call show_each(tower, level_members(g), depth + 1)

  end subroutine show

  recursive subroutine show_each(tower, members, depth)

    type(expansion), intent(in) :: tower
    type(branch)   , intent(in) :: members
    integer        , intent(in) :: depth

    type(graph), pointer :: first

    if (sequence_empty(members)) return

    first => sequence_first(members)
    call show(tower, first, depth)
    call show_each(tower, sequence_rest(members), depth)

  end subroutine show_each

  function status(tower, g) result(text)

    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: g
    character(len=:), allocatable :: text

    real(dp), allocatable :: x(:)

    select case (tower % status_of(g))
    case (VALUE_KNOWN)
       call tower % value_of(g, x)
       text = '   holds ' // count_of(size(x))
    case (VALUE_UNKNOWN)
       text = '   not yet known'
    case (VALUE_UNATTACHED)
       text = ''
    case default
       error stop 'assembled_tower: a value status is one of the three'
    end select

  end function status

  function extent(tower, g) result(text)

    type(expansion), intent(in) :: tower
    type(graph)    , intent(in) :: g
    character(len=:), allocatable :: text

    text = ''
    if (tower % extent_of(g) > 0) text = '   extent ' // count_of(tower % extent_of(g))

  end function extent

  function count_of(n) result(text)

    integer, intent(in) :: n
    character(len=:), allocatable :: text

    character(len=12) :: buffer

    write(buffer,'(i0)') n
    text = trim(buffer)

  end function count_of

  !-------------------------------------------------------------------!
  ! The grid's partial in the design, against a central difference,
  ! with the sum of the steps printed on either side.
  !-------------------------------------------------------------------!

  subroutine grid_partials()

    integer , parameter :: num_instants = 6
    real(dp), parameter :: delta = 1.0e-6_dp

    type(designed_grid) :: steps
    type(stored_directed_graph) :: instants
    type(stored_field) :: knobs, direction
    class(field), allocatable :: out
    real(dp) :: design(num_instants - 1), v(num_instants - 1)
    real(dp), allocatable :: dt(:), exact(:), plus(:), minus(:)

    design = [1.0_dp, 2.0_dp, 1.5_dp, 0.5_dp, 3.0_dp]
    v      = 0.0_dp
    v(2)   = 1.0_dp

    steps    = designed_grid(duration)
    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    knobs     = stored_field('design', instants % vertex_set(), size(design))
    direction = stored_field('v', instants % vertex_set(), size(design))
    call knobs     % set_real_vector(design)
    call direction % set_real_vector(v)

    call steps % apply(instants, [knobs], out)
    call out % real_vector(dt)

    call steps % partial_action(instants, [knobs], &
         & [variation(steps % argument(1), direction)], out)
    call out % real_vector(exact)

    call differenced(steps, instants, knobs, design, v, delta, plus, minus)

    call show_grid(dt, exact, (plus - minus) / (2.0_dp * delta))

  end subroutine grid_partials

  subroutine show_grid(dt, exact, differenced)

    real(dp), intent(in) :: dt(:), exact(:), differenced(:)

    write(*,'(a)')          ' '
    write(*,'(a)')          ' a designed grid of five steps'
    write(*,'(a,6f10.5)')   '   steps                       ', dt
    write(*,'(a,f10.5)')    '   their sum                   ', sum(dt)
    write(*,'(a,6f10.5)')   '   partial in design 2         ', exact
    write(*,'(a,6f10.5)')   '   central difference          ', differenced
    write(*,'(a,es10.2)')   '   the sum is unchanged along it', sum(exact)

  end subroutine show_grid

  !-------------------------------------------------------------------!
  ! The steps either side of the design, the design left as it was
  ! found.
  !-------------------------------------------------------------------!

  subroutine differenced(steps, instants, knobs, design, v, delta, plus, minus)

    type(designed_grid)        , intent(in)    :: steps
    type(stored_directed_graph), intent(in)    :: instants
    type(stored_field)         , intent(inout) :: knobs
    real(dp)                   , intent(in)    :: design(:), v(:), delta
    real(dp), allocatable      , intent(out)   :: plus(:), minus(:)

    class(field), allocatable :: out

    call knobs % set_real_vector(design + delta * v)
    call steps % apply(instants, [knobs], out)
    call out % real_vector(plus)

    call knobs % set_real_vector(design - delta * v)
    call steps % apply(instants, [knobs], out)
    call out % real_vector(minus)

    call knobs % set_real_vector(design)

  end subroutine differenced

  !-------------------------------------------------------------------!
  ! A block marched by stages. Its slices hold the stages of the step
  ! arriving at them and then the instant itself, and the block's
  ! first instant holds no stages because no step arrives there.
  !
  ! The weights of one step are printed against the tableau they come
  ! from: the rows below the highest degree carry a step, so they are
  ! the step times an entry, and the recovery of the highest degree
  ! carries none and is the entry itself.
  !-------------------------------------------------------------------!

  subroutine stage_block()

    integer , parameter :: num_instants = 3
    real(dp), parameter :: gamma = (3.0_dp + sqrt(3.0_dp)) / 6.0_dp

    type(expansion) :: staged
    type(family_holder) :: schemes(1)
    type(graph), pointer :: g, coupling
    real(dp), allocatable :: w(:)
    real(dp) :: step

    allocate(schemes(1) % scheme, source=crouzeix_two_stage())
    call staged % build(van_der_pol(2), schemes, [num_instants], &
         & uniform_grid(duration), 0, [real(dp) ::])

    write(*,'(a)') ' '
    write(*,'(a)') ' a block marched by a two-stage tableau'
    call show(staged, staged % node(staged % root()), 0)

    step = duration / real(num_instants - 1, dp)
    g => second_slice(staged)
    coupling => level_coupling(g)
    call staged % value_of(coupling, w)

    write(*,'(a)')        ' '
    write(*,'(a,l1)')     ' every level and coupling valid     ', &
         & staged % consistent(staged % node(staged % root()))
    write(*,'(a,i0)')     ' weights on one step                ', size(w)
    write(*,'(a,5f10.5)') '   value row, a_11 a_21 a_22 b_1 b_2', w(1:5)
    write(*,'(a,5f10.5)') '   step times the tableau           ', &
         & step * [gamma, 1.0_dp - 2.0_dp * gamma, gamma, 0.5_dp, 0.5_dp]
    write(*,'(a,2f10.5)') '   recovery of the highest degree   ', w(size(w) - 1:)
    write(*,'(a,2f10.5)') '   the tableau weights themselves   ', [0.5_dp, 0.5_dp]

  end subroutine stage_block

  !-------------------------------------------------------------------!
  ! The second slice of the only block of the only sweep, which is
  ! the first one a step arrives at.
  !-------------------------------------------------------------------!

  function second_slice(tower) result(g)

    type(expansion), intent(in) :: tower
    type(graph), pointer :: g

    type(branch) :: members

    members = level_members(tower % node(tower % root()))   ! sweeps
    g => sequence_first(members)
    g => sequence_first(level_members(g))                   ! the horizon's blocks
    g => sequence_first(level_members(g))                   ! the block
    g => sequence_first(sequence_rest(level_members(g)))    ! its second slice

  end function second_slice

end program assembled_tower
