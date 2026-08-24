! What each level carries beyond its structure.
!
! A node is a graph and holds only its two branches. Everything else
! about it is kept in a map keyed on its identity, and there are
! three: what it is called, what values it holds, and how many
! members the set it denotes has. This builds a complete tower and
! attaches all three at every level, then reads it back with one
! traversal that never asks which level it is at.
!
!    expansion    label  the physics and what it is differentiated in
!                 value  the design, known; one term per sweep, not
!                        yet known
!                 set    the extent of the design
!    sweep        label  which derivative of the functional it makes
!    horizon      label  the duration it covers
!    block        label  the family and order marching it
!                 value  the steps, known
!    slice        label  the instant it sits at
!    component    label  the degree it holds
!                 value  the freedoms, known at the initial instants
!                        and not yet known after them
!                 set    the extent of the freedoms
!
! The value map separates three states: no row at all, a row whose
! value is not yet determined, and a row that holds one. That is the
! same distinction the kernel draws between a branch that is NULL,
! UNKNOWN and KNOWN, and it is what a march turns over as it goes:
! every component after the initial instants starts not yet known and
! becomes known when its block is solved.
!
! The map enforces that order. A row is attached before it can be
! marked, so a value can only become known where one was expected,
! and a graph that was never given a row reads as holding no value
! rather than as holding zero.
!
! The couplings are left out here. They are columns three and four
! and are built in coupling_relation.f90; this is about what hangs on
! an identity, not about what reads what.
program level_maps

  use iso_fortran_env       , only : dp => REAL64
  use graph_fractal         , only : graph, branch
  use view_sequence         , only : sequence_empty, sequence_first, sequence_rest
  use view_level            , only : level_storage, level_is_leaf, &
       & level_num_members, level_members
  use map_value             , only : value_map, VALUE_UNATTACHED, VALUE_UNKNOWN, &
       & VALUE_KNOWN
  use map_label             , only : label_map
  use map_set               , only : set_map
  use map_set_representation, only : counted_set_representation

  implicit none

  integer , parameter :: max_derivative_degree = 1     ! the primal and one tangent
  integer , parameter :: max_state_degree      = 2     ! q, q', q''
  integer , parameter :: num_instants          = 3     ! per block
  integer , parameter :: num_freedoms          = 1     ! an ordinary differential equation
  real(dp), parameter :: duration              = 1.15_dp

  type(level_storage) :: store
  type(value_map)     :: values
  type(label_map)     :: labels
  type(set_map)       :: extents

  integer, allocatable :: sweeps(:)
  integer :: expansion, s

  sweeps = [(one_sweep(s), s = 0, max_derivative_degree)]
  expansion = store % assemble(sweeps, 0)

  call labels % bind(store % node(expansion), 'expansion of the van der pol functional in nu')
  call extents % bind(store % node(expansion), counted_set_representation(1))
  call attach_known(store % node(expansion), [1.0_dp])

  write(*,'(a)') ' the tower, and what each level carries'
  call show(store % node(expansion), 0)

  write(*,'(a)')    ' '
  write(*,'(a,i0)') ' nodes owned by the storage      ', store % num_nodes()
  write(*,'(a,i0)') ' components not yet known        ', not_yet_known(store % node(expansion))

  call determine(store % node(expansion))
  write(*,'(a,i0)') ' after every block is solved     ', not_yet_known(store % node(expansion))

contains

  !-------------------------------------------------------------------!
  ! A row attached and marked in one step, for a value that is known
  ! as soon as the node exists.
  !-------------------------------------------------------------------!

  subroutine attach_known(g, x)

    type(graph), intent(in) :: g
    real(dp)   , intent(in) :: x(:)

    call values % attach_unknown(g)
    call values % mark_known(g, x)

  end subroutine attach_known

  !-------------------------------------------------------------------!
  ! One component: a leaf holding its freedoms. The first two
  ! instants of the first block carry the initial conditions and are
  ! known; everything after them waits on a march.
  !-------------------------------------------------------------------!

  integer function one_component(degree, instant, block_index) result(at)

    integer, intent(in) :: degree, instant, block_index

    character(len=1) :: d

    at = store % assemble([integer ::], 0)
    write(d,'(i1)') degree

    call labels  % bind(store % node(at), 'component of degree ' // d)
    call extents % bind(store % node(at), counted_set_representation(num_freedoms))

    if (block_index == 1 .and. instant <= 2) then
       call attach_known(store % node(at), spread(0.0_dp, 1, num_freedoms))
    else
       call values % attach_unknown(store % node(at))
    end if

  end function one_component

  !-------------------------------------------------------------------!
  ! One slice: the components of every degree at one instant.
  !-------------------------------------------------------------------!

  integer function one_slice(instant, block_index) result(at)

    integer, intent(in) :: instant, block_index

    character(len=2) :: k
    integer :: d

    at = store % assemble([(one_component(d, instant, block_index), &
         & d = 0, max_state_degree)], 0)

    write(k,'(i2)') instant
    call labels % bind(store % node(at), 'slice at instant' // k)

  end function one_slice

  !-------------------------------------------------------------------!
  ! One block: the instants it marches, the family that marches them,
  ! and the steps it takes.
  !-------------------------------------------------------------------!

  integer function one_block(block_index, family) result(at)

    integer         , intent(in) :: block_index
    character(len=*), intent(in) :: family

    real(dp) :: steps(num_instants)
    integer  :: k

    at = store % assemble([(one_slice(k, block_index), k = 1, num_instants)], 0)

    call labels % bind(store % node(at), family)

    steps    = duration / real(2 * num_instants, dp)
    steps(1) = 0.0_dp
    call attach_known(store % node(at), steps)

  end function one_block

  !-------------------------------------------------------------------!
  ! One horizon: the blocks that partition the duration.
  !-------------------------------------------------------------------!

  integer function one_horizon() result(at)

    character(len=8) :: t

    at = store % assemble([one_block(1, 'bdf of order 2'), &
         &                 one_block(2, 'adams-moulton of order 3')], 0)

    write(t,'(f8.4)') duration
    call labels % bind(store % node(at), 'horizon of duration' // t)

  end function one_horizon

  !-------------------------------------------------------------------!
  ! One sweep: a traversal of its own horizon. Sweep zero is the
  ! primal and the sweeps above it are the tangents, so each owns a
  ! horizon of its own and the value rows of one never stand for
  ! another's.
  !-------------------------------------------------------------------!

  integer function one_sweep(sensitivity) result(at)

    integer, intent(in) :: sensitivity

    character(len=1) :: s

    at = store % assemble([one_horizon()], 0)

    write(s,'(i1)') sensitivity
    if (sensitivity == 0) then
       call labels % bind(store % node(at), 'sweep 0, the functional itself')
    else
       call labels % bind(store % node(at), 'sweep ' // s // ', derivative ' // s // ' in nu')
    end if
    call values % attach_unknown(store % node(at))

  end function one_sweep

  !-------------------------------------------------------------------!
  ! One traversal, printing what the three maps hold. It never asks
  ! which level it is at: the label is data, and a map that holds
  ! nothing for a node says so.
  !-------------------------------------------------------------------!

  recursive subroutine show(g, depth)

    type(graph), intent(in) :: g
    integer    , intent(in) :: depth

    character(len=:), allocatable :: name

    name = ' '
    if (labels % labelled(g)) name = labels % label_of(g)

    write(*,'(a,a,a,a)') repeat('   ', depth + 1), name, &
         & '   [' // status_name(values % status_of(g)) // ']', extent_of(g)

    if (level_is_leaf(g)) return
    call show_each(level_members(g), depth + 1)

  end subroutine show

  recursive subroutine show_each(members, depth)

    type(branch), intent(in) :: members
    integer     , intent(in) :: depth

    type(graph), pointer :: first

    if (sequence_empty(members)) return

    first => sequence_first(members)
    call show(first, depth)
    call show_each(sequence_rest(members), depth)

  end subroutine show_each

  pure function status_name(status) result(name)

    integer, intent(in) :: status
    character(len=:), allocatable :: name

    select case (status)
    case (VALUE_KNOWN)
       name = 'known'
    case (VALUE_UNKNOWN)
       name = 'not yet known'
    case (VALUE_UNATTACHED)
       name = 'no value'
    case default
       error stop 'level_maps: a value status is one of the three'
    end select

  end function status_name

  function extent_of(g) result(text)

    type(graph), intent(in) :: g
    character(len=:), allocatable :: text

    character(len=3) :: n

    text = ''
    if (.not. extents % describes(g)) return

    write(n,'(i3)') extents % num_members_of(g)
    text = '   extent' // n

  end function extent_of

  !-------------------------------------------------------------------!
  ! What a march turns over: every component that was waiting is
  ! marked with the value its block determined. The rows were
  ! attached when the tower was built, so nothing new is keyed here.
  !-------------------------------------------------------------------!

  recursive subroutine determine(g)

    type(graph), intent(in) :: g

    if (level_is_leaf(g)) then
       if (values % status_of(g) == VALUE_UNKNOWN) then
          call values % mark_known(g, spread(1.0_dp, 1, num_freedoms))
       end if
       return
    end if

    call determine_each(level_members(g))

  end subroutine determine

  recursive subroutine determine_each(members)

    type(branch), intent(in) :: members

    type(graph), pointer :: first

    if (sequence_empty(members)) return

    first => sequence_first(members)
    call determine(first)
    call determine_each(sequence_rest(members))

  end subroutine determine_each

  !-------------------------------------------------------------------!
  ! How many components are still waiting on a march.
  !-------------------------------------------------------------------!

  recursive integer function not_yet_known(g) result(n)

    type(graph), intent(in) :: g

    n = 0
    if (level_is_leaf(g)) then
       if (values % status_of(g) == VALUE_UNKNOWN) n = 1
       return
    end if

    n = counted(level_members(g))

  end function not_yet_known

  recursive integer function counted(members) result(n)

    type(branch), intent(in) :: members

    type(graph), pointer :: first

    n = 0
    if (sequence_empty(members)) return

    first => sequence_first(members)
    n = not_yet_known(first) + counted(sequence_rest(members))

  end function counted

end program level_maps
