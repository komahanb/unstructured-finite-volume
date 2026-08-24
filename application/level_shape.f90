! Two blocks built through level_storage and read back through the
! level view, one per family shape.
!
! A multistep block: slices hold their components directly. A
! multistage block: slices hold stages, and stages hold components.
! The traversal is the same procedure at every level and never asks
! which kind of level it is at; the two depths are what differ.
!
! Every traversal follows the spine once, through sequence_first and
! sequence_rest. Nothing is indexed by position.
!
! The last part shows what the consistency check refuses: two
! couplings of equal size, one beginning with this level's own
! members and one naming a member of no level, so that counting
! cannot separate them and only identity can.
program level_shape

  use graph_fractal, only : graph, branch, known_branch
  use view_sequence, only : sequence_empty, sequence_first, sequence_rest
  use view_level   , only : level_storage, level_is_leaf, level_num_members, &
       & level_members, level_couples, level_consistent

  implicit none

  integer, parameter :: max_instants     = 2
  integer, parameter :: max_stages       = 2
  integer, parameter :: max_state_degree = 1

  type(level_storage) :: store
  type(graph), pointer :: root
  integer :: multistep, multistage

  multistep  = one_block(one_multistep_slice)
  multistage = one_block(one_multistage_slice)

  write(*,'(a)') ' a multistep block: slices hold components'
  root => store % node(multistep)
  call show(root, 1)

  write(*,'(a)') ' '
  write(*,'(a)') ' a multistage block: slices hold stages'
  root => store % node(multistage)
  call show(root, 1)

  write(*,'(a)')    ' '
  write(*,'(a,i0)') ' nodes owned by the storage   ', store % num_nodes()

  call the_stranger_is_refused()

contains

  !-------------------------------------------------------------------!
  ! The indices of n members, each built by one call of make: none,
  ! or one member followed by the rest. The first is built before
  ! the rest so that the storage receives them in order.
  !-------------------------------------------------------------------!

  recursive function members_of(n, make) result(members)

    integer, intent(in) :: n
    interface
       integer function make()
       end function make
    end interface
    integer, allocatable :: members(:)

    integer :: first

    if (n == 0) then
       allocate(members(0))
       return
    end if

    first   = make()
    members = [first, members_of(n - 1, make)]

  end function members_of

  integer function one_leaf() result(leaf)

    leaf = store % assemble([integer ::], 0)

  end function one_leaf

  !-------------------------------------------------------------------!
  ! A level whose members couple: the coupling's carriers begin with
  ! those members.
  !-------------------------------------------------------------------!

  integer function coupled(members) result(level)

    integer, intent(in) :: members(:)

    level = store % assemble(members, store % assemble(members, 0))

  end function coupled

  !-------------------------------------------------------------------!
  ! The components of one bundle, degree 0 .. max_state_degree.
  !-------------------------------------------------------------------!

  function components() result(members)

    integer, allocatable :: members(:)

    members = members_of(max_state_degree + 1, one_leaf)

  end function components

  integer function one_multistep_slice() result(slice)

    slice = coupled(components())

  end function one_multistep_slice

  integer function one_stage() result(stage)

    stage = coupled(components())

  end function one_stage

  integer function one_multistage_slice() result(slice)

    slice = coupled(members_of(max_stages, one_stage))

  end function one_multistage_slice

  !-------------------------------------------------------------------!
  ! One block: instants 0 .. max_instants, each built by make, coupled.
  !-------------------------------------------------------------------!

  integer function one_block(make) result(block)

    interface
       integer function make()
       end function make
    end interface

    block = coupled(members_of(max_instants + 1, make))

  end function one_block

  !-------------------------------------------------------------------!
  ! Print one level, then each of its members in turn.
  !-------------------------------------------------------------------!

  recursive subroutine show(g, depth)

    type(graph), intent(in) :: g
    integer    , intent(in) :: depth

    if (level_is_leaf(g)) then
       write(*,'(a,a)') repeat('   ', depth), 'leaf'
       return
    end if

    write(*,'(a,a,i0,a,l1,a,l1)') repeat('   ', depth), 'members ', &
         & level_num_members(g), '   couples ', level_couples(g), &
         & '   consistent ', level_consistent(g)

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

  !-------------------------------------------------------------------!
  ! Two couplings of two carriers each. The first begins with this
  ! level's own members; the second names a member of no level.
  !-------------------------------------------------------------------!

  subroutine the_stranger_is_refused()

    type(graph), pointer :: subject_node, other
    integer, allocatable :: mine(:)
    integer :: stranger, alien, subject

    mine     = members_of(2, one_leaf)
    stranger = one_leaf()
    alien    = store % assemble([mine(1), stranger], 0)

    subject = coupled(mine)
    subject_node => store % node(subject)
    write(*,'(a)')    ' '
    write(*,'(a,l1)') ' coupling beginning with its own members is consistent ', &
         & level_consistent(subject_node)

    other => store % node(alien)
    subject_node % branch(2) = known_branch(other)
    write(*,'(a,i0)') ' the other coupling carries the same count             ', &
         & level_num_members(subject_node)
    write(*,'(a,l1)') ' and is refused by identity                            ', &
         & level_consistent(subject_node)

  end subroutine the_stranger_is_refused

end program level_shape
