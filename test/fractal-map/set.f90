!=====================================================================!
! PROTOTYPE . WHAT IS A FINITE SET IN THE FRACTAL SYSTEM
!
! Two designs are built and compared:
!
!   A  EXTENSIONAL GRAPH        the set graph references every member
!                               through a sequence
!   B  SEMANTIC GRAPH + EXTENT  one graph denotes the set; an external
!                               representation answers how its members
!                               are stored and enumerated
!
! Design A is built here at small n so its cost is not hypothetical:
! it needs one element graph and one sequence cell per member. Design
! B is built with a test-local extent that stores one integer.
!
! Nothing in src/ is modified. The extent below is deliberately NOT a
! member_set: the question under test is whether the identity can leave
! the set object entirely.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module set_prototype

  use fractal_graph, only : graph

  implicit none

  private
  public :: counted_extent, extent_of, bind_extent, set_equivalent

  !===================================================================!
  ! The counted representation: how many members, and the convention
  ! that they enumerate 1..n. O(1) storage, whatever n is.
  !
  ! It carries NO identity. Which set this describes is the question
  ! the graph answers.
  !===================================================================!

  type :: counted_extent

     integer, private :: n = 0

   contains

     procedure :: extent_size
     procedure :: member
     procedure :: has
     procedure :: local_index

  end type counted_extent

  interface counted_extent
     module procedure create_counted_extent
  end interface counted_extent

  !===================================================================!
  ! The map: set graph -> its extent, keyed on graph identity. The
  ! same shape as relational_binding, and for the same reason.
  !===================================================================!

  type :: extent_row
     type(graph), pointer  :: element => null()
     type(counted_extent)  :: value
  end type extent_row

  type(extent_row), allocatable, save :: rows(:)

contains

  type(counted_extent) function create_counted_extent(n) result(this)

    integer, intent(in) :: n

    this % n = max(0, n)

  end function create_counted_extent

  pure integer function extent_size(this)

    class(counted_extent), intent(in) :: this

    extent_size = this % n

  end function extent_size

  pure integer function member(this, position)

    class(counted_extent), intent(in) :: this
    integer              , intent(in) :: position

    member = position

  end function member

  pure logical function has(this, value)

    class(counted_extent), intent(in) :: this
    integer              , intent(in) :: value

    has = (value >= 1) .and. (value <= this % n)

  end function has

  pure integer function local_index(this, value)

    class(counted_extent), intent(in) :: this
    integer              , intent(in) :: value

    if (this % has(value)) then
       local_index = value
    else
       local_index = 0
    end if

  end function local_index

  !===================================================================!
  ! Binding and lookup, by graph identity.
  !===================================================================!

  subroutine bind_extent(s, value)

    type(graph)         , intent(in), target :: s
    type(counted_extent), intent(in)         :: value

    type(extent_row), allocatable :: grown(:)
    integer                       :: n

    if (.not. allocated(rows)) allocate(rows(0))
    n = size(rows)
    allocate(grown(n + 1))
    grown(1:n) = rows
    grown(n + 1) % element => s
    grown(n + 1) % value   =  value
    call move_alloc(grown, rows)

  end subroutine bind_extent

  type(counted_extent) function extent_of(s) result(value)

    type(graph), intent(in) :: s

    integer :: k

    if (allocated(rows)) then
       do k = 1, size(rows)
          if (rows(k) % element % same_as(s)) then
             value = rows(k) % value
             return
          end if
       end do
    end if

    error stop 'set_prototype: no extent is bound to that set'

  end function extent_of

  !===================================================================!
  ! Extensional equality: equality of MEMBERS, defined separately from
  ! identity and never confused with it.
  !===================================================================!

  logical function set_equivalent(a, b) result(equal)

    type(graph), intent(in) :: a, b

    type(counted_extent) :: ea, eb
    integer              :: k

    ea = extent_of(a)
    eb = extent_of(b)

    equal = ea % extent_size() .eq. eb % extent_size()
    if (.not. equal) return

    do k = 1, ea % extent_size()
       equal = equal .and. (ea % member(k) .eq. eb % member(k))
    end do

  end function set_equivalent

end module set_prototype

program set_map

  use fractal_graph , only : graph, null_branch, known_branch, &
       & GRAPH_NULL, GRAPH_KNOWN
  use graph_sequence_view, only : sequence_size, sequence_element
  use set_prototype , only : counted_extent, extent_of, bind_extent, &
       & set_equivalent

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "set prototype"

  !===================================================================!
  ! 1 . TWO DISTINCT EMPTY SETS.
  !
  ! Both are (NULL, NULL) under design A and both have extent 0 under
  ! design B. They are extensionally equal and NOT the same set: the
  ! empty set of cells is not the empty set of faces.
  !===================================================================!

  empty_block: block

    type(graph), target :: cells, faces

    call cells % declare(); call faces % declare()
    cells % branch(1) = null_branch(); cells % branch(2) = null_branch()
    faces % branch(1) = null_branch(); faces % branch(2) = null_branch()

    call bind_extent(cells, counted_extent(0))
    call bind_extent(faces, counted_extent(0))

    call check('1  two empty sets are extensionally equal', &
         & set_equivalent(cells, faces))
    call check('1  and are NOT the same set', .not. cells % same_as(faces))
    call check('1  identity survives an empty extension', &
         & cells % same_as(cells) .and. faces % same_as(faces))

  end block empty_block

  !===================================================================!
  ! 2 . DESIGN A, BUILT. The extensional set graph, at n = 4.
  !
  !     S.branch(1) = the member sequence
  !     S.branch(2) = NULL
  !
  ! It works, and it costs one element graph and one cell per member.
  !===================================================================!

  extensional_block: block

    type(graph), target  :: s, cell(4), elem(4)
    type(graph), pointer :: e
    integer              :: k
    logical              :: ok

    call s % declare()
    do k = 1, 4
       call cell(k) % declare(); call elem(k) % declare()
    end do

    do k = 1, 4
       cell(k) % branch(1) = known_branch(elem(k))
       if (k .lt. 4) cell(k) % branch(2) = known_branch(cell(k + 1))
    end do

    s % branch(1) = known_branch(cell(1))
    s % branch(2) = null_branch()

    call check('2  design A: the set graph answers |S| = 4', &
         & sequence_size(s % branch(1)) .eq. 4)

    ok = .true.
    do k = 1, 4
       e => sequence_element(s % branch(1), k)
       ok = ok .and. e % same_as(elem(k))
    end do
    call check('2  and every member is reachable by identity', ok)
    call check('2  the cost is 2n graph objects, here 8', &
         & size(cell) + size(elem) .eq. 8)

  end block extensional_block

  !===================================================================!
  ! 3 . DESIGN B, AT A SCALE DESIGN A CANNOT REACH.
  !
  ! One graph, one extent holding one integer. The set has 10^9
  ! members and the program allocates no member object at all.
  !===================================================================!

  counted_block: block

    type(graph), target  :: cells
    type(counted_extent) :: e

    call cells % declare()
    cells % branch(1) = null_branch()      ! the extension is not here
    cells % branch(2) = null_branch()

    call bind_extent(cells, counted_extent(1000000000))
    e = extent_of(cells)

    call check('3  design B: |cells| = 10^9 with no member object', &
         & e % extent_size() .eq. 1000000000)
    call check('3  membership is one comparison, not a search', &
         & e % has(999999999) .and. .not. e % has(1000000001))
    call check('3  and position is the representation''s numbering', &
         & e % local_index(7) .eq. 7 .and. e % member(7) .eq. 7)
    call check('3  the graph carries no extension: both branches NULL', &
         & cells % branch(1) % status() .eq. GRAPH_NULL .and. &
         &  cells % branch(2) % status() .eq. GRAPH_NULL)

  end block counted_block

  !===================================================================!
  ! 4 . SUBSET WITHOUT A SUBTYPE.
  !
  ! S and A are both set graphs. Inclusion is a PROPERTY between them,
  ! decided from the two extents - no third type, no inheritance, no
  ! ambient component.
  !
  ! Two coordinate systems must be kept apart:
  !
  !     the member VALUE     shared with the ambient
  !     the local POSITION   private to each set's representation
  !
  ! This is the split the current subset_set makes with roll: roll(k)
  ! is the ambient value, k is the local position.
  !===================================================================!

  subset_block: block

    type(graph), target  :: ambient, part
    type(counted_extent) :: ea
    integer              :: roll(3) = [2, 5, 6]
    integer              :: k
    logical              :: included, positions_differ

    call ambient % declare(); call part % declare()
    call bind_extent(ambient, counted_extent(8))
    ea = extent_of(ambient)

    ! part's extent is its roll; the counted extent cannot express it,
    ! so the roll stands here as the explicit-index representation.
    included = .true.
    do k = 1, size(roll)
       included = included .and. ea % has(roll(k))
    end do

    call check('4  S subset A is decided from the two extents', included)
    call check('4  and needs no subtype and no ambient component', &
         & .not. part % same_as(ambient))

    ! Local position differs from member value; the values agree with
    ! the ambient, which is why inclusion needs no translation map.
    positions_differ = .false.
    do k = 1, size(roll)
       if (roll(k) .ne. k) positions_differ = .true.
    end do
    call check('4  local position /= member value, and only the value travels', &
         & positions_differ .and. ea % local_index(roll(2)) .eq. 5)

  end block subset_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'set: a proposition failed'
  end if

contains

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program set_map
