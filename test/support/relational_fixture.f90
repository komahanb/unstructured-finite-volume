!=====================================================================!
! RELATIONAL FIXTURE . TEST SCAFFOLDING
!
! Builds the pair
!
!     type(graph)              read as (S, P)
!     type(relational_binding) element graph -> legacy object
!
! from a legacy relational_graph, so that a suite which already
! constructs its member sets and relations need not construct them
! twice while graph_profile is cut over.
!
! MIGRATION SCAFFOLDING, NOT PRODUCTION. This module reads
! relational_graph on purpose and dies with it: when the remaining
! suites are retargeted by capability and graph_structure is deleted,
! this file goes with it. It is under test/ and is never compiled into
! libufvm.
!
! STORAGE. A graph references its cells and a view borrows pointers
! into a binding, so both must outlive what reads them. The fixture
! owns the graphs and the bindings and never moves them: each call
! allocates its own, records it, and reallocates nothing afterwards.
! Callers therefore receive POINTERS and may reuse the same two
! variables at every site, because the objects behind earlier ones
! survive. A fixture must outlive every view built from it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relational_fixture

  use fractal_graph        , only : graph, null_branch, known_branch
  use graph_carrier        , only : member_set
  use graph_relation       , only : relation
  use graph_structure      , only : relational_graph
  use graph_relational_view, only : relational_binding

  implicit none

  private
  public :: fractal_fixture

  !===================================================================!
  ! One allocation of cell or element graphs. Recorded so the fixture
  ! can free it; never reallocated, so nothing it holds ever moves.
  !===================================================================!

  type :: storage_block
     type(graph), pointer :: cells(:) => null()
  end type storage_block

  type :: binding_block
     type(relational_binding), pointer :: object => null()
  end type binding_block

  type :: fractal_fixture

     type(storage_block), allocatable, private :: blocks(:)
     type(binding_block), allocatable, private :: bindings(:)

   contains

     procedure :: to_fractal
     procedure, private :: reserve

     final :: release_fixture

  end type fractal_fixture

contains

  !===================================================================!
  ! Allocate n declared graphs, record the block, and answer it.
  !===================================================================!

  function reserve(this, n) result(cells)

    class(fractal_fixture), intent(inout) :: this
    integer               , intent(in)    :: n
    type(graph), pointer                  :: cells(:)

    type(storage_block), allocatable :: grown(:)
    integer                          :: m, k

    allocate(cells(max(n, 1)))
    do k = 1, n
       call cells(k) % declare()
    end do

    if (.not. allocated(this % blocks)) allocate(this % blocks(0))
    m = size(this % blocks)
    allocate(grown(m + 1))
    grown(1:m) = this % blocks
    grown(m + 1) % cells => cells
    call move_alloc(grown, this % blocks)

  end function reserve

  !===================================================================!
  ! Read a relational_graph and answer the same structure as a graph
  ! plus the binding that resolves its elements.
  !
  !     branch(1) = the sequence of member sets
  !     branch(2) = the sequence of relations
  !===================================================================!

  subroutine to_fractal(this, source, g, b)

    class(fractal_fixture)           , intent(inout) :: this
    class(relational_graph)          , intent(in)    :: source
    type(graph)             , pointer, intent(out)   :: g
    type(relational_binding), pointer, intent(out)   :: b

    type(graph), pointer             :: scell(:), selem(:), rcell(:), relem(:)
    type(graph), pointer             :: holder(:)
    type(binding_block), allocatable :: grown(:)
    integer                          :: ns, nr, k, m

    ns = source % num_member_sets()
    nr = source % num_relations()

    scell => this % reserve(ns)
    selem => this % reserve(ns)
    rcell => this % reserve(nr)
    relem => this % reserve(nr)

    holder => this % reserve(1)
    g => holder(1)

    allocate(b)
    if (.not. allocated(this % bindings)) allocate(this % bindings(0))
    m = size(this % bindings)
    allocate(grown(m + 1))
    grown(1:m) = this % bindings
    grown(m + 1) % object => b
    call move_alloc(grown, this % bindings)

    do k = 1, ns
       call b % bind_set(selem(k), source % member_set_at(k))
       scell(k) % branch(1) = known_branch(selem(k))
       if (k .lt. ns) then
          scell(k) % branch(2) = known_branch(scell(k + 1))
       else
          scell(k) % branch(2) = null_branch()
       end if
    end do

    do k = 1, nr
       call b % bind_relation(relem(k), source % relation_at(k))
       rcell(k) % branch(1) = known_branch(relem(k))
       if (k .lt. nr) then
          rcell(k) % branch(2) = known_branch(rcell(k + 1))
       else
          rcell(k) % branch(2) = null_branch()
       end if
    end do

    if (ns .gt. 0) then
       g % branch(1) = known_branch(scell(1))
    else
       g % branch(1) = null_branch()
    end if

    if (nr .gt. 0) then
       g % branch(2) = known_branch(rcell(1))
    else
       g % branch(2) = null_branch()
    end if

  end subroutine to_fractal

  subroutine release_fixture(this)

    type(fractal_fixture), intent(inout) :: this

    integer :: k

    if (allocated(this % bindings)) then
       do k = 1, size(this % bindings)
          if (associated(this % bindings(k) % object)) then
             deallocate(this % bindings(k) % object)
          end if
       end do
    end if

    if (allocated(this % blocks)) then
       do k = 1, size(this % blocks)
          if (associated(this % blocks(k) % cells)) then
             deallocate(this % blocks(k) % cells)
          end if
       end do
    end if

  end subroutine release_fixture

end module relational_fixture
