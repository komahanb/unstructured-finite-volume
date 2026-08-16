!=====================================================================!
! RELATIONAL VIEW
!
! One view of a graph, in which the two branches are interpreted as the
! pair (S, P):
!
!     branch(1) = the sequence of member sets
!     branch(2) = the sequence of relations
!
! (S, P) is a view. The same graph remains readable under any other,
! and the kernel learns none of these words.
!
! THE BINDING. A branch references a GRAPH, never an arbitrary object.
! The member sets and relations of this repository are not yet graphs,
! so each element graph stands for one of them and a binding maps the
! element to the legacy object it denotes. The binding OWNS its
! objects, because a borrowed pointer handed to a caller outlives
! whatever lent it.
!
! THE STORAGE LAW.
!
!     relational_binding is not assignable; bind_* preserves every
!     outstanding object pointer until the binding is destroyed.
!
! That law is not free. A row holds a POINTER to an individually
! allocated object, never the object itself: when the row array grows,
! the rows are copied and the objects do not move. Holding the object
! in the row - as an allocatable component - was measured and rejected,
! because growth relocates the array and every borrowed pointer then
! reads freed storage: silently wrong first, fatal next. See
! test/graph-relational/lifetime.f90, which holds this law on every run.
!
! ASSIGNMENT IS REFUSED, at run time, because no Fortran mechanism
! prohibits it at compile time; four were compiled and measured in
! test/graph-relational/fortran-assignment. Extension and replacement
! are different operations: bind_* extends a binding and preserves what
! it lent, and there is no replacement that can.
!
! The refusing procedure takes its left-hand side INTENT(INOUT), not
! INTENT(OUT), because an INTENT(OUT) dummy of a finalizable type is
! finalized on entry - which would destroy the lender before refusing
! to destroy it.
!
! Because a row holds a pointer rather than the object, the pointer
! this module returns does not point into the binding. The binding
! therefore needs no TARGET attribute at any call site.
!
! The binding is storage keyed on identity, not ontology. It is where
! the wrappers that used to sit inside relational_graph belong: they
! existed only because a Fortran array carries one dynamic type, and
! that was always a storage fact.
!
! Sequence behaviour is delegated to graph_sequence_view. Nothing here
! traverses a spine.
!
! TWO FAILURES, STRUCTURALLY APART.
!
!     malformed sequence   refused, by graph_sequence_view
!     relationally invalid answered .false. by relational_valid
!
! A view over an existing graph reports invalidity; it does not refuse
! a graph it did not construct. The old create_graph refused at
! construction because it was a constructor; this is not one.
!
! THE VALIDITY LAW, unchanged in content:
!
!     G is relationally valid iff every domain of every relation in the
!     relation sequence is a member set in the member-set sequence.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_relational_view

  use fractal_graph      , only : graph
  use graph_carrier      , only : member_set
  use graph_relation     , only : relation
  use graph_sequence_view, only : sequence_size, sequence_element, &
       & sequence_contains

  implicit none

  private
  public :: relational_binding
  public :: num_member_sets, member_set_at
  public :: num_relations, relation_at
  public :: holds_set, relational_valid

  !===================================================================!
  ! Owned storage. One row per bound element.
  !===================================================================!

  type :: bound_set
     type(graph)      , pointer :: element => null()
     class(member_set), pointer :: object  => null()
  end type bound_set

  type :: bound_relation
     type(graph)    , pointer :: element => null()
     class(relation), pointer :: object  => null()
  end type bound_relation

  type :: relational_binding

     type(bound_set)     , allocatable, private :: sets(:)
     type(bound_relation), allocatable, private :: relations(:)

   contains

     procedure :: bind_set
     procedure :: bind_relation
     procedure :: set_for
     procedure :: relation_for

     procedure, private :: refuse_assignment
     generic :: assignment(=) => refuse_assignment

     final :: release_binding

  end type relational_binding

contains

  !===================================================================!
  ! Binding. The object is copied into owned storage; the element is
  ! referenced.
  !===================================================================!

  subroutine bind_set(this, element, object)

    class(relational_binding), intent(inout)        :: this
    type(graph)              , intent(in) , target  :: element
    class(member_set)        , intent(in)           :: object

    type(bound_set), allocatable :: grown(:)
    integer                      :: n

    if (.not. allocated(this % sets)) allocate(this % sets(0))
    n = size(this % sets)
    allocate(grown(n + 1))
    grown(1:n) = this % sets
    grown(n + 1) % element => element
    allocate(grown(n + 1) % object, source=object)
    call move_alloc(grown, this % sets)

  end subroutine bind_set

  subroutine bind_relation(this, element, object)

    class(relational_binding), intent(inout)       :: this
    type(graph)              , intent(in), target  :: element
    class(relation)          , intent(in)          :: object

    type(bound_relation), allocatable :: grown(:)
    integer                           :: n

    if (.not. allocated(this % relations)) allocate(this % relations(0))
    n = size(this % relations)
    allocate(grown(n + 1))
    grown(1:n) = this % relations
    grown(n + 1) % element => element
    allocate(grown(n + 1) % object, source=object)
    call move_alloc(grown, this % relations)

  end subroutine bind_relation

  !===================================================================!
  ! Lookup by element identity, answering a reference into owned
  ! storage.
  !===================================================================!

  function set_for(this, element) result(s)

    class(relational_binding), intent(in) :: this
    type(graph)              , intent(in) :: element
    class(member_set), pointer            :: s

    integer :: k

    if (allocated(this % sets)) then
       do k = 1, size(this % sets)
          if (this % sets(k) % element % same_as(element)) then
             s => this % sets(k) % object
             return
          end if
       end do
    end if

    error stop 'graph_relational_view: no member set is bound to that element'

  end function set_for

  function relation_for(this, element) result(r)

    class(relational_binding), intent(in) :: this
    type(graph)              , intent(in) :: element
    class(relation), pointer              :: r

    integer :: k

    if (allocated(this % relations)) then
       do k = 1, size(this % relations)
          if (this % relations(k) % element % same_as(element)) then
             r => this % relations(k) % object
             return
          end if
       end do
    end if

    error stop 'graph_relational_view: no relation is bound to that element'

  end function relation_for

  !===================================================================!
  ! Refusal. Replacing a binding cannot preserve what it lent, so it is
  ! not an operation. INTENT(INOUT): an INTENT(OUT) dummy would be
  ! finalized before this body ran.
  !===================================================================!

  subroutine refuse_assignment(lhs, rhs)

    class(relational_binding), intent(inout) :: lhs
    type(relational_binding) , intent(in)    :: rhs

    error stop 'graph_relational_view: a relational_binding is not assignable'

  end subroutine refuse_assignment

  !===================================================================!
  ! Release. Individually allocated objects are individually freed.
  !===================================================================!

  subroutine release_binding(this)

    type(relational_binding), intent(inout) :: this

    integer :: k

    if (allocated(this % sets)) then
       do k = 1, size(this % sets)
          if (associated(this % sets(k) % object)) then
             deallocate(this % sets(k) % object)
          end if
       end do
    end if

    if (allocated(this % relations)) then
       do k = 1, size(this % relations)
          if (associated(this % relations(k) % object)) then
             deallocate(this % relations(k) % object)
          end if
       end do
    end if

  end subroutine release_binding

  !===================================================================!
  ! The view. Counting and indexing are the sequence view's; this
  ! module only names the two branches and resolves the binding.
  !===================================================================!

  integer function num_member_sets(g) result(n)

    type(graph), intent(in) :: g

    n = sequence_size(g % branch(1))

  end function num_member_sets

  integer function num_relations(g) result(n)

    type(graph), intent(in) :: g

    n = sequence_size(g % branch(2))

  end function num_relations

  function member_set_at(g, b, k) result(s)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    integer                 , intent(in) :: k
    class(member_set), pointer           :: s

    type(graph), pointer :: element

    element => sequence_element(g % branch(1), k)
    s => b % set_for(element)

  end function member_set_at

  function relation_at(g, b, k) result(r)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    integer                 , intent(in) :: k
    class(relation), pointer             :: r

    type(graph), pointer :: element

    element => sequence_element(g % branch(2), k)
    r => b % relation_for(element)

  end function relation_at

  !===================================================================!
  ! Does this graph hold that member set. One scan of the binding to
  ! find the element that denotes it, then one traversal of the
  ! sequence: O(m + n), never O(m*n).
  !===================================================================!

  logical function holds_set(g, b, s) result(held)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(member_set)       , intent(in) :: s

    integer :: k

    held = .false.
    if (.not. allocated(b % sets)) return

    do k = 1, size(b % sets)
       if (b % sets(k) % object % same_as(s)) then
          held = sequence_contains(g % branch(1), b % sets(k) % element)
          return
       end if
    end do

  end function holds_set

  !===================================================================!
  ! The validity law. Every domain of every relation is a member set
  ! of this graph.
  !===================================================================!

  logical function relational_valid(g, b) result(ok)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b

    class(relation), pointer       :: r
    class(member_set), allocatable :: d
    integer                        :: k, j

    ok = .true.
    do k = 1, num_relations(g)
       r => relation_at(g, b, k)
       do j = 1, r % arity()
          d = r % domain(j)
          if (.not. holds_set(g, b, d)) then
             ok = .false.
             return
          end if
       end do
    end do

  end function relational_valid

end module graph_relational_view
