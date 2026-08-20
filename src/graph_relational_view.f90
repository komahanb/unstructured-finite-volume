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
! the wrappers that used to sit inside the retired container belong:
! they existed only because a Fortran array carries one dynamic type,
! and that was always a storage fact.
!
! Sequence behaviour is delegated to graph_sequence_view. Nothing here
! traverses a spine.
!
! THREE FAILURES, STRUCTURALLY APART.
!
!     malformed sequence   refused, by graph_sequence_view
!     unstorable object    refused, by bind_set / bind_relation
!     relationally invalid answered .false. by relational_valid
!
! A view over an existing graph reports invalidity; it does not refuse
! a graph it did not construct. The retired constructor refused at
! construction because it was a constructor; this is not one. What the
! binding refuses is not the view but the storage: an object with no
! identity cannot be compared, and a borrowing view cannot be owned,
! because copying it into owned storage copies a reference to a base
! the binding does not keep alive.
!
! THE VALIDITY LAW:
!
!     G is relationally valid iff
!       (i)   no member set occurs twice in the member-set sequence,
!       (ii)  no relation occurs twice in the relation sequence, and
!       (iii) every domain of every relation is a member set of G.
!
! S and P are SETS; the branches represent them as sequences, and a
! sequence may repeat what a set cannot. (i) and (ii) are that gap,
! answered rather than refused, because repetition is a property of a
! graph already built.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_relational_view

  use graph_fractal      , only : graph
  use graph_fractal      , only : set_graph => graph
  use relation_finitary     , only : relation
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
     type(set_graph)  , pointer :: object  => null()
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
  !
  ! What may be stored: an object with an assigned identity, because
  ! the view compares objects and nothing else; and, for a relation,
  ! one that is materialized, because a borrowing view copied into
  ! owned storage carries a reference to a base the binding does not
  ! keep alive. A view rides above a bound relation, never inside one.
  !===================================================================!

  subroutine bind_set(this, element, object)

    class(relational_binding), intent(inout)        :: this
    type(graph)              , intent(in) , target  :: element
    type(set_graph)          , intent(in)           :: object

    type(bound_set), allocatable :: grown(:)
    integer                      :: n

    ! An undeclared token does not match itself.
    if (.not. object % same_as(object)) then
       error stop 'graph_relational_view: a binding stores identified objects'
    end if

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

    if (.not. object % same_as(object)) then
       error stop 'graph_relational_view: a binding stores identified objects'
    end if

    if (.not. object % materialized()) then
       error stop 'graph_relational_view: a binding owns whole relations; a view cannot be bound'
    end if

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
    type(set_graph), pointer              :: s

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
    type(set_graph), pointer             :: s

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
    type(set_graph)         , intent(in) :: s

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
  ! The validity law. S and P are sets, and every domain of every
  ! relation is a member set of this graph.
  !===================================================================!

  logical function relational_valid(g, b) result(ok)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b

    type(set_graph)  , pointer     :: s, s_earlier
    class(relation)  , pointer     :: r, r_earlier
    type(set_graph)                :: d
    integer                        :: k, j

    ok = .false.

    do k = 1, num_member_sets(g)                     ! (i) S is a set
       s => member_set_at(g, b, k)
       do j = 1, k - 1
          s_earlier => member_set_at(g, b, j)
          if (s % same_as(s_earlier)) return
       end do
    end do

    do k = 1, num_relations(g)                       ! (ii) P is a set
       r => relation_at(g, b, k)
       do j = 1, k - 1
          r_earlier => relation_at(g, b, j)
          if (r % same_as(r_earlier)) return
       end do
    end do

    do k = 1, num_relations(g)                       ! (iii) closure
       r => relation_at(g, b, k)
       do j = 1, r % arity()
          d = r % domain(j)
          if (.not. holds_set(g, b, d)) return
       end do
    end do

    ok = .true.

  end function relational_valid

end module graph_relational_view
