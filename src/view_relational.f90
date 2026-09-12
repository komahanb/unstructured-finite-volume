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
! and the kernel defines none of these terms.
!
! THE BINDING. A branch references a GRAPH, never an arbitrary object.
! The member sets and relations of this repository are not yet graphs,
! so each element graph denotes one of them and a binding maps the
! element to the legacy object it denotes. The binding OWNS its
! objects, because a non-owning pointer passed to a caller can outlive
! the object it references.
!
! THE STORAGE LAW.
!
!     a binding is a counted reference to its objects; bind_* preserves
!     every outstanding object pointer while any owner of the binding
!     remains.
!
! That law has a cost. A row stores a POINTER to an individually
! allocated object, never the object itself: when the row array grows,
! the rows are copied and the objects do not move. Storing the object
! in the row - as an allocatable component - was measured and rejected,
! because growth relocates the array and every non-owning pointer then
! reads freed storage: incorrect values first, a fatal error next. See
! test/graph-relational/lifetime.f90, which checks this law on every run.
!
! ASSIGNMENT BINDS ONE MORE OWNER of the same objects
! (util_counted_storage); the last owner's finalization deallocates
! them. No Fortran mechanism prohibits assignment at compile time;
! four were compiled and measured in
! test/graph-relational/fortran-assignment, and a deep copy was
! measured and rejected (test/graph-relational/lifetime.f90).
! Extension and replacement are different operations: bind_* extends
! a binding with one owner and preserves every pointer it has
! returned; a binding with more than one owner is immutable and
! bind_* refuses. A copy made without defined assignment (allocate
! with source=, a structure constructor, polymorphic assignment) is
! not an owner: the first of the two to be finalized releases the
! binding and the other is refused at its next access.
!
! Because a row stores a pointer rather than the object, the pointer
! this module returns does not point into the binding. The binding
! therefore needs no TARGET attribute at any call site.
!
! The binding is storage keyed on identity, not ontology: its rows are
! an identity-row table of copied element tokens (map_token_rows), so
! it references no element graph and bind_* needs no TARGET either.
! The wrappers that were stored inside the retired container belong
! here: they existed only because a Fortran array has one dynamic
! type, and that is a storage fact.
!
! Sequence behaviour is delegated to view_sequence. Nothing here
! reads a cell.
!
! THREE FAILURES, STRUCTURALLY APART.
!
!     malformed sequence   rejected, by view_sequence
!     unstorable object    rejected, by bind_set / bind_relation
!     relationally invalid returned .false. by relational_valid
!
! A view over an existing graph reports invalidity; it does not reject
! a graph it did not construct. The retired constructor rejected at
! construction because it was a constructor; this is not one. What the
! binding rejects is not the view but the storage: an object with no
! identity cannot be compared, and a borrowing view cannot be owned,
! because copying it into owned storage copies a reference to a base
! the binding does not keep allocated.
!
! THE VALIDITY LAW:
!
!     G is relationally valid iff
!       (i)   no member set occurs twice in the member-set sequence,
!       (ii)  no relation occurs twice in the relation sequence, and
!       (iii) every domain of every relation is a member set of G.
!
! S and P are SETS; the branches represent them as sequences, and a
! sequence may repeat what a set cannot. (i) and (ii) are that
! difference, reported rather than rejected, because repetition is a
! property of a constructed graph.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_relational

  use graph_fractal    , only : graph, branch
  use token_identity   , only : token
  use util_counted_storage, only : counted_storage, counted_reference
  use map_token_rows   , only : identity_rows
  use relation_finitary, only : relation
  use view_sequence    , only : sequence_num_elements, sequence_element, &
       & sequence_empty, sequence_first, sequence_rest

  implicit none

  private
  public :: relational_binding
  public :: num_member_sets, member_set_at
  public :: num_relations, relation_at
  public :: has_set, relational_valid

  !===================================================================!
  ! Owned storage. The elements are the keys of two identity-row
  ! tables; the objects run parallel to the rows, each allocated
  ! separately and referenced by pointer. The cell is jointly owned
  ! by every binding value assigned from it.
  !===================================================================!

  type :: bound_set
     type(graph)    , pointer :: object => null()
  end type bound_set

  type :: bound_relation
     class(relation), pointer :: object => null()
  end type bound_relation

  type, extends(counted_storage) :: bound_objects

     type(identity_rows)              :: set_rows, relation_rows
     type(bound_set)     , allocatable :: sets(:)
     type(bound_relation), allocatable :: relations(:)

   contains

     procedure :: clear => release_objects

  end type bound_objects

  type :: relational_binding

     type(counted_reference), private :: reference

   contains

     procedure :: bind_set
     procedure :: bind_relation
     procedure :: set_for
     procedure :: relation_for
     procedure :: num_owners

     procedure, private :: cell
     procedure, private :: extended

  end type relational_binding

contains

  !===================================================================!
  ! Binding. The object is copied into owned storage; the element's
  ! identity is the key, so an element is declared and bound once.
  !
  ! What may be stored: an object with an assigned identity, because
  ! the view compares objects and nothing else; and, for a relation,
  ! one that is materialized, because a borrowing view copied into
  ! owned storage stores a reference to a base the binding does not
  ! keep allocated. A view is defined over a bound relation, never
  ! stored inside one.
  !===================================================================!

  subroutine bind_set(this, element, object)

    class(relational_binding), intent(inout) :: this
    type(graph)              , intent(in)    :: element
    type(graph)              , intent(in)    :: object

    type(bound_objects), pointer :: objects
    integer :: at

    ! An undeclared token does not match itself.
    if (.not. object % same_as(object)) then
       error stop 'view_relational: a binding stores identified objects'
    end if

    objects => this % extended()
    at = objects % set_rows % append(element % id(), &
         & 'view_relational: a binding is keyed on assigned identity', &
         & 'view_relational: an element is bound once')

    if (.not. allocated(objects % sets)) allocate(objects % sets(0))
    objects % sets = [objects % sets, bound_set()]
    allocate(objects % sets(at) % object, source=object)

  end subroutine bind_set

  subroutine bind_relation(this, element, object)

    class(relational_binding), intent(inout) :: this
    type(graph)              , intent(in)    :: element
    class(relation)          , intent(in)    :: object

    type(bound_objects), pointer :: objects
    integer :: at

    if (.not. object % same_as(object)) then
       error stop 'view_relational: a binding stores identified objects'
    end if

    if (.not. object % materialized()) then
       error stop 'view_relational: a binding owns whole relations; a view cannot be bound'
    end if

    objects => this % extended()
    at = objects % relation_rows % append(element % id(), &
         & 'view_relational: a binding is keyed on assigned identity', &
         & 'view_relational: an element is bound once')

    if (.not. allocated(objects % relations)) allocate(objects % relations(0))
    objects % relations = [objects % relations, bound_relation()]
    allocate(objects % relations(at) % object, source=object)

  end subroutine bind_relation

  !===================================================================!
  ! The cell this binding is bound to. A binding whose objects were
  ! released stops the program: the last owner deallocated them, or
  ! this value is a bitwise copy whose twin released the binding.
  ! A binding never bound has no cell.
  !===================================================================!

  function cell(this) result(objects)

    class(relational_binding), intent(in) :: this
    type(bound_objects), pointer          :: objects

    class(counted_storage), pointer :: storage

    objects => null()
    if (this % reference % released()) then
       error stop 'view_relational: this binding''s objects have been released'
    end if
    storage => this % reference % storage()
    if (.not. associated(storage)) return
    select type (storage)
    type is (bound_objects)
       objects => storage
    end select

  end function cell

  !===================================================================!
  ! The cell to extend: acquired by the first bind_*; a binding with
  ! more than one owner is immutable and stops the program.
  !===================================================================!

  function extended(this) result(objects)

    class(relational_binding), intent(inout) :: this
    type(bound_objects), pointer             :: objects

    type(bound_objects) :: template

    if (this % reference % num_owners() > 1) then
       error stop 'view_relational: a binding is extended by its sole owner'
    end if
    if (.not. this % reference % live()) call this % reference % acquire(template)
    objects => this % cell()

  end function extended

  ! The number of binding values bound to these objects; zero for a
  ! binding without any.
  pure integer function num_owners(this)

    class(relational_binding), intent(in) :: this

    num_owners = this % reference % num_owners()

  end function num_owners

  !===================================================================!
  ! Lookup by element identity, returning a reference into owned
  ! storage.
  !===================================================================!

  function set_for(this, element) result(s)

    class(relational_binding), intent(in) :: this
    type(graph)              , intent(in) :: element
    type(graph), pointer                  :: s

    type(bound_objects), pointer :: objects

    objects => this % cell()
    if (.not. associated(objects)) then
       error stop 'view_relational: no member set is bound to that element'
    end if
    s => objects % sets(objects % set_rows % row(element % id(), &
         & 'view_relational: no member set is bound to that element')) % object

  end function set_for

  function relation_for(this, element) result(r)

    class(relational_binding), intent(in) :: this
    type(graph)              , intent(in) :: element
    class(relation), pointer              :: r

    type(bound_objects), pointer :: objects

    objects => this % cell()
    if (.not. associated(objects)) then
       error stop 'view_relational: no relation is bound to that element'
    end if
    r => objects % relations(objects % relation_rows % row(element % id(), &
         & 'view_relational: no relation is bound to that element')) % object

  end function relation_for

  !===================================================================!
  ! Rejection. Replacing a binding cannot preserve the pointers it has
  ! returned, so replacement is not an operation. INTENT(INOUT): an
  ! INTENT(OUT) dummy would be finalized before this body executed.
  !===================================================================!

  !===================================================================!
  ! Release every object of a cell, when its last owner is finalized.
  !===================================================================!

  subroutine release_objects(this)

    class(bound_objects), intent(inout) :: this

    type(identity_rows) :: no_rows
    integer :: k

    do k = 1, this % set_rows % num_rows()
       if (associated(this % sets(k) % object)) deallocate(this % sets(k) % object)
    end do

    do k = 1, this % relation_rows % num_rows()
       if (associated(this % relations(k) % object)) deallocate(this % relations(k) % object)
    end do

    if (allocated(this % sets)) deallocate(this % sets)
    if (allocated(this % relations)) deallocate(this % relations)
    ! the cell is recycled by the next acquisition: its key tables are
    ! emptied with its objects
    this % set_rows = no_rows
    this % relation_rows = no_rows

  end subroutine release_objects


  integer function num_member_sets(g) result(n)

    type(graph), intent(in) :: g

    n = sequence_num_elements(g % branch(1))

  end function num_member_sets

  integer function num_relations(g) result(n)

    type(graph), intent(in) :: g

    n = sequence_num_elements(g % branch(2))

  end function num_relations

  function member_set_at(g, b, k) result(s)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    integer                 , intent(in) :: k
    type(graph), pointer             :: s

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
  ! Whether this graph contains that member set. One scan of the binding to
  ! find the element that denotes it, then one traversal of the
  ! sequence: O(m + n), never O(m*n).
  !===================================================================!

  logical function has_set(g, b, s) result(stored)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    type(graph)             , intent(in) :: s

    type(bound_objects), pointer :: objects
    integer :: k

    stored = .false.

    objects => b % cell()
    if (.not. associated(objects)) return

    do k = 1, objects % set_rows % num_rows()
       if (objects % sets(k) % object % same_as(s)) then
          stored = sequence_has_key(g % branch(1), objects % set_rows % key(k))
          return
       end if
    end do

  end function has_set

  !===================================================================!
  ! Does the sequence contain the element with this identity: the
  ! element is stored as a key, so membership is read by token. A
  ! chain of cells that reaches UNKNOWN is refused by sequence_first.
  !===================================================================!

  recursive logical function sequence_has_key(b, key) result(found)

    type(branch), intent(in) :: b
    type(token) , intent(in) :: key

    type(graph), pointer :: element

    found = .false.
    if (sequence_empty(b)) return

    element => sequence_first(b)
    found   =  key % matches(element % id())
    if (.not. found) found = sequence_has_key(sequence_rest(b), key)

  end function sequence_has_key

  !===================================================================!
  ! The validity law. S and P are sets, and every domain of every
  ! relation is a member set of this graph.
  !===================================================================!

  logical function relational_valid(g, b) result(valid)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b

    type(graph)  , pointer     :: s, s_earlier
    class(relation)  , pointer     :: r, r_earlier
    type(graph)                :: d
    integer                        :: k, j

    valid = .false.

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
          if (.not. has_set(g, b, d)) return
       end do
    end do

    valid = .true.

  end function relational_valid

end module view_relational
