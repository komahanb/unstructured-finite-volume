!=====================================================================!
! LABEL MAP
!
! What a set is CALLED. Keyed on graph identity, stored outside the
! graph:
!
!     set graph  ->  character label
!
! The third association of the family, and the smallest. A graph
! answers WHICH object; a representation answers HOW ITS MEMBERS ARE
! STORED; this map answers WHAT IT IS CALLED, and nothing here is ever
! consulted to decide the first question.
!
!                    A LABEL IS NOT AN IDENTITY
!
! Two graphs may carry the same label and remain two sets, exactly as
! two graphs with equal extensions remain two sets. So this map is
! searched BY IDENTITY and never by name: there is no lookup from a
! label back to a graph, deliberately, because such a lookup would
! have to answer for a string that names two objects. Naming is not
! addressing.
!
!     same label, different graph   ->  still different
!     copied token, same graph      ->  same label
!
!                    WHY THIS EXISTS AT ALL
!
! The domain capability audit found one question the re-rooted domain
! law could not answer. type(graph) deliberately carries no name - the
! kernel type is branch(2) and a private token - yet production code
! names a derived subset after the domain it came from:
!
!     sg = subset_set(dom % name(), global_carrier, kept(1:n))
!
! That is metadata, and metadata belongs beside the mathematics rather
! than inside it. Putting a label in fractal_graph would teach the core
! a word it does not need; putting one in set_representation would make
! two descriptions of one set disagree about its name. So it lives
! here, orthogonal to both, and a consumer that does not name anything
! never carries this map at all.
!
! This is a LABEL map and not an attribute system. One noun, one role.
! When a second kind of metadata earns its place it will be measured
! and named on its own, not smuggled in behind a generalization.
!
!                        THE UNNAMED ANSWER
!
! An unbound graph is not an error - it is a set nobody named. label_of
! answers the empty string, matching the metadata convention the
! carriers already keep (name() answers '' when nobody said). Binding
! is the strict half: a set is named ONCE, and an unsigned token is
! refused, because both would leave two answers to one question.
!
!                          THE STORAGE LAW
!
! Rows key on type(token), copied at bind. This map borrows no graph
! object merely to recognize it, so it may outlive every variable that
! built it, and bind demands no TARGET.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_label_map

  use fractal_graph , only : graph
  use graph_identity, only : token

  implicit none

  private
  public :: label_map

  !===================================================================!
  ! One row: WHICH set - by value - and what it is called.
  !===================================================================!

  type :: label_row
     type(token)                   :: identity
     character(len=:), allocatable :: label
  end type label_row

  type :: label_map

     type(label_row), allocatable, private :: rows(:)

   contains

     procedure :: bind    => bind_label
     procedure :: labelled
     procedure :: label_of

  end type label_map

contains

  !===================================================================!
  ! Name a set. A set is named once: a second binding would leave two
  ! answers to one question, and the later one would win silently.
  !===================================================================!

  subroutine bind_label(this, element, label)

    class(label_map), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: label

    type(label_row), allocatable :: grown(:)
    type(token)                  :: key
    integer                      :: n

    ! An undeclared token does not match itself.
    key = element % id()
    if (.not. key % matches(key)) then
       error stop 'graph_label_map: a label map is keyed on assigned identity'
    end if

    if (row_at(this, key) /= 0) then
       error stop 'graph_label_map: a set is named once'
    end if

    if (.not. allocated(this % rows)) allocate(this % rows(0))

    n = size(this % rows)
    allocate(grown(n + 1))
    grown(1:n) = this % rows
    grown(n + 1) % identity = key
    grown(n + 1) % label    = label
    call move_alloc(grown, this % rows)

  end subroutine bind_label

  !===================================================================!
  ! Where the row for an identity is, or zero. The one comparison, and
  ! it reads nothing outside this map.
  !===================================================================!

  pure integer function row_at(this, key) result(at)

    class(label_map), intent(in) :: this
    type(token)     , intent(in) :: key

    integer :: k

    at = 0
    if (.not. allocated(this % rows)) return

    do k = 1, size(this % rows)
       if (this % rows(k) % identity % matches(key)) then
          at = k
          return
       end if
    end do

  end function row_at

  pure logical function labelled(this, element)

    class(label_map), intent(in) :: this
    type(graph)     , intent(in) :: element

    type(token) :: key

    key = element % id()
    labelled = row_at(this, key) /= 0

  end function labelled

  !===================================================================!
  ! What the set is called, as a value; '' when nobody named it. The
  ! unnamed answer is not a refusal, because being unnamed is not a
  ! contradiction - metadata is allowed to be absent.
  !===================================================================!

  pure function label_of(this, element) result(label)

    class(label_map), intent(in)  :: this
    type(graph)     , intent(in)  :: element
    character(len=:), allocatable :: label

    type(token) :: key
    integer     :: at

    key = element % id()
    at  = row_at(this, key)

    if (at == 0) then
       label = ''
    else
       label = this % rows(at) % label
    end if

  end function label_of

end module graph_label_map
