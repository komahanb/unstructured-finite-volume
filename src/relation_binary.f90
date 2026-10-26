!=====================================================================!
! Module that states level 1 of the new tower, the binary
! specialization.
!
! Arity two has its own contract (AGENTS.md 5.3) because arity two has
! its own canonical queries. For
!
!      P  <=  A x B
!
! the specialization adds what no general arity can provide:
!
!      source, target      the two ends of the signature, named
!      image(a)            the members of B that a relates to
!      preimage(b)         the members of A that relate to b
!
! and the transpose is the canonical slot permutation, produced here
! as a view.
!
! Members in, members out. Every query uses member values, never
! storage rows. A sparse carrier may store { 10 20 30 }; image(20) is
! called with 20 and returns members of the other domain. The map
! between a member and its storage row is the carrier's own
! local_index, the inverse enumeration the carriers guarantee, so the
! indexed lookup below never assumes a domain is 1..n. The image of a
! non-member is the empty set: relating to nothing is a result, not an
! error.
!
!      carrier members    10    20    30
!      rows                1     2     3        local_index(20) = 2
!
! Two tiers of traversal. The deferred primitives are the views:
! image_view and preimage_view return a fibre with a private reference
! to the stored index, no allocation, no copy, the hot-loop path
! (AGENTS.md 33). The allocating image and preimage are defined above
! them as conveniences, written once for the whole family as copies of
! the views. Scalar reads and read callbacks cannot modify the stored
! index. A caller storing a view stores a non-owning reference: the
! relation must have TARGET and must not be reassigned or deallocated
! while the view is used. An owned copy is obtained through values
! instead.
!
! The CSR representation. csr_relation stores both directions, built
! once at construction:
!
!      xfwd, tgt      row a-local  ->  its B members     image
!      xbwd, src      row b-local  ->  its A members     preimage
!
!      P = { (1,1) (1,2) (2,2) (3,1) },  A = { 1 2 3 },  B = { 1 2 }
!
!      forward     row      1          2       3
!                  xfwd     1          3       4       5
!                  tgt      1   2   |  2   |   1             image(1) = { 1 2 }
!
!      backward    row      1             2
!                  xbwd     1             3            5
!                  src      1   3     |   1   2              preimage(2) = { 1 2 }
!
! so each fibre is one row slice, has([a,b]) is one row scan, and after
! member-to-row lookup, duplicate collapse and both index builds are
! linear in members plus tuples. Set semantics hold here exactly as in
! the stored relation: a tuple passed in twice is in the relation once,
! first appearance retaining its position.
!
! Complexity, parameterized exactly. Every fibre first reads the
! member's position from the carrier, so the total cost is
!
!      T_image(a)  =  T_local_index(a) + O(deg a)
!
! and the slice alone is O(deg). A counted carrier computes
! local_index in O(1), so the bound reduces to O(deg) there, the mesh
! path's case. A listed carrier uses a sorted inverse index, so its
! lookup adds O(log n) for arbitrary sparse member values.
!
! The view, and its lifetime requirement. transpose_of(r) returns a
! view: O(1) to construct, no topology copied, image and preimage
! swapped, the signature read in reverse. A view is non-owning: the
! view stores its base by pointer, and the base must outlive the view;
! the caller's base must have the target attribute. That is the whole
! cost of an O(1) transpose in a language of value semantics, and it
! is stated explicitly.
!
!      P     <=  A x B          image_view(a)  = row a of xfwd
!      P^T   <=  B x A          image_view(b)  = P's preimage_view(b)
!                               a pointer to P; nothing copied
!
! The ownership policy, declared for the levels above. When the graph
! is constructed and contains relations (AGENTS.md 14), the law is:
!
!      the graph owns stable relations;
!      views and fibre references may reference them.
!
! A graph accessor must therefore return its relations by reference to
! owned, stable storage, never as temporary copies that a view or
! fibre could reference after deallocation. The owner of the base
! decides its lifetime; every borrower's lifetime is strictly
! contained in it.
!
! Identity is not equality. same_as compares assigned identity: a view
! has its own token, so a view is never same_as its base, and the
! involution
!
!      (P^T)^T = P
!
! is a statement about extension, the same tuples over the same
! domains, not about tokens. Test it by comparing tuples and comparing
! domains slot against slot; only a deliberate canonicalization could
! guarantee it by identity, and none is guaranteed here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_binary

  use graph_fractal         , only : graph
  use relation_finitary     , only : relation
  use map_set               , only : set_map
  use map_set_representation, only : set_representation

  implicit none

  private

  public :: binary_relation, csr_relation, transposed_relation, transpose_of
  public :: group_by_key
  public :: ragged
  public :: transpose_padded
  public :: integer_fibre

  !-------------------------------------------------------------------!
  ! Type that models a finite ordered fibre in compiled integer
  ! coordinates. Its private reference grants read access only: neither
  ! a member nor a callback argument is a writable alias. Construction
  ! and assignment copy only the reference. The target array must
  ! outlive every use of the fibre.
  !-------------------------------------------------------------------!

  type :: integer_fibre
     private
     integer, pointer :: entries(:) => null()
   contains
     procedure, non_overridable :: num_members => fibre_num_members
     procedure, non_overridable :: member => fibre_member
     procedure, non_overridable :: values => fibre_values
     procedure, private, non_overridable :: read_members => fibre_read
     procedure, private, non_overridable :: read_context => fibre_read_context
     generic :: read => read_members, read_context
  end type integer_fibre

  interface integer_fibre
     module procedure create_fibre
  end interface integer_fibre

  abstract interface

     !----------------------------------------------------------------!
     ! Interface for a reader of a fibre's members
     !----------------------------------------------------------------!

     subroutine integer_reader_interface(members)
       integer, intent(in) :: members(:)
     end subroutine integer_reader_interface

     !----------------------------------------------------------------!
     ! Interface for a reader of a fibre's members with state of the
     ! caller's own
     !----------------------------------------------------------------!

     subroutine integer_context_reader_interface(members, context)
       integer , intent(in)    :: members(:)
       class(*), intent(inout) :: context
     end subroutine integer_context_reader_interface

  end interface

  !-------------------------------------------------------------------!
  ! Type that models the abstract binary relation: the general
  ! contract, plus the queries only arity two defines. Every descendant
  ! declares a signature of two domains, so arity is the root's.
  !-------------------------------------------------------------------!

  type, abstract, extends(relation) :: binary_relation

   contains

     !----------------------------------------------------------------!
     ! The deferred primitives: fibres as non-owning references, no
     ! allocation
     !----------------------------------------------------------------!

     procedure(binary_fibre_view_interface), deferred :: image_view
     procedure(binary_fibre_view_interface), deferred :: preimage_view

     !----------------------------------------------------------------!
     ! The conveniences, written once for the family: copies of the
     ! views, for callers that require an owned copy rather than a
     ! non-owning reference
     !----------------------------------------------------------------!

     procedure :: image
     procedure :: preimage

     procedure :: source
     procedure :: target

  end type binary_relation

  abstract interface

     !----------------------------------------------------------------!
     ! Interface for the fibre of a member, as a view
     !----------------------------------------------------------------!

     function binary_fibre_view_interface(this, member) result(fibre)
       import binary_relation, integer_fibre
       class(binary_relation), target, intent(in) :: this
       integer               , intent(in)         :: member
       type(integer_fibre)                        :: fibre
     end function binary_fibre_view_interface

  end interface

  !-------------------------------------------------------------------!
  ! Type that models the CSR representation: both directions
  ! materialized once, every query an O(degree) slice.
  !
  ! Two queries, two stores, and they are disjoint:
  !
  !     signature      which domains        semantic, identities
  !     coordinates    which row a member   compiled, numbering only
  !
  ! The signature, stored by the root relation, returns domain(k) and
  ! nothing else; the coordinates return local_index and nothing else.
  ! A coordinate representation stores no identity: it cannot state
  ! which set it numbers, and does not need to, because the signature
  ! beside it already does.
  !
  ! The coordinates are stored by value, copied out of the caller's set
  ! map at construction. That is deliberate and differs from the
  ! field's case: many fields share one domain, so copying an extent
  ! per field was measured to cost too much; a CSR relation's row
  ! numbering is part of its own compiled execution contract, and the
  ! hot path may not search for it. So the coordinates are stored here,
  ! and image/preimage/has read them directly: no map row scan, no
  ! graph traversal, no label lookup.
  !-------------------------------------------------------------------!

  type, extends(binary_relation) :: csr_relation

     ! Compiled: how a member value becomes a row, each direction. The
     ! signature, which two domains, is the root's.
     class(set_representation), allocatable, private :: source_coords
     class(set_representation), allocatable, private :: target_coords

     integer             , private :: nnz = 0
     integer, allocatable, private :: xfwd(:), tgt(:)
     integer, allocatable, private :: xbwd(:), src(:)

   contains

     procedure :: has           => csr_has
     procedure :: num_tuples    => csr_num_tuples
     procedure :: tuples        => csr_tuples
     procedure :: image_view    => csr_image_view
     procedure :: preimage_view => csr_preimage_view
     procedure :: materialized  => csr_materialized

  end type csr_relation

  interface csr_relation
     module procedure create_csr
  end interface csr_relation

  !-------------------------------------------------------------------!
  ! Type that models the transpose view: a borrower. The view stores
  ! its base by pointer and evaluates every tuple query through the
  ! base, ends swapped; its signature is the base's two domains
  ! swapped, copied at construction. The base must outlive the view.
  !-------------------------------------------------------------------!

  type, extends(binary_relation) :: transposed_relation

     class(binary_relation), pointer, private :: base => null()

   contains

     procedure :: has           => view_has
     procedure :: num_tuples    => view_num_tuples
     procedure :: tuples        => view_tuples
     procedure :: image_view    => view_image_view
     procedure :: preimage_view => view_preimage_view

     ! No materialized binding, deliberately: the root's default
     ! already returns false, and a borrower (copying it copies a
     ! pointer to a base it does not keep allocated) is exactly what
     ! the default rejects. Views are defined over graph-owned
     ! relations, never stored inside them.

  end type transposed_relation

  !-------------------------------------------------------------------!
  ! Type that models a list of lists, compressed: the entries of list
  ! k are entries(first(k) : first(k+1) - 1). This is the shape
  ! group_by_key produces and every compressed traversal reads. The
  ! padded shape, one fixed width with a count per list, is the same
  ! relation stored with unused capacity, and the two are converted
  ! here and nowhere else.
  !
  !      padded, width 3             ragged
  !      x(:, 1)   4  7  0   n = 2   first     1  3  4  7
  !      x(:, 2)   5  0  0   n = 1   entries   4  7 | 5 | 6  8  9
  !      x(:, 3)   6  8  9   n = 3
  !-------------------------------------------------------------------!

  type :: ragged
     integer, allocatable :: first(:)
     integer, allocatable :: entries(:)
   contains
     procedure :: num_lists => ragged_num_lists
     procedure :: length    => ragged_length
     procedure :: list      => ragged_list
     procedure :: padded    => ragged_padded
  end type ragged

  interface ragged
     module procedure ragged_of_padded
     module procedure ragged_of_lists
  end interface ragged

contains

  !===================================================================!
  ! Creates a fibre referencing the members given
  !===================================================================!

  function create_fibre(members) result(fibre)

    integer, target, intent(in) :: members(:)
    type(integer_fibre)         :: fibre

    fibre % entries => members

  end function create_fibre

  !===================================================================!
  ! Returns the number of members of the fibre
  !===================================================================!

  pure integer function fibre_num_members(this) result(n)

    class(integer_fibre), intent(in) :: this

    n = 0
    if (associated(this % entries)) n = size(this % entries)

  end function fibre_num_members

  !===================================================================!
  ! Returns the member at a position of the fibre; a position outside
  ! it stops the program
  !===================================================================!

  pure integer function fibre_member(this, index) result(member)

    class(integer_fibre), intent(in) :: this
    integer             , intent(in) :: index

    if (.not. associated(this % entries)) then
       error stop 'relation_binary: fibre_member was called on a fibre with no entries associated'
    end if
    if (index < 1 .or. index > size(this % entries)) then
       error stop 'relation_binary: the index argument lies outside the fibre''s extent'
    end if
    member = this % entries(index)

  end function fibre_member

  !===================================================================!
  ! Returns an owned copy of the fibre's members
  !===================================================================!

  function fibre_values(this) result(members)

    class(integer_fibre), intent(in) :: this
    integer, allocatable             :: members(:)

    call copy_fibre(this, members)

  end function fibre_values

  !===================================================================!
  ! Copies the fibre's members, none for a fibre of no entries
  !===================================================================!

  subroutine copy_fibre(this, members)

    class(integer_fibre), intent(in)  :: this
    integer, allocatable, intent(out) :: members(:)

    if (associated(this % entries)) then
       members = this % entries
    else
       allocate(members(0))
    end if

  end subroutine copy_fibre

  !===================================================================!
  ! Passes the fibre's members to a reader, without a copy
  !===================================================================!

  subroutine fibre_read(this, reader)

    class(integer_fibre), intent(in)    :: this
    procedure(integer_reader_interface) :: reader

    if (associated(this % entries)) then
       call reader(this % entries)
    else
       call reader([integer ::])
    end if

  end subroutine fibre_read

  !===================================================================!
  ! Passes the fibre's members to a reader with the caller's state.
  ! Explicit per-call state allows a module procedure to read a fibre
  ! without a compiler-generated closure or module-global mutable
  ! state.
  !===================================================================!

  subroutine fibre_read_context(this, reader, context)

    class(integer_fibre), intent(in)            :: this
    procedure(integer_context_reader_interface) :: reader
    class(*)            , intent(inout)         :: context

    if (associated(this % entries)) then
       call reader(this % entries, context)
    else
       call reader([integer ::], context)
    end if

  end subroutine fibre_read_context

  !===================================================================!
  ! Returns the source, the first end of the signature, as arity two
  ! names it
  !===================================================================!

  type(graph) function source(this) result(domain)

    class(binary_relation), intent(in) :: this

    domain = this % domain(1)

  end function source

  !===================================================================!
  ! Returns the target, the second end of the signature
  !===================================================================!

  type(graph) function target(this) result(domain)

    class(binary_relation), intent(in) :: this

    domain = this % domain(2)

  end function target

  !===================================================================!
  ! Returns the image of a member, an owned copy of what the view
  ! references. Written once, here, for every binary implementation
  ! present and future.
  !===================================================================!

  subroutine image(this, member, indices)

    class(binary_relation), target, intent(in)  :: this
    integer                       , intent(in)  :: member
    integer, allocatable          , intent(out) :: indices(:)
    type(integer_fibre)                         :: fibre

    fibre = this % image_view(member)
    call copy_fibre(fibre, indices)

  end subroutine image

  !===================================================================!
  ! Returns the preimage of a member, an owned copy of what the view
  ! references
  !===================================================================!

  subroutine preimage(this, member, indices)

    class(binary_relation), target, intent(in)  :: this
    integer                       , intent(in)  :: member
    integer, allocatable          , intent(out) :: indices(:)
    type(integer_fibre)                         :: fibre

    fibre = this % preimage_view(member)
    call copy_fibre(fibre, indices)

  end subroutine preimage

  !===================================================================!
  ! Creates a CSR relation: a name, the two carriers (any concretions,
  ! each validated by its own checks) and the tuple table, one column
  ! per tuple, members throughout. Input checks first, as at every
  ! constructor of the level; then the duplicate collapse and both
  ! index builds, linear apart from carrier lookup:
  !
  !      count rows        one pass with local_index
  !      place tuples      counting sort by source row
  !      collapse          one marker array over target rows
  !      backward build    the same, mirrored
  !===================================================================!

  type(csr_relation) function create_csr(name, source, target, table, sets) &
       & result(this)

    character(len=*), intent(in) :: name
    type(graph)     , intent(in) :: source
    type(graph)     , intent(in) :: target
    integer         , intent(in) :: table(:,:)
    type(set_map)   , intent(in) :: sets

    integer, allocatable :: aloc(:), bloc(:), order(:), marker(:)
    integer, allocatable :: source_members(:), target_rows(:)
    integer              :: na, nb, nt
    integer              :: j, p, q, row, col, num_retained
    character(len=250)   :: message

    call this % declare(name, [source, target])

    if (size(table, 1) /= 2) then
       write(message,'(a,i0)') 'relation_binary: each tuple must have exactly one part per &
            &position (source, target); size(table,1) = ', size(table, 1)
       error stop trim(message)
    end if

    !-----------------------------------------------------------------!
    ! Compilation. The map is read here and only here: the two extents
    ! are copied in, and from this line on the relation numbers its own
    ! rows. Nothing below, and nothing in image, preimage or has, reads
    ! the map.
    !-----------------------------------------------------------------!

    call sets % extent_of(source, this % source_coords)
    call sets % extent_of(target, this % target_coords)

    na = this % source_coords % num_members()
    nb = this % target_coords % num_members()
    nt = size(table, 2)

    ! Validate through the coordinates' own membership, and read every
    ! member's row through their own inverse enumeration
    allocate(aloc(nt), bloc(nt))
    do j = 1, nt
       if (.not. this % source_coords % has(table(1, j))) then
          write(message,'(a,i0,a,i0)') 'relation_binary: tuple ', j, ' names a source member the &
               &source domain does not contain; table(1,j) = ', table(1, j)
          error stop trim(message)
       end if
       if (.not. this % target_coords % has(table(2, j))) then
          write(message,'(a,i0,a,i0)') 'relation_binary: tuple ', j, ' names a target member the &
               &target domain does not contain; table(2,j) = ', table(2, j)
          error stop trim(message)
       end if
       aloc(j) = this % source_coords % local_index(table(1, j))
       bloc(j) = this % target_coords % local_index(table(2, j))
    end do

    ! Forward: group the tuples by source row, duplicates included, then
    ! collapse each row with one pass over the marker array: marker(col)
    ! is the last row that retained col
    block
      integer, allocatable :: ptr(:), identity(:)
      allocate(identity(nt))
      identity = [(j, j = 1, nt)]
      call group_by_key(na, aloc, identity, ptr, order)

      allocate(this % xfwd(na + 1))
      allocate(this % tgt(nt), marker(max(nb, 1)))
      allocate(source_members(nt), target_rows(nt))
      marker = 0
      num_retained  = 0
      do row = 1, na
         p = ptr(row)
         q = ptr(row + 1) - 1
         this % xfwd(row) = num_retained + 1
         do j = p, q
            col = bloc(order(j))
            if (marker(col) /= row) then
               marker(col)                  = row
               num_retained                 = num_retained + 1
               this % tgt(num_retained)     = table(2, order(j))
               source_members(num_retained) = table(1, order(j))
               target_rows(num_retained)    = col
            end if
         end do
      end do
      this % xfwd(na + 1) = num_retained + 1
      this % nnz          = num_retained
    end block

    ! Backward: the retained tuples grouped by target row
    call group_by_key(nb, target_rows(1:num_retained), source_members(1:num_retained), &
         & this % xbwd, this % src)

  end function create_csr

  !===================================================================!
  ! Returns the image of a member as a view: one local_index, one
  ! slice, zero allocation. Members in, members out; a non-member's
  ! fibre is the empty slice. Cost: T_local_index(member) + O(1) to
  ! construct, O(degree) to read.
  !===================================================================!

  function csr_image_view(this, member) result(fibre)

    class(csr_relation), target, intent(in) :: this
    integer                    , intent(in) :: member
    type(integer_fibre)                     :: fibre

    integer :: row

    row = this % source_coords % local_index(member)
    if (row == 0) then
       fibre % entries => this % tgt(1:0)
       return
    end if

    fibre % entries => this % tgt(this % xfwd(row) : this % xfwd(row + 1) - 1)

  end function csr_image_view

  !===================================================================!
  ! Returns the preimage of a member as a view, the same way through
  ! the backward index
  !===================================================================!

  function csr_preimage_view(this, member) result(fibre)

    class(csr_relation), target, intent(in) :: this
    integer                    , intent(in) :: member
    type(integer_fibre)                     :: fibre

    integer :: row

    row = this % target_coords % local_index(member)
    if (row == 0) then
       fibre % entries => this % src(1:0)
       return
    end if

    fibre % entries => this % src(this % xbwd(row) : this % xbwd(row + 1) - 1)

  end function csr_preimage_view

  !===================================================================!
  ! Returns whether a tuple is in the relation: one row, one scan,
  ! O(degree), as stated above
  !===================================================================!

  pure logical function csr_has(this, tuple)

    class(csr_relation), intent(in) :: this
    integer            , intent(in) :: tuple(:)

    integer :: row, j

    csr_has = .false.

    if (size(tuple) /= 2) return

    row = this % source_coords % local_index(tuple(1))
    if (row == 0) return

    do j = this % xfwd(row), this % xfwd(row + 1) - 1
       if (this % tgt(j) == tuple(2)) then
          csr_has = .true.
          return
       end if
    end do

  end function csr_has

  !===================================================================!
  ! Returns the number of tuples
  !===================================================================!

  pure integer function csr_num_tuples(this)

    class(csr_relation), intent(in) :: this

    csr_num_tuples = this % nnz

  end function csr_num_tuples

  !===================================================================!
  ! Returns the tuple table, rebuilt row by row: non-hot generic
  ! access
  !===================================================================!

  pure subroutine csr_tuples(this, table)

    class(csr_relation) , intent(in)  :: this
    integer, allocatable, intent(out) :: table(:,:)

    integer :: row, j, a

    allocate(table(2, this % nnz))
    do row = 1, size(this % xfwd) - 1
       a = this % source_coords % member(row)
       do j = this % xfwd(row), this % xfwd(row + 1) - 1
          table(1, j) = a
          table(2, j) = this % tgt(j)
       end do
    end do

  end subroutine csr_tuples

  !===================================================================!
  ! Returns the transpose view: O(1), nothing copied, a new identity of
  ! its own. The base must have the target attribute, and must outlive
  ! the view: a view references, it never owns.
  !===================================================================!

  function transpose_of(base) result(view)

    class(binary_relation), target, intent(in) :: base
    type(transposed_relation)                  :: view

    view % base => base
    call view % declare(base % name() // '^T', [base % domain(2), base % domain(1)])

  end function transpose_of

  !===================================================================!
  ! Returns whether a tuple is in the view: through the base, ends
  ! swapped
  !===================================================================!

  pure logical function view_has(this, tuple)

    class(transposed_relation), intent(in) :: this
    integer                   , intent(in) :: tuple(:)

    view_has = .false.
    if (size(tuple) /= 2) return

    view_has = this % base % has([tuple(2), tuple(1)])

  end function view_has

  !===================================================================!
  ! Returns the number of tuples: the base's
  !===================================================================!

  pure integer function view_num_tuples(this)

    class(transposed_relation), intent(in) :: this

    view_num_tuples = this % base % num_tuples()

  end function view_num_tuples

  !===================================================================!
  ! Returns the tuple table: the base's, its two rows swapped
  !===================================================================!

  pure subroutine view_tuples(this, table)

    class(transposed_relation), intent(in)  :: this
    integer, allocatable      , intent(out) :: table(:,:)

    integer, allocatable :: forward(:,:)

    call this % base % tuples(forward)
    allocate(table(2, size(forward, 2)))
    table(1, :) = forward(2, :)
    table(2, :) = forward(1, :)

  end subroutine view_tuples

  !===================================================================!
  ! Returns the image of a member in the view: the base's preimage
  !===================================================================!

  function view_image_view(this, member) result(fibre)

    class(transposed_relation), target, intent(in) :: this
    integer                           , intent(in) :: member
    type(integer_fibre)                            :: fibre

    fibre = this % base % preimage_view(member)

  end function view_image_view

  !===================================================================!
  ! Returns the preimage of a member in the view: the base's image
  !===================================================================!

  function view_preimage_view(this, member) result(fibre)

    class(transposed_relation), target, intent(in) :: this
    integer                           , intent(in) :: member
    type(integer_fibre)                            :: fibre

    fibre = this % base % image_view(member)

  end function view_preimage_view

  !===================================================================!
  ! Returns true: a CSR relation is materialized
  !===================================================================!

  pure logical function csr_materialized(this)

    class(csr_relation), intent(in) :: this

    csr_materialized = .true.

  end function csr_materialized

  !===================================================================!
  ! Groups a finite family of (key, value) pairs by key: the fibres of
  ! a stored binary relation over one slot, as the compressed rows
  ! ptr(k) .. ptr(k+1)-1 into grouped(:). One counting pass, one prefix
  ! sum, one scatter, stable, so input order is retained within each
  ! key. A key outside 1..nkeys is skipped: a pair with no key belongs
  ! to no fibre. This is the one grouping kernel in the codebase; CSR
  ! builds, incidence lists, padded transposes, and triple combination
  ! are its callers.
  !
  !      keys      2  1  2  3            ptr       1   2      4   5
  !      values    a  b  c  d    -->     grouped   b | a  c | d
  !                                      key       1   2      3
  !===================================================================!

  pure subroutine group_by_key(nkeys, keys, values, ptr, grouped)

    integer             , intent(in)  :: nkeys
    integer             , intent(in)  :: keys(:)
    integer             , intent(in)  :: values(:)
    integer, allocatable, intent(out) :: ptr(:)
    integer, allocatable, intent(out) :: grouped(:)

    integer, allocatable :: cursor(:)
    integer              :: j, k, n

    allocate(ptr(nkeys + 1))
    ptr = 0
    do j = 1, size(keys)
       if (keys(j) >= 1 .and. keys(j) <= nkeys) then
          ptr(keys(j) + 1) = ptr(keys(j) + 1) + 1
       end if
    end do

    ptr(1) = 1
    do k = 1, nkeys
       ptr(k + 1) = ptr(k + 1) + ptr(k)
    end do

    n = ptr(nkeys + 1) - 1
    allocate(grouped(max(n, 0)))
    allocate(cursor(nkeys))
    cursor = ptr(1:nkeys)
    do j = 1, size(keys)
       if (keys(j) >= 1 .and. keys(j) <= nkeys) then
          grouped(cursor(keys(j))) = values(j)
          cursor(keys(j)) = cursor(keys(j)) + 1
       end if
    end do

  end subroutine group_by_key

  !===================================================================!
  ! Transposes a binary relation stored as padded lists: forward(k,
  ! key) lists the values of each key with per-key counts; the result
  ! lists, for each value 1..n_values, the keys that list it, in the
  ! same padded shape. One grouping, then the pad.
  !
  !      forward   key 1: 2 3     key 2: 3        (cells to faces)
  !      reverse   value 2: 1     value 3: 1 2    (faces to cells)
  !===================================================================!

  pure subroutine transpose_padded(forward, num_forward, n_values, &
       & reverse, num_reverse)

    integer             , intent(in)  :: forward(:,:)
    integer             , intent(in)  :: num_forward(:)
    integer             , intent(in)  :: n_values
    integer, allocatable, intent(out) :: reverse(:,:)
    integer, allocatable, intent(out) :: num_reverse(:)

    integer, allocatable :: keys(:), values(:), ptr(:), grouped(:)
    integer              :: key, k, n

    allocate(keys(sum(num_forward)), values(sum(num_forward)))
    n = 0
    do key = 1, size(num_forward)
       do k = 1, num_forward(key)
          n = n + 1
          keys(n)   = forward(k, key)
          values(n) = key
       end do
    end do

    call group_by_key(n_values, keys, values, ptr, grouped)
    call ragged_padded(ragged(ptr, grouped), reverse, num_reverse)

  end subroutine transpose_padded

  !===================================================================!
  ! Creates the ragged list from its padded shape: the lists copied
  ! once, in order. Counts outside 0..the width stop the program.
  !===================================================================!

  pure function ragged_of_padded(x, num_x) result(this)

    integer, intent(in) :: x(:,:)
    integer, intent(in) :: num_x(:)
    type(ragged)        :: this

    integer :: k, at

    if (size(x, 2) /= size(num_x)) error stop 'ragged: x and num_x must name the same number of lists'
    if (any(num_x < 0) .or. any(num_x > size(x, 1))) then
       error stop 'ragged: every count in num_x must lie within 0..size(x,1)'
    end if

    allocate(this % first(size(num_x) + 1), this % entries(sum(num_x)))
    at = 1
    do k = 1, size(num_x)
       this % first(k) = at
       this % entries(at : at + num_x(k) - 1) = x(1:num_x(k), k)
       at = at + num_x(k)
    end do
    this % first(size(num_x) + 1) = at

  end function ragged_of_padded

  !===================================================================!
  ! Creates the ragged list from its two arrays, copied once. first
  ! must start at 1 and end one past the last entry.
  !===================================================================!

  pure function ragged_of_lists(first, entries) result(this)

    integer, intent(in) :: first(:)
    integer, intent(in) :: entries(:)
    type(ragged)        :: this

    if (size(first) < 1) error stop 'ragged: first must have at least one entry, the one past the last list'
    if (first(1) /= 1 .or. first(size(first)) /= size(entries) + 1) then
       error stop 'ragged: first must start at 1 and end one past the last entry'
    end if

    this % first   = first
    this % entries = entries

  end function ragged_of_lists

  !===================================================================!
  ! Returns the number of lists
  !===================================================================!

  pure integer function ragged_num_lists(this)

    class(ragged), intent(in) :: this

    ragged_num_lists = size(this % first) - 1

  end function ragged_num_lists

  !===================================================================!
  ! Returns the length of list k
  !===================================================================!

  pure integer function ragged_length(this, k)

    class(ragged), intent(in) :: this
    integer      , intent(in) :: k

    ragged_length = this % first(k + 1) - this % first(k)

  end function ragged_length

  !===================================================================!
  ! Returns list k, copied
  !===================================================================!

  pure function ragged_list(this, k) result(members)

    class(ragged), intent(in) :: this
    integer      , intent(in) :: k
    integer, allocatable      :: members(:)

    members = this % entries(this % first(k) : this % first(k + 1) - 1)

  end function ragged_list

  !===================================================================!
  ! Returns the padded shape: the widest list's width, a count per
  ! list, and zeros past each count
  !===================================================================!

  pure subroutine ragged_padded(this, x, num_x)

    class(ragged)       , intent(in)  :: this
    integer, allocatable, intent(out) :: x(:,:)
    integer, allocatable, intent(out) :: num_x(:)

    integer :: k, n

    n = this % num_lists()
    allocate(num_x(n))
    do k = 1, n
       num_x(k) = this % length(k)
    end do
    allocate(x(max(maxval(num_x), 0), n))
    x = 0
    do k = 1, n
       x(1:num_x(k), k) = this % entries(this % first(k) : this % first(k + 1) - 1)
    end do

  end subroutine ragged_padded

end module relation_binary
