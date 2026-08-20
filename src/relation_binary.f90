!=====================================================================!
! LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION
!
! Arity two earns its own contract (AGENTS.md 5.3) because arity two
! has its own canonical questions. For
!
!      P  <=  A x B
!
! the specialization adds what no general arity can promise:
!
!      source, target      the two ends of the signature, named
!      image(a)            the members of B that a relates to
!      preimage(b)         the members of A that relate to b
!
! and the transpose is the canonical slot permutation, delivered
! here as a VIEW.
!
!                    MEMBERS IN, MEMBERS OUT
!
! Every query speaks MEMBER VALUES, never storage rows. A sparse
! carrier may hold { 10 20 30 }; image(20) is asked with 20 and
! answers members of the far side. The bridge between a member and
! its storage row is the carrier's own local_index - the inverse
! enumeration the carriers guarantee - so the indexed lookup below
! never assumes a domain is 1..n. The image of an outsider is the
! empty set: relating to nothing is an answer, not an error.
!
!                    TWO TIERS OF TRAVERSAL
!
! The deferred primitives are the VIEWS: image_view and
! preimage_view answer a fibre as a pointer into the stored index -
! no allocation, no copy, the hot-loop road (AGENTS.md 33). The
! allocating image and preimage stand above them as conveniences,
! written once for the whole family as copies of the views. A
! caller holding a view holds a borrow: it lives while the relation
! lives, and no longer.
!
!                        THE CSR CITIZEN
!
! csr_relation stores both directions, built once at construction:
!
!      xfwd, tgt      row a-local  ->  its B members     image
!      xbwd, src      row b-local  ->  its A members     preimage
!
! so each fibre is one row slice, has([a,b]) is one row scan, and
! construction - validation, duplicate collapse, both index builds -
! is linear in members plus tuples. Set semantics hold here exactly
! as in the stored relation: a tuple handed in twice is in the
! relation once, first appearance keeping its place.
!
! COMPLEXITY, PARAMETERIZED HONESTLY. Every fibre first asks the
! carrier where the member stands, so the true cost is
!
!      T_image(a)  =  T_local_index(a) + O(deg a)
!
! and the slice alone is O(deg). A counted carrier answers
! local_index in O(1), so the promise collapses to O(deg) there -
! the mesh path's case. A carrier whose local_index scans (the
! listed fixture does) pays its scan on every query; if such a
! carrier ever matters at scale, it owes itself an index.
!
!                     THE VIEW, AND ITS DEBT
!
! transpose_of(r) answers a lightweight view: O(1) to make, no
! topology copied, image and preimage swapped, the signature read
! in reverse. A view BORROWS - it holds its base by pointer, and
! the base must outlive it; the caller's base must carry the target
! attribute. That is the whole cost of an O(1) transpose in a
! language of value semantics, and it is stated rather than hidden.
!
! THE OWNERSHIP POLICY, DECLARED FOR THE LEVELS ABOVE. When the
! graph arrives and contains relations (AGENTS.md 14), the law is:
!
!      the graph OWNS stable relations;
!      views and fibre borrows may BORROW them.
!
! A graph accessor must therefore hand out its relations by
! reference to owned, stable storage - never as temporary copies
! that a view or fibre could dangle into. Whoever owns the base
! decides its lifetime; every borrower lives strictly inside it.
!
!                  IDENTITY IS NOT EQUALITY
!
! same_as answers minted identity: a view signs its own token, so a
! view is never same_as its base, and the involution
!
!      (P^T)^T = P
!
! is a statement about EXTENSION - the same tuples over the same
! domains - not about stamps. Test it by comparing tuples and
! judging domains slot against slot; only a deliberate
! canonicalization could ever promise it by identity, and none is
! promised here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module relation_binary

  use graph_fractal           , only : set_graph => graph
  use relation_finitary          , only : relation
  use map_set           , only : set_map
  use map_set_representation, only : set_representation
  use map_label         , only : label_map

  implicit none

  private
  public :: binary_relation, csr_relation, transposed_view, transpose_of
  public :: inclusion_of
  public :: group_by_key
  public :: transpose_padded

  !===================================================================!
  ! The abstract binary relation: the general contract, plus the
  ! questions only arity two can ask. Arity is answered here, once,
  ! for every descendant: two.
  !===================================================================!

  type, abstract, extends(relation) :: binary_relation

   contains

     procedure :: arity => binary_arity

     !----------------------------------------------------------------!
     ! The deferred primitives: fibres as borrows, no allocation.
     !----------------------------------------------------------------!

     procedure(binary_fibre_view_interface), deferred :: image_view
     procedure(binary_fibre_view_interface), deferred :: preimage_view

     !----------------------------------------------------------------!
     ! The conveniences, written once for the family: copies of the
     ! views, for callers who would rather own than borrow.
     !----------------------------------------------------------------!

     procedure :: image
     procedure :: preimage

     procedure :: source
     procedure :: target

  end type binary_relation

  abstract interface

     function binary_fibre_view_interface(this, member) result(fibre)
       import binary_relation
       class(binary_relation), target, intent(in) :: this
       integer               , intent(in)         :: member
       integer, pointer                           :: fibre(:)
     end function binary_fibre_view_interface

  end interface

  !===================================================================!
  ! The CSR citizen: both directions materialized once, every query
  ! an O(degree) slice.
  !===================================================================!

  !===================================================================!
  ! TWO QUESTIONS, TWO STORES, AND THEY DO NOT MEET.
  !
  !     signature      WHICH domains        semantic, identities
  !     coordinates    WHICH ROW a member   compiled, numbering only
  !
  ! The signature answers domain(k) and nothing else; the coordinates
  ! answer local_index and nothing else. A coordinate representation
  ! carries NO identity - it cannot say which set it numbers, and does
  ! not need to, because the signature beside it already did.
  !
  ! The coordinates are held BY VALUE, copied out of the caller's set
  ! map at construction. That is deliberate and is not the field's
  ! situation: many fields share one domain, so copying an extent per
  ! field was measured harmful; a CSR relation's row numbering is part
  ! of its own compiled execution contract, and the hot path may not go
  ! looking for it. So it is here, and image/preimage/has reach it
  ! directly - no map row scan, no graph traversal, no label lookup.
  !===================================================================!

  type, extends(binary_relation) :: csr_relation

     ! Semantic: which two domains.
     type(set_graph), allocatable, private :: signature(:)

     ! Compiled: how a member value becomes a row, each direction.
     class(set_representation), allocatable, private :: source_coords
     class(set_representation), allocatable, private :: target_coords

     integer             , private :: nnz = 0
     integer, allocatable, private :: xfwd(:), tgt(:)
     integer, allocatable, private :: xbwd(:), src(:)

   contains

     procedure :: domain        => csr_domain
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

  !===================================================================!
  ! The transpose view: a borrower. It holds its base by pointer
  ! and answers every question through it, ends swapped. The base
  ! must outlive the view.
  !===================================================================!

  type, extends(binary_relation) :: transposed_view

     class(binary_relation), pointer, private :: base => null()

   contains

     procedure :: domain        => view_domain
     procedure :: has           => view_has
     procedure :: num_tuples    => view_num_tuples
     procedure :: tuples        => view_tuples
     procedure :: image_view    => view_image_view
     procedure :: preimage_view => view_preimage_view

     ! No materialized binding, on purpose: the root's fail-closed
     ! default already answers false, and a borrower - copying it
     ! copies a pointer to a base it does not keep alive - is
     ! exactly what the default guards against. Views live OVER
     ! graph-owned relations, never inside them.

  end type transposed_view

contains

  !===================================================================!
  ! Arity two, for every binary citizen, forever.
  !===================================================================!

  pure integer function binary_arity(this)

    class(binary_relation), intent(in) :: this

    binary_arity = 2

  end function binary_arity

  !===================================================================!
  ! The two ends of the signature, named as arity two names them.
  !===================================================================!

  type(set_graph) function source(this) result(domain)

    class(binary_relation), intent(in) :: this

    domain = this % domain(1)

  end function source

  type(set_graph) function target(this) result(domain)

    class(binary_relation), intent(in) :: this

    domain = this % domain(2)

  end function target

  !===================================================================!
  ! The conveniences: own a copy of what the view borrows. Written
  ! once, here, for every binary citizen present and future.
  !===================================================================!

  subroutine image(this, member, indices)

    class(binary_relation), target, intent(in)  :: this
    integer                       , intent(in)  :: member
    integer, allocatable          , intent(out) :: indices(:)

    indices = this % image_view(member)

  end subroutine image

  subroutine preimage(this, member, indices)

    class(binary_relation), target, intent(in)  :: this
    integer                       , intent(in)  :: member
    integer, allocatable          , intent(out) :: indices(:)

    indices = this % preimage_view(member)

  end subroutine preimage

  !===================================================================!
  ! Declare a CSR relation: a name, the two carriers - any
  ! concretions, each judged by its own laws - and the tuple table,
  ! one column per tuple, members throughout. Refusals first, as at
  ! every gate of the level; then the duplicate collapse and both
  ! index builds, all linear:
  !
  !      count rows        one pass with local_index
  !      place tuples      counting sort by source row
  !      collapse          one stamp array over target rows
  !      backward build    the same, mirrored
  !===================================================================!

  type(csr_relation) function create_csr(name, source, target, table, sets) &
       & result(this)

    character(len=*), intent(in) :: name
    type(set_graph) , intent(in) :: source
    type(set_graph) , intent(in) :: target
    integer         , intent(in) :: table(:,:)
    type(set_map)   , intent(in) :: sets

    integer, allocatable :: aloc(:), bloc(:), order(:), stamp(:)
    integer, allocatable :: keepa(:), keepb(:)
    integer              :: na, nb, nt
    integer              :: j, p, q, row, col, kept

    if (.not. source % same_as(source) .or. &
         & .not. target % same_as(target)) then
       error stop 'relation_binary: a signature refers to declared domains only'
    end if

    if (size(table, 1) /= 2) then
       error stop 'relation_binary: each tuple has exactly one part per position'
    end if

    !----------------------------------------------------------------!
    ! COMPILATION. The map is read here and only here: the two extents
    ! are copied in, and from this line on the relation numbers its own
    ! rows. Nothing below, and nothing in image, preimage or has, asks
    ! the map anything.
    !----------------------------------------------------------------!

    call sets % extent_of(source, this % source_coords)
    call sets % extent_of(target, this % target_coords)

    na = this % source_coords % size()
    nb = this % target_coords % size()
    nt = size(table, 2)

    ! Validate through the coordinates' own membership, and read every
    ! member's row through their own inverse enumeration.
    allocate(aloc(nt), bloc(nt))
    do j = 1, nt
       if (.not. this % source_coords % has(table(1, j)) .or. &
            & .not. this % target_coords % has(table(2, j))) then
          error stop 'relation_binary: a tuple names a member its domain does not hold'
       end if
       aloc(j) = this % source_coords % local_index(table(1, j))
       bloc(j) = this % target_coords % local_index(table(2, j))
    end do

    ! Forward: group the tuples by source row - duplicates and all -
    ! then collapse each row with one stamp sweep.
    block
      integer, allocatable :: ptr(:), identity(:)
      allocate(identity(nt))
      identity = [(j, j = 1, nt)]
      call group_by_key(na, aloc, identity, ptr, order)

      allocate(this % xfwd(na + 1))
      allocate(this % tgt(nt), stamp(max(nb, 1)))
      allocate(keepa(nt), keepb(nt))
      stamp = 0
      kept  = 0
      do row = 1, na
         p = ptr(row)
         q = ptr(row + 1) - 1
         this % xfwd(row) = kept + 1
         do j = p, q
            col = bloc(order(j))
            if (stamp(col) /= row) then
               stamp(col)       = row
               kept             = kept + 1
               this % tgt(kept) = table(2, order(j))
               keepa(kept)      = table(1, order(j))
               keepb(kept)      = col
            end if
         end do
      end do
      this % xfwd(na + 1) = kept + 1
      this % nnz          = kept
    end block

    ! Backward: the kept tuples grouped by target row.
    call group_by_key(nb, keepb(1:kept), keepa(1:kept), &
         & this % xbwd, this % src)

    allocate(this % signature(2))
    this % signature(1) = source
    this % signature(2) = target

    call this % declare(name)

  end function create_csr

  type(set_graph) function csr_domain(this, position) result(domain)

    class(csr_relation), intent(in) :: this
    integer            , intent(in) :: position

    domain = this % signature(position)

  end function csr_domain

  !===================================================================!
  ! The fibre views: one local_index, one slice, zero allocation.
  ! Members in, members out; an outsider's fibre is the empty
  ! borrow. Cost, honestly: T_local_index(member) + O(1) to make,
  ! O(degree) to read.
  !===================================================================!

  function csr_image_view(this, member) result(fibre)

    class(csr_relation), target, intent(in) :: this
    integer                    , intent(in) :: member
    integer, pointer                        :: fibre(:)

    integer :: row

    row = this % source_coords % local_index(member)
    if (row == 0) then
       fibre => this % tgt(1:0)
       return
    end if

    fibre => this % tgt(this % xfwd(row) : this % xfwd(row + 1) - 1)

  end function csr_image_view

  function csr_preimage_view(this, member) result(fibre)

    class(csr_relation), target, intent(in) :: this
    integer                    , intent(in) :: member
    integer, pointer                        :: fibre(:)

    integer :: row

    row = this % target_coords % local_index(member)
    if (row == 0) then
       fibre => this % src(1:0)
       return
    end if

    fibre => this % src(this % xbwd(row) : this % xbwd(row + 1) - 1)

  end function csr_preimage_view

  !===================================================================!
  ! Membership: one row, one scan - O(degree), as promised.
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

  pure integer function csr_num_tuples(this)

    class(csr_relation), intent(in) :: this

    csr_num_tuples = this % nnz

  end function csr_num_tuples

  !===================================================================!
  ! The tuple table, rebuilt row by row - non-hot generic access.
  !===================================================================!

  pure subroutine csr_tuples(this, table)

    class(csr_relation), intent(in)   :: this
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
  ! Make the transpose view: O(1), nothing copied, a fresh identity
  ! of its own. The base must be a target, and must outlive the
  ! view - a view borrows, it never owns.
  !===================================================================!

  function transpose_of(base) result(view)

    class(binary_relation), target, intent(in) :: base
    type(transposed_view)                      :: view

    view % base => base
    call view % declare(base % name() // '^T')

  end function transpose_of

  !===================================================================!
  ! The view's answers: everything through the base, ends swapped.
  !===================================================================!

  type(set_graph) function view_domain(this, position) result(domain)

    class(transposed_view), intent(in) :: this
    integer               , intent(in) :: position

    domain = this % base % domain(3 - position)

  end function view_domain

  pure logical function view_has(this, tuple)

    class(transposed_view), intent(in) :: this
    integer               , intent(in) :: tuple(:)

    view_has = .false.
    if (size(tuple) /= 2) return

    view_has = this % base % has([tuple(2), tuple(1)])

  end function view_has

  pure integer function view_num_tuples(this)

    class(transposed_view), intent(in) :: this

    view_num_tuples = this % base % num_tuples()

  end function view_num_tuples

  pure subroutine view_tuples(this, table)

    class(transposed_view), intent(in)  :: this
    integer, allocatable  , intent(out) :: table(:,:)

    integer, allocatable :: forward(:,:)

    call this % base % tuples(forward)
    allocate(table(2, size(forward, 2)))
    table(1, :) = forward(2, :)
    table(2, :) = forward(1, :)

  end subroutine view_tuples

  function view_image_view(this, member) result(fibre)

    class(transposed_view), target, intent(in) :: this
    integer                       , intent(in) :: member
    integer, pointer                           :: fibre(:)

    fibre => this % base % preimage_view(member)

  end function view_image_view

  function view_preimage_view(this, member) result(fibre)

    class(transposed_view), target, intent(in) :: this
    integer                       , intent(in) :: member
    integer, pointer                           :: fibre(:)

    fibre => this % base % image_view(member)

  end function view_preimage_view

  pure logical function csr_materialized(this)

    class(csr_relation), intent(in) :: this

    csr_materialized = .true.

  end function csr_materialized

  !===================================================================!
  ! The relational face of a subobject: the inclusion
  !
  !      I_S  <=  S x A ,       (s, s) for every s in S
  !
  ! total, functional and injective by construction - each member of
  ! S relates to exactly its own image in the ambient. What the
  ! subset states as membership, this states as a relation, ready
  ! for the algebra: restriction of a relation to a subset is a
  ! composition with its inclusion.
  !===================================================================!

  type(csr_relation) function inclusion_of(s, host, sets, labels) &
       & result(inclusion)

    type(set_graph) , intent(in) :: s
    type(set_graph) , intent(in) :: host
    type(set_map)   , intent(in) :: sets
    type(label_map) , intent(in) :: labels

    integer, allocatable :: table(:,:)
    integer              :: k, n

    n = sets % size_of(s)

    allocate(table(2, n))
    do k = 1, n
       table(:, k) = [sets % member_of(s, k), sets % member_of(s, k)]
    end do

    inclusion = csr_relation(labels % label_of(s) // ' in ' // &
         &                   labels % label_of(host), s, host, table, sets)

  end function inclusion_of

  !===================================================================!
  ! Group a finite family of (key, value) pairs by key: the fibers
  ! of a stored binary relation over one slot, as the compressed
  ! rows ptr(k) .. ptr(k+1)-1 into grouped(:). One counting pass,
  ! one prefix sum, one scatter - stable, so arrival order is kept
  ! within each key. A key outside 1..nkeys is skipped: a pair with
  ! no key belongs to no fiber. This is the one grouping kernel in
  ! the codebase; CSR builds, incidence lists, padded transposes,
  ! and triple combination are its callers.
  !===================================================================!

  pure subroutine group_by_key(nkeys, keys, values, ptr, grouped)

    integer             , intent(in)  :: nkeys
    integer             , intent(in)  :: keys(:)
    integer             , intent(in)  :: values(:)
    integer, allocatable, intent(out) :: ptr(:)
    integer, allocatable, intent(out) :: grouped(:)

    integer, allocatable :: cursor(:)
    integer :: j, k, n

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
  ! Transpose a binary relation stored as padded lists: forward(k,
  ! key) lists the values of each key with per-key counts; the
  ! result lists, for each value 1..n_values, the keys that touch
  ! it, in the same padded shape. One grouping, then the pad.
  !===================================================================!

  pure subroutine transpose_padded(forward, num_forward, n_values, &
       & reverse, num_reverse)

    integer             , intent(in)  :: forward(:,:)
    integer             , intent(in)  :: num_forward(:)
    integer             , intent(in)  :: n_values
    integer, allocatable, intent(out) :: reverse(:,:)
    integer, allocatable, intent(out) :: num_reverse(:)

    integer, allocatable :: keys(:), values(:), ptr(:), grouped(:)
    integer :: key, k, n, v

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

    allocate(num_reverse(n_values))
    do v = 1, n_values
       num_reverse(v) = ptr(v + 1) - ptr(v)
    end do
    allocate(reverse(maxval(num_reverse), n_values))
    reverse = 0
    do v = 1, n_values
       reverse(1:num_reverse(v), v) = grouped(ptr(v) : ptr(v + 1) - 1)
    end do

  end subroutine transpose_padded

end module relation_binary
