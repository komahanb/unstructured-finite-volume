!=====================================================================!
! LEVEL 1 OF THE NEW TOWER . THE BINARY SPECIALIZATION
!
! Arity two earns its own contract (AGENTS.md 5.3) because arity two
! has its own canonical questions. For
!
!      R  <=  A x B
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
!      (R^T)^T = R
!
! is a statement about EXTENSION - the same tuples over the same
! domains - not about stamps. Test it by comparing tuples and
! judging domains slot against slot; only a deliberate
! canonicalization could ever promise it by identity, and none is
! promised here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_binary_relation

  use graph_carrier , only : member_set
  use graph_relation, only : relation, slot

  implicit none

  private
  public :: binary_relation, csr_relation, transposed_view, transpose_of

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

  type, extends(binary_relation) :: csr_relation

     type(slot), allocatable, private :: signature(:)

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

  function source(this) result(domain)

    class(binary_relation), intent(in) :: this
    class(member_set), allocatable     :: domain

    domain = this % domain(1)

  end function source

  function target(this) result(domain)

    class(binary_relation), intent(in) :: this
    class(member_set), allocatable     :: domain

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

  type(csr_relation) function create_csr(name, source, target, table) &
       & result(this)

    character(len=*) , intent(in) :: name
    class(member_set), intent(in) :: source
    class(member_set), intent(in) :: target
    integer          , intent(in) :: table(:,:)

    integer, allocatable :: aloc(:), bloc(:), order(:), fill(:), stamp(:)
    integer, allocatable :: keepa(:), keepb(:)
    integer              :: na, nb, nt
    integer              :: j, p, q, row, col, kept

    if (.not. source % same_as(source) .or. &
         & .not. target % same_as(target)) then
       error stop 'graph_binary_relation: a signature refers to declared domains only'
    end if

    if (size(table, 1) /= 2) then
       error stop 'graph_binary_relation: each tuple has exactly one part per slot'
    end if

    na = source % size()
    nb = target % size()
    nt = size(table, 2)

    ! Validate through the carriers' own membership, and read every
    ! member's row through the carriers' own inverse enumeration.
    allocate(aloc(nt), bloc(nt))
    do j = 1, nt
       if (.not. source % has(table(1, j)) .or. &
            & .not. target % has(table(2, j))) then
          error stop 'graph_binary_relation: a tuple names a member its domain does not hold'
       end if
       aloc(j) = source % local_index(table(1, j))
       bloc(j) = target % local_index(table(2, j))
    end do

    ! Forward: count per source row, prefix-sum, place - duplicates
    ! and all - then collapse each row with one stamp sweep.
    allocate(this % xfwd(na + 1), fill(na), order(nt))
    this % xfwd = 0
    do j = 1, nt
       this % xfwd(aloc(j) + 1) = this % xfwd(aloc(j) + 1) + 1
    end do
    this % xfwd(1) = 1
    do row = 1, na
       this % xfwd(row + 1) = this % xfwd(row + 1) + this % xfwd(row)
    end do
    fill = this % xfwd(1:na) - 1
    do j = 1, nt
       fill(aloc(j)) = fill(aloc(j)) + 1
       order(fill(aloc(j))) = j
    end do

    allocate(this % tgt(nt), stamp(max(nb, 1)))
    allocate(keepa(nt), keepb(nt))
    stamp = 0
    kept  = 0
    do row = 1, na
       p = this % xfwd(row)
       q = this % xfwd(row + 1) - 1
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

    ! Backward: the same construction, mirrored, over the kept
    ! tuples alone.
    allocate(this % xbwd(nb + 1), this % src(kept))
    deallocate(fill); allocate(fill(max(nb, 1)))
    this % xbwd = 0
    do j = 1, kept
       this % xbwd(keepb(j) + 1) = this % xbwd(keepb(j) + 1) + 1
    end do
    this % xbwd(1) = 1
    do row = 1, nb
       this % xbwd(row + 1) = this % xbwd(row + 1) + this % xbwd(row)
    end do
    fill(1:nb) = this % xbwd(1:nb) - 1
    do j = 1, kept
       fill(keepb(j)) = fill(keepb(j)) + 1
       this % src(fill(keepb(j))) = keepa(j)
    end do

    allocate(this % signature(2))
    allocate(this % signature(1) % carrier, source=source)
    allocate(this % signature(2) % carrier, source=target)

    call this % declare(name)

  end function create_csr

  function csr_domain(this, slot_index) result(domain)

    class(csr_relation), intent(in) :: this
    integer            , intent(in) :: slot_index
    class(member_set), allocatable  :: domain

    allocate(domain, source=this % signature(slot_index) % carrier)

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

    row = this % signature(1) % carrier % local_index(member)
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

    row = this % signature(2) % carrier % local_index(member)
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

    row = this % signature(1) % carrier % local_index(tuple(1))
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
       a = this % signature(1) % carrier % member(row)
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

  function view_domain(this, slot_index) result(domain)

    class(transposed_view), intent(in) :: this
    integer               , intent(in) :: slot_index
    class(member_set), allocatable     :: domain

    domain = this % base % domain(3 - slot_index)

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

end module graph_binary_relation
