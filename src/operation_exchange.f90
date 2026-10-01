!=====================================================================!
! THE EXCHANGE: the distribution of the unknowns over the images and
! the operators it induces on vectors.
!
! Of the whole, an exchange serves vectors of the whole length: every
! image has every entry; those it owns are its own, and every other
! is stale until a gather. Of a part, it serves vectors of the part's
! length: the part the partitioner cuts from the coupling of the
! blocks, in the part's own numbering, the owned unknowns first and
! the halo after them. The exchange is
!
!      gather(x)       of the whole: x(i) := x(i) of the owner of i,
!                      for every i; of a part: every halo entry takes
!                      its owner's value
!      total(s)        the sum of s over the images
!      maximum(s)      the largest s over the images
!      part(coupling)  the part of this image, of the whole's
!                      ownership
!      restriction(x)  the part's entries of a whole vector
!      assembled(x)    the whole vector of the parts' owned entries
!      add_summed_rows(x, y)
!                      a summed row's owner adds to y the other parts'
!                      terms of the row at x
!
! The ownership is given per unknown and is the same for every
! unknown of a block, so the part is cut from the coupling of the
! blocks and consists of whole blocks. The neighbours of a part are
! the owners of its halo; the relation is symmetric, the adjacency of
! the blocks being so. The part sends to each neighbour the unknowns
! of its owned blocks adjacent to that neighbour's, one segment per
! neighbour in ascending image order, each segment in the order of
! the colouring of the whole blocks, then of the block; the halo it
! receives from a neighbour is that neighbour's segment for it, in
! the same order. Both ends compute the segments from the same
! ownership and the same adjacency, so nothing is exchanged to form
! them, and a gather of a part is one write of the sent segments to
! this image's outbox, one synchronisation with the neighbours, and
! one contiguous read of each neighbour's segment: the volume per
! image is that of its own halo, not of the interface of the whole,
! and no image waits for an image that is not its neighbour. The
! part also stores the colouring of the whole blocks, the same on
! every image, so that a relaxation sweeping the colours in turn
! gathers each colour's segments alone, and relaxes in the order one
! image does. A summed row reads every point - the integral of an
! unknown over the region - so it is not a coupling of the cut: its
! owner states the terms of its own columns, every other part keeps
! the terms of its owned columns here, and one sum over the images
! of their partial sums completes the row. Under one image nothing is
! exchanged.
!
! The outbox is a coarray allocated in this module, by every image
! at once with the same bounds. The coarray runtime in use
! (OpenCoarrays 2.9.2 under gfortran 14) mis-sizes a remote read
! whose source or destination is vector-subscripted, and rejects one
! of zero length: every remote read here is a contiguous section of
! positive length into a contiguous vector. The reductions need no
! coarray variable: the collectives are sound.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_exchange

  use util_precision       , only : dp
  use view_directed        , only : directed_graph
  use view_directed_stored , only : stored_directed_graph
  use relation_partition   , only : partition_relation
  use transform_partitioner, only : partitioner, PARTITION_ADOPTED
  use operation_traversal  , only : traversal, TRAVERSAL_COLOURING
  use field_calculus       , only : field

  implicit none

  private
  public :: exchange, blocked, state_summed_rows

  type :: exchange

     integer :: image = 1
     integer :: images = 1
     ! the unknowns of a block, which share their owner
     integer :: width = 1
     ! the whole's ownership, one owning image per unknown
     integer, allocatable :: owner(:)
     ! over the entries of the vectors this exchange serves
     logical, allocatable :: is_owned(:)
     integer, allocatable :: owned(:)
     ! of a part: the whole's index of each entry, owned first
     integer, allocatable :: whole_index(:)
     ! of a part: the colour of each of its blocks in the colouring of
     ! the whole blocks, and the number of colours
     integer, allocatable :: block_colour(:)
     integer :: num_colours = 0
     ! of a part: the neighbouring images, ascending; the entries sent,
     ! in the order of the outbox, those of neighbour k and colour c
     ! at sent_first(c, k):sent_first(c + 1, k) - 1; the halo entries
     ! received, in the order of the neighbours' segments, those of
     ! neighbour k and colour c at received_first(c, k):
     ! received_first(c + 1, k) - 1, read from the neighbour's outbox
     ! from remote_first(c, k)
     integer, allocatable :: neighbours(:)
     integer, allocatable :: sent(:), received(:)
     integer, allocatable :: sent_first(:,:), received_first(:,:), remote_first(:,:)
     ! of a part: the summed rows, each one's row in this part when it
     ! owns it (zero otherwise), and the terms of this part's owned
     ! columns in the rows it does not own
     integer , allocatable :: summed_row(:), term_row(:), term_column(:)
     real(dp), allocatable :: term_weight(:)

   contains

     procedure :: gather
     procedure :: total
     procedure :: maximum
     procedure :: part
     procedure :: restriction
     procedure :: assembled
     procedure :: restricted_triples
     procedure :: add_summed_rows

  end type exchange

  interface exchange
     module procedure create
  end interface exchange

  ! THE OUTBOX: the segments this image sends, in two buffers
  ! alternated by every gather of a part. A neighbour reads a buffer
  ! after the synchronisation of its gather; the next gather writes
  ! the other buffer, and the one after it passes a synchronisation
  ! the neighbour reaches only after its read, so no read is overtaken
  ! by a write. One vector, the second buffer after the first, each
  ! of the largest segment total over the images, the bounds being
  ! the same on every image.
  real(dp), allocatable :: outbox(:)[:]
  integer :: buffer = 2

contains

  !===================================================================!
  ! The exchange of the whole of the ownership given, one owning image
  ! per unknown, the unknowns in blocks of the width given. An owner
  ! outside 1..num_images(), or an ownership not in whole blocks, is
  ! invalid input.
  !===================================================================!

  function create(owner, width) result(this)

    integer, intent(in) :: owner(:)
    integer, intent(in), optional :: width
    type(exchange) :: this

    integer :: i
    character(len=250) :: message

    this % image  = this_image()
    this % images = num_images()
    if (present(width)) this % width = width
    if (any(owner < 1) .or. any(owner > this % images)) then
       write(message,'(a,i0,a,i0)') 'operation_exchange: every unknown is owned by one of the images &
            &1..', this % images, '; owners range to ', maxval(owner)
       error stop trim(message)
    end if
    if (mod(size(owner), this % width) /= 0) then
       write(message,'(a,i0,a,i0)') 'operation_exchange: the unknowns come in whole blocks; unknowns = ', &
            & size(owner), ', width = ', this % width
       error stop trim(message)
    end if
    this % owner    = owner
    this % is_owned = owner == this % image
    this % owned    = pack([(i, i = 1, size(owner))], this % is_owned)

  end function create

  !===================================================================!
  ! THE PART of this image, cut by the partitioner from the coupling
  ! of the blocks with the whole's ownership adopted: whole blocks,
  ! the owned first; and the segments sent to and received from each
  ! neighbour, from one pass over the adjacency of the blocks. Every
  ! image calls this at once: the outbox is allocated here. A
  ! coupling of another vertex count than the ownership, or a block
  ! whose unknowns have several owners, is invalid input.
  !===================================================================!

  function part(this, coupling) result(local)

    class(exchange)      , intent(in) :: this
    class(directed_graph), intent(in) :: coupling
    type(exchange) :: local

    type(stored_directed_graph) :: blocks
    type(partitioner) :: cut
    type(traversal) :: colouring
    class(field), allocatable :: coloured
    type(partition_relation) :: rel
    class(directed_graph), allocatable :: cut_part
    integer, allocatable :: block_owner(:), local_blocks(:), colours(:)
    integer :: nb, w, i, b, j, owned_blocks, length
    character(len=250) :: message

    w  = this % width
    nb = size(this % owner) / w
    if (coupling % num_vertices() /= size(this % owner)) then
       write(message,'(a,i0,a,i0)') 'operation_exchange: the coupling has one vertex per unknown; &
            &vertices = ', coupling % num_vertices(), ', unknowns = ', size(this % owner)
       error stop trim(message)
    end if
    block_owner = this % owner(1:nb * w:w)
    do b = 1, nb
       if (any(this % owner((b - 1) * w + 1:b * w) /= block_owner(b))) then
          write(message,'(a,i0,a)') 'operation_exchange: the unknowns of block ', b, ' have several owners'
          error stop trim(message)
       end if
    end do
    blocks = blocked(coupling, w)
    colouring = traversal(TRAVERSAL_COLOURING)
    call colouring % apply(blocks, output=coloured)
    call coloured % integer_vector(colours)

    ! this image's part alone
    cut = partitioner(PARTITION_ADOPTED, this % images, part=this % image, adopted=block_owner)
    call cut % partition_graph(blocks, cut_part, rel)
    owned_blocks = count(block_owner == this % image)
    local_blocks = [(rel % global_vertex_index(i), i = 1, cut_part % num_vertices())]

    local % image  = this % image
    local % images = this % images
    local % width  = w
    local % owner  = this % owner
    length = size(local_blocks) * w
    allocate(local % whole_index(length))
    do i = 1, size(local_blocks)
       local % whole_index((i - 1) * w + 1:i * w) = [((local_blocks(i) - 1) * w + j, j = 1, w)]
    end do
    local % is_owned     = this % owner(local % whole_index) == this % image
    local % block_colour = colours(local_blocks)
    local % num_colours  = maxval(colours)
    local % owned        = pack([(i, i = 1, length)], local % is_owned)

    call segments(local, blocks, block_owner, colours, local_blocks, owned_blocks)
    call allocate_outbox(size(local % sent))

  end function part

  !===================================================================!
  ! THE SEGMENTS of a part. The neighbours are the owners of the halo
  ! blocks. One pass over the blocks owned by this image or by a
  ! neighbour reads each block's adjacent owners: a block of this
  ! image adjacent to neighbour k is sent to k; a block of neighbour k
  ! adjacent to this image is received from k; a block of neighbour k
  ! adjacent to an image before this one lies before k's segment for
  ! this image in k's outbox. A halo block received from no neighbour,
  ! or a received block outside the part, is invalid input: the
  ! partitioner's halo is the set of blocks adjacent to the owned
  ! ones, which is what the pass reads.
  !===================================================================!

  subroutine segments(local, blocks, block_owner, colours, local_blocks, owned_blocks)

    type(exchange)             , intent(inout) :: local
    type(stored_directed_graph), intent(in)    :: blocks
    integer                    , intent(in)    :: block_owner(:), colours(:), local_blocks(:)
    integer                    , intent(in)    :: owned_blocks

    integer, allocatable :: neighbour_index(:), marked(:), adjacent(:)
    integer, allocatable :: send_count(:,:), receive_count(:,:), before(:)
    integer, allocatable :: send_at(:,:), receive_at(:,:), place_of(:)
    integer :: me, w, nc, nn, nb, b, o, i, k, c, m, pass, p, j, s, r

    me = local % image
    w  = local % width
    nc = local % num_colours
    nb = size(block_owner)

    ! the neighbours, ascending, from the owners of the halo blocks
    allocate(neighbour_index(local % images), source=0)
    do i = owned_blocks + 1, size(local_blocks)
       neighbour_index(block_owner(local_blocks(i))) = 1
    end do
    nn = 0
    do k = 1, local % images
       if (neighbour_index(k) == 0) cycle
       nn = nn + 1
       neighbour_index(k) = nn
    end do
    allocate(local % neighbours(nn))
    do k = 1, local % images
       if (neighbour_index(k) > 0) local % neighbours(neighbour_index(k)) = k
    end do

    allocate(send_count(nc, nn), receive_count(nc, nn), before(nn), marked(local % images))
    allocate(local % sent_first(nc + 1, nn), local % received_first(nc + 1, nn), local % remote_first(nc + 1, nn))
    allocate(send_at(nc, nn), receive_at(nc, nn))
    allocate(place_of(nb * w), source=0)
    place_of(local % whole_index) = [(i, i = 1, size(local % whole_index))]
    marked = 0

    ! the first pass counts the blocks of each segment, the second
    ! places their unknowns in the order of the counts
    do pass = 1, 2
       if (pass == 1) then
          send_count    = 0
          receive_count = 0
          before        = 0
       else
          call segment_bounds(send_count, w, local % sent_first)
          call segment_bounds(receive_count, w, local % received_first)
          do k = 1, nn
             local % remote_first(:, k) = before(k) + local % received_first(:, k) - local % received_first(1, k) + 1
          end do
          allocate(local % sent(w * sum(send_count)))
          allocate(local % received(w * sum(receive_count)))
          send_at    = local % sent_first(1:nc, :)
          receive_at = local % received_first(1:nc, :)
       end if
       do b = 1, nb
          o = block_owner(b)
          if (o /= me) then
             if (neighbour_index(o) == 0) cycle
          end if
          c = colours(b)
          call blocks % adjacent_vertices(b, adjacent)
          do i = 1, size(adjacent)
             m = block_owner(adjacent(i))
             ! each adjacent image once per block and pass
             if (m == o .or. marked(m) == (pass - 1) * nb + b) cycle
             marked(m) = (pass - 1) * nb + b
             if (o == me) then
                k = neighbour_index(m)
                if (pass == 1) then
                   send_count(c, k) = send_count(c, k) + 1
                else
                   do j = 1, w
                      s = send_at(c, k)
                      local % sent(s) = place_of((b - 1) * w + j)
                      send_at(c, k) = s + 1
                   end do
                end if
             else if (m == me) then
                k = neighbour_index(o)
                if (pass == 1) then
                   receive_count(c, k) = receive_count(c, k) + 1
                else
                   do j = 1, w
                      r = receive_at(c, k)
                      p = place_of((b - 1) * w + j)
                      if (p == 0) error stop 'operation_exchange: a block sent by a neighbour lies outside the part'
                      local % received(r) = p
                      receive_at(c, k) = r + 1
                   end do
                end if
             else if (m < me .and. pass == 1) then
                k = neighbour_index(o)
                before(k) = before(k) + w
             end if
          end do
       end do
    end do

    if (size(local % received) /= (size(local_blocks) - owned_blocks) * w) then
       error stop 'operation_exchange: the halo of the part is not the union of the neighbours'' segments'
    end if

  end subroutine segments

  ! the first entry of each (colour, neighbour) segment, the segments
  ! of one neighbour consecutive by colour and the neighbours
  ! consecutive, each block counted contributing w entries
  pure subroutine segment_bounds(count, w, first)
    integer, intent(in)  :: count(:,:), w
    integer, intent(out) :: first(:,:)
    integer :: c, k, n
    n = 1
    do k = 1, size(count, 2)
       do c = 1, size(count, 1)
          first(c, k) = n
          n = n + w * count(c, k)
       end do
       first(size(count, 1) + 1, k) = n
    end do
  end subroutine segment_bounds

  !===================================================================!
  ! The outbox of the length given on this image, allocated by every
  ! image at once with the largest length over the images; an outbox
  ! already of that length is left as it is.
  !===================================================================!

  subroutine allocate_outbox(length)

    integer, intent(in) :: length

    integer :: needed

    needed = max(length, 1)
    call co_max(needed)
    if (allocated(outbox)) then
       if (size(outbox) == 2 * needed) return
       deallocate(outbox)
    end if
    allocate(outbox(2 * needed)[*])
    outbox = 0.0_dp

  end subroutine allocate_outbox

  !===================================================================!
  ! THE GATHER. Of the whole, by one sum over the images of the vector
  ! with every entry not owned set to zero: each entry is owned by one
  ! image, so the sum is the owner's value everywhere. Of a part, by
  ! the neighbour exchange: the sent segments written to this image's
  ! outbox, one synchronisation with the neighbours, and each
  ! neighbour's segment for this image read from its outbox into the
  ! halo. Every image calls it, with the same colour.
  !===================================================================!

  subroutine gather(this, x, colour, summed_into, factor)

    class(exchange), intent(in)    :: this
    real(dp)       , intent(inout) :: x(:)
    ! of a part, the segments of this colour alone
    integer        , intent(in), optional :: colour
    ! of a part, the other parts' terms of the summed rows at x, times
    ! the factor, added to this vector
    real(dp)       , intent(inout), optional :: summed_into(:)
    real(dp)       , intent(in)   , optional :: factor

    real(dp), allocatable :: segment(:)
    integer :: c1, c2, k, a, b, ra, rb, shift

    if (this % images == 1) return
    if (.not. allocated(this % whole_index)) then
       where (.not. this % is_owned) x = 0.0_dp
       call co_sum(x)
       return
    end if
    c1 = 1
    c2 = this % num_colours
    if (present(colour)) then
       c1 = colour
       c2 = colour
    end if
    buffer = 3 - buffer
    shift  = (buffer - 1) * (size(outbox) / 2)
    allocate(segment(size(this % received)))
    do k = 1, size(this % neighbours)
       a = this % sent_first(c1, k)
       b = this % sent_first(c2 + 1, k) - 1
       if (b < a) cycle
       outbox(shift + a:shift + b) = x(this % sent(a:b))
    end do
    if (size(this % neighbours) > 0) sync images(this % neighbours)
    do k = 1, size(this % neighbours)
       a  = this % received_first(c1, k)
       b  = this % received_first(c2 + 1, k) - 1
       ! a section of no entries is not read: the runtime rejects
       ! a remote read of zero length
       if (b < a) cycle
       ra = shift + this % remote_first(c1, k)
       rb = ra + b - a
       ! into a contiguous vector first: the runtime mis-sizes a
       ! remote read whose destination is vector-subscripted
       segment(1:b - a + 1) = outbox(ra:rb)[this % neighbours(k)]
       x(this % received(a:b)) = segment(1:b - a + 1)
    end do
    if (present(summed_into)) call this % add_summed_rows(x, summed_into, factor)

  end subroutine gather

  !===================================================================!
  ! The part's entries of a whole vector, and the whole vector of the
  ! parts' owned entries, by one sum over the images.
  !===================================================================!

  function restriction(this, x) result(y)

    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: x(:)
    real(dp), allocatable :: y(:)

    y = x(this % whole_index)

  end function restriction

  function assembled(this, x) result(y)

    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: x(:)
    real(dp), allocatable :: y(:)

    allocate(y(size(this % owner)), source=0.0_dp)
    y(this % whole_index(this % owned)) = x(this % owned)
    if (this % images > 1) call co_sum(y)

  end function assembled

  !===================================================================!
  ! The triples of an operator on the whole restricted to the part:
  ! those of the rows it owns, in the part's numbering. A column of an
  ! owned row outside the part is invalid input: the part contains every
  ! block adjacent to its own.
  !===================================================================!

  subroutine restricted_triples(this, rows, columns, weights, part_rows, part_columns, part_weights, summed)

    class(exchange)      , intent(in)  :: this
    integer              , intent(in)  :: rows(:), columns(:)
    real(dp)             , intent(in)  :: weights(:)
    integer , allocatable, intent(out) :: part_rows(:), part_columns(:)
    real(dp), allocatable, intent(out) :: part_weights(:)
    ! the whole's summed rows: of those the terms of the owned columns
    ! alone, the others' being added by add_summed_rows
    integer              , intent(in), optional :: summed(:)

    integer, allocatable :: place_of(:)
    logical, allocatable :: owned_row(:), is_summed(:)
    integer :: i

    allocate(place_of(size(this % owner)), source=0)
    place_of(this % whole_index) = [(i, i = 1, size(this % whole_index))]
    owned_row = this % owner(rows) == this % image
    if (present(summed)) then
       allocate(is_summed(size(this % owner)), source=.false.)
       is_summed(summed) = .true.
       owned_row = owned_row .and. (.not. is_summed(rows) .or. this % owner(columns) == this % image)
    end if
    part_rows    = place_of(pack(rows, owned_row))
    part_columns = place_of(pack(columns, owned_row))
    part_weights = pack(weights, owned_row)
    if (any(part_columns == 0)) then
       error stop 'operation_exchange: a column of an owned row lies outside the part'
    end if

  end subroutine restricted_triples

  !===================================================================!
  ! THE SUMMED ROWS of a part, from the whole's rows and the whole
  ! operator's terms (rows, columns, weights): the terms this part owns
  ! the column of, in the rows it does not own.
  !===================================================================!

  subroutine state_summed_rows(this, summed, rows, columns, weights)

    type(exchange), intent(inout) :: this
    integer       , intent(in)    :: summed(:), rows(:), columns(:)
    real(dp)      , intent(in)    :: weights(:)

    integer, allocatable :: place_of(:), index_of(:)
    logical, allocatable :: from_here(:)
    integer :: i

    allocate(place_of(size(this % owner)), index_of(size(this % owner)), source=0)
    place_of(this % whole_index) = [(i, i = 1, size(this % whole_index))]
    index_of(summed) = [(i, i = 1, size(summed))]
    this % summed_row = merge(place_of(summed), 0, this % owner(summed) == this % image)
    from_here = index_of(rows) > 0 .and. this % owner(rows) /= this % image .and. this % owner(columns) == this % image
    this % term_row    = index_of(pack(rows, from_here))
    this % term_column = place_of(pack(columns, from_here))
    this % term_weight = pack(weights, from_here)

  end subroutine state_summed_rows

  !===================================================================!
  ! The other parts' terms of the summed rows at x added to y at the
  ! rows this part owns: one sum over the images of every part's
  ! partial sums. Every image calls it, the owner of no summed row too.
  !===================================================================!

  subroutine add_summed_rows(this, x, y, factor)

    class(exchange), intent(in)    :: this
    real(dp)       , intent(in)    :: x(:)
    real(dp)       , intent(inout) :: y(:)
    ! the terms times this factor, one when absent
    real(dp)       , intent(in), optional :: factor

    real(dp), allocatable :: partial(:)
    integer :: e

    if (.not. allocated(this % summed_row)) return
    if (size(this % summed_row) == 0) return
    allocate(partial(size(this % summed_row)), source=0.0_dp)
    do e = 1, size(this % term_row)
       partial(this % term_row(e)) = partial(this % term_row(e)) + this % term_weight(e) * x(this % term_column(e))
    end do
    call co_sum(partial)
    if (present(factor)) partial = factor * partial
    do e = 1, size(this % summed_row)
       if (this % summed_row(e) > 0) y(this % summed_row(e)) = y(this % summed_row(e)) + partial(e)
    end do

  end subroutine add_summed_rows

  !===================================================================!
  ! The coupling of the blocks of `width` consecutive unknowns: an
  ! edge between two blocks for every edge of the unknowns' coupling
  ! that joins them.
  !===================================================================!

  function blocked(coupling, width) result(blocks)

    class(directed_graph), intent(in) :: coupling
    integer              , intent(in) :: width
    type(stored_directed_graph) :: blocks

    integer, allocatable :: tails(:), heads(:)
    integer :: e, ne, t, h, n

    ne = coupling % num_edges()
    allocate(tails(ne), heads(ne))
    n = 0
    do e = 1, ne
       if (.not. coupling % edge_has_head(e)) cycle
       t = (coupling % edge_tail(e) - 1) / width + 1
       h = (coupling % edge_head(e) - 1) / width + 1
       if (t == h) cycle
       n = n + 1
       tails(n) = t
       heads(n) = h
    end do
    blocks = stored_directed_graph(coupling % num_vertices() / width, tails=tails(1:n), heads=heads(1:n))

  end function blocked

  real(dp) function total(this, s)
    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: s
    total = s
    if (this % images > 1) call co_sum(total)
  end function total

  real(dp) function maximum(this, s)
    class(exchange), intent(in) :: this
    real(dp)       , intent(in) :: s
    maximum = s
    if (this % images > 1) call co_max(maximum)
  end function maximum

end module operation_exchange
