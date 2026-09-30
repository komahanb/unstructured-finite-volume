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
!
! The ownership is given per unknown and is the same for every
! unknown of a block, so the part is cut from the coupling of the
! blocks and consists of whole blocks. The interface is the union of the
! halos of every part: each image cuts every part from the same
! ownership, so no exchange forms it. A gather of a part sums the
! interface alone over the images: the halo a product or a relaxation
! over owned rows reads, at a length that falls against the whole as
! the unknowns grow. The part also stores the colouring of the whole
! blocks, the same on every image, and its interface in colour order,
! so that a relaxation sweeping the colours in turn gathers each
! colour's interface alone, and relaxes in the order one image does.
! Under one image nothing is exchanged.
!
! A reduction needs no coarray variable: the coarray runtime in use
! (OpenCoarrays 2.9.2 under gfortran 14) registers a module coarray
! before the runtime is initialised and mis-sizes a vector-subscript
! access through a component, while its collectives are sound.
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
  public :: exchange, blocked

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
     ! of a part: the interface entries this part owns and those of its
     ! halo, each with its place in the interface
     integer, allocatable :: sent(:), sent_at(:), received(:), received_at(:)
     integer :: interface_length = 0
     ! of a part: the colour of each of its blocks in the colouring of
     ! the whole blocks; the interface of colour c at
     ! colour_first(c):colour_first(c+1)-1, and the entries sent to it
     ! and received from it at sent_first(c) and received_first(c)
     integer, allocatable :: block_colour(:), colour_first(:), sent_first(:), received_first(:)

   contains

     procedure :: gather
     procedure :: total
     procedure :: maximum
     procedure :: part
     procedure :: restriction
     procedure :: assembled
     procedure :: restricted_triples

  end type exchange

  interface exchange
     module procedure create
  end interface exchange

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
  ! the owned first; and the interface, from the halo of every part.
  ! A coupling of another vertex count than the ownership, or a block
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
    integer, allocatable :: block_owner(:), interface_of(:), local_blocks(:), position(:), colours(:)
    logical, allocatable :: in_interface(:)
    integer, allocatable :: order(:)
    integer :: nb, w, k, i, b, j, c, owned_blocks, length
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

    ! every part's halo, from the same ownership on every image
    allocate(in_interface(nb), source=.false.)
    do k = 1, this % images
       cut = partitioner(PARTITION_ADOPTED, this % images, part=k, adopted=block_owner)
       call cut % partition_graph(blocks, cut_part, rel)
       owned_blocks = count(block_owner == k)
       do i = owned_blocks + 1, cut_part % num_vertices()
          in_interface(rel % global_vertex_index(i)) = .true.
       end do
       if (k == this % image) local_blocks = [(rel % global_vertex_index(i), i = 1, cut_part % num_vertices())]
    end do

    local % image  = this % image
    local % images = this % images
    local % width  = w
    local % owner  = this % owner
    length = size(local_blocks) * w
    allocate(local % whole_index(length))
    do i = 1, size(local_blocks)
       local % whole_index((i - 1) * w + 1:i * w) = [((local_blocks(i) - 1) * w + j, j = 1, w)]
    end do
    local % is_owned = this % owner(local % whole_index) == this % image
    local % block_colour = colours(local_blocks)
    local % owned    = pack([(i, i = 1, length)], local % is_owned)

    ! each interface unknown's place in the interface, and the part's
    ! entries that send to it and receive from it
    allocate(interface_of(nb), source=0)
    allocate(local % colour_first(maxval(colours) + 1))
    local % interface_length = 0
    do c = 1, maxval(colours)
       local % colour_first(c) = local % interface_length + 1
       do b = 1, nb
          if (.not. in_interface(b) .or. colours(b) /= c) cycle
          interface_of(b) = local % interface_length
          local % interface_length = local % interface_length + w
       end do
    end do
    local % colour_first(maxval(colours) + 1) = local % interface_length + 1
    allocate(position(length), source=0)
    do i = 1, length
       b = (local % whole_index(i) - 1) / w + 1
       if (in_interface(b)) position(i) = interface_of(b) + local % whole_index(i) - (b - 1) * w
    end do
    ! in the order of the interface, so each colour's entries are
    ! consecutive
    order = interface_order(position)
    local % sent        = pack(order, local % is_owned(order))
    local % sent_at     = position(local % sent)
    local % received    = pack(order, .not. local % is_owned(order))
    local % received_at = position(local % received)
    local % sent_first     = colour_bounds(local % sent_at, local % colour_first)
    local % received_first = colour_bounds(local % received_at, local % colour_first)

  end function part

  ! the part's entries that lie in the interface, by ascending place
  function interface_order(position) result(order)
    integer, intent(in) :: position(:)
    integer, allocatable :: order(:)
    integer, allocatable :: at(:)
    integer :: i
    allocate(at(max(0, maxval(position))), source=0)
    do i = 1, size(position)
       if (position(i) > 0) at(position(i)) = i
    end do
    order = pack(at, at > 0)
  end function interface_order

  ! for places ascending, the first of each colour's range of places
  function colour_bounds(places, colour_first) result(first)
    integer, intent(in) :: places(:), colour_first(:)
    integer, allocatable :: first(:)
    integer :: c, k
    allocate(first(size(colour_first)))
    k = 1
    do c = 1, size(colour_first)
       do while (k <= size(places))
          if (places(k) >= colour_first(c)) exit
          k = k + 1
       end do
       first(c) = k
    end do
  end function colour_bounds

  !===================================================================!
  ! THE GATHER. Of the whole, by one sum over the images of the vector
  ! with every entry not owned set to zero: each entry is owned by one
  ! image, so the sum is the owner's value everywhere. Of a part, by
  ! one sum of the interface, each image placing the interface entries
  ! it owns: every halo entry then takes its owner's value.
  !===================================================================!

  subroutine gather(this, x, colour)

    class(exchange), intent(in)    :: this
    real(dp)       , intent(inout) :: x(:)
    ! of a part, the interface of this colour alone
    integer        , intent(in), optional :: colour

    real(dp), allocatable :: shared(:)
    integer :: first, last, s1, s2, r1, r2

    if (this % images == 1) return
    if (.not. allocated(this % whole_index)) then
       where (.not. this % is_owned) x = 0.0_dp
       call co_sum(x)
       return
    end if
    first = 1
    last  = this % interface_length
    s1 = 1
    s2 = size(this % sent)
    r1 = 1
    r2 = size(this % received)
    if (present(colour)) then
       first = this % colour_first(colour)
       last  = this % colour_first(colour + 1) - 1
       s1 = this % sent_first(colour)
       s2 = this % sent_first(colour + 1) - 1
       r1 = this % received_first(colour)
       r2 = this % received_first(colour + 1) - 1
    end if
    allocate(shared(first:last), source=0.0_dp)
    shared(this % sent_at(s1:s2)) = x(this % sent(s1:s2))
    call co_sum(shared)
    x(this % received(r1:r2)) = shared(this % received_at(r1:r2))

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

  subroutine restricted_triples(this, rows, columns, weights, part_rows, part_columns, part_weights)

    class(exchange)      , intent(in)  :: this
    integer              , intent(in)  :: rows(:), columns(:)
    real(dp)             , intent(in)  :: weights(:)
    integer , allocatable, intent(out) :: part_rows(:), part_columns(:)
    real(dp), allocatable, intent(out) :: part_weights(:)

    integer, allocatable :: place_of(:)
    logical, allocatable :: owned_row(:)
    integer :: i

    allocate(place_of(size(this % owner)), source=0)
    place_of(this % whole_index) = [(i, i = 1, size(this % whole_index))]
    owned_row = this % owner(rows) == this % image
    part_rows    = place_of(pack(rows, owned_row))
    part_columns = place_of(pack(columns, owned_row))
    part_weights = pack(weights, owned_row)
    if (any(part_columns == 0)) then
       error stop 'operation_exchange: a column of an owned row lies outside the part'
    end if

  end subroutine restricted_triples

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
