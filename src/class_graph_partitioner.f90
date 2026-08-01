!=====================================================================!
! Concrete graph partitioners.
!
! P cuts a graph into parts. One concrete type carries the rule, so a
! caller can hold partitioners in a plain array and a new rule costs a
! case rather than a class.
!
!         o---o---o---o---o---o
!                     :                cut where few edges cross, so
!         o---o---o   :   o---o---o    the parts have little to say
!                part 1     part 2     to each other
!
! What comes out is one part, in its own numbering, and it is still a
! graph. It also remembers how it relates to the whole - which cells
! it owns, which it only borrows, and what each of its own numbers
! was called in the full graph.
!
!         full graph    1   2   3   4   5   6   7   8
!                                   |   |   |
!         part 2                    1   2   3
!
!                       full_vertex_index(2) = 4
!
!=====================================================================!
!
!                      OWNED, BORROWED, OVERLAP
!
! A part owns the cells it must produce answers for. It borrows the
! neighbouring cells that other parts own, because a face term needs
! the value on both sides. Together those are the overlap - what
! this part must be able to see to finish what it owns.
!
!            part 1                        part 2
!       +---------------+            +---------------+
!       |  o    o    o  |            |  o    o    o  |
!       |  o    o    o--|------------|--b    o    o  |
!       +---------------+            +---------------+
!                    \______________/
!                       part 1 borrows one cell from part 2
!
! Every cell is owned by exactly one part. That is what stops a
! conserved quantity being counted twice when the parts are added back
! together.
!
!=====================================================================!
!
!                        WHAT IS NOT HERE
!
! No physics, no geometry, no solver behaviour. A partitioner works
! out which cells go where and carries the data across the same cut.
! It does not evaluate anything.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_partitioner

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_partitioner, graph, graph_data
  use abstract_graph_types, only : graph_field, graph_support
  use class_graph         , only : stored_graph
  use class_graph_support , only : vertex_support, edge_support
  use class_graph_field   , only : vertex_field, edge_field

  implicit none

  private
  public :: partitioner
  public :: PARTITION_LINEAR, PARTITION_BREADTH_FIRST, PARTITION_ADOPTED

  !-------------------------------------------------------------------!
  ! Slice the vertex numbering into equal blocks. Cheap, deterministic,
  ! and blind to how the cells are actually joined - useful as a
  ! reference the other rules are judged against.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_LINEAR = 1

  !-------------------------------------------------------------------!
  ! Grow each part outward from a seed cell, one ring at a time, until
  ! it has its share. Follows the connections, so it cuts fewer edges.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_BREADTH_FIRST = 2

  !-------------------------------------------------------------------!
  ! Take a map computed elsewhere - a mesh file, an outside library,
  ! a previous run.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_ADOPTED = 3

  !===================================================================!
  ! One partitioner: how to cut, into how many, and which part to hand
  ! back.
  !===================================================================!

  type, extends(graph_partitioner) :: partitioner

     integer :: rule   = PARTITION_LINEAR
     integer :: nparts = 1
     integer :: part   = 1

     !----------------------------------------------------------------!
     ! The map an adopted partition was handed, one owning part per
     ! whole-graph vertex.
     !----------------------------------------------------------------!

     integer, allocatable :: adopted(:)

   contains

     procedure :: defined_on_graph => p_defined_on_graph
     procedure :: defined_on_data  => p_defined_on_data
     procedure :: partition_graph  => p_partition_graph
     procedure :: partition_data   => p_partition_data

  end type partitioner

  interface partitioner
     module procedure create
  end interface partitioner

contains

  !===================================================================!
  ! Build a partitioner. Say how to cut, into how many pieces, and
  ! which piece you want back.
  !===================================================================!

  pure type(partitioner) function create(rule, nparts, part, adopted) result(this)

    integer, intent(in)           :: rule
    integer, intent(in)           :: nparts
    integer, intent(in), optional :: part
    integer, intent(in), optional :: adopted(:)

    this % rule   = rule
    this % nparts = nparts

    if (present(part))    this % part = part
    if (present(adopted)) allocate(this % adopted, source=adopted)

  end function create

  !===================================================================!
  ! A graph can be cut when it has cells to cut up, and when an
  ! adopted map - if that is the rule - actually covers them.
  !===================================================================!

  pure logical function p_defined_on_graph(this, input_graph)

    class(partitioner), intent(in) :: this
    class(graph)      , intent(in) :: input_graph

    p_defined_on_graph = input_graph % num_vertices() > 0 .and. this % nparts >= 1

    if (this % rule == PARTITION_ADOPTED) then
       if (.not. allocated(this % adopted)) then
          p_defined_on_graph = .false.
       else if (size(this % adopted) < input_graph % num_vertices()) then
          p_defined_on_graph = .false.
       end if
    end if

  end function p_defined_on_graph

  !===================================================================!
  ! Data can be carried across when the graph can be cut and the data
  ! sits on that graph.
  !===================================================================!

  pure logical function p_defined_on_data(this, input_graph, input_data)

    class(partitioner), intent(in) :: this
    class(graph)      , intent(in) :: input_graph
    class(graph_data) , intent(in) :: input_data

    p_defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (graph_field)
       p_defined_on_data = p_defined_on_data .and. input_data % num_entries() >= 0
    class default
       p_defined_on_data = .false.
    end select

  end function p_defined_on_data

  !===================================================================!
  ! P. Work out who owns what, gather this part's cells, and rebuild
  ! the piece as a graph in its own numbering.
  !===================================================================!

  subroutine p_partition_graph(this, full_graph, part_graph)

    class(partitioner), intent(in)              :: this
    class(graph)      , intent(in)              :: full_graph
    class(graph)      , allocatable, intent(out) :: part_graph

    integer, allocatable :: owner(:), mine(:), whereis(:)
    integer, allocatable :: ltail(:), lhead(:), efull(:), eowner(:), vowner(:)
    integer :: nv, ne, e, t, h, k, nkeep

    nv = full_graph % num_vertices()
    ne = full_graph % num_edges()

    call stamp_owners(this, full_graph, owner)
    call gather_part(full_graph, owner, this % part, mine, whereis)

    ! Keep an edge when both its ends are in this part and at least one
    ! of them is owned here. An edge with neither end owned belongs
    ! entirely to another part; keeping it here would add its flux to
    ! the balance twice.
    allocate(ltail(ne), lhead(ne), efull(ne), eowner(ne))
    nkeep = 0
    do e = 1, ne
       t = full_graph % edge_tail(e)
       h = full_graph % edge_head(e)

       if (whereis(t) == 0) cycle
       if (h >= 1) then
          if (whereis(h) == 0) cycle
       end if
       if (owner(t) /= this % part) then
          if (h < 1) cycle
          if (owner(h) /= this % part) cycle
       end if

       nkeep = nkeep + 1
       ltail(nkeep) = whereis(t)
       if (h >= 1) then
          lhead(nkeep) = whereis(h)
       else
          lhead(nkeep) = 0
       end if
       efull(nkeep) = e

       ! An edge is owned by the part that owns its tail, unless the
       ! tail is borrowed, in which case the head's owner answers for it.
       if (owner(t) == this % part) then
          eowner(nkeep) = owner(t)
       else
          eowner(nkeep) = owner(t)
       end if
    end do

    allocate(vowner(size(mine)))
    do k = 1, size(mine)
       vowner(k) = owner(mine(k))
    end do

    allocate(part_graph, source = &
         & stored_graph(size(mine), tails=ltail(1:nkeep), heads=lhead(1:nkeep), &
         &              number=this % part))

    ! Stamp the relation back to the whole onto the piece. Without it
    ! the assembler cannot restore whole-graph order, and says so
    ! rather than guessing.
    select type (part_graph)
    type is (stored_graph)
       part_graph % cut    = .true.
       part_graph % nparts = this % nparts
       part_graph % me     = this % part
       allocate(part_graph % vfull , source=mine)
       allocate(part_graph % vowner, source=vowner)
       allocate(part_graph % efull , source=efull(1:nkeep))
       allocate(part_graph % eowner, source=eowner(1:nkeep))
    end select

  end subroutine p_partition_graph

  !===================================================================!
  ! Decide which part owns each cell of the whole graph.
  !===================================================================!

  subroutine stamp_owners(this, full_graph, owner)

    class(partitioner)  , intent(in)  :: this
    class(graph)        , intent(in)  :: full_graph
    integer, allocatable, intent(out) :: owner(:)

    integer :: nv

    nv = full_graph % num_vertices()
    allocate(owner(nv))

    select case (this % rule)

    case (PARTITION_ADOPTED)
       owner = this % adopted(1:nv)

    case (PARTITION_BREADTH_FIRST)
       call stamp_breadth_first(full_graph, this % nparts, owner)

    case default
       call stamp_linear(nv, this % nparts, owner)

    end select

  end subroutine stamp_owners

  !===================================================================!
  ! Equal blocks of the numbering. The first few parts take one extra
  ! cell when the count does not divide evenly.
  !===================================================================!

  pure subroutine stamp_linear(nv, nparts, owner)

    integer, intent(in)    :: nv, nparts
    integer, intent(inout) :: owner(:)

    integer :: v, base, extra, lo, hi, k

    base  = nv / nparts
    extra = mod(nv, nparts)

    hi = 0
    do k = 1, nparts
       lo = hi + 1
       hi = lo + base - 1
       if (k <= extra) hi = hi + 1
       do v = lo, min(hi, nv)
          owner(v) = k
       end do
    end do

  end subroutine stamp_linear

  !===================================================================!
  ! Grow every part outward from a seed, one ring at a time, so each
  ! part comes out connected and few edges are left crossing.
  !===================================================================!

  subroutine stamp_breadth_first(full_graph, nparts, owner)

    class(graph), intent(in)    :: full_graph
    integer     , intent(in)    :: nparts
    integer     , intent(inout) :: owner(:)

    integer, allocatable :: queue(:), nbrs(:)
    integer :: nv, share, k, v, seed, head_of_queue, tail_of_queue, taken, i

    nv    = full_graph % num_vertices()
    owner = 0
    share = (nv + nparts - 1) / nparts

    allocate(queue(nv))

    do k = 1, nparts

       ! Start from the first cell nobody has claimed.
       seed = 0
       do v = 1, nv
          if (owner(v) == 0) then
             seed = v
             exit
          end if
       end do
       if (seed == 0) exit

       owner(seed)   = k
       queue(1)      = seed
       head_of_queue = 1
       tail_of_queue = 1
       taken         = 1

       ! Take neighbours, then their neighbours, until this part has
       ! its share or the ring runs out.
       do while (head_of_queue <= tail_of_queue .and. taken < share)
          v = queue(head_of_queue)
          head_of_queue = head_of_queue + 1
          call full_graph % adjacent_vertices(v, nbrs)
          do i = 1, size(nbrs)
             if (taken >= share) exit
             if (owner(nbrs(i)) == 0) then
                owner(nbrs(i))       = k
                taken                = taken + 1
                tail_of_queue        = tail_of_queue + 1
                queue(tail_of_queue) = nbrs(i)
             end if
          end do
       end do

    end do

    ! Anything the rings never reached goes to the last part, so every
    ! cell ends up owned exactly once.
    do v = 1, nv
       if (owner(v) == 0) owner(v) = nparts
    end do

  end subroutine stamp_breadth_first

  !===================================================================!
  ! Collect one part's cells: the ones it owns first, then the ones it
  ! must borrow to work out its own answers.
  !
  ! The owned-first order has a practical consequence: a part's owned
  ! values sit at the front of every vector, so the piece a solver
  ! reduces over is a contiguous slice.
  !===================================================================!

  subroutine gather_part(full_graph, owner, part, mine, whereis)

    class(graph)        , intent(in)  :: full_graph
    integer             , intent(in)  :: owner(:)
    integer             , intent(in)  :: part
    integer, allocatable, intent(out) :: mine(:)
    integer, allocatable, intent(out) :: whereis(:)

    integer, allocatable :: nbrs(:)
    integer :: nv, v, i, n

    nv = full_graph % num_vertices()

    allocate(whereis(nv))
    whereis = 0
    allocate(mine(nv))
    n = 0

    do v = 1, nv
       if (owner(v) == part) then
          n = n + 1
          mine(n)    = v
          whereis(v) = n
       end if
    end do

    do v = 1, nv
       if (owner(v) /= part) cycle
       call full_graph % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (owner(nbrs(i)) /= part .and. whereis(nbrs(i)) == 0) then
             n = n + 1
             mine(n)          = nbrs(i)
             whereis(nbrs(i)) = n
          end if
       end do
    end do

    mine = mine(1:n)

  end subroutine gather_part

  !===================================================================!
  ! Carry the data across the very same cut, by the map the part graph
  ! already holds. Nothing here recomputes the cut, so the values
  ! cannot drift out of step with the structure.
  !===================================================================!

  subroutine p_partition_data(this, full_graph, full_data, part_graph, part_data)

    class(partitioner), intent(in)               :: this
    class(graph)      , intent(in)               :: full_graph
    class(graph_data) , intent(in)               :: full_data
    class(graph)      , intent(in)               :: part_graph
    class(graph_data) , allocatable, intent(out) :: part_data

    select type (full_data)

    class is (vertex_field)
       call carry_vertex_field(full_data, part_graph, part_data)

    class is (edge_field)
       call carry_edge_field(full_data, part_graph, part_data)

    end select

  end subroutine p_partition_data

  !===================================================================!
  ! A cell field follows the vertex map: the part's cell l was called
  ! full_vertex_index(l) in the whole graph, so its values come from
  ! there.
  !===================================================================!

  subroutine carry_vertex_field(full_data, part_graph, part_data)

    type(vertex_field), intent(in)               :: full_data
    class(graph)      , intent(in)               :: part_graph
    class(graph_data) , allocatable, intent(out) :: part_data

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    real(dp), allocatable :: fv(:), lv(:)
    integer , allocatable :: locals(:)
    integer :: nlocal, ncomp, l, c

    nlocal = part_graph % num_vertices()
    ncomp  = full_data % num_components()

    allocate(locals(nlocal))
    do l = 1, nlocal
       locals(l) = l
    end do

    on  = vertex_support(locals)
    out = vertex_field(full_data % name(), on, ncomp=ncomp, unit_name=full_data % units())

    call full_data % get_real_vector(fv)
    allocate(lv(nlocal * ncomp))
    lv = 0.0_dp

    do l = 1, nlocal
       do c = 1, ncomp
          associate (from => (part_graph % full_vertex_index(l) - 1) * ncomp + c)
            if (from >= 1 .and. from <= size(fv)) lv((l - 1) * ncomp + c) = fv(from)
          end associate
       end do
    end do

    call out % set_real_vector(lv)
    allocate(part_data, source=out)

  end subroutine carry_vertex_field

  !===================================================================!
  ! A face field follows the edge map in the same way.
  !===================================================================!

  subroutine carry_edge_field(full_data, part_graph, part_data)

    type(edge_field), intent(in)               :: full_data
    class(graph)    , intent(in)               :: part_graph
    class(graph_data), allocatable, intent(out) :: part_data

    type(edge_field)      :: out
    type(edge_support)    :: on
    real(dp), allocatable :: fv(:), lv(:)
    integer , allocatable :: locals(:)
    integer :: nlocal, ncomp, l, c

    nlocal = part_graph % num_edges()
    ncomp  = full_data % num_components()

    allocate(locals(nlocal))
    do l = 1, nlocal
       locals(l) = l
    end do

    on  = edge_support(locals)
    out = edge_field(full_data % name(), on, ncomp=ncomp, unit_name=full_data % units())

    call full_data % get_real_vector(fv)
    allocate(lv(nlocal * ncomp))
    lv = 0.0_dp

    do l = 1, nlocal
       do c = 1, ncomp
          associate (from => (part_graph % full_edge_index(l) - 1) * ncomp + c)
            if (from >= 1 .and. from <= size(fv)) lv((l - 1) * ncomp + c) = fv(from)
          end associate
       end do
    end do

    call out % set_real_vector(lv)
    allocate(part_data, source=out)

  end subroutine carry_edge_field

end module class_graph_partitioner
