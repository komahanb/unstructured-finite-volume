!=====================================================================!
! Concrete graph partitioners.
!
! P cuts a graph into parts. One concrete type holds the rule, so a
! caller can hold partitioners in a plain array and a new rule costs a
! case rather than a class.
!
!         o---o---o---o---o---o
!                     :                cut where few edges cross, so
!         o---o---o   :   o---o---o    the parts have little to say
!                part 1     part 2     to each other
!
! What comes out is one part, in its own numbering, and it is still a
! graph. It also records how it relates to the whole - which cells
! it owns, which it only borrows, and what each of its own numbers
! was called in the global graph.
!
!         global graph    1   2   3   4   5   6   7   8
!                                   |   |   |
!         part 2                    1   2   3
!
!                       global_vertex_index(2) = 4
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
! out which cells go where and holds the data across the same cut.
! It does not evaluate anything.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_partitioner

  use iso_fortran_env     , only : dp => REAL64
  use graph_grammar       , only : ordinary_graph, graph_field
  use graph_set       , only : set, index_set, subset
  use graph_calculus      , only : graph_partitioner
  use class_graph         , only : stored_graph
  use class_graph_field   , only : field
  use class_graph_walk    , only : walk, WALK_VISIT_ORDER

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

     procedure :: defined_on_graph
     procedure :: defined_on_data
     procedure :: partition_graph
     procedure :: partition_data

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

  pure logical function defined_on_graph(this, input_graph)

    class(partitioner), intent(in) :: this
    class(ordinary_graph)      , intent(in) :: input_graph

    defined_on_graph = input_graph % num_vertices() > 0 .and. this % nparts >= 1

    if (this % rule == PARTITION_ADOPTED) then
       if (.not. allocated(this % adopted)) then
          defined_on_graph = .false.
       else if (size(this % adopted) < input_graph % num_vertices()) then
          defined_on_graph = .false.
       end if
    end if

  end function defined_on_graph

  !===================================================================!
  ! Data can be carried across when the graph can be cut and the data
  ! sits on that graph.
  !===================================================================!

  logical function defined_on_data(this, input_graph, input_data)

    class(partitioner), intent(in) :: this
    class(ordinary_graph)      , intent(in) :: input_graph
    class(graph_field) , intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (graph_field)
       defined_on_data = defined_on_data .and. input_data % num_entries() >= 0
    class default
       defined_on_data = .false.
    end select

  end function defined_on_data

  !===================================================================!
  ! P. Work out who owns what, gather this part's cells, and rebuild
  ! the piece as a graph in its own numbering.
  !===================================================================!

  subroutine partition_graph(this, global_graph, part_graph)

    class(partitioner), intent(in)              :: this
    class(ordinary_graph)      , intent(in)              :: global_graph
    class(ordinary_graph)      , allocatable, intent(out) :: part_graph

    integer, allocatable :: owner(:), mine(:), whereis(:)
    integer, allocatable :: ltail(:), lhead(:), eglobal(:), eowner(:), vowner(:)
    integer :: nv, ne, e, t, h, k, nkeep

    nv = global_graph % num_vertices()
    ne = global_graph % num_edges()

    call assign_owners(this, global_graph, owner)
    call gather_part(global_graph, owner, this % part, mine, whereis)

    ! Keep an edge when both its ends are in this part and at least one
    ! of them is owned here. An edge with neither end owned belongs
    ! entirely to another part; keeping it here would add its flux to
    ! the balance twice.
    allocate(ltail(ne), lhead(ne), eglobal(ne), eowner(ne))
    nkeep = 0
    do e = 1, ne
       t = global_graph % edge_tail(e)
       h = global_graph % edge_head(e)

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
       eglobal(nkeep) = e

       ! An edge is owned by the part that owns its TAIL - always,
       ! whether this part holds that tail or borrows it. One global
       ! edge therefore has exactly one owner across all parts, which
       ! is what makes assembly reconstruct a global edge field
       ! exactly once. (The branch below is vestigial: both arms
       ! assign the same thing. An earlier design let the head's
       ! owner answer for a borrowed tail; that rule was never
       ! implemented, and the uniqueness law does not need it.)
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

    ! The piece is born knowing how it relates to the whole. Without
    ! that record the assembler cannot restore whole-graph order and
    ! defined_on_graph answers .false. rather than assuming a map -
    ! and the record goes in through the door, because a graph told
    ! its frame after birth would answer one question two ways in
    ! one lifetime.
    allocate(part_graph, source = &
         & stored_graph(size(mine), tails=ltail(1:nkeep), heads=lhead(1:nkeep), &
         &              number  = this % part,   &
         &              nparts  = this % nparts, &
         &              vglobal = mine,          &
         &              vowner  = vowner,        &
         &              eglobal = eglobal(1:nkeep), &
         &              eowner  = eowner(1:nkeep)))

  end subroutine partition_graph

  !===================================================================!
  ! Decide which part owns each cell of the whole graph.
  !===================================================================!

  subroutine assign_owners(this, global_graph, owner)

    class(partitioner)  , intent(in)  :: this
    class(ordinary_graph)        , intent(in)  :: global_graph
    integer, allocatable, intent(out) :: owner(:)

    integer :: nv

    nv = global_graph % num_vertices()
    allocate(owner(nv))

    select case (this % rule)

    case (PARTITION_ADOPTED)
       owner = this % adopted(1:nv)

    case (PARTITION_BREADTH_FIRST)
       call assign_owners_breadth_first(global_graph, this % nparts, owner)

    case default
       call assign_owners_linear(nv, this % nparts, owner)

    end select

  end subroutine assign_owners

  !===================================================================!
  ! Equal blocks of the numbering. The first few parts take one extra
  ! cell when the count does not divide evenly.
  !===================================================================!

  pure subroutine assign_owners_linear(nv, nparts, owner)

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

  end subroutine assign_owners_linear

  !===================================================================!
  ! Grow every part outward from a seed, one ring at a time, so each
  ! part comes out connected and few edges are left crossing.
  !===================================================================!

  subroutine assign_owners_breadth_first(global_graph, nparts, owner)

    class(ordinary_graph), intent(in)    :: global_graph
    integer     , intent(in)    :: nparts
    integer     , intent(inout) :: owner(:)

    type(stored_graph) :: untaken
    type(walk)         :: visit
    class(graph_field), allocatable :: reached
    integer, allocatable :: locals(:), whereis(:), tails(:), heads(:), order(:)
    integer :: nv, ne, share, k, v, e, t, h, n, m

    nv    = global_graph % num_vertices()
    ne    = global_graph % num_edges()
    owner = 0
    share = (nv + nparts - 1) / nparts

    allocate(locals(nv), whereis(nv), tails(ne), heads(ne))

    do k = 1, nparts

       ! The unclaimed remainder, as the graph it is. The walk owns
       ! breadth-first; this routine only asks and reads.
       n = 0
       whereis = 0
       do v = 1, nv
          if (owner(v) == 0) then
             n = n + 1
             locals(n)  = v
             whereis(v) = n
          end if
       end do
       if (n == 0) exit

       m = 0
       do e = 1, ne
          t = global_graph % edge_tail(e)
          if (.not. global_graph % edge_has_head(e)) cycle
          h = global_graph % edge_head(e)
          if (whereis(t) > 0 .and. whereis(h) > 0) then
             m = m + 1
             tails(m) = whereis(t)
             heads(m) = whereis(h)
          end if
       end do

       untaken = stored_graph(n, tails=tails(1:m), heads=heads(1:m))

       ! The first unclaimed cell seeds the part; the visit order
       ! says who its share of the ring is.
       visit = walk(WALK_VISIT_ORDER, seed=1)
       call visit % apply(untaken, output=reached)
       call reached % get_integer_vector(order)

       do v = 1, n
          if (order(v) >= 1 .and. order(v) <= share) owner(locals(v)) = k
       end do

    end do

    ! Anything the rings never reached goes to the last part, so every
    ! cell ends up owned exactly once.
    do v = 1, nv
       if (owner(v) == 0) owner(v) = nparts
    end do

  end subroutine assign_owners_breadth_first

  !===================================================================!
  ! Collect one part's cells: the ones it owns first, then the ones it
  ! must borrow to work out its own answers.
  !
  ! The owned-first order has a practical consequence: a part's owned
  ! values sit at the front of every vector, so the piece a solver
  ! reduces over is a contiguous slice.
  !===================================================================!

  subroutine gather_part(global_graph, owner, part, mine, whereis)

    class(ordinary_graph)        , intent(in)  :: global_graph
    integer             , intent(in)  :: owner(:)
    integer             , intent(in)  :: part
    integer, allocatable, intent(out) :: mine(:)
    integer, allocatable, intent(out) :: whereis(:)

    integer, allocatable :: nbrs(:)
    integer :: nv, v, i, n

    nv = global_graph % num_vertices()

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
       call global_graph % adjacent_vertices(v, nbrs)
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

  subroutine partition_data(this, global_graph, global_data, part_graph, part_data)

    class(partitioner), intent(in)               :: this
    class(ordinary_graph)      , intent(in)               :: global_graph
    class(graph_field) , intent(in)               :: global_data
    class(ordinary_graph)      , intent(in)               :: part_graph
    class(graph_field) , allocatable, intent(out) :: part_data

    class(set), allocatable :: dom

    associate (u1 => this); end associate

    select type (global_data)

    class is (field)
       call global_data % domain(dom)
       ! Classify by embedding; coverage decides the carry inside.
       if (dom % is_subobject_of(global_graph % vertex_set())) then
          call carry_field(global_data, dom, &
               & global_graph % vertex_set(), part_graph, .true., part_data)
       else if (dom % is_subobject_of(global_graph % edge_set())) then
          call carry_field(global_data, dom, &
               & global_graph % edge_set(), part_graph, .false., part_data)
       else
          error stop 'partition: this field does not live on this graph''s domains'
       end if

    class default
       error stop 'partition: this data does not ride on this transform'
    end select

  end subroutine partition_data

  !===================================================================!
  ! One carry for both families and both coverages. A FULL field -
  ! domain equals the global set - lands on the part's own
  ! set, every part member valued through the global map. A
  ! PROPER SUBSET travels as a subset: the part-local members whose
  ! global names the subset holds, each value fetched through the
  ! GLOBAL DOMAIN'S local_index - never by raw member arithmetic -
  ! and seated on a new subobject of the part's set. A new
  ! ambient means a new declared subset: identity is not preserved
  ! across transport, extension and values are.
  !===================================================================!

  subroutine carry_field(global_data, dom, global_set, part_graph, &
       &                 on_vertices, part_data)

    type(field)       , intent(in)               :: global_data
    class(set) , intent(in)               :: dom
    type(index_set) , intent(in)               :: global_set
    class(ordinary_graph)      , intent(in)               :: part_graph
    logical           , intent(in)               :: on_vertices
    class(graph_field), allocatable, intent(out) :: part_data

    type(field)           :: out
    type(index_set)     :: part_set
    type(subset)      :: sp
    real(dp), allocatable :: fv(:), lv(:)
    integer , allocatable :: kept(:)
    integer :: nlocal, ncomp, l, c, g, n, at

    if (on_vertices) then
       nlocal       = part_graph % num_vertices()
       part_set = part_graph % vertex_set()
    else
       nlocal       = part_graph % num_edges()
       part_set = part_graph % edge_set()
    end if
    ncomp = global_data % num_components()

    call global_data % get_real_vector(fv)

    if (dom % equals(global_set)) then

       ! Full coverage: the part field lives on the part's set.
       allocate(lv(nlocal * ncomp))
       lv = 0.0_dp
       do l = 1, nlocal
          g = global_of(part_graph, l, on_vertices)
          at = dom % local_index(g)
          if (at >= 1) then
             do c = 1, ncomp
                lv((l - 1) * ncomp + c) = fv((at - 1) * ncomp + c)
             end do
          end if
       end do
       out = field(global_data % name(), part_set, ncomp=ncomp, &
            &      unit_name=global_data % units())
       call out % set_real_vector(lv)

    else

       ! Proper subset: gather the part members the subset names.
       allocate(kept(nlocal))
       n = 0
       do l = 1, nlocal
          g = global_of(part_graph, l, on_vertices)
          if (dom % has(g)) then
             n = n + 1
             kept(n) = l
          end if
       end do
       sp = subset(dom % name(), part_set, kept(1:n))

       allocate(lv(n * ncomp))
       do l = 1, n
          g  = global_of(part_graph, kept(l), on_vertices)
          at = dom % local_index(g)
          do c = 1, ncomp
             lv((l - 1) * ncomp + c) = fv((at - 1) * ncomp + c)
          end do
       end do
       out = field(global_data % name(), sp, ncomp=ncomp, &
            &      unit_name=global_data % units())
       call out % set_real_vector(lv)

    end if

    allocate(part_data, source=out)

  end subroutine carry_field

  !===================================================================!
  ! What the part's l-th member was called in the whole.
  !===================================================================!

  pure integer function global_of(part_graph, l, on_vertices)

    class(ordinary_graph), intent(in) :: part_graph
    integer     , intent(in) :: l
    logical     , intent(in) :: on_vertices

    if (on_vertices) then
       global_of = part_graph % global_vertex_index(l)
    else
       global_of = part_graph % global_edge_index(l)
    end if

  end function global_of

end module class_graph_partitioner
