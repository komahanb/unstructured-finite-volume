!=====================================================================!
! The concrete graph coarsener.
!
! C maps a larger graph to a smaller one by merging cells into
! blocks. Same structure, lower resolution:
!
!      o o o o                 O   O
!      o o o o     ------>
!      o o o o                 O   O
!
!      twelve cells            four blocks
!
! Two fine cells joined by a face produce a face between their
! blocks, unless both are in the same block - then that face is
! interior to the block and is removed. That is the complete rule.
!
! Coarsening and partitioning are both transforms and they are not the
! same family. Partitioning changes WHICH PART STORES WHICH CELLS:
! whole to parts, same resolution. Coarsening changes the resolution:
! fine to coarse, same whole.
!
! Uses: a multigrid level, where the smooth part of an error is
! damped at lower cost; an initial iterate for a fine solve; and a
! reduced view of a mesh too large to render.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_coarsener

  use util_precision  , only : dp
  use view_directed , only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use transform_structure, only : transform, through_blocks
  use view_directed_stored         , only : stored_directed_graph
  use field_stored   , only : stored_field

  implicit none

  private
  public :: coarsener
  public :: COARSEN_PAIRWISE, COARSEN_ADOPTED

  !===================================================================!
  ! COARSENER. The transform to fewer, larger cells - a
  ! multigrid level. coarsen_data states how several fine values
  ! merge onto one coarse cell: added for a residual, averaged for
  ! a state, volume-weighted for unequal cells.
  !===================================================================!

  !-------------------------------------------------------------------!
  ! Traverse the cells and merge each unassigned one with a neighbour
  ! that is also unassigned. Low cost, deterministic, and it follows
  ! the adjacency rather than the numbering.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_PAIRWISE = 1

  !-------------------------------------------------------------------!
  ! Take a map computed elsewhere: which block each cell is assigned to.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_ADOPTED = 2

  !===================================================================!
  ! One coarsener, storing its merge rule and the resulting map.
  !===================================================================!

  type, extends(transform) :: coarsener

     integer :: rule = COARSEN_PAIRWISE

     !----------------------------------------------------------------!
     ! Which block each fine cell belongs to, and how many blocks there
     ! are. The map is the coarsening; everything else follows from it.
     !----------------------------------------------------------------!

     integer, allocatable :: block_of(:)
     integer              :: nblocks = 0

     !----------------------------------------------------------------!
     ! Add the fine values or average them. A residual is summed,
     ! because it is an extensive quantity. A state is averaged,
     ! because it is an intensive quantity.
     !----------------------------------------------------------------!

     logical :: average = .true.

   contains

     procedure :: defined_on_graph
     procedure :: defined_on_data
     procedure :: coarsen_graph
     procedure :: coarsen_data
     procedure :: blocks

  end type coarsener

  interface coarsener
     module procedure create
  end interface coarsener

contains

  !===================================================================!
  ! Build a coarsener that follows one rule. Pairwise requires no
  ! further input; adopted supplies its own map, cell -> block. The
  ! average flag selects a block's value: the mean of its cells, or
  ! their sum.
  !===================================================================!

  pure type(coarsener) function create(rule, block_of, nblocks, average) result(this)

    integer, intent(in)           :: rule
    integer, intent(in), optional :: block_of(:)
    integer, intent(in), optional :: nblocks
    logical, intent(in), optional :: average

    this % rule = rule

    if (present(block_of)) allocate(this % block_of, source=block_of)
    if (present(nblocks))  this % nblocks = nblocks
    if (present(average))  this % average = average

  end function create

  !===================================================================!
  ! Whether this coarsener is defined on the graph.
  !
  ! A graph with one cell is accepted: no coarser graph exists, so
  ! coarsening returns it unchanged, and returning the input
  ! unchanged is a valid result. Rejecting it would require a special
  ! case in every caller - including the smallest level of a
  ! multigrid hierarchy, where a special case is least acceptable.
  !
  ! What is rejected is an adopted map that does not cover the graph,
  ! because no valid block assignment exists then.
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(coarsener), intent(in) :: this
    class(directed_graph)    , intent(in) :: input_graph

    defined_on_graph = input_graph % num_vertices() > 0

    if (this % rule == COARSEN_ADOPTED) then
       if (.not. allocated(this % block_of)) then
          defined_on_graph = .false.
       else if (size(this % block_of) < input_graph % num_vertices()) then
          defined_on_graph = .false.
       end if
    end if

  end function defined_on_graph

  !===================================================================!
  ! Whether this coarsener is defined on the data: true for a field
  ! whose entries match the graph; the graph check is evaluated first.
  !===================================================================!

  logical function defined_on_data(this, input_graph, input_data)

    class(coarsener) , intent(in) :: this
    class(directed_graph)     , intent(in) :: input_graph
    class(field), intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (stored_field)
       block
         type(graph) :: dom
         integer         :: n_dom
         dom   = input_data % domain()
         n_dom = input_data % num_entries()
         ! Full coverage, not merely family: this kernel indexes
         ! every vertex densely (AGENTS.md 5B: routing is not
         ! admissibility).
         defined_on_data = defined_on_data &
              & .and. dom % same_as(input_graph % vertex_set()) &
              & .and. input_data % num_entries() >= 0
       end block
    class default
       defined_on_data = .false.
    end select

  end function defined_on_data

  !===================================================================!
  ! C. Compute the blocks, then create a face between two blocks
  ! wherever a fine face joined them.
  !===================================================================!

  subroutine coarsen_graph(this, fine_graph, coarse_graph)

    class(coarsener), intent(in)               :: this
    class(directed_graph)    , intent(in)               :: fine_graph
    class(directed_graph)    , allocatable, intent(out) :: coarse_graph

    integer, allocatable :: blk(:), tails(:), heads(:)
    integer :: nb, ne, e, t, h, bt, bh, n
    logical, allocatable :: adjacent_blocks(:,:)

    call blocks_of(this, fine_graph, blk, nb)

    ne = fine_graph % num_edges()
    allocate(tails(ne), heads(ne))
    allocate(adjacent_blocks(nb, nb))
    adjacent_blocks = .false.
    n = 0

    do e = 1, ne

       t  = fine_graph % edge_tail(e)
       bt = blk(t)

       if (.not. fine_graph % edge_has_head(e)) then
          ! A boundary face remains a boundary face. The block
          ! containing its tail cell inherits it.
          n = n + 1
          tails(n) = bt
          heads(n) = 0
          cycle
       end if

       h  = fine_graph % edge_head(e)
       bh = blk(h)

       ! A face whose two cells are in the same block is removed. A
       ! face between two blocks is created once, whatever the number
       ! of fine faces joining them.
       if (bt == bh) cycle
       if (adjacent_blocks(bt, bh)) cycle

       adjacent_blocks(bt, bh) = .true.
       n = n + 1
       tails(n) = bt
       heads(n) = bh

    end do

    allocate(coarse_graph, source = &
         & stored_directed_graph(nb, tails=tails(1:n), heads=heads(1:n), &
         &              number=fine_graph % id()))

  end subroutine coarsen_graph

  !===================================================================!
  ! The aggregate map, returned publicly: which block each fine cell
  ! belongs to, and how many blocks there are. A multigrid reads
  ! this to build its coarse operator; the coarsener owns the map and
  ! every other caller reads it.
  !===================================================================!

  subroutine blocks(this, fine_graph, assignment, nblocks)

    class(coarsener)    , intent(in)  :: this
    class(directed_graph)        , intent(in)  :: fine_graph
    integer, allocatable, intent(out) :: assignment(:)
    integer             , intent(out) :: nblocks

    call blocks_of(this, fine_graph, assignment, nblocks)

  end subroutine blocks

  !===================================================================!
  ! Which block each fine cell belongs to.
  !===================================================================!

  subroutine blocks_of(this, fine_graph, blk, nb)

    class(coarsener)    , intent(in)  :: this
    class(directed_graph)        , intent(in)  :: fine_graph
    integer, allocatable, intent(out) :: blk(:)
    integer             , intent(out) :: nb

    integer, allocatable :: nbrs(:)
    integer :: nv, v, i, paired_vertex

    nv = fine_graph % num_vertices()
    allocate(blk(nv))

    if (this % rule == COARSEN_ADOPTED) then
       blk = this % block_of(1:nv)
       nb  = max(maxval(blk), this % nblocks)
       return
    end if

    ! Pairwise. Traverse the cells; each unassigned one starts a block
    ! and adds the first unassigned neighbour found. A cell with no
    ! unassigned neighbour is a block of one.
    blk = 0
    nb  = 0

    do v = 1, nv
       if (blk(v) /= 0) cycle
       nb      = nb + 1
       blk(v)  = nb
       call fine_graph % adjacent_vertices(v, nbrs)
       paired_vertex = 0
       do i = 1, size(nbrs)
          if (blk(nbrs(i)) == 0) then
             paired_vertex = nbrs(i)
             exit
          end if
       end do
       if (paired_vertex /= 0) blk(paired_vertex) = nb
    end do

  end subroutine blocks_of

  !===================================================================!
  ! Restrict the values onto the blocks. Several fine cells map to
  ! one coarse cell, so the rule specifies how they merge - summed
  ! when the value is extensive, averaged when it is intensive.
  !
  ! A wrong choice raises no error: a residual that is averaged
  ! instead of summed is too small by the block size, and a
  ! multigrid cycle converges more slowly than the correct one.
  !===================================================================!

  subroutine coarsen_data(this, fine_graph, fine_data, coarse_graph, coarse_data)

    class(coarsener) , intent(in)               :: this
    class(directed_graph)     , intent(in)               :: fine_graph
    class(field), intent(in)               :: fine_data
    class(directed_graph)     , intent(in)               :: coarse_graph
    class(field), allocatable, intent(out) :: coarse_data

    type(stored_field)    :: out
    integer , allocatable :: blk(:)
    real(dp), allocatable :: fv(:), cv(:)
    integer :: nb, nv, num_components

    select type (fine_data)
    class is (stored_field)

       call blocks_of(this, fine_graph, blk, nb)

       nv    = fine_graph % num_vertices()
       num_components = fine_data % num_components()

       out = stored_field(fine_data % name(), coarse_graph % vertex_set(), coarse_graph % num_vertices(), &
            &             num_components=num_components, unit_name=fine_data % units())

       call fine_data % real_vector(fv)
       call through_blocks(blk(1:nv), nb, num_components, fv, cv, transposed=.false., &
            & average=this % average)

       call out % set_real_vector(cv)
       allocate(coarse_data, source=out)

    class default
       error stop 'coarsen: the transform is defined on stored_field data only'
    end select

  end subroutine coarsen_data

end module transform_coarsener
