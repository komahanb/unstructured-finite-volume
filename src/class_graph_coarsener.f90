!=====================================================================!
! The concrete graph coarsener.
!
! C makes a smaller graph out of a bigger one by gluing cells into
! blocks. Same shape, less detail:
!
!      o o o o                 O   O
!      o o o o     ------>
!      o o o o                 O   O
!
!      twelve cells            four blocks
!
! Two fine cells joined by a face put a face between their blocks,
! unless they landed in the same block - then that face disappears
! inside it. That is the whole rule.
!
! Coarsening and partitioning are both transforms and they are not the
! same family. Partitioning changes WHO HOLDS WHAT: whole to parts,
! same detail. Coarsening changes the detail: fine to coarse, same
! whole.
!
! What it is for: a multigrid level, where the slow smooth part of an
! error is cheap to damp; a first guess to start a fine solve from;
! and a quick look at a mesh too big to draw.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_coarsener

  use iso_fortran_env     , only : dp => REAL64
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use graph_fractal      , only : set_graph => graph
  use transform_structure, only : graph_transform
  use class_graph         , only : directed_stored_graph
  use class_graph_field   , only : field

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
  ! Walk the cells and glue each unclaimed one to a neighbour that is
  ! also unclaimed. Cheap, deterministic, and it follows the
  ! connections rather than the numbering.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_PAIRWISE = 1

  !-------------------------------------------------------------------!
  ! Take a map computed elsewhere: which block each cell goes to.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_ADOPTED = 2

  !===================================================================!
  ! One coarsener, holding how it glues and the map it ends up with.
  !===================================================================!

  type, extends(graph_transform) :: coarsener

     integer :: rule = COARSEN_PAIRWISE

     !----------------------------------------------------------------!
     ! Which block each fine cell belongs to, and how many blocks there
     ! are. The map is the coarsening; everything else follows from it.
     !----------------------------------------------------------------!

     integer, allocatable :: block_of(:)
     integer              :: nblocks = 0

     !----------------------------------------------------------------!
     ! Add the fine values or average them. A residual adds, because
     ! it is a total. A state averages, because it is a level.
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
  ! Build a coarsener that follows one rule. Pairwise needs nothing
  ! more; adopted brings its own map, cell -> block. The average flag
  ! picks what a block's value is: the mean of its cells, or their
  ! sum.
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
  ! Can this coarsener say anything about that graph?
  !
  ! A graph already down to one cell is accepted: it is as coarse as
  ! a graph gets, so coarsening returns it unchanged, and returning
  ! the input unchanged is a valid result. Rejecting it would push a
  ! special case onto every caller - including the smallest level of
  ! a multigrid hierarchy, the worst place for one.
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
  ! Can this coarsener say anything about that data? Yes for a field
  ! whose entries match the graph; the graph gate answers first.
  !===================================================================!

  logical function defined_on_data(this, input_graph, input_data)

    class(coarsener) , intent(in) :: this
    class(directed_graph)     , intent(in) :: input_graph
    class(graph_field), intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (field)
       block
         type(set_graph) :: dom
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
  ! C. Work out the blocks, then draw a face between two blocks
  ! wherever a fine face ran between them.
  !===================================================================!

  subroutine coarsen_graph(this, fine_graph, coarse_graph)

    class(coarsener), intent(in)               :: this
    class(directed_graph)    , intent(in)               :: fine_graph
    class(directed_graph)    , allocatable, intent(out) :: coarse_graph

    integer, allocatable :: blk(:), tails(:), heads(:)
    integer :: nb, ne, e, t, h, bt, bh, n
    logical, allocatable :: drawn(:,:)

    call blocks_of(this, fine_graph, blk, nb)

    ne = fine_graph % num_edges()
    allocate(tails(ne), heads(ne))
    allocate(drawn(nb, nb))
    drawn = .false.
    n = 0

    do e = 1, ne

       t  = fine_graph % edge_tail(e)
       bt = blk(t)

       if (.not. fine_graph % edge_has_head(e)) then
          ! A wall stays a wall. The block behind it inherits it.
          n = n + 1
          tails(n) = bt
          heads(n) = 0
          cycle
       end if

       h  = fine_graph % edge_head(e)
       bh = blk(h)

       ! A face whose two cells landed in the same block vanishes
       ! inside it. A face between two blocks is drawn once, however
       ! many fine faces ran between them.
       if (bt == bh) cycle
       if (drawn(bt, bh)) cycle

       drawn(bt, bh) = .true.
       n = n + 1
       tails(n) = bt
       heads(n) = bh

    end do

    allocate(coarse_graph, source = &
         & directed_stored_graph(nb, tails=tails(1:n), heads=heads(1:n), &
         &              number=fine_graph % id()))

  end subroutine coarsen_graph

  !===================================================================!
  ! The aggregate map, answered publicly: which block each fine cell
  ! belongs to, and how many blocks there are. A multigrid reads
  ! this to build its coarse operator; the map is the coarsener's to
  ! own and everyone else's to consume.
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
    integer :: nv, v, i, mate

    nv = fine_graph % num_vertices()
    allocate(blk(nv))

    if (this % rule == COARSEN_ADOPTED) then
       blk = this % block_of(1:nv)
       nb  = max(maxval(blk), this % nblocks)
       return
    end if

    ! Pairwise. Walk the cells; each unclaimed one starts a block and
    ! pulls in the first unclaimed neighbour it can find. A cell with
    ! no free neighbour is a block of one.
    blk = 0
    nb  = 0

    do v = 1, nv
       if (blk(v) /= 0) cycle
       nb      = nb + 1
       blk(v)  = nb
       call fine_graph % adjacent_vertices(v, nbrs)
       mate = 0
       do i = 1, size(nbrs)
          if (blk(nbrs(i)) == 0) then
             mate = nbrs(i)
             exit
          end if
       end do
       if (mate /= 0) blk(mate) = nb
    end do

  end subroutine blocks_of

  !===================================================================!
  ! Lift the values onto the blocks. Several fine cells land on one
  ! coarse cell, so this has to say how they merge - added when the
  ! value is a total, averaged when it is a level.
  !
  ! Getting that choice wrong is quiet: a residual that gets averaged
  ! instead of added comes out too small by the block size, and a
  ! multigrid cycle simply converges more slowly than it should.
  !===================================================================!

  subroutine coarsen_data(this, fine_graph, fine_data, coarse_graph, coarse_data)

    class(coarsener) , intent(in)               :: this
    class(directed_graph)     , intent(in)               :: fine_graph
    class(graph_field), intent(in)               :: fine_data
    class(directed_graph)     , intent(in)               :: coarse_graph
    class(graph_field), allocatable, intent(out) :: coarse_data

    type(field)    :: out
    integer , allocatable :: blk(:), tally(:)
    real(dp), allocatable :: fv(:), cv(:)
    integer :: nb, nv, ncomp, v, c, b

    select type (fine_data)
    class is (field)

       call blocks_of(this, fine_graph, blk, nb)

       nv    = fine_graph % num_vertices()
       ncomp = fine_data % num_components()

       out = field(fine_data % name(), coarse_graph % vertex_set(), coarse_graph % num_vertices(), &
            &             ncomp=ncomp, unit_name=fine_data % units())

       call fine_data % get_real_vector(fv)
       allocate(cv(nb * ncomp), tally(nb))
       cv    = 0.0_dp
       tally = 0

       do v = 1, nv
          b = blk(v)
          tally(b) = tally(b) + 1
          do c = 1, ncomp
             associate (to => (b - 1) * ncomp + c, from => (v - 1) * ncomp + c)
               if (from <= size(fv)) cv(to) = cv(to) + fv(from)
             end associate
          end do
       end do

       if (this % average) then
          do b = 1, nb
             if (tally(b) > 0) then
                do c = 1, ncomp
                   cv((b - 1) * ncomp + c) = cv((b - 1) * ncomp + c) / real(tally(b), dp)
                end do
             end if
          end do
       end if

       call out % set_real_vector(cv)
       allocate(coarse_data, source=out)

    class default
       error stop 'coarsen: this data does not ride on this transform'
    end select

  end subroutine coarsen_data

end module class_graph_coarsener
