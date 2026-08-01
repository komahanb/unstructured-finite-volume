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
! same family. Partitioning changes the frame: whole to parts, same
! detail. Coarsening changes the detail: fine to coarse, same whole.
!
! What it is for: a multigrid level, where the slow smooth part of an
! error is cheap to kill; a first guess to start a fine solve from;
! and a quick look at a mesh too big to draw.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_coarsener

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_coarsener, graph, graph_data, graph_field
  use class_graph         , only : stored_graph
  use class_graph_support , only : vertex_support
  use class_graph_field   , only : vertex_field

  implicit none

  private
  public :: coarsener
  public :: COARSEN_PAIRWISE, COARSEN_ADOPTED

  !-------------------------------------------------------------------!
  ! Walk the cells and glue each unclaimed one to a neighbour that is
  ! also unclaimed. Cheap, deterministic, and it follows the
  ! connections rather than the numbering.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_PAIRWISE = 1

  !-------------------------------------------------------------------!
  ! Take a map somebody else worked out: which block each cell goes to.
  !-------------------------------------------------------------------!

  integer, parameter :: COARSEN_ADOPTED = 2

  !===================================================================!
  ! One coarsener, carrying how it glues and the map it ends up with.
  !===================================================================!

  type, extends(graph_coarsener) :: coarsener

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

     procedure :: defined_on_graph => c_defined_on_graph
     procedure :: defined_on_data  => c_defined_on_data
     procedure :: coarsen_graph    => c_coarsen_graph
     procedure :: coarsen_data     => c_coarsen_data

  end type coarsener

  interface coarsener
     module procedure create
  end interface coarsener

contains

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
  ! A graph already down to one cell is not refused. It is as coarse
  ! as a graph gets, so coarsening it hands the same graph back, and
  ! doing nothing is a perfectly good answer. Refusing would push the
  ! decision onto every caller, who would then have to special-case
  ! the smallest level of a multigrid hierarchy - which is exactly the
  ! place nobody wants a special case.
  !
  ! What is refused is a coarsening that was handed a map and cannot
  ! read it, because then it genuinely does not know what to do.
  !===================================================================!

  pure logical function c_defined_on_graph(this, input_graph)

    class(coarsener), intent(in) :: this
    class(graph)    , intent(in) :: input_graph

    c_defined_on_graph = input_graph % num_vertices() > 0

    if (this % rule == COARSEN_ADOPTED) then
       if (.not. allocated(this % block_of)) then
          c_defined_on_graph = .false.
       else if (size(this % block_of) < input_graph % num_vertices()) then
          c_defined_on_graph = .false.
       end if
    end if

  end function c_defined_on_graph

  pure logical function c_defined_on_data(this, input_graph, input_data)

    class(coarsener) , intent(in) :: this
    class(graph)     , intent(in) :: input_graph
    class(graph_data), intent(in) :: input_data

    c_defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (vertex_field)
       c_defined_on_data = c_defined_on_data .and. input_data % num_entries() >= 0
    class default
       c_defined_on_data = .false.
    end select

  end function c_defined_on_data

  !===================================================================!
  ! C. Work out the blocks, then draw a face between two blocks
  ! wherever a fine face ran between them.
  !===================================================================!

  subroutine c_coarsen_graph(this, fine_graph, coarse_graph)

    class(coarsener), intent(in)               :: this
    class(graph)    , intent(in)               :: fine_graph
    class(graph)    , allocatable, intent(out) :: coarse_graph

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
         & stored_graph(nb, tails=tails(1:n), heads=heads(1:n), &
         &              number=fine_graph % id()))

  end subroutine c_coarsen_graph

  !===================================================================!
  ! Which block each fine cell belongs to.
  !===================================================================!

  subroutine blocks_of(this, fine_graph, blk, nb)

    class(coarsener)    , intent(in)  :: this
    class(graph)        , intent(in)  :: fine_graph
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

  subroutine c_coarsen_data(this, fine_graph, fine_data, coarse_graph, coarse_data)

    class(coarsener) , intent(in)               :: this
    class(graph)     , intent(in)               :: fine_graph
    class(graph_data), intent(in)               :: fine_data
    class(graph)     , intent(in)               :: coarse_graph
    class(graph_data), allocatable, intent(out) :: coarse_data

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    integer , allocatable :: blk(:), indices(:), tally(:)
    real(dp), allocatable :: fv(:), cv(:)
    integer :: nb, nv, ncomp, v, c, b

    select type (fine_data)
    class is (vertex_field)

       call blocks_of(this, fine_graph, blk, nb)

       nv    = fine_graph % num_vertices()
       ncomp = fine_data % num_components()

       allocate(indices(nb))
       do b = 1, nb
          indices(b) = b
       end do

       on  = vertex_support(indices)
       out = vertex_field(fine_data % name(), on, ncomp=ncomp, &
            &             unit_name=fine_data % units())

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

    end select

  end subroutine c_coarsen_data

end module class_graph_coarsener
