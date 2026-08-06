!=====================================================================!
! The concrete graph assembler.
!
! P inverse, and only that. It puts a piece back into whole-graph
! order and brings its data with it.
!
!         part 2   1   2   3
!                  |   |   |
!         whole    3   4   5           by the map the part already
!                                      holds - the assembler reads
!                                      the relation, it never invents
!                                      one
!
! The law it has to satisfy:
!
!         assemble( partition( G ) )     ==  G
!         assemble( partition( G, D ) )  ==  ( G, D )
!
!=====================================================================!
!
!                   ONLY OWNED VALUES ARE COLLECTED
!
! A part borrows the cells around its edge so it can work out its own
! answers. A borrowed value is a copy of a value another part owns.
! Collecting both copies counts a conserved quantity twice - mass
! appears from nowhere, and it appears only in parallel, only near a
! partition boundary, where such an error is hardest to locate.
!
!            part 1                        part 2
!       +---------------+            +---------------+
!       |  o    o    o  |            |  o    o    o  |
!       |  o    o    O--|------------|--b    o    o  |
!       +---------------+            +---------------+
!                       \____________/
!                    part 1 borrows this cell.
!                    part 2 owns it and answers for it.
!                    exactly one of them is collected.
!
!=====================================================================!
!
!                     WHAT ONE PART CAN AND CANNOT DO
!
! The contract hands the assembler a single part. So it can restore
! everything that part owns, and it cannot invent what it never saw.
!
! With one part, that is the whole graph and the round trip is exact.
! With several, each call fills in that part's own share and leaves
! the rest alone, so summing the answers over all the parts rebuilds
! the whole. The union of the owned sets is the whole graph and the
! sets do not overlap, which is what makes that sum right.
!
! ASSEMBLER MEANS THIS AND NOTHING ELSE. No physics, no boundary
! conditions, no residual, no matrix, no file, no solver behaviour.
! Any of those added here obscures the one-line law the type exists
! to keep visible.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_assembler

  use iso_fortran_env     , only : dp => REAL64
  use structure_graph, only : graph
  use data_graph_field, only : graph_field
  use transform_graph_transform, only : graph_assembler
  use structure_graph, only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use structure_stored_graph         , only : stored_graph
  use structure_support , only : support
  use data_field   , only : field

  implicit none

  private
  public :: assembler

  !===================================================================!
  ! The assembler holds nothing. Everything it needs - the index maps
  ! and the ownership - is already recorded on the part it is handed.
  !===================================================================!

  type, extends(graph_assembler) :: assembler

   contains

     procedure :: defined_on_graph
     procedure :: assemble_graph
     procedure :: assemble_data

  end type assembler

  interface assembler
     module procedure create
  end interface assembler

contains

  !===================================================================!
  ! Build an assembler. There is nothing to choose - everything
  ! assembly needs is recorded on the part it will be handed.
  !===================================================================!

  pure type(assembler) function create() result(this)

    type(assembler) :: blank

    this = blank

  end function create

  !===================================================================!
  ! A part can be assembled back only if it holds its relation to
  ! the whole. A graph straight off a mesh file holds no relation,
  ! and this check returns .false. rather than assuming an identity
  ! map and being wrong in silence.
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(assembler), intent(in) :: this
    class(graph)    , intent(in) :: input_graph

    associate (u1 => this); end associate

    defined_on_graph = input_graph % has_part_relation() .or. &
         &               input_graph % num_parts() == 1

  end function defined_on_graph

  !===================================================================!
  ! Put the piece back in whole-graph order.
  !
  ! Every cell of the part is renamed to what the whole graph called
  ! it, and every edge with it. What comes out is a graph again,
  ! and holds no partition record - because a whole graph is not a
  ! part of anything.
  !===================================================================!

  subroutine assemble_graph(this, part_graph, global_graph)

    class(assembler), intent(in)               :: this
    class(graph)    , intent(in)               :: part_graph
    class(graph)    , allocatable, intent(out) :: global_graph

    integer, allocatable :: tails(:), heads(:)
    integer :: ne, e, nv_global, l, biggest

    associate (u1 => this); end associate

    ne = part_graph % num_edges()

    ! The whole graph is at least as big as the largest whole-graph
    ! index recorded on this part.
    biggest = 0
    do l = 1, part_graph % num_vertices()
       biggest = max(biggest, part_graph % global_vertex_index(l))
    end do
    nv_global = biggest

    allocate(tails(ne), heads(ne))
    do e = 1, ne
       tails(e) = part_graph % global_vertex_index(part_graph % edge_tail(e))
       if (part_graph % edge_has_head(e)) then
          heads(e) = part_graph % global_vertex_index(part_graph % edge_head(e))
       else
          heads(e) = 0
       end if
    end do

    allocate(global_graph, source = &
         & stored_graph(nv_global, tails=tails, heads=heads, number=part_graph % id()))

  end subroutine assemble_graph

  !===================================================================!
  ! Assemble the data back onto the whole graph.
  !
  ! The answer is laid out on the whole graph. Only the entries this
  ! part owns are written; everything else is left at zero, so adding
  ! the answers from every part rebuilds the whole field exactly once.
  !===================================================================!

  subroutine assemble_data(this, part_graph, part_data, global_graph, global_data)

    class(assembler) , intent(in)               :: this
    class(graph)     , intent(in)               :: part_graph
    class(graph_field), intent(in)               :: part_data
    class(graph)     , intent(in)               :: global_graph
    class(graph_field), allocatable, intent(out) :: global_data

    associate (u1 => this); end associate

    select type (part_data)

    class is (field)
       if (part_data % on % side() == GRAPH_SIDE_VERTEX) then
          call assemble_vertex_field(part_data, part_graph, global_graph, global_data)
       else
          call assemble_edge_field(part_data, part_graph, global_graph, global_data)
       end if

    class default
       error stop 'assemble: this data does not ride on this transform'
    end select

  end subroutine assemble_data

  !===================================================================!
  ! Cell values go back to the cells they came from - the owned ones
  ! only.
  !===================================================================!

  subroutine assemble_vertex_field(part_data, part_graph, global_graph, global_data)

    type(field), intent(in)               :: part_data
    class(graph)      , intent(in)               :: part_graph
    class(graph)      , intent(in)               :: global_graph
    class(graph_field) , allocatable, intent(out) :: global_data

    type(field)    :: out
    type(support)  :: on
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: indices(:)
    integer :: nglobal, ncomp, l, c, f, me

    nglobal = global_graph % num_vertices()
    ncomp = part_data % num_components()
    me    = part_graph % id()

    indices = [(f, f = 1, nglobal)]

    on  = support(GRAPH_SIDE_VERTEX, indices)
    out = field(part_data % name(), on, ncomp=ncomp, unit_name=part_data % units())

    call part_data % get_real_vector(lv)
    allocate(fv(nglobal * ncomp))
    fv = 0.0_dp

    do l = 1, part_graph % num_vertices()

       ! Borrowed cells belong to another part. Skipping them keeps
       ! each value out of the answer exactly once.
       if (part_graph % has_part_relation()) then
          if (part_graph % vertex_owner_part(l) /= me) cycle
       end if

       f = part_graph % global_vertex_index(l)
       do c = 1, ncomp
          associate (to => (f - 1) * ncomp + c, from => (l - 1) * ncomp + c)
            if (to >= 1 .and. to <= size(fv) .and. from <= size(lv)) fv(to) = lv(from)
          end associate
       end do

    end do

    call out % set_real_vector(fv)
    allocate(global_data, source=out)

  end subroutine assemble_vertex_field

  !===================================================================!
  ! Face values, the same way, by the edge map.
  !===================================================================!

  subroutine assemble_edge_field(part_data, part_graph, global_graph, global_data)

    type(field), intent(in)                :: part_data
    class(graph)    , intent(in)                :: part_graph
    class(graph)    , intent(in)                :: global_graph
    class(graph_field), allocatable, intent(out) :: global_data

    type(field)      :: out
    type(support)    :: on
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: indices(:)
    integer :: nglobal, ncomp, l, c, f, me

    nglobal = global_graph % num_edges()
    ncomp = part_data % num_components()
    me    = part_graph % id()

    indices = [(f, f = 1, nglobal)]

    on  = support(GRAPH_SIDE_EDGE, indices)
    out = field(part_data % name(), on, ncomp=ncomp, unit_name=part_data % units())

    call part_data % get_real_vector(lv)
    allocate(fv(nglobal * ncomp))
    fv = 0.0_dp

    do l = 1, part_graph % num_edges()

       if (part_graph % has_part_relation()) then
          if (part_graph % edge_owner_part(l) /= me) cycle
       end if

       f = part_graph % global_edge_index(l)
       do c = 1, ncomp
          associate (to => (f - 1) * ncomp + c, from => (l - 1) * ncomp + c)
            if (to >= 1 .and. to <= size(fv) .and. from <= size(lv)) fv(to) = lv(from)
          end associate
       end do

    end do

    call out % set_real_vector(fv)
    allocate(global_data, source=out)

  end subroutine assemble_edge_field

end module transform_assembler
