!=====================================================================!
! The concrete graph assembler.
!
! P inverse, and only that. It puts a piece back into whole-graph
! order and brings its data with it.
!
!         part 2   1   2   3
!                  |   |   |
!         whole    3   4   5           by the map the part already
!                                      carries - the assembler reads
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
! answers. Those borrowed values are somebody else's copy. Collect
! them too and a conserved quantity is counted twice - mass appears
! from nowhere, and it appears only in parallel, only near a partition
! boundary, which is about the worst place to go looking for it.
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
! The moment any of those move in here, the one-line law stops being
! readable and this becomes the old assembler again under a new name.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_assembler

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_assembler, graph, graph_data
  use abstract_graph_types, only : graph_field
  use class_graph         , only : stored_graph
  use class_graph_support , only : vertex_support, edge_support
  use class_graph_field   , only : vertex_field, edge_field

  implicit none

  private
  public :: assembler

  !===================================================================!
  ! The assembler holds nothing. Everything it needs to find the way
  ! home is already stamped on the part it is handed.
  !===================================================================!

  type, extends(graph_assembler) :: assembler

   contains

     procedure :: defined_on_graph => a_defined_on_graph
     procedure :: defined_on_data  => a_defined_on_data
     procedure :: assemble_graph   => a_assemble_graph
     procedure :: assemble_data    => a_assemble_data

  end type assembler

  interface assembler
     module procedure create
  end interface assembler

contains

  pure type(assembler) function create() result(this)

    ! Nothing to set up. Everything the assembler needs to find the
    ! way home is stamped on the part it will be handed.

  end function create

  !===================================================================!
  ! A part can go home only if it remembers where home was. A graph
  ! that came straight off a mesh file carries no relation, and this
  ! says so rather than guessing an identity map and being wrong in
  ! silence.
  !===================================================================!

  pure logical function a_defined_on_graph(this, input_graph)

    class(assembler), intent(in) :: this
    class(graph)    , intent(in) :: input_graph

    a_defined_on_graph = input_graph % has_part_relation() .or. &
         &               input_graph % num_parts() == 1

  end function a_defined_on_graph

  pure logical function a_defined_on_data(this, input_graph, input_data)

    class(assembler) , intent(in) :: this
    class(graph)     , intent(in) :: input_graph
    class(graph_data), intent(in) :: input_data

    a_defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (graph_field)
       a_defined_on_data = a_defined_on_data .and. input_data % num_entries() >= 0
    class default
       a_defined_on_data = .false.
    end select

  end function a_defined_on_data

  !===================================================================!
  ! Put the piece back in whole-graph order.
  !
  ! Every local cell is renamed to what it was called in the whole
  ! graph, and every edge with it. What comes out is a graph again,
  ! with no partition stamped on it - because a whole graph is not a
  ! part of anything.
  !===================================================================!

  subroutine a_assemble_graph(this, part_graph, full_graph)

    class(assembler), intent(in)               :: this
    class(graph)    , intent(in)               :: part_graph
    class(graph)    , allocatable, intent(out) :: full_graph

    integer, allocatable :: tails(:), heads(:)
    integer :: ne, e, nv_full, l, biggest

    ne = part_graph % num_edges()

    ! The whole graph is at least as big as the largest whole-graph
    ! name this part knows.
    biggest = 0
    do l = 1, part_graph % num_vertices()
       biggest = max(biggest, part_graph % full_vertex_index(l))
    end do
    nv_full = biggest

    allocate(tails(ne), heads(ne))
    do e = 1, ne
       tails(e) = part_graph % full_vertex_index(part_graph % edge_tail(e))
       if (part_graph % edge_has_head(e)) then
          heads(e) = part_graph % full_vertex_index(part_graph % edge_head(e))
       else
          heads(e) = 0
       end if
    end do

    allocate(full_graph, source = &
         & stored_graph(nv_full, tails=tails, heads=heads, number=part_graph % id()))

  end subroutine a_assemble_graph

  !===================================================================!
  ! Bring the data home.
  !
  ! The answer is laid out on the whole graph. Only the entries this
  ! part owns are written; everything else is left at zero, so adding
  ! the answers from every part rebuilds the whole field exactly once.
  !===================================================================!

  subroutine a_assemble_data(this, part_graph, part_data, full_graph, full_data)

    class(assembler) , intent(in)               :: this
    class(graph)     , intent(in)               :: part_graph
    class(graph_data), intent(in)               :: part_data
    class(graph)     , intent(in)               :: full_graph
    class(graph_data), allocatable, intent(out) :: full_data

    select type (part_data)

    class is (vertex_field)
       call bring_vertex_field_home(part_data, part_graph, full_graph, full_data)

    class is (edge_field)
       call bring_edge_field_home(part_data, part_graph, full_graph, full_data)

    end select

  end subroutine a_assemble_data

  !===================================================================!
  ! Cell values go back to the cells they came from - the owned ones
  ! only.
  !===================================================================!

  subroutine bring_vertex_field_home(part_data, part_graph, full_graph, full_data)

    type(vertex_field), intent(in)               :: part_data
    class(graph)      , intent(in)               :: part_graph
    class(graph)      , intent(in)               :: full_graph
    class(graph_data) , allocatable, intent(out) :: full_data

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: indices(:)
    integer :: nfull, ncomp, l, c, f, me

    nfull = full_graph % num_vertices()
    ncomp = part_data % num_components()
    me    = part_graph % id()

    allocate(indices(nfull))
    do f = 1, nfull
       indices(f) = f
    end do

    on  = vertex_support(indices)
    out = vertex_field(part_data % name(), on, ncomp=ncomp, unit_name=part_data % units())

    call part_data % get_real_vector(lv)
    allocate(fv(nfull * ncomp))
    fv = 0.0_dp

    do l = 1, part_graph % num_vertices()

       ! Borrowed cells belong to somebody else. Skip them, or the
       ! same value lands in the answer twice.
       if (part_graph % has_part_relation()) then
          if (part_graph % vertex_owner_part(l) /= me) cycle
       end if

       f = part_graph % full_vertex_index(l)
       do c = 1, ncomp
          associate (to => (f - 1) * ncomp + c, from => (l - 1) * ncomp + c)
            if (to >= 1 .and. to <= size(fv) .and. from <= size(lv)) fv(to) = lv(from)
          end associate
       end do

    end do

    call out % set_real_vector(fv)
    allocate(full_data, source=out)

  end subroutine bring_vertex_field_home

  !===================================================================!
  ! Face values, the same way, by the edge map.
  !===================================================================!

  subroutine bring_edge_field_home(part_data, part_graph, full_graph, full_data)

    type(edge_field), intent(in)                :: part_data
    class(graph)    , intent(in)                :: part_graph
    class(graph)    , intent(in)                :: full_graph
    class(graph_data), allocatable, intent(out) :: full_data

    type(edge_field)      :: out
    type(edge_support)    :: on
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: indices(:)
    integer :: nfull, ncomp, l, c, f, me

    nfull = full_graph % num_edges()
    ncomp = part_data % num_components()
    me    = part_graph % id()

    allocate(indices(nfull))
    do f = 1, nfull
       indices(f) = f
    end do

    on  = edge_support(indices)
    out = edge_field(part_data % name(), on, ncomp=ncomp, unit_name=part_data % units())

    call part_data % get_real_vector(lv)
    allocate(fv(nfull * ncomp))
    fv = 0.0_dp

    do l = 1, part_graph % num_edges()

       if (part_graph % has_part_relation()) then
          if (part_graph % edge_owner_part(l) /= me) cycle
       end if

       f = part_graph % full_edge_index(l)
       do c = 1, ncomp
          associate (to => (f - 1) * ncomp + c, from => (l - 1) * ncomp + c)
            if (to >= 1 .and. to <= size(fv) .and. from <= size(lv)) fv(to) = lv(from)
          end associate
       end do

    end do

    call out % set_real_vector(fv)
    allocate(full_data, source=out)

  end subroutine bring_edge_field_home

end module class_graph_assembler
