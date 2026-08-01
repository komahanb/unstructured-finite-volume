!=====================================================================!
! The concrete graph refiner.
!
! R is the other way from C. One cell becomes several:
!
!      O   O                   o o o o
!               ------>        o o o o
!      O   O                   o o o o
!
!      four blocks             twelve cells
!
! Every child of one parent is joined to its siblings, and a face
! between two parents becomes a face between one child of each. That
! second rule is a choice: a real mesh refiner would know the geometry
! and join the children that actually touch. This one keeps the shape
! and the connectivity honest without pretending to know where
! anything is.
!
! What it is for: sharpening a mesh where an error measure says it is
! needed, and carrying a coarse multigrid correction back down.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_refiner

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_refiner, graph, graph_data, graph_field
  use class_graph         , only : stored_graph
  use class_graph_support , only : vertex_support
  use class_graph_field   , only : vertex_field

  implicit none

  private
  public :: refiner

  !===================================================================!
  ! One refiner: how many children each cell splits into.
  !
  ! The children of coarse cell v are numbered
  !
  !      (v - 1) * split + 1  ...  v * split
  !
  ! so a child's parent is (child - 1) / split + 1, and neither
  ! direction needs a stored map.
  !===================================================================!

  type, extends(graph_refiner) :: refiner

     integer :: split = 2

   contains

     procedure :: defined_on_graph => r_defined_on_graph
     procedure :: defined_on_data  => r_defined_on_data
     procedure :: refine_graph     => r_refine_graph
     procedure :: refine_data      => r_refine_data

  end type refiner

  interface refiner
     module procedure create
  end interface refiner

contains

  pure type(refiner) function create(split) result(this)

    integer, intent(in) :: split

    this % split = max(split, 1)

  end function create

  pure logical function r_defined_on_graph(this, input_graph)

    class(refiner), intent(in) :: this
    class(graph)  , intent(in) :: input_graph

    r_defined_on_graph = input_graph % num_vertices() > 0 .and. this % split >= 1

  end function r_defined_on_graph

  pure logical function r_defined_on_data(this, input_graph, input_data)

    class(refiner)   , intent(in) :: this
    class(graph)     , intent(in) :: input_graph
    class(graph_data), intent(in) :: input_data

    r_defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (vertex_field)
       r_defined_on_data = r_defined_on_data .and. input_data % num_entries() >= 0
    class default
       r_defined_on_data = .false.
    end select

  end function r_defined_on_data

  !===================================================================!
  ! R. Split every cell, join each family together, and carry every
  ! face down to one child on each side.
  !===================================================================!

  subroutine r_refine_graph(this, coarse_graph, fine_graph)

    class(refiner), intent(in)               :: this
    class(graph)  , intent(in)               :: coarse_graph
    class(graph)  , allocatable, intent(out) :: fine_graph

    integer, allocatable :: tails(:), heads(:)
    integer :: nv, ne, v, e, i, j, n, room

    nv = coarse_graph % num_vertices()
    ne = coarse_graph % num_edges()

    ! Room for every sibling pair plus one child face per coarse face.
    room = nv * this % split * this % split + ne
    allocate(tails(room), heads(room))
    n = 0

    ! Siblings are joined to each other, so a family stays connected.
    do v = 1, nv
       do i = 1, this % split - 1
          do j = i + 1, this % split
             n = n + 1
             tails(n) = child_of(v, i, this % split)
             heads(n) = child_of(v, j, this % split)
          end do
       end do
    end do

    ! A coarse face lands between the last child of its tail and the
    ! first child of its head. A wall stays a wall, on the first child.
    do e = 1, ne
       n = n + 1
       tails(n) = child_of(coarse_graph % edge_tail(e), this % split, this % split)
       if (coarse_graph % edge_has_head(e)) then
          heads(n) = child_of(coarse_graph % edge_head(e), 1, this % split)
       else
          heads(n) = 0
       end if
    end do

    allocate(fine_graph, source = &
         & stored_graph(nv * this % split, tails=tails(1:n), heads=heads(1:n), &
         &              number=coarse_graph % id()))

  end subroutine r_refine_graph

  !===================================================================!
  ! The i-th child of coarse cell v.
  !===================================================================!

  pure integer function child_of(v, i, split)

    integer, intent(in) :: v, i, split

    child_of = (v - 1) * split + i

  end function child_of

  !===================================================================!
  ! Carry the values down. Every child starts from its parent's value.
  !
  ! This is injection, the simplest honest choice: a child of a cell
  ! holding 3.0 holds 3.0. A geometric refiner would interpolate so
  ! the result stayed smooth across the new faces, but that needs to
  ! know where the children are, and this one does not.
  !===================================================================!

  subroutine r_refine_data(this, coarse_graph, coarse_data, fine_graph, fine_data)

    class(refiner)   , intent(in)               :: this
    class(graph)     , intent(in)               :: coarse_graph
    class(graph_data), intent(in)               :: coarse_data
    class(graph)     , intent(in)               :: fine_graph
    class(graph_data), allocatable, intent(out) :: fine_data

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    integer , allocatable :: indices(:)
    real(dp), allocatable :: cv(:), fv(:)
    integer :: nfine, ncomp, v, i, c, child

    select type (coarse_data)
    class is (vertex_field)

       nfine = fine_graph % num_vertices()
       ncomp = coarse_data % num_components()

       allocate(indices(nfine))
       do v = 1, nfine
          indices(v) = v
       end do

       on  = vertex_support(indices)
       out = vertex_field(coarse_data % name(), on, ncomp=ncomp, &
            &             unit_name=coarse_data % units())

       call coarse_data % get_real_vector(cv)
       allocate(fv(nfine * ncomp))
       fv = 0.0_dp

       do v = 1, coarse_graph % num_vertices()
          do i = 1, this % split
             child = child_of(v, i, this % split)
             do c = 1, ncomp
                associate (to => (child - 1) * ncomp + c, from => (v - 1) * ncomp + c)
                  if (to <= size(fv) .and. from <= size(cv)) fv(to) = cv(from)
                end associate
             end do
          end do
       end do

       call out % set_real_vector(fv)
       allocate(fine_data, source=out)

    end select

  end subroutine r_refine_data

end module class_graph_refiner
