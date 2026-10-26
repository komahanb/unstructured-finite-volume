!=====================================================================!
! Module that states the concrete graph refiner.
!
! R is the inverse direction of C. One cell becomes several:
!
!      O   O                   o o o o
!               ------>        o o o o
!      O   O                   o o o o
!
!      four blocks             twelve cells, split = 3
!
! Every child of one parent is joined to its siblings, and a face
! between two parents becomes a face between one child of each. The
! second rule is the strongest statement the connectivity supports: a
! geometric refiner, with positions available, joins exactly the
! children that are adjacent; this refiner has no positions, so it
! preserves the shape and asserts nothing beyond it.
!
! Its purpose is to refine a mesh where an error measure requires it,
! and to prolongate a coarse multigrid correction to the fine level.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_refiner

  use util_precision      , only : dp
  use view_directed       , only : directed_graph
  use field_calculus      , only : field
  use transform_structure , only : transform, through_blocks
  use view_directed_stored, only : stored_directed_graph
  use field_stored        , only : stored_field

  implicit none

  private

  public :: refiner

  !-------------------------------------------------------------------!
  ! Type that models the refiner: the transform from one cell to
  ! split cells.
  !
  ! The pair's identity is one-sided, coarsen(refine(G)) = G, and only
  ! in that direction, because refinement invents detail that
  ! coarsening cannot recover. refine_data states how one coarse value
  ! maps onto the new cells: copied, or interpolated.
  !
  ! The children of coarse cell v are numbered
  !
  !      (v - 1) * split + 1  ...  v * split
  !
  ! so a child's parent is (child - 1) / split + 1, and neither
  ! direction needs a stored map:
  !
  !   coarse cell     1          2          3
  !   children     1 2 3      4 5 6      7 8 9       split = 3
  !-------------------------------------------------------------------!

  type, extends(transform) :: refiner

     integer :: split = 2

   contains

     procedure :: defined_on_graph
     procedure :: refine_graph
     procedure :: refine_data

  end type refiner

  interface refiner
     module procedure create
  end interface refiner

contains

  !===================================================================!
  ! Creates a refiner that splits every cell into split cells. A split
  ! below one is taken as one: the identity refinement.
  !===================================================================!

  pure type(refiner) function create(split) result(this)

    integer, intent(in) :: split

    this % split = max(split, 1)

  end function create

  !===================================================================!
  ! Returns whether this refiner is defined on that graph: any graph
  ! with at least one cell qualifies
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(refiner)       , intent(in) :: this
    class(directed_graph), intent(in) :: input_graph

    defined_on_graph = input_graph % num_vertices() > 0 .and. this % split >= 1

  end function defined_on_graph

  !===================================================================!
  ! Refines the graph, R: splits every cell, joins the children of each
  ! cell to one another, and maps every face to one child on each side:
  !
  !   coarse      1 ------------ e ------------> 2          split = 3
  !
  !   fine        1   2   3                      4   5   6
  !               siblings joined:  1-2, 1-3, 2-3     4-5, 4-6, 5-6
  !               the face e:       3 -----------> 4
  !                                 the last child of the tail to
  !                                 the first child of the head
  !===================================================================!

  subroutine refine_graph(this, coarse_graph, fine_graph)

    class(refiner)                    , intent(in)  :: this
    class(directed_graph)             , intent(in)  :: coarse_graph
    class(directed_graph), allocatable, intent(out) :: fine_graph

    integer, allocatable :: tails(:), heads(:)
    integer              :: nv, ne, v, e, i, j, n, capacity

    nv = coarse_graph % num_vertices()
    ne = coarse_graph % num_edges()

    ! Space for every sibling pair, and one child face per coarse face
    capacity = nv * this % split * this % split + ne
    allocate(tails(capacity), heads(capacity))
    n = 0

    ! Siblings are joined to each other, so the children of one cell
    ! stay connected
    do v = 1, nv
       do i = 1, this % split - 1
          do j = i + 1, this % split
             n = n + 1
             tails(n) = child_of(v, i, this % split)
             heads(n) = child_of(v, j, this % split)
          end do
       end do
    end do

    ! A coarse face maps to an edge between the last child of its tail
    ! and the first child of its head. A boundary face stays a boundary
    ! face, on the tail's last child.
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
         & stored_directed_graph(nv * this % split, tails=tails(1:n), heads=heads(1:n), &
         &              number=coarse_graph % id()))

  end subroutine refine_graph

  !===================================================================!
  ! Returns the i-th child of coarse cell v
  !===================================================================!

  pure integer function child_of(v, i, split)

    integer, intent(in) :: v, i, split

    child_of = (v - 1) * split + i

  end function child_of

  !===================================================================!
  ! Prolongates the values: every child takes its parent's value.
  !
  ! This is injection, the transpose of the block sum the coarsener
  ! takes, read through the same map: a child of a cell storing 3.0
  ! stores 3.0. A geometric refiner interpolates so that the result
  ! stays smooth across the new faces, but interpolation requires the
  ! children's positions, and this refiner has none.
  !===================================================================!

  subroutine refine_data(this, coarse_graph, coarse_data, fine_graph, fine_data)

    class(refiner)           , intent(in)  :: this
    class(directed_graph)    , intent(in)  :: coarse_graph
    class(field)             , intent(in)  :: coarse_data
    class(directed_graph)    , intent(in)  :: fine_graph
    class(field), allocatable, intent(out) :: fine_data

    type(stored_field)    :: out
    real(dp), allocatable :: cv(:), fv(:)
    integer               :: nfine, num_components, v

    select type (coarse_data)
    class is (stored_field)

       nfine = fine_graph % num_vertices()
       num_components = coarse_data % num_components()

       out = stored_field(coarse_data % name(), fine_graph % vertex_set(), fine_graph % num_vertices(), &
            &             num_components=num_components, unit_name=coarse_data % units())

       ! The block map read backward: every child its parent's value
       call coarse_data % real_vector(cv)
       call through_blocks([((v - 1) / this % split + 1, v = 1, nfine)], &
            & coarse_graph % num_vertices(), num_components, fv, cv, transposed=.true.)

       call out % set_real_vector(fv)
       allocate(fine_data, source=out)

    class default
       error stop 'refine: coarse_data''s dynamic type is not class(stored_field), &
            &which refine_data requires'
    end select

  end subroutine refine_data

end module transform_refiner
