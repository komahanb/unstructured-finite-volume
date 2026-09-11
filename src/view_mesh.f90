!=====================================================================!
! The mesh: a stored graph embedded in space.
!
! LEVEL 1 OF THE STRATIFICATION - measurements are calculus
! content: values defined on structure, with no objective in them.
! This is the tower's one use of inheritance: the mesh IS a graph
! - cells are the vertices, interior faces the edges, boundary
! faces the edges without heads - so it extends the stored graph
! rather than containing one.
!
!     +-----+-----+
!     |  1  |  2  |            (1)---(2)
!     +-----+-----+             |     |        the same mesh, as
!     |  3  |  4  |            (3)---(4)       its graph
!     +-----+-----+
!
! Everything else is measurement defined on that graph, stored as
! typed fields with fixed names:
!
!      cell_volume()    one number per cell
!      cell_centre()    three per cell
!      face_area()      one per face
!      face_delta()     the centre-to-centre distance, one per face
!      face_normal()    three per face
!      face_centre()    three per face
!      face_weights()   the interpolation weight, one per face
!
! values_of(field, values) reads any of them as a real array in one
! statement.
!
! No string names any of these. The dictionary in
! geometry-to-operator-mapping.md specifies which operator argument
! each one supplies; an operator receives those numbers at
! construction and never reads the mesh.
!
! THE CHECK AT CONSTRUCTION. Every geometry array must match the
! structure it measures - one volume per cell, one area per face -
! and the constructor stops the program on a mismatch rather than
! storing inconsistent data. A mesh that is constructed is a mesh
! whose measurements match its structure.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh

  use iso_fortran_env, only : error_unit
  use util_precision  , only : dp
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph

  implicit none

  private
  public :: mesh, require, values_of

  !===================================================================!
  ! One mesh: the inherited structure, plus seven measurements.
  !===================================================================!

  type, extends(stored_directed_graph) :: mesh

     integer :: dimension = 3

     type(stored_field) :: volumes
     type(stored_field) :: cell_centres
     type(stored_field) :: areas
     type(stored_field) :: deltas
     type(stored_field) :: normals
     type(stored_field) :: face_centers_
     type(stored_field) :: weights
     ! the translation of a face's head cell into the face's frame: a
     ! face identifying two sides of a periodic box stores the
     ! period, every other face zero
     type(stored_field) :: shifts

   contains

     procedure :: cell_volume
     procedure :: cell_centre
     procedure :: face_area
     procedure :: face_delta
     procedure :: face_normal
     procedure :: face_centre
     procedure :: face_weights
     procedure :: face_shift
     procedure :: neighbourhood

  end type mesh

  interface mesh
     module procedure create
  end interface mesh

contains

  !===================================================================!
  ! Build a mesh from its structure and its measurements. The
  ! structure arguments are the stored graph's own; the geometry
  ! is passed as plain arrays, one entry per cell or per face, vector
  ! quantities three wide in entry order. The check: every array must
  ! match the structure, or the constructor stops.
  !===================================================================!

  impure type(mesh) function create(nv, tails, heads, volumes, &
       & cell_centres, areas, deltas, normals, face_centres, weights, &
       & vtags, etags, number, dimension, shifts) result(this)

    integer         , intent(in)           :: nv
    integer         , intent(in)           :: tails(:)
    integer         , intent(in)           :: heads(:)
    real(dp)        , intent(in)           :: volumes(:)
    real(dp)        , intent(in)           :: cell_centres(:)
    real(dp)        , intent(in)           :: areas(:)
    real(dp)        , intent(in)           :: deltas(:)
    real(dp)        , intent(in)           :: normals(:)
    real(dp)        , intent(in)           :: face_centres(:)
    real(dp)        , intent(in)           :: weights(:)
    character(len=*), intent(in), optional :: vtags(:)
    character(len=*), intent(in), optional :: etags(:)
    integer         , intent(in), optional :: number
    integer         , intent(in), optional :: dimension
    real(dp)        , intent(in), optional :: shifts(:)

    type(graph) :: cells, faces
    integer :: ne, d

    ! The structure first, through the parent's own constructor.
    this % stored_directed_graph = stored_directed_graph(nv, tails=tails, heads=heads, &
         & vtags=vtags, etags=etags, number=number)

    ne = this % num_edges()

    d = 3
    if (present(dimension)) d = dimension
    this % dimension = d

    ! The check: measurement sizes must match the structure they measure.
    call require(size(volumes)      == nv    , 'one volume per cell')
    call require(size(cell_centres) == d * nv, 'd centre parts per cell')
    call require(size(areas)        == ne    , 'one area per face')
    call require(size(deltas)       == ne    , 'one delta per face')
    call require(size(normals)      == d * ne, 'd normal parts per face')
    call require(size(face_centres) == d * ne, 'd centre parts per face')
    call require(size(weights)      == ne    , 'one weight per face')
    if (present(shifts)) call require(size(shifts) == d * ne, 'd shift parts per face')

    ! Geometry is defined on the graph's OWN carriers, so a field's
    ! domain returns the mesh identity every consumer compares against.
    cells = this % vertex_set()
    faces = this % edge_set()

    this % volumes       = stored_field('cell_volume' , cells, nv, unit_name='m3')
    this % cell_centres  = stored_field('cell_centre' , cells, nv, num_components=d, unit_name='m')
    this % areas         = stored_field('face_area'   , faces, ne, unit_name='m2')
    this % deltas        = stored_field('face_delta'  , faces, ne, unit_name='m')
    this % normals       = stored_field('face_normal' , faces, ne, num_components=d, unit_name='-')
    this % face_centers_ = stored_field('face_centre' , faces, ne, num_components=d, unit_name='m')
    this % weights       = stored_field('face_weights', faces, ne, unit_name='-')
    this % shifts        = stored_field('face_shift'  , faces, ne, num_components=d, unit_name='m')

    call this % volumes       % set_real_vector(volumes)
    call this % cell_centres  % set_real_vector(cell_centres)
    call this % areas         % set_real_vector(areas)
    call this % deltas        % set_real_vector(deltas)
    call this % normals       % set_real_vector(normals)
    call this % face_centers_ % set_real_vector(face_centres)
    call this % weights       % set_real_vector(weights)
    if (present(shifts)) then
       call this % shifts % set_real_vector(shifts)
    else
       call this % shifts % set_real_vector(spread(0.0_dp, 1, d * ne))
    end if

  end function create

  !===================================================================!
  ! The check itself: report which condition failed, then stop. A
  ! mesh with inconsistent measurements must not exist, and a file of
  ! incorrect cells must not be written: the writer checks through
  ! this procedure too.
  !===================================================================!

  subroutine require(satisfied, condition_description)

    logical         , intent(in) :: satisfied
    character(len=*), intent(in) :: condition_description

    if (satisfied) return

    write(error_unit, *) 'mesh check: expected ', condition_description
    error stop 'mesh: a precondition on the structure or its measurements failed'

  end subroutine require

  !===================================================================!
  ! The seven accessors. Each returns a copy of the stored field, so
  ! a caller may read it, pass its values to an operator, and never
  ! reference the mesh again.
  !===================================================================!

  type(stored_field) function cell_volume(this)

    class(mesh), intent(in) :: this

    cell_volume = this % volumes

  end function cell_volume

  type(stored_field) function cell_centre(this)

    class(mesh), intent(in) :: this

    cell_centre = this % cell_centres

  end function cell_centre

  type(stored_field) function face_area(this)

    class(mesh), intent(in) :: this

    face_area = this % areas

  end function face_area

  type(stored_field) function face_delta(this)

    class(mesh), intent(in) :: this

    face_delta = this % deltas

  end function face_delta

  type(stored_field) function face_normal(this)

    class(mesh), intent(in) :: this

    face_normal = this % normals

  end function face_normal

  type(stored_field) function face_centre(this)

    class(mesh), intent(in) :: this

    face_centre = this % face_centers_

  end function face_centre

  type(stored_field) function face_shift(this)

    class(mesh), intent(in) :: this

    face_shift = this % shifts

  end function face_shift

  !===================================================================!
  ! The cells around the seeds, ring by ring, each once, with the
  ! translation that places each in the first seed's frame: crossing
  ! a face from its tail to its head adds the face's shift, the other
  ! way subtracts it. As many rings as given, or until at least the
  ! count required is reached, or until no new cell is found. A cell
  ! reached again with another translation stops the program: the
  ! periodic box has too few cells for the reach.
  !===================================================================!

  subroutine neighbourhood(this, seeds, rings, at_least, members, offsets)

    class(mesh), intent(in) :: this
    integer    , intent(in) :: seeds(:), rings, at_least
    integer , allocatable, intent(out) :: members(:)
    real(dp), allocatable, intent(out) :: offsets(:,:)

    real(dp), allocatable :: shift(:), translation(:)
    integer , allocatable :: edges(:)
    integer :: d, r, k, e, f, other, before, at
    real(dp) :: sign

    d = this % dimension
    call values_of(this % shifts, shift)
    members = seeds
    allocate(offsets(d, size(seeds)), source=0.0_dp)
    ! the seeds share one frame: a second seed is the head or the
    ! tail of a face the first seed is on
    do k = 2, size(seeds)
       call this % incident_edges(seeds(1), edges)
       do e = 1, size(edges)
          f = edges(e)
          if (.not. this % edge_has_head(f)) cycle
          if (this % edge_tail(f) == seeds(1) .and. this % edge_head(f) == seeds(k)) then
             offsets(:, k) = shift(d * f - d + 1:d * f)
          else if (this % edge_head(f) == seeds(1) .and. this % edge_tail(f) == seeds(k)) then
             offsets(:, k) = -shift(d * f - d + 1:d * f)
          end if
       end do
    end do
    r = 0
    do
       if (rings > 0) then
          if (r >= rings) exit
       else
          if (r >= 1 .and. size(members) >= at_least) exit
       end if
       before = size(members)
       do k = 1, before
          call this % incident_edges(members(k), edges)
          do e = 1, size(edges)
             f = edges(e)
             if (.not. this % edge_has_head(f)) cycle
             if (this % edge_tail(f) == members(k)) then
                other = this % edge_head(f)
                sign  = 1.0_dp
             else
                other = this % edge_tail(f)
                sign  = -1.0_dp
             end if
             translation = offsets(:, k) + sign * shift(d * f - d + 1:d * f)
             at = findloc(members, other, dim=1)
             if (at > 0) then
                if (any(abs(offsets(:, at) - translation) > 0.0_dp)) then
                   error stop 'mesh: a periodic box stores more cells than the neighbourhood reaches around it'
                end if
                cycle
             end if
             members = [members, other]
             offsets = reshape([offsets, translation], [d, size(members)])
          end do
       end do
       r = r + 1
       if (size(members) == before) exit
    end do

  end subroutine neighbourhood

  type(stored_field) function face_weights(this)

    class(mesh), intent(in) :: this

    face_weights = this % weights

  end function face_weights

  !===================================================================!
  ! The real values of one measurement as an array, so that a caller
  ! reads them in one statement: call values_of(m % face_area(), a).
  !===================================================================!

  pure subroutine values_of(f, values)

    type(stored_field)   , intent(in)  :: f
    real(dp), allocatable, intent(out) :: values(:)

    call f % real_vector(values)

  end subroutine values_of

end module view_mesh
