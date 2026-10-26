!=====================================================================!
! Module that builds a mesh. A gmsh file becomes a mesh in one call,
!
!      m = mesh_from_gmsh('square-10.msh')
!
! The pipeline is the prime decomposition stated in
! view_mesh_geometry. The loader parses the file into member sets, the
! cell-to-vertex relation, the vertex coordinates and the tag names;
! the geometry module derives the face set, the incidence relations and
! every measurement; and this builder assembles the results as the
! mesh, the directed view whose vertices are cells and whose edges are
! the two-cell faces, with the measurements as fields and the tag names
! on the boundary edges:
!
!   file.msh --> loader --> cells -> vertices, coordinates, tags
!                              |
!                              v
!                           geometry --> faces, incidence, measurements
!                              |
!                              v
!                           builder --> the mesh
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_builder

  use util_precision   , only : dp
  use view_mesh        , only : mesh, values_of
  use util_string      , only : string
  use view_mesh_loader , only : mesh_loader
  use view_gmsh_loader , only : gmsh_loader
  use view_mesh_geometry, only : element_dimension
  use relation_binary  , only : transpose_padded
  use view_mesh_geometry, only : derive_faces, derive_face_cells, &
       & mesh_from_incidence

  implicit none

  private

  public :: mesh_from_gmsh
  public :: mesh_from_points
  public :: mesh_of_corners

contains

  !===================================================================!
  ! Returns the mesh of a gmsh file: load, derive, assemble. Tails and
  ! heads come from the face-to-cell relation: a face with one cell is
  ! a boundary face, an edge without a head. The per-face normal is the
  ! normal oriented from its tail cell; the per-face weight is the tail
  ! cell's interpolation share. Tag names are assigned to the boundary
  ! faces, read from the file's tag table by tag number.
  !
  ! A periodic link of the file identifies an entity with its master
  ! under a translation T, master + T = image. Each boundary face whose
  ! vertices are all images under one link is merged into the master
  ! face with those vertices' masters: the master face keeps its
  ! geometry, takes the image face's cell as its head, and stores the
  ! shift -T, the translation of the head to the master's side; the
  ! image face is removed. A link under which an image face has no
  ! master face is invalid input.
  !
  !   master face                          image face, + T
  !   cell p |                             | cell q
  !
  !   merged:  the master face, tail p, head q, shift -T;
  !            the image face removed
  !===================================================================!

  impure type(mesh) function mesh_from_gmsh(filename) result(m)

    character(len=*), intent(in) :: filename

    class(mesh_loader), allocatable :: gl

    ! The raw data the file supplies
    integer                   :: num_vertices, num_edges, bnum_faces, num_cells, num_tags
    integer , allocatable     :: vertex_numbers(:), vertex_tags(:)
    real(dp), allocatable     :: vertices(:,:)
    integer , allocatable     :: edge_numbers(:), edge_tags(:)
    integer , allocatable     :: edge_vertices(:,:), num_edge_vertices(:)
    integer , allocatable     :: edge_types(:)
    integer , allocatable     :: bface_numbers(:), bface_tags(:)
    integer , allocatable     :: bface_vertices(:,:), bnum_face_vertices(:)
    integer , allocatable     :: bface_types(:)
    integer , allocatable     :: cell_numbers(:), cell_tags(:)
    integer , allocatable     :: cell_vertices(:,:), num_cell_vertices(:)
    integer , allocatable     :: cell_types(:)
    integer , allocatable     :: tag_numbers(:), tag_physical_dimensions(:)
    type(string), allocatable :: tag_info(:)

    ! The derived relations and measurements
    integer               :: spatial_dim, num_faces
    integer , allocatable :: face_vertices(:,:), num_face_vertices(:)
    integer , allocatable :: face_tags(:), face_types(:)
    integer , allocatable :: vertex_cells(:,:), num_vertex_cells(:)
    integer , allocatable :: face_cells(:,:), num_face_cells(:)

    character(len=64), allocatable :: etags(:)
    integer                        :: f, row

    ! The periodic links and the shift of every face
    integer , allocatable :: first(:), last(:), image(:), master(:), remaining(:)
    real(dp), allocatable :: translation(:,:), face_shift(:,:)

    allocate(gl, source=gmsh_loader(trim(filename)))

    call gl % mesh_data( &
         & num_vertices, vertex_numbers, vertex_tags, vertices, &
         & num_edges, edge_numbers, edge_tags, edge_vertices, &
         & num_edge_vertices, &
         & bnum_faces, bface_numbers, bface_tags, bface_vertices, &
         & bnum_face_vertices, &
         & num_cells, cell_numbers, cell_tags, cell_vertices, &
         & num_cell_vertices, &
         & cell_types, bface_types, edge_types, &
         & num_tags, tag_numbers, tag_physical_dimensions, tag_info)

    spatial_dim = maxval(element_dimension(cell_types))

    ! The face member set and its relations
    call derive_faces(spatial_dim, num_vertices, num_cells, &
         & cell_vertices, num_cell_vertices, cell_types, cell_tags, &
         & bnum_faces, bface_vertices, bnum_face_vertices, bface_tags, &
         & bface_types, num_faces, face_vertices, num_face_vertices, &
         & face_tags, face_types)

    call transpose_padded(cell_vertices, num_cell_vertices, &
         & num_vertices, vertex_cells, num_vertex_cells)

    call derive_face_cells(num_faces, face_vertices, num_face_vertices, &
         & cell_vertices, num_cell_vertices, vertex_cells, &
         & num_vertex_cells, face_cells, num_face_cells)

    ! The periodic faces merged, the image faces removed
    select type (gl)
    type is (gmsh_loader)
       call gl % periodic_links(vertex_numbers, first, last, image, master, translation)
    end select
    allocate(face_shift(spatial_dim, num_faces), source=0.0_dp)
    call merge_periodic_faces(first, last, image, master, translation, face_vertices, &
         & num_face_vertices, face_cells, num_face_cells, face_shift, remaining)
    num_faces         = size(remaining)
    face_vertices     = face_vertices(:, remaining)
    num_face_vertices = num_face_vertices(remaining)
    face_tags         = face_tags(remaining)
    face_types        = face_types(remaining)
    face_cells        = face_cells(:, remaining)
    num_face_cells    = num_face_cells(remaining)
    face_shift        = face_shift(:, remaining)

    ! Tag names on the boundary faces, blank inside; a name is found by
    ! its physical number and dimension, not by its line in the file
    allocate(etags(num_faces))
    etags = ''
    do f = 1, num_faces
       if (num_face_cells(f) < 2) then
          row = findloc(tag_numbers, face_tags(f), dim=1, &
               & mask=tag_physical_dimensions == spatial_dim - 1)
          if (row > 0) etags(f) = tag_info(row) % str
       end if
    end do

    ! The measurements, and the mesh
    m = mesh_from_incidence(spatial_dim, vertices(1:spatial_dim, :), cell_vertices, num_cell_vertices, &
         & face_vertices, num_face_vertices, face_cells, num_face_cells, etags, &
         & face_shift=face_shift, cell_types=cell_types)

  end function mesh_from_gmsh

  !===================================================================!
  ! Merges the periodic faces. For each link, the boundary faces whose
  ! vertices are all images: the master face is the boundary face whose
  ! vertex set is the set of masters. face_cells and num_face_cells of
  ! the master face take the image's cell as head; face_shift of the
  ! master face is -T. remaining lists the faces that remain, in order.
  !===================================================================!

  pure subroutine merge_periodic_faces(first, last, image, master, translation, face_vertices, &
       & num_face_vertices, face_cells, num_face_cells, face_shift, remaining)

    integer , intent(in)               :: first(:), last(:), image(:), master(:)
    real(dp), intent(in)               :: translation(:,:)
    integer , intent(in)               :: face_vertices(:,:), num_face_vertices(:)
    integer , intent(inout)            :: face_cells(:,:), num_face_cells(:)
    real(dp), intent(inout)            :: face_shift(:,:)
    integer , allocatable, intent(out) :: remaining(:)

    logical, allocatable :: removed(:)
    integer, allocatable :: masters(:)
    integer              :: num_faces, f, g, l, i, k, at, d

    num_faces = size(num_face_vertices)
    d         = size(face_shift, 1)
    allocate(removed(num_faces), source=.false.)

    do l = 1, size(first)
       do f = 1, num_faces
          if (num_face_cells(f) /= 1 .or. removed(f)) cycle

          ! The masters of the face's vertices under this link
          masters = spread(0, 1, num_face_vertices(f))
          do k = 1, num_face_vertices(f)
             do i = first(l), last(l)
                if (image(i) == face_vertices(k, f)) masters(k) = master(i)
             end do
          end do
          if (any(masters == 0)) cycle

          ! The boundary face with that vertex set
          at = 0
          do g = 1, num_faces
             if (g == f .or. num_face_cells(g) /= 1 .or. removed(g)) cycle
             if (num_face_vertices(g) /= num_face_vertices(f)) cycle
             if (same_set(face_vertices(1:num_face_vertices(g), g), masters)) at = g
          end do
          if (at == 0) then
             error stop 'view_mesh_builder: a periodic image face has no master face with the &
                  &masters of its vertices'
          end if

          face_cells(2, at)    = face_cells(1, f)
          num_face_cells(at)   = 2
          face_shift(:, at)    = -translation(1:d, l)
          removed(f)           = .true.
       end do
    end do

    remaining = pack([(f, f = 1, num_faces)], .not. removed)

  end subroutine merge_periodic_faces

  !===================================================================!
  ! Returns whether two lists hold the same set of vertices
  !===================================================================!

  pure logical function same_set(a, b)
    integer, intent(in) :: a(:), b(:)
    integer             :: k
    same_set = size(a) == size(b)
    if (.not. same_set) return
    do k = 1, size(a)
       if (all(b /= a(k))) same_set = .false.
    end do
  end function same_set

  !===================================================================!
  ! Returns the vertex-centred mesh of points on a line, ascending.
  ! Each point is the centre of one cell, the cell between the
  ! midpoints to its neighbours, so that the cell lengths are the
  ! trapezoidal weights of the points. An interior face lies at each
  ! midpoint; a boundary face lies at each end point, tagged lower and
  ! upper, and is at the centre of its cell, so its delta is zero:
  !
  !   points      x1          x2          x3          x4
  !   cells    |--x1--|-----x2-----|-----x3-----|--x4--|
  !   faces   lower  (x1+x2)/2   (x2+x3)/2   (x3+x4)/2  upper
  !   lengths   (x2-x1)/2  (x3-x1)/2   (x4-x2)/2  (x4-x3)/2
  !
  ! Invalid input: fewer than two points, or points that do not
  ! ascend.
  !===================================================================!

  impure type(mesh) function mesh_from_points(positions) result(m)

    real(dp), intent(in) :: positions(:)

    integer , allocatable         :: tails(:), heads(:)
    real(dp), allocatable         :: volumes(:), areas(:), deltas(:), normals(:), face_centres(:), weights(:)
    character(len=5), allocatable :: tags(:)
    integer                       :: n, f

    n = size(positions)
    if (n < 2) error stop 'view_mesh_builder: a line of points requires two points at least'
    if (any(positions(2:) <= positions(:n - 1))) then
       error stop 'view_mesh_builder: the points of a line ascend'
    end if

    ! The interior faces, then the boundary faces at the two ends
    allocate(tails(n + 1), heads(n + 1), areas(n + 1), deltas(n + 1), normals(n + 1), &
         & face_centres(n + 1), weights(n + 1), tags(n + 1))
    do f = 1, n - 1
       tails(f)        = f
       heads(f)        = f + 1
       face_centres(f) = (positions(f) + positions(f + 1)) / 2.0_dp
       deltas(f)       = positions(f + 1) - positions(f)
       normals(f)      = 1.0_dp
       weights(f)      = 0.5_dp
       tags(f)         = ''
    end do
    tails(n : n + 1)        = [1, n]
    heads(n : n + 1)        = 0
    face_centres(n : n + 1) = [positions(1), positions(n)]
    deltas(n : n + 1)       = 0.0_dp
    normals(n : n + 1)      = [-1.0_dp, 1.0_dp]
    weights(n : n + 1)      = 1.0_dp
    tags(n : n + 1)         = ['lower', 'upper']
    areas = 1.0_dp

    ! The cells' lengths, from midpoint to midpoint
    allocate(volumes(n))
    volumes(1)         = (positions(2) - positions(1)) / 2.0_dp
    volumes(2 : n - 1) = (positions(3 : n) - positions(1 : n - 2)) / 2.0_dp
    volumes(n)         = (positions(n) - positions(n - 1)) / 2.0_dp

    m = mesh(n, tails=tails, heads=heads, volumes=volumes, cell_centres=positions, areas=areas, &
         & deltas=deltas, normals=normals, face_centres=face_centres, weights=weights, etags=tags, &
         & dimension=1)

  end function mesh_from_points

  !===================================================================!
  ! Returns the corners of a mesh as its points: the mesh whose cells
  ! are the corners some cell of m reads, numbered in the order of the
  ! corners, each at its corner's position and with a share of the
  ! volumes of the cells around it, a cell's volume over its corner
  ! count. Its faces join every two corners of a cell, once. The cells
  ! of m are its elements, their corners renumbered:
  !
  !     4 ----- 3
  !     |       |       a cell of m, corners 1 2 3 4: each corner a
  !     |   e   |       point of the result, with a quarter of e's
  !     |       |       volume; the faces 1-2 1-3 1-4 2-3 2-4 3-4
  !     1 ----- 2
  !
  ! Invalid input: a mesh that keeps no corners, and a periodic mesh,
  ! whose corners on the two sides of a periodic face are not
  ! identified.
  !===================================================================!

  impure type(mesh) function mesh_of_corners(m) result(c)

    type(mesh), intent(in) :: m

    integer , allocatable :: number(:), count_of(:), start(:), other(:), tails(:), heads(:)
    real(dp), allocatable :: volume(:), centre(:), shift(:), cell_volume(:)
    real(dp), allocatable :: areas(:), deltas(:), normals(:), face_centres(:), weights(:)
    real(dp)              :: span(3)
    integer               :: d, n, ncell, e, i, j, a, b, k, f, nf, low, high

    if (.not. m % has_corners()) then
       error stop 'view_mesh_builder: the corners of a mesh are its points where the mesh keeps them, as a mesh &
            &read from a file does'
    end if
    call values_of(m % face_shift(), shift)
    if (any(shift /= 0.0_dp)) then
       error stop 'view_mesh_builder: the corners on the two sides of a periodic face are not identified'
    end if
    d     = m % dimension
    ncell = size(m % num_cell_corners)

    ! The corners some cell reads, numbered in order
    allocate(number(size(m % corners, 2)), source=0)
    do e = 1, ncell
       number(m % cell_corners(1:m % num_cell_corners(e), e)) = 1
    end do
    n = 0
    do k = 1, size(number)
       if (number(k) == 0) cycle
       n = n + 1
       number(k) = n
    end do

    ! The positions, and each cell's volume shared among its corners
    allocate(centre(d * n), volume(n), source=0.0_dp)
    do k = 1, size(number)
       if (number(k) > 0) centre(d * (number(k) - 1) + 1:d * number(k)) = m % corners(1:d, k)
    end do
    call values_of(m % cell_volume(), cell_volume)
    do e = 1, ncell
       do i = 1, m % num_cell_corners(e)
          k = number(m % cell_corners(i, e))
          volume(k) = volume(k) + cell_volume(e) / real(m % num_cell_corners(e), dp)
       end do
    end do

    ! The faces: for each corner the corners above it that share a cell
    ! with it, each once, gathered by rows as start and other
    allocate(count_of(n), source=0)
    do e = 1, ncell
       do i = 1, m % num_cell_corners(e)
          do j = 1, m % num_cell_corners(e)
             a = number(m % cell_corners(i, e))
             b = number(m % cell_corners(j, e))
             if (a < b) count_of(a) = count_of(a) + 1
          end do
       end do
    end do
    allocate(start(n + 1))
    start(1) = 1
    do k = 1, n
       start(k + 1) = start(k) + count_of(k)
    end do
    allocate(other(start(n + 1) - 1))
    count_of = 0
    do e = 1, ncell
       do i = 1, m % num_cell_corners(e)
          do j = 1, m % num_cell_corners(e)
             a = number(m % cell_corners(i, e))
             b = number(m % cell_corners(j, e))
             if (a >= b) cycle
             other(start(a) + count_of(a)) = b
             count_of(a) = count_of(a) + 1
          end do
       end do
    end do
    allocate(tails(size(other)), heads(size(other)))
    nf = 0
    do a = 1, n
       low  = start(a)
       high = start(a + 1) - 1
       do f = low, high
          if (any(other(low:f - 1) == other(f))) cycle
          nf = nf + 1
          tails(nf) = a
          heads(nf) = other(f)
       end do
    end do

    ! Each face's measurements, from the corners it joins
    allocate(areas(nf), deltas(nf), normals(d * nf), face_centres(d * nf), weights(nf))
    do f = 1, nf
       span(1:d) = centre(d * (heads(f) - 1) + 1:d * heads(f)) - centre(d * (tails(f) - 1) + 1:d * tails(f))
       areas(f)  = 1.0_dp
       deltas(f) = norm2(span(1:d))
       normals(d * (f - 1) + 1:d * f)      = span(1:d) / deltas(f)
       face_centres(d * (f - 1) + 1:d * f) = centre(d * (tails(f) - 1) + 1:d * tails(f)) + span(1:d) / 2.0_dp
       weights(f) = 0.5_dp
    end do

    c = mesh(n, tails=tails(1:nf), heads=heads(1:nf), volumes=volume, cell_centres=centre, areas=areas, &
         & deltas=deltas, normals=normals, face_centres=face_centres, weights=weights, dimension=d)

    ! The cells of m, the elements of their corners
    allocate(c % corners(3, n), source=0.0_dp)
    do k = 1, size(number)
       if (number(k) > 0) c % corners(:, number(k)) = m % corners(:, k)
    end do
    c % cell_corners = m % cell_corners
    do e = 1, ncell
       do i = 1, m % num_cell_corners(e)
          c % cell_corners(i, e) = number(m % cell_corners(i, e))
       end do
    end do
    c % num_cell_corners = m % num_cell_corners
    c % cell_kinds       = m % cell_kinds

  end function mesh_of_corners

end module view_mesh_builder
