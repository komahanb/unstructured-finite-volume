!=====================================================================!
! The mesh builder: a gmsh file becomes a mesh in one call,
!
!      m = mesh_from_gmsh('square-10.msh')
!
! The pipeline is the prime decomposition stated in
! view_mesh_geometry: the loader parses the file into member sets,
! the cell-to-vertex relation, the vertex coordinates, and the tag
! names; the geometry module derives the face set, the incidence
! relations, and every measurement; this builder assembles the
! results as the mesh - the directed view whose vertices are cells and
! whose edges are the two-cell faces, with the measurements as
! fields and the tag names on the boundary edges.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_builder

  use util_precision  , only : dp
  use view_mesh     , only : mesh
  use util_string         , only : string
  use view_mesh_loader, only : mesh_loader
  use view_gmsh_loader    , only : gmsh_loader
  use view_mesh_geometry  , only : element_dimension
  use relation_binary, only : transpose_padded
  use view_mesh_geometry  , only : derive_faces, derive_face_cells, &
       & mesh_from_incidence

  implicit none

  private
  public :: mesh_from_gmsh

contains

  !===================================================================!
  ! Load, derive, assemble. Tails and heads come from the face-to-cell
  ! relation - a face with one cell is a boundary face, an edge
  ! without a head. The per-face normal is the normal oriented from
  ! its tail cell; the per-face weight is the tail cell's
  ! interpolation share. Tag names are assigned to the boundary faces,
  ! read from the file's tag table by tag number.
  !
  ! A periodic link of the file identifies an entity with its master
  ! under a translation T, master + T = image. Each boundary face
  ! whose vertices are all images under one link is merged into the
  ! master face with those vertices' masters: the master face keeps
  ! its geometry, takes the image face's cell as its head, and stores
  ! the shift -T, the translation of the head to the master's side;
  ! the image face is removed. A link under which an image face has
  ! no master face is invalid input.
  !===================================================================!

  impure type(mesh) function mesh_from_gmsh(filename) result(m)

    character(len=*), intent(in) :: filename

    class(mesh_loader), allocatable :: gl

    ! the raw data the file supplies
    integer :: num_vertices, num_edges, bnum_faces, num_cells, num_tags
    integer , allocatable :: vertex_numbers(:), vertex_tags(:)
    real(dp), allocatable :: vertices(:,:)
    integer , allocatable :: edge_numbers(:), edge_tags(:)
    integer , allocatable :: edge_vertices(:,:), num_edge_vertices(:)
    integer , allocatable :: edge_types(:)
    integer , allocatable :: bface_numbers(:), bface_tags(:)
    integer , allocatable :: bface_vertices(:,:), bnum_face_vertices(:)
    integer , allocatable :: bface_types(:)
    integer , allocatable :: cell_numbers(:), cell_tags(:)
    integer , allocatable :: cell_vertices(:,:), num_cell_vertices(:)
    integer , allocatable :: cell_types(:)
    integer , allocatable :: tag_numbers(:), tag_physical_dimensions(:)
    type(string), allocatable :: tag_info(:)

    ! the derived relations and measurements
    integer :: spatial_dim, num_faces
    integer , allocatable :: face_vertices(:,:), num_face_vertices(:)
    integer , allocatable :: face_tags(:), face_types(:)
    integer , allocatable :: vertex_cells(:,:), num_vertex_cells(:)
    integer , allocatable :: face_cells(:,:), num_face_cells(:)

    character(len=64), allocatable :: etags(:)
    integer :: f, tag

    ! the periodic links and the shift of every face
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

    ! the face member set and its relations
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


    ! the periodic faces merged, the image faces removed
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

    ! tag names on the boundary faces, blank inside
    allocate(etags(num_faces))
    etags = ''
    do f = 1, num_faces
       if (num_face_cells(f) < 2) then
          tag = face_tags(f)
          if (tag >= 1 .and. tag <= size(tag_info)) then
             etags(f) = tag_info(tag) % str
          end if
       end if
    end do

    ! the measurements, and the mesh
    m = mesh_from_incidence(spatial_dim, vertices(1:spatial_dim, :), cell_vertices, num_cell_vertices, &
         & face_vertices, num_face_vertices, face_cells, num_face_cells, etags, &
         & face_shift=face_shift, cell_types=cell_types)

  end function mesh_from_gmsh

  !===================================================================!
  ! THE MERGE. For each link, the boundary faces whose vertices are
  ! all images: the master face is the boundary face whose vertex set
  ! is the set of masters. face_cells and num_face_cells of the master
  ! face take the image's cell as head; face_shift of the master face
  ! is -T. remaining lists the faces that remain, in order.
  !===================================================================!

  pure subroutine merge_periodic_faces(first, last, image, master, translation, face_vertices, &
       & num_face_vertices, face_cells, num_face_cells, face_shift, remaining)

    integer , intent(in)    :: first(:), last(:), image(:), master(:)
    real(dp), intent(in)    :: translation(:,:)
    integer , intent(in)    :: face_vertices(:,:), num_face_vertices(:)
    integer , intent(inout) :: face_cells(:,:), num_face_cells(:)
    real(dp), intent(inout) :: face_shift(:,:)
    integer , allocatable, intent(out) :: remaining(:)

    logical, allocatable :: removed(:)
    integer, allocatable :: masters(:)
    integer :: num_faces, f, g, l, i, k, at, d

    num_faces = size(num_face_vertices)
    d         = size(face_shift, 1)
    allocate(removed(num_faces), source=.false.)

    do l = 1, size(first)
       do f = 1, num_faces
          if (num_face_cells(f) /= 1 .or. removed(f)) cycle

          ! the masters of the face's vertices under this link
          masters = spread(0, 1, num_face_vertices(f))
          do k = 1, num_face_vertices(f)
             do i = first(l), last(l)
                if (image(i) == face_vertices(k, f)) masters(k) = master(i)
             end do
          end do
          if (any(masters == 0)) cycle

          ! the boundary face with that vertex set
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

  pure logical function same_set(a, b)
    integer, intent(in) :: a(:), b(:)
    integer :: k
    same_set = size(a) == size(b)
    if (.not. same_set) return
    do k = 1, size(a)
       if (all(b /= a(k))) same_set = .false.
    end do
  end function same_set

end module view_mesh_builder
