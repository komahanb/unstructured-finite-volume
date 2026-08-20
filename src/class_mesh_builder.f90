!=====================================================================!
! The mesh builder: a gmsh file becomes a mesh in one call,
!
!      m = mesh_from_gmsh('square-10.msh')
!
! The pipeline is the prime decomposition stated in
! graph_mesh_geometry: the loader parses the file into member sets,
! the cell-to-vertex relation, the vertex coordinates, and the tag
! names; the geometry module derives the face set, the incidence
! relations, and every measurement; this builder seats the results
! as the mesh - the directed view whose vertices are cells and
! whose edges are the two-cell faces, with the measurements as
! fields and the tag names on the boundary edges.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_mesh_builder

  use iso_fortran_env      , only : dp => REAL64
  use class_graph_mesh     , only : mesh
  use util_string         , only : string
  use interface_mesh_loader, only : mesh_loader
  use class_gmsh_loader    , only : gmsh_loader
  use graph_mesh_geometry  , only : elem_type_dimension
  use relation_binary, only : transpose_padded
  use graph_mesh_geometry  , only : derive_faces, derive_face_cells, &
       & cell_centers_of, face_centers_areas_of, &
       & cell_face_normals_of, cell_volumes_of, centroidal_vectors_of, &
       & face_deltas_of, face_weights_of

  implicit none

  private
  public :: mesh_from_gmsh

contains

  !===================================================================!
  ! Load, derive, seat. Tails and heads come from the face-to-cell
  ! relation - a face with one cell is a boundary face, an edge
  ! without a head. The per-face normal is the one its tail cell
  ! sees; the per-face weight is the tail cell's interpolation
  ! share. Tag names are stamped on the boundary faces, read from
  ! the file's tag table by tag number.
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
    integer , allocatable :: cell_faces(:,:), num_cell_faces(:)
    real(dp), allocatable :: cell_centers(:,:), face_centers(:,:)
    real(dp), allocatable :: face_areas(:), cell_volumes(:)
    real(dp), allocatable :: cell_face_normals(:,:,:), lvec(:,:)
    real(dp), allocatable :: face_deltas(:), face_cell_weights(:,:)

    ! the seat's per-face arrays
    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: normals(:)
    real(dp), allocatable :: weights(:)
    character(len=64), allocatable :: etags(:)
    integer :: f, c, l, tag

    allocate(gl, source=gmsh_loader(trim(filename)))

    call gl % get_mesh_data( &
         & num_vertices, vertex_numbers, vertex_tags, vertices, &
         & num_edges, edge_numbers, edge_tags, edge_vertices, &
         & num_edge_vertices, &
         & bnum_faces, bface_numbers, bface_tags, bface_vertices, &
         & bnum_face_vertices, &
         & num_cells, cell_numbers, cell_tags, cell_vertices, &
         & num_cell_vertices, &
         & cell_types, bface_types, edge_types, &
         & num_tags, tag_numbers, tag_physical_dimensions, tag_info)

    spatial_dim = maxval(elem_type_dimension(cell_types))

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

    call transpose_padded(face_cells, num_face_cells, num_cells, &
         & cell_faces, num_cell_faces)

    ! the measurements
    call cell_centers_of(vertices, cell_vertices, num_cell_vertices, &
         & cell_centers)
    call face_centers_areas_of(spatial_dim, vertices, face_vertices, &
         & num_face_vertices, face_centers, face_areas)
    call cell_face_normals_of(spatial_dim, vertices, face_vertices, &
         & num_face_vertices, cell_faces, num_cell_faces, face_centers, &
         & cell_centers, cell_face_normals)
    call cell_volumes_of(spatial_dim, face_centers, face_areas, &
         & cell_face_normals, cell_faces, num_cell_faces, cell_volumes)
    call centroidal_vectors_of(face_cells, num_face_cells, cell_centers, &
         & face_centers, lvec)
    call face_deltas_of(num_faces, lvec, cell_face_normals, cell_faces, &
         & num_cell_faces, face_deltas)
    call face_weights_of(face_cells, num_face_cells, cell_centers, &
         & face_centers, face_cell_weights)

    ! tails and heads from the face-to-cell relation; a face with
    ! one cell is a boundary face, an edge without a head
    allocate(tails(num_faces), heads(num_faces))
    do f = 1, num_faces
       tails(f) = face_cells(1, f)
       heads(f) = 0
       if (num_face_cells(f) >= 2) heads(f) = face_cells(2, f)
    end do

    ! one outward normal per face, read from its tail cell's side
    allocate(normals(3 * num_faces))
    normals = 0.0_dp
    do c = 1, num_cells
       do l = 1, num_cell_faces(c)
          f = cell_faces(l, c)
          if (face_cells(1, f) == c) then
             normals(3 * f - 2 : 3 * f) = cell_face_normals(:, l, c)
          end if
       end do
    end do

    ! one interpolation weight per face: the tail cell's share
    allocate(weights(num_faces))
    do f = 1, num_faces
       weights(f) = face_cell_weights(1, f)
    end do

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

    m = mesh(num_cells, tails=tails, heads=heads, &
         & volumes      = cell_volumes, &
         & cell_centers = reshape(cell_centers, [3 * num_cells]), &
         & areas        = face_areas, &
         & deltas       = face_deltas, &
         & normals      = normals, &
         & face_centers = reshape(face_centers, [3 * num_faces]), &
         & weights      = weights, &
         & etags        = etags)

  end function mesh_from_gmsh

end module class_mesh_builder
