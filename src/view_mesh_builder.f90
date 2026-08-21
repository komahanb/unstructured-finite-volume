!=====================================================================!
! The mesh builder: a gmsh file becomes a mesh in one call,
!
!      m = mesh_from_gmsh('square-10.msh')
!
! The pipeline is the prime decomposition stated in
! view_mesh_geometry: the loader parses the file into member sets,
! the cell-to-vertex relation, the vertex coordinates, and the tag
! names; the geometry module derives the face set, the incidence
! relations, and every measurement; this builder places the results
! as the mesh - the directed view whose vertices are cells and
! whose edges are the two-cell faces, with the measurements as
! fields and the tag names on the boundary edges.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_builder

  use iso_fortran_env      , only : dp => REAL64
  use view_mesh     , only : mesh
  use util_string         , only : string
  use view_mesh_loader, only : mesh_loader
  use view_gmsh_loader    , only : gmsh_loader
  use view_mesh_geometry  , only : element_dimension
  use relation_binary, only : transpose_padded
  use view_mesh_geometry  , only : derive_faces, derive_face_cells, &
       & derive_cell_centres, &
       & derive_face_vectors, outward_sign, derive_cell_volumes, derive_centroidal_vectors, &
       & derive_face_deltas, derive_face_weights

  implicit none

  private
  public :: mesh_from_gmsh

contains

  !===================================================================!
  ! Load, derive, place. Tails and heads come from the face-to-cell
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
    real(dp), allocatable :: cell_centres(:,:), face_centres(:,:)
    real(dp), allocatable :: face_areas(:), cell_volumes(:)
    real(dp), allocatable :: face_vectors(:,:), lvec(:,:)
    real(dp), allocatable :: face_deltas(:), face_cell_weights(:,:)

    ! the view's per-face arrays
    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: normals(:)
    real(dp), allocatable :: weights(:)
    character(len=64), allocatable :: etags(:)
    integer :: f, tag

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

    call transpose_padded(face_cells, num_face_cells, num_cells, &
         & cell_faces, num_cell_faces)

    ! the measurements
    call derive_cell_centres(vertices, cell_vertices, num_cell_vertices, &
         & cell_centres)
    call derive_face_vectors(spatial_dim, vertices, face_vertices, &
         & num_face_vertices, face_centres, face_vectors, face_areas)
    call derive_cell_volumes(spatial_dim, face_centres, face_vectors, &
         & cell_centres, cell_faces, num_cell_faces, cell_volumes)
    call derive_centroidal_vectors(face_cells, num_face_cells, cell_centres, &
         & face_centres, lvec)
    call derive_face_deltas(lvec, face_vectors, face_deltas)
    call derive_face_weights(face_cells, num_face_cells, cell_centres, &
         & face_centres, face_cell_weights)

    ! tails and heads from the face-to-cell relation; a face with
    ! one cell is a boundary face, an edge without a head
    allocate(tails(num_faces), heads(num_faces))
    do f = 1, num_faces
       tails(f) = face_cells(1, f)
       heads(f) = 0
       if (num_face_cells(f) >= 2) heads(f) = face_cells(2, f)
    end do

    ! one unit normal per face, pointing out of its tail cell
    allocate(normals(3 * num_faces))
    do f = 1, num_faces
       normals(3 * f - 2 : 3 * f) = &
            & outward_sign(face_vectors(:, f), face_centres(:, f), &
            &              cell_centres(:, face_cells(1, f))) &
            & * face_vectors(:, f) / face_areas(f)
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
         & cell_centres = reshape(cell_centres, [3 * num_cells]), &
         & areas        = face_areas, &
         & deltas       = face_deltas, &
         & normals      = normals, &
         & face_centres = reshape(face_centres, [3 * num_faces]), &
         & weights      = weights, &
         & etags        = etags)

  end function mesh_from_gmsh

end module view_mesh_builder
