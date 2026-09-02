!=====================================================================!
! MESH DERIVATION: the mesh expressed in the eight primes.
!
! A mesh file supplies three primitive pieces of data:
!
!      member sets    vertices, cells, and the tagged boundary
!                     faces                        (prime: graph)
!      one relation   C2V subset of cells x vertices
!      one field      coordinates : vertices -> R^3
!      labels         physical tag names           (prime: map)
!
! Everything else a finite-volume solver reads is DERIVED here by
! relation algebra and geometry, as pure functions of those inputs:
!
!      faces          the interior faces are the shared-vertex
!                     intersections of cell pairs; with the file's
!                     boundary faces they form the face member set
!                     and the relation F2V subset of faces x vertices
!      V2C            transpose(C2V), by relation_binary's
!                     transpose_padded
!      F2C            {(f, c) : F2V(f) subset of C2V(c)}, at most
!                     two cells per face; one cell = boundary face
!      C2F            transpose(F2C)
!      measurements   fields from coordinates and the relations:
!                     cell centres and volumes, face centres, areas,
!                     normals, centre-to-centre vectors, deltas,
!                     and interpolation weights
!      directed view  cells as graph vertices, the two-cell faces
!                     as edges with tail/head = F2C(f); assembled
!                     by view_mesh_builder into view_mesh
!
! Every routine takes plain arrays and returns plain arrays; no
! graph type is needed to measure geometry. Incidence tables are
! padded rectangles with per-entry counts, entries numbered from 1.
!
! Checks that stop the program: a face vertex count outside 2..4,
! a face incident to no cell, a nonpositive face area, and a negative
! cell volume (an inside-out cell).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_geometry

  use util_precision  , only : dp
  use view_mesh       , only : mesh
  use relation_binary , only : transpose_padded

  implicit none

  private
  public :: derive_faces
  public :: derive_face_cells
  public :: derive_cell_centres
  public :: derive_face_vectors
  public :: outward_sign
  public :: derive_cell_volumes
  public :: derive_centroidal_vectors
  public :: derive_face_deltas
  public :: derive_face_weights
  public :: mesh_from_incidence
  public :: find
  public :: element_dimension
  public :: element_num_vertices
  public :: element_kind, elements, gmsh_kinds, widest_element, face_kind

  !===================================================================!
  ! THE ELEMENT TABLE: the meaning of each gmsh element number, stated
  ! once. Each kind stores its dimension, its vertex count, its faces
  ! - which of its vertices form each face, in the order that makes
  ! the face normal outward - and the paraview type that renders it.
  ! Every query about an element kind is a read of this table: the
  ! loader's widths, the face count the algebraic face total needs,
  ! the ordering of a shared face, the writer's cell type. A gmsh
  ! number absent from the table has dimension -1 and no other data.
  !===================================================================!

  type :: element_kind
     integer :: dimension    = -1
     integer :: num_vertices = 0
     integer :: num_faces    = 0
     integer :: vtk_type     = 0
     integer :: num_face_vertices(6) = 0
     integer :: face_vertices(4, 6)  = 0
  end type element_kind

  integer, parameter :: gmsh_kinds = 15

  type(element_kind), parameter :: elements(gmsh_kinds) = [ &
       ! 1: a 2-node line; its faces are its two vertices
       & element_kind(1, 2, 2, 3, [1, 1, 0, 0, 0, 0], &
       &   reshape([1,0,0,0, 2,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0], [4, 6])), &
       ! 2: a 3-node triangle; its faces are its three edges
       & element_kind(2, 3, 3, 5, [2, 2, 2, 0, 0, 0], &
       &   reshape([1,2,0,0, 2,3,0,0, 3,1,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0], [4, 6])), &
       ! 3: a 4-node quadrangle; its faces are its four edges
       & element_kind(2, 4, 4, 9, [2, 2, 2, 2, 0, 0], &
       &   reshape([1,2,0,0, 2,3,0,0, 3,4,0,0, 4,1,0,0, 0,0,0,0, 0,0,0,0], [4, 6])), &
       ! 4: a 4-node tetrahedron
       & element_kind(3, 4, 4, 10, [3, 3, 3, 3, 0, 0], &
       &   reshape([1,2,3,0, 3,2,4,0, 1,3,4,0, 1,4,2,0, 0,0,0,0, 0,0,0,0], [4, 6])), &
       ! 5: an 8-node hexahedron
       & element_kind(3, 8, 6, 12, [4, 4, 4, 4, 4, 4], &
       &   reshape([1,5,8,4, 2,3,7,6, 5,6,7,8, 1,4,3,2, 7,3,4,8, 1,2,6,5], [4, 6])), &
       ! 6: a 6-node prism, three quadrangles and two triangles
       & element_kind(3, 6, 5, 13, [4, 4, 4, 3, 3, 0], &
       &   reshape([1,4,5,6, 1,3,6,4, 2,5,6,3, 1,2,3,0, 4,6,5,0, 0,0,0,0], [4, 6])), &
       ! 7: a 5-node pyramid, one quadrangle and four triangles
       & element_kind(3, 5, 5, 14, [4, 3, 3, 3, 3, 0], &
       &   reshape([1,2,3,4, 1,2,5,0, 2,3,5,0, 3,4,5,0, 1,5,4,0, 0,0,0,0], [4, 6])), &
       ! 8 to 14: second-order elements, which this table does not define
       & element_kind(), element_kind(), element_kind(), element_kind(), &
       & element_kind(), element_kind(), element_kind(), &
       ! 15: a 1-node point
       & element_kind(0, 1, 0, 1, [0, 0, 0, 0, 0, 0], reshape([0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0], [4, 6])) ]


contains

  !===================================================================!
  ! The face member set and F2V. The file supplies only the tagged
  ! boundary faces; the interior faces are derived: two cells share
  ! an interior face when they share 2 vertices in 2d (a segment),
  ! or 3 or 4 vertices in 3d (a triangle or quadrilateral). The
  ! total count is computed algebraically before the search,
  !
  !      num_faces = (sum of per-cell face counts - nbfaces)/2
  !                  + nbfaces,
  !
  ! because every interior face is counted by exactly two cells.
  ! Boundary faces retain the file's order and tags, positions
  ! 1..nbfaces; interior faces follow in cell-pair order. A shared
  ! vertex count above 4 stops the program.
  !===================================================================!

  impure subroutine derive_faces(spatial_dim, num_points, num_cells, &
       & cell_vertices, num_cell_vertices, cell_types, cell_tags, &
       & nbfaces, bface_vertices, bnum_face_vertices, bface_tags, &
       & bface_types, num_faces, face_vertices, num_face_vertices, &
       & face_tags, face_types)

    integer, intent(in)  :: spatial_dim
    integer, intent(in)  :: num_points
    integer, intent(in)  :: num_cells
    integer, intent(in)  :: cell_vertices(:,:)
    integer, intent(in)  :: num_cell_vertices(:)
    integer, intent(in)  :: cell_types(:)
    integer, intent(in)  :: cell_tags(:)
    integer, intent(in)  :: nbfaces
    integer, intent(in)  :: bface_vertices(:,:)
    integer, intent(in)  :: bnum_face_vertices(:)
    integer, intent(in)  :: bface_tags(:)
    integer, intent(in)  :: bface_types(:)
    integer, intent(out) :: num_faces
    integer, allocatable, intent(out) :: face_vertices(:,:)
    integer, allocatable, intent(out) :: num_face_vertices(:)
    integer, allocatable, intent(out) :: face_tags(:)
    integer, allocatable, intent(out) :: face_types(:)

    integer :: iface, icell, jcell, ivertex
    integer :: shared_node_count, shared_kind
    integer, allocatable :: shared_vertices(:)
    integer :: nf

    num_faces = (sum(element_num_faces(cell_types)) - nbfaces)/2 + nbfaces

    allocate(face_vertices(widest_element(spatial_dim - 1), num_faces))
    allocate(shared_vertices(maxval(num_cell_vertices)))
    allocate(num_face_vertices(num_faces))
    allocate(face_tags(num_faces))
    allocate(face_types(num_faces))
    face_vertices     = 0
    num_face_vertices = 0
    face_tags         = 0
    face_types        = 0

    ! the file's boundary faces first, order and tags preserved
    do iface = 1, nbfaces
       face_tags(iface)         = bface_tags(iface)
       face_types(iface)        = bface_types(iface)
       num_face_vertices(iface) = bnum_face_vertices(iface)
       face_vertices(1:bnum_face_vertices(iface), iface) = &
            & bface_vertices(1:bnum_face_vertices(iface), iface)
    end do

    ! interior faces from shared vertices, upper-triangle cell pairs
    nf = nbfaces

    do icell = 1, num_cells

       shared_vertices = 0

       do jcell = icell + 1, num_cells

          shared_node_count = 0
          do ivertex = 1, num_cell_vertices(icell)
             if (any(cell_vertices(1:num_cell_vertices(jcell), jcell) .eq. &
                  & cell_vertices(ivertex, icell))) then
                shared_node_count = shared_node_count + 1
                shared_vertices(shared_node_count) = &
                     & cell_vertices(ivertex, icell)
             end if
          end do

          ! a shared vertex is never a face, and in 3d neither is a
          ! shared segment: a face has as many vertices as the space
          ! has dimensions, at the least, and is a kind the table defines
          if (shared_node_count >= spatial_dim) then
             shared_kind = face_kind(spatial_dim - 1, shared_node_count)
             if (shared_kind == 0) then
                error stop 'view_mesh_geometry: two cells share the vertices of a face the table defines'
             end if
             nf = nf + 1
             face_tags(nf)         = cell_tags(icell)
             face_types(nf)        = shared_kind
             num_face_vertices(nf) = shared_node_count
             call order_face_vertices(cell_types(jcell), cell_vertices(:, jcell), &
                  & shared_vertices(1:shared_node_count))
             face_vertices(1:shared_node_count, nf) = shared_vertices(1:shared_node_count)
          end if


       end do

    end do

    if (nf .ne. num_faces) then
       error stop 'view_mesh_geometry: the derived face count matches the algebraic count'
    end if

    if (maxval(face_vertices) .ne. num_points) then
       error stop 'view_mesh_geometry: every vertex belongs to a face'
    end if

  end subroutine derive_faces

  !===================================================================!
  ! F2C: a cell contains a face when the cell contains every one of
  ! the face's vertices; the candidates are the cells incident to the
  ! face's first vertex (V2C supplies them), so the search is
  ! local. At most two cells contain a face; exactly one means a
  ! boundary face. A face with no cell stops the program.
  !===================================================================!

  pure subroutine derive_face_cells(num_faces, face_vertices, &
       & num_face_vertices, cell_vertices, num_cell_vertices, &
       & vertex_cells, num_vertex_cells, face_cells, num_face_cells)

    integer, intent(in)  :: num_faces
    integer, intent(in)  :: face_vertices(:,:)
    integer, intent(in)  :: num_face_vertices(:)
    integer, intent(in)  :: cell_vertices(:,:)
    integer, intent(in)  :: num_cell_vertices(:)
    integer, intent(in)  :: vertex_cells(:,:)
    integer, intent(in)  :: num_vertex_cells(:)
    integer, allocatable, intent(out) :: face_cells(:,:)
    integer, allocatable, intent(out) :: num_face_cells(:)

    integer :: icell, iface, k, iv
    logical :: contains_all

    allocate(num_face_cells(num_faces))
    allocate(face_cells(2, num_faces))
    num_face_cells = 0
    face_cells     = 0

    do iface = 1, num_faces
       do k = 1, num_vertex_cells(face_vertices(1, iface))
          icell = vertex_cells(k, face_vertices(1, iface))
          contains_all = .true.
          do iv = 1, num_face_vertices(iface)
             if (.not. any(cell_vertices(1:num_cell_vertices(icell), icell) &
                  &        .eq. face_vertices(iv, iface))) then
                contains_all = .false.
                exit
             end if
          end do
          if (.not. contains_all) cycle
          num_face_cells(iface) = num_face_cells(iface) + 1
          face_cells(num_face_cells(iface), iface) = icell
       end do
    end do

    if (minval(num_face_cells) .lt. 1) then
       error stop 'view_mesh_geometry: every face is incident to a cell'
    end if

  end subroutine derive_face_cells

  !===================================================================!
  ! Cell centroids by the divergence theorem: since div(x_i x) is
  ! (d + 1) x_i, and x . n is constant on a flat face,
  !
  !      x_c = sum_f sigma_cf (x_f . S_f) x_f / ((d + 1) V)
  !
  ! with sigma_cf the sign that turns S_f out of the cell and x_f the
  ! face's centroid - the dual of the volume sum, and exact for any
  ! cell with flat faces. It is the true centroid, about which the
  ! first moment vanishes, so a cell average evaluated there is second
  ! order; the vertex mean is not, except by symmetry.
  !===================================================================!

  pure subroutine derive_cell_centres(spatial_dim, face_centres, face_vectors, &
       & interior, cell_faces, num_cell_faces, cell_volumes, cell_centres)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: face_centres(:,:)
    real(dp), intent(in)  :: face_vectors(:,:)
    real(dp), intent(in)  :: interior(:,:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), intent(in)  :: cell_volumes(:)
    real(dp), allocatable, intent(out) :: cell_centres(:,:)

    integer :: lcell, lface, gface

    allocate(cell_centres(size(face_centres, 1), size(num_cell_faces)))
    cell_centres = 0.0_dp

    ! the sign orients each face's area vector out of this cell,
    ! evaluated against any interior point; the vertex mean is one
    ! such point
    do lcell = 1, size(num_cell_faces)
       do lface = 1, num_cell_faces(lcell)
          gface = cell_faces(lface, lcell)
          cell_centres(:, lcell) = cell_centres(:, lcell) &
               & + outward_sign(face_vectors(:, gface), face_centres(:, gface), interior(:, lcell)) &
               & * dot_product(face_vectors(:, gface), face_centres(:, gface)) &
               & * face_centres(:, gface)
       end do
       cell_centres(:, lcell) = cell_centres(:, lcell) &
            & / (real(spatial_dim + 1, dp) * cell_volumes(lcell))
    end do

  end subroutine derive_cell_centres

  !===================================================================!
  ! Face centres and area vectors. The one measurement of a face is
  ! its area vector S_f: the corner fan from the first vertex, each
  ! simplex of d-1 edge vectors taken through the dual and summed,
  ! over (d-1)!. In 2d that is the segment turned a quarter turn in
  ! plane, (t_y, -t_x, 0), so |S_f| is the length; in 3d half the
  ! sum of the fan's cross products, one triangle plus a second for
  ! a quadrilateral, so |S_f| is the vector area - for a planar face
  ! the scalar area, for a non-planar quadrilateral the area the
  ! divergence theorem integrates. Neither case is written out here:
  ! one formula serves every dimension. The scalar area is its norm;
  ! the unit normal is S_f over its norm and is computed where it is
  ! needed. The centre is the vertex mean. A face whose area vector
  ! vanishes - coincident points, a degenerate fan - stops the program
  ! before any division by it.
  !===================================================================!

  pure subroutine derive_face_vectors(spatial_dim, coordinates, &
       & face_vertices, num_face_vertices, face_centres, face_vectors, &
       & face_areas)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: coordinates(:,:)
    integer , intent(in)  :: face_vertices(:,:)
    integer , intent(in)  :: num_face_vertices(:)
    real(dp), allocatable, intent(out) :: face_centres(:,:)
    real(dp), allocatable, intent(out) :: face_vectors(:,:)
    real(dp), allocatable, intent(out) :: face_areas(:)

    integer  :: iface, num_faces, i, k, rows
    real(dp) :: s(size(coordinates, 1)), piece(size(coordinates, 1))
    real(dp) :: fan(size(coordinates, 1), spatial_dim - 1), mean(size(coordinates, 1))
    real(dp) :: measure, whole

    rows      = size(coordinates, 1)
    num_faces = size(num_face_vertices)
    allocate(face_centres(rows, num_faces))
    allocate(face_vectors(rows, num_faces))
    allocate(face_areas(num_faces))
    face_centres = 0.0_dp
    face_vectors = 0.0_dp

    do iface = 1, num_faces
       associate(facenodes => face_vertices(1:num_face_vertices(iface), iface))

         ! the fan from the first vertex: each simplex's dual summed is
         ! the area vector, and each simplex's mean weighted by its
         ! measure is the face's centroid - exact for any flat face,
         ! the midpoint of a segment, the vertex mean of a triangle
         s     = 0.0_dp
         mean  = 0.0_dp
         whole = 0.0_dp
         do i = 2, num_face_vertices(iface) - spatial_dim + 2
            do k = 1, spatial_dim - 1
               fan(:, k) = coordinates(:, facenodes(i + k - 1)) - coordinates(:, facenodes(1))
            end do
            piece   = dual(spatial_dim, fan)
            measure = norm2(piece)
            s       = s + piece
            mean    = mean + measure * (coordinates(:, facenodes(1)) &
                 & + sum(fan, dim=2) / real(spatial_dim, dp))
            whole   = whole + measure
         end do
         s = s / real(factorial(spatial_dim - 1), dp)

         face_vectors(:, iface) = s
         face_areas(iface)      = norm2(s)
         if (whole > 0.0_dp) then
            face_centres(:, iface) = mean / whole
         else
            face_centres(:, iface) = sum(coordinates(:, facenodes), dim=2) &
                 & / real(num_face_vertices(iface), kind=dp)
         end if

       end associate
    end do

    ! the lower bound is the build's own real kind's, not one fixed
    ! kind's: a quadruple-precision build has a far smaller tiny()
    ! than a double-precision build, and the check scales with it
    if (minval(face_areas) < 10.0_dp * tiny(1.0_dp)) then
       error stop 'view_mesh_geometry: a face has a nonzero area vector'
    end if

  end subroutine derive_face_vectors

  !===================================================================!
  ! The sign that points a face's area vector out of a cell: +1 when
  ! S_f points from the cell centre toward the face centre, -1 when
  ! S_f must be reversed. The same interior face has opposite signs
  ! relative to its two cells, which makes the assembled operator
  ! conservative.
  !===================================================================!

  pure real(dp) function outward_sign(face_vector, face_centre, cell_centre)

    real(dp), intent(in) :: face_vector(:), face_centre(:), cell_centre(:)

    outward_sign = 1.0_dp
    if (dot_product(face_vector, face_centre - cell_centre) .lt. 0.0_dp) then
       outward_sign = -1.0_dp
    end if

  end function outward_sign


  !===================================================================!
  ! Cell volumes by the divergence theorem,
  !
  !      V = (1/d) sum over the cell's faces of  sigma_cf (S_f . x_f),
  !
  ! with d the spatial dimension, S_f the face's area vector and
  ! sigma_cf the sign that points it out of the cell. A negative
  ! volume is an inside-out cell and stops the program.
  !===================================================================!

  pure subroutine derive_cell_volumes(spatial_dim, face_centres, face_vectors, &
       & cell_centres, cell_faces, num_cell_faces, cell_volumes)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: face_centres(:,:)
    real(dp), intent(in)  :: face_vectors(:,:)
    real(dp), intent(in)  :: cell_centres(:,:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), allocatable, intent(out) :: cell_volumes(:)

    integer :: lcell, lface, gface

    allocate(cell_volumes(size(num_cell_faces)))
    cell_volumes = 0.0_dp

    do lcell = 1, size(num_cell_faces)
       do lface = 1, num_cell_faces(lcell)
          gface = cell_faces(lface, lcell)
          cell_volumes(lcell) = cell_volumes(lcell) + &
               & outward_sign(face_vectors(:, gface), face_centres(:, gface), &
               &              cell_centres(:, lcell)) &
               & * dot_product(face_vectors(:, gface), face_centres(:, gface)) &
               & / real(spatial_dim, dp)
       end do
    end do

    if (minval(cell_volumes) .lt. 0.0_dp) then
       error stop 'view_mesh_geometry: a cell volume is nonnegative'
    end if

  end subroutine derive_cell_volumes

  !===================================================================!
  ! The centre-to-centre vector of each face: cell centre to cell
  ! centre across an interior face, cell centre to face centre on
  ! a boundary face. This is the segment a two-point gradient
  ! differences along.
  !===================================================================!

  pure subroutine derive_centroidal_vectors(face_cells, num_face_cells, &
       & cell_centres, face_centres, lvec)

    integer , intent(in)  :: face_cells(:,:)
    integer , intent(in)  :: num_face_cells(:)
    real(dp), intent(in)  :: cell_centres(:,:)
    real(dp), intent(in)  :: face_centres(:,:)
    real(dp), allocatable, intent(out) :: lvec(:,:)

    integer :: iface

    allocate(lvec(size(face_centres, 1), size(num_face_cells)))
    lvec = real(0, dp)

    do iface = 1, size(num_face_cells)
       if (num_face_cells(iface) .eq. 1) then
          lvec(:, iface) = face_centres(:, iface) &
               & - cell_centres(:, face_cells(1, iface))
       else
          lvec(:, iface) = cell_centres(:, face_cells(2, iface)) &
               & - cell_centres(:, face_cells(1, iface))
       end if
    end do

  end subroutine derive_centroidal_vectors

  !===================================================================!
  ! Face deltas: the centre-to-centre vector projected on the unit
  ! normal, |l_f . S_f| / |S_f|, one per face - the sign of the
  ! normal cancels under the modulus. On a skewed mesh the segment
  ! crosses its face obliquely, so this normal distance is shorter
  ! than the segment; it is the denominator of every two-point face
  ! gradient.
  !===================================================================!

  pure subroutine derive_face_deltas(lvec, face_vectors, face_deltas)

    real(dp), intent(in)  :: lvec(:,:)
    real(dp), intent(in)  :: face_vectors(:,:)
    real(dp), allocatable, intent(out) :: face_deltas(:)

    integer :: f

    allocate(face_deltas(size(face_vectors, 2)))

    do f = 1, size(face_vectors, 2)
       face_deltas(f) = abs(dot_product(lvec(:, f), face_vectors(:, f))) &
            & / norm2(face_vectors(:, f))
    end do

  end subroutine derive_face_deltas

  !===================================================================!
  ! Interpolation weights per face: the two cells' fractions by
  ! inverse distance to the face centre, stored as (w, 1-w). A
  ! boundary face has one cell, which receives the whole weight.
  !===================================================================!

  pure subroutine derive_face_weights(face_cells, num_face_cells, &
       & cell_centres, face_centres, face_cell_weights)

    integer , intent(in)  :: face_cells(:,:)
    integer , intent(in)  :: num_face_cells(:)
    real(dp), intent(in)  :: cell_centres(:,:)
    real(dp), intent(in)  :: face_centres(:,:)
    real(dp), allocatable, intent(out) :: face_cell_weights(:,:)

    integer  :: iface
    real(dp) :: d1, d2, dinv1, dinv2, weight

    allocate(face_cell_weights(2, size(num_face_cells)))

    do iface = 1, size(num_face_cells)

       d1    = distance(cell_centres(:, face_cells(1, iface)), &
            &           face_centres(:, iface))
       dinv1 = 1.0_dp/d1

       if (num_face_cells(iface) .ne. 1) then
          d2    = distance(cell_centres(:, face_cells(2, iface)), &
               &           face_centres(:, iface))
          dinv2 = 1.0_dp/d2
       else
          dinv2 = 0.0_dp
       end if

       weight = dinv1/(dinv1 + dinv2)

       face_cell_weights(1:2, iface) = [weight, 1.0_dp - weight]

    end do

  end subroutine derive_face_weights


  !===================================================================!
  ! The dual of d-1 vectors in d dimensions: the vector normal to all
  ! of them whose length is the volume they span - the cross product
  ! at d = 3, the quarter turn (t_y, -t_x) at d = 2, and in general
  ! the cofactors of the matrix whose columns they are, component i
  ! being (-1)^(i+1) times the minor with row i struck. One formula
  ! serves every dimension, so no dimension is written out. The
  ! vectors have as many rows as the coordinates have; only the
  ! first d rows are read, and the components past d are zero.
  !===================================================================!

  pure function dual(d, vectors) result(s)

    integer , intent(in) :: d
    real(dp), intent(in) :: vectors(:,:)
    real(dp) :: s(size(vectors, 1))

    real(dp) :: minor(d - 1, d - 1)
    integer  :: i, k

    s = 0.0_dp
    do i = 1, d
       minor = vectors([(k, k = 1, i - 1), (k, k = i + 1, d)], 1:d - 1)
       s(i) = real((-1) ** (i + 1), dp) * determinant(minor)
    end do

  end function dual

  !===================================================================!
  ! The determinant, by expansion along the first row. The minors
  ! here are at most (d-1) square, so the expansion's operation count
  ! is small and a factorisation is not warranted.
  !===================================================================!

  pure recursive function determinant(a) result(det)

    real(dp), intent(in) :: a(:,:)
    real(dp) :: det

    integer :: n, j, k

    n = size(a, 1)
    if (n == 0) then
       det = 1.0_dp
    else if (n == 1) then
       det = a(1, 1)
    else
       det = 0.0_dp
       do j = 1, n
          det = det + real((-1) ** (j + 1), dp) * a(1, j) &
               & * determinant(a(2:n, [(k, k = 1, j - 1), (k, k = j + 1, n)]))
       end do
    end if

  end function determinant

  pure integer function factorial(n)

    integer, intent(in) :: n

    integer :: k

    factorial = 1
    do k = 2, n
       factorial = factorial * k
    end do

  end function factorial

  !===================================================================!
  ! Compute the geometric distance between two points.
  !===================================================================!

  pure real(dp) function distance(x, y)

    real(dp), intent(in)  :: X(:), y(:) ! The shape is [[x,y,z], [1:2]].

    distance = sqrt(sum((x-y)**2))

  end function distance

  !===================================================================!
  ! THE MESH FROM ITS INCIDENCES. Every measurement a mesh stores is
  ! a function of the coordinates and three relations - each cell's
  ! vertices, each face's vertices, each face's cells - and of nothing
  ! else: not of element types, not of how the faces were found. So
  ! this is the one terminal step of every mesh construction: the gmsh
  ! builder reaches it after deriving its faces, the spatial level
  ! after enumerating its own, and the measurements are computed once.
  !
  ! The coordinates may have more rows than the space has dimensions
  ! (a two-dimensional mesh read from a three-dimensional file); the
  ! mesh retains the first spatial_dim of them. A face with one cell
  ! is a boundary face, an edge without a head, and stores the tag
  ! given; the unit normal points out of the tail cell; the weight is
  ! the tail cell's interpolation fraction.
  !===================================================================!

  impure type(mesh) function mesh_from_incidence(spatial_dim, coordinates, &
       & cell_vertices, num_cell_vertices, face_vertices, num_face_vertices, &
       & face_cells, num_face_cells, etags) result(m)

    integer         , intent(in) :: spatial_dim
    real(dp)        , intent(in) :: coordinates(:,:)
    integer         , intent(in) :: cell_vertices(:,:)
    integer         , intent(in) :: num_cell_vertices(:)
    integer         , intent(in) :: face_vertices(:,:)
    integer         , intent(in) :: num_face_vertices(:)
    integer         , intent(in) :: face_cells(:,:)
    integer         , intent(in) :: num_face_cells(:)
    character(len=*), intent(in), optional :: etags(:)

    integer , allocatable :: cell_faces(:,:), num_cell_faces(:), tails(:), heads(:)
    real(dp), allocatable :: interior(:,:), cell_centres(:,:), face_centres(:,:)
    real(dp), allocatable :: face_vectors(:,:), face_areas(:), cell_volumes(:), lvec(:,:)
    real(dp), allocatable :: face_deltas(:), face_cell_weights(:,:), normals(:), weights(:)
    integer :: num_cells, num_faces, c, f, d

    d         = spatial_dim
    num_cells = size(num_cell_vertices)
    num_faces = size(num_face_vertices)

    if (size(coordinates, 1) < d) then
       error stop 'view_mesh_geometry: the coordinates have a row for every dimension of the space'
    end if

    call transpose_padded(face_cells, num_face_cells, num_cells, cell_faces, num_cell_faces)

    call derive_face_vectors(d, coordinates, face_vertices, num_face_vertices, &
         & face_centres, face_vectors, face_areas)

    ! an interior point of each cell for the outward signs: the mean
    ! of its vertices, inside any convex cell
    allocate(interior(size(coordinates, 1), num_cells))
    do c = 1, num_cells
       interior(:, c) = sum(coordinates(:, cell_vertices(1:num_cell_vertices(c), c)), dim=2) &
            & / real(num_cell_vertices(c), dp)
    end do

    call derive_cell_volumes(d, face_centres, face_vectors, interior, cell_faces, &
         & num_cell_faces, cell_volumes)
    call derive_cell_centres(d, face_centres, face_vectors, interior, cell_faces, &
         & num_cell_faces, cell_volumes, cell_centres)
    call derive_centroidal_vectors(face_cells, num_face_cells, cell_centres, face_centres, lvec)
    call derive_face_deltas(lvec, face_vectors, face_deltas)
    call derive_face_weights(face_cells, num_face_cells, cell_centres, face_centres, &
         & face_cell_weights)

    allocate(tails(num_faces), heads(num_faces), normals(d * num_faces), weights(num_faces))
    do f = 1, num_faces
       tails(f) = face_cells(1, f)
       heads(f) = 0
       if (num_face_cells(f) >= 2) heads(f) = face_cells(2, f)
       normals(d * (f - 1) + 1 : d * f) = &
            & outward_sign(face_vectors(:, f), face_centres(:, f), cell_centres(:, tails(f))) &
            & * face_vectors(1:d, f) / face_areas(f)
       weights(f) = face_cell_weights(1, f)
    end do

    m = mesh(num_cells, tails=tails, heads=heads, &
         & volumes      = cell_volumes, &
         & cell_centres = reshape(cell_centres(1:d, :), [d * num_cells]), &
         & areas        = face_areas, &
         & deltas       = face_deltas, &
         & normals      = normals, &
         & face_centres = reshape(face_centres(1:d, :), [d * num_faces]), &
         & weights      = weights, &
         & etags        = etags, &
         & dimension    = d)

  end function mesh_from_incidence

  !===================================================================!
  ! Return the index of a target value if it is present in the array;
  ! return -1 otherwise.
  !===================================================================!

  pure type(integer) function find(array, target_value)

    integer, intent(in) :: array(:)
    integer, intent(in) :: target_value
    integer :: i, num_entries

    num_entries = size(array, dim=1)

    do i = 1, num_entries
       if (array(i) .eq. target_value) then
          find = i
          return
       endif
    end do

    find = -1

  end function find

  !===================================================================!
  ! Return the number of faces a gmsh cell type has. The name could
  ! generalize to the number of lower dimensional entities.
  !===================================================================!


  !===================================================================!
  ! Return the spatial dimension of a gmsh element type: a point is
  ! 0, a line is 1, a triangle or quadrangle is 2, and a tet, hex,
  ! prism, or pyramid is 3. The mesh classifies elements into cells,
  ! faces, and edges by dimension rather than by type.
  !===================================================================!


  !===================================================================!
  ! Return the number of vertices a gmsh element type has. The name
  ! could generalize to the number of lower dimensional entities.
  !===================================================================!


  !===================================================================!
  ! Put a face's vertices into the cell's own vertex order. Each
  ! cell type stores a face table that specifies which vertices bound
  ! which face, ordered so that the normal points outward. Match the
  ! unordered face against the table and return the face in the
  ! table's order.
  !===================================================================!



  !===================================================================!
  ! Three elemental reads of the table; a number outside the table
  ! reads as dimension -1, zero vertices, zero faces.
  !===================================================================!

  pure elemental integer function element_dimension(elem_type) result (dim)

    integer, intent(in) :: elem_type

    dim = -1
    if (elem_type >= 1 .and. elem_type <= gmsh_kinds) dim = elements(elem_type) % dimension

  end function element_dimension

  pure elemental integer function element_num_vertices(elem_type) result (num_vertices)

    integer, intent(in) :: elem_type

    num_vertices = 0
    if (elem_type >= 1 .and. elem_type <= gmsh_kinds) num_vertices = elements(elem_type) % num_vertices

  end function element_num_vertices

  pure elemental integer function element_num_faces(elem_type) result (num_faces)

    integer, intent(in) :: elem_type

    num_faces = 0
    if (elem_type >= 1 .and. elem_type <= gmsh_kinds) num_faces = elements(elem_type) % num_faces

  end function element_num_faces

  !===================================================================!
  ! The most vertices any element of a dimension has - the width a
  ! padded list of such elements needs - and the kind of a face of a
  ! dimension with this many vertices, or zero where the table
  ! defines none.
  !===================================================================!

  pure integer function widest_element(dimension)

    integer, intent(in) :: dimension

    widest_element = max(0, maxval(elements % num_vertices, mask=elements % dimension == dimension))

  end function widest_element

  pure integer function face_kind(dimension, num_vertices)

    integer, intent(in) :: dimension, num_vertices

    integer :: t

    face_kind = 0
    do t = 1, gmsh_kinds
       if (elements(t) % dimension == dimension .and. elements(t) % num_vertices == num_vertices) then
          face_kind = t
          return
       end if
    end do

  end function face_kind

  !===================================================================!
  ! The vertices of a face a cell shares, put in the order the table
  ! gives that face of that cell, which is the order that orients its
  ! area vector outward. The face is found among the cell's faces by
  ! its vertex set; a vertex set that is no face of the cell stops the
  ! program, because the mesh is then nonconforming.
  !===================================================================!

  impure subroutine order_face_vertices(cell_type, cell_vertices, face_vertices_unordered)

    integer, intent(in)    :: cell_type
    integer, intent(in)    :: cell_vertices(:)
    integer, intent(inout) :: face_vertices_unordered(:)

    type(element_kind) :: cell_shape
    integer :: iface, ivertex, n, match_count

    if (element_num_faces(cell_type) == 0) then
       error stop 'view_mesh_geometry: the element table lists the faces of this cell type'
    end if

    cell_shape = elements(cell_type)

    do iface = 1, cell_shape % num_faces
       n = cell_shape % num_face_vertices(iface)
       if (n /= size(face_vertices_unordered)) cycle
       match_count = 0
       do ivertex = 1, n
          if (any(face_vertices_unordered == cell_vertices(cell_shape % face_vertices(ivertex, iface)))) then
             match_count = match_count + 1
          end if
       end do
       if (match_count == n) then
          face_vertices_unordered = cell_vertices(cell_shape % face_vertices(1:n, iface))
          return
       end if
    end do

    error stop 'view_mesh_geometry: the vertices two cells share are a face of the cell'

  end subroutine order_face_vertices

end module view_mesh_geometry
