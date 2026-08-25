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
! a face touching no cell, a nonpositive face area, and a negative
! cell volume (an inside-out cell).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_geometry

  use util_precision  , only : dp

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
  ! THE ELEMENT TABLE: what gmsh's element numbers mean, once. Each
  ! kind carries its dimension, its vertex count, its faces - which of
  ! its vertices make each face, in the order that turns the face
  ! outward - and the paraview type that draws it. Every question
  ! about an element kind is a read of this table: the loader's
  ! widths, the face count the algebraic face total needs, the
  ! ordering of a shared face, the writer's cell type. A gmsh number
  ! the table does not carry has dimension -1 and nothing else.
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
       ! 8 to 14: second-order elements, which this tower does not carry
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
  ! total count is known algebraically before the search,
  !
  !      num_faces = (sum of per-cell face counts - nbfaces)/2
  !                  + nbfaces,
  !
  ! because every interior face is counted by exactly two cells.
  ! Boundary faces keep the file's order and tags, positions
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
          ! has dimensions, at the least, and is a kind the table has
          if (shared_node_count >= spatial_dim) then
             shared_kind = face_kind(spatial_dim - 1, shared_node_count)
             if (shared_kind == 0) then
                error stop 'view_mesh_geometry: two cells share the vertices of a face the table has'
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
       error stop 'view_mesh_geometry: every last vertex belongs to a face'
    end if

  end subroutine derive_faces

  !===================================================================!
  ! F2C: a cell contains a face when it holds every one of the
  ! face's vertices; the candidates are the cells touching the
  ! face's first vertex (V2C supplies them), so the search is
  ! local. At most two cells hold a face; exactly one means a
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
    logical :: holds_all

    allocate(num_face_cells(num_faces))
    allocate(face_cells(2, num_faces))
    num_face_cells = 0
    face_cells     = 0

    do iface = 1, num_faces
       do k = 1, num_vertex_cells(face_vertices(1, iface))
          icell = vertex_cells(k, face_vertices(1, iface))
          holds_all = .true.
          do iv = 1, num_face_vertices(iface)
             if (.not. any(cell_vertices(1:num_cell_vertices(icell), icell) &
                  &        .eq. face_vertices(iv, iface))) then
                holds_all = .false.
                exit
             end if
          end do
          if (.not. holds_all) cycle
          num_face_cells(iface) = num_face_cells(iface) + 1
          face_cells(num_face_cells(iface), iface) = icell
       end do
    end do

    if (minval(num_face_cells) .lt. 1) then
       error stop 'view_mesh_geometry: every face touches a cell'
    end if

  end subroutine derive_face_cells

  !===================================================================!
  ! Cell centroids by the divergence theorem: x_c = sum over the
  ! cell's faces of (x_f . S_f) x_f / ((d+1) V), where x_f is the
  ! face centre, S_f is the face's area vector, and V is the cell
  ! volume. This is the true centroid (first moment vanishes), not
  ! the vertex mean. On a skewed mesh this differs from the vertex
  ! mean; on orthogonal gmsh meshes the two coincide.
  !===================================================================!

  pure subroutine derive_cell_centres(spatial_dim, face_centres, face_vectors, &
       & cell_faces, num_cell_faces, cell_volumes, cell_centres)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: face_centres(:,:)
    real(dp), intent(in)  :: face_vectors(:,:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), intent(in)  :: cell_volumes(:)
    real(dp), allocatable, intent(out) :: cell_centres(:,:)

    integer :: lcell, lface, gface
    real(dp) :: fc_dot_fv, denom

    allocate(cell_centres(3, size(num_cell_faces)))
    cell_centres = 0.0_dp

    denom = real(spatial_dim + 1, dp)

    do lcell = 1, size(num_cell_faces)
       do lface = 1, num_cell_faces(lcell)
          gface = cell_faces(lface, lcell)
          fc_dot_fv = dot_product(face_vectors(:, gface), face_centres(:, gface))
          cell_centres(:, lcell) = cell_centres(:, lcell) + &
               & fc_dot_fv * face_centres(:, gface)
       end do
       if (cell_volumes(lcell) > 0.0_dp) then
          cell_centres(:, lcell) = cell_centres(:, lcell) / (denom * cell_volumes(lcell))
       end if
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
  ! divergence theorem sees. Neither is written here: one formula
  ! serves every dimension. The scalar area is its norm; the unit
  ! normal is S_f over its norm and is taken where it is needed. The
  ! centre is the vertex mean. A face whose area vector vanishes -
  ! coincident points, a degenerate fan - stops the program before
  ! anything divides by it.
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

    integer  :: iface, num_faces, i, k
    real(dp) :: s(3), fan(3, spatial_dim - 1)

    num_faces = size(num_face_vertices)
    allocate(face_centres(3, num_faces))
    allocate(face_vectors(3, num_faces))
    allocate(face_areas(num_faces))
    face_centres = 0.0_dp
    face_vectors = 0.0_dp

    do iface = 1, num_faces
       associate(facenodes => face_vertices(1:num_face_vertices(iface), iface))

         face_centres(:, iface) = sum(coordinates(:, facenodes), dim=2) &
              & / real(num_face_vertices(iface), kind=dp)

         s = 0.0_dp
         do i = 2, num_face_vertices(iface) - spatial_dim + 2
            do k = 1, spatial_dim - 1
               fan(:, k) = coordinates(:, facenodes(i + k - 1)) - coordinates(:, facenodes(1))
            end do
            s = s + dual(spatial_dim, fan)
         end do
         s = s / real(factorial(spatial_dim - 1), dp)

         face_vectors(:, iface) = s
         face_areas(iface)      = norm2(s)

       end associate
    end do

    if (minval(face_areas) .lt. 10.0d0 * tiny(1.0d0)) then
       error stop 'view_mesh_geometry: a face has a nonzero area vector'
    end if

  end subroutine derive_face_vectors

  !===================================================================!
  ! The sign that points a face's area vector out of a cell: +1 when
  ! S_f leaves the cell centre toward the face centre, -1 when it
  ! must be turned. The same interior face carries opposite signs
  ! seen from its two cells, which is what makes the assembled
  ! operator conservative.
  !===================================================================!

  pure real(dp) function outward_sign(face_vector, face_centre, cell_centre)

    real(dp), intent(in) :: face_vector(3), face_centre(3), cell_centre(3)

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

    allocate(lvec(3, size(num_face_cells)))
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
  ! normal drops out under the modulus. On a skewed mesh the segment
  ! crosses its face at a slant, so this normal distance is shorter
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
  ! Interpolation weights per face: the two cells' shares by
  ! inverse distance to the face centre, stored as (w, 1-w). A
  ! boundary face has one cell, which takes the whole weight.
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
  ! vectors arrive as many rows as the coordinates carry; only the
  ! first d are read, and the components past d are zero.
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
  ! here are at most (d-1) square, so the expansion costs nothing
  ! worth a factorisation.
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
  ! Compute mesh measurements from incidence relations and coordinates.
  ! A pure factoring of the tail of mesh_from_gmsh: takes the spatial
  ! dimension, coordinates, and the five incidence relations that
  ! topology dictates, and derives all measurements - cell centres and
  ! volumes, face centres, areas, vectors, centre-to-centre separations,
  ! normal projections, and interpolation weights. Both view_mesh_builder
  ! (after derive_faces) and gti_space (to replace framed()) call this.
  !===================================================================!

  impure subroutine mesh_from_incidence(spatial_dim, coordinates, &
       & cell_vertices, num_cell_vertices, &
       & face_vertices, num_face_vertices, &
       & cell_faces, num_cell_faces, &
       & face_cells, num_face_cells, &
       & cell_centres, face_centres, face_vectors, cell_volumes, &
       & lvec, face_deltas, face_cell_weights)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: coordinates(:,:)
    integer , intent(in)  :: cell_vertices(:,:)
    integer , intent(in)  :: num_cell_vertices(:)
    integer , intent(in)  :: face_vertices(:,:)
    integer , intent(in)  :: num_face_vertices(:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    integer , intent(in)  :: face_cells(:,:)
    integer , intent(in)  :: num_face_cells(:)
    real(dp), allocatable, intent(out) :: cell_centres(:,:)
    real(dp), allocatable, intent(out) :: face_centres(:,:)
    real(dp), allocatable, intent(out) :: face_vectors(:,:)
    real(dp), allocatable, intent(out) :: cell_volumes(:)
    real(dp), allocatable, intent(out) :: lvec(:,:)
    real(dp), allocatable, intent(out) :: face_deltas(:)
    real(dp), allocatable, intent(out) :: face_cell_weights(:,:)

    real(dp), allocatable :: temp_centres(:,:)
    integer :: ncells, c

    ncells = size(num_cell_faces)

    ! Compute face geometry: centres, area vectors (face_vectors), and scalar areas
    call derive_face_vectors(spatial_dim, coordinates, face_vertices, &
         & num_face_vertices, face_centres, face_vectors, face_deltas)

    ! Temporary vertex mean for outward sign checks in volume calculation
    allocate(temp_centres(3, ncells))
    do c = 1, ncells
       associate(vids => cell_vertices(1:num_cell_vertices(c), c))
         temp_centres(:, c) = sum(coordinates(:, vids), dim=2) &
              & / real(num_cell_vertices(c), kind=dp)
       end associate
    end do

    ! Cell volumes using the divergence theorem and the temporary centres
    call derive_cell_volumes(spatial_dim, face_centres, face_vectors, &
         & temp_centres, cell_faces, num_cell_faces, cell_volumes)

    ! True centroids using volumes
    call derive_cell_centres(spatial_dim, face_centres, face_vectors, &
         & cell_faces, num_cell_faces, cell_volumes, cell_centres)

    ! Centre-to-centre vectors for each face
    call derive_centroidal_vectors(face_cells, num_face_cells, &
         & cell_centres, face_centres, lvec)

    ! Recompute face deltas using the true centres (not temporary)
    call derive_face_deltas(lvec, face_vectors, face_deltas)

    ! Interpolation weights for each face
    call derive_face_weights(face_cells, num_face_cells, &
         & cell_centres, face_centres, face_cell_weights)

  end subroutine mesh_from_incidence

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
  ! Return the number of faces a gmsh cell type owns. The name could
  ! generalize to the number of lower dimensional entities.
  !===================================================================!


  !===================================================================!
  ! Return the spatial dimension of a gmsh element type: a point is
  ! 0, a line is 1, a triangle or quadrangle is 2, and a tet, hex,
  ! prism, or pyramid is 3. The mesh classifies elements into cells,
  ! faces, and edges by dimension rather than by type.
  !===================================================================!


  !===================================================================!
  ! Return the number of vertices a gmsh element type owns. The name
  ! could generalize to the number of lower dimensional entities.
  !===================================================================!


  !===================================================================!
  ! Put a face's vertices into the cell's own winding order. Each
  ! cell type carries a wiring table that says which corners bound
  ! which face, wound so that the normal points outward. Match the
  ! unordered face against the table and hand it back in the table's
  ! order.
  !===================================================================!



  !===================================================================!
  ! The table read three ways, elementally, a number outside the table
  ! reading as no dimension, no vertices, no faces.
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
  ! The most vertices any element of a dimension owns - the width a
  ! padded list of such elements needs - and the kind of a face of a
  ! dimension with this many vertices, or zero where the table has
  ! none.
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
  ! gives that face of that cell, which is the order that turns its
  ! area vector outward. The face is found among the cell's faces by
  ! its vertex set; a set that is no face of the cell stops the
  ! program, the mesh being nonconforming.
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
