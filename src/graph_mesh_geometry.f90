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
!                     cell centers and volumes, face centers, areas,
!                     normals, center-to-center vectors, deltas,
!                     and interpolation weights
!      directed view  cells as graph vertices, the two-cell faces
!                     as edges with tail/head = F2C(f); assembled
!                     by class_mesh_builder into class_graph_mesh
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

module graph_mesh_geometry

  use iso_fortran_env  , only : dp => REAL64

  implicit none

  private
  public :: derive_faces
  public :: derive_face_cells
  public :: cell_centers_of
  public :: face_centers_areas_of
  public :: cell_face_normals_of
  public :: cell_volumes_of
  public :: centroidal_vectors_of
  public :: face_deltas_of
  public :: face_weights_of
  public :: find
  public :: elem_type_dimension
  public :: elem_type_vertex_count

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
    integer :: shared_node_count
    integer :: shared_vertices(4)
    integer :: nf

    num_faces = (sum(elem_type_face_count(cell_types)) - nbfaces)/2 + nbfaces

    allocate(face_vertices(4, num_faces))
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

          select case (shared_node_count)
          case (0)
             ! the cells do not touch
          case (1)
             ! a shared vertex is never a face
          case (2)
             ! a shared segment is a face in 2d; in 3d two cells can
             ! touch along a segment without sharing a face
             if (spatial_dim .eq. 2) then
                nf = nf + 1
                face_tags(nf)         = cell_tags(icell)
                face_types(nf)        = 1
                num_face_vertices(nf) = 2
                call order_face_vertices(cell_types(jcell), &
                     & cell_vertices(:, jcell), shared_vertices(1:2))
                face_vertices(1:2, nf) = shared_vertices(1:2)
             end if
          case (3)
             nf = nf + 1
             face_tags(nf)         = cell_tags(icell)
             face_types(nf)        = 2
             num_face_vertices(nf) = 3
             call order_face_vertices(cell_types(jcell), &
                  & cell_vertices(:, jcell), shared_vertices(1:3))
             face_vertices(1:3, nf) = shared_vertices(1:3)
          case (4)
             nf = nf + 1
             face_tags(nf)         = cell_tags(icell)
             face_types(nf)        = 3
             num_face_vertices(nf) = 4
             call order_face_vertices(cell_types(jcell), &
                  & cell_vertices(:, jcell), shared_vertices(1:4))
             face_vertices(1:4, nf) = shared_vertices(1:4)
          case default
             error stop 'graph_mesh_geometry: two cells share at most four vertices'
          end select

       end do

    end do

    if (nf .ne. num_faces) then
       error stop 'graph_mesh_geometry: the derived face count matches the algebraic count'
    end if

    if (maxval(face_vertices) .ne. num_points) then
       error stop 'graph_mesh_geometry: every last vertex belongs to a face'
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
       error stop 'graph_mesh_geometry: every face touches a cell'
    end if

  end subroutine derive_face_cells

  !===================================================================!
  ! Cell centers: the mean of each cell's vertex coordinates.
  !===================================================================!

  pure subroutine cell_centers_of(coordinates, cell_vertices, &
       & num_cell_vertices, cell_centers)

    real(dp), intent(in)  :: coordinates(:,:)
    integer , intent(in)  :: cell_vertices(:,:)
    integer , intent(in)  :: num_cell_vertices(:)
    real(dp), allocatable, intent(out) :: cell_centers(:,:)

    integer :: icell

    allocate(cell_centers(3, size(num_cell_vertices)))

    do icell = 1, size(num_cell_vertices)
       associate(&
            & num_vertices => real(num_cell_vertices(icell), kind=dp), &
            & vids => cell_vertices(1:num_cell_vertices(icell), icell))
         cell_centers(:, icell) = sum(coordinates(:, vids), dim=2)/num_vertices
       end associate
    end do

  end subroutine cell_centers_of

  !===================================================================!
  ! Face centers and areas. In 2d a face is a segment: the center
  ! is its midpoint and the area its length; a zero length stops
  ! the program (two coincident points). In 3d the center is the
  ! vertex mean and the area sums the half cross-products of the
  ! corner fan (one triangle, plus a second for a quadrilateral);
  ! a negative area stops the program.
  !===================================================================!

  pure subroutine face_centers_areas_of(spatial_dim, coordinates, &
       & face_vertices, num_face_vertices, face_centers, face_areas)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: coordinates(:,:)
    integer , intent(in)  :: face_vertices(:,:)
    integer , intent(in)  :: num_face_vertices(:)
    real(dp), allocatable, intent(out) :: face_centers(:,:)
    real(dp), allocatable, intent(out) :: face_areas(:)

    integer  :: iface, num_faces
    real(dp) :: n(3)

    num_faces = size(num_face_vertices)
    allocate(face_areas(num_faces))
    allocate(face_centers(3, num_faces))
    face_areas   = 0.0_dp
    face_centers = 0.0_dp

    if (spatial_dim .eq. 2) then

       do iface = 1, num_faces
          associate(facenodes => face_vertices(1:num_face_vertices(iface), iface))
            face_centers(1:3, iface) = &
                 & sum(coordinates(1:3, facenodes), dim=2)/real(2, kind=dp)
            associate(v1 => coordinates(:, facenodes(1)), &
                 &    v2 => coordinates(:, facenodes(2)))
              face_areas(iface) = distance(v1, v2)
            end associate
          end associate
       end do

       if (abs(minval(face_areas)) .lt. 10.0d0*tiny(1.0d0)) then
          error stop 'graph_mesh_geometry: a face has two distinct end points'
       end if

    else

       do iface = 1, num_faces
          associate(facenodes => face_vertices(1:num_face_vertices(iface), iface))

            associate(num_vertices => real(num_face_vertices(iface), kind=dp))
              face_centers(:, iface) = &
                   & sum(coordinates(:, facenodes), dim=2)/num_vertices
            end associate

            associate(&
                 & t12 => coordinates(:, facenodes(2)) - coordinates(:, facenodes(1)), &
                 & t13 => coordinates(:, facenodes(3)) - coordinates(:, facenodes(1)))

              call cross_product(t12, t13, n)
              face_areas(iface) = norm2(n)/real(2, dp)

              if (num_face_vertices(iface) .gt. 3) then
                 associate(t14 => coordinates(:, facenodes(4)) - coordinates(:, facenodes(1)))
                   call cross_product(t13, t14, n)
                   face_areas(iface) = face_areas(iface) + norm2(n)/real(2, dp)
                 end associate
              end if

            end associate

          end associate
       end do

       if (minval(face_areas) .lt. 0.0_dp) then
          error stop 'graph_mesh_geometry: a face area is nonnegative'
       end if

    end if

  end subroutine face_centers_areas_of

  !===================================================================!
  ! Outward unit normals, one per (cell, local face) pair. In 3d
  ! the normal comes from the corner-fan cross products; in 2d it
  ! is the segment tangent turned a quarter turn in plane. Either
  ! way it is flipped, if needed, to point from the cell center
  ! toward the face center - out of the cell - so the same
  ! interior face carries opposite normals seen from its two
  ! cells, which is what makes the assembled operator
  ! conservative.
  !===================================================================!

  pure subroutine cell_face_normals_of(spatial_dim, coordinates, &
       & face_vertices, num_face_vertices, cell_faces, num_cell_faces, &
       & face_centers, cell_centers, cell_face_normals)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: coordinates(:,:)
    integer , intent(in)  :: face_vertices(:,:)
    integer , intent(in)  :: num_face_vertices(:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), intent(in)  :: face_centers(:,:)
    real(dp), intent(in)  :: cell_centers(:,:)
    real(dp), allocatable, intent(out) :: cell_face_normals(:,:,:)

    integer  :: icell, iface, gface, num_cells
    integer  :: fv1, fv2
    real(dp) :: t(3), n(3), tmp(3)

    num_cells = size(num_cell_faces)
    allocate(cell_face_normals(3, maxval(num_cell_faces), num_cells))
    cell_face_normals = 0.0_dp

    do icell = 1, num_cells
       do iface = 1, num_cell_faces(icell)

          gface = cell_faces(iface, icell)

          if (spatial_dim .eq. 2) then

             fv1 = face_vertices(1, gface)
             fv2 = face_vertices(2, gface)

             t = coordinates(:, fv2) - coordinates(:, fv1)
             t = t/norm2(t)

             n(1) =  t(2)
             n(2) = -t(1)
             n(3) =  0.0d0

          else

             associate(ifv => face_vertices(1:num_face_vertices(gface), gface))
               associate(&
                    & t12 => coordinates(:, ifv(2)) - coordinates(:, ifv(1)), &
                    & t13 => coordinates(:, ifv(3)) - coordinates(:, ifv(1)))

                 call cross_product(t12, t13, n)

                 if (num_face_vertices(gface) .gt. 3) then
                    associate(t14 => coordinates(:, ifv(4)) - coordinates(:, ifv(1)))
                      call cross_product(t13, t14, tmp)
                      n = n + tmp
                    end associate
                 end if

               end associate
             end associate

             n = n/norm2(n)

          end if

          if (dot_product(n, face_centers(:, gface) - cell_centers(:, icell)) &
               & .lt. real(0, dp)) then
             n = -n
          end if

          cell_face_normals(:, iface, icell) = n

       end do
    end do

  end subroutine cell_face_normals_of

  !===================================================================!
  ! Cell volumes by the divergence theorem,
  !
  !      V = (1/d) sum over the cell's faces of  A_f (n_f . x_f),
  !
  ! with d the spatial dimension. A negative volume is an
  ! inside-out cell and stops the program.
  !===================================================================!

  pure subroutine cell_volumes_of(spatial_dim, face_centers, face_areas, &
       & cell_face_normals, cell_faces, num_cell_faces, cell_volumes)

    integer , intent(in)  :: spatial_dim
    real(dp), intent(in)  :: face_centers(:,:)
    real(dp), intent(in)  :: face_areas(:)
    real(dp), intent(in)  :: cell_face_normals(:,:,:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), allocatable, intent(out) :: cell_volumes(:)

    integer :: lcell, lface, gface

    allocate(cell_volumes(size(num_cell_faces)))
    cell_volumes = 0.0_dp

    do lcell = 1, size(num_cell_faces)
       do lface = 1, num_cell_faces(lcell)
          gface = cell_faces(lface, lcell)
          associate(&
               & face_center => face_centers(:, gface)/real(spatial_dim, dp), &
               & face_normal => cell_face_normals(:, lface, lcell), &
               & face_area   => face_areas(gface))
            cell_volumes(lcell) = cell_volumes(lcell) + &
                 & face_area*dot_product(face_normal, face_center)
          end associate
       end do
    end do

    if (minval(cell_volumes) .lt. 0.0_dp) then
       error stop 'graph_mesh_geometry: a cell volume is nonnegative'
    end if

  end subroutine cell_volumes_of

  !===================================================================!
  ! The center-to-center vector of each face: cell center to cell
  ! center across an interior face, cell center to face center on
  ! a boundary face. This is the segment a two-point gradient
  ! differences along.
  !===================================================================!

  pure subroutine centroidal_vectors_of(face_cells, num_face_cells, &
       & cell_centers, face_centers, lvec)

    integer , intent(in)  :: face_cells(:,:)
    integer , intent(in)  :: num_face_cells(:)
    real(dp), intent(in)  :: cell_centers(:,:)
    real(dp), intent(in)  :: face_centers(:,:)
    real(dp), allocatable, intent(out) :: lvec(:,:)

    integer :: iface

    allocate(lvec(3, size(num_face_cells)))
    lvec = real(0, dp)

    do iface = 1, size(num_face_cells)
       if (num_face_cells(iface) .eq. 1) then
          lvec(:, iface) = face_centers(:, iface) &
               & - cell_centers(:, face_cells(1, iface))
       else
          lvec(:, iface) = cell_centers(:, face_cells(2, iface)) &
               & - cell_centers(:, face_cells(1, iface))
       end if
    end do

  end subroutine centroidal_vectors_of

  !===================================================================!
  ! Face deltas: the center-to-center vector projected on the face
  ! normal. On a skewed mesh the segment crosses its face at a
  ! slant, so this normal distance is shorter than the segment; it
  ! is the denominator of every two-point face gradient.
  !===================================================================!

  pure subroutine face_deltas_of(num_faces, lvec, cell_face_normals, &
       & cell_faces, num_cell_faces, face_deltas)

    integer , intent(in)  :: num_faces
    real(dp), intent(in)  :: lvec(:,:)
    real(dp), intent(in)  :: cell_face_normals(:,:,:)
    integer , intent(in)  :: cell_faces(:,:)
    integer , intent(in)  :: num_cell_faces(:)
    real(dp), allocatable, intent(out) :: face_deltas(:)

    integer :: lcell, lface, gface

    allocate(face_deltas(num_faces))
    face_deltas = real(0, dp)

    do lcell = 1, size(num_cell_faces)
       do lface = 1, num_cell_faces(lcell)
          gface = cell_faces(lface, lcell)
          face_deltas(gface) = &
               & abs(dot_product(lvec(:, gface), &
               &                 cell_face_normals(:, lface, lcell)))
       end do
    end do

  end subroutine face_deltas_of

  !===================================================================!
  ! Interpolation weights per face: the two cells' shares by
  ! inverse distance to the face center, stored as (w, 1-w). A
  ! boundary face has one cell, which takes the whole weight.
  !===================================================================!

  pure subroutine face_weights_of(face_cells, num_face_cells, &
       & cell_centers, face_centers, face_cell_weights)

    integer , intent(in)  :: face_cells(:,:)
    integer , intent(in)  :: num_face_cells(:)
    real(dp), intent(in)  :: cell_centers(:,:)
    real(dp), intent(in)  :: face_centers(:,:)
    real(dp), allocatable, intent(out) :: face_cell_weights(:,:)

    integer  :: iface
    real(dp) :: d1, d2, dinv1, dinv2, weight

    allocate(face_cell_weights(2, size(num_face_cells)))

    do iface = 1, size(num_face_cells)

       d1    = distance(cell_centers(:, face_cells(1, iface)), &
            &           face_centers(:, iface))
       dinv1 = 1.0_dp/d1

       if (num_face_cells(iface) .ne. 1) then
          d2    = distance(cell_centers(:, face_cells(2, iface)), &
               &           face_centers(:, iface))
          dinv2 = 1.0_dp/d2
       else
          dinv2 = 0.0_dp
       end if

       weight = dinv1/(dinv1 + dinv2)

       face_cell_weights(1:2, iface) = [weight, 1.0_dp - weight]

    end do

  end subroutine face_weights_of


  !===================================================================!
  ! Compute the cross product for area computations.
  !===================================================================!

  pure subroutine cross_product(a, b, pdt)

    real(dp), intent(in)  :: a(3), b(3)
    real(dp), intent(out) :: pdt(3)

    ! The skew-symmetric form would generalize this product to n dimensions.
    pdt(1) = a(2) * b(3) - a(3) * b(2)
    pdt(2) = a(3) * b(1) - a(1) * b(3)
    pdt(3) = a(1) * b(2) - a(2) * b(1)

  end subroutine cross_product

  !===================================================================!
  ! Compute the geometric distance between two points.
  !===================================================================!

  pure real(dp) function distance(x, y)

    real(dp), intent(in)  :: X(:), y(:) ! The shape is [[x,y,z], [1:2]].

    distance = sqrt(sum((x-y)**2))

  end function distance

  !===================================================================!
  ! Return the index of a target value if it is present in the array;
  ! return -1 otherwise.
  !===================================================================!

  pure type(integer) function find(array, target_value)

    integer, intent(in) :: array(:)
    integer, intent(in) :: target_value
    integer :: i, nentries

    nentries = size(array, dim=1)

    do i = 1, nentries
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

  pure elemental integer function elem_type_face_count(cell_type) result (num_faces)

    integer, intent(in) :: cell_type

    select case (cell_type)
    case (1)
       ! A 2-node line.
       num_faces = 2 ! The faces are its two vertices.
    case (2)
       ! A 3-node triangle.
       num_faces = 3  ! The faces are its three edges.
    case (3)
       ! A 4-node quadrangle.
       num_faces = 4 ! The faces are its four edges.
    case (4)
       ! A 4-node tetrahedron.
       num_faces = 4
    case (5)
       ! An 8-node hexahedron.
       num_faces = 6
    case (6)
       ! A 6-node prism.
       num_faces = 5
    case(7)
       ! A 5-node prism (a pyramid).
       num_faces = 5
    case default
       num_faces = 0
    end select

  end function elem_type_face_count

  !===================================================================!
  ! Return the spatial dimension of a gmsh element type: a point is
  ! 0, a line is 1, a triangle or quadrangle is 2, and a tet, hex,
  ! prism, or pyramid is 3. The mesh classifies elements into cells,
  ! faces, and edges by dimension rather than by type.
  !===================================================================!

  pure elemental integer function elem_type_dimension(elem_type) result (dim)

    integer, intent(in) :: elem_type

    select case (elem_type)
    case (15)
       ! A 1-node point.
       dim = 0
    case (1)
       ! A 2-node line.
       dim = 1
    case (2:3)
       ! A triangle or a quadrangle.
       dim = 2
    case (4:7)
       ! A tet, hex, prism, or pyramid.
       dim = 3
    case default
       dim = -1
    end select

  end function elem_type_dimension

  !===================================================================!
  ! Return the number of vertices a gmsh element type owns. The name
  ! could generalize to the number of lower dimensional entities.
  !===================================================================!

  pure elemental integer function elem_type_vertex_count(elem_type) result (num_vertices)

    integer, intent(in) :: elem_type

    select case (elem_type)
    case (1)
       ! A 2-node line.
       num_vertices = 2
    case (2)
       ! A 3-node triangle.
       num_vertices = 3
    case (3)
       ! A 4-node quadrangle.
       num_vertices = 4
    case (4)
       ! A 4-node tetrahedron.
       num_vertices = 4
    case (5)
       ! An 8-node hexahedron.
       num_vertices = 8
    case (6)
       ! A 6-node prism.
       num_vertices = 6
    case(7)
       ! A 5-node prism (a pyramid).
       num_vertices = 5
    case default
       num_vertices = 0
    end select

  end function elem_type_vertex_count

  !===================================================================!
  ! Put a face's vertices into the cell's own winding order. Each
  ! cell type carries a wiring table that says which corners bound
  ! which face, wound so that the normal points outward. Match the
  ! unordered face against the table and hand it back in the table's
  ! order.
  !===================================================================!

  impure subroutine order_face_vertices(cell_type, cell_vertices, face_vertices_unordered)

    integer, intent(in)    :: cell_type
    integer, intent(in)    :: cell_vertices(:)
    integer, intent(inout) :: face_vertices_unordered(:)

    integer, allocatable :: face_vertices(:,:)
    integer :: num_face_vertices, num_cell_faces
    integer :: iface, ivertex
    integer :: match_count

    select case (cell_type)
    case (2)
       ! A 3-node triangle (a 2d cell); its faces are its three edges.
       num_face_vertices = 2
       num_cell_faces = 3

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3)]
       face_vertices(:,3) = [cell_vertices(3), cell_vertices(1)]

    case (3)
       ! A 4-node quadrangle (a 2d cell); its faces are its four edges.
       num_face_vertices = 2
       num_cell_faces = 4

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3)]
       face_vertices(:,3) = [cell_vertices(3), cell_vertices(4)]
       face_vertices(:,4) = [cell_vertices(4), cell_vertices(1)]

    case (4)
       ! A 4-node tetrahedron.
       num_face_vertices = 3
       num_cell_faces = 4

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       ! The normals must point outward.
       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2), cell_vertices(3)]
       face_vertices(:,2) = [cell_vertices(3), cell_vertices(2), cell_vertices(4)]
       face_vertices(:,3) = [cell_vertices(1), cell_vertices(3), cell_vertices(4)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(4), cell_vertices(2)]

    case (5)
       ! An 8-node hexahedron.
       num_face_vertices = 4
       num_cell_faces = 6

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(5), cell_vertices(8), cell_vertices(4)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3), cell_vertices(7), cell_vertices(6)]
       face_vertices(:,3) = [cell_vertices(5), cell_vertices(6), cell_vertices(7), cell_vertices(8)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(4), cell_vertices(3), cell_vertices(2)]
       face_vertices(:,5) = [cell_vertices(7), cell_vertices(3), cell_vertices(4), cell_vertices(8)]
       face_vertices(:,6) = [cell_vertices(1), cell_vertices(2), cell_vertices(6), cell_vertices(5)]

    case (6)
       ! A 6-node prism with five faces.
       num_cell_faces = 5
       num_face_vertices = 4 ! Take the maximum.

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(4), cell_vertices(5), cell_vertices(6)]
       face_vertices(:,2) = [cell_vertices(1), cell_vertices(3), cell_vertices(6), cell_vertices(4)]
       face_vertices(:,3) = [cell_vertices(2), cell_vertices(5), cell_vertices(6), cell_vertices(3)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(2), cell_vertices(3)]
       face_vertices(:,5) = [cell_vertices(4), cell_vertices(6), cell_vertices(5)]

    case(7)
       ! A 5-node prism (a pyramid) with five faces.
       num_cell_faces = 5
       num_face_vertices = 4 ! Take the maximum.

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2), cell_vertices(3), cell_vertices(4)]
       face_vertices(:,2) = [cell_vertices(1), cell_vertices(2), cell_vertices(5)]
       face_vertices(:,3) = [cell_vertices(2), cell_vertices(3), cell_vertices(5)]
       face_vertices(:,4) = [cell_vertices(3), cell_vertices(4), cell_vertices(5)]
       face_vertices(:,5) = [cell_vertices(1), cell_vertices(5), cell_vertices(4)]

    case default
       print *, "unknown type"
       error stop
    end select

    !-----------------------------------------------------------------!
    ! Compare the input with every candidate face and find the face
    ! whose matching count equals the input vertex count.
    !-----------------------------------------------------------------!

    match_face: do iface = 1, num_cell_faces

       match_count = 0
       match_vertex: do ivertex = 1, num_face_vertices
          if (any (face_vertices_unordered(:) .eq. face_vertices(ivertex,iface) ) .eqv. .true.) then
             ! exit match_vertex ! skip and go to next face
             ! else
             match_count = match_count + 1
          end if
       end do match_vertex

       ! This is the face with the correct ordering of vertices.
       if (match_count .eq. num_face_vertices) then
          face_vertices_unordered = face_vertices(:, iface)
          exit match_face
       end if

    end do match_face

    if (allocated(face_vertices)) deallocate(face_vertices)

  end subroutine order_face_vertices


end module graph_mesh_geometry
