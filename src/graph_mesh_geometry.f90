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
!      V2C            transpose(C2V), by graph_binary_relation's
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
  use module_mesh_utils, only : distance, cross_product, &
       & elem_type_face_count, order_face_vertices

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

end module graph_mesh_geometry
