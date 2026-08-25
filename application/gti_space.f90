! THE SPATIAL LEVEL: the framework's mesh over a two-dimensional
! domain, and the framework's operator on it.
!
! The mesh is structured in two parametric coordinates on the unit
! square and mapped to the domain by its geometry:
!
!      cartesian     x = a xi,          y = b eta
!      circular      x = a xi cos 2 pi eta, y = a xi sin 2 pi eta
!      elliptical    x = a xi cos 2 pi eta, y = b xi sin 2 pi eta
!
! so one indexing serves every shape and the shape is a mapping. The
! spacing along each coordinate IS the time grid's, uniform or drawn
! from a seed through the same partition, so a seed means the same
! thing in space as in time - and every cell is the
! polygon of its mapped corners: its area, centroid and face geometry
! are read off those corners, whatever the shape. A polar mapping
! collapses the inner corners onto the origin, so the innermost ring
! is one cell, a polygon of the first ring's corners.
!
! What is built is the mesh view_mesh already defines - cells as
! vertices, faces as edges, a boundary face an edge without a head
! and tagged - and the operator is operation_diffusion's: the fitted
! balance, a polynomial form of the order asked for fitted on each
! face's neighbourhood and aimed along the face normal, so a skewed
! face is measured in the right direction and the order is the
! form's. The outer boundary holds no flux. Nothing spatial is
! discretised here; this module maps a shape and hands the mesh on.
!
! The operator is the flux balance, integrated over each cell; a
! caller wanting the laplacian divides each row by its cell's area.
module gti_space

  use util_precision            , only : dp
  use view_mesh                 , only : mesh
  use view_mesh_geometry        , only : mesh_from_incidence
  use field_stored              , only : stored_field
  use operation_stencil         , only : stencil
  use operation_diffusion       , only : diffusion_stencil
  use operation_conduction      , only : conduction
  use operation_robin_condition , only : robin_condition, neumann
  use field_forms               , only : polynomial_form
  use view_paraview_writer      , only : paraview_writer, polygon_cell
  use relation_binary           , only : ragged
  use util_string               , only : string
  use operation_grid            , only : uniform_grid, random_grid
  use gti_configuration         , only : chosen_from
  use gti_march                 , only : partitioned

  implicit none

  private
  public :: room, spatial_mesh, spatial_operator, written_paraview, coarse_cells
  public :: cartesian, circular, elliptical, geometry_of

  integer, parameter :: cartesian  = 1
  integer, parameter :: circular   = 2
  integer, parameter :: elliptical = 3

  type :: room

     integer :: geometry  = cartesian
     integer :: num_cells = 0
     integer :: num_faces = 0

     ! the framework's mesh
     type(mesh) :: m

     ! corners, and each cell's corners in order, ragged - kept for
     ! writing the cells out as polygons
     real(dp), allocatable :: corner(:,:)
     integer , allocatable :: first_corner(:)
     integer , allocatable :: cell_corner(:)

     real(dp), allocatable :: centre(:,:)
     real(dp), allocatable :: volume(:)

     ! each cell's place on the parametric grid, for coarsening
     integer , allocatable :: cell_ij(:,:)
     integer :: n1 = 0, n2 = 0

  end type room

  ! a face while the mesh is being built: its cells, its corners
  type :: face_record
     integer :: tail = 0, head = 0, corner_a = 0, corner_b = 0
  end type face_record

contains

  !-------------------------------------------------------------------!
  ! The geometry a name denotes; an unknown name stops the program.
  !-------------------------------------------------------------------!

  integer function geometry_of(name) result(geometry)

    character(len=*), intent(in) :: name

    ! the constants are the places in this list
    geometry = chosen_from(name, ['cartesian ', 'circular  ', 'elliptical'], 'spatial_geometry')

  end function geometry_of

  pure function mapped(geometry, a, b, xi, eta) result(x)

    integer , intent(in) :: geometry
    real(dp), intent(in) :: a, b, xi, eta
    real(dp) :: x(2)

    real(dp) :: theta

    select case (geometry)
    case (cartesian)
       x = [a * xi, b * eta]
    case (circular)
       theta = 2.0_dp * acos(-1.0_dp) * eta
       x = [a * xi * cos(theta), a * xi * sin(theta)]
    case default
       theta = 2.0_dp * acos(-1.0_dp) * eta
       x = [a * xi * cos(theta), b * xi * sin(theta)]
    end select

  end function mapped

  !-------------------------------------------------------------------!
  ! The mesh: counts along the two coordinates, the extents a and b,
  ! the spacing, and the geometry. A count below two along either
  ! coordinate, or an extent that is not positive, stops the program.
  !-------------------------------------------------------------------!

  function spatial_mesh(geometry, a, b, n1, n2, drawn, seed) result(this)

    integer , intent(in) :: geometry, n1, n2, seed
    real(dp), intent(in) :: a, b
    logical , intent(in) :: drawn
    type(room) :: this

    real(dp), allocatable :: xi(:), eta(:), dxi(:), deta(:)
    type(face_record), allocatable :: faces(:)
    integer :: i, j, c, f, polar, cells, ring

    if (n1 < 2 .or. n2 < 2) then
       error stop 'gti_space: at least two cells along each coordinate'
    end if
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
       error stop 'gti_space: an extent is positive'
    end if

    this % geometry = geometry
    polar = merge(1, 0, geometry /= cartesian)

    ! the spacing along each coordinate is the time grid's own draw,
    ! scaled to the unit interval; the second coordinate continues
    ! the draw past the first, as a seed offset by the first's count
    if (drawn) then
       call partitioned(random_grid(1.0_dp, seed),      n1 + 1, dxi,  xi)
       call partitioned(random_grid(1.0_dp, seed + n1), n2 + 1, deta, eta)
    else
       call partitioned(uniform_grid(1.0_dp), n1 + 1, dxi,  xi)
       call partitioned(uniform_grid(1.0_dp), n2 + 1, deta, eta)
    end if

    allocate(this % corner(2, (n1 + 1 - polar) * (n2 + 1)))
    do j = polar, n1
       do i = 0, n2
          this % corner(:, corner_index(i, j, n2, polar)) = mapped(geometry, a, b, xi(j + 1), eta(i + 1))
       end do
    end do

    if (polar == 1) then
       cells = 1 + (n1 - 1) * n2
    else
       cells = n1 * n2
    end if
    this % num_cells = cells

    allocate(this % first_corner(cells + 1))
    allocate(this % cell_corner(merge(n2 + 4 * (n1 - 1) * n2, 4 * n1 * n2, polar == 1)))
    allocate(this % centre(2, cells), this % volume(cells), this % cell_ij(2, cells))
    this % n1 = n1
    this % n2 = n2

    c = 0
    this % first_corner(1) = 1

    if (polar == 1) then
       c = 1
       do i = 0, n2 - 1
          this % cell_corner(i + 1) = corner_index(i, 1, n2, polar)
       end do
       this % first_corner(2) = n2 + 1
       this % cell_ij(:, 1) = [0, 1]
    end if

    do j = 1 + polar, n1
       do i = 1, n2
          c = c + 1
          call quad(this, c, i, j, n2)
          this % cell_ij(:, c) = [i, j]
       end do
    end do

    !----------------------------------------------------------------!
    ! Faces. Interior ones along the first coordinate between rings
    ! or columns and along the second between neighbours, periodic
    ! in the second for a polar mapping; boundary ones on the outer
    ! edge of the domain, headless.
    !----------------------------------------------------------------!

    allocate(faces(2 * n1 * n2 + 2 * (n1 + n2) + n2))
    f = 0

    do j = 1 + polar, n1 - 1
       do i = 1, n2
          call face_between(faces, f, cell_index(i, j, n2, polar), &
               & cell_index(i, j + 1, n2, polar), corner_index(i - 1, j, n2, polar), &
               & corner_index(i, j, n2, polar))
       end do
    end do

    if (polar == 1) then
       do i = 1, n2
          call face_between(faces, f, 1, cell_index(i, 2, n2, polar), &
               & corner_index(i - 1, 1, n2, polar), corner_index(i, 1, n2, polar))
       end do
    end if

    do j = 1 + polar, n1
       do i = 1, n2 - 1 + polar
          ring = i + 1
          if (ring > n2) ring = 1
          call face_between(faces, f, cell_index(i, j, n2, polar), &
               & cell_index(ring, j, n2, polar), corner_index(i, j - 1, n2, polar), &
               & corner_index(i, j, n2, polar))
       end do
    end do

    do i = 1, n2
       call face_between(faces, f, cell_index(i, n1, n2, polar), 0, &
            & corner_index(i - 1, n1, n2, polar), corner_index(i, n1, n2, polar))
    end do

    if (polar == 0) then
       do i = 1, n2
          call face_between(faces, f, cell_index(i, 1, n2, polar), 0, &
               & corner_index(i - 1, 0, n2, polar), corner_index(i, 0, n2, polar))
       end do
       do j = 1, n1
          call face_between(faces, f, cell_index(1, j, n2, polar), 0, &
               & corner_index(0, j - 1, n2, polar), corner_index(0, j, n2, polar))
          call face_between(faces, f, cell_index(n2, j, n2, polar), 0, &
               & corner_index(n2, j - 1, n2, polar), corner_index(n2, j, n2, polar))
       end do
    end if

    this % num_faces = f
    call measured(this, faces(1:f))

  end function spatial_mesh

  pure integer function corner_index(i, j, n2, polar) result(c)

    integer, intent(in) :: i, j, n2, polar

    ! a polar mesh's first row of corners is its innermost ring's;
    ! the row below it, every point the origin, is never made
    c = (j - polar) * (n2 + 1) + i + 1

  end function corner_index

  pure integer function cell_index(i, j, n2, polar) result(c)

    integer, intent(in) :: i, j, n2, polar

    if (polar == 1) then
       c = 1 + (j - 2) * n2 + i
    else
       c = (j - 1) * n2 + i
    end if

  end function cell_index

  subroutine quad(this, c, i, j, n2)

    type(room), intent(inout) :: this
    integer   , intent(in)    :: c, i, j, n2

    integer :: at, polar

    polar = merge(1, 0, this % geometry /= cartesian)
    at = this % first_corner(c)
    this % cell_corner(at)     = corner_index(i - 1, j - 1, n2, polar)
    this % cell_corner(at + 1) = corner_index(i,     j - 1, n2, polar)
    this % cell_corner(at + 2) = corner_index(i,     j,     n2, polar)
    this % cell_corner(at + 3) = corner_index(i - 1, j,     n2, polar)
    this % first_corner(c + 1) = at + 4

  end subroutine quad

  subroutine face_between(faces, f, tail, head, corner_a, corner_b)

    type(face_record), intent(inout) :: faces(:)
    integer          , intent(inout) :: f
    integer          , intent(in)    :: tail, head, corner_a, corner_b

    f = f + 1
    faces(f) = face_record(tail, head, corner_a, corner_b)

  end subroutine face_between

  !-------------------------------------------------------------------!
  ! The mesh from the corners, the cells and the faces enumerated
  ! above, through the framework's one ending of every mesh pipeline:
  ! areas, centroids, normals, deltas and weights are its, computed
  ! as for a mesh read from a file. The outer boundary is the wall.
  ! The cell centres and areas the level reads are then the mesh's.
  !-------------------------------------------------------------------!

  subroutine measured(this, faces)

    type(room)       , intent(inout) :: this
    type(face_record), intent(in)    :: faces(:)

    integer , allocatable :: cell_vertices(:,:), num_cell_vertices(:)
    integer , allocatable :: face_vertices(:,:), num_face_vertices(:), face_cells(:,:), num_face_cells(:)
    character(len=4), allocatable :: tags(:)
    type(ragged) :: corners
    type(stored_field) :: measure
    real(dp), allocatable :: values(:)
    integer :: f, nf

    nf = size(faces)
    corners = ragged(this % first_corner, this % cell_corner)
    call corners % padded(cell_vertices, num_cell_vertices)

    allocate(face_vertices(2, nf), num_face_vertices(nf), face_cells(2, nf), num_face_cells(nf), tags(nf))
    do f = 1, nf
       face_vertices(:, f)  = [faces(f) % corner_a, faces(f) % corner_b]
       num_face_vertices(f) = 2
       face_cells(:, f)     = [faces(f) % tail, faces(f) % head]
       num_face_cells(f)    = merge(2, 1, faces(f) % head > 0)
       tags(f)              = merge('    ', 'wall', faces(f) % head > 0)
    end do

    this % m = mesh_from_incidence(2, this % corner, cell_vertices, num_cell_vertices, &
         & face_vertices, num_face_vertices, face_cells, num_face_cells, tags)

    measure = this % m % cell_centre()
    call measure % real_vector(values)
    this % centre = reshape(values, [2, this % num_cells])
    measure = this % m % cell_volume()
    call measure % real_vector(this % volume)

  end subroutine measured

  !-------------------------------------------------------------------!
  ! The spatial operator: the diffusion statement on the mesh, the
  ! conductivity kappa through every face, no flux at the wall, and
  ! the polynomial form of the degree asked for, each fit taking as
  ! many rings as its members need. What comes back is the flux
  ! balance per cell, integrated.
  !-------------------------------------------------------------------!

  function spatial_operator(this, kappa, degree) result(op)

    type(room), intent(in) :: this
    real(dp)  , intent(in) :: kappa
    integer   , intent(in) :: degree
    type(stencil) :: op

    type(robin_condition) :: wall(1)

    if (degree < 1) then
       error stop 'gti_space: a form of degree below one fits no gradient'
    end if

    wall(1) = neumann('wall', 0.0_dp)
    op = diffusion_stencil(this % m, conduction(kappa), wall, polynomial_form(degree, this % m % dimension))

  end function spatial_operator

  !-------------------------------------------------------------------!
  ! The cells coarsened by pairs along each parametric coordinate:
  ! one aggregate per block of two by two, the innermost polar cell
  ! its own. The aggregate of every cell, numbered from one.
  !-------------------------------------------------------------------!

  function coarse_cells(this) result(aggregate)

    type(room), intent(in) :: this
    integer, allocatable :: aggregate(:)

    integer :: c, i, j, n2c, polar

    polar = merge(1, 0, this % geometry /= cartesian)
    n2c   = (this % n2 + 1) / 2
    allocate(aggregate(this % num_cells))

    do c = 1, this % num_cells
       i = this % cell_ij(1, c)
       j = this % cell_ij(2, c)
       if (polar == 1 .and. c == 1) then
          aggregate(c) = 1
       else
          aggregate(c) = polar + ((j - 1 - polar) / 2) * n2c + (i - 1) / 2 + 1
       end if
    end do

  end function coarse_cells

  !-------------------------------------------------------------------!
  ! One instant written for paraview by the framework's writer: the
  ! cells as polygons in the order their corners lie, one scalar per
  ! cell for each name. A numbered series of these is read as steps
  ! in time.
  !-------------------------------------------------------------------!

  subroutine written_paraview(this, path, names, values)

    type(room)      , intent(in) :: this
    character(len=*), intent(in) :: path, names(:)
    real(dp)        , intent(in) :: values(:,:)

    type(paraview_writer) :: writer

    if (size(values, 1) /= this % num_cells .or. size(values, 2) /= size(names)) then
       error stop 'gti_space: one value per cell per name'
    end if

    writer = paraview_writer(this % m, this % corner, &
         & ragged(this % first_corner, this % cell_corner), &
         & spread(polygon_cell, 1, this % num_cells))
    call writer % write(path, values, string(names))

  end subroutine written_paraview

end module gti_space
