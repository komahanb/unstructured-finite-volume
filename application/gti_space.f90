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
! spacing along each coordinate is given as the time grid's is - by
! counts, and uniform or drawn from a seed - and every cell is the
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
  use operation_stencil         , only : stencil
  use operation_diffusion       , only : diffusion_stencil
  use operation_conduction      , only : conduction
  use operation_robin_condition , only : robin_condition, neumann
  use field_forms               , only : polynomial_form

  implicit none

  private
  public :: room, spatial_mesh, spatial_operator, written_vtk, coarse_cells
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

    select case (trim(name))
    case ('cartesian')
       geometry = cartesian
    case ('circular')
       geometry = circular
    case ('elliptical')
       geometry = elliptical
    case default
       write(*,'(a)') ' spatial_geometry names ' // trim(name) // &
            & ', which this program has nothing for.'
       error stop 'gti_space: a geometry is cartesian, circular or elliptical'
    end select

  end function geometry_of

  !-------------------------------------------------------------------!
  ! The spacing along one coordinate: n widths, uniform or drawn, that
  ! sum to one. The draw is the time grid's, so a seed means the same
  ! thing in space as in time.
  !-------------------------------------------------------------------!

  pure subroutine spacing(n, drawn, seed, salt, u)

    integer, intent(in) :: n
    logical, intent(in) :: drawn
    integer, intent(in) :: seed, salt
    real(dp), allocatable, intent(out) :: u(:)

    integer, parameter :: modulus = 2147483647
    real(dp) :: w(n)
    integer  :: k, x

    ! allocated here with its own lower bound, which an assignment
    ! from a function result would not keep
    allocate(u(0:n))

    w = 1.0_dp
    if (drawn) then
       do k = 1, n
          x = modulo(seed * 40503 + (k + salt) * 65537, modulus)
          x = modulo(int(mod(1103515245_8 * int(x, 8) + 12345_8, int(modulus, 8))), modulus)
          w(k) = 0.5_dp + real(x, dp) / real(modulus, dp)
       end do
    end if

    u(0) = 0.0_dp
    do k = 1, n
       u(k) = u(k - 1) + w(k)
    end do
    u = u / u(n)

  end subroutine spacing

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

    real(dp), allocatable :: xi(:), eta(:)
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

    call spacing(n1, drawn, seed, 0, xi)
    call spacing(n2, drawn, seed, n1, eta)

    allocate(this % corner(2, (n1 + 1) * (n2 + 1)))
    do j = 0, n1
       do i = 0, n2
          this % corner(:, corner_index(i, j, n2)) = mapped(geometry, a, b, xi(j), eta(i))
       end do
    end do

    if (polar == 1) then
       cells = 1 + (n1 - 1) * n2
    else
       cells = n1 * n2
    end if
    this % num_cells = cells

    allocate(this % first_corner(cells + 1))
    allocate(this % cell_corner(4 * n1 * n2 + n2))
    allocate(this % centre(2, cells), this % volume(cells), this % cell_ij(2, cells))
    this % n1 = n1
    this % n2 = n2

    c = 0
    this % first_corner(1) = 1

    if (polar == 1) then
       c = 1
       do i = 0, n2 - 1
          this % cell_corner(i + 1) = corner_index(i, 1, n2)
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

    call polygon_geometry(this)

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
               & cell_index(i, j + 1, n2, polar), corner_index(i - 1, j, n2), &
               & corner_index(i, j, n2))
       end do
    end do

    if (polar == 1) then
       do i = 1, n2
          call face_between(faces, f, 1, cell_index(i, 2, n2, polar), &
               & corner_index(i - 1, 1, n2), corner_index(i, 1, n2))
       end do
    end if

    do j = 1 + polar, n1
       do i = 1, n2 - 1 + polar
          ring = i + 1
          if (ring > n2) ring = 1
          call face_between(faces, f, cell_index(i, j, n2, polar), &
               & cell_index(ring, j, n2, polar), corner_index(i, j - 1, n2), &
               & corner_index(i, j, n2))
       end do
    end do

    do i = 1, n2
       call face_between(faces, f, cell_index(i, n1, n2, polar), 0, &
            & corner_index(i - 1, n1, n2), corner_index(i, n1, n2))
    end do

    if (polar == 0) then
       do i = 1, n2
          call face_between(faces, f, cell_index(i, 1, n2, polar), 0, &
               & corner_index(i - 1, 0, n2), corner_index(i, 0, n2))
       end do
       do j = 1, n1
          call face_between(faces, f, cell_index(1, j, n2, polar), 0, &
               & corner_index(0, j - 1, n2), corner_index(0, j, n2))
          call face_between(faces, f, cell_index(n2, j, n2, polar), 0, &
               & corner_index(n2, j - 1, n2), corner_index(n2, j, n2))
       end do
    end if

    this % num_faces = f
    call framed(this, faces(1:f))

  end function spatial_mesh

  pure integer function corner_index(i, j, n2) result(c)

    integer, intent(in) :: i, j, n2

    c = j * (n2 + 1) + i + 1

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

    integer :: at

    at = this % first_corner(c)
    this % cell_corner(at)     = corner_index(i - 1, j - 1, n2)
    this % cell_corner(at + 1) = corner_index(i,     j - 1, n2)
    this % cell_corner(at + 2) = corner_index(i,     j,     n2)
    this % cell_corner(at + 3) = corner_index(i - 1, j,     n2)
    this % first_corner(c + 1) = at + 4

  end subroutine quad

  !-------------------------------------------------------------------!
  ! Area and centroid of every cell from its corners, by the shoelace
  ! sums. A cell of no area stops the program.
  !-------------------------------------------------------------------!

  subroutine polygon_geometry(this)

    type(room), intent(inout) :: this

    real(dp) :: area, cx, cy, cross
    real(dp) :: p(2), q(2)
    integer  :: c, k, m, n

    do c = 1, this % num_cells
       n    = this % first_corner(c + 1) - this % first_corner(c)
       area = 0.0_dp
       cx   = 0.0_dp
       cy   = 0.0_dp
       do k = 0, n - 1
          m     = mod(k + 1, n)
          p     = this % corner(:, this % cell_corner(this % first_corner(c) + k))
          q     = this % corner(:, this % cell_corner(this % first_corner(c) + m))
          cross = p(1) * q(2) - q(1) * p(2)
          area  = area + cross
          cx    = cx + (p(1) + q(1)) * cross
          cy    = cy + (p(2) + q(2)) * cross
       end do
       area = 0.5_dp * area
       if (abs(area) <= tiny(1.0_dp)) then
          error stop 'gti_space: every cell has area'
       end if
       this % volume(c)    = abs(area)
       this % centre(:, c) = [cx, cy] / (6.0_dp * area)
    end do

  end subroutine polygon_geometry

  subroutine face_between(faces, f, tail, head, corner_a, corner_b)

    type(face_record), intent(inout) :: faces(:)
    integer          , intent(inout) :: f
    integer          , intent(in)    :: tail, head, corner_a, corner_b

    f = f + 1
    faces(f) = face_record(tail, head, corner_a, corner_b)

  end subroutine face_between

  !-------------------------------------------------------------------!
  ! The framework's mesh from the polygons: one face per record with
  ! its length, its unit normal out of the tail, its centre, the
  ! distance between the centroids projected on the normal, the
  ! tail's inverse-distance share, and a tag on the headless ones.
  ! The conventions are view_mesh_builder's, read there.
  !-------------------------------------------------------------------!

  subroutine framed(this, faces)

    type(room)       , intent(inout) :: this
    type(face_record), intent(in)    :: faces(:)

    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: areas(:), deltas(:), normals(:), centres(:), weights(:), &
         & cell_centres(:)
    character(len=4), allocatable :: tags(:)
    real(dp) :: a(2), b(2), n(2), xf(2), ct(2), ch(2), d1, d2
    integer  :: f, nf, c

    nf = size(faces)
    allocate(tails(nf), heads(nf), areas(nf), deltas(nf), normals(3 * nf), &
         & centres(3 * nf), weights(nf), tags(nf), cell_centres(3 * this % num_cells))

    do c = 1, this % num_cells
       cell_centres(3 * c - 2:3 * c) = [this % centre(1, c), this % centre(2, c), 0.0_dp]
    end do

    do f = 1, nf
       tails(f) = faces(f) % tail
       heads(f) = faces(f) % head

       a  = this % corner(:, faces(f) % corner_a)
       b  = this % corner(:, faces(f) % corner_b)
       xf = 0.5_dp * (a + b)
       ct = this % centre(:, faces(f) % tail)

       areas(f) = norm2(b - a)
       n = [b(2) - a(2), a(1) - b(1)] / areas(f)
       if (dot_product(n, xf - ct) < 0.0_dp) n = -n

       normals(3 * f - 2:3 * f) = [n(1), n(2), 0.0_dp]
       centres(3 * f - 2:3 * f) = [xf(1), xf(2), 0.0_dp]

       d1 = norm2(ct - xf)
       if (faces(f) % head > 0) then
          ch         = this % centre(:, faces(f) % head)
          deltas(f)  = abs(dot_product(ch - ct, n))
          d2         = norm2(ch - xf)
          weights(f) = (1.0_dp / d1) / (1.0_dp / d1 + 1.0_dp / d2)
          tags(f)    = ''
       else
          deltas(f)  = abs(dot_product(xf - ct, n))
          weights(f) = 1.0_dp
          tags(f)    = 'wall'
       end if
    end do

    this % m = mesh(this % num_cells, tails=tails, heads=heads, &
         & volumes=this % volume, cell_centres=cell_centres, areas=areas, &
         & deltas=deltas, normals=normals, face_centres=centres, weights=weights, &
         & etags=tags)

  end subroutine framed

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
    op = diffusion_stencil(this % m, conduction(kappa), wall, polynomial_form(degree))

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
  ! One legacy vtk file: the cells as polygons, and one scalar per
  ! cell for each name given. A numbered series of these is read by
  ! paraview as steps in time.
  !-------------------------------------------------------------------!

  subroutine written_vtk(this, path, names, values)

    type(room)      , intent(in) :: this
    character(len=*), intent(in) :: path, names(:)
    real(dp)        , intent(in) :: values(:,:)

    integer :: u, c, k, n, total

    if (size(values, 1) /= this % num_cells .or. size(values, 2) /= size(names)) then
       error stop 'gti_space: one value per cell per name'
    end if

    open(newunit=u, file=path, action='write', status='replace')
    write(u,'(a)') '# vtk DataFile Version 3.0'
    write(u,'(a)') 'graph time integrator, one instant'
    write(u,'(a)') 'ASCII'
    write(u,'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(u,'(a,i0,a)') 'POINTS ', size(this % corner, 2), ' double'
    do c = 1, size(this % corner, 2)
       write(u,'(3es24.15)') this % corner(1, c), this % corner(2, c), 0.0_dp
    end do

    total = this % first_corner(this % num_cells + 1) - 1 + this % num_cells
    write(u,'(a,i0,a,i0)') 'CELLS ', this % num_cells, ' ', total
    do c = 1, this % num_cells
       n = this % first_corner(c + 1) - this % first_corner(c)
       write(u,'(i0,*(1x,i0))') n, (this % cell_corner(k) - 1, &
            & k = this % first_corner(c), this % first_corner(c + 1) - 1)
    end do

    write(u,'(a,i0)') 'CELL_TYPES ', this % num_cells
    do c = 1, this % num_cells
       n = this % first_corner(c + 1) - this % first_corner(c)
       write(u,'(i0)') merge(9, 7, n == 4)
    end do

    write(u,'(a,i0)') 'CELL_DATA ', this % num_cells
    do k = 1, size(names)
       write(u,'(a)') 'SCALARS ' // trim(names(k)) // ' double 1'
       write(u,'(a)') 'LOOKUP_TABLE default'
       do c = 1, this % num_cells
          write(u,'(es24.15)') values(c, k)
       end do
    end do

    close(u)

  end subroutine written_vtk

end module gti_space
