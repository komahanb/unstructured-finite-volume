!=====================================================================!
! The paraview writer: a mesh and its cell fields as one .vtu file,
! the UnstructuredGrid format paraview reads.
!
! WHAT THE CALLER SUPPLIES. The mesh is cell-centred - cells as
! vertices, faces as edges, measurements as fields - and holds no
! corners, so the corners arrive beside it:
!
!      m                    the mesh: its cell count and cell_volume()
!      coordinates(:,:)     (ndim, num_points), any ndim
!      cell_vertices(:,:)   (corner, cell), padded; 1-based points
!      num_cell_vertices(:) the corners each cell owns
!      cell_types(:)        gmsh's numbers (1 line, 2 triangle, 3 quad,
!                           4 tet, 5 hex, 6 prism, 7 pyramid), or
!                           polygon_cell, or hypercube_cell
!
! and the fields to write ride the call to write, one column per
! label; the cell volume is always written first, under 'volume'.
!
!      use view_paraview_writer, only : paraview_writer, hypercube_cell
!      w = paraview_writer(m, coordinates, cell_vertices, &
!                          num_cell_vertices, cell_types)
!      call w % write('state_0001.vtu', phic, labels)
!
! phic(cell, label) is real(dp) and labels(:) is type(string) from
! util_string; both are optional and omitted together. A numbered
! series of files is read by paraview as steps in time.
!
! FROM A GMSH FILE. mesh_from_gmsh discards the corners it reads, so
! a caller wanting the picture reads them once more through the
! loader and passes them on:
!
!      allocate(gl, source=gmsh_loader(file))
!      call gl % mesh_data(num_vertices, ..., vertices, ..., &
!                          cell_vertices, num_cell_vertices, &
!                          cell_types, ...)
!      m = mesh_from_gmsh(file)
!      w = paraview_writer(m, vertices, cell_vertices, &
!                          num_cell_vertices, cell_types)
!
! THE HYPERCUBE. Our own convention for a cell in any dimension: a
! box of 2^ndim corners in tensor order, the first coordinate
! varying fastest - corner k (from 0) is at the upper end of axis a
! when bit a-1 of k is set. In two dimensions
!
!      3 ---- 4          corners  1 = (lo, lo)   2 = (hi, lo)
!      |      |                   3 = (lo, hi)   4 = (hi, hi)
!      1 ---- 2
!
! which is NOT paraview's quadrangle cycle; the writer reorders. A
! mapped (curved) box is still a hypercube here: the order is the
! parametric one, the coordinates are wherever the map put them.
!
! THE DIMENSION. Paraview draws three coordinates, so the writer
! draws the axes asked for
!
!      axes  = [1, 2, 3]     the default: the first three, or fewer
!      axes  = [1, 3]        the xz view of a three-dimensional mesh
!      axes  = [3, 1]        x from the third coordinate, y from the first
!
! and the coordinates not drawn are either dropped - a projection,
! every cell drawn, flat - or fixed - a slice
!
!      slice = [0.5]                one value per coordinate not
!      slice = [0.25, 0.75]         drawn, in coordinate order
!
! A slice keeps the cells whose extent along each fixed coordinate
! contains its value (the extent is closed, so a slice on a face
! keeps the cells on both sides) and draws each one's cross-section
! there: for a four-dimensional mesh, axes = [1,2,3] with one value
! draws hexahedra, axes = [1,2] with two values draws quadrangles,
! one axis with three values draws lines. The cross-section is
! computed for hypercube cells only - the hypercube one dimension
! down, its corners interpolated between the two opposite faces at
! the parameter the value has in the cell's extent, exact when the
! cell is a box and first order when it is mapped. Any other cell
! type has a clipped polytope for a section, which this writer does
! not compute, so a slice through it is refused; a hypercube above
! three dimensions is drawn only through a slice.
!
!      w = paraview_writer(m, coordinates, cell_vertices, &
!                          num_cell_vertices, cell_types, &
!                          axes=[1, 2, 3], slice=[0.5_dp])
!
! THE GATE. Every array must fit the mesh's cells and the points,
! the drawn axes must be distinct coordinates, a slice must fix
! every coordinate not drawn, and a hypercube must own 2^ndim
! corners; the constructor stops the program otherwise.
!
! THE FILE. Ascii .vtu, points always three wide, Float64 cell data
! (the cell volume of the mesh, then the labelled fields), the
! drawn cells only - a slice writes fewer cells than the mesh has,
! each with its own points.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_paraview_writer

  use iso_fortran_env, only : error_unit, int32
  use util_precision , only : dp
  use util_string    , only : string
  use field_stored   , only : stored_field
  use view_mesh      , only : mesh

  implicit none

  private
  public :: paraview_writer, linear_cell_type
  public :: polygon_cell, hypercube_cell

  ! our own conventions, gmsh having no such elements
  integer, parameter :: polygon_cell   = -1  ! an agglomerated polygon
  integer, parameter :: hypercube_cell = -2  ! 2^ndim corners in tensor order

  !===================================================================!
  ! This datatype enumerates the paraview cell types.
  !===================================================================!

  type :: linear_cell_type

     integer(kind=int32) :: VTK_EMPTY_CELL        = 0
     integer(kind=int32) :: VTK_VERTEX            = 1
     integer(kind=int32) :: VTK_POLY_VERTEX       = 2
     integer(kind=int32) :: VTK_LINE              = 3
     integer(kind=int32) :: VTK_POLY_LINE         = 4
     integer(kind=int32) :: VTK_TRIANGLE          = 5
     integer(kind=int32) :: VTK_TRIANGLE_STRIP    = 6
     integer(kind=int32) :: VTK_POLYGON           = 7
     integer(kind=int32) :: VTK_PIXEL             = 8
     integer(kind=int32) :: VTK_QUAD              = 9
     integer(kind=int32) :: VTK_TETRA             = 10
     integer(kind=int32) :: VTK_VOXEL             = 11
     integer(kind=int32) :: VTK_HEXAHEDRON        = 12
     integer(kind=int32) :: VTK_WEDGE             = 13
     integer(kind=int32) :: VTK_PYRAMID           = 14
     integer(kind=int32) :: VTK_PENTAGONAL_PRISM  = 15
     integer(kind=int32) :: VTK_HEXAGONAL_PRISM   = 16
     integer(kind=int32) :: VTK_POLYHEDRON        = 42

   contains

     procedure :: element_type
     procedure :: hypercube_type

  end type linear_cell_type

  !===================================================================!
  ! This datatype writes a mesh and its cell fields to paraview. It
  ! holds what is drawn: three coordinates per point, each drawn
  ! cell's points in order, ragged, its paraview type, and the mesh
  ! cell it draws - a slice draws fewer cells than the mesh has.
  !===================================================================!

  type :: paraview_writer

     real(dp), allocatable :: points(:,:)
     integer , allocatable :: first_point(:)
     integer , allocatable :: cell_points(:)
     integer , allocatable :: types(:)
     integer , allocatable :: cells(:)

     ! the mesh's own count and volumes, one per mesh cell
     integer               :: num_cells = 0
     real(dp), allocatable :: volumes(:)

     type(linear_cell_type) :: cell_type

   contains

     procedure :: write

  end type paraview_writer

  !===================================================================!
  ! This interface admits multiple constructors.
  !===================================================================!

  interface paraview_writer
     module procedure construct
  end interface paraview_writer

contains

  !===================================================================!
  ! This function maps gmsh element numbers to paraview cell types.
  !===================================================================!

  pure elemental type(integer) function element_type(this, gmsh_type) &
       & result (paraview_type)

    class(linear_cell_type), intent(in) :: this
    integer                , intent(in) :: gmsh_type

    select case (gmsh_type)
    case (1) ! A 2-node line.
       paraview_type = this % VTK_LINE
    case (2) ! A 3-node triangle.
       paraview_type = this % VTK_TRIANGLE
    case (3) ! A 4-node quadrangle.
       paraview_type = this % VTK_QUAD
    case (4) ! A 4-node tetrahedron.
       paraview_type = this % VTK_TETRA
    case (5) ! An 8-node hexahedron.
       paraview_type = this % VTK_HEXAHEDRON
    case (6) ! A 6-node prism.
       paraview_type = this % VTK_WEDGE
    case (7) ! A 5-node prism (a pyramid).
       paraview_type = this % VTK_PYRAMID
    case (polygon_cell)
       paraview_type = this % VTK_POLYGON
    case default
       paraview_type = this % VTK_POLYHEDRON
    end select

  end function element_type

  !===================================================================!
  ! The paraview type of a hypercube of n axes: a point, a line, a
  ! quadrangle, a hexahedron. Paraview draws nothing above three.
  !===================================================================!

  type(integer) function hypercube_type(this, n) result (paraview_type)

    class(linear_cell_type), intent(in) :: this
    integer                , intent(in) :: n

    select case (n)
    case (0)
       paraview_type = this % VTK_VERTEX
    case (1)
       paraview_type = this % VTK_LINE
    case (2)
       paraview_type = this % VTK_QUAD
    case (3)
       paraview_type = this % VTK_HEXAHEDRON
    case default
       call gate(.false., 'a hypercube of at most three drawn axes')
       paraview_type = this % VTK_EMPTY_CELL
    end select

  end function hypercube_type

  !===================================================================!
  ! This is the constructor for the paraview writer. The gate: the
  ! corner arrays must fit the mesh's cells, the drawn axes must be
  ! distinct coordinates, and a slice must fix every other one.
  !===================================================================!

  impure type(paraview_writer) function construct(m, coordinates, &
       & cell_vertices, num_cell_vertices, cell_types, axes, slice) &
       & result (this)

    type(mesh), intent(in)           :: m
    real(dp)  , intent(in)           :: coordinates(:,:)    ! (ndim, point)
    integer   , intent(in)           :: cell_vertices(:,:)  ! (corner, cell), padded
    integer   , intent(in)           :: num_cell_vertices(:)
    integer   , intent(in)           :: cell_types(:)
    integer   , intent(in), optional :: axes(:)
    real(dp)  , intent(in), optional :: slice(:)

    type(stored_field)   :: volume
    integer              :: ndim, num_points, c, k, n
    integer, allocatable :: drawn(:), fixed(:)
    logical              :: sliced

    ndim             = size(coordinates, 1)
    num_points       = size(coordinates, 2)
    this % num_cells = m % num_vertices()

    call gate(size(num_cell_vertices) == this % num_cells, 'one corner count per cell')
    call gate(size(cell_types)        == this % num_cells, 'one type per cell')
    call gate(size(cell_vertices, 2)  == this % num_cells, 'one corner list per cell')
    call gate(all(num_cell_vertices >= 1) .and. &
         &    all(num_cell_vertices <= size(cell_vertices, 1)), &
         &    'corner counts within the corner lists')
    do c = 1, this % num_cells
       n = num_cell_vertices(c)
       call gate(all(cell_vertices(1:n, c) >= 1) .and. &
            &    all(cell_vertices(1:n, c) <= num_points), 'corners among the points')
       if (cell_types(c) == hypercube_cell) then
          call gate(n == 2 ** ndim, 'a hypercube of 2^ndim corners')
       end if
    end do

    ! the axes drawn, and the ones not drawn in coordinate order
    if (present(axes)) then
       drawn = axes
    else
       drawn = [(k, k = 1, min(ndim, 3))]
    end if
    call gate(size(drawn) >= 1 .and. size(drawn) <= 3, 'one to three axes drawn')
    call gate(all(drawn >= 1) .and. all(drawn <= ndim), 'drawn axes among the coordinates')
    do k = 2, size(drawn)
       call gate(all(drawn(1:k-1) /= drawn(k)), 'drawn axes distinct')
    end do
    fixed = pack([(k, k = 1, ndim)], [(all(drawn /= k), k = 1, ndim)])

    sliced = .false.
    if (present(slice)) then
       call gate(size(slice) == size(fixed), 'one slice value per coordinate not drawn')
       sliced = size(fixed) >= 1
    end if

    volume = m % cell_volume()
    call volume % real_vector(this % volumes)
    call gate(size(this % volumes) == this % num_cells, 'one volume per cell from the mesh')

    if (sliced) then
       call gate(all(cell_types == hypercube_cell), 'hypercube cells for a slice')
       call sectioned(this, coordinates, cell_vertices, drawn, fixed, slice)
    else
       call projected(this, coordinates, cell_vertices, num_cell_vertices, &
            & cell_types, drawn)
    end if

  end function construct

  !===================================================================!
  ! The projection: every point, its drawn coordinates first and
  ! zeros after; every cell, its corners as given, a hypercube's
  ! reordered to paraview's order along the drawn axes.
  !===================================================================!

  subroutine projected(this, coordinates, cell_vertices, num_cell_vertices, &
       & cell_types, drawn)

    type(paraview_writer), intent(inout) :: this
    real(dp)             , intent(in)    :: coordinates(:,:)
    integer              , intent(in)    :: cell_vertices(:,:)
    integer              , intent(in)    :: num_cell_vertices(:)
    integer              , intent(in)    :: cell_types(:)
    integer              , intent(in)    :: drawn(:)

    integer :: ndim, num_points, c, k, n, at

    ndim       = size(coordinates, 1)
    num_points = size(coordinates, 2)

    allocate(this % points(3, num_points))
    this % points = 0.0_dp
    do k = 1, size(drawn)
       this % points(k, :) = coordinates(drawn(k), :)
    end do

    this % cells = [(c, c = 1, this % num_cells)]
    allocate(this % first_point(this % num_cells + 1))
    allocate(this % types(this % num_cells))
    allocate(this % cell_points(sum(num_cell_vertices)))

    this % first_point(1) = 1
    do c = 1, this % num_cells
       n  = num_cell_vertices(c)
       at = this % first_point(c)
       if (cell_types(c) == hypercube_cell) then
          call gate(ndim <= 3, 'a hypercube above three dimensions drawn through a slice')
          this % cell_points(at : at + n - 1) = &
               & cell_vertices(hypercube_order(ndim, [(k, k = 1, ndim)], drawn), c)
          this % types(c) = this % cell_type % hypercube_type(ndim)
       else
          this % cell_points(at : at + n - 1) = cell_vertices(1:n, c)
          this % types(c) = this % cell_type % element_type(cell_types(c))
       end if
       this % first_point(c + 1) = at + n
    end do

  end subroutine projected

  !===================================================================!
  ! The slice: the cells whose extent along every fixed axis contains
  ! the value there, each drawn as its cross-section - a hypercube of
  ! the drawn axes with its own points, shared with no other cell.
  !===================================================================!

  subroutine sectioned(this, coordinates, cell_vertices, drawn, fixed, slice)

    type(paraview_writer), intent(inout) :: this
    real(dp)             , intent(in)    :: coordinates(:,:)
    integer              , intent(in)    :: cell_vertices(:,:)
    integer              , intent(in)    :: drawn(:)
    integer              , intent(in)    :: fixed(:)
    real(dp)             , intent(in)    :: slice(:)

    integer :: ndim, num_corners, num_drawn_corners, num_kept, c, i, j, k
    integer , allocatable :: kept(:)
    real(dp), allocatable :: corners(:,:), section(:,:)

    ndim              = size(coordinates, 1)
    num_corners       = 2 ** ndim
    num_drawn_corners = 2 ** size(drawn)

    ! the kept cells first, so the arrays are sized once
    allocate(kept(this % num_cells))
    num_kept = 0
    do c = 1, this % num_cells
       corners = coordinates(:, cell_vertices(1:num_corners, c))
       if (within(corners, fixed, slice)) then
          num_kept       = num_kept + 1
          kept(num_kept) = c
       end if
    end do

    allocate(this % points(3, num_drawn_corners * num_kept))
    allocate(this % cell_points(num_drawn_corners * num_kept))
    allocate(this % first_point(num_kept + 1))
    allocate(this % types(num_kept))
    this % points = 0.0_dp
    this % cells  = kept(1:num_kept)

    do i = 1, num_kept
       c       = kept(i)
       corners = coordinates(:, cell_vertices(1:num_corners, c))
       section = sectioned_corners(corners, fixed, slice, drawn)
       do j = 1, num_drawn_corners
          k = (i - 1) * num_drawn_corners + j
          this % points(1:size(drawn), k) = section(drawn, j)
          this % cell_points(k)           = k
       end do
       this % first_point(i) = (i - 1) * num_drawn_corners + 1
       this % types(i)       = this % cell_type % hypercube_type(size(drawn))
    end do
    this % first_point(num_kept + 1) = num_drawn_corners * num_kept + 1

  end subroutine sectioned

  !===================================================================!
  ! Whether a cell's extent along every fixed axis contains the slice
  ! value there; the extent is closed, so a slice on a face keeps the
  ! cells on both sides.
  !===================================================================!

  pure logical function within(corners, fixed, slice)

    real(dp), intent(in) :: corners(:,:)
    integer , intent(in) :: fixed(:)
    real(dp), intent(in) :: slice(:)

    integer :: f

    within = .true.
    do f = 1, size(fixed)
       within = within .and. &
            & minval(corners(fixed(f), :)) <= slice(f) .and. &
            & slice(f) <= maxval(corners(fixed(f), :))
    end do

  end function within

  !===================================================================!
  ! The cross-section of one hypercube: each fixed axis is cut in
  ! turn, and the corners that remain are put in paraview's order
  ! along the drawn axes.
  !===================================================================!

  pure function sectioned_corners(corners, fixed, slice, drawn) result (section)

    real(dp), intent(in)  :: corners(:,:)
    integer , intent(in)  :: fixed(:)
    real(dp), intent(in)  :: slice(:)
    integer , intent(in)  :: drawn(:)
    real(dp), allocatable :: section(:,:)

    integer, allocatable :: carried(:)
    integer :: f, j, k

    section = corners
    carried = [(k, k = 1, size(corners, 1))]

    do f = 1, size(fixed)
       j       = findloc(carried, fixed(f), dim=1)
       section = cut(section, fixed(f), j, slice(f))
       carried = pack(carried, carried /= fixed(f))
    end do

    section = section(:, hypercube_order(size(carried), carried, drawn))

  end function sectioned_corners

  !===================================================================!
  ! One cut: the hypercube's corners, in tensor order with bit j-1
  ! running along the axis cut, become the hypercube one dimension
  ! down, each corner interpolated between the two faces at the
  ! parameter the slice value has in the cell's extent. A cell flat
  ! along the axis is its own section.
  !===================================================================!

  pure function cut(corners, axis, j, value) result (section)

    real(dp), intent(in)  :: corners(:,:)
    integer , intent(in)  :: axis, j
    real(dp), intent(in)  :: value
    real(dp), allocatable :: section(:,:)

    real(dp) :: lo, hi, t
    integer  :: bit, k, low, k0, k1

    lo = minval(corners(axis, :))
    hi = maxval(corners(axis, :))
    t  = 0.0_dp
    if (hi > lo) t = (value - lo) / (hi - lo)

    bit = 2 ** (j - 1)
    allocate(section(size(corners, 1), size(corners, 2) / 2))
    do k = 0, size(section, 2) - 1
       low = iand(k, bit - 1)
       k0  = low + 2 * (k - low)
       k1  = k0 + bit
       section(:, k + 1) = (1.0_dp - t) * corners(:, k0 + 1) + t * corners(:, k1 + 1)
    end do

  end function cut

  !===================================================================!
  ! The permutation from tensor order over the carried axes to
  ! paraview's order: the corners re-indexed with the drawn axes as
  ! the leading bits, in the order drawn, then paraview's cycle of
  ! the quadrangle and of the hexahedron's two faces.
  !===================================================================!

  pure function hypercube_order(n, carried, drawn) result (perm)

    integer, intent(in)  :: n
    integer, intent(in)  :: carried(:), drawn(:)
    integer, allocatable :: perm(:)

    integer, allocatable :: order(:), tensor(:), ordering(:)
    integer :: c, q, r, k

    ! the drawn axes that are carried, then the rest of the carried
    order = pack(drawn, [(any(carried == drawn(k)), k = 1, size(drawn))])
    order = [order, pack(carried, [(all(drawn /= carried(k)), k = 1, size(carried))])]

    allocate(tensor(2 ** n))
    do q = 0, 2 ** n - 1
       r = 0
       do c = 1, n
          if (btest(q, c - 1)) r = ibset(r, findloc(carried, order(c), dim=1) - 1)
       end do
       tensor(q + 1) = r + 1
    end do

    select case (n)
    case (2)
       ordering = [1, 2, 4, 3]
    case (3)
       ordering = [1, 2, 4, 3, 5, 6, 8, 7]
    case default
       ordering = [(k, k = 1, 2 ** n)]
    end select

    perm = tensor(ordering)

  end function hypercube_order

  !===================================================================!
  ! Write the mesh and the optional cell fields to a .vtu file in the
  ! UnstructuredGrid format.
  !===================================================================!

  impure subroutine write(this, filename, phic, solution_labels)

    ! These are the arguments.
    class(paraview_writer) , intent(in)           :: this
    character(len=*)       , intent(in)           :: filename
    real(dp)               , intent(in), optional :: phic(:,:) ! (icell, ivar)
    type(string)           , intent(in), optional :: solution_labels(:)

    ! These are the locals.
    integer :: ierr
    integer :: fhandle
    integer :: iresult
    integer :: num_drawn

    num_drawn = size(this % cells)

    if (present(solution_labels) .and. present(phic)) then
       call gate(size(phic, 1) == this % num_cells, 'one value per mesh cell per label')
       call gate(size(phic, 2) == size(solution_labels), 'one label per field')
    end if

    ! Open the output file for formatted writing.
    open(newunit=fhandle, file=trim(filename), iostat=ierr, action='write', &
         & form='formatted', status='replace')
    if (ierr .ne. 0) then
       write(error_unit, '(a)') 'paraview writer: opening ' // trim(filename) // ' failed'
       error stop 'paraview writer: the file cannot be written'
    end if

    !-----------------------------------------------------------------!
    ! Write the basic header information.
    !-----------------------------------------------------------------!

    write(fhandle, '(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(fhandle, '(a)') '<UnstructuredGrid>'
    write(fhandle, '(a,i0,a,i0,a)') '<Piece NumberOfPoints="', size(this % points, 2), &
         & '" NumberOfCells="', num_drawn, '">'

    !-----------------------------------------------------------------!
    ! Write the vertices.
    !-----------------------------------------------------------------!

    write_points: block

      integer :: ivertex

      write(fhandle, '(a)') '<Points>'
      write(fhandle, '(a)') '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do ivertex = 1, size(this % points, 2)
         write(fhandle, '(3(1x,es24.16))') this % points(:, ivertex)
      end do
      write(fhandle, '(a)') '</DataArray>'
      write(fhandle, '(a)') '</Points>'

    end block write_points

    write_cells: block

      integer :: icell, jvertex

      write(fhandle, '(a)') '<Cells>'

      !---------------------------------------------------------------!
      ! Write the cell-to-vertex connectivities.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
      do icell = 1, num_drawn
         ! Correct for paraview's 0-based numbering.
         write(fhandle, '(*(1x,i0))') (this % cell_points(jvertex) - 1, &
              & jvertex = this % first_point(icell), this % first_point(icell + 1) - 1)
      end do
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! Write the cell-to-vertex connectivity offsets.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
      do icell = 1, num_drawn
         write(fhandle, '(i0)') this % first_point(icell + 1) - 1
      end do
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! Write the cell types.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
      do icell = 1, num_drawn
         write(fhandle, '(i0)') this % types(icell)
      end do
      write(fhandle, '(a)') '</DataArray>'

      write(fhandle, '(a)') '</Cells>'

      write(fhandle, '(a)') '<PointData></PointData>'

      !---------------------------------------------------------------!
      ! Write the cell data.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<CellData>'

      ! Export the cell volumes.
      write(fhandle, '(a)') '<DataArray type="Float64" Name="volume" format="ascii">'
      do icell = 1, num_drawn
         write(fhandle, '(es24.16)') this % volumes(this % cells(icell))
      end do
      write(fhandle, '(a)') '</DataArray>'

      if (present(solution_labels) .and. present(phic)) then
         do iresult = 1, size(solution_labels)
            write(fhandle, '(a,a,a)') '<DataArray type="Float64" Name="', &
                 & solution_labels(iresult) % str, '" format="ascii">'
            do icell = 1, num_drawn
               write(fhandle, '(es24.16)') phic(this % cells(icell), iresult)
            end do
            write(fhandle, '(a)') '</DataArray>'
         end do
      end if

      write(fhandle, '(a)') '</CellData>'

    end block write_cells

    ! Close the opened tags.
    write(fhandle, '(a)') '</Piece>'
    write(fhandle, '(a)') '</UnstructuredGrid>'
    write(fhandle, '(a)') '</VTKFile>'

    close(unit=fhandle)

  end subroutine write

  !===================================================================!
  ! The gate itself: state what failed, then stop. A picture of the
  ! wrong cells must not be written.
  !===================================================================!

  subroutine gate(fits, what)

    logical         , intent(in) :: fits
    character(len=*), intent(in) :: what

    if (fits) return

    write(error_unit, *) 'paraview writer gate: expected ', what
    error stop 'paraview writer: the corners do not fit the mesh'

  end subroutine gate

end module view_paraview_writer
