!=====================================================================!
! Module that reads the surface of a Nastran bulk data file, for
! drawing: the points of its GRID cards, of short fields or of long
! ones (GRID*, the third coordinate on the continuation line), and the
! corners of its quadrangles and triangles (CQUAD4, CQUAD8, CQUAD9,
! CTRIA3, CTRIA6), the corners first on each card:
!
!   short fields, 8 columns each
!   |GRID    |ID      |CP      |X1      |X2      |X3      |
!    1        9        17       25       33       41       49
!
!   long fields, 16 columns each, X3 on the continuation line
!   |GRID*   |ID              |CP              |X1              |X2              |
!    1        9                25               41               57               73
!   |*       |X3              |
!    1        9                25
!
!   a cell's corners, from the fourth field
!   |CQUAD4  |EID     |PID     |G1      |G2      |G3      |G4      |
!    1        9        17       25       33       41       49       57
!
! Each cell is drawn by its corners. A point no corner names is
! dropped, and the points are renumbered in the order the cells'
! corners first name them. Any other card is passed over.
!
!      call bdf_surface('hub4b.bdf', points, corners, count, kinds)
!
! points(:, p) stores the three coordinates of point p; corners(:, c)
! the points of cell c, padded to four; count(c) its number of
! corners; kinds(c) its gmsh element number, 3 the quadrangle and 2
! the triangle. A number in the short fields may omit the exponent's
! letter, 1.5-3 for 1.5e-3, as Nastran admits.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_bdf_loader

  use util_precision, only : dp

  implicit none

  private

  public :: bdf_surface

contains

  !===================================================================!
  ! Reads the points and the cells of a bulk data file, and renumbers
  ! the points the cells name:
  !
  !   ids read          101   102   205   300
  !   place(id)          -1    -2    -3    -4    minus the read position
  !   cells name        205, 101, 102
  !   place(id)           2     3     1    -4    the order first named;
  !                                              300 is dropped
  !===================================================================!

  subroutine bdf_surface(file, points, corners, count, kinds)

    character(len=*)     , intent(in)  :: file
    real(dp), allocatable, intent(out) :: points(:,:)
    integer , allocatable, intent(out) :: corners(:,:), count(:), kinds(:)

    character(len=128)    :: line, next
    real(dp), allocatable :: read_points(:,:)
    integer , allocatable :: ids(:), cells(:,:), cell_count(:), cell_kind(:), place(:)
    integer               :: unit, status, np, nc, id, c, k, n

    ! The cards, read into arrays that double as they fill
    allocate(read_points(3, 1024), ids(1024), cells(4, 1024), cell_count(1024), cell_kind(1024))
    np = 0
    nc = 0
    open(newunit=unit, file=file, status='old', action='read', iostat=status)
    if (status /= 0) error stop 'view_bdf_loader: the bulk data file cannot be read: ' // trim(file)
    do
       read(unit, '(a)', iostat=status) line
       if (status /= 0) exit
       if (line(1:5) == 'GRID*') then
          read(unit, '(a)', iostat=status) next
          if (status /= 0) error stop 'view_bdf_loader: a GRID* card without its continuation in ' // trim(file)
          call grow_points()
          np = np + 1
          read(line(9:24), *) ids(np)
          read_points(1, np) = number(line(41:56))
          read_points(2, np) = number(line(57:72))
          read_points(3, np) = number(next(9:24))
       else if (line(1:5) == 'GRID ' .or. line(1:8) == 'GRID    ') then
          call grow_points()
          np = np + 1
          read(line(9:16), *) ids(np)
          read_points(1, np) = number(line(25:32))
          read_points(2, np) = number(line(33:40))
          read_points(3, np) = number(line(41:48))
       else if (line(1:8) == 'CQUAD4  ' .or. line(1:8) == 'CQUAD8  ' .or. line(1:8) == 'CQUAD9  ') then
          call grow_cells()
          nc = nc + 1
          do k = 1, 4
             read(line(17 + 8 * k:24 + 8 * k), *) cells(k, nc)
          end do
          cell_count(nc) = 4
          cell_kind(nc)  = 3
       else if (line(1:8) == 'CTRIA3  ' .or. line(1:8) == 'CTRIA6  ') then
          call grow_cells()
          nc = nc + 1
          do k = 1, 3
             read(line(17 + 8 * k:24 + 8 * k), *) cells(k, nc)
          end do
          cells(4, nc)   = 0
          cell_count(nc) = 3
          cell_kind(nc)  = 2
       end if
    end do
    close(unit)

    ! The points the corners name, renumbered in the order they are
    ! first named
    allocate(place(minval(ids(1:np)):maxval(ids(1:np))), source=0)
    do k = 1, np
       place(ids(k)) = -k
    end do
    do c = 1, nc
       do k = 1, cell_count(c)
          if (cells(k, c) < lbound(place, 1) .or. cells(k, c) > ubound(place, 1)) then
             error stop 'view_bdf_loader: a cell names a point the file does not give: ' // trim(file)
          end if
          if (place(cells(k, c)) == 0) error stop 'view_bdf_loader: a cell names a point the file does not give: ' &
               & // trim(file)
       end do
    end do
    n = 0
    allocate(points(3, np))
    do c = 1, nc
       do k = 1, cell_count(c)
          id = cells(k, c)
          if (place(id) < 0) then
             n = n + 1
             points(:, n) = read_points(:, -place(id))
             place(id) = n
          end if
          cells(k, c) = place(id)
       end do
    end do
    points  = points(:, 1:n)
    corners = cells(:, 1:nc)
    count   = cell_count(1:nc)
    kinds   = cell_kind(1:nc)

  contains

    ! The points' arrays doubled when full
    subroutine grow_points()
      real(dp), allocatable :: p(:,:)
      integer , allocatable :: i(:)
      if (np < size(ids)) return
      allocate(p(3, 2 * size(ids)), i(2 * size(ids)))
      p(:, 1:np) = read_points(:, 1:np)
      i(1:np)    = ids(1:np)
      call move_alloc(p, read_points)
      call move_alloc(i, ids)
    end subroutine grow_points

    ! The cells' arrays doubled when full
    subroutine grow_cells()
      integer, allocatable :: q(:,:), m(:), t(:)
      if (nc < size(cell_count)) return
      allocate(q(4, 2 * size(cell_count)), m(2 * size(cell_count)), t(2 * size(cell_count)))
      q(:, 1:nc) = cells(:, 1:nc)
      m(1:nc)    = cell_count(1:nc)
      t(1:nc)    = cell_kind(1:nc)
      call move_alloc(q, cells)
      call move_alloc(m, cell_count)
      call move_alloc(t, cell_kind)
    end subroutine grow_cells

  end subroutine bdf_surface

  !===================================================================!
  ! Returns the number in a field, its exponent's letter supplied where
  ! Nastran omits it: a sign after the first digit begins the exponent,
  !
  !   '1.5-3'   ->  '1.5e-3'
  !   '-2.+4'   ->  '-2.e+4'
  !===================================================================!

  real(dp) function number(field)

    character(len=*), intent(in) :: field

    character(len=len(field) + 1) :: text
    integer                       :: k, first

    text  = adjustl(field)
    first = 1
    if (text(1:1) == '+' .or. text(1:1) == '-') first = 2
    if (scan(text, 'eEdD') == 0) then
       do k = first + 1, len_trim(text)
          if (text(k:k) == '+' .or. text(k:k) == '-') then
             text = text(1:k - 1) // 'e' // text(k:)
             exit
          end if
       end do
    end if
    read(text, *) number

  end function number

end module view_bdf_loader
