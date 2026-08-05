!=====================================================================!
! The mesh builder: the gmsh path onto the tower.
!
! A gmsh file becomes a tower mesh in one call:
!
!      m = mesh_from_gmsh('square-10.msh')
!
! THE BRIDGE, STATED PLAINLY. Parsing the file and computing the
! measurements - volumes from coordinates, areas by triangulation,
! normals, deltas, weights - still live in the old world
! (class_gmsh_loader, class_mesh). This builder runs that machinery
! once, at load, and carries its answers into the tower's seat: the
! new mesh, whose measurements travel as typed fields from there on.
! Nothing downstream of this call touches the old world.
!
! The builder dies when the measurement machinery is ported into the
! tower; its call site will not change when that happens, which is
! the point of building it now.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_mesh_builder

  use iso_fortran_env      , only : dp => REAL64
  use class_graph_mesh     , only : mesh
  use class_mesh           , only : legacy => mesh
  use interface_mesh_loader, only : mesh_loader
  use class_gmsh_loader    , only : gmsh_loader

  implicit none

  private
  public :: mesh_from_gmsh

contains

  !===================================================================!
  ! Load, measure, seat. The old machinery parses and measures; the
  ! arrays it produces are rearranged into the seat's order: one
  ! entry per cell or face, vectors three wide, tails and heads from
  ! the face-to-cell map, tag names stamped on the boundary faces.
  !===================================================================!

  impure type(mesh) function mesh_from_gmsh(filename) result(m)

    character(len=*), intent(in) :: filename

    type(legacy) :: old
    class(mesh_loader), allocatable :: gl

    integer , allocatable :: tails(:), heads(:)
    real(dp), allocatable :: normals(:)
    real(dp), allocatable :: weights(:)
    character(len=64), allocatable :: etags(:)
    integer :: nc, nf, f, c, l, tag

    allocate(gl, source=gmsh_loader(trim(filename)))
    old = legacy(gl)

    nc = old % num_cells
    nf = old % num_faces

    ! Tails and heads from the face-to-cell map. A face with one
    ! cell is a boundary face: an edge without a head.
    allocate(tails(nf), heads(nf))
    do f = 1, nf
       tails(f) = old % face_cells(1, f)
       heads(f) = 0
       if (old % num_face_cells(f) >= 2) heads(f) = old % face_cells(2, f)
    end do

    ! One outward normal per face, read from its tail cell's side.
    allocate(normals(3 * nf))
    normals = 0.0_dp
    do c = 1, nc
       do l = 1, old % num_cell_faces(c)
          f = old % cell_faces(l, c)
          if (old % face_cells(1, f) == c) then
             normals(3 * f - 2 : 3 * f) = old % cell_face_normals(:, l, c)
          end if
       end do
    end do

    ! One interpolation weight per face: the tail cell's share.
    allocate(weights(nf))
    do f = 1, nf
       weights(f) = old % face_cell_weights(1, f)
    end do

    ! Tag names on the boundary faces, blank inside: the strings the
    ! mesh file itself declared, entering here once.
    allocate(etags(nf))
    etags = ''
    do f = 1, nf
       if (old % num_face_cells(f) < 2) then
          tag = old % face_tags(f)
          if (tag >= 1 .and. tag <= size(old % tag_info)) then
             etags(f) = old % tag_info(tag) % str
          end if
       end if
    end do

    m = mesh(nc, tails=tails, heads=heads, &
         & volumes      = old % cell_volumes, &
         & cell_centers = reshape(old % cell_centers, [3 * nc]), &
         & areas        = old % face_areas, &
         & deltas       = old % face_deltas, &
         & normals      = normals, &
         & face_centers = reshape(old % face_centers, [3 * nf]), &
         & weights      = weights, &
         & etags        = etags)

  end function mesh_from_gmsh

end module class_mesh_builder
