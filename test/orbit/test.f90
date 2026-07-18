!=====================================================================!
! The escape-time painting, on a CFD mesh, through the graph's orbit.
!
! The idea is simple. Every cell's centroid is a complex number once
! the mesh is stretched over a window of the complex plane. One step
! of the julia map z -> z**2 + c moves that number somewhere: either
! off the window (the orbit escapes) or into some other cell. So the
! map, quantized to the mesh, gives every cell one arrow out - a
! stored successor array, built once. That array IS the functional
! graph; after that there is no arithmetic left, only following
! arrows: each cell's orbit is walked by the graph's orbit method,
! and the cell is painted by how long its orbit lasted. Cells whose
! orbits escape get their escape time; cells still home when patience
! runs out (or trapped on a cycle) are painted the step limit - the
! dark interior of the set. The painting goes out as a paraview file.
!
! The julia constant is the douady rabbit's, which sits inside the
! mandelbrot set, so the painting is guaranteed both light (escapers
! near the corners) and dark (the rabbit's body).
!
! Checks:
!   1. the painting matches an independent recount of the arrows (a
!      plain loop, no orbit machinery, repaints the same picture)
!   2. the painting has light and dark (escaped and home both nonzero)
!   3. the window's corner escapes at once (its first step leaves)
!   4. the painting is on disk
!=====================================================================!

program test_orbit_painting

  use iso_fortran_env      , only : dp => REAL64
  use class_mesh           , only : mesh
  use class_gmsh_loader    , only : gmsh_loader
  use class_graph          , only : mesh_graph
  use class_paraview_writer, only : paraview_writer
  use class_string         , only : string

  implicit none

  character(len=*), parameter :: mesh_file  = '../square-80.msh'
  character(len=*), parameter :: paint_file = 'julia.vtu'

  ! the julia constant (douady rabbit) and the window of the complex
  ! plane the mesh is stretched over
  complex(dp), parameter :: julia_constant    = (-0.123_dp, 0.745_dp)
  real(dp)   , parameter :: window_half_width = 1.6_dp
  real(dp)   , parameter :: escape_radius     = 2.0_dp
  integer    , parameter :: step_limit        = 40   ! the patience

  class(gmsh_loader)    , allocatable :: loader
  class(mesh)           , allocatable :: grid
  class(paraview_writer), allocatable :: painter
  type(mesh_graph) :: g

  complex(dp), allocatable :: zc(:)     ! cell centroids on the complex plane
  integer    , allocatable :: succ(:)   ! the one arrow out of each cell (0 = off the window)
  real(dp)   , allocatable :: escape_time(:,:)
  integer    , allocatable :: visited(:)
  type(string) :: field_labels(1)

  integer :: ncells, icell, n_escaped, n_home, nfail

  nfail = 0

  meshing: block

    allocate(loader, source = gmsh_loader(mesh_file))
    allocate(grid  , source = mesh(loader))
    deallocate(loader)

    g      = mesh_graph(grid)
    ncells = g % num_vertices

  end block meshing

  !===================================================================!
  ! Stretch the mesh over the window: an affine map from the bounding
  ! box of the centroids onto the square of half-width 1.6 centred at
  ! the origin, where the rabbit lives.
  !===================================================================!

  stretch_over_the_plane: block

    real(dp) :: xlo, xhi, ylo, yhi

    xlo = minval(grid % cell_centers(1, :)); xhi = maxval(grid % cell_centers(1, :))
    ylo = minval(grid % cell_centers(2, :)); yhi = maxval(grid % cell_centers(2, :))

    allocate(zc(ncells))
    do icell = 1, ncells
       zc(icell) = cmplx( &
            & window_half_width*(2.0_dp*(grid % cell_centers(1, icell) - xlo)/(xhi - xlo) - 1.0_dp), &
            & window_half_width*(2.0_dp*(grid % cell_centers(2, icell) - ylo)/(yhi - ylo) - 1.0_dp), &
            & kind = dp)
    end do

  end block stretch_over_the_plane

  !===================================================================!
  ! Quantize the julia map to the mesh: one step from each centroid,
  ! then the nearest centroid claims the landing point. Past the
  ! escape radius nobody claims it - the arrow leaves the window and
  ! the successor is 0. Built once; the dynamics is now pure graph.
  !===================================================================!

  build_the_arrows: block

    complex(dp) :: z
    real(dp)    :: best, d2
    integer     :: v, w, nearest

    allocate(succ(ncells))

    do v = 1, ncells

       z = zc(v)*zc(v) + julia_constant

       if (abs(z) .gt. escape_radius) then
          succ(v) = 0
          cycle
       end if

       nearest = 1
       best    = huge(1.0_dp)
       do w = 1, ncells
          d2 = (real(z, dp) - real(zc(w), dp))**2 + (aimag(z) - aimag(zc(w)))**2
          if (d2 .lt. best) then
             best    = d2
             nearest = w
          end if
       end do
       succ(v) = nearest

    end do

  end block build_the_arrows

  !===================================================================!
  ! Paint: walk every cell's orbit and colour by how long it lasted.
  ! Escaped (the final cell's arrow points off the window): the escape
  ! time. Still home at the limit, or trapped on a cycle: the full
  ! step limit - the dark interior.
  !===================================================================!

  paint_by_patience: block

    allocate(escape_time(ncells, 1))

    n_escaped = 0
    n_home    = 0

    do icell = 1, ncells
       visited = g % orbit(icell, successor, step_limit)
       if (succ(visited(size(visited))) .eq. 0) then
          escape_time(icell, 1) = real(size(visited), dp)
          n_escaped = n_escaped + 1
       else
          escape_time(icell, 1) = real(step_limit, dp)
          n_home = n_home + 1
       end if
    end do

  end block paint_by_patience

  field_labels(1) = string('escape_time')
  allocate(painter, source = paraview_writer(grid))
  call painter % write(paint_file, escape_time, field_labels)

  write(*,'(1x,a,i0,a,i0,a,i0,a)') "painted ", ncells, " cells: ", &
       & n_escaped, " escaped, ", n_home, " home"

  call check_the_painting(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all orbit painting checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " orbit painting check(s)"
     error stop
  end if

contains

  ! the single successor of cell v, read from the stored arrows
  pure integer function successor(v)
    integer, intent(in) :: v
    successor = succ(v)
  end function successor

  subroutine report(ok, label, nfail)
    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

  subroutine check_the_painting(nfail)

    integer, intent(inout) :: nfail
    logical :: on_disk

    ! 1: repaint every cell by following the stored arrows in a plain
    ! loop - no orbit machinery - and demand the same picture. A cycle
    ! just spins here until the limit, which is the same "home" the
    ! orbit's early cycle exit paints.
    recount_the_painting: block

      real(dp) :: repaint
      logical  :: agrees
      integer  :: v, steps

      agrees = .true.
      do icell = 1, ncells
         v     = icell
         steps = 0
         do
            steps = steps + 1
            if (succ(v) .eq. 0) then
               repaint = real(steps, dp)          ! escaped at this length
               exit
            end if
            if (steps .eq. step_limit) then
               repaint = real(step_limit, dp)     ! still home
               exit
            end if
            v = succ(v)
         end do
         if (repaint .ne. escape_time(icell, 1)) agrees = .false.
      end do

      call report(agrees, &
           & "painting matches an independent recount of the arrows", nfail)

    end block recount_the_painting

    call report(n_escaped .gt. 0 .and. n_home .gt. 0, &
         & "the painting has light and dark", nfail)

    ! 3: the cell nearest the window's corner sits far outside the
    ! set - its very first step must leave, so it is painted 1
    corner_escapes: block

      complex(dp) :: corner_z
      real(dp)    :: d2, best
      integer     :: v, corner

      corner_z = cmplx(-window_half_width, -window_half_width, kind = dp)
      corner   = 1
      best     = huge(1.0_dp)
      do v = 1, ncells
         d2 = abs(zc(v) - corner_z)**2
         if (d2 .lt. best) then
            best   = d2
            corner = v
         end if
      end do

      call report(nint(escape_time(corner, 1)) .eq. 1, &
           & "the window's corner escapes at once", nfail)

    end block corner_escapes

    inquire(file = paint_file, exist = on_disk)
    call report(on_disk, "the painting is on disk", nfail)

  end subroutine check_the_painting

end program test_orbit_painting
