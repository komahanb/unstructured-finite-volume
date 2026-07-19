!=====================================================================!
! The escape-time painting, on a CFD mesh, through the graph's orbit -
! now with both fast lanes the squint provides.
!
! The idea is unchanged. Every cell's centroid is a complex number
! once the mesh is stretched over a window of the complex plane. One
! step of the julia map z -> z**2 + c moves that number somewhere:
! either off the window (the orbit escapes) or into some other cell.
! The map, quantized to the mesh, gives every cell one arrow out - a
! stored successor array. After that the dynamics is pure graph, and
! the cell is painted by how long its orbit lasts: escapers get their
! escape time, cells still home at the limit get the full patience.
!
! What is new is where the work went.
!
! Squint in space: finding which cell claims a landing point used to
! scan every centroid - the whole painting cost cells x cells
! distance checks. Now the mesh graph is partitioned into compact
! parts, the parts squint into a quotient graph, and a landing point
! is claimed by descending: nearest part first, then only that part's
! cells and its quotient-neighbours' cells. Same arrows, a fraction
! of the checks - and a check below proves the "same".
!
! Squint in time: orbits that merge share their tails, so the whole
! painting is resolved by escape_times in one pass - one visit per
! cell - instead of one orbit per cell. Two recounts prove it paints
! the same picture: the orbit method cell by cell, and a plain loop
! with no machinery at all.
!
! Checks:
!   1. the squint finds the same arrows as the full scan
!   2. the one-pass painting matches orbit-by-orbit painting
!   3. and matches a plain-loop recount with no machinery at all
!   4. the painting has light and dark (escaped and home both nonzero)
!   5. the window's corner escapes at once (its first step leaves)
!   6. the painting is on disk
!=====================================================================!

program test_orbit_painting

  use iso_fortran_env      , only : dp => REAL64
  use class_mesh           , only : mesh
  use class_gmsh_loader    , only : gmsh_loader
  use class_graph          , only : mesh_graph
  use class_stored_graph   , only : stored_graph
  use class_paraview_writer, only : paraview_writer
  use class_string         , only : string

  implicit none

  character(len=*), parameter :: mesh_file  = '../square-tri-40.msh'
  character(len=*), parameter :: paint_file = 'julia.vtu'

  ! the julia constant (douady rabbit) and the window of the complex
  ! plane the mesh is stretched over
  complex(dp), parameter :: julia_constant    = (-0.123_dp, 0.745_dp)
  real(dp)   , parameter :: window_half_width = 1.6_dp
  real(dp)   , parameter :: escape_radius     = 2.0_dp
  integer    , parameter :: step_limit        = 40   ! the patience
  integer    , parameter :: cells_per_part    = 64   ! spatial squint factor

  class(gmsh_loader)    , allocatable :: loader
  class(mesh)           , allocatable :: grid
  class(paraview_writer), allocatable :: painter
  type(mesh_graph)   :: g
  type(stored_graph) :: coarse

  complex(dp), allocatable :: zc(:)          ! cell centroids on the complex plane
  complex(dp), allocatable :: zpart(:)       ! part centroids on the complex plane
  real(dp)   , allocatable :: part_radius(:) ! how far a part's cells stray from its centroid
  integer    , allocatable :: succ(:)        ! the one arrow out of each cell (0 = off the window)
  integer    , allocatable :: times(:)       ! the painting, resolved in one pass
  real(dp)   , allocatable :: escape_time(:,:)
  type(string) :: field_labels(1)

  integer  :: ncells, nparts, icell, n_home
  integer  :: checks_by_descent, full_scans, nfail
  real(dp) :: t0, t1

  nfail = 0

  ! a stale painting from a previous run would satisfy the on-disk
  ! check without the writer doing anything - burn the old canvas first
  fresh_canvas: block
    integer :: unit, ierr
    open(newunit = unit, file = paint_file, status = 'old', iostat = ierr)
    if (ierr .eq. 0) close(unit, status = 'delete')
  end block fresh_canvas

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
  ! The spatial squint: partition the mesh graph into compact parts by
  ! coordinate bisection, squint the parts into the quotient graph,
  ! and place each part's centroid on the window - the coarse level a
  ! landing point descends through.
  !===================================================================!

  squint_the_mesh: block

    integer, allocatable :: members(:)
    integer :: k

    nparts = max(1, ncells/cells_per_part)
    call g % partition_rcb(grid % cell_centers, nparts)
    coarse = stored_graph(g)

    allocate(zpart(nparts), part_radius(nparts))
    do k = 1, nparts
       members        = g % owned(k)
       zpart(k)       = sum(zc(members))/real(size(members), dp)
       part_radius(k) = sqrt(maxval(abs(zc(members) - zpart(k))**2))
    end do

  end block squint_the_mesh

  !===================================================================!
  ! Quantize the julia map to the mesh: one step from each centroid,
  ! then the nearest centroid claims the landing point - found by
  ! descending the squint: nearest part first, then only that part's
  ! cells and its quotient-neighbours' cells. Past the escape radius
  ! nobody claims it - the arrow leaves the window, successor 0.
  !===================================================================!

  build_the_arrows: block

    complex(dp) :: z
    integer :: v

    allocate(succ(ncells))
    checks_by_descent = 0
    full_scans        = 0

    call cpu_time(t0)
    do v = 1, ncells
       z = zc(v)*zc(v) + julia_constant
       if (abs(z) .gt. escape_radius) then
          succ(v) = 0
       else
          succ(v) = nearest_cell_by_descent(z, checks_by_descent)
       end if
    end do
    call cpu_time(t1)

  end block build_the_arrows

  !===================================================================!
  ! Paint in one pass: escape_times resolves every cell at once -
  ! orbits that merge share their tails. Escapers carry their escape
  ! time; cells still home at the limit carry the full patience - the
  ! dark interior of the set.
  !===================================================================!

  times = g % escape_times(successor, step_limit)

  allocate(escape_time(ncells, 1))
  escape_time(:, 1) = real(times, dp)
  n_home = count(times .eq. step_limit)

  field_labels(1) = string('escape_time')
  allocate(painter, source = paraview_writer(grid))
  call painter % write(paint_file, escape_time, field_labels)

  write(*,'(1x,a,i0,a,i0,a,i0,a)') "painted ", ncells, " cells: ", &
       & ncells - n_home, " escaped, ", n_home, " home"
  write(*,'(1x,a,i0,a,i0,a,f6.3,a)') "the squint claimed the arrows with ", &
       & checks_by_descent, " distance checks where the full scan needs ", &
       & ncells*ncells, " (", t1 - t0, " s)"
  write(*,'(1x,a,i0,a)') "the certificate sent ", full_scans, &
       & " landing(s) back to the full scan"

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

  ! which cell claims the landing point z: descend the squint - the
  ! nearest part's cells and its quotient-neighbours' cells - then
  ! certify the answer by the triangle inequality: no cell of a far
  ! part can beat what we hold, because every one of its cells sits
  ! at least its centroid distance minus its radius away. On the rare
  ! landing the certificate cannot vouch for, fall back to the full
  ! scan - counted, so the suite reports how rare.
  integer function nearest_cell_by_descent(z, checks) result(nearest)

    complex(dp), intent(in)    :: z
    integer    , intent(inout) :: checks

    integer, allocatable :: parts(:)
    logical, allocatable :: candidate(:)
    real(dp), allocatable :: dpart(:)
    real(dp) :: best, d2
    logical  :: certified
    integer  :: k, k0, i, j, v

    ! nearest part centroid (every part's distance kept for the
    ! certificate below)
    allocate(dpart(nparts))
    k0   = 1
    best = huge(1.0_dp)
    do k = 1, nparts
       dpart(k) = abs(z - zpart(k))
       checks   = checks + 1
       if (dpart(k)**2 .lt. best) then
          best = dpart(k)**2
          k0   = k
       end if
    end do

    ! its cells and its quotient-neighbours' cells
    parts = [k0, coarse % neighbours(k0)]
    allocate(candidate(nparts))
    candidate = .false.
    candidate(parts) = .true.

    nearest = 0
    best    = huge(1.0_dp)
    do i = 1, size(parts)
       associate(members => g % owned(parts(i)))
         do j = 1, size(members)
            v  = members(j)
            d2 = abs(z - zc(v))**2
            checks = checks + 1
            if (d2 .lt. best) then
               best    = d2
               nearest = v
            end if
         end do
       end associate
    end do

    ! the certificate: every unsearched part must sit entirely
    ! farther away than the cell we hold
    certified = .true.
    do k = 1, nparts
       if (candidate(k)) cycle
       if (dpart(k) - part_radius(k) .lt. sqrt(best)) then
          certified = .false.
          exit
       end if
    end do

    if (.not. certified) then
       full_scans = full_scans + 1
       do v = 1, ncells
          d2 = abs(z - zc(v))**2
          checks = checks + 1
          if (d2 .lt. best) then
             best    = d2
             nearest = v
          end if
       end do
    end if

  end function nearest_cell_by_descent

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

    ! 1: the full scan - every centroid consulted for every landing
    ! point - must claim exactly the arrows the descent claimed
    full_scan_agrees: block

      complex(dp) :: z
      real(dp)    :: best, d2
      logical     :: agrees
      integer     :: v, w, nearest

      agrees = .true.
      do v = 1, ncells
         z = zc(v)*zc(v) + julia_constant
         if (abs(z) .gt. escape_radius) then
            nearest = 0
         else
            nearest = 1
            best    = huge(1.0_dp)
            do w = 1, ncells
               d2 = abs(z - zc(w))**2
               if (d2 .lt. best) then
                  best    = d2
                  nearest = w
               end if
            end do
         end if
         if (nearest .ne. succ(v)) agrees = .false.
      end do

      call report(agrees, "the squint finds the same arrows", nfail)

    end block full_scan_agrees

    ! 2: orbit-by-orbit painting - the orbit method walked from every
    ! cell, classified by whether its final arrow leaves - must paint
    ! the same picture the one-pass resolution painted
    orbit_recount: block

      integer, allocatable :: visited(:)
      logical :: agrees
      integer :: v, painted

      agrees = .true.
      do v = 1, ncells
         visited = g % orbit(v, successor, step_limit)
         if (succ(visited(size(visited))) .eq. 0) then
            painted = size(visited)
         else
            painted = step_limit
         end if
         if (painted .ne. times(v)) agrees = .false.
      end do

      call report(agrees, "one-pass painting matches orbit-by-orbit painting", nfail)

    end block orbit_recount

    ! 3: a plain loop with no machinery - follow the stored arrows,
    ! count, cap - must repaint the same picture too
    recount_the_painting: block

      logical :: agrees
      integer :: v, w, steps, repaint

      agrees = .true.
      do v = 1, ncells
         w     = v
         steps = 0
         do
            steps = steps + 1
            if (succ(w) .eq. 0) then
               repaint = steps               ! escaped at this length
               exit
            end if
            if (steps .eq. step_limit) then
               repaint = step_limit          ! still home
               exit
            end if
            w = succ(w)
         end do
         if (repaint .ne. times(v)) agrees = .false.
      end do

      call report(agrees, &
           & "painting matches an independent recount of the arrows", nfail)

    end block recount_the_painting

    call report(n_home .gt. 0 .and. n_home .lt. ncells, &
         & "the painting has light and dark", nfail)

    ! 5: the cell nearest the window's corner sits far outside the
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

      call report(times(corner) .eq. 1, &
           & "the window's corner escapes at once", nfail)

    end block corner_escapes

    inquire(file = paint_file, exist = on_disk)
    call report(on_disk, "the painting is on disk", nfail)

  end subroutine check_the_painting

end program test_orbit_painting
