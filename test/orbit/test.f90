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
!   7. the ladder painted, every level its own mesh derived from the
!      base by the mesh's own machinery: julia-1 and julia-2 on real
!      refinements (every triangle split into four, and again), and
!      julia+1 and julia+2 on agglomerated meshes (each part of the
!      quotient fused into one polygon); each file fresh on disk with
!      light and dark, cell counts thinning as the squint deepens
!=====================================================================!

program test_orbit_painting

  use iso_fortran_env      , only : dp => REAL64
  use class_mesh           , only : mesh
  use class_gmsh_loader    , only : gmsh_loader
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
  class(mesh)           , allocatable :: grid, base_grid
  class(paraview_writer), allocatable :: writer
  type(mesh)   :: g
  type(stored_graph) :: coarse

  complex(dp), allocatable :: zc(:)          ! cell centroids on the complex plane
  complex(dp), allocatable :: zpart(:)       ! part centroids on the complex plane
  real(dp)   , allocatable :: part_radius(:) ! how far a part's cells stray from its centroid
  integer    , allocatable :: succ(:)        ! the one arrow out of each cell (0 = off the window)
  integer    , allocatable :: times(:)       ! the painting, resolved in one pass
  real(dp)   , allocatable :: escape_time(:,:)
  type(string) :: field_labels(1)

  ! arrows of the coarse levels' paintings, host-carried for the rules
  integer, allocatable :: arrows_coarse1(:), arrows_coarse2(:)

  integer  :: ncells, nparts, n_home
  integer  :: ncells_base, ncells_fine1, ncells_fine2
  integer  :: ncells_coarse1, ncells_coarse2
  integer  :: checks_by_descent, full_scans, nfail
  real(dp) :: t0, t1
  real(dp) :: xlo, xhi, ylo, yhi   ! the window map, shared with the ladder

  nfail = 0

  field_labels(1) = string('escape_time')

  ! the base canvas comes off disk once; every other level of the
  ! ladder is derived from it by the mesh's own machinery
  load_the_base: block
    allocate(loader, source = gmsh_loader(mesh_file))
    allocate(base_grid, source = mesh(loader))
    deallocate(loader)
  end block load_the_base

  ! the base painting, with the deep check battery on its heels
  call paint_mesh_level(base_grid, paint_file, nfail)
  ncells_base = ncells
  write(*,'(1x,a,i0,a,i0)') "the squint claimed the base arrows with ", &
       & checks_by_descent, " distance checks where the full scan needs ", ncells*ncells
  write(*,'(1x,a,i0,a)') "the certificate sent ", full_scans, &
       & " landing(s) back to the full scan"
  call check_the_painting(nfail)

  ! the squint side: the base cells huddle into parts, the quotient
  ! graph carries the coarse solution, and the agglomerated mesh -
  ! every part fused into one polygon - carries it to paraview
  call paint_coarse_levels(nfail)

  ! the zoom side: the mesh refined by its own machinery, the whole
  ! pipeline rerun on each refined mesh
  zoom: block
    type(mesh) :: refined_once, refined_twice
    refined_once = base_grid % refined()
    call paint_mesh_level(refined_once, 'julia-1.vtu', nfail)
    ncells_fine1 = ncells
    refined_twice = refined_once % refined()
    call paint_mesh_level(refined_twice, 'julia-2.vtu', nfail)
    ncells_fine2 = ncells
  end block zoom

  call report(ncells_coarse2 .lt. ncells_coarse1 .and. &
       &      ncells_coarse1 .lt. ncells_base    .and. &
       &      ncells_base    .lt. ncells_fine1   .and. &
       &      ncells_fine1   .lt. ncells_fine2, &
       & "the ladder thins as the squint deepens", nfail)

  write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a)') "the ladder, finest to coarsest: ", &
       & ncells_fine2, " / ", ncells_fine1, " / ", ncells_base, " / ", &
       & ncells_coarse1, " / ", ncells_coarse2, " cells, each level its own mesh"

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all orbit painting checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " orbit painting check(s)"
     error stop
  end if

contains

  ! the affine map from mesh coordinates onto the window
  pure complex(dp) function to_window(x, y) result(zw)
    real(dp), intent(in) :: x, y
    zw = cmplx(window_half_width*(2.0_dp*(x - xlo)/(xhi - xlo) - 1.0_dp), &
         &     window_half_width*(2.0_dp*(y - ylo)/(yhi - ylo) - 1.0_dp), kind = dp)
  end function to_window

  !===================================================================!
  ! Paint one level of the ladder: take the mesh handed in - loaded
  ! or derived - and run the whole pipeline on it: the window map,
  ! the squint the descent reads, the arrows, the one-pass escape
  ! times, the painting written on that mesh's own cells. Light
  ! checks here; the deep battery runs once, on the base painting.
  !===================================================================!

  subroutine paint_mesh_level(level_grid, level_file, nfail)

    class(mesh)     , intent(in)    :: level_grid
    character(len=*), intent(in)    :: level_file
    integer         , intent(inout) :: nfail

    integer, allocatable :: members(:)

    complex(dp) :: z
    logical     :: on_disk
    integer     :: v, k, unit, ierr

    ! a stale painting would satisfy the on-disk check - burn it
    open(newunit = unit, file = level_file, status = 'old', iostat = ierr)
    if (ierr .eq. 0) close(unit, status = 'delete')

    ! this level's mesh becomes the pipeline's canvas
    if (allocated(grid))   deallocate(grid)
    if (allocated(writer)) deallocate(writer)
    allocate(grid, source = level_grid)

    g      = grid
    ncells = g % num_vertices

    ! stretch this mesh over the same window
    xlo = minval(grid % cell_centers(1, :)); xhi = maxval(grid % cell_centers(1, :))
    ylo = minval(grid % cell_centers(2, :)); yhi = maxval(grid % cell_centers(2, :))
    zc  = [(to_window(grid % cell_centers(1, v), grid % cell_centers(2, v)), v = 1, ncells)]

    ! the squint the descent reads
    nparts = max(1, ncells/cells_per_part)
    call g % partition_rcb(grid % cell_centers, nparts)
    coarse = stored_graph(g)

    if (allocated(zpart))       deallocate(zpart)
    if (allocated(part_radius)) deallocate(part_radius)
    allocate(zpart(nparts), part_radius(nparts))
    do k = 1, nparts
       members        = g % owned(k)
       zpart(k)       = sum(zc(members))/real(size(members), dp)
       part_radius(k) = sqrt(maxval(abs(zc(members) - zpart(k))**2))
    end do

    ! the arrows, claimed by descent
    if (allocated(succ)) deallocate(succ)
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

    ! paint in one pass and write on this mesh's own cells
    times = g % escape_times(successor, step_limit)
    if (allocated(escape_time)) deallocate(escape_time)
    allocate(escape_time(ncells, 1))
    escape_time(:, 1) = real(times, dp)
    n_home = count(times .eq. step_limit)
    call cpu_time(t1)

    allocate(writer, source = paraview_writer(grid))
    call writer % write(level_file, escape_time, field_labels)

    write(*,'(1x,a,i0,a,i0,a,i0,a,f6.3,a,a)') "painted ", ncells, " cells: ", &
         & ncells - n_home, " escaped, ", n_home, " home (", t1 - t0, " s) -> ", level_file

    inquire(file = level_file, exist = on_disk)
    call report(on_disk .and. n_home .gt. 0 .and. n_home .lt. ncells, &
         & level_file // " on its own mesh, light and dark", nfail)

  end subroutine paint_mesh_level

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


  !===================================================================!
  ! The squint side of the ladder. The base cells huddle into parts
  ! (about four cells each), the quotient graph gets its own arrows
  ! and escape times, and the agglomerated mesh - each part fused
  ! into one polygon - carries the coarse solution to paraview.
  ! Twice over for julia+2, squinting the quotient again.
  !===================================================================!

  subroutine paint_coarse_levels(nfail)

    integer, intent(inout) :: nfail

    type(stored_graph) :: quotient1, quotient2
    type(mesh)         :: agg1, agg2

    complex(dp), allocatable :: zq1(:), zq2(:)
    integer    , allocatable :: level_times(:), members(:)
    integer    , allocatable :: anc1(:), anc2(:)

    integer :: v, k, n1, n2

    ! squint once, by aggregation - a seed huddles with its
    ! neighbours, so every part is connected by construction and its
    ! boundary traces into one loop (coordinate bisection cuts by
    ! count and can scatter a part into disconnected islands)
    call g % partition_aggregate()
    quotient1 = stored_graph(g)
    n1        = quotient1 % num_vertices
    anc1      = [(g % part_of(v), v = 1, ncells)]

    allocate(zq1(n1))
    do k = 1, n1
       members = g % owned(k)
       zq1(k)  = sum(zc(members))/real(size(members), dp)
    end do

    arrows_coarse1 = quantized_arrows(zq1)
    level_times    = quotient1 % escape_times(coarse1_successor, step_limit)

    agg1 = grid % agglomerated(anc1, n1)
    call write_coarse_painting(agg1, real(level_times, dp), 'julia+1.vtu', nfail)
    ncells_coarse1 = n1

    ! squint twice, the same way
    call quotient1 % partition_aggregate()
    quotient2 = stored_graph(quotient1)
    n2        = quotient2 % num_vertices
    anc2      = [(quotient1 % part_of(anc1(v)), v = 1, ncells)]

    allocate(zq2(n2))
    do k = 1, n2
       members = quotient1 % owned(k)
       zq2(k)  = sum(zq1(members))/real(size(members), dp)
    end do

    arrows_coarse2 = quantized_arrows(zq2)
    level_times    = quotient2 % escape_times(coarse2_successor, step_limit)

    agg2 = grid % agglomerated(anc2, n2)
    call write_coarse_painting(agg2, real(level_times, dp), 'julia+2.vtu', nfail)
    ncells_coarse2 = n2

  end subroutine paint_coarse_levels

  ! burn any stale copy, write the coarse painting on its own
  ! agglomerated mesh, and demand it lands with light and dark
  subroutine write_coarse_painting(agg, field, filename, nfail)

    type(mesh)      , intent(in)    :: agg
    real(dp)        , intent(in)    :: field(:)
    character(len=*), intent(in)    :: filename
    integer         , intent(inout) :: nfail

    class(paraview_writer), allocatable :: coarse_writer

    real(dp), allocatable :: sheet(:,:)
    logical               :: on_disk
    integer               :: unit, ierr

    open(newunit = unit, file = filename, status = 'old', iostat = ierr)
    if (ierr .eq. 0) close(unit, status = 'delete')

    allocate(sheet(size(field), 1))
    sheet(:, 1) = field
    allocate(coarse_writer, source = paraview_writer(agg))
    call coarse_writer % write(filename, sheet, field_labels)

    inquire(file = filename, exist = on_disk)
    call report(on_disk .and. any(field .ge. real(step_limit, dp)) .and. &
         &      any(field .lt. real(step_limit, dp)), &
         & filename // " on its own agglomerated mesh, light and dark", nfail)

  end subroutine write_coarse_painting

  ! the coarse levels' successors, reading the host-carried arrows
  pure integer function coarse1_successor(v)
    integer, intent(in) :: v
    coarse1_successor = arrows_coarse1(v)
  end function coarse1_successor

  pure integer function coarse2_successor(v)
    integer, intent(in) :: v
    coarse2_successor = arrows_coarse2(v)
  end function coarse2_successor

  ! the arrows of a small centroid set, by full scan: one step of the
  ! map from each centroid, the nearest centroid claims the landing
  pure function quantized_arrows(zs) result(arrows)

    complex(dp), intent(in) :: zs(:)

    integer, allocatable :: arrows(:)

    complex(dp) :: z
    real(dp)    :: best, d2
    integer     :: v, w, nearest

    allocate(arrows(size(zs)))
    do v = 1, size(zs)
       z = zs(v)*zs(v) + julia_constant
       if (abs(z) .gt. escape_radius) then
          arrows(v) = 0
          cycle
       end if
       nearest = 1
       best    = huge(1.0_dp)
       do w = 1, size(zs)
          d2 = abs(z - zs(w))**2
          if (d2 .lt. best) then
             best    = d2
             nearest = w
          end if
       end do
       arrows(v) = nearest
    end do

  end function quantized_arrows

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
