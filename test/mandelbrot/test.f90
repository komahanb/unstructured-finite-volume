!=====================================================================!
! The escape-time fractal, computed as PHYSICS: the julia set of the
! douady rabbit painted by the ordinary solver stack -
!
!    (mandelbrot law) ──▶ (assembler) ──▶ (forward euler) ──▶ (writer)
!     z² + c - z as a       no flux,        h = 1: each step     one
!     reacting source       source only     IS the map           picture
!
! and then refined by its own adjoint, cycle after cycle:
!
!    march ──▶ adjoint ──▶ flag ──▶ refine ──▶ march ──▶ ...
!
! The run is read from paint.cfg (keyword per line, the class_config
! pattern); defaults reproduce the douady rabbit on square-tri-40.
!
! Checks:
!   1. the law's coupling graph is real: two variables, ONE edge
!      (the complex square mixes them) - the first law with an edge
!   2. the law's jacobian dS/dq matches finite differences - this is
!      the matrix the adjoint will ride
!   3. the marched escape counts agree with a plain z -> z² + c loop
!      cell by cell (identical almost everywhere; the marched
!      arithmetic may round differently at an escape tie, so a
!      whisker of ±1 on a fraction of a percent of cells is honest)
!   4. the adjoint rides the marcher's own chain: seed dJ/dz_T at the
!      last step, pull back through (I - dS/dq)^T edge by edge, and
!      what lands on the first vertex is dJ/dz0 for every cell at
!      once - certified against central differences of a re-marched J
!   5. the adjoint chooses where the mesh deserves more cells: flag
!      the most sensitive fraction (the fractal's boundary lights
!      up), refine once, march only the flagged cells' children,
!      inherit the rest - and the flagged cells must hold at least
!      twice their blind share of the energy error
!   6. cycle the whole loop and the sensitivity recedes: each refined
!      canvas answers with a smaller bulk |dJ/dz0| and a smaller
!      dual-weighted error estimate - while the peak may sharpen,
!      because finer cells land nearer the chaotic boundary
!
! julia-physics.vtu holds the rabbit, julia-sensitivity-<k>.vtu the
! adjoint's view of cycle k's canvas, julia-refined.vtu the sharpened
! painting of cycle 1's selective march.
!
! A nonzero exit (error stop) means a check failed.
!=====================================================================!

program test_mandelbrot

  use iso_fortran_env      , only : dp => REAL64
  use class_gmsh_loader    , only : gmsh_loader
  use class_mesh           , only : mesh
  use class_assembler      , only : assembler
  use class_diffusion_flux , only : diffusion_flux
  use class_mandelbrot     , only : mandelbrot_source
  use class_forward_euler  , only : forward_euler
  use class_paraview_writer, only : paraview_writer
  use class_string         , only : string
  use class_file           , only : file
  use class_config         , only : split_words
  use interface_physics    , only : point_state

  implicit none

  ! the run - paint.cfg overrides, these defaults ARE the douady
  ! rabbit drawn with the same eyes as the orbit painting
  character(len=256) :: mesh_file         = '../square-tri-40.msh'
  real(dp)           :: c(2)              = [-0.123_dp, 0.745_dp]
  real(dp)           :: window_half_width = 1.6_dp
  integer            :: patience          = 40
  real(dp)           :: flag_fraction     = 0.25_dp
  integer            :: num_refinement_cycles = 2

  real(dp), parameter :: escape_radius_sq = 4.0_dp

  class(gmsh_loader), allocatable :: loader
  class(mesh)       , allocatable :: grid
  class(assembler)  , allocatable :: fvm
  type(forward_euler)             :: ti
  type(mandelbrot_source)         :: law_g

  ! shared by the checks, refreshed each cycle: the window map, every
  ! cell's start, the painting, the adjoint and its sensitivity
  real(dp), allocatable :: z0(:,:), eta(:), lam(:,:)
  integer , allocatable :: marched(:)
  real(dp), allocatable :: peak_eta(:), bulk_eta(:), est_error(:)
  real(dp)              :: xlo, xhi, ylo, yhi

  integer :: cycle, nfail

  nfail = 0

  call read_config('paint.cfg')
  allocate(peak_eta(num_refinement_cycles), bulk_eta(num_refinement_cycles), &
       &   est_error(num_refinement_cycles))

  call check_coupling_graph(nfail)
  call check_jacobian(nfail)

  allocate(loader, source = gmsh_loader(trim(mesh_file)))
  allocate(grid  , source = mesh(loader))

  ! the loop the adjoint drives: march the canvas, ask it where it is
  ! tender, refine, march again
  cycles: do cycle = 1, num_refinement_cycles

     call march_and_paint(cycle, nfail)
     call compute_sensitivity(cycle)
     peak_eta(cycle) = maxval(eta)

     if (cycle .eq. 1) then
        call check_adjoint_gradient(nfail)
        call check_adjoint_refinement(nfail)
     end if

     if (cycle .lt. num_refinement_cycles) call refine_canvas()

  end do cycles

  call check_receding_sensitivity(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all mandelbrot-physics checks passed"
  else
     write(*,'(1x,i0,a)') nfail, " mandelbrot-physics checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! The run reader: keyword per line, the class_config pattern.
  ! A missing file means the defaults paint.
  !===================================================================!

  subroutine read_config(filename)

    character(len=*), intent(in) :: filename

    type(file)                :: cfile
    type(string), allocatable :: lines(:), tok(:)
    integer                   :: iline, ntok
    logical                   :: found

    inquire(file = filename, exist = found)
    if (.not. found) return

    cfile = file(filename, 256)
    call cfile % read_lines(lines)

    do iline = 1, size(lines)

       call split_words(lines(iline) % str, ntok, tok)
       if (ntok .lt. 1) cycle
       if (tok(1) % str(1:1) .eq. "#") cycle

       select case (trim(tok(1) % str))
       case ("mesh")
          mesh_file = trim(tok(2) % str)
       case ("constant")                       ! constant  c_re c_im
          c = [tok(2) % asreal(), tok(3) % asreal()]
       case ("window")
          window_half_width = tok(2) % asreal()
       case ("patience")
          patience = tok(2) % asinteger()
       case ("flag_fraction")
          flag_fraction = tok(2) % asreal()
       case ("refinement_cycles")
          num_refinement_cycles = tok(2) % asinteger()
       case default
          write(*,'(1x,3a)') "paint.cfg: unknown keyword '", trim(tok(1) % str), "'"
          error stop
       end select

    end do

    write(*,'(1x,a,i0,a,f5.2)') "paint.cfg read: ", num_refinement_cycles, &
         & " cycle(s), flag fraction ", flag_fraction

  end subroutine read_config

  !===================================================================!
  ! 1: the law IS a graph, and its one edge is the complex square
  !===================================================================!

  subroutine check_coupling_graph(nfail)

    integer, intent(inout) :: nfail
    type(mandelbrot_source) :: law

    law = mandelbrot_source(c)

    if (law % num_vertices .eq. 2 .and. law % num_edges .eq. 1 .and. &
         & all(law % neighbours(1) .eq. [2]) .and. law % degree(2) .eq. 1) then
       write(*,'(1x,a)') "PASS : the law couples its two variables - one edge, (u)---(v)"
    else
       write(*,'(1x,a)') "FAIL : the coupling graph is wrong"
       nfail = nfail + 1
    end if

  end subroutine check_coupling_graph

  !===================================================================!
  ! 2: dS/dq (the edge as numbers, the adjoint's engine) vs central
  ! finite differences at a nowhere-special state
  !===================================================================!

  subroutine check_jacobian(nfail)

    integer, intent(inout) :: nfail

    type(mandelbrot_source) :: law
    type(point_state)       :: st
    real(dp) :: dS(2,2), fd(2,2), Sp(2), Sm(2), err
    real(dp), parameter :: eps = 1.0e-7_dp
    integer  :: j

    law = mandelbrot_source(c)

    allocate(st % q(2), st % gradq(3,2))
    st % nv = 2
    st % q  = [0.3_dp, -0.2_dp]
    st % gradq = 0.0_dp

    dS = real(law % dsource_dq(st), dp)

    do j = 1, 2
       st % q(j) = st % q(j) + eps
       Sp = real(law % value(st), dp)
       st % q(j) = st % q(j) - 2.0_dp*eps
       Sm = real(law % value(st), dp)
       st % q(j) = st % q(j) + eps
       fd(:,j) = (Sp - Sm)/(2.0_dp*eps)
    end do

    err = maxval(abs(dS - fd))
    if (err .lt. 1.0e-6_dp) then
       write(*,'(1x,a,es10.2)') "PASS : dS/dq matches finite differences, err ", err
    else
       write(*,'(1x,a,es10.2)') "FAIL : dS/dq vs finite differences, err ", err
       nfail = nfail + 1
    end if

  end subroutine check_jacobian

  !===================================================================!
  ! March the current canvas: every cell starts at its own point of
  ! the window (the map is fixed on cycle 1 so every level paints the
  ! SAME window), forty steps of size one, escape counts read off the
  ! trajectory. On cycle 1 the counts are certified against a plain
  ! z -> z² + c loop with no framework at all, and the rabbit lands
  ! on disk.
  !===================================================================!

  subroutine march_and_paint(cycle, nfail)

    integer, intent(in)    :: cycle
    integer, intent(inout) :: nfail

    type(paraview_writer), allocatable :: writer
    type(string)          :: labels(1)
    real(dp), allocatable :: paint(:,:)
    integer , allocatable :: plain(:)
    real(dp) :: u, v, zr, zi
    integer  :: ncells, icell, k, mismatches, nsteps

    if (allocated(fvm)) deallocate(fvm)
    allocate(fvm, source = assembler(grid))

    ! no flux, only the reacting law; walls that ask nothing
    call fvm % set_equation(diffusion_flux(0.0_dp, 2), mandelbrot_source(c))
    call fvm % set_neumann("BoundaryLeft"  , 0.0_dp)
    call fvm % set_neumann("BoundaryRight" , 0.0_dp)
    call fvm % set_neumann("BoundaryTop"   , 0.0_dp)
    call fvm % set_neumann("BoundaryBottom", 0.0_dp)

    ncells = fvm % grid % num_cells

    ! the window map, fixed once on the first canvas
    if (cycle .eq. 1) then
       xlo = minval(fvm % grid % cell_centers(1, :)); xhi = maxval(fvm % grid % cell_centers(1, :))
       ylo = minval(fvm % grid % cell_centers(2, :)); yhi = maxval(fvm % grid % cell_centers(2, :))
    end if

    if (allocated(z0)) deallocate(z0)
    allocate(z0(2, ncells))
    do icell = 1, ncells
       z0(1, icell) = window_half_width*(2.0_dp*(fvm % grid % cell_centers(1, icell) - xlo)/(xhi - xlo) - 1.0_dp)
       z0(2, icell) = window_half_width*(2.0_dp*(fvm % grid % cell_centers(2, icell) - ylo)/(yhi - ylo) - 1.0_dp)
       fvm % phi(fvm % grid % dof(icell, 1)) = z0(1, icell)
       fvm % phi(fvm % grid % dof(icell, 2)) = z0(2, icell)
    end do

    ! steps of size one: the map, marched
    ti = forward_euler(fvm, 0.0_dp, real(patience, dp), 1.0_dp)
    call ti % integrate()

    ! escape counts from the trajectory (step k holds z after k-1
    ! applications of the map; patient cells keep the full patience)
    if (allocated(marched)) deallocate(marched)
    allocate(marched(ncells))
    do icell = 1, ncells
       marched(icell) = patience
       do k = 2, ti % num_vertices
          u = real(ti % U(k, fvm % grid % dof(icell, 1), 1), dp)
          v = real(ti % U(k, fvm % grid % dof(icell, 2), 1), dp)
          if (u*u + v*v .gt. escape_radius_sq) then
             marched(icell) = k - 1
             exit
          end if
       end do
    end do

    first_canvas_only: if (cycle .eq. 1) then

       ! the plain map, no framework at all
       allocate(plain(ncells))
       do icell = 1, ncells
          call march_cell(z0(1, icell), z0(2, icell), zr, zi, nsteps)
          plain(icell) = nsteps
       end do

       mismatches = count(marched .ne. plain)
       write(*,'(1x,a,i0,a,i0,a)') "escape counts: ", ncells - mismatches, " of ", ncells, &
            & " cells identical to the plain map"
       if (real(mismatches, dp)/real(ncells, dp) .le. 0.005_dp .and. &
            & maxval(abs(marched - plain)) .le. 1) then
          write(*,'(1x,a)') "PASS : the marched map IS the map (ties may round, nothing more)"
       else
          write(*,'(1x,a,i0,a,i0)') "FAIL : mismatches ", mismatches, "  max drift ", &
               & maxval(abs(marched - plain))
          nfail = nfail + 1
       end if

       ! the picture: the rabbit, by simulation
       allocate(paint(ncells, 1))
       paint(:,1) = real(marched, dp)
       labels(1)  = string("escape")
       allocate(writer, source = paraview_writer(fvm % grid))
       call writer % write('julia-physics.vtu', paint, labels)
       write(*,'(1x,a)') "painted julia-physics.vtu (the douady rabbit, by physics)"

    end if first_canvas_only

  end subroutine march_and_paint

  !===================================================================!
  ! The adjoint of this cycle's painting. The functional is the final
  ! energy  J = 1/2 sum_c vol_c |z_c(T)|²; its derivative seeds the
  ! last vertex of the step chain, accumulate_adjoint walks the chain
  ! back, and every edge pulls the adjoint through (I - dS/dq)^T -
  ! the law's matrix from check 2 doing the adjoint's work. What
  ! lands on the first vertex is dJ/dz0 for every cell at once; its
  ! magnitude is the sensitivity the flags read, written to disk as
  ! this cycle's julia-sensitivity-<k>.vtu.
  !===================================================================!

  subroutine compute_sensitivity(cycle)

    integer, intent(in) :: cycle

    type(paraview_writer), allocatable :: writer
    type(string)          :: labels(1)
    real(dp), allocatable :: seed(:), mass(:), paint(:,:)
    character(len=64)     :: fname
    integer               :: icell

    law_g = mandelbrot_source(c)

    associate(n => fvm % num_state_vars, nt => ti % num_vertices, &
         &    ncells => fvm % grid % num_cells)

      allocate(mass(n), seed(n))
      if (allocated(lam)) deallocate(lam)
      allocate(lam(n, nt))

      call fvm % get_lumped_mass(mass)
      seed = mass * real(ti % U(nt, :, 1), dp)

      call ti % accumulate_adjoint(seed, step_pullback, nt, lam)

      ! how tender the painting is to where each cell starts
      if (allocated(eta)) deallocate(eta)
      allocate(eta(ncells))
      do icell = 1, ncells
         eta(icell) = hypot(lam(fvm % grid % dof(icell,1), 1), &
              &             lam(fvm % grid % dof(icell,2), 1))
      end do

      write(*,'(1x,a,i0,a,i0,a,es9.2,a,es9.2)') "cycle ", cycle, ": ", &
           & ncells, " cells, |dJ/dz0| spans ", minval(eta), " to ", maxval(eta)

      ! the receding measures: the bulk sensitivity (geometric mean -
      ! the log-average a 39-decade span calls for) and the dual-
      ! weighted error estimate, each cell's sensitivity times its
      ! reach (sum eta * sqrt(vol) - how far a cell's start can stray)
      bulk_eta(cycle)  = exp(sum(log(max(eta, 1.0e-30_dp)))/real(ncells, dp))
      est_error(cycle) = sum(eta * sqrt(fvm % grid % cell_volumes))

      write(fname,'(a,i0,a)') 'julia-sensitivity-', cycle, '.vtu'
      allocate(paint(ncells, 1))
      paint(:,1) = log10(max(eta, 1.0e-30_dp))
      labels(1)  = string("log10_sensitivity")
      allocate(writer, source = paraview_writer(fvm % grid))
      call writer % write(trim(fname), paint, labels)

    end associate

  end subroutine compute_sensitivity

  !===================================================================!
  ! 4: the adjoint against central differences of a re-marched J, on
  ! two probes chosen where a difference can breathe: a chaotic cell
  ! (huge eta) leaves the linear regime within machine epsilon, a
  ! dead cell (tiny eta) drowns in roundoff - so probe near eta ~ 1
  ! and eta ~ 0.01, both clear of the freeze threshold, the one kink
  ! a difference cannot straddle.
  !===================================================================!

  subroutine check_adjoint_gradient(nfail)

    integer, intent(inout) :: nfail

    real(dp), parameter :: eps = 1.0e-6_dp
    real(dp) :: jp, jm, fd(2), ad(2), relerr, worst
    real(dp) :: zr, zi
    integer  :: icell, j, p1, p2, probe, picked(2), nsteps

    picked = [pick_probe(1.0_dp), pick_probe(1.0e-2_dp)]

    worst = 0.0_dp
    do probe = 1, 2
       icell = picked(probe)
       p1 = fvm % grid % dof(icell, 1)
       p2 = fvm % grid % dof(icell, 2)
       ad = [lam(p1,1), lam(p2,1)]
       do j = 1, 2
          call march_cell(z0(1,icell) + merge(eps,0.0_dp,j.eq.1), &
               &          z0(2,icell) + merge(eps,0.0_dp,j.eq.2), zr, zi, nsteps)
          jp = 0.5_dp*fvm % grid % cell_volumes(icell)*(zr*zr + zi*zi)
          call march_cell(z0(1,icell) - merge(eps,0.0_dp,j.eq.1), &
               &          z0(2,icell) - merge(eps,0.0_dp,j.eq.2), zr, zi, nsteps)
          jm = 0.5_dp*fvm % grid % cell_volumes(icell)*(zr*zr + zi*zi)
          fd(j) = (jp - jm)/(2.0_dp*eps)
       end do
       relerr = maxval(abs(fd - ad))/maxval(abs(ad))
       worst  = max(worst, relerr)
    end do

    if (worst .lt. 1.0e-5_dp) then
       write(*,'(1x,a,es10.2)') &
            & "PASS : the chain's reverse walk IS dJ/dz0 (vs differences), err ", worst
    else
       write(*,'(1x,a,es10.2)') "FAIL : adjoint vs finite differences, err ", worst
       nfail = nfail + 1
    end if

  end subroutine check_adjoint_gradient

  !===================================================================!
  ! One edge of the step chain, pulled back: per cell, apply the
  ! transpose of the step's jacobian I - dS/dq - read at the TAIL's
  ! state, the state the step left from - to the head's adjoint.
  !===================================================================!

  subroutine step_pullback(tail, head, adjoint_head, contribution)

    integer , intent(in)  :: tail, head
    real(dp), intent(in)  :: adjoint_head(:)
    real(dp), intent(out) :: contribution(:)

    type(point_state) :: st
    real(dp) :: dS(2,2), a1, a2
    integer  :: icell, p1, p2

    allocate(st % q(2), st % gradq(3,2))
    st % nv    = 2
    st % gradq = 0.0_dp

    do icell = 1, fvm % grid % num_cells
       p1 = fvm % grid % dof(icell, 1)
       p2 = fvm % grid % dof(icell, 2)
       st % q = [real(ti % U(tail, p1, 1), dp), real(ti % U(tail, p2, 1), dp)]
       st % x = fvm % grid % cell_centers(:, icell)
       dS = real(law_g % dsource_dq(st), dp)
       a1 = adjoint_head(p1)
       a2 = adjoint_head(p2)
       contribution(p1) = (1.0_dp - dS(1,1))*a1 - dS(2,1)*a2
       contribution(p2) = -dS(1,2)*a1 + (1.0_dp - dS(2,2))*a2
    end do

  end subroutine step_pullback

  !===================================================================!
  ! The closest any step of a cell's marched orbit comes to the
  ! freeze threshold |z|² = 4 - the kink the probes must stay off.
  !===================================================================!

  pure real(dp) function margin_of(icell) result(margin)

    integer, intent(in) :: icell

    real(dp) :: rr
    integer  :: k

    margin = huge(1.0_dp)
    do k = 1, ti % num_vertices
       rr = real(ti % U(k, fvm % grid % dof(icell,1), 1), dp)**2 &
          + real(ti % U(k, fvm % grid % dof(icell,2), 1), dp)**2
       margin = min(margin, abs(rr - 4.0_dp))
    end do

  end function margin_of

  !===================================================================!
  ! The map with the law's freeze, run to the full patience - the
  ! same dynamics the framework marches. Hands back the final state
  ! (for energies and difference probes) and the escape count (for
  ! the painting), one walk for both.
  !===================================================================!

  pure subroutine march_cell(zr0, zi0, zr, zi, n)

    real(dp), intent(in)  :: zr0, zi0
    real(dp), intent(out) :: zr, zi
    integer , intent(out) :: n

    real(dp) :: t
    integer  :: k

    zr = zr0; zi = zi0
    n  = patience
    do k = 1, patience
       if (zr*zr + zi*zi .le. escape_radius_sq) then
          t  = zr*zr - zi*zi + c(1)
          zi = 2.0_dp*zr*zi  + c(2)
          zr = t
       end if
       ! outside after this step - mapped there, or born there and
       ! held - is escaped at this length, exactly as the framework's
       ! frozen trajectory reads
       if (n .eq. patience .and. zr*zr + zi*zi .gt. escape_radius_sq) n = k
    end do

  end subroutine march_cell

  !===================================================================!
  ! The margin-safe cell whose sensitivity sits closest to the given
  ! target - how the difference probes are chosen.
  !===================================================================!

  pure integer function pick_probe(target) result(best)

    real(dp), intent(in) :: target

    real(dp) :: gap, d
    integer  :: icell

    best = 0
    gap  = huge(1.0_dp)
    do icell = 1, fvm % grid % num_cells
       if (margin_of(icell) .lt. 0.05_dp) cycle
       d = abs(log(max(eta(icell), 1.0e-30_dp)) - log(target))
       if (d .lt. gap) then
          gap  = d
          best = icell
       end if
    end do

  end function pick_probe

  !===================================================================!
  ! The flag threshold: bisect tau until the configured fraction of
  ! the cells carry more sensitivity than it.
  !===================================================================!

  pure real(dp) function flag_threshold() result(tau)

    real(dp) :: lo, hi
    integer  :: k, want

    want = int(flag_fraction*real(fvm % grid % num_cells, dp))
    lo = minval(eta); hi = maxval(eta)
    do k = 1, 60
       tau = 0.5_dp*(lo + hi)
       if (count(eta .ge. tau) .gt. want) then
          lo = tau
       else
          hi = tau
       end if
    end do

  end function flag_threshold

  !===================================================================!
  ! 5: the adjoint chooses where the mesh deserves more cells. Flag
  ! the most sensitive fraction of the cells - the fractal's boundary,
  ! by the adjoint's numbers - refine the mesh once, march only the
  ! flagged cells' children, and let every other child inherit its
  ! parent's answer down the refinement tree (children of cell v sit
  ! at 4(v-1)+1..4v). The fully marched fine level is the oracle, and
  ! the verdict is measured in the functional the adjoint serves: the
  ! flagged cells must hold at least TWICE their blind share of the
  ! energy error inheritance would commit (a random fraction holds
  ! its own share; the adjoint's here holds well over twice - what a
  ! derivative cannot see is the jump across an escape ring, so total
  ! capture is not on offer from any first-order indicator). The
  ! paint mismatch count is reported alongside as the picture's own
  ! tally.
  !===================================================================!

  subroutine check_adjoint_refinement(nfail)

    integer, intent(inout) :: nfail

    type(mesh)                         :: fine
    type(paraview_writer), allocatable :: writer
    type(string)          :: labels(1)
    logical , allocatable :: flagged(:)
    integer , allocatable :: adaptive(:)
    real(dp), allocatable :: paint(:,:)
    real(dp) :: tau, x0f, y0f, recovered
    real(dp) :: zr, zi, e2_child, e2_parent, err, e_inherit, e_adaptive
    integer  :: v, child, nfine, nfull, miss_inherit, miss_adaptive

    tau = flag_threshold()
    allocate(flagged(fvm % grid % num_cells))
    flagged = eta .ge. tau

    ! one refinement of the whole canvas; the adjoint decides which
    ! children are worth marching
    fine  = fvm % grid % refined()
    nfine = fine % num_cells

    allocate(adaptive(nfine))

    e_inherit     = 0.0_dp
    e_adaptive    = 0.0_dp
    miss_inherit  = 0
    miss_adaptive = 0
    do v = 1, fvm % grid % num_cells

       ! the parent's final energy, read off the marched trajectory
       e2_parent = real(ti % U(ti % num_vertices, fvm % grid % dof(v,1), 1), dp)**2 &
            &    + real(ti % U(ti % num_vertices, fvm % grid % dof(v,2), 1), dp)**2

       do child = 4*(v-1) + 1, 4*v

          ! the child's own start, through the SAME window map, and
          ! the oracle's answer for it
          x0f = window_half_width*(2.0_dp*(fine % cell_centers(1, child) - xlo)/(xhi - xlo) - 1.0_dp)
          y0f = window_half_width*(2.0_dp*(fine % cell_centers(2, child) - ylo)/(yhi - ylo) - 1.0_dp)
          call march_cell(x0f, y0f, zr, zi, nfull)
          e2_child = zr*zr + zi*zi

          ! what inheriting the parent's answer would cost, in the
          ! functional's own currency
          err       = 0.5_dp*fine % cell_volumes(child)*abs(e2_child - e2_parent)
          e_inherit = e_inherit + err
          if (.not. flagged(v)) e_adaptive = e_adaptive + err

          ! the adaptive painting marches only what the adjoint chose
          if (flagged(v)) then
             adaptive(child) = nfull
          else
             adaptive(child) = marched(v)
          end if

          if (marched(v)      .ne. nfull) miss_inherit  = miss_inherit  + 1
          if (adaptive(child) .ne. nfull) miss_adaptive = miss_adaptive + 1

       end do
    end do

    recovered = 1.0_dp
    if (e_inherit .gt. tiny(1.0_dp)) recovered = 1.0_dp - e_adaptive/e_inherit

    write(*,'(1x,a,i0,a,f5.1,a)') "the adjoint flagged ", count(flagged), &
         & " of the cells; their children hold ", 100.0_dp*recovered, &
         & "% of the energy error inheritance would commit"
    write(*,'(1x,a,i0,a,i0,a,i0,a)') "the picture's tally: ", miss_inherit, &
         & " of ", nfine, " fine cells wrong by inheritance, ", miss_adaptive, &
         & " after the adjoint's marching"

    associate(share => real(count(flagged), dp)/real(fvm % grid % num_cells, dp))
      if (recovered .ge. 2.0_dp*share .and. share .le. 1.35_dp*flag_fraction) then
         write(*,'(1x,a,f4.1,a)') &
              & "PASS : the adjoint's flags target ", recovered/share, &
              & "x their blind share of the error"
      else
         write(*,'(1x,a,f6.2,a,f6.2,a)') "FAIL : flagged cells hold ", &
              & 100.0_dp*recovered, "% of the error against a blind ", &
              & 100.0_dp*share, "%"
         nfail = nfail + 1
      end if
    end associate

    ! the sharpened painting the flags bought
    allocate(paint(nfine, 1))
    paint(:,1) = real(adaptive, dp)
    labels(1)  = string("escape")
    allocate(writer, source = paraview_writer(fine))
    call writer % write('julia-refined.vtu', paint, labels)
    write(*,'(1x,a)') "painted julia-refined.vtu"

  end subroutine check_adjoint_refinement

  !===================================================================!
  ! Trade the canvas for its refinement: every triangle into four,
  ! by the mesh's own machinery. The next cycle marches this one.
  !===================================================================!

  subroutine refine_canvas()

    type(mesh) :: fine

    fine = grid % refined()
    deallocate(grid)
    allocate(grid, source = fine)

  end subroutine refine_canvas

  !===================================================================!
  ! 6: run the loop and the sensitivity recedes - in the measures
  ! that ought to recede. The bulk sensitivity falls (a finer cell
  ! claims less of the objective) and so does the dual-weighted error
  ! estimate sum eta*sqrt(vol) (less error left to buy). The PEAK is
  ! deliberately not asserted: a finer canvas lands centroids nearer
  ! the chaotic boundary, and chaos amplifies faster than volume
  ! shrinks - the sensitivity ridge SHARPENS as it localizes. That is
  ! the fractal being a fractal, and the run reports it honestly.
  !===================================================================!

  subroutine check_receding_sensitivity(nfail)

    integer, intent(inout) :: nfail

    integer :: k
    logical :: recedes

    if (num_refinement_cycles .lt. 2) return

    recedes = .true.
    do k = 2, num_refinement_cycles
       write(*,'(1x,a,i0,a,i0,a,f7.2,a,f7.2,a,f7.2)') "cycle ", k-1, " -> ", k, &
            & ": bulk sensitivity fell x", bulk_eta(k-1)/bulk_eta(k), &
            & ", error estimate fell x", est_error(k-1)/est_error(k), &
            & ", peak moved x", peak_eta(k)/peak_eta(k-1)
       if (bulk_eta(k) .ge. bulk_eta(k-1))   recedes = .false.
       if (est_error(k) .ge. est_error(k-1)) recedes = .false.
    end do

    if (recedes) then
       write(*,'(1x,a)') "PASS : the sensitivity recedes as the canvas refines"
    else
       write(*,'(1x,a)') "FAIL : a refined canvas answered with MORE estimated error"
       nfail = nfail + 1
    end if

  end subroutine check_receding_sensitivity

end program test_mandelbrot
