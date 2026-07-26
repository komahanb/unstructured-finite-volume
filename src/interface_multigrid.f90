!=====================================================================!
! This module is the abstract multigrid method: the squint made
! mechanical. Build a coarse hierarchy from a fine operator (setup),
! then use it either as a
! preconditioner (apply: one V-cycle, z = M^-1 r - the deferred apply
! inherited from preconditioner) or as a standalone solver (solve:
! cycle until converged).
!
! Everything mechanical lives here, written once: the level storage,
! the setup driver (coarsen - tentative prolongation over the parts -
! smoothed prolongation - galerkin coarse operator, recursing until the
! coarsest level is small, then a dense LU there), the cycle engine,
! and the coarse direct solve.
!
! The cycle itself is not code - it is data, and the data is a graph.
! A schedule is a directed graph of stations and moves: each vertex
! numbered by the level it visits, each arrow the next leg of the
! trip. The graph reads its own route (source_path); apply_cycle only
! interprets the numbers, which count squints from the given mesh:
! level 0 is the mesh itself, level +k is k squints coarser (the
! deepest solved exactly), level -k is k zooms finer - refine(n)
! builds those by splitting the mesh operator's own graph into
! children, the operator lifted by interpolation alone. A repeated
! station smooths again. setup builds the classical V as one such
! schedule and apply follows it; nothing about the V is wired in.
!
! A concrete kind answers exactly one deferred question: how does this
! level squint - which fine vertices huddle into which coarse part.
! Algebraic multigrid answers by strength (the strength graph's own
! aggregation); geometric multigrid answers by coordinates (recursive
! bisection of the mesh graph). Everything after the squint is shared.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_multigrid

  use iso_fortran_env        , only : dp => REAL64
  use class_csr              , only : csr_matrix
  use interface_graph        , only : digraph, power_iteration
  use class_stored_graph     , only : stored_graph, stored_digraph
  use interface_linear_solver, only : preconditioner

  implicit none

  private
  public :: multigrid

  !-------------------------------------------------------------------!
  ! This type holds one multigrid level.
  !-------------------------------------------------------------------!

  type :: multigrid_level
     type(csr_matrix)      :: A           ! This is the operator on this level.
     type(csr_matrix)      :: P           ! This is the prolongation (fine x coarse);
                                          ! the restriction is NOT stored -
                                          ! it is P's edges walked against
                                          ! their arrows (matvec_transpose).
     real(dp), allocatable :: Dinv(:)     ! This is the smoother diagonal 1/diag(A).
     integer , allocatable :: parts(:)    ! The squint's answer is kept: which coarse vertex each vertex joins.
  end type multigrid_level

  !-------------------------------------------------------------------!
  ! This type carries the travelling state of a cycle: each level's
  ! right-hand side and correction while the trip is underway.
  !-------------------------------------------------------------------!

  type :: level_vectors
     real(dp), allocatable :: b(:), x(:)
  end type level_vectors

  !-------------------------------------------------------------------!
  ! A multigrid is a preconditioner (apply = one cycle) that can also
  ! solve standalone. The hierarchy and every mechanical part live on
  ! this abstract type; a concrete kind only answers how to squint.
  !-------------------------------------------------------------------!

  type, abstract, extends(preconditioner) :: multigrid

     type(multigrid_level), allocatable :: levels(:)   ! Level 1 is the finest; nlevel is the coarsest.
     integer                            :: nlevel = 0

     ! The coarsest level holds a dense LU factor.
     real(dp), allocatable :: Ac_dense(:,:)
     integer , allocatable :: ipiv(:)
     integer               :: ncoarse = 0

     ! The refined levels below the given mesh are the zoom side of
     ! the ladder, built on request by refine; finer_levels(k) is
     ! level -k.
     type(multigrid_level), allocatable :: finer_levels(:)
     integer                            :: n_finer = 0

     ! These are the parameters of the mechanism.
     real(dp) :: jacobi_w            = 0.6667_dp   ! This is the smoother weight (~2/3).
     integer  :: max_levels          = 25
     integer  :: coarse_size         = 50
     integer  :: npre                = 1
     integer  :: npost               = 1
     integer  :: children_per_vertex = 4           ! This counts the children a refined vertex splits into.

     ! These levels keep their ancestry for export, configured before
     ! setup; storage is on demand, and nothing is kept unasked.
     integer, allocatable :: exported_levels(:)

     ! The classical V is built by setup as an ordinary schedule.
     type(stored_digraph) :: v_cycle_schedule

   contains

     ! The one deferred question asks which fine vertex huddles into
     ! which coarse part on level lev.
     procedure(coarsen_interface), deferred :: coarsen

     ! The mechanism is shared.
     procedure :: setup       ! Build the coarse hierarchy from A.
     procedure :: refine      ! Build n refined levels below the given mesh.
     procedure :: apply       ! Apply z = M^-1 r, one trip along the V.
     procedure :: apply_cycle ! Produce z from one trip along any injected schedule.
     procedure :: solve       ! Solve A x = b, cycling to the tolerance.
     procedure :: num_levels    ! Count the levels in the coarse hierarchy.
     procedure :: export_levels ! Configure which levels keep their ancestry.
     procedure :: ancestors     ! Report the ancestry map of a configured signed level.

     procedure, private :: keeps
     procedure, private :: smooth
     procedure, private :: spectral_radius

  end type multigrid

  !-------------------------------------------------------------------!
  ! This is the deferred interface of the squint.
  !-------------------------------------------------------------------!

  abstract interface

     !================================================================!
     ! Stamp each fine vertex of level lev with its coarse part
     ! (agg(v) in 1..naggr). The level's operator is available as
     ! this % levels(lev) % A.
     !================================================================!

     subroutine coarsen_interface(this, lev, agg, naggr)
       import :: multigrid
       class(multigrid)    , intent(inout) :: this
       integer             , intent(in)    :: lev
       integer, allocatable, intent(out)   :: agg(:)
       integer             , intent(out)   :: naggr
     end subroutine coarsen_interface

  end interface

contains

  !===================================================================!
  ! Build the hierarchy from the fine operator A: at each level ask
  ! the concrete kind for its squint (coarsen), then build the
  ! tentative prolongation over the parts, smooth it, and form the
  ! galerkin coarse operator Ac = R A P - recursing until the coarsest
  ! level is small, where a dense LU waits.
  !===================================================================!

  impure subroutine setup(this, A)

    class(multigrid), intent(inout) :: this
    type(csr_matrix), intent(in)    :: A

    integer               :: lev, n, naggr, info, i
    integer , allocatable :: agg(:)
    real(dp), allocatable :: dinv(:)
    real(dp)              :: omega, rho
    type(csr_matrix)      :: P0, AP0, Pt

    !-----------------------------------------------------------------!
    ! Setup is re-entrant: drop any previously-built hierarchy,
    ! refined levels included - they hang off the old level-1
    ! operator.
    !-----------------------------------------------------------------!

    if (allocated(this % levels))       deallocate(this % levels)
    if (allocated(this % Ac_dense))     deallocate(this % Ac_dense)
    if (allocated(this % ipiv))         deallocate(this % ipiv)
    if (allocated(this % finer_levels)) deallocate(this % finer_levels)
    this % n_finer = 0

    allocate(this % levels(this % max_levels))
    this % levels(1) % A = A

    lev = 1
    coarsening: do

       n = this % levels(lev) % A % num_vertices
       if (n .le. this % coarse_size .or. lev .ge. this % max_levels) exit coarsening

       ! Build the smoother data on this level.
       dinv = 1.0_dp/this % levels(lev) % A % get_diagonal()
       this % levels(lev) % Dinv = dinv

       !--------------------------------------------------------------!
       ! Ask the concrete kind for its squint, then build the
       ! tentative prolongation. The answer is kept only when a
       ! configured export needs it for its composition - storage is
       ! on demand.
       !--------------------------------------------------------------!

       call this % coarsen(lev, agg, naggr)
       if (naggr .ge. n .or. naggr .le. 1) exit coarsening   ! There is no useful coarsening.
       if (this % keeps(lev)) this % levels(lev) % parts = agg

       P0 = tentative_prolongation(agg, naggr)

       rho   = this % spectral_radius(lev)
       omega = (4.0_dp/3.0_dp)/rho

       !--------------------------------------------------------------!
       ! Build the smoothed prolongation
       !
       !    P = (I - omega Dinv A) P0 = P0 - omega Dinv (A P0).
       !
       ! The recipe is multigrid's; it is composed from the general
       ! csr operations.
       !--------------------------------------------------------------!

       AP0 = this % levels(lev) % A % matmat(P0)
       call AP0 % scale_rows(dinv)
       this % levels(lev) % P = P0 % add(1.0_dp, -omega, AP0)

       !--------------------------------------------------------------!
       ! Form the galerkin coarse operator  Ac = P^T (A P). The
       ! restriction is not a thing we keep - the cycle walks P's
       ! edges backward via matvec_transpose - so a scratch transpose
       ! serves this one product and evaporates.
       !--------------------------------------------------------------!

       Pt = this % levels(lev) % P % transpose()
       this % levels(lev+1) % A = Pt % matmat( &
            & this % levels(lev) % A % matmat(this % levels(lev) % P))

       lev = lev + 1

    end do coarsening

    this % nlevel  = lev
    this % ncoarse = this % levels(lev) % A % num_vertices

    !-----------------------------------------------------------------!
    ! Factor the coarsest level with a dense LU, applying a tiny
    ! diagonal shift when singular (e.g. a pure-neumann constant null
    ! space). The LU is self-contained, so the library carries no
    ! LAPACK dependency.
    !-----------------------------------------------------------------!

    call this % levels(lev) % A % to_dense(this % Ac_dense)
    allocate(this % ipiv(this % ncoarse))
    call dense_lu_factor(this % Ac_dense, this % ncoarse, this % ipiv, info)
    if (info .ne. 0) then
       call this % levels(lev) % A % to_dense(this % Ac_dense)
       do i = 1, this % ncoarse
          this % Ac_dense(i,i) = this % Ac_dense(i,i) + 1.0e-10_dp
       end do
       call dense_lu_factor(this % Ac_dense, this % ncoarse, this % ipiv, info)
       if (info .ne. 0) error stop "multigrid: coarse factorization failed"
    end if

    !-----------------------------------------------------------------!
    ! Record the classical V as an ordinary schedule: down the levels
    ! and back, stations numbered 0..nlevel-1..0, with one arrow
    ! between each.
    !
    !    (0) --> (1) --> ... --> (nlevel-1) --> ... --> (1) --> (0)
    !-----------------------------------------------------------------!

    build_v_schedule: block
      integer, allocatable :: level_visits(:)
      integer :: ns
      ns = 2*this % nlevel - 1
      allocate(level_visits(ns))
      do i = 1, this % nlevel
         level_visits(i) = i - 1
      end do
      do i = this % nlevel + 1, ns
         level_visits(i) = 2*this % nlevel - i - 1
      end do
      this % v_cycle_schedule = stored_digraph(ns, &
           & tails = [(i, i = 1, ns-1)], heads = [(i+1, i = 1, ns-1)], &
           & numbers = level_visits)
    end block build_v_schedule

  end subroutine setup

  !===================================================================!
  ! Build n refined levels below the given mesh - the zoom. The
  ! level-1 operator hands over its own graph (simple_graph), the
  ! graph splits into children (the stored_graph refinement, which
  ! adopts the parent map as its partition), and the wheels already
  ! on the shelf do the rest: tentative_prolongation over the parents
  ! builds the interpolation (walked backward, it restricts), and the
  ! lifted operator is A_finer = P A P^T - interpolation alone, no
  ! new discretization, said plainly: a refined level smooths the
  ! same problem through a finer lens. Repeat n times, each level
  ! refining the one before it.
  !===================================================================!

  impure subroutine refine(this, n_finer)

    class(multigrid), intent(inout) :: this
    integer         , intent(in)    :: n_finer

    type(stored_graph)    :: coarser, refined
    type(csr_matrix)      :: A_coarser, PA, Pt
    integer , allocatable :: parent(:)
    integer               :: k, v

    if (n_finer .lt. 1) then
       error stop "multigrid: refine needs at least one level"
    end if
    if (this % nlevel .lt. 1) then
       error stop "multigrid: refine needs a built hierarchy - call setup first"
    end if

    if (allocated(this % finer_levels)) deallocate(this % finer_levels)
    allocate(this % finer_levels(n_finer))
    this % n_finer = n_finer

    coarser   = this % levels(1) % A % simple_graph()
    A_coarser = this % levels(1) % A

    do k = 1, n_finer

       refined = stored_graph(coarser, this % children_per_vertex)

       ! The parent map came in as the refined graph's partition.
       allocate(parent(refined % num_vertices))
       do v = 1, refined % num_vertices
          parent(v) = refined % part_of(v)
       end do

       associate(L => this % finer_levels(k))
         L % P    = tentative_prolongation(parent, coarser % num_vertices)

         !-----------------------------------------------------------!
         ! The lifted operator is A_finer = P A P^T - a scratch
         ! transpose serves the one product and evaporates.
         !-----------------------------------------------------------!

         PA       = L % P % matmat(A_coarser)
         Pt       = L % P % transpose()
         L % A    = PA % matmat(Pt)
         L % Dinv = 1.0_dp/L % A % get_diagonal()
         ! The parent map is kept only when a configured export needs it.
         if (this % keeps(-k)) L % parts = parent
       end associate

       deallocate(parent)
       coarser   = refined
       A_coarser = this % finer_levels(k) % A

    end do

  end subroutine refine

  !===================================================================!
  ! Apply z = M^-1 r: one trip along the classical V, which is just
  ! the schedule setup built; nothing about the V is wired in.
  !===================================================================!

  pure subroutine apply(this, r, z)
    class(multigrid), intent(in)  :: this
    real(dp)        , intent(in)  :: r(:)
    real(dp)        , intent(out) :: z(:)
    call this % apply_cycle(this % v_cycle_schedule, r, z)
  end subroutine apply

  !===================================================================!
  ! Solve A x = b standalone by cycling until ||b - A x||/||b|| falls
  ! below max_tol or max_it cycles pass. iters returns the cycles
  ! taken.
  !===================================================================!

  pure subroutine solve(this, b, x, max_tol, max_it, iters)

    class(multigrid), intent(in)    :: this
    real(dp)        , intent(in)    :: b(:)
    real(dp)        , intent(inout) :: x(:)
    real(dp)        , intent(in)    :: max_tol
    integer         , intent(in)    :: max_it
    integer         , intent(out)   :: iters

    real(dp), allocatable :: r(:), z(:), Ax(:)
    real(dp)              :: bnorm

    if (this % nlevel .lt. 1) then
       error stop "multigrid: solve needs a built hierarchy - call setup first"
    end if

    associate(A => this % levels(1) % A)

      allocate(r(A % num_vertices), z(A % num_vertices), Ax(A % num_vertices))

      bnorm = norm2(b)
      if (bnorm .eq. 0.0_dp) bnorm = 1.0_dp

      do iters = 1, max_it
         call A % matvec(x, Ax)
         r = b - Ax
         if (norm2(r)/bnorm .le. max_tol) exit
         call this % apply(r, z)
         x = x + z
      end do

    end associate

  end subroutine solve

  !===================================================================!
  ! Return the number of levels in the hierarchy; one means a single
  ! grid.
  !===================================================================!

  pure integer function num_levels(this)
    class(multigrid), intent(in) :: this
    num_levels = this % nlevel
  end function num_levels

  !===================================================================!
  ! Follow an injected cycle. The schedule is a directed graph of
  ! stations and moves; the graph reads its own route (source_path,
  ! which audits the shape: one source, one arrow out, no revisit,
  ! full coverage) and hands back the vertex numbers in trip order.
  ! This subroutine only interprets the numbers, which count squints
  ! from the given mesh: level 0 is the mesh, +k is k squints coarser
  ! (nlevel-1 the deepest, solved exactly), -k is k zooms finer.
  !
  !   - At a station, smooth: npre sweeps while the trip is heading
  !     coarser or staying, npost once it is heading finer or ending.
  !     At the deepest coarse level, solve exactly instead.
  !   - A repeated station is the level visited again: two stations
  !     in a row at the same level smooth twice.
  !   - Moving coarser inside the hierarchy restricts the residual
  !     and zeroes the coarser correction; moving back toward the
  !     mesh prolongs the correction and adds it.
  !   - Moving into refined territory lifts the problem through the
  !     interpolation: b and x both ride up. Moving back projects the
  !     improved iterate down. A refined station smooths the lifted
  !     problem - the same equation through a finer lens.
  !
  ! The trip must start and end at level 0, where the residual lives.
  ! An ill-formed schedule stops with a report - the structural
  ! audit, not an extension hook.
  !===================================================================!

  pure subroutine apply_cycle(this, schedule, r, z)

    class(multigrid), intent(in)  :: this
    class(digraph)  , intent(in)  :: schedule
    real(dp)        , intent(in)  :: r(:)
    real(dp)        , intent(out) :: z(:)

    type(level_vectors), allocatable :: state(:)
    logical            , allocatable :: wanted(:)
    integer            , allocatable :: visits(:)
    real(dp)           , allocatable :: res(:), correction(:)
    integer                          :: i, n, l, l_next

    if (this % nlevel .lt. 1) then
       error stop "multigrid: apply needs a built hierarchy - call setup first"
    end if

    ! The graph reads the route; the numbers come back in trip order.
    visits = schedule % source_path()
    n      = size(visits)

    if (visits(1) .ne. 0 .or. visits(n) .ne. 0) then
       error stop "multigrid: the trip must start and end at level 0 - the given mesh"
    end if

    !-----------------------------------------------------------------!
    ! Audit every station against the ladder before anything is
    ! touched, and note which levels the trip actually visits - only
    ! those get travelling state, so a deep refined ladder does not
    ! tax a trip that never goes there.
    !-----------------------------------------------------------------!

    allocate(wanted(-this % n_finer : this % nlevel - 1))
    wanted = .false.
    do i = 1, n
       if (visits(i) .gt. this % nlevel - 1 .or. visits(i) .lt. -this % n_finer) then
          error stop "multigrid: a station names a level outside the ladder"
       end if
       wanted(visits(i)) = .true.
    end do

    allocate(state(-this % n_finer : this % nlevel - 1))
    do l = -this % n_finer, this % nlevel - 1
       if (.not. wanted(l)) cycle
       allocate(state(l) % b(rows_at(this, l)))
       allocate(state(l) % x(rows_at(this, l)))
    end do
    state(0) % b = r
    state(0) % x = 0.0_dp

    do i = 1, n

       l = visits(i)

       l_next = l
       if (i .lt. n) l_next = visits(i+1)

       !--------------------------------------------------------------!
       ! At the station, solve exactly at the deepest coarse level
       ! and smooth anywhere else.
       !--------------------------------------------------------------!

       if (l .eq. this % nlevel - 1) then
          state(l) % x = state(l) % b
          call dense_lu_solve(this % Ac_dense, this % ncoarse, this % ipiv, state(l) % x)
       else if (i .lt. n .and. l_next .ge. l) then
          call this % smooth(l, state(l) % b, state(l) % x, this % npre)
       else
          call this % smooth(l, state(l) % b, state(l) % x, this % npost)
       end if

       if (i .eq. n) exit

       ! Make the move; every destination was audited up front.
       if (l_next .eq. l + 1) then

          if (l .ge. 0) then
             !--------------------------------------------------------!
             ! Heading coarser inside the hierarchy, restrict the
             ! residual and zero the coarser correction.
             !--------------------------------------------------------!

             if (allocated(res)) deallocate(res)
             allocate(res(rows_at(this, l)))
             call this % levels(l+1) % A % matvec(state(l) % x, res)
             res = state(l) % b - res
             ! To restrict is to walk the prolongation's edges backward.
             call this % levels(l+1) % P % matvec_transpose(res, state(l_next) % b)
             state(l_next) % x = 0.0_dp
          else
             !--------------------------------------------------------!
             ! Leaving refined territory, project the improved iterate
             ! back down through the interpolation, against its
             ! arrows.
             !--------------------------------------------------------!

             call this % finer_levels(-l) % P % matvec_transpose(state(l) % x, state(l_next) % x)
          end if

       else if (l_next .eq. l - 1) then

          if (l_next .ge. 0) then
             !--------------------------------------------------------!
             ! Heading back toward the mesh inside the hierarchy,
             ! prolong the correction and add it.
             !--------------------------------------------------------!

             if (allocated(correction)) deallocate(correction)
             allocate(correction(rows_at(this, l_next)))
             call this % levels(l_next+1) % P % matvec(state(l) % x, correction)
             state(l_next) % x = state(l_next) % x + correction
          else
             !--------------------------------------------------------!
             ! Heading into refined territory, lift the whole problem
             ! through the interpolation - b and x both ride up.
             !--------------------------------------------------------!

             call this % finer_levels(-l_next) % P % matvec(state(l) % b, state(l_next) % b)
             call this % finer_levels(-l_next) % P % matvec(state(l) % x, state(l_next) % x)
          end if

       else if (l_next .ne. l) then
          error stop "multigrid: a move goes one level along the ladder, or stays"
       end if

    end do

    z = state(0) % x

  end subroutine apply_cycle

  !===================================================================!
  ! Configure which levels keep their ancestry for export; call this
  ! before setup and refine, with levels counted in squints from the
  ! mesh (positive coarser, negative finer). Storage is on demand: a
  ! step's map is kept only when some requested level needs it for
  ! its composition, and an unconfigured hierarchy keeps nothing.
  !===================================================================!

  pure subroutine export_levels(this, levels)

    class(multigrid), intent(inout) :: this
    integer         , intent(in)    :: levels(:)

    this % exported_levels = levels

  end subroutine export_levels

  !===================================================================!
  ! Report whether the coarsening (+k) or refinement (-k) answer at
  ! step k is needed by any configured export.
  !===================================================================!

  pure logical function keeps(this, k) result(kept)

    class(multigrid), intent(in) :: this
    integer         , intent(in) :: k

    kept = .false.
    if (.not. allocated(this % exported_levels)) return

    if (k .gt. 0) then
       kept = any(this % exported_levels .ge. k)
    else
       kept = any(this % exported_levels .le. k)
    end if

  end function keeps

  !===================================================================!
  ! Return the ancestry map of a configured level. Positive l gives
  ! each mesh cell's coarse vertex after l squints - a per-cell field
  ! on the mesh, ready for the mesh writer, because a quotient IS a
  ! partition of the mesh it came from. Negative l gives each vertex
  ! of the refined level its mesh-cell ancestor - the map sub-cell
  ! work needs to find its way home. Level 0 is the mesh itself:
  ! every cell is its own ancestor. A level whose ancestry was not
  ! configured for export stops with a report.
  !===================================================================!

  pure function ancestors(this, l) result(anc)

    class(multigrid), intent(in) :: this
    integer         , intent(in) :: l

    integer, allocatable :: anc(:)
    integer              :: k, i

    if (this % nlevel .lt. 1) then
       error stop "multigrid: ancestors needs a built hierarchy - call setup first"
    end if

    if (l .eq. 0) then
       anc = [(i, i = 1, this % levels(1) % A % num_vertices)]
       return
    end if

    if (l .gt. this % nlevel - 1 .or. l .lt. -this % n_finer) then
       error stop "multigrid: ancestors of a level outside the ladder"
    end if

    if (l .gt. 0) then

       do k = 1, l
          if (.not. allocated(this % levels(k) % parts)) then
             error stop "multigrid: ancestry was not kept - configure export_levels before setup"
          end if
       end do

       anc = this % levels(1) % parts
       do k = 2, l
          anc = this % levels(k) % parts(anc)
       end do

    else

       do k = 1, -l
          if (.not. allocated(this % finer_levels(k) % parts)) then
             error stop "multigrid: ancestry was not kept - configure export_levels before refine"
          end if
       end do

       anc = this % finer_levels(-l) % parts
       do k = -l - 1, 1, -1
          anc = this % finer_levels(k) % parts(anc)
       end do

    end if

  end function ancestors

  !===================================================================!
  ! Return the rows of the operator at a level counted in squints
  ! from the mesh: 0 is the mesh, +k is coarser via levels(k+1), and
  ! -k is finer.
  !===================================================================!

  pure integer function rows_at(this, l)
    class(multigrid), intent(in) :: this
    integer         , intent(in) :: l
    if (l .ge. 0) then
       rows_at = this % levels(l+1) % A % num_vertices
    else
       rows_at = this % finer_levels(-l) % A % num_vertices
    end if
  end function rows_at

  !===================================================================!
  ! Apply the weighted-jacobi smoother on level lev; it is symmetric,
  ! so it keeps M^-1 SPD:
  !
  !   x <- x + w Dinv (b - A x)
  !===================================================================!

  pure subroutine smooth(this, lev, b, x, nsweep)
    class(multigrid), intent(in)    :: this
    integer         , intent(in)    :: lev   ! Squints from the mesh: 0 is the mesh, +k coarser, -k finer.
    real(dp)        , intent(in)    :: b(:)
    real(dp)        , intent(inout) :: x(:)
    integer         , intent(in)    :: nsweep
    if (lev .ge. 0) then
       call jacobi_sweeps(this % levels(lev+1), this % jacobi_w, b, x, nsweep)
    else
       call jacobi_sweeps(this % finer_levels(-lev), this % jacobi_w, b, x, nsweep)
    end if
  end subroutine smooth

  !===================================================================!
  ! This is the sweep kernel both level directions share: nsweep damped
  ! passes of x <- x + w Dinv (b - A x), each one a matvec over the
  ! level's edges followed by a pointwise (per-vertex) correction.
  !===================================================================!

  pure subroutine jacobi_sweeps(L, w, b, x, nsweep)
    type(multigrid_level), intent(in)    :: L
    real(dp)             , intent(in)    :: w
    real(dp)             , intent(in)    :: b(:)
    real(dp)             , intent(inout) :: x(:)
    integer              , intent(in)    :: nsweep
    real(dp), allocatable :: Ax(:)
    integer :: s
    allocate(Ax(L % A % num_vertices))
    do s = 1, nsweep
       call L % A % matvec(x, Ax)
       x = x + w*L % Dinv*(b - Ax)
    end do
  end subroutine jacobi_sweeps

  !===================================================================!
  ! Estimate the spectral radius of (Dinv A) on level lev - the graph
  ! module's power kernel iterates the damped action, and a
  ! gershgorin floor guards the short iteration from underestimating.
  ! It serves the smoothed-prolongation weight omega = (4/3)/rho.
  !===================================================================!

  impure real(dp) function spectral_radius(this, lev) result(rho)
    class(multigrid), intent(in) :: this
    integer         , intent(in) :: lev
    real(dp), allocatable :: v(:)
    real(dp) :: gers, rowsum
    integer  :: n, i, k

    n = this % levels(lev) % A % num_vertices

    ! Start non-constant and deterministic; the constant is near-null.
    allocate(v(n))
    do i = 1, n
       v(i) = sin(real(i,dp)*0.9_dp) + 0.1_dp
    end do

    call power_iteration(damped_action, n, 15, rho, v0 = v)

    ! Apply the gershgorin floor.
    associate(A => this % levels(lev) % A, Dinv => this % levels(lev) % Dinv)
      gers = 0.0_dp
      do i = 1, n
         rowsum = 0.0_dp
         do k = A % out_xadj(i), A % out_xadj(i+1) - 1
            rowsum = rowsum + abs(Dinv(i)*A % vals(k))
         end do
         gers = max(gers, rowsum)
      end do
    end associate
    rho = max(rho, gers)
    if (rho .le. 0.0_dp) rho = 1.0_dp

  contains

    !=================================================================!
    ! This action applies the damped operator Dinv A; the level is
    ! captured from the host.
    !=================================================================!

    subroutine damped_action(x, y)
      real(dp), intent(in)  :: x(:)
      real(dp), intent(out) :: y(:)
      call this % levels(lev) % A % matvec(x, y)
      y = this % levels(lev) % Dinv * y
    end subroutine damped_action

  end function spectral_radius

  !===================================================================!
  ! Build the tentative prolongation P0 (n x naggr): piecewise
  ! constant over the parts, and column-normalized so the constant
  ! near-null space is represented exactly (orthonormal columns). It
  ! is built from the stamped parts alone - no input matrix - so it
  ! is a plain builder.
  !===================================================================!

  pure function tentative_prolongation(agg, naggr) result(P0)
    integer         , intent(in) :: agg(:), naggr
    type(csr_matrix)             :: P0
    integer , allocatable :: counts(:), row_ptr(:), col_idx(:)
    real(dp), allocatable :: vals(:)
    integer :: n, i

    n = size(agg)

    ! Count the members of each part.
    allocate(counts(naggr))
    counts = 0
    do i = 1, n
       counts(agg(i)) = counts(agg(i)) + 1
    end do

    ! Each row carries one entry: its part's normalized column.
    allocate(row_ptr(n+1), col_idx(n), vals(n))
    do i = 1, n
       row_ptr(i) = i
       col_idx(i) = agg(i)
       vals(i)    = 1.0_dp/sqrt(real(counts(agg(i)), dp))
    end do
    row_ptr(n+1) = n + 1

    P0 = csr_matrix(n, naggr, row_ptr, col_idx, vals)

  end function tentative_prolongation

  !===================================================================!
  ! Factor a dense matrix with a self-contained LU with partial
  ! pivoting for the coarse direct solve; the coarsest level is small,
  ! so a plain O(n^3) factor is fine and keeps the library free of
  ! any LAPACK dependency.
  !===================================================================!

  pure subroutine dense_lu_factor(A, n, ipiv, info)
    integer , intent(in)    :: n
    real(dp), intent(inout) :: A(n,n)
    integer , intent(out)   :: ipiv(n)
    integer , intent(out)   :: info
    integer  :: i, j, k, p
    real(dp) :: amax, t
    info = 0
    do k = 1, n
       ! Pick the pivot row: the largest magnitude in the column.
       p    = k
       amax = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) .gt. amax) then
             amax = abs(A(i,k))
             p = i
          end if
       end do
       ipiv(k) = p
       if (amax .eq. 0.0_dp) then
          info = k
          return
       end if
       ! Swap the pivot row into place.
       if (p .ne. k) then
          do j = 1, n
             t = A(k,j)
             A(k,j) = A(p,j)
             A(p,j) = t
          end do
       end if
       ! Eliminate below the pivot.
       do i = k+1, n
          A(i,k) = A(i,k)/A(k,k)
          do j = k+1, n
             A(i,j) = A(i,j) - A(i,k)*A(k,j)
          end do
       end do
    end do
  end subroutine dense_lu_factor

  !===================================================================!
  ! Solve with the dense LU factor and its recorded pivots.
  !===================================================================!

  pure subroutine dense_lu_solve(A, n, ipiv, b)
    integer , intent(in)    :: n
    real(dp), intent(in)    :: A(n,n)
    integer , intent(in)    :: ipiv(n)
    real(dp), intent(inout) :: b(n)
    integer  :: i, k, p
    real(dp) :: t
    ! Apply the recorded row swaps.
    do k = 1, n
       p = ipiv(k)
       if (p .ne. k) then
          t = b(k)
          b(k) = b(p)
          b(p) = t
       end if
    end do
    ! Solve forward with the unit lower triangle.
    do i = 2, n
       do k = 1, i-1
          b(i) = b(i) - A(i,k)*b(k)
       end do
    end do
    ! Solve backward with the upper triangle.
    do i = n, 1, -1
       do k = i+1, n
          b(i) = b(i) - A(i,k)*b(k)
       end do
       b(i) = b(i)/A(i,i)
    end do
  end subroutine dense_lu_solve

end module interface_multigrid
