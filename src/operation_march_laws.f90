!=====================================================================!
! THE LAWS OF A MARCH: the constraints of the solver Lagrangian of the
! time discretization,
!
!    L_s = W(s) + sum_j mu_j c_j(s),    c_j(s) <= 0,
!
! each evaluated at the design s of a march solved: its instants
! t_1 < .. < t_N and the family of each step. W is the work, the
! solves of the steps' stages.
!
! One tolerance eps is given. Each law is stated as a ratio to its
! bound, so that it is met where the ratio is at most one:
!
!   accuracy   the largest over the cells c and instants k of the
!              energy norm of the transported error, per unit measure
!              of the cell, N_c(e(t_k)) / |c|^1/2, over eps times the
!              root of the largest energy of the free motion per unit
!              measure of the region, (max_k E(t_k) / |Omega|)^1/2;
!   phase      sum over the steps of |arg R(i w_hi h_k) - w_hi h_k|,
!              over eps: the phase lost by the vibration;
!   amplitude  sum over the steps of -ln|R(i w_hi h_k)|, over eps: the
!              amplitude the families take from the vibration;
!   stability  sum over the steps of ln|R(i w h_k)|, at w_hi and at the
!              highest frequency of the free motions, over eps: the
!              growth the families give a mode;
!   ratio      the largest h_k / (r h_(k-1)), r the family's bound on
!              the ratio of its steps: its zero-stability; zero where
!              no family of the march bounds it.
!
! R is each step's amplification on y' = lambda y (family), w_hi the
! frequency of the vibration and E the energy of the free motion
! (measure_free_motion), e the transported error (transported_errors)
! and N its energy norm (energy_norms).
! A law reads the step or instant where it binds.
! The first that fails, in the order phase, amplitude, stability,
! ratio, accuracy, is named: a frequency the instants cannot hold
! fails the phase law before any estimate of the error is read, and
! the accuracy is then not evaluated, its ratio -1.
! The frequency, the highest frequency and the energy may be given:
! designs compared with each other are measured by the same ones, the
! largest measured over them, so that a family that damps the
! vibration out does not lower the bounds it is measured against.
!
! THE SELECTION (select_march) chooses the instants and the family of
! a march to meet these laws with the least work, by rounds of
! marches; the solve law fails at a design whose march does not
! converge, and is then named first.
!
! Invalid input: as energy_norms; a chain of other than one block per
! instant; a family whose amplification is not stated.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_march_laws

  use util_precision      , only : dp
  use operation_field     , only : discrete_field
  use operation_manifold  , only : discrete_manifold
  use operation_family    , only : family, chain
  use operation_lagrangian, only : discrete_lagrangian
  use operation_minimization, only : discrete_minimizer
  use util_verbosity      , only : verbosity, set_verbosity

  implicit none

  private
  public :: march_laws, select_march

  type :: march_laws
     real(dp) :: tolerance = 0.0_dp
     ! each law's ratio to its bound, and where it binds
     real(dp) :: accuracy  = 0.0_dp, phase = 0.0_dp, amplitude = 0.0_dp, stability = 0.0_dp, ratio = 0.0_dp
     ! the solve, huge where the march does not converge
     real(dp) :: solve = 0.0_dp
     integer  :: accuracy_instant = 0, accuracy_cell = 0, phase_step = 0, amplitude_step = 0
     integer  :: stability_step = 0, ratio_step = 0
     ! the frequency of the vibration, the highest of the free motions,
     ! the largest energy of the free motion, and the work
     real(dp) :: frequency = 0.0_dp, highest = 0.0_dp, energy = 0.0_dp
     integer  :: work = 0
     ! the largest over the cells of the energy norm of each step's
     ! local error, at the instant ending it, where the accuracy is read
     real(dp), allocatable :: local_norm(:)
   contains
     procedure :: failed
     procedure :: lines
  end type march_laws

  interface march_laws
     module procedure create
  end interface march_laws

contains

  function create(L_h, solution, tolerance, frequency, highest, energy) result(this)

    type(discrete_lagrangian), intent(in) :: L_h
    type(discrete_field)     , intent(in) :: solution
    real(dp)                 , intent(in) :: tolerance
    real(dp)                 , intent(in), optional :: frequency, highest, energy
    type(march_laws) :: this

    type(discrete_manifold) :: points
    type(family) :: f
    real(dp), allocatable :: local(:,:,:), global(:,:,:), norm(:,:), energies(:), measure(:)
    real(dp) :: h, h_before, step_phase, step_amplitude, grown, top_grown, region
    real(dp) :: phase_largest, amplitude_largest, stability_largest
    integer :: n, b, k, c

    if (tolerance <= 0.0_dp) error stop 'operation_march_laws: the tolerance is positive'
    select type (on => solution % on)
    type is (discrete_manifold)
       points = on
    class default
       error stop 'operation_march_laws: the laws are of a solution at the points of a manifold'
    end select
    n = points % num_instants()
    if (L_h % schemes % num_blocks() /= n) then
       error stop 'operation_march_laws: the laws are of a chain of one block per instant'
    end if
    this % tolerance = tolerance

    ! the free motion, measured or given, and the families' laws step
    ! by step
    if (present(frequency) .and. present(highest) .and. present(energy)) then
       this % frequency = frequency
       this % highest   = highest
       this % energy    = energy
    else
       call L_h % measure_free_motion(solution, frequency=this % frequency, highest=this % highest, energy=energies)
       this % energy = maxval(energies)
       if (present(frequency)) this % frequency = frequency
       if (present(highest))   this % highest   = highest
       if (present(energy))    this % energy    = energy
    end if
    phase_largest     = -1.0_dp
    amplitude_largest = -huge(1.0_dp)
    stability_largest = -huge(1.0_dp)
    top_grown         = 0.0_dp
    h_before = 0.0_dp
    do b = 2, n
       f = L_h % schemes % scheme(b)
       h = points % instant(b) - points % instant(b - 1)
       this % work = this % work + f % num_stages()
       step_phase     = abs(atan2(aimag(f % amplification(cmplx(0.0_dp, this % frequency * h, dp))), &
            & real(f % amplification(cmplx(0.0_dp, this % frequency * h, dp)))) - this % frequency * h)
       step_amplitude = -log(abs(f % amplification(cmplx(0.0_dp, this % frequency * h, dp))))
       grown = max(log(abs(f % amplification(cmplx(0.0_dp, this % frequency * h, dp)))), &
            &      log(abs(f % amplification(cmplx(0.0_dp, this % highest * h, dp)))))
       this % phase     = this % phase + step_phase
       this % amplitude = this % amplitude + step_amplitude
       top_grown        = top_grown + grown
       if (step_phase > phase_largest) then
          phase_largest   = step_phase
          this % phase_step = b
       end if
       if (step_amplitude > amplitude_largest) then
          amplitude_largest   = step_amplitude
          this % amplitude_step = b
       end if
       if (grown > stability_largest) then
          stability_largest   = grown
          this % stability_step = b
       end if
       if (h_before > 0.0_dp .and. f % step_ratio_bound() < huge(1.0_dp)) then
          if (h / (f % step_ratio_bound() * h_before) > this % ratio) then
             this % ratio      = h / (f % step_ratio_bound() * h_before)
             this % ratio_step = b
          end if
       end if
       h_before = h
    end do
    this % phase     = this % phase / tolerance
    this % amplitude = this % amplitude / tolerance
    this % stability = top_grown / tolerance

    ! the accuracy: each cell's norm per unit measure against the
    ! energy of the free motion per unit measure of the region; not
    ! evaluated where the phase or the amplitude fails, the estimate
    ! being blind to the frequency those laws cannot represent
    if (this % phase > 1.0_dp .or. this % amplitude > 1.0_dp) then
       this % accuracy = -1.0_dp
       return
    end if
    global = L_h % transported_errors(solution, local)
    norm   = L_h % energy_norms(solution, global, highest=this % frequency)
    this % local_norm = maxval(L_h % energy_norms(solution, local, highest=this % frequency), dim=1)
    measure = cell_measures(L_h, points, size(norm, 1))
    region  = sum(measure)
    if (this % energy <= 0.0_dp) then
       this % accuracy = huge(1.0_dp)
    else
       this % accuracy = 0.0_dp
       do k = 1, size(norm, 2)
          do c = 1, size(norm, 1)
             if (norm(c, k) / sqrt(measure(c)) > this % accuracy) then
                this % accuracy         = norm(c, k) / sqrt(measure(c))
                this % accuracy_instant = k
                this % accuracy_cell    = c
             end if
          end do
       end do
       this % accuracy = this % accuracy / (tolerance * sqrt(this % energy / region))
    end if

  end function create

  ! the measure of each cell the norms are of: of each element of a law
  ! projected on the cells, the sum of its quadrature weights; one for
  ! a law at the points of one cell
  function cell_measures(L_h, points, cells) result(measure)
    type(discrete_lagrangian), intent(in) :: L_h
    type(discrete_manifold)  , intent(in) :: points
    integer                  , intent(in) :: cells
    real(dp), allocatable :: measure(:)
    integer :: e
    allocate(measure(cells), source=1.0_dp)
    if (allocated(L_h % elements)) then
       do e = 1, cells
          measure(e) = sum(L_h % elements % weight(:, e))
       end do
    else if (points % num_cells() > 1) then
       error stop 'operation_march_laws: the accuracy of a law at the points of several cells is not stated'
    end if
  end function cell_measures

  ! the name of the first law that fails, in the order solve, phase,
  ! amplitude, stability, ratio, accuracy; empty where every law is met
  function failed(this) result(name)
    class(march_laws), intent(in) :: this
    character(len=:), allocatable :: name
    name = ''
    if (this % accuracy > 1.0_dp) name = 'accuracy'
    if (this % ratio > 1.0_dp) name = 'ratio'
    if (this % stability > 1.0_dp) name = 'stability'
    if (this % amplitude > 1.0_dp) name = 'amplitude'
    if (this % phase > 1.0_dp) name = 'phase'
    if (this % solve > 1.0_dp) name = 'solve'
  end function failed

  ! the laws, one to a line, each its ratio to its bound and where it
  ! binds
  function lines(this) result(text)
    class(march_laws), intent(in) :: this
    character(len=:), allocatable :: text
    character(len=160) :: line
    character(len=1), parameter :: newline = achar(10)
    text = ''
    if (this % solve > 1.0_dp) then
       text = 'solve      the march does not converge'
       return
    end if
    write(line, '(a, es10.3, a, es10.3, a, es10.3, a, i0)') 'frequency of the vibration ', this % frequency, &
         & ', highest ', this % highest, ', largest energy of the free motion ', this % energy, ', work ', this % work
    text = text // trim(line) // newline
    write(line, '(a, es10.3, a, i0)') 'phase      ', this % phase, '  largest at the step ending at instant ', &
         & this % phase_step
    text = text // trim(line) // newline
    write(line, '(a, es10.3, a, i0)') 'amplitude  ', this % amplitude, '  largest at the step ending at instant ', &
         & this % amplitude_step
    text = text // trim(line) // newline
    write(line, '(a, es10.3, a, i0)') 'stability  ', this % stability, '  largest at the step ending at instant ', &
         & this % stability_step
    text = text // trim(line) // newline
    write(line, '(a, es10.3, a, i0)') 'ratio      ', this % ratio, '  at the step ending at instant ', this % ratio_step
    text = text // trim(line) // newline
    if (this % accuracy < 0.0_dp) then
       line = 'accuracy   not evaluated: the phase or the amplitude fails'
    else
       write(line, '(a, es10.3, a, i0, a, i0)') 'accuracy   ', this % accuracy, '  at instant ', &
            & this % accuracy_instant, ', cell ', this % accuracy_cell
    end if
    text = text // trim(line)
  end function lines

  !===================================================================!
  ! THE SELECTION OF A MARCH: the instants and the family of a march
  ! chosen to meet the laws at the tolerance with the least work, by
  ! rounds.
  !
  ! Each round marches the design and reads its laws against the
  ! largest frequency, highest frequency and energy measured over the
  ! rounds (march_laws).
  ! Where every law is met, the design is selected.
  ! Where the phase or the amplitude fails, the count each family of the
  ! pool needs on steps of one length follows from its amplification
  ! alone (uniform_count); the family of least work, the count times its
  ! stages, is taken at that count.
  ! Where the accuracy fails, the instants are placed by equidistributing
  ! the energy norms of the steps' local errors, eta_k ~ C_k h_k^(p+1):
  ! the density of the steps is proportional to C^(1/(p+1)), and the
  ! count is the present one times the accuracy's ratio to the power
  ! 1/p.
  ! p is the family's order until two rounds of the family give the
  ! order observed between them, log(A_1 / A_2) / log(N_2 / N_1), A the
  ! accuracy's ratio and N the count.
  ! Where the observed order is not positive, refining does not help,
  ! and the family of next least work is taken.
  ! Where the march does not converge, the solve law fails for that
  ! design, and the family of next least work is taken.
  !
  ! With count given, the count is fixed: the instants are placed once
  ! by equidistribution, and the laws of that design are read.
  !
  ! The selection ends where every law is met, where a round would
  ! measure a design measured before, or where no family of the pool is
  ! left; laws % failed() then names the law that fails.
  ! The march's minimizer by must report a singular block rather than
  ! stop (dense_direct, singular_reported), for a design that fails to
  ! end its round.
  ! Invalid input: as march_laws; an empty pool, or one with a family of
  ! several steps.
  !===================================================================!

  subroutine select_march(L_h, estimate, solution, laws, tolerance, pool, count, by)

    type(discrete_lagrangian), intent(inout) :: L_h
    type(discrete_field)     , intent(in)    :: estimate
    type(discrete_field)     , intent(out)   :: solution
    type(march_laws)         , intent(out)   :: laws
    real(dp)                 , intent(in)    :: tolerance
    type(family)             , intent(in)    :: pool(:)
    integer                  , intent(in), optional :: count
    type(discrete_minimizer) , intent(in), optional :: by

    type(discrete_field) :: placed
    type(march_laws) :: unsolved
    real(dp), allocatable :: t(:), energies(:), eta(:), t_measured(:), eta_measured(:)
    integer , allocatable :: tried(:,:)
    real(dp) :: frequency, highest, energy, f, hi, accuracy_before, p, longest
    integer :: n, n_before, g, g_before, g_measured, round, k, level
    logical :: met

    if (size(pool) == 0) error stop 'operation_march_laws: a selection reads a pool of one family at least'
    if (any([(pool(k) % history_depth(2) > 1, k = 1, size(pool))])) then
       error stop 'operation_march_laws: a selection gives every step one family, so its pool is of families of &
            &one step; a family of several steps would read history before the first instant'
    end if
    t = L_h % points % instant
    g = 1
    frequency  = 0.0_dp
    highest    = 0.0_dp
    energy     = 0.0_dp
    n_before   = 0
    g_before   = 0
    g_measured = 0
    accuracy_before = 0.0_dp
    allocate(tried(2, 0))
    level = verbosity
    call set_verbosity(0)
    round = 0
    do
       round = round + 1
       if (any(tried(1, :) == size(t) .and. tried(2, :) == g)) exit
       tried = reshape([tried, size(t), g], [2, size(tried, 2) + 1])
       call design_set(L_h, t, pool(g))
       placed = initial_estimate(L_h, estimate)
       call L_h % minimize(placed, solution, continued=.true., by=by, converged=met)
       if (.not. met) then
          ! the solve law fails at this design
          laws = unsolved
          laws % tolerance = tolerance
          laws % solve     = huge(1.0_dp)
          if (level >= 1) print '(a, i0, a, i0, a, i0, a)', 'selection round ', round, ': family ', g, &
               & '  instants ', size(t), '  the march does not converge'
          g = next_family(g)
          if (g == 0) exit
          cycle
       end if
       call L_h % measure_free_motion(solution, frequency=f, highest=hi, energy=energies)
       frequency = max(frequency, f)
       highest   = max(highest, hi)
       energy    = max(energy, maxval(energies))
       laws = march_laws(L_h, solution, tolerance, frequency=frequency, highest=highest, energy=energy)
       if (level >= 1) print '(a, i0, a, i0, a, i0, a, 5es10.2, a, a)', 'selection round ', round, ': family ', g, &
            & '  instants ', size(t), '  phase, amplitude, stability, ratio, accuracy ', laws % phase, &
            & laws % amplitude, laws % stability, laws % ratio, laws % accuracy, '  failed ', laws % failed()
       if (len(laws % failed()) == 0) exit
       if (present(count) .and. round > 1) exit
       if (laws % accuracy < 0.0_dp .or. laws % stability > 1.0_dp .or. laws % ratio > 1.0_dp) then
          ! the phase, the amplitude or the stability: the family of
          ! least work at the step their bound allows, its instants
          ! placed again where its accuracy was measured, else of one
          ! length
          g = least_work()
          if (g == 0) exit
          longest = longest_step(pool(g))
          if (g == g_measured) then
             t = equidistributed(t_measured, eta_measured, real(pool(g) % accuracy(), dp), size(t_measured), &
                  & longest, count)
          else
             n = uniform_count(pool(g), frequency, highest, t(1), t(size(t)), tolerance)
             if (present(count)) n = count
             t = [(t(1) + (t(size(t)) - t(1)) * real(k - 1, dp) / real(n - 1, dp), k = 1, n)]
          end if
          cycle
       end if
       ! the accuracy: the order observed between two rounds of the
       ! family, the count to meet it, and the instants equidistributed,
       ! no step longer than the phase, the amplitude and the stability
       ! allow
       p = real(pool(g) % accuracy(), dp)
       if (g == g_before .and. n_before > 0 .and. n_before /= size(t)) then
          p = log(accuracy_before / laws % accuracy) / log(real(size(t) - 1, dp) / real(n_before - 1, dp))
          if (p <= 0.0_dp) then
             g = next_family(g)
             if (g == 0) exit
             cycle
          end if
       end if
       g_before        = g
       n_before        = size(t)
       accuracy_before = laws % accuracy
       eta = laws % local_norm
       g_measured   = g
       t_measured   = t
       eta_measured = eta
       n = ceiling((size(t) - 1) * laws % accuracy**(1.0_dp / p)) + 1
       t = equidistributed(t, eta, p, n, longest_step(pool(g)), count)
    end do
    call set_verbosity(level)

  contains

    ! the longest step of the family on which the phase, the amplitude
    ! and the stability laws are met on steps of one length
    real(dp) function longest_step(f)
      type(family), intent(in) :: f
      longest_step = (t(size(t)) - t(1)) / real(uniform_count(f, frequency, highest, t(1), t(size(t)), tolerance) &
           & - 1, dp)
    end function longest_step

    ! the next family of the pool not tried, after g; zero where none is
    ! left
    integer function next_family(g) result(next)
      integer, intent(in) :: g
      integer :: i
      next = 0
      do i = 1, size(pool)
         if (i == g) cycle
         if (any(tried(2, :) == i)) cycle
         next = i
         return
      end do
    end function next_family

    ! the family of least work, the count its amplification requires times its
    ! stages; zero where every family has been tried at that count
    integer function least_work() result(chosen)
      integer :: i, m
      real(dp) :: w, least
      chosen = 0
      least  = huge(1.0_dp)
      do i = 1, size(pool)
         m = uniform_count(pool(i), frequency, highest, t(1), t(size(t)), tolerance)
         if (i /= g_measured .and. any(tried(1, :) == m .and. tried(2, :) == i)) cycle
         w = real(m - 1, dp) * pool(i) % num_stages()
         if (w < least) then
            least  = w
            chosen = i
         end if
      end do
    end function least_work

  end subroutine select_march

  !===================================================================!
  ! THE COUNT OF INSTANTS a family needs on steps of one length over
  ! [a, b] for the phase, the amplitude and the stability laws at the
  ! frequency w and the highest frequency: the least n at which the
  ! sums over the n - 1 steps of |arg R(i w h) - w h|, of -ln|R(i w h)|
  ! and of ln|R(i w h)| at w and at the highest are at most the
  ! tolerance. The count doubles until they are, then is bisected.
  !===================================================================!

  integer function uniform_count(f, w, top, a, b, tolerance) result(n)
    type(family), intent(in) :: f
    real(dp)    , intent(in) :: w, top, a, b, tolerance
    integer :: low, high, middle
    high = 2
    do while (.not. met_at(high))
       high = 2 * high
    end do
    low = max(2, high / 2)
    do while (high - low > 1)
       middle = (low + high) / 2
       if (met_at(middle)) then
          high = middle
       else
          low = middle
       end if
    end do
    n = high
  contains
    logical function met_at(m)
      integer, intent(in) :: m
      real(dp) :: h
      complex(dp) :: r, r_top
      h = (b - a) / real(m - 1, dp)
      r     = f % amplification(cmplx(0.0_dp, w * h, dp))
      r_top = f % amplification(cmplx(0.0_dp, top * h, dp))
      met_at = (m - 1) * abs(atan2(aimag(r), real(r)) - w * h) <= tolerance .and. &
           &   (m - 1) * (-log(abs(r))) <= tolerance .and. &
           &   (m - 1) * max(log(abs(r)), log(abs(r_top))) <= tolerance
    end function met_at
  end function uniform_count

  ! the design: the instants t, the family f at every step
  subroutine design_set(L_h, t, f)
    type(discrete_lagrangian), intent(inout) :: L_h
    real(dp)                 , intent(in)    :: t(:)
    type(family)             , intent(in)    :: f
    integer :: k
    L_h % points % instant = t
    L_h % schemes = chain([(f, k = 1, size(t))], from=[(k, k = 1, size(t))])
  end subroutine design_set

  ! the estimate on the instants of L_h: the estimate's first instant at
  ! every instant, which a march continued from instant to instant reads
  ! alone
  function initial_estimate(L_h, estimate) result(placed)
    type(discrete_lagrangian), intent(in) :: L_h
    type(discrete_field)     , intent(in) :: estimate
    type(discrete_field) :: placed
    real(dp), allocatable :: values(:,:)
    integer :: k, c, ncells
    ncells = L_h % points % num_cells()
    allocate(values(L_h % points % num_points(), size(estimate % value, 2)))
    do k = 1, L_h % points % num_instants()
       do c = 1, ncells
          values((k - 1) * ncells + c, :) = estimate % value(c, :)
       end do
    end do
    placed = discrete_field(L_h % points, values, estimate % name)
  end function initial_estimate

  !===================================================================!
  ! THE INSTANTS EQUIDISTRIBUTED over [t_1, t_N]: their density on step
  ! k is proportional to C_k^(1/(p+1)), C_k = eta_k / h_k^(p+1) the
  ! error density of the step whose local error is eta_k at its instant
  ! k, scaled to n instants, and at least 1 / longest: no step is longer
  ! than longest. A step whose local error is zero takes the least
  ! positive density of the others. The instants are the quantiles of
  ! the integral of the density, linear within each step, and their
  ! count is that integral's, at least n; count given, it is count.
  !===================================================================!

  function equidistributed(t, eta, p, n, longest, count) result(placed)
    real(dp), intent(in) :: t(:), eta(:), p, longest
    integer , intent(in) :: n
    integer , intent(in), optional :: count
    real(dp), allocatable :: placed(:)
    real(dp), allocatable :: density(:), integral(:)
    real(dp) :: target, least_density
    integer :: k, j, m, total
    m = size(t)
    allocate(density(m), source=0.0_dp)
    do k = 2, m
       density(k) = (eta(k) / (t(k) - t(k - 1))**(p + 1.0_dp))**(1.0_dp / (p + 1.0_dp))
    end do
    if (any(density(2:) > 0.0_dp)) then
       least_density = minval(density(2:), mask=density(2:) > 0.0_dp)
    else
       least_density = 1.0_dp
    end if
    where (density(2:) <= 0.0_dp) density(2:) = least_density
    density(2:) = density(2:) * real(n - 1, dp) / sum(density(2:) * (t(2:) - t(:m - 1)))
    density(2:) = max(density(2:), 1.0_dp / longest)
    allocate(integral(m), source=0.0_dp)
    do k = 2, m
       integral(k) = integral(k - 1) + density(k) * (t(k) - t(k - 1))
    end do
    total = max(n, ceiling(integral(m)) + 1)
    if (present(count)) total = count
    allocate(placed(total))
    placed(1)     = t(1)
    placed(total) = t(m)
    k = 2
    do j = 2, total - 1
       target = integral(m) * real(j - 1, dp) / real(total - 1, dp)
       do while (integral(k) < target)
          k = k + 1
       end do
       placed(j) = t(k - 1) + (target - integral(k - 1)) / density(k)
    end do
  end function equidistributed

end module operation_march_laws
