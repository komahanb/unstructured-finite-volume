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
! Invalid input: as energy_norms; a chain of other than one block per
! instant; a family whose amplification is not stated.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_march_laws

  use util_precision      , only : dp
  use operation_field     , only : discrete_field
  use operation_manifold  , only : discrete_manifold
  use operation_family    , only : family
  use operation_lagrangian, only : discrete_lagrangian

  implicit none

  private
  public :: march_laws

  type :: march_laws
     real(dp) :: tolerance = 0.0_dp
     ! each law's ratio to its bound, and where it binds
     real(dp) :: accuracy  = 0.0_dp, phase = 0.0_dp, amplitude = 0.0_dp, stability = 0.0_dp, ratio = 0.0_dp
     integer  :: accuracy_instant = 0, accuracy_cell = 0, phase_step = 0, amplitude_step = 0
     integer  :: stability_step = 0, ratio_step = 0
     ! the frequency of the vibration, the highest of the free motions,
     ! the largest energy of the free motion, and the work
     real(dp) :: frequency = 0.0_dp, highest = 0.0_dp, energy = 0.0_dp
     integer  :: work = 0
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
    local  = L_h % local_errors(solution)
    global = L_h % transported_errors(solution, local)
    norm   = L_h % energy_norms(solution, global, highest=this % frequency)
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

  ! the name of the first law that fails, in the order phase,
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
  end function failed

  ! the laws, one to a line, each its ratio to its bound and where it
  ! binds
  function lines(this) result(text)
    class(march_laws), intent(in) :: this
    character(len=:), allocatable :: text
    character(len=160) :: line
    character(len=1), parameter :: newline = achar(10)
    text = ''
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

end module operation_march_laws
