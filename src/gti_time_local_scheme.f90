!=====================================================================!
! GTI TIME-LOCAL SCHEME MOTIF (LOCAL SEAT)
!
! A time-local motif is a coefficient rule that builds the state
! tuple of one evaluation point from a finite set of time samples:
!
!      q      =  sum_i a_i q_i
!      qdot   =  sum_i b_i q_i
!      qddot  =  sum_i c_i q_i
!
! one rule per state component, one weight per sample. The motif
! knows coefficient rows and nothing else: not BDF, not DIRK, not
! ABM, not RK, not adaptivity - those are later phases' ways of
! CHOOSING rows, while this seat only applies them. No time graph
! is built, and nothing is solved: the evaluator combines samples
! into a point and asks the residual form its value through the
! existing form evaluator.
!
! The laws of the input, refused loudly before any construction:
!
!      at least one sample
!      every rule names a state component from the GTI_STATE_*
!        vocabulary, no component twice
!      every rule carries weights, one per sample
!      every sample provides q, and all sample q vectors agree
!        in shape
!
! The built point:
!
!      state    seats up to the largest ruled component; only
!               ruled components have occupants
!      design   copied in as given
!      time     the requested time; window, step, and stage stay
!               zero in this phase, for no time graph names them
!
! Output component fields are minted by copying a sample q field
! and replacing its values - the copy carries the domain identity,
! so q, qdot, and qddot all live where the samples live. One
! concrete field type is supported; any other dynamic type is
! refused loudly rather than guessed at.
!
! The evaluator carries nothing: no form, no graph, no scheme
! table, no solver, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_local_schemes

  use iso_fortran_env      , only : dp => REAL64
  use class_graph_field    , only : field
  use gti_state_bundles    , only : gti_state_bundle, gti_field_slot, &
       & GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  use gti_design_bundles   , only : gti_design_bundle
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_value_buffers    , only : gti_value_buffer
  use gti_form_interface   , only : gti_differentiable_form
  use gti_form_evaluators  , only : gti_form_evaluator

  implicit none

  private
  public :: gti_time_sample
  public :: gti_time_component_rule
  public :: gti_time_motif
  public :: gti_time_local_residual_evaluator

  !===================================================================!
  ! One time sample: a state bundle holding this sample's q, and
  ! the time it was taken.
  !===================================================================!

  type :: gti_time_sample

     type(gti_state_bundle) :: state
     real(dp)               :: time = 0.0_dp

  end type gti_time_sample

  !===================================================================!
  ! One coefficient row: which state component it builds, and the
  ! weight each sample contributes.
  !===================================================================!

  type :: gti_time_component_rule

     integer :: state_component = GTI_STATE_Q
     real(dp), allocatable :: weight(:)

  end type gti_time_component_rule

  !===================================================================!
  ! The motif: the coefficient rows of one time-local scheme.
  ! The types keep their public singular names; Fortran denies a
  ! type its host module's name, so the module speaks in the
  ! plural.
  !===================================================================!

  type :: gti_time_motif

     type(gti_time_component_rule), allocatable :: rule(:)

   contains

     procedure :: size          => motif_size
     procedure :: has_component => motif_has_component

  end type gti_time_motif

  !===================================================================!
  ! The stateless verb-pair: build a point from motif and samples,
  ! or go on and ask the form its value there.
  !===================================================================!

  type :: gti_time_local_residual_evaluator

   contains

     procedure :: build_point
     procedure :: value

  end type gti_time_local_residual_evaluator

contains

  !===================================================================!
  ! How many coefficient rows the motif carries.
  !===================================================================!

  pure function motif_size(this) result(nrules)

    class(gti_time_motif), intent(in) :: this
    integer :: nrules

    if (allocated(this % rule)) then
       nrules = size(this % rule)
    else
       nrules = 0
    end if

  end function motif_size

  !===================================================================!
  ! Does some row build this state component?
  !===================================================================!

  pure function motif_has_component(this, state_component) result(has)

    class(gti_time_motif), intent(in) :: this
    integer              , intent(in) :: state_component
    logical :: has

    integer :: k

    has = .false.
    if (.not. allocated(this % rule)) return

    do k = 1, size(this % rule)
       if (this % rule(k) % state_component == state_component) then
          has = .true.
          return
       end if
    end do

  end function motif_has_component

  !===================================================================!
  ! Build the local evaluation point: admit the motif and the
  ! samples, combine sample q vectors through the coefficient
  ! rows, seat the built components, carry the design, stamp the
  ! time. Window, step, and stage stay zero - no time graph names
  ! them in this phase.
  !===================================================================!

  subroutine build_point(this, motif, samples, design, time, point)

    class(gti_time_local_residual_evaluator), intent(in)    :: this
    type(gti_time_motif)                    , intent(in)    :: motif
    type(gti_time_sample)                   , intent(in)    :: samples(:)
    type(gti_design_bundle)                 , intent(in)    :: design
    real(dp)                                , intent(in)    :: time
    type(gti_evaluation_point)              , intent(inout) :: point

    type(gti_state_bundle) :: built
    real(dp), allocatable  :: qs(:,:), vec(:)
    integer :: k, seats, seat

    call require_lawful_motif(motif, samples)

    call gather_sample_q(samples, qs)

    !----------------------------------------------------------------!
    ! Seats up to the largest ruled component; only ruled
    ! components gain occupants.
    !----------------------------------------------------------------!

    seats = 0
    do k = 1, motif % size()
       seats = max(seats, motif % rule(k) % state_component + 1)
    end do

    allocate(built % component(seats))

    do k = 1, motif % size()
       vec  = matmul(qs, motif % rule(k) % weight)
       seat = 1 + motif % rule(k) % state_component
       call fill_component_slot(samples(1) % state % component(1 + GTI_STATE_Q), &
            & vec, component_label(motif % rule(k) % state_component), &
            & built % component(seat))
    end do

    point % state     = built
    point % design    = design
    point % time      = time
    point % window_id = 0
    point % step_id   = 0
    point % stage_id  = 0

  end subroutine build_point

  !===================================================================!
  ! Build the point and ask the residual form its value there,
  ! held to its declared output shape by the form evaluator.
  !===================================================================!

  subroutine value(this, form, motif, samples, design, time, output)

    class(gti_time_local_residual_evaluator), intent(in)    :: this
    class(gti_differentiable_form)          , intent(in)    :: form
    type(gti_time_motif)                    , intent(in)    :: motif
    type(gti_time_sample)                   , intent(in)    :: samples(:)
    type(gti_design_bundle)                 , intent(in)    :: design
    real(dp)                                , intent(in)    :: time
    type(gti_value_buffer)                  , intent(inout) :: output

    type(gti_evaluation_point) :: point
    type(gti_form_evaluator)   :: evaluator

    call this % build_point(motif, samples, design, time, point)

    call evaluator % value(form, point, output)

  end subroutine value

  !===================================================================!
  ! The motif laws, refused in order: samples exist; every row
  ! speaks the GTI_STATE_* vocabulary; no component is built
  ! twice; every row carries weights, one per sample.
  !===================================================================!

  pure subroutine require_lawful_motif(motif, samples)

    type(gti_time_motif) , intent(in) :: motif
    type(gti_time_sample), intent(in) :: samples(:)

    integer :: k, j

    if (size(samples) < 1) then
       error stop 'gti_time_local_scheme: at least one sample is required'
    end if

    do k = 1, motif % size()
       select case (motif % rule(k) % state_component)
       case (GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT)
          continue
       case default
          error stop 'gti_time_local_scheme: unknown state component is refused'
       end select
    end do

    do k = 1, motif % size()
       do j = k + 1, motif % size()
          if (motif % rule(k) % state_component == &
               & motif % rule(j) % state_component) then
             error stop 'gti_time_local_scheme: duplicate state component rule is refused'
          end if
       end do
    end do

    do k = 1, motif % size()
       if (.not. allocated(motif % rule(k) % weight)) then
          error stop 'gti_time_local_scheme: a rule has weights'
       end if
       if (size(motif % rule(k) % weight) /= size(samples)) then
          error stop 'gti_time_local_scheme: rule weight count matches sample count'
       end if
    end do

  end subroutine require_lawful_motif

  !===================================================================!
  ! The sample laws, proven while gathering: every sample provides
  ! q, and every q vector has the shape of the first. The gathered
  ! matrix holds sample i's values in column i.
  !===================================================================!

  pure subroutine gather_sample_q(samples, qs)

    type(gti_time_sample) , intent(in)  :: samples(:)
    real(dp), allocatable , intent(out) :: qs(:,:)

    real(dp), allocatable :: v(:)
    integer :: i, n

    do i = 1, size(samples)
       if (.not. samples(i) % state % has_component(GTI_STATE_Q)) then
          error stop 'gti_time_local_scheme: every sample provides q'
       end if
    end do

    call samples(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(v)
    n = size(v)

    allocate(qs(n, size(samples)))
    qs(:, 1) = v

    do i = 2, size(samples)
       call samples(i) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(v)
       if (size(v) /= n) then
          error stop 'gti_time_local_scheme: sample q shapes match'
       end if
       qs(:, i) = v
    end do

  end subroutine gather_sample_q

  !===================================================================!
  ! Mint one built component: copy the prototype q field - the
  ! copy carries the domain identity - rename it for the component
  ! it now holds, and replace its values. One concrete field type
  ! is supported; any other dynamic type is refused, never guessed.
  !===================================================================!

  subroutine fill_component_slot(prototype, vec, label, slot)

    type(gti_field_slot), intent(in)    :: prototype
    real(dp)            , intent(in)    :: vec(:)
    character(len=*)    , intent(in)    :: label
    type(gti_field_slot), intent(inout) :: slot

    type(field) :: minted

    select type (proto => prototype % value)
    type is (field)
       minted = proto
       minted % label = label
       call minted % set_real_vector(vec)
       slot % value = minted
    class default
       error stop 'gti_time_local_scheme: output field dynamic type is supported'
    end select

  end subroutine fill_component_slot

  !===================================================================!
  ! The component's name, in the established vocabulary.
  !===================================================================!

  pure function component_label(state_component) result(label)

    integer, intent(in) :: state_component
    character(len=:), allocatable :: label

    select case (state_component)
    case (GTI_STATE_Q)
       label = 'q'
    case (GTI_STATE_QDOT)
       label = 'qdot'
    case default
       label = 'qddot'
    end select

  end function component_label

end module gti_time_local_schemes
