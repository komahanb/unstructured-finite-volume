!=====================================================================!
! GTI TIME-LOCAL UNKNOWN RESIDUAL PROBLEM
!
! The first unknown seat: one sample of a local time motif is
! declared unknown, a trial q is injected into it, and the
! residual answers at that trial -
!
!      residual(q_trial) = R( U(q_trial), xi, t )
!
! through the existing time-local seat, with U built by the motif
! rows exactly as for any sample set. This is the local state of
! the residual problem at one trial, and nothing more: no Newton,
! no Jacobian, no linear solve, no derivative approximation - a
! later phase will DRIVE this function; this phase only exposes it.
!
! The injection law:
!
!      the trial replaces VALUES, never identity - the injected q
!      field is a copy of the unknown sample's own q field with
!      its real vector replaced, so the domain identity rides
!      unchanged;
!
!      the trial replaces q ONLY - every other sample, every
!      sample time, and every occupied non-q seat of every sample,
!      the unknown included, is preserved in the working copy;
!
!      the inputs are never mutated - the working sample set is a
!      deep copy, and the caller's samples stay exactly as given.
!
! Refused loudly before any copy: an empty sample set, an index
! outside the samples, an unknown sample without q, a trial with
! no real values, a trial whose shape disagrees with the unknown
! q field, and any q field of unsupported dynamic type.
!
! The problem carries nothing: no form, no motif, no samples, no
! graph, no solver, no map - two verbs over arguments.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_local_unknown_problems

  use iso_fortran_env      , only : dp => REAL64
  use class_graph_field    , only : field
  use gti_state_bundles    , only : gti_field_slot, GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_value_buffers    , only : gti_value_buffer
  use gti_form_interface   , only : gti_differentiable_form
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator

  implicit none

  private
  public :: gti_time_local_unknown_problem

  !===================================================================!
  ! The stateless verb-pair: inject a trial into one unknown seat,
  ! or go on and ask the residual its value there. The type keeps
  ! the public singular name; Fortran denies a type its host
  ! module's name, so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_local_unknown_problem

   contains

     procedure :: inject_trial_q
     procedure :: value

  end type gti_time_local_unknown_problem

contains

  !===================================================================!
  ! Build the working sample set with the trial q seated in the
  ! unknown sample. Not a solve: only the construction of the
  ! local state of the residual problem at one trial.
  !===================================================================!

  subroutine inject_trial_q(this, samples, unknown_index, trial_q, work_samples)

    class(gti_time_local_unknown_problem), intent(in)  :: this
    type(gti_time_sample)                , intent(in)  :: samples(:)
    integer                              , intent(in)  :: unknown_index
    type(gti_value_buffer)               , intent(in)  :: trial_q
    type(gti_time_sample), allocatable   , intent(out) :: work_samples(:)

    real(dp), allocatable :: trial_values(:), unknown_values(:)

    if (size(samples) < 1) then
       error stop 'gti_time_local_unknown_problem: at least one sample is required'
    end if

    if (unknown_index < 1 .or. unknown_index > size(samples)) then
       error stop 'gti_time_local_unknown_problem: unknown sample index is in range'
    end if

    if (.not. samples(unknown_index) % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_time_local_unknown_problem: unknown sample provides q'
    end if

    call trial_q % get_real(trial_values)
    if (size(trial_values) == 0) then
       error stop 'gti_time_local_unknown_problem: trial q has values'
    end if

    !----------------------------------------------------------------!
    ! The shape law: the trial agrees with the unknown q field in
    ! entries, components, and stored scalars.
    !----------------------------------------------------------------!

    associate(unknown_q => samples(unknown_index) % state % component(1 + GTI_STATE_Q) % value)

      call unknown_q % get_real_vector(unknown_values)

      if (trial_q % nentries /= unknown_q % num_entries() .or. &
           & trial_q % ncomp /= unknown_q % num_components() .or. &
           & size(trial_values) /= size(unknown_values)) then
         error stop 'gti_time_local_unknown_problem: trial q shape matches unknown q'
      end if

    end associate

    !----------------------------------------------------------------!
    ! Deep-copy the samples, then replace exactly one seat: the
    ! unknown's q, minted from its own field so the domain
    ! identity rides unchanged.
    !----------------------------------------------------------------!

    work_samples = samples

    call mint_trial_slot(samples(unknown_index) % state % component(1 + GTI_STATE_Q), &
         & trial_values, &
         & work_samples(unknown_index) % state % component(1 + GTI_STATE_Q))

  end subroutine inject_trial_q

  !===================================================================!
  ! residual(q_trial): inject the trial, then ask the form its
  ! value through the existing time-local seat - the motif builds
  ! U from the working samples, and the form evaluator holds the
  ! answer to its declared shape.
  !===================================================================!

  subroutine value(this, form, motif, samples, unknown_index, trial_q, &
       & design, time, output)

    class(gti_time_local_unknown_problem), intent(in)    :: this
    class(gti_differentiable_form)       , intent(in)    :: form
    type(gti_time_motif)                 , intent(in)    :: motif
    type(gti_time_sample)                , intent(in)    :: samples(:)
    integer                              , intent(in)    :: unknown_index
    type(gti_value_buffer)               , intent(in)    :: trial_q
    type(gti_design_bundle)              , intent(in)    :: design
    real(dp)                             , intent(in)    :: time
    type(gti_value_buffer)               , intent(inout) :: output

    type(gti_time_sample), allocatable      :: work_samples(:)
    type(gti_time_local_residual_evaluator) :: evaluator

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the injection replaces it wholesale
    allocate(work_samples(0))

    call this % inject_trial_q(samples, unknown_index, trial_q, work_samples)

    call evaluator % value(form, motif, work_samples, design, time, output)

  end subroutine value

  !===================================================================!
  ! Mint the injected q: copy the unknown sample's own q field -
  ! the copy carries the domain identity - and replace its values
  ! with the trial. One concrete field type is supported; any
  ! other dynamic type is refused, never guessed.
  !===================================================================!

  subroutine mint_trial_slot(prototype, values, slot)

    type(gti_field_slot), intent(in)    :: prototype
    real(dp)            , intent(in)    :: values(:)
    type(gti_field_slot), intent(inout) :: slot

    type(field) :: minted

    select type (proto => prototype % value)
    type is (field)
       minted = proto
       call minted % set_real_vector(values)
       slot % value = minted
    class default
       error stop 'gti_time_local_unknown_problem: q field dynamic type is supported'
    end select

  end subroutine mint_trial_slot

end module gti_time_local_unknown_problems
