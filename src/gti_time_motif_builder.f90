!=====================================================================!
! GTI NAMED TIME MOTIF BUILDERS
!
! The seam this module respects: builders CHOOSE coefficient rows,
! the time-local seat APPLIES them, and nothing solves anything.
! A named scheme formula enters here and leaves as a plain
! gti_time_motif - two rows, one weight per sample - that the
! existing gti_time_local_residual_evaluator consumes unchanged.
!
! Three named families, each over samples ordered oldest -> newest:
!
!   BDF, uniform step, orders 1 and 2
!
!      order 1, samples [q_{n-1}, q_n]:
!         q    = [ 0, 1 ]
!         qdot = [ -1/h, 1/h ]
!
!      order 2, samples [q_{n-2}, q_{n-1}, q_n]:
!         q    = [ 0, 0, 1 ]
!         qdot = [ 1/(2h), -2/h, 3/(2h) ]
!
!   DIRK, one diagonal stage, samples [q_base, q_stage] where
!   q_base = q_n + h sum_{j<i} a_ij k_j is assembled OUTSIDE:
!
!         q    = [ 0, 1 ]
!         qdot = [ -1/(h gamma), 1/(h gamma) ]
!
!   ABM, implicit corrector, samples [q_base, q_new] where q_base
!   carries the Adams predictor/history terms, assembled OUTSIDE:
!
!         q    = [ 0, 1 ]
!         qdot = [ -1/(h beta0), 1/(h beta0) ],
!         beta0 = 1 at order 1, 1/2 at order 2
!
! What a builder does not do: no Butcher table, no stage storage,
! no Adams predictor, no history integrals, no stage or corrector
! solve, no time graph, no adaptivity. Unsupported orders,
! nonpositive steps, and a zero DIRK diagonal are refused loudly
! before any row is minted.
!
! The builder carries nothing: no map, no graph, no solver, no
! scheme table, no form, no sample, no state - a name and a step
! size in, coefficient rows out.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_motif_builders

  use iso_fortran_env      , only : dp => REAL64
  use gti_form_interface   , only : GTI_STATE_Q, GTI_STATE_QDOT
  use gti_time_local_schemes, only : gti_time_motif, gti_time_component_rule

  implicit none

  private
  public :: gti_time_motif_builder

  !===================================================================!
  ! The stateless builder verb-set. The type keeps the public
  ! singular name; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_time_motif_builder

   contains

     procedure :: bdf_uniform
     procedure :: dirk_stage
     procedure :: abm_corrector

  end type gti_time_motif_builder

contains

  !===================================================================!
  ! BDF over a uniform step: the newest sample is the state, and
  ! the backward difference of the ordered samples is its rate.
  !===================================================================!

  pure subroutine bdf_uniform(this, order, h, motif)

    class(gti_time_motif_builder), intent(in)  :: this
    integer                      , intent(in)  :: order
    real(dp)                     , intent(in)  :: h
    type(gti_time_motif)         , intent(out) :: motif

    if (order /= 1 .and. order /= 2) then
       error stop 'gti_time_motif_builder: BDF order is supported'
    end if

    if (h <= 0.0_dp) then
       error stop 'gti_time_motif_builder: step size is positive'
    end if

    select case (order)
    case (1)
       call emit_rows([0.0_dp, 1.0_dp], &
            &         [-1.0_dp / h, 1.0_dp / h], motif)
    case default
       call emit_rows([0.0_dp, 0.0_dp, 1.0_dp], &
            &         [1.0_dp / (2.0_dp * h), -2.0_dp / h, &
            &          3.0_dp / (2.0_dp * h)], motif)
    end select

  end subroutine bdf_uniform

  !===================================================================!
  ! One DIRK diagonal stage: the stage sample is the state, and
  ! its distance from the externally assembled base, over h gamma,
  ! is its rate. The Butcher table lives elsewhere; only the
  ! diagonal reaches this seat.
  !===================================================================!

  pure subroutine dirk_stage(this, gamma, h, motif)

    class(gti_time_motif_builder), intent(in)  :: this
    real(dp)                     , intent(in)  :: gamma
    real(dp)                     , intent(in)  :: h
    type(gti_time_motif)         , intent(out) :: motif

    if (h <= 0.0_dp) then
       error stop 'gti_time_motif_builder: step size is positive'
    end if

    if (gamma == 0.0_dp) then
       error stop 'gti_time_motif_builder: DIRK diagonal is nonzero'
    end if

    call emit_rows([0.0_dp, 1.0_dp], &
         &         [-1.0_dp / (h * gamma), 1.0_dp / (h * gamma)], motif)

  end subroutine dirk_stage

  !===================================================================!
  ! The ABM implicit corrector: the new sample is the state, and
  ! its distance from the externally assembled predictor/history
  ! base, over h beta0, is its rate. beta0 is 1 at order 1 and
  ! 1/2 at order 2; the predictor itself never enters this seat.
  !===================================================================!

  pure subroutine abm_corrector(this, order, h, motif)

    class(gti_time_motif_builder), intent(in)  :: this
    integer                      , intent(in)  :: order
    real(dp)                     , intent(in)  :: h
    type(gti_time_motif)         , intent(out) :: motif

    real(dp) :: beta0

    if (order /= 1 .and. order /= 2) then
       error stop 'gti_time_motif_builder: ABM order is supported'
    end if

    if (h <= 0.0_dp) then
       error stop 'gti_time_motif_builder: step size is positive'
    end if

    if (order == 1) then
       beta0 = 1.0_dp
    else
       beta0 = 0.5_dp
    end if

    call emit_rows([0.0_dp, 1.0_dp], &
         &         [-1.0_dp / (h * beta0), 1.0_dp / (h * beta0)], motif)

  end subroutine abm_corrector

  !===================================================================!
  ! Mint the two-row motif every builder in this phase emits: one
  ! q row, one qdot row.
  !===================================================================!

  pure subroutine emit_rows(q_weight, qdot_weight, motif)

    real(dp)            , intent(in)  :: q_weight(:), qdot_weight(:)
    type(gti_time_motif), intent(out) :: motif

    allocate(motif % rule(2))

    motif % rule(1) % state_component = GTI_STATE_Q
    motif % rule(1) % weight          = q_weight

    motif % rule(2) % state_component = GTI_STATE_QDOT
    motif % rule(2) % weight          = qdot_weight

  end subroutine emit_rows

end module gti_time_motif_builders
