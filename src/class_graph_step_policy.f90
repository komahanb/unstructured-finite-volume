!=====================================================================!
! The step policy: the walk's governance, and nothing else.
!
! An adaptive walk separates cleanly into mechanism and decision.
! The marcher owns the mechanism - trial steps, evidence measured
! as the distance from the extrapolating predictor, commitment of
! an accepted edge. This family owns every DECISION as three
! deferred verbs:
!
!      propose ····· the step to try first at a new edge
!      judge ······· accept or reject, from the evidence alone
!      retry ······· the step to try after a rejection
!
! A policy sees scalars - the estimate, the step, the attempt -
! and never a state, a graph, or a statement: policies decide,
! mechanisms execute, and a different adaptivity is a different
! policy, never a different marcher.
!
! The reference member halves: start at the caller's step, accept
! within the caller's tolerance, halve on every rejection. A
! policy example - not the policy.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_step_policy

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: step_policy
  public :: halving_policy

  type, abstract :: step_policy

   contains

     procedure(policy_propose), deferred :: propose
     procedure(policy_judge)  , deferred :: judge
     procedure(policy_retry)  , deferred :: retry

  end type step_policy

  abstract interface

     subroutine policy_propose(this, h)
       import :: step_policy, dp
       class(step_policy), intent(inout) :: this
       real(dp)          , intent(out)   :: h
     end subroutine policy_propose

     subroutine policy_judge(this, estimate, h, attempt, accept)
       import :: step_policy, dp
       class(step_policy), intent(inout) :: this
       real(dp)          , intent(in)    :: estimate
       real(dp)          , intent(in)    :: h
       integer           , intent(in)    :: attempt
       logical           , intent(out)   :: accept
     end subroutine policy_judge

     subroutine policy_retry(this, estimate, h)
       import :: step_policy, dp
       class(step_policy), intent(inout) :: this
       real(dp)          , intent(in)    :: estimate
       real(dp)          , intent(inout) :: h
     end subroutine policy_retry

  end interface

  type, extends(step_policy) :: halving_policy

     real(dp) :: first_step = 1.0_dp
     real(dp) :: tolerance  = 1.0e-2_dp

   contains

     procedure :: propose => halving_propose
     procedure :: judge   => halving_judge
     procedure :: retry   => halving_retry

  end type halving_policy

contains

  subroutine halving_propose(this, h)

    class(halving_policy), intent(inout) :: this
    real(dp)             , intent(out)   :: h

    if (this % first_step <= 0.0_dp) then
       error stop 'step_policy: the first step is positive'
    end if

    if (this % tolerance <= 0.0_dp) then
       error stop 'step_policy: the tolerance is positive'
    end if

    h = this % first_step

  end subroutine halving_propose

  subroutine halving_judge(this, estimate, h, attempt, accept)

    class(halving_policy), intent(inout) :: this
    real(dp)             , intent(in)    :: estimate
    real(dp)             , intent(in)    :: h
    integer              , intent(in)    :: attempt
    logical              , intent(out)   :: accept

    associate(unread_h => h, unread_attempt => attempt)
    end associate

    accept = estimate <= this % tolerance

  end subroutine halving_judge

  subroutine halving_retry(this, estimate, h)

    class(halving_policy), intent(inout) :: this
    real(dp)             , intent(in)    :: estimate
    real(dp)             , intent(inout) :: h

    associate(unread_this => this, unread_estimate => estimate)
    end associate

    h = 0.5_dp * h

  end subroutine halving_retry

end module class_graph_step_policy
