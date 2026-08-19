!=====================================================================!
! The refusals that must die at the controller seat, one per
! invocation:
!
!      attempts0    an attempt budget of zero
!      novertex     advancing a graph with no initial vertex
!      badstep      a policy proposing a zero step
!      badorder     a policy proposing order three
!      earlyorder2  a policy proposing order two against one
!                   history vertex
!      polinit      the halving policy with a zero initial step
!      poltol       the halving policy with a zero tolerance
!      polorder     the halving policy preferring order three
!
! A policy that contradicts itself dies at the deciding seat, not
! in a motif builder's message far from it. A spent budget is
! deliberately absent from this list - failure to advance is a
! lawful answer, proven positively in test.f90, never a refusal.
!
! Every case must error stop; a case that returns is a failure of
! the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_liar_policies

  use iso_fortran_env, only : dp => REAL64
  use gti_time_adaptive_controllers, only : gti_time_step_policy, &
       & gti_time_step_proposal, gti_time_attempt_evidence

  implicit none

  private
  public :: liar_policy

  !===================================================================!
  ! One configurable liar: it proposes whatever its knobs say,
  ! accepts everything, and never learns.
  !===================================================================!

  type, extends(gti_time_step_policy) :: liar_policy
     real(dp) :: proposes_step  = 0.5_dp
     integer  :: proposes_order = 1
   contains
     procedure :: first_proposal => liar_first_proposal
     procedure :: judge          => liar_judge
     procedure :: retry_proposal => liar_retry_proposal
  end type liar_policy

contains

  subroutine liar_first_proposal(this, history_count, proposal)
    class(liar_policy)          , intent(inout) :: this
    integer                     , intent(in)    :: history_count
    type(gti_time_step_proposal), intent(out)   :: proposal
    associate(unread => history_count)
    end associate
    proposal % step  = this % proposes_step
    proposal % order = this % proposes_order
  end subroutine liar_first_proposal

  subroutine liar_judge(this, evidence, accept)
    class(liar_policy)             , intent(inout) :: this
    type(gti_time_attempt_evidence), intent(in)    :: evidence
    logical                        , intent(out)   :: accept
    associate(unread => this)
    end associate
    associate(unread => evidence)
    end associate
    accept = .true.
  end subroutine liar_judge

  subroutine liar_retry_proposal(this, evidence, proposal)
    class(liar_policy)             , intent(inout) :: this
    type(gti_time_attempt_evidence), intent(in)    :: evidence
    type(gti_time_step_proposal)   , intent(out)   :: proposal
    associate(unread => evidence)
    end associate
    proposal % step  = this % proposes_step
    proposal % order = this % proposes_order
  end subroutine liar_retry_proposal

end module gti_liar_policies

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_options
  use gti_time_adaptive_controllers, only : &
       & gti_time_adaptive_controller, gti_time_adaptive_controller_options, &
       & gti_time_adaptive_advance_result, gti_time_halving_step_policy
  use gti_liar_policies    , only : liar_policy
  use gti_toy_forms        , only : toy_qdot_square_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, xi_field

  type(gti_time_graph)    :: time_graph, empty_graph
  type(gti_design_bundle) :: design

  type(gti_time_adaptive_controller)         :: controller
  type(gti_time_adaptive_controller_options) :: options
  type(gti_time_forward_options)             :: forward_options
  type(gti_time_halving_step_policy)         :: halving
  type(liar_policy)                          :: liar
  type(gti_time_adaptive_advance_result)     :: result

  type(toy_qdot_square_form) :: r_form
  character(len=32) :: which

  call get_command_argument(1, which)

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 1)
  call q0_field % set_real_vector([2.0_dp])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([1.0_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  allocate(time_graph % vertex(1))
  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q0_field
  time_graph % vertex(1) % sample % time = 0.0_dp
  time_graph % vertex(1) % has_solution = .true.

  select case (trim(which))

  case ('attempts0')

     options % max_attempts = 0
     call controller % advance(r_form, time_graph, design, halving, &
          & forward_options, options, result)

  case ('novertex')

     call controller % advance(r_form, empty_graph, design, halving, &
          & forward_options, options, result)

  case ('badstep')

     liar % proposes_step = 0.0_dp
     call controller % advance(r_form, time_graph, design, liar, &
          & forward_options, options, result)

  case ('badorder')

     liar % proposes_order = 3
     call controller % advance(r_form, time_graph, design, liar, &
          & forward_options, options, result)

  case ('earlyorder2')

     liar % proposes_order = 2
     call controller % advance(r_form, time_graph, design, liar, &
          & forward_options, options, result)

  case ('polinit')

     halving % initial_step = 0.0_dp
     call controller % advance(r_form, time_graph, design, halving, &
          & forward_options, options, result)

  case ('poltol')

     halving % error_tolerance = 0.0_dp
     call controller % advance(r_form, time_graph, design, halving, &
          & forward_options, options, result)

  case ('polorder')

     halving % preferred_order = 3
     call controller % advance(r_form, time_graph, design, halving, &
          & forward_options, options, result)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
