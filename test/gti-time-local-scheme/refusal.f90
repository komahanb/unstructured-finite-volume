!=====================================================================!
! The refusals that must die at the time-local seat, one per
! invocation:
!
!      nosamples    a motif applied to no samples at all
!      unknowncomp  a rule building a component outside the three
!                   GTI_STATE_* orders
!      dupq         two rules both building q
!      noweights    a rule with no weights allocated
!      wrongcount   a rule with three weights over two samples
!      missingq     a sample whose bundle holds no q
!      qshape       two samples whose q vectors disagree in shape
!
! Every case must error stop before any point is built; a case
! that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT
  use gti_design_bundles   , only : gti_design_bundle
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator

  implicit none

  type(graph) :: states, short
  type(field) :: qa_field, qb_field, qc_field

  type(gti_time_sample)   :: samples(2)
  type(gti_time_sample)   :: none(0)
  type(gti_time_motif)    :: motif
  type(gti_design_bundle) :: design

  type(gti_time_local_residual_evaluator) :: evaluator
  type(gti_evaluation_point)              :: point
  character(len=32) :: which

  call get_command_argument(1, which)

  call states % declare()
  call short  % declare()

  qa_field = field('qa', states, 3)
  call qa_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  qb_field = field('qb', states, 3)
  call qb_field % set_real_vector([3.0_dp, 5.0_dp, 8.0_dp])
  qc_field = field('qc', short, 2)
  call qc_field % set_real_vector([1.0_dp, 2.0_dp])

  allocate(samples(1) % state % component(1))
  samples(1) % state % component(1) % value = qa_field
  allocate(samples(2) % state % component(1))
  samples(2) % state % component(1) % value = qb_field

  allocate(motif % rule(2))
  motif % rule(1) % state_component = GTI_STATE_Q
  motif % rule(1) % weight = [0.0_dp, 1.0_dp]
  motif % rule(2) % state_component = GTI_STATE_QDOT
  motif % rule(2) % weight = [-1.0_dp, 1.0_dp]

  select case (trim(which))

  case ('nosamples')

     call evaluator % build_point(motif, none, design, 0.0_dp, point)

  case ('unknowncomp')

     motif % rule(2) % state_component = 999
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case ('dupq')

     motif % rule(2) % state_component = GTI_STATE_Q
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case ('noweights')

     deallocate(motif % rule(1) % weight)
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case ('wrongcount')

     motif % rule(1) % weight = [1.0_dp, 0.0_dp, 0.0_dp]
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case ('missingq')

     deallocate(samples(2) % state % component)
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case ('qshape')

     samples(2) % state % component(1) % value = qc_field
     call evaluator % build_point(motif, samples, design, 0.0_dp, point)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
