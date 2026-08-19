!=====================================================================!
! The adaptive growth suite: the time graph grows one candidate
! at a time, transactionally -
!
!      append -> solve -> accept keeps / reject rolls back -
!
! and the three structural laws hold: an accepted candidate
! becomes part of G_time, a rejected candidate leaves no trace,
! and a later accepted relation reads an earlier accepted vertex
! as history. Candidates carry their own motifs, so a DIRK step
! is accepted after a BDF step with the same verbs, and a grown
! graph is a graph like any other: the UNCHANGED solve-gradient
! driver runs on it and lands on the known closed forms
! value = 43/9, total gradient = 4/3. A non-converged candidate
! reports, rolls back whole, and writes nothing.
!
! No estimator, no step-size law: the accept decision is a
! logical from an external controller, and the transaction is
! the whole novelty.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_adaptive_growth_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_options
  use gti_time_adaptive_growth_drivers, only : gti_time_growth_candidate, &
       & gti_time_growth_step_result, gti_time_adaptive_growth_driver
  use gti_time_functional_seed_drivers, only : gti_time_functional
  use gti_time_solve_gradient_drivers, only : gti_time_solve_gradient_driver, &
       & gti_time_solve_gradient_options, gti_time_solve_gradient_result
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_sum_time_functional, toy_square_plus_form

  implicit none

  type(graph) :: states, scalars, designs
  type(field) :: q1_field, zero_field, one_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph, hetero_graph, flat_graph
  type(gti_design_bundle)      :: design, no_design

  type(gti_time_adaptive_growth_driver) :: grower
  type(gti_time_growth_candidate)       :: cand_a, cand_b, cand_c, cand_flat
  type(gti_time_growth_step_result)     :: step
  type(gti_time_forward_options)        :: forward_options, starved

  type(gti_time_functional)             :: terminal
  type(gti_time_solve_gradient_driver)  :: composer
  type(gti_time_solve_gradient_options) :: compose_options
  type(gti_time_solve_gradient_result)  :: composed

  type(toy_time_residual_form)  :: r_form
  type(toy_sum_time_functional) :: f_form
  type(toy_square_plus_form)    :: rootless

  type(gti_value_buffer) :: eta

  real(dp), allocatable :: rv(:)
  real(dp) :: root1(3), root2(3), root_dirk(3)
  integer :: av, ar, nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti adaptive graph growth suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call scalars % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  one_field = field('one', scalars, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  root1     = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root2     = (2.0_dp * root1 - 0.5_dp) / 3.0_dp
  root_dirk = (root1 - 0.5_dp) / 2.0_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! The base graph: one solved vertex, no relations. Candidates
  ! name their own vertex by the sentinel -1.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph)

  call make_candidate(cand_a, 0.5_dp, [1, -1], 0.5_dp)
  call builder % bdf_uniform(1, 0.5_dp, cand_a % relation % motif)

  call make_candidate(cand_b, 1.0_dp, [2, -1], 1.0_dp)
  call builder % bdf_uniform(1, 0.5_dp, cand_b % relation % motif)

  call make_candidate(cand_c, 2.0_dp, [2, -1], 2.0_dp)
  call builder % dirk_stage(0.5_dp, 2.0_dp, cand_c % relation % motif)

  !-------------------------------------------------------------------!
  ! Append and rollback, bare.
  !-------------------------------------------------------------------!

  call grower % append_candidate(time_graph, cand_a, av, ar)

  call report(time_graph % num_vertices() == 2 .and. av == 2, &
       & "append grows the vertex set from 1 to 2", nfail)

  call report(time_graph % num_relations() == 1 .and. ar == 1, &
       & "append grows the relation set from 0 to 1", nfail)

  call report(all(time_graph % relation(1) % sample_vertex == [1, 2]), &
       & "the sentinel -1 resolved to the appended vertex 2", nfail)

  call report(time_graph % relation(1) % unknown_vertex() == 2, &
       & "the appended relation's unknown is the appended vertex", nfail)

  call grower % discard_last_candidate(time_graph, av, ar)

  call report(time_graph % num_vertices() == 1 .and. &
       & time_graph % num_relations() == 0, &
       & "rollback restores one vertex and no relations", nfail)

  call time_graph % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp) .and. &
       & time_graph % vertex(1) % has_solution, &
       & "rollback preserves the surviving vertex's q and flag", nfail)

  !-------------------------------------------------------------------!
  ! Accepted: the candidate becomes part of the graph.
  !-------------------------------------------------------------------!

  call grower % try_candidate(r_form, time_graph, cand_a, design, &
       & forward_options, .true., step)

  call report(step % accepted .and. step % solved .and. step % converged, &
       & "an accepted converged candidate reports all three flags", nfail)

  call report(step % num_vertices_after == 2 .and. &
       & step % num_relations_after == 1, &
       & "the graph keeps the accepted vertex and relation", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "the accepted vertex carries the solved q = (2 q1 - xi)/3", nfail)

  call report(time_graph % vertex(2) % has_solution, &
       & "the accepted vertex is marked solved", nfail)

  !-------------------------------------------------------------------!
  ! Rejected after convergence: no trace remains.
  !-------------------------------------------------------------------!

  call grower % try_candidate(r_form, time_graph, cand_b, design, &
       & forward_options, .false., step)

  call report((.not. step % accepted) .and. step % solved .and. &
       & step % converged, &
       & "a rejected converged candidate solved, converged, and was refused", nfail)

  call report(step % num_vertices_after == 2 .and. &
       & step % num_relations_after == 1 .and. &
       & time_graph % num_vertices() == 2, &
       & "the rejected candidate left no vertex and no relation behind", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "the earlier accepted vertex keeps its value through the rejection", nfail)

  !-------------------------------------------------------------------!
  ! Accepted later: the same candidate, kept this time, reads the
  ! earlier accepted vertex as history.
  !-------------------------------------------------------------------!

  call grower % try_candidate(r_form, time_graph, cand_b, design, &
       & forward_options, .true., step)

  call report(step % accepted .and. time_graph % num_vertices() == 3 .and. &
       & time_graph % num_relations() == 2, &
       & "the re-proposed candidate is accepted: three vertices, two relations", nfail)

  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-9_dp), &
       & "the new vertex solved FROM the earlier accepted vertex as history", nfail)

  !-------------------------------------------------------------------!
  ! Heterogeneous growth: a DIRK candidate after a BDF one.
  !-------------------------------------------------------------------!

  call seed_base_graph(hetero_graph)

  call grower % try_candidate(r_form, hetero_graph, cand_a, design, &
       & forward_options, .true., step)
  call grower % try_candidate(r_form, hetero_graph, cand_c, design, &
       & forward_options, .true., step)

  call report(step % accepted .and. hetero_graph % num_vertices() == 3 .and. &
       & hetero_graph % num_relations() == 2, &
       & "a DIRK candidate is accepted after a BDF one", nfail)

  call hetero_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root_dirk, 1.0e-9_dp), &
       & "the DIRK vertex solved to (v2 - xi)/2 = [0, 1/3, 1]", nfail)

  !-------------------------------------------------------------------!
  ! A grown graph is a graph: the unchanged solve-gradient driver
  ! runs on it and lands on the known closed forms.
  !-------------------------------------------------------------------!

  allocate(terminal % term(1))
  terminal % term(1) % vertex_index    = 3
  terminal % term(1) % evaluation_time = 1.0_dp

  call composer % solve(r_form, f_form, time_graph, terminal, design, eta, &
       & compose_options, composed)

  call report(composed % completed .and. &
       & abs(composed % value - 43.0_dp / 9.0_dp) < 1.0e-12_dp, &
       & "the grown graph carries the full pipeline: value = 43/9", nfail)

  call composed % total_design_gradient_action % get_real(rv)
  call report(matches(rv, [4.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the grown graph's total design gradient is 4/3", nfail)

  !-------------------------------------------------------------------!
  ! A non-converged candidate rolls back whole and writes nothing.
  !-------------------------------------------------------------------!

  allocate(flat_graph % vertex(1))
  allocate(flat_graph % vertex(1) % sample % state % component(1))
  flat_graph % vertex(1) % sample % state % component(1) % value = one_field
  flat_graph % vertex(1) % has_solution = .true.

  cand_flat % vertex % sample % time = 1.0_dp
  allocate(cand_flat % vertex % sample % state % component(1))
  cand_flat % vertex % sample % state % component(1) % value = one_field
  allocate(cand_flat % relation % motif % rule(1))
  cand_flat % relation % motif % rule(1) % state_component = GTI_STATE_Q
  cand_flat % relation % motif % rule(1) % weight = [1.0_dp]
  cand_flat % relation % sample_vertex  = [-1]
  cand_flat % relation % unknown_sample = 1

  starved % newton % max_iterations = 1

  call grower % try_candidate(rootless, flat_graph, cand_flat, no_design, &
       & starved, .true., step)

  call report(step % solved .and. (.not. step % converged) .and. &
       & (.not. step % accepted), &
       & "a rootless candidate reports solved, unconverged, rejected", nfail)

  call report(step % num_vertices_after == 1 .and. &
       & step % num_relations_after == 0 .and. &
       & flat_graph % num_vertices() == 1, &
       & "the non-converged candidate rolled back whole: nothing written", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all adaptive growth checks passed"
  else
     error stop
  end if

contains

  subroutine seed_base_graph(g)

    type(gti_time_graph), intent(out) :: g

    allocate(g % vertex(1))
    allocate(g % vertex(1) % sample % state % component(1))
    g % vertex(1) % sample % state % component(1) % value = q1_field
    g % vertex(1) % sample % time = 0.0_dp
    g % vertex(1) % has_solution = .true.

  end subroutine seed_base_graph

  subroutine make_candidate(candidate, vertex_time, tuple, evaluation_time)

    type(gti_time_growth_candidate), intent(out) :: candidate
    real(dp)                       , intent(in)  :: vertex_time
    integer                        , intent(in)  :: tuple(:)
    real(dp)                       , intent(in)  :: evaluation_time

    candidate % vertex % sample % time = vertex_time
    allocate(candidate % vertex % sample % state % component(1))
    candidate % vertex % sample % state % component(1) % value = zero_field

    candidate % relation % sample_vertex   = tuple
    candidate % relation % unknown_sample  = 2
    candidate % relation % evaluation_time = evaluation_time

  end subroutine make_candidate

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < tolerance)

  end function matches

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_time_adaptive_growth_driver
