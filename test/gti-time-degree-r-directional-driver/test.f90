!=====================================================================!
! The degree-r directional traversal suite: over the solved
! nonlinear BDF1 chain
!
!      R = (q_u - q_h) + q_u^2 - xi,   q = 1 everywhere, xi = 1,
!
! ONE driver walks the graph at r = 1, 2, 3, 4 and computes every
! directional derivative along eta = 1 from the same J_u = 3 per
! relation, built once and eliminated r times:
!
!      r1:  q2^(1) =  1/3          (rhs      1)
!           q2^(2) = -2/27         (rhs  -2/9)
!           q2^(3) =  4/81         (rhs   4/27)
!           q2^(4) = -40/729       (rhs -40/243)
!      r2:  q3^(1) =  4/9          (rhs      4/3)
!           q3^(2) = -38/243       (rhs    -38/81)
!           q3^(3) =  340/2187     (rhs    340/729)
!           q3^(4) = -14848/59049  (rhs -14848/19683)
!
! Though the residual's own calculus stops at the quadratic, the
! third and fourth derivatives are alive: the chain-rule
! right-hand side transports products of lower-degree curvature,
! assembled by gti_chain_rule_assemblies at every degree. Degree
! 3 and degree 4 are verification cases of the SAME driver, not
! drivers of their own; at r = 2 the general loop must reproduce
! the dedicated degree-2 seat exactly. The graph's primal content
! survives every traversal untouched, and every number arrives
! through partial actions alone.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_degree_r_directional_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_degree_r_directional_drivers, only : &
       & gti_time_degree_r_directional_driver, &
       & gti_time_degree_r_directional_options, &
       & gti_time_degree_r_directional_result
  use gti_time_degree2_directional_drivers, only : &
       & gti_time_degree2_directional_driver, &
       & gti_time_degree2_directional_options, &
       & gti_time_degree2_directional_result
  use gti_toy_forms        , only : toy_qdot_square_form

  implicit none

  type(graph) :: states, designs
  type(field) :: one_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_degree_r_directional_driver)  :: driver
  type(gti_time_degree_r_directional_options) :: options
  type(gti_time_degree_r_directional_result)  :: result

  type(gti_time_degree2_directional_driver)  :: old_driver
  type(gti_time_degree2_directional_options) :: old_options
  type(gti_time_degree2_directional_result)  :: old_result

  type(toy_qdot_square_form) :: r_form
  type(gti_value_buffer)     :: eta

  real(dp), allocatable :: rv(:), ov(:)
  logical :: norms_ok, agrees
  integer :: v, nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti degree-r directional traversal suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call designs % declare()

  one_field = field('one', states, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([1.0_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp])

  !-------------------------------------------------------------------!
  ! The solved chain: q = 1 at every vertex IS the root of
  ! (q_u - q_h) + q_u^2 - xi at xi = 1.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(3))
  do v = 1, 3
     allocate(time_graph % vertex(v) % sample % state % component(1))
     time_graph % vertex(v) % sample % state % component(1) % value = one_field
     time_graph % vertex(v) % sample % time = real(v - 1, dp)
     time_graph % vertex(v) % has_solution = .true.
  end do

  allocate(time_graph % relation(2))
  call builder % bdf_uniform(1, 1.0_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 1.0_dp
  call builder % bdf_uniform(1, 1.0_dp, time_graph % relation(2) % motif)
  time_graph % relation(2) % sample_vertex   = [2, 3]
  time_graph % relation(2) % unknown_sample  = 2
  time_graph % relation(2) % evaluation_time = 2.0_dp

  norms_ok = .true.

  !-------------------------------------------------------------------!
  ! r = 1: the traversal is the tangent march.
  !-------------------------------------------------------------------!

  options % max_degree = 1
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call report(result % completed, &
       & "r=1: solve_all completes", nfail)

  call report(result % max_degree == 1, &
       & "r=1: the result records max_degree = 1", nfail)

  call report(result % completed_relations == 2, &
       & "r=1: both relations completed, vertices 2 then 3", nfail)

  call result % vertex_derivative(1, 2) % get_real(rv)
  call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "r=1: q2^(1) = 1/3", nfail)

  call result % vertex_derivative(1, 3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "r=1: q3^(1) = 4/9, v2's derivative fed r2 as history", nfail)

  call report(size(result % vertex_derivative, 1) == 1 .and. &
       & size(result % step(1) % rhs) == 1 .and. &
       & size(result % step(1) % derivative) == 1 .and. &
       & size(result % step(1) % linear_residual_norm) == 1, &
       & "r=1: no degree-2/3/4 storage exists in the result", nfail)

  !-------------------------------------------------------------------!
  ! r = 2: the general loop must BE the dedicated degree-2 seat.
  !-------------------------------------------------------------------!

  options % max_degree = 2
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(1, 2) % get_real(rv)
  call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "r=2: q2^(1) = 1/3", nfail)

  call result % vertex_derivative(2, 2) % get_real(rv)
  call report(matches(rv, [-2.0_dp / 27.0_dp], 1.0e-12_dp), &
       & "r=2: q2^(2) = -2/27", nfail)

  call result % vertex_derivative(1, 3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "r=2: q3^(1) = 4/9", nfail)

  call result % vertex_derivative(2, 3) % get_real(rv)
  call report(matches(rv, [-38.0_dp / 243.0_dp], 1.0e-12_dp), &
       & "r=2: q3^(2) = -38/243", nfail)

  call old_driver % solve_all(r_form, time_graph, design, eta, &
       & old_options, old_result)

  agrees = .true.
  call result % vertex_derivative(1, 2) % get_real(rv)
  call old_result % vertex_first(2) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(2, 2) % get_real(rv)
  call old_result % vertex_second(2) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(1, 3) % get_real(rv)
  call old_result % vertex_first(3) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(2, 3) % get_real(rv)
  call old_result % vertex_second(3) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(1) % rhs(1) % get_real(rv)
  call old_result % step(1) % first_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(1) % rhs(2) % get_real(rv)
  call old_result % step(1) % second_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(2) % rhs(1) % get_real(rv)
  call old_result % step(2) % first_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(2) % rhs(2) % get_real(rv)
  call old_result % step(2) % second_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)

  call report(agrees, &
       & "r=2: the general loop matches the dedicated degree-2 driver", nfail)

  !-------------------------------------------------------------------!
  ! r = 3: alive though the residual's third partials are zero.
  !-------------------------------------------------------------------!

  options % max_degree = 3
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(3, 2) % get_real(rv)
  call report(matches(rv, [4.0_dp / 81.0_dp], 1.0e-12_dp), &
       & "r=3: q2^(3) = 4/81, lower-degree curvature transported", nfail)

  call result % vertex_derivative(3, 3) % get_real(rv)
  call report(matches(rv, [340.0_dp / 2187.0_dp], 1.0e-12_dp), &
       & "r=3: q3^(3) = 340/2187", nfail)

  call result % step(1) % rhs(3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 27.0_dp], 1.0e-12_dp), &
       & "r=3: r1's degree-3 right-hand side is -B3 = [4/27]", nfail)

  call result % step(2) % rhs(3) % get_real(rv)
  call report(matches(rv, [340.0_dp / 729.0_dp], 1.0e-12_dp), &
       & "r=3: r2's degree-3 right-hand side is -B3 = [340/729]", nfail)

  !-------------------------------------------------------------------!
  ! r = 4: the top of the assembler's vocabulary.
  !-------------------------------------------------------------------!

  options % max_degree = 4
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(4, 2) % get_real(rv)
  call report(matches(rv, [-40.0_dp / 729.0_dp], 1.0e-12_dp), &
       & "r=4: q2^(4) = -40/729", nfail)

  call result % vertex_derivative(4, 3) % get_real(rv)
  call report(matches(rv, [-14848.0_dp / 59049.0_dp], 1.0e-12_dp), &
       & "r=4: q3^(4) = -14848/59049", nfail)

  call result % step(1) % rhs(4) % get_real(rv)
  call report(matches(rv, [-40.0_dp / 243.0_dp], 1.0e-12_dp), &
       & "r=4: r1's degree-4 right-hand side is -B4 = [-40/243]", nfail)

  call result % step(2) % rhs(4) % get_real(rv)
  call report(matches(rv, [-14848.0_dp / 19683.0_dp], 1.0e-12_dp), &
       & "r=4: r2's degree-4 right-hand side is -B4 = [-14848/19683]", nfail)

  !-------------------------------------------------------------------!
  ! Shared laws: exact eliminations, and a graph that never felt
  ! the traversals.
  !-------------------------------------------------------------------!

  call report(norms_ok, &
       & "every linear residual norm at every degree is below 1e-12", nfail)

  do v = 1, 3
     call time_graph % vertex(v) % sample % state % component(1 + GTI_STATE_Q) % &
          & value % get_real_vector(rv)
     if (.not. matches(rv, [1.0_dp], 1.0e-14_dp)) exit
  end do
  call report(v == 4, &
       & "the graph's q values survive all four traversals untouched", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "the graph's solution flags survive all four traversals untouched", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all degree-r traversal checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < tolerance)

  end function matches

  subroutine fold_norms(result, norms_ok)
    type(gti_time_degree_r_directional_result), intent(in)    :: result
    logical                                    , intent(inout) :: norms_ok
    integer :: r, s
    do r = 1, size(result % step)
       do s = 1, size(result % step(r) % linear_residual_norm)
          norms_ok = norms_ok .and. &
               & result % step(r) % linear_residual_norm(s) <= 1.0e-12_dp
       end do
    end do
  end subroutine fold_norms

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

end program test_gti_time_degree_r_directional_driver
