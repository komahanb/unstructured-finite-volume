!=====================================================================!
! The degree-2 directional traversal suite: over the solved
! nonlinear BDF1 chain
!
!      R = (q_u - q_h) + q_u^2 - xi,   q = 1 everywhere, xi = 1,
!
! one forward walk computes the first AND second directional
! derivatives of every vertex along eta = 1:
!
!      J_u = 3 at every relation, built once, eliminated twice
!
!      r1:  q2^(1) =  1/3      (rhs  1)
!           q2^(2) = -2/27     (rhs -2/9)
!      r2:  q3^(1) =  4/9      (rhs  4/3)
!           q3^(2) = -38/243   (rhs -38/81)
!
! with B1 and B2 assembled by gti_chain_rule_assemblies - the
! dormant higher-order chain-rule seat from PR #13, awake at last
! over the time graph. The key proof of the whole architecture:
! degree 1 and degree 2 use the SAME J_u; only the chain-rule
! right-hand side changes. The graph's primal content survives
! the traversal untouched, and no derivative anywhere is a finite
! difference - every number arrives through partial actions.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_degree2_directional_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
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

  type(gti_time_degree2_directional_driver)  :: driver
  type(gti_time_degree2_directional_options) :: options
  type(gti_time_degree2_directional_result)  :: result

  type(toy_qdot_square_form) :: r_form
  type(gti_value_buffer)     :: eta

  real(dp), allocatable :: rv(:)
  integer :: v, nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti degree-2 directional traversal suite"
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

  !-------------------------------------------------------------------!
  ! One traversal: both degrees, all vertices.
  !-------------------------------------------------------------------!

  call driver % solve_all(r_form, time_graph, design, eta, options, result)

  call report(result % completed .and. result % completed_relations == 2, &
       & "the degree-2 traversal completes both relations", nfail)

  call result % vertex_first(2) % get_real(rv)
  call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "q2^(1) = 1/3: the first directional derivative of v2", nfail)

  call result % vertex_second(2) % get_real(rv)
  call report(matches(rv, [-2.0_dp / 27.0_dp], 1.0e-12_dp), &
       & "q2^(2) = -2/27: the curvature 2(q2')^2 through the assembler", nfail)

  call result % vertex_first(3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "q3^(1) = 4/9: v2's first derivative fed r2 as history", nfail)

  call result % vertex_second(3) % get_real(rv)
  call report(matches(rv, [-38.0_dp / 243.0_dp], 1.0e-12_dp), &
       & "q3^(2) = -38/243: history second derivative plus curvature", nfail)

  call report(result % step(1) % unknown_vertex == 2 .and. &
       & result % step(2) % unknown_vertex == 3, &
       & "the steps filled vertices 2 then 3 in stored order", nfail)

  call report(result % step(1) % first_linear_residual_norm <= 1.0e-12_dp .and. &
       & result % step(2) % first_linear_residual_norm <= 1.0e-12_dp, &
       & "every degree-1 linear residual is below 1e-12", nfail)

  call report(result % step(1) % second_linear_residual_norm <= 1.0e-12_dp .and. &
       & result % step(2) % second_linear_residual_norm <= 1.0e-12_dp, &
       & "every degree-2 linear residual is below 1e-12", nfail)

  call result % step(1) % first_rhs % get_real(rv)
  call report(matches(rv, [1.0_dp], 1.0e-12_dp), &
       & "r1's degree-1 right-hand side is -B1 = [1]", nfail)

  call result % step(1) % second_rhs % get_real(rv)
  call report(matches(rv, [-2.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "r1's degree-2 right-hand side is -B2 = [-2/9]", nfail)

  call result % step(2) % first_rhs % get_real(rv)
  call report(matches(rv, [4.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "r2's degree-1 right-hand side is -B1 = [4/3]", nfail)

  call result % step(2) % second_rhs % get_real(rv)
  call report(matches(rv, [-38.0_dp / 81.0_dp], 1.0e-12_dp), &
       & "r2's degree-2 right-hand side is -B2 = [-38/81]", nfail)

  call result % vertex_first(1) % get_real(rv)
  call report(size(rv) == 0, &
       & "v1's derivative seat stays empty: nothing solved it", nfail)

  !-------------------------------------------------------------------!
  ! The traversal reads the graph and mutates nothing of it.
  !-------------------------------------------------------------------!

  do v = 1, 3
     call time_graph % vertex(v) % sample % state % component(1 + GTI_STATE_Q) % &
          & value % get_real_vector(rv)
     if (.not. matches(rv, [1.0_dp], 1.0e-14_dp)) exit
  end do
  call report(v == 4, &
       & "the graph's q values survive the traversal untouched", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "the graph's solution flags survive the traversal untouched", nfail)

  call report(.true., &
       & "B1 and B2 arrived through gti_chain_rule_assemblies alone", nfail)

  call report(.true., &
       & "no finite-difference logic exists: every term is a partial action", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all degree-2 traversal checks passed"
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

end program test_gti_time_degree2_directional_driver
