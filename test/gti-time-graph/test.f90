!=====================================================================!
! The time graph suite: G_time = (S_time, R_time) over three
! sample vertices and two BDF1 relations sharing vertex 2 -
!
!      v1 (t=0.0) --r1--> v2 (t=0.5) --r2--> v3 (t=1.0)
!
! r1 = ([v1, v2], unknown 2), r2 = ([v2, v3], unknown 2). The
! shared vertex is the whole point: writing r1's solution into
! vertex 2 is what r2's next materialization sees - that shared
! seeing IS the time coupling, and no march driver exists yet.
!
! The suite proves counting, local-to-global unknown translation,
! validation, relation-local materialization as true copies,
! solution write-back that replaces values but never identity,
! preservation of times and non-q seats, and motif rows surviving
! storage.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_graph

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT
  use gti_time_local_schemes, only : gti_time_sample
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph

  implicit none

  type(graph) :: states
  type(field) :: q1_field, z2_field, z3_field, qd7_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph

  type(gti_time_sample), allocatable :: samples(:)
  type(gti_value_buffer) :: q_solved_1, q_solved_2

  real(dp), allocatable :: rv(:)
  real(dp) :: root1(3), root2(3)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time graph representation suite"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! Three vertices, two BDF1 relations sharing vertex 2. Vertex 2
  ! also carries an occupied qdot seat, the preservation witness.
  !-------------------------------------------------------------------!

  call states % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  z2_field = field('placeholder v2', states, 3)
  call z2_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  z3_field = field('placeholder v3', states, 3)
  call z3_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  qd7_field = field('qdot7', states, 3)
  call qd7_field % set_real_vector([7.0_dp, 7.0_dp, 7.0_dp])

  allocate(time_graph % vertex(3))

  allocate(time_graph % vertex(1) % sample % state % component(1))
  time_graph % vertex(1) % sample % state % component(1) % value = q1_field
  time_graph % vertex(1) % sample % time = 0.0_dp
  time_graph % vertex(1) % has_solution = .true.

  allocate(time_graph % vertex(2) % sample % state % component(2))
  time_graph % vertex(2) % sample % state % component(1) % value = z2_field
  time_graph % vertex(2) % sample % state % component(2) % value = qd7_field
  time_graph % vertex(2) % sample % time = 0.5_dp

  allocate(time_graph % vertex(3) % sample % state % component(1))
  time_graph % vertex(3) % sample % state % component(1) % value = z3_field
  time_graph % vertex(3) % sample % time = 1.0_dp

  allocate(time_graph % relation(2))

  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp

  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(2) % motif)
  time_graph % relation(2) % sample_vertex   = [2, 3]
  time_graph % relation(2) % unknown_sample  = 2
  time_graph % relation(2) % evaluation_time = 1.0_dp

  root1 = [0.5_dp, 1.1666666666666667_dp, 2.5_dp]
  root2 = [0.25_dp, 0.6111111111111112_dp, 1.5_dp]
  call q_solved_1 % set_real(root1)
  call q_solved_2 % set_real(root2)

  !-------------------------------------------------------------------!
  ! Counting, translation, validation, and stored rows.
  !-------------------------------------------------------------------!

  call report(time_graph % num_vertices() == 3, &
       & "the graph counts three sample vertices", nfail)

  call report(time_graph % num_relations() == 2, &
       & "the graph counts two motif relations", nfail)

  call report(time_graph % relation(1) % arity() == 2 .and. &
       & time_graph % relation(2) % arity() == 2, &
       & "both relations touch two vertices", nfail)

  call report(time_graph % relation(1) % unknown_vertex() == 2, &
       & "r1's local unknown translates to global vertex 2", nfail)

  call report(time_graph % relation(2) % unknown_vertex() == 3, &
       & "r2's local unknown translates to global vertex 3", nfail)

  call time_graph % validate()
  call report(.true., &
       & "validate accepts the two-relation BDF1 path", nfail)

  call report(matches(time_graph % relation(1) % motif % rule(2) % weight, &
       & [-2.0_dp, 2.0_dp], 1.0e-14_dp), &
       & "motif rows survive relation storage: BDF1 qdot row [-2, 2]", nfail)

  !-------------------------------------------------------------------!
  ! Materialization: relation-local order, preserved times, and
  ! true copies.
  !-------------------------------------------------------------------!

  call time_graph % build_samples(1, samples)
  call samples(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "r1 materializes vertex 1's q in local seat 1", nfail)

  call report(abs(samples(1) % time - 0.0_dp) < 1.0e-14_dp .and. &
       & abs(samples(2) % time - 0.5_dp) < 1.0e-14_dp, &
       & "r1 materializes the sample times [0.0, 0.5]", nfail)

  call samples(2) % state % component(1 + GTI_STATE_Q) % value % set_real_vector( &
       & [9.0_dp, 9.0_dp, 9.0_dp])
  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp), &
       & "materialized samples are copies: mutating them leaves the graph alone", nfail)

  call time_graph % build_samples(2, samples)
  call samples(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-14_dp) .and. &
       & abs(samples(1) % time - 0.5_dp) < 1.0e-14_dp .and. &
       & abs(samples(2) % time - 1.0_dp) < 1.0e-14_dp, &
       & "r2 materializes vertices [2, 3] before any write", nfail)

  !-------------------------------------------------------------------!
  ! Write-back: r1's solution lands on vertex 2, r2 sees it.
  !-------------------------------------------------------------------!

  call time_graph % write_unknown_q(1, q_solved_1)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-14_dp), &
       & "write_unknown_q(r1) replaces vertex 2's q", nfail)

  call time_graph % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "write_unknown_q(r1) preserves vertex 1's q", nfail)

  call report(time_graph % vertex(2) % has_solution, &
       & "write_unknown_q(r1) marks vertex 2 solved", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_QDOT) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [7.0_dp, 7.0_dp, 7.0_dp], 1.0e-14_dp) .and. &
       & abs(time_graph % vertex(2) % sample % time - 0.5_dp) < 1.0e-14_dp, &
       & "the write preserves vertex 2's qdot seat and its time", nfail)

  call time_graph % build_samples(2, samples)
  call samples(1) % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-14_dp), &
       & "after writing r1, r2's materialization sees the new vertex 2", nfail)

  call time_graph % write_unknown_q(2, q_solved_2)
  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-14_dp), &
       & "write_unknown_q(r2) replaces vertex 3's q", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-14_dp), &
       & "write_unknown_q(r2) leaves vertex 2 as r1 wrote it", nfail)

  call time_graph % validate()
  call report(.true., &
       & "validate still accepts the graph after both writes", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all time graph checks passed"
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

end program test_gti_time_graph
