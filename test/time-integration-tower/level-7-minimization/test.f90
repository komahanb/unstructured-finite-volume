!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 7 . MINIMIZATION
!
! The level answers one question:
!
!      CAN THE ALREADY DOMAIN-EXPLICIT MINIMIZER SOLVE THE TEMPORAL
!      STEP ON Q WHILE THE GRAPH HOST HAS FIVE UNRELATED VERTICES?
!
! It can, and the interesting part is how little had to happen for
! that to be true.
!
!                    THE MINIMIZER WAS ALREADY RIGHT
!
! graph_minimization takes its unknown domain as an EXPLICIT
! argument - `attach(action, on, unknown_domain, ncomp)` - and asks
! the ACTION for the residual domain rather than the host. Its own
! comment says so: "no hidden fallback to the host's vertices; a
! caller that means vertices says so at its own call site."
!
! So once Level 6's correction made the step operator preserve its
! action's domain, this level needed no production change at all.
! That contrast is the evidence: the same seam that was open in the
! discretization was already closed in the minimization, by a
! contract written the other way round. Nothing here "improves"
! graph_minimization, because nothing here found anything to
! improve.
!
!                    THE TWO SOLVES
!
! The solver's matvec is the attached operation minus its affine
! part, so an implicit step becomes a linear system with the affine
! constant moved to the right-hand side:
!
!   backward euler   c = R(0) = -q0 = [-2, 0]
!                    (I + hM) q1 = q0,   M = [[1,0],[-1,1]]
!                    rhs = [2, 0]        ->  q1 = [4/3, 4/9]
!
!   bdf-2            c = R(0) = -2 q1 + (1/2) q0
!                    rhs = -c = [5/3, 8/9]
!                                        ->  q2 = [5/6, 47/72]
!
! Both right-hand sides are FIELDS ON Q, both solutions come back
! on Q, and the host has five vertices throughout.
!
! NO MARCHER APPEARS HERE. The marcher stamps a step along a chain;
! this level solves one step, twice, by hand. Whether the marcher
! can carry a state domain that is not its host's vertex set is a
! question for a level that has not been built, and the import gate
! refuses the module until then.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_7

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, verdict
  use time_assert           , only : NQ, NT, TOL
  use time_assert           , only : H_STEP, Q0, Q_BE1, Q_BDF2
  use graph_fractal        , only : set_graph => graph
  use map_set        , only : set_map
  use graph_directed_view   , only : directed_graph
  use graph_field_calculus  , only : graph_field
  use class_graph           , only : directed_stored_graph
  use class_graph_field     , only : field
  use class_graph_step      , only : step_operator, backward_euler, bdf
  use class_graph_gmres     , only : gmres
  use time_carriers_fixture , only : time_carriers
  use triangular_decay_fixture, only : triangular_decay

  implicit none

  type(set_graph)      :: q, t, e
  type(set_map)      :: sets
  type(directed_stored_graph)     :: ht
  type(triangular_decay) :: decay
  integer                :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 7 . solve"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)

  ! The same compatibility host as Level 6: five vertices, and not Q.
  ht = directed_stored_graph(NT, tails=[1,2,3,4], heads=[2,3,4,5])

  decay = triangular_decay(q, NQ)

  call check_the_host_is_still_not_q(nfail)
  call check_backward_euler_solve(nfail)
  call check_bdf2_solve(nfail)
  call check_unknown_domain_is_the_caller_s_word(nfail)

  call verdict(nfail, "level 7")

contains

  !===================================================================!
  ! Restated at this level rather than inherited, because every
  ! assertion below is only interesting while it holds.
  !===================================================================!

  subroutine check_the_host_is_still_not_q(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: hv

    hv = ht % vertex_set()

    call report(ht % num_vertices() .eq. NT .and. sets % size_of(q) .eq. NQ &
         &      .and. .not. hv % same_as(q), &
         & "the host still has five vertices and Q still has two: " // &
         & "the solve below is not a coincidence of sizes", nfail)

  end subroutine check_the_host_is_still_not_q

  !===================================================================!
  ! THE first complete implicit temporal step.
  !===================================================================!

  subroutine check_backward_euler_solve(nfail)

    integer, intent(inout) :: nfail

    type(step_operator)             :: step
    type(gmres)                     :: solver
    type(field)                     :: rhs
    class(graph_field), allocatable :: answer
    type(set_graph)  :: d
    integer         :: n_d
    real(dp), allocatable           :: c(:), v(:)

    step = backward_euler(decay, H_STEP)
    step % qold = Q0

    call solver % attach(step, ht, q, NQ, ncomp=1)

    call solver % domain(ht, d, n_d)
    call report(d % same_as(q), &
         & "the SOLVER answers Q when asked its domain, on a " // &
         & "five-vertex host", nfail)

    ! The affine part: R(0) = -q0, because S(0) = 0.
    call solver % constant(c)
    call report(size(c) .eq. NQ .and. &
         &      maxval(abs(c + Q0)) .lt. TOL, &
         & "the backward-euler affine constant is -q0 = [-2, 0], " // &
         & "measured rather than assumed", nfail)

    ! ... so the linear statement is (I + hM) q1 = q0.
    rhs = field('rhs', q, NQ, ncomp=1)
    call rhs % set_real_vector(Q0)

    call solver % apply(ht, [rhs], answer)

    d = answer % domain()
    call report(d % same_as(q), &
         & "and the SOLUTION lands on Q, through the operation " // &
         & "face, with the host present and unread", nfail)

    call answer % get_real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v - Q_BE1)) .lt. TOL, &
         & "q1 = [4/3, 4/9] - THE FIRST COMPLETE IMPLICIT TEMPORAL " // &
         & "STEP, solved on the state domain", nfail)

  end subroutine check_backward_euler_solve

  !===================================================================!
  ! The two-step scheme, solved the same way - a second instance,
  ! attached to a different statement on the same domain.
  !===================================================================!

  subroutine check_bdf2_solve(nfail)

    integer, intent(inout) :: nfail

    type(step_operator)             :: step
    type(gmres)                     :: solver
    type(field)                     :: rhs
    class(graph_field), allocatable :: answer
    type(set_graph)  :: d
    real(dp), allocatable           :: c(:), v(:)
    real(dp)                        :: expected_rhs(NQ)

    step = bdf(2, decay, H_STEP)
    step % qold   = Q_BE1
    step % qolder = Q0

    call solver % attach(step, ht, q, NQ, ncomp=1)

    ! c = -2 q1 + (1/2) q0, so the linear right-hand side is -c.
    expected_rhs = 2.0_dp * Q_BE1 - 0.5_dp * Q0

    call solver % constant(c)
    call report(size(c) .eq. NQ .and. &
         &      maxval(abs(c + expected_rhs)) .lt. TOL, &
         & "the bdf-2 affine constant is -2q1 + q0/2, so the linear " // &
         & "right-hand side is [5/3, 8/9]", nfail)

    call report(maxval(abs(expected_rhs - &
         &      [5.0_dp/3.0_dp, 8.0_dp/9.0_dp])) .lt. TOL, &
         & "which is the oracle, arrived at from the history rather " // &
         & "than copied in", nfail)

    rhs = field('rhs', q, NQ, ncomp=1)
    call rhs % set_real_vector(expected_rhs)

    call solver % apply(ht, [rhs], answer)

    d = answer % domain()
    call report(d % same_as(q), &
         & "the bdf-2 solution lands on Q as well", nfail)

    call answer % get_real_vector(v)
    call report(size(v) .eq. NQ .and. maxval(abs(v - Q_BDF2)) .lt. TOL, &
         & "q2 = [5/6, 47/72] - structural reach at Level 2, scheme " // &
         & "coefficients at Level 6, implicit solve here", nfail)

  end subroutine check_bdf2_solve

  !===================================================================!
  ! THE level's evidence about the minimizer's contract, stated as
  ! a refusal rather than as praise: the unknown domain is the
  ! CALLER'S word, and a caller that hands over the host's vertices
  ! gets a different - and here, an impossible - problem.
  !
  ! This is what "already domain-explicit" means operationally. The
  ! minimizer never guesses; it takes what it is given and holds the
  ! action to it.
  !===================================================================!

  subroutine check_unknown_domain_is_the_caller_s_word(nfail)

    integer, intent(inout) :: nfail

    type(step_operator)            :: step
    type(gmres)                    :: solver
    type(set_graph) :: d, hv
    integer         :: n_d

    step = backward_euler(decay, H_STEP)
    step % qold = Q0

    call solver % attach(step, ht, q, NQ, ncomp=1)
    call solver % domain(ht, d, n_d)
    hv = ht % vertex_set()

    call report(d % same_as(q) .and. .not. d % same_as(hv), &
         & "the unknown domain is exactly what the call site said - " // &
         & "Q - and never the host's vertices", nfail)

    ! And the residual domain the minimizer holds the action to is
    ! the action's own, which is why attach agreed to exist at all:
    ! a five-member residual against a two-member unknown would have
    ! been refused as a rectangular problem.
    call step % domain(ht, d, n_d)
    call report(d % same_as(q), &
         & "and the residual domain it validates against is the " // &
         & "step's, which is the action's: 2 unknowns, 2 residuals, " // &
         & "a square problem the solver family accepts", nfail)

  end subroutine check_unknown_domain_is_the_caller_s_word

end program time_level_7
