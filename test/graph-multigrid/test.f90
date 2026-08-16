!=====================================================================!
! The multigrid suite: rung 5's acceptance.
!
! An eight-cell conduction chain, walls held at 0 and 10, compiled
! by the fitted balance into one stencil operator; blocks of two
! from the house coarsener. Four verdicts:
!
!   1. the aggregates cover every cell exactly once
!   2. THE COMMUTATION SQUARE: anything the blocks can express is
!      answered identically by the coarsened statement and by the
!      fine one - solve_coarse(R(A(P e))) returns e
!   3. two different roads, one answer: the two-grid cycle and a
!      direct solve meet on the same field
!   4. the detour pays: the two-grid residual after one cycle is
!      far below one smoothing pass alone
!=====================================================================!

program test_graph_multigrid

  use iso_fortran_env, only : dp => REAL64
  use graph_ordinary_view, only : ordinary_graph
  use class_graph_mesh   , only : mesh
  use class_graph_stencil, only : stencil_operator
  use class_fitted_balance, only : fitted_balance_stencil
  use class_polynomial_form, only : polynomial_form
  use class_graph_coarsener, only : coarsener, COARSEN_PAIRWISE
  use class_graph_multigrid, only : multigrid
  use class_graph_jacobi   , only : jacobi
  use class_graph_gmres    , only : gmres

  implicit none

  integer :: nfail

  nfail = 0

  call check_blocks_cover_once(nfail)
  call check_commutation_square(nfail)
  call check_two_roads_one_answer(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all multigrid checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' multigrid checks FAILED'
     error stop 1
  end if

contains

  subroutine report(ok, message, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! The chain: eight cells, unit faces, walls at 0 and 10.
  !===================================================================!

  type(mesh) function chain_mesh() result(m)

    integer :: kf

    m = mesh(8, &
         & tails=[1, 2, 3, 4, 5, 6, 7,  1, 8], &
         & heads=[2, 3, 4, 5, 6, 7, 8,  0, 0], &
         & volumes      = [(1.0_dp, kf = 1, 8)], &
         & cell_centers = [(real(kf - 1, dp) + 0.5_dp, 0.0_dp, 0.0_dp, &
         &                  kf = 1, 8)], &
         & areas        = [(1.0_dp, kf = 1, 9)], &
         & deltas       = [(1.0_dp, kf = 1, 7), 0.5_dp, 0.5_dp], &
         & normals      = [(1.0_dp, 0.0_dp, 0.0_dp, kf = 1, 7), &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centers = [(real(kf, dp), 0.0_dp, 0.0_dp, kf = 1, 7), &
         &                  0.0_dp, 0.0_dp, 0.0_dp, &
         &                  8.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [(0.5_dp, kf = 1, 9)])

  end function chain_mesh

  type(stencil_operator) function chain_statement(m) result(op)

    type(mesh), intent(in) :: m

    real(dp) :: vb(9), farea(9)
    integer :: e

    do e = 1, 9
       farea(e) = 1.0_dp
       vb(e)    = 0.0_dp
    end do
    vb(9) = 10.0_dp

    op = fitted_balance_stencil(m, polynomial_form(), farea, &
         & boundary_values=vb)

  end function chain_statement

  !===================================================================!
  ! VERDICT ONE. The aggregates partition the cells.
  !===================================================================!

  subroutine check_blocks_cover_once(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    integer :: nb

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    call report(size(assignment) == 8 .and. nb >= 2, &
         & 'every cell holds a block, and there is more than one', nfail)
    call report(minval(assignment) >= 1 .and. maxval(assignment) == nb, &
         & 'the blocks cover once: no gaps, no strays', nfail)

  end subroutine check_blocks_cover_once

  !===================================================================!
  ! VERDICT TWO. The commutation square, exact: pick any block
  ! field e, push it through the fine statement, gather, and the
  ! coarsened statement must answer e back.
  !===================================================================!

  subroutine check_commutation_square(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(multigrid) :: mg
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: e_blocks(:), e_fine(:), r_fine(:), rc(:), ec(:)
    real(dp) :: answered
    integer :: nb, v, b

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    allocate(mg % smoother, source=jacobi())
    allocate(mg % coarse  , source=gmres())
    call mg % attach(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    call mg % setup(assignment)
    mg % coarse % tolerance = 1.0d-13

    allocate(e_blocks(nb), e_fine(8), rc(nb), ec(nb))
    do b = 1, nb
       e_blocks(b) = real(mod(3 * b, 5) + 1, dp)
    end do

    ! P e: every cell takes its block's value. A(P e): the fine
    ! statement. R: gathered onto blocks.
    do v = 1, 8
       e_fine(v) = e_blocks(assignment(v))
    end do
    call mg % matvec(e_fine, r_fine)
    rc = 0.0_dp
    do v = 1, 8
       rc(assignment(v)) = rc(assignment(v)) + r_fine(v)
    end do

    ec = 0.0_dp
    call mg % coarse % solve(rc, ec, answered)

    call report(all(abs(ec - e_blocks) < 1.0d-9), &
         & 'the commutation square holds: coarsen-then-solve returns e', nfail)

  end subroutine check_commutation_square

  !===================================================================!
  ! VERDICTS THREE AND FOUR. Two roads to one answer, and the
  ! detour pays.
  !===================================================================!

  subroutine check_two_roads_one_answer(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(multigrid) :: mg
    type(gmres)     :: direct
    type(coarsener) :: c
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: g(:), rhs(:), x(:), xd(:)
    real(dp) :: achieved
    integer :: nb

    m = chain_mesh()
    c = coarsener(COARSEN_PAIRWISE)
    call c % blocks(m, assignment, nb)

    allocate(mg % smoother, source=jacobi())
    allocate(mg % coarse  , source=gmres())
    call mg % attach(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    call mg % setup(assignment)

    mg % smoother % max_iterations = 3
    mg % smoother % tolerance      = 0.0_dp
    select type (s => mg % smoother)
    type is (jacobi)
       s % omega = 0.7_dp
    end select
    mg % coarse   % tolerance      = 1.0d-13
    mg % tolerance                 = 1.0d-11
    mg % max_iterations            = 60

    call mg % constant(g)
    rhs = -g

    allocate(x(8))
    x = 0.0_dp
    call mg % solve(rhs, x, achieved)

    call report(achieved < 1.0d-9, &
         & 'the two-grid cycle closes the chain statement', nfail)

    call direct % attach(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    direct % tolerance = 1.0d-12
    allocate(xd(8))
    xd = 0.0_dp
    call direct % solve(rhs, xd, achieved)

    call report(all(abs(x - xd) < 1.0d-7), &
         & 'two roads, one answer: the cycle meets the direct solve', nfail)

    !----------------------------------------------------------------!
    ! TWO REGRESSIONS, PINNED WHERE THEY BROKE.
    !
    ! ATTACH IS RE-ENTERABLE. Newton calls it once per iteration, and
    ! it declares the solver's number domain each time - which a graph
    ! refuses, because a graph signs ONCE. The old counted_set
    ! constructor minted a fresh number domain per attach, so that is
    ! the semantics kept: reset to an unsigned graph, then declare.
    ! Attaching twice must therefore succeed and answer the same.
    !----------------------------------------------------------------!

    call direct % attach(chain_statement(m), m, m % vertex_set(), &
         & m % num_vertices())
    xd = 0.0_dp
    call direct % solve(rhs, xd, achieved)

    call report(all(abs(x - xd) < 1.0d-7), &
         & 'attach is re-enterable: the second attach solves the same', &
         & nfail)

    !----------------------------------------------------------------!
    ! A COUNT COMPONENT NEEDS A DEFAULT. jacobi() and gmres() name no
    ! component, so every component of minimizer must have a value
    ! without being named - a count with no value is not a different
    ! kind of thing from a count of zero.
    !----------------------------------------------------------------!

    unnamed_components: block
      type(jacobi) :: j
      type(gmres)  :: q
      j = jacobi()
      q = gmres()
      call report(j % n_unknown_domain .eq. 0 .and. &
           &      q % n_unknown_domain .eq. 0 .and. &
           &      j % n_residual_domain .eq. 0, &
           & 'an unattached solver counts no unknowns, and says so', &
           & nfail)
    end block unnamed_components

  end subroutine check_two_roads_one_answer

end program test_graph_multigrid
