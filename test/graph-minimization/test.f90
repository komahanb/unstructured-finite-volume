!=====================================================================!
! The minimization suite: rung 4's acceptance.
!
! The whole tower stands in one problem: heat conduction on a chain
! of three cells, walls held at 0 and 10,
!
!       0 |--- (1) --- (2) --- (3) ---| 10        k = 1, A = 1
!         d 0.5     1       1     0.5
!
! whose exact answer is q = [5/3, 5, 25/3]. The mesh seats the
! measurements, the constitution supplies the coefficients, the
! calculus builds the rows, and the minimization drives the residual
! down - every level doing only its own job, and the solver speaking
! only its own words: matvec, inner product, norm, diagonal, sweep
! order, each a delegation to an engine seat.
!=====================================================================!

program test_graph_minimization

  use iso_fortran_env, only : dp => REAL64
  use graph_grammar  , only : graph
  use class_graph_mesh   , only : mesh
  use class_mesh_builder , only : mesh_from_gmsh
  use class_robin_condition, only : robin_condition, dirichlet
  use class_conduction     , only : conduction
  use class_graph_differential_operator, only : differential_operator
  use class_graph_differential_operator, only : vertex_differential_operator
  use class_graph_jacobi   , only : jacobi
  use class_graph_conjugate_gradient, only : conjugate_gradient

  implicit none

  integer :: nfail

  nfail = 0

  call check_solver_words(nfail)
  call check_hand_problem(nfail)
  call check_real_mesh(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all minimization checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' minimization checks FAILED'
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
  ! The three-cell conduction problem, assembled from the levels.
  !===================================================================!

  type(mesh) function chain_mesh() result(m)

    m = mesh(3, tails=[1, 2, 1, 3], heads=[2, 3, 0, 0], &
         & volumes      = [1.0_dp, 1.0_dp, 1.0_dp], &
         & cell_centers = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp, &
         &                 2.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
         & deltas       = [1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centers = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp, &
         &                 3.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 0.5_dp, 1.0_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', '', 'west', 'east'])

  end function chain_mesh

  !===================================================================!
  ! The assembly: the constitution's coefficients placed into the
  ! operator's argument list - conduction inside, one condition per
  ! wall - and nothing else.
  !===================================================================!

  type(differential_operator) function chain_operator(m) result(op)

    type(mesh), intent(in) :: m

    type(robin_condition) :: west, east
    real(dp), allocatable :: c(:), b(:)

    call assemble(m, conduction(1.0_dp), &
         & [dirichlet('west', 0.0_dp), dirichlet('east', 10.0_dp)], &
         & 1.0_dp, c, b)

    op = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=[1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp], boundary_values=b)

  end function chain_operator

  subroutine assemble(m, law, conditions, kappa, c, b)

    type(mesh), intent(in)            :: m
    type(conduction), intent(in)      :: law
    type(robin_condition), intent(in) :: conditions(:)
    real(dp), intent(in)              :: kappa
    real(dp), allocatable, intent(out) :: c(:), b(:)

    class(graph), allocatable :: members
    real(dp), allocatable :: cw(:), bw(:)
    integer :: k, f, e

    call law % edge_coefficients(m, c)
    allocate(b(size(c)))
    b = 0.0_dp

    do k = 1, size(conditions)
       call conditions(k) % faces(m, members)
       call conditions(k) % operator_coefficients(m, kappa, cw)
       call conditions(k) % boundary_values(m, bw)
       do f = 1, members % num_vertices()
          e = members % global_vertex_index(f)
          c(e) = cw(f)
          b(e) = bw(f)
       end do
    end do

  end subroutine assemble

  !===================================================================!
  ! The solver's words, checked against hand values on the chain:
  ! the probed diagonal is the stencil's own, and the colouring
  ! never gives neighbours one colour.
  !===================================================================!

  subroutine check_solver_words(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: d(:)
    integer , allocatable :: colours(:), nbrs(:)
    logical :: proper
    integer :: v, i

    m = chain_mesh()
    call js % attach(chain_operator(m), m)

    call js % diagonal(d)
    call report(all(abs(d - [-3.0_dp, -2.0_dp, -3.0_dp]) < 1.0d-12), &
         & 'the probed diagonal is the stencil diagonal, by hand', nfail)

    call js % sweep_order(colours)
    proper = .true.
    do v = 1, m % num_vertices()
       call m % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (colours(nbrs(i)) == colours(v)) proper = .false.
       end do
    end do
    call report(proper, 'the sweep order never gives neighbours one colour', nfail)

    call report(abs(js % inner_product([1.0_dp, 2.0_dp, 3.0_dp], &
         &                             [2.0_dp, 1.0_dp, 0.0_dp]) - 4.0_dp) &
         & < 1.0d-14, 'the inner product is the measured sum', nfail)

    call report(abs(js % norm([3.0_dp, 4.0_dp, 0.0_dp]) - 5.0_dp) < 1.0d-14, &
         & 'the norm is the norm', nfail)

  end subroutine check_solver_words

  !===================================================================!
  ! The hand problem, driven to its exact answer.
  !===================================================================!

  subroutine check_hand_problem(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(jacobi) :: js
    real(dp), allocatable :: g(:), rhs(:), x(:), y(:)
    real(dp) :: achieved
    real(dp), parameter :: exact(3) = [5.0_dp/3.0_dp, 5.0_dp, 25.0_dp/3.0_dp]

    m = chain_mesh()
    call js % attach(chain_operator(m), m)
    js % max_iterations = 5000
    js % tolerance      = 1.0d-12

    ! The equation action(q) = 0 reads matvec(q) = -constant.
    call js % constant(g)
    rhs = -g

    allocate(x(3))
    x = 0.0_dp
    call js % solve(rhs, x, achieved)

    call report(achieved < 1.0d-10, 'jacobi drives the residual down', nfail)
    call report(all(abs(x - exact) < 1.0d-8), &
         & 'and lands on the exact answer: 5/3, 5, 25/3', nfail)

    call js % matvec(x, y)
    call report(js % norm(rhs - y) < 1.0d-10, &
         & 'the answer satisfies the assembled equation', nfail)

  end subroutine check_hand_problem

  !===================================================================!
  ! The real mesh. A symmetric positive operator - the conduction
  ! coefficients negated - and a consistent right hand side
  ! manufactured from a known state; conjugate gradient must return
  ! a state the operator cannot tell apart from it.
  !===================================================================!

  subroutine check_real_mesh(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(conjugate_gradient) :: cg
    type(conduction) :: law
    real(dp), allocatable :: c(:), deltas(:), xref(:), b(:), x(:), y(:)
    real(dp) :: achieved
    integer :: v

    m = mesh_from_gmsh('../square-10.msh')

    law = conduction(1.0_dp)
    call law % edge_coefficients(m, c)
    c = -c

    block
      use class_graph_field, only : field
      type(field) :: fd
      fd = m % face_delta()
      call fd % get_real_vector(deltas)
    end block

    call cg % attach(vertex_differential_operator(order=2, &
         & coefficients=c, spacings=deltas), m)
    cg % max_iterations = 2000
    cg % tolerance      = 1.0d-10

    allocate(xref(m % num_vertices()))
    do v = 1, size(xref)
       xref(v) = real(mod(3 * v, 17), dp)
    end do

    call cg % matvec(xref, b)

    allocate(x(size(xref)))
    x = 0.0_dp
    call cg % solve(b, x, achieved)

    call cg % matvec(x, y)
    call report(cg % norm(b - y) < 1.0d-7 * (1.0_dp + cg % norm(b)), &
         & 'conjugate gradient closes the manufactured equation on a real mesh', &
         & nfail)

  end subroutine check_real_mesh

end program test_graph_minimization
