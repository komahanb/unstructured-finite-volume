! The implicit function theorem on a space residual, solved once with
! no chain: the same implicit_forward and implicit_dual the marcher
! makes at every instant, here with a conjugate-gradient solver.
!
!     R_i(q, xi) = q_i^2 - xi = 0   on three points,  q = sqrt(xi)
!
program law2
  use iso_fortran_env, only : dp => REAL64
  use view_directed_stored, only : stored_directed_graph
  use graph_fractal, only : graph
  use field_stored, only : stored_field
  use operation_chain_rule, only : argument_path
  use operation_linearization, only : implicit_forward, implicit_dual
  use operation_newton, only : newton
  use operation_conjugate_gradient, only : conjugate_gradient
  use operation_gmres, only : gmres
  use toy_differentiable_forms, only : equilibrium_law, fill_path
  implicit none

  type(stored_directed_graph) :: three, lone
  type(graph) :: points, cells
  type(equilibrium_law) :: equil
  type(newton) :: solver
  type(conjugate_gradient) :: cg
  type(stored_field) :: qf, xif
  type(stored_field), allocatable :: inputs(:), duals(:)
  type(argument_path) :: xipath(1)
  real(dp), allocatable :: derivatives(:,:), lambda(:), g(:)
  real(dp) :: q(3), achieved, xi
  integer :: nfail

  nfail = 0
  xi    = 4.0_dp

  three  = stored_directed_graph(3, tails=[integer ::], heads=[integer ::])
  lone   = stored_directed_graph(1, tails=[integer ::], heads=[integer ::])
  points = three % vertex_set()
  cells  = lone % vertex_set()
  equil  = equilibrium_law(xi)

  xif = stored_field('xi', cells, 1, num_components=1)
  call xif % set_real_vector([xi])

  ! 1. local invertibility: newton = minimizer o tangent_of, xi held
  allocate(solver % inner, source=gmres())
  solver % inner % tolerance = 1.0d-14
  solver % tolerance = 1.0d-12
  call solver % attach(equil, three, points, 3, held_inputs=[xif])
  q = [3.0_dp, 1.0_dp, 5.0_dp]
  call solver % solve([0.0_dp, 0.0_dp, 0.0_dp], q, achieved)
  call check(maxval(abs(q - sqrt(xi))) < 1.0d-9, 'newton solves q^2 = xi on three points: q = 2')

  qf = stored_field('q', points, 3, num_components=1)
  call qf % set_real_vector(q)
  inputs = [qf, xif]

  ! 2. forward: dq/dxi to order 3 along xi' = 1, through CG
  call fill_path(xipath(1), equil % argument(2), [1.0_dp, 0.0_dp, 0.0_dp], cells)
  cg % tolerance = 1.0d-14
  cg % max_iterations = 50
  call implicit_forward(equil, three, inputs, equil % argument(1), 3, xipath, &
       & derivatives, solver=cg)
  call check(maxval(abs(derivatives(:, 1) - 0.5_dp / sqrt(xi)))        < 1.0d-9, &
       & "q'   = 1/(2 sqrt xi)        = 0.25")
  call check(maxval(abs(derivatives(:, 2) + 0.25_dp / xi ** 1.5_dp))   < 1.0d-9, &
       & "q''  = -1/(4 xi^3/2)        = -0.03125")
  call check(maxval(abs(derivatives(:, 3) - 0.375_dp / xi ** 2.5_dp))  < 1.0d-9, &
       & "q''' = 3/(8 xi^5/2)         = 0.01171875")

  ! 3. dual: J = sum q, seed = dJ/dq = 1; lambda = A^-T seed; g = (D_xi R)^T lambda
  call implicit_dual(equil, three, inputs, equil % argument(1), [1.0_dp, 1.0_dp, 1.0_dp], &
       & lambda, wrt=[equil % argument(2)], duals=duals, solver=cg)
  call duals(1) % real_vector(g)
  call check(maxval(abs(lambda - 0.5_dp / sqrt(xi))) < 1.0d-9, 'lambda = seed / (2q) = 0.25')
  call check(size(g) == 1 .and. abs(g(1) + 3.0_dp * 0.5_dp / sqrt(xi)) < 1.0d-9, &
       & 'the dual in xi lives on the one-entry xi domain: (D_xi R)^T lambda = -0.75')

  ! 4. the pairing law between the two roads: <seed, q'> = <-g, xi'>
  call check(abs(sum(derivatives(:, 1)) + g(1)) < 1.0d-12, &
       & "<dJ/dq, dq/dxi> = -<(D_xi R)^T lambda, 1>: forward and dual agree")

  if (nfail == 0) then
     print '(a)', ' PASS : the implicit function theorem holds on a space residual through the marcher''s two solves'
  else
     print '(a,i0,a)', ' FAIL : ', nfail, ' checks'
     error stop
  end if

contains

  subroutine check(ok, what)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: what
    if (ok) then
       print '(a)', ' PASS : ' // what
    else
       print '(a)', ' FAIL : ' // what
       nfail = nfail + 1
    end if
  end subroutine check

end program law2
