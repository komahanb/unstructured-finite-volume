!=====================================================================!
! Solver restriction and residual domain refusals, EXPECTED TO
! TERMINATE: a selection outside the whole domain, a repeated or
! empty selection, a split block, a selection beyond a solver's
! stated metadata, an elimination retaining nothing, a partition the
! selection misses, and a solve before the restricted operator is
! stated; then, on the residual boundary, a state, a design, a
! direction or a stored input defined on a graph of equal extent but
! another identity, a design with a value count other than the point
! count, and a host graph other than the residual's own. One case per
! argument.
!=====================================================================!
program restriction_refusal

  use iso_fortran_env, only : dp => REAL64
  use operation_stencil, only : stencil
  use operation_dense_direct, only : dense_direct
  use operation_gauss_seidel, only : gauss_seidel
  use operation_multigrid, only : multigrid
  use operation_elimination, only : elimination
  use operation_temporal_minimization, only : temporal_minimizer
  use operation_residual, only : residual_operator
  use operation_expression, only : unknown, derivative, constant, stated, operator(+), operator(*), operator(**)
  use operation_action, only : variation
  use view_directed_stored, only : stored_directed_graph
  use field_stored, only : stored_field, typed_field_domain
  use field_calculus, only : field

  implicit none

  type(stencil) :: a, tying
  type(dense_direct) :: factorisation
  type(gauss_seidel) :: sweeps
  type(multigrid) :: levels
  type(elimination) :: schur
  type(temporal_minimizer) :: solver
  type(residual_operator) :: residual
  type(stored_directed_graph) :: other
  type(typed_field_domain) :: fields
  type(stored_field) :: frozen(2), moved
  class(field), allocatable :: image
  real(dp) :: x(2), achieved, q(4)
  character(len=32) :: case_name

  call get_command_argument(1, case_name)
  a = stencil([1, 2, 3], [1, 2, 3], [2.0_dp, 2.0_dp, 2.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], 'diagonal')

  ! the implicit march q' = -q^3 over two instants, the first fixed
  tying    = stencil([4, 4, 4], [4, 3, 1], [1.0_dp, -10.0_dp, 10.0_dp], [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 'tying rows')
  residual = residual_operator(tying, &
       & stated(derivative(unknown(), 1) + constant(1.0_dp) * derivative(unknown(), 0) ** 3, 1, 'cubic decay'), &
       & [0, 2], 4, 2, [0], [1, 2], [1.0_dp, -1.0_dp])
  q      = [1.0_dp, -1.0_dp, 0.9_dp, -0.7_dp]
  frozen = residual % frozen_tuple(q, [0.0_dp, 0.0_dp])

  select case (trim(case_name))
  case ('state_domain')
     other  = stored_directed_graph(4, tails=[integer ::], heads=[integer ::])
     fields = typed_field_domain(other % vertex_set(), 4)
     moved  = fields % state(q)
     call residual % apply(residual % unknown_graph(), residual % bind([moved, frozen(2)]), image)
  case ('design_domain')
     other  = stored_directed_graph(2, tails=[integer ::], heads=[integer ::])
     fields = typed_field_domain(other % vertex_set(), 2)
     moved  = fields % design([0.0_dp, 0.0_dp])
     call residual % apply(residual % unknown_graph(), residual % bind([frozen(1), moved]), image)
  case ('design_count')
     fields = typed_field_domain(residual % design_domain(), 3)
     moved  = fields % design([0.0_dp, 0.0_dp, 0.0_dp])
     call residual % apply(residual % unknown_graph(), residual % bind([frozen(1), moved]), image)
  case ('host')
     other = stored_directed_graph(4, tails=[integer ::], heads=[integer ::])
     call residual % apply(other, residual % bind(frozen), image)
  case ('direction_domain')
     other  = stored_directed_graph(4, tails=[integer ::], heads=[integer ::])
     fields = typed_field_domain(other % vertex_set(), 4)
     moved  = fields % direction(q)
     call residual % partial_action(residual % unknown_graph(), residual % bind(frozen), &
          & [variation(residual % argument(1), moved)], image)
  case ('stored_domain')
     other  = stored_directed_graph(2, tails=[integer ::], heads=[integer ::])
     fields = typed_field_domain(other % vertex_set(), 2)
     moved  = fields % design([0.0_dp, 0.0_dp])
     allocate(solver % inner, source=factorisation)
     call solver % state(residual, residual % unknown_graph(), residual % unknown_domain(), 4, stored_inputs=[moved])
     call solver % partition([1, 1, 2, 2], [1, 2])
     call solver % solve([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], q, achieved)
  case ('outside')
     call factorisation % state(a, a % pattern, a % pattern % vertex_set(), 3)
     call factorisation % restrict([1, 4])
  case ('repeated')
     call factorisation % restrict([2, 2])
  case ('empty')
     call factorisation % restrict([integer ::])
  case ('split_block')
     sweeps % block_width = 2
     call sweeps % restrict([2, 3])
  case ('aggregates')
     levels % aggregates = [1, 1, 2, 2]
     call levels % restrict([5])
  case ('flags')
     schur % eliminated = [.false., .true., .false., .true.]
     call schur % restrict([1, 6])
  case ('no_retained')
     schur % eliminated = [.false., .true., .true., .true.]
     call schur % restrict([2, 3])
  case ('partition')
     call solver % partition([1, 1, 2, 2], [2, 1])
     call solver % restrict([6])
  case ('unstated')
     call factorisation % state(a, a % pattern, a % pattern % vertex_set(), 3)
     call factorisation % restrict([1, 3])
     x = 0.0_dp
     call factorisation % solve([1.0_dp, 1.0_dp], x, achieved)
  case default
     error stop 'unknown refusal case'
  end select
  print *, 'Invalid restriction was admitted.'

end program restriction_refusal
