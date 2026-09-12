!=====================================================================!
! Solver restriction refusals, EXPECTED TO TERMINATE: a selection
! outside the whole domain, a repeated or empty selection, a split
! block, a selection beyond a solver's stated metadata, an elimination
! retaining nothing, a partition the selection misses, and a solve
! before the restricted operator is stated. One case per argument.
!=====================================================================!
program restriction_refusal

  use iso_fortran_env, only : dp => REAL64
  use operation_stencil, only : stencil
  use operation_dense_direct, only : dense_direct
  use operation_gauss_seidel, only : gauss_seidel
  use operation_multigrid, only : multigrid
  use operation_elimination, only : elimination
  use operation_temporal_minimization, only : temporal_minimizer

  implicit none

  type(stencil) :: a
  type(dense_direct) :: factorisation
  type(gauss_seidel) :: sweeps
  type(multigrid) :: levels
  type(elimination) :: schur
  type(temporal_minimizer) :: solver
  real(dp) :: x(2), achieved
  character(len=32) :: case_name

  call get_command_argument(1, case_name)
  a = stencil([1, 2, 3], [1, 2, 3], [2.0_dp, 2.0_dp, 2.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], 'diagonal')

  select case (trim(case_name))
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
