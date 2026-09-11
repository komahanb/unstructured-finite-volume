program context_contract
  use util_precision, only : dp
  use iso_fortran_env, only : real128
  use operation_minimization, only : minimizer, absolute, by_count
  use operation_gmres, only : gmres
  use operation_dense_direct, only : dense_direct
  use gti_march, only : march_context, precision_needed
  implicit none
  type(march_context) :: a, b, independent
  class(minimizer), allocatable :: inner
  integer :: version_a, version_b
  real(real128) :: spacing_a, spacing_b
  character(len=:), allocatable :: kind_name
  call a % set_linear_solver('iterative')
  call a % set_storage('sparse')
  call a % set_linear_limits(4, 2, 7)
  call a % set_stopping(1.0e-7_dp, absolute, by_count, 9)
  call a % set_predictor_order(2)
  call a % set_newton_order(3)
  call a % set_coarse_nodes([1, 1, 2, 2])
  call a % set_space_coupling('sequential')
  call a % set_time_coupling('coupled')
  call a % read_inner(inner, 8, 2, spread(.false., 1, 8))
  select type(inner)
  type is(gmres)
     if (inner % restart /= 4 .or. inner % max_iterations /= 7) error stop 'initial linear configuration'
  class default
     error stop 'context A uses GMRES'
  end select
  inner % tolerance = 0.314_dp
  call a % store_inner(inner)
  call b % set_stopping(1.0e-3_dp, absolute, by_count, 3)
  call b % read_inner(inner, 8, 2, spread(.false., 1, 8))
  select type(inner)
  type is(dense_direct)
  class default
     error stop 'context B uses direct factorisation'
  end select
  call b % store_inner(inner)
  call a % read_inner(inner, 8, 2, spread(.false., 1, 8))
  if (inner % tolerance /= 0.314_dp) error stop 'B changed A retained minimizer'
  call a % store_inner(inner)
  version_a = a % next_version()
  version_b = b % next_version()
  if (version_a /= 1 .or. version_b /= 1) error stop 'independent versions'
  independent = a % configuration()
  if (independent % next_version() /= 1 .or. a % next_version() /= 2) error stop 'fresh version sequence'
  if (independent % predictor_order() /= 2 .or. independent % newton_order() /= 3) error stop 'fresh orders'
  if (any(independent % coarse_nodes(4) /= [1,1,2,2])) error stop 'fresh coarse members'
  if (independent % coupling_named() /= a % coupling_named()) error stop 'fresh coupling'
  call independent % read_inner(inner, 8, 2, spread(.false., 1, 8))
  if (inner % tolerance /= 1.0e-7_dp) error stop 'fresh copied a retained solver'
  call a % read_inner(inner, 8, 2, spread(.false., 1, 8))
  if (inner % tolerance /= 0.314_dp) error stop 'fresh changed source solver'
  call a % store_inner(inner)
  call a % set_linear_limits(2, 1, 3)
  call a % read_inner(inner, 8, 2, spread(.false., 1, 8))
  select type(inner)
  type is(gmres)
     if (inner % restart /= 2 .or. inner % max_iterations /= 3) error stop 'linear limit invalidation'
  end select
  call a % store_inner(inner)
  call a % set_linear_stopping(1.0e-6_dp, absolute, by_count)
  call a % read_inner(inner, 8, 2, spread(.false., 1, 8))
  if (inner % tolerance /= 1.0e-6_dp) error stop 'linear tolerance invalidation'
  call a % store_inner(inner)
  call a % read_inner(inner, 1, 1, [.false.])
  select type(inner)
  type is(gmres)
     if (inner % restart /= 1) error stop 'dimension invalidation'
  end select
  call precision_needed(1.0_dp, 1.0_dp, 1.0_dp, spacing_a, kind_name, context=a)
  call precision_needed(1.0_dp, 1.0_dp, 1.0_dp, spacing_b, kind_name, context=b)
  if (abs(real(spacing_a, dp) / 1.0e-7_dp - 1.0_dp) > 1.0e-14_dp) error stop 'A nonlinear tolerance'
  if (abs(real(spacing_b, dp) / 1.0e-3_dp - 1.0_dp) > 1.0e-14_dp) error stop 'B nonlinear tolerance'
  call precision_needed(1.0_dp, 1.0_dp, 1.0_dp, spacing_a, kind_name)
  if (abs(real(spacing_a, dp) / 1.0e-12_dp - 1.0_dp) > 1.0e-14_dp) error stop 'standalone default tolerance'
  print *, 'PASS: independent settings, retained solvers, fresh configuration, versions and invalidation'
end program context_contract
