!=====================================================================!
! THE TAYLOR-GREEN VORTEX under chosen parameters: the same statement
! as taylor_green_vortex.f90, with the duration, the number of
! instants, the order of the finite differences, the orders of the
! three time families and the instants their blocks begin at read
! from the command line, and the convergence of every solve printed.
!
!      taylor_green_study [T] [instants] [order] [dirk] [bdf] [adams] [from2] [from3] [mesh] [length] [method]
!                         [adjoint] [checkpoint] [snapshots] [stages]
!
! Every argument has the value of taylor_green_vortex.f90 when
! absent: 1.0 10 2 2 2 2 5 8 box.msh 0 difference. The tenth, when
! positive, is the number of instants per Adams block: the Adams
! family is repeated in blocks of that length from the third block's
! first instant to the end, so that a long duration is solved block
! by block rather than as one system; when negative, the chain is
! the DIRK family alone, one block per instant, with the number of
! stages of the fifteenth argument when that is positive. The
! eleventh names the spatial method, difference or volume. The
! twelfth, when one, pairs the equations with an adjoint and states
! the kinetic energy integrated over the manifold as the objective,
! so that the reverse sweep runs; the thirteenth and the fourteenth
! are the schedule of the sweeps, a snapshot every checkpoint blocks
! or snapshots placed by the binomial schedule. 16 x 16 cells, finite
! volumes of order 2, SDIRK-4 (five stages), 100 instants, the
! adjoint: 118 s and a peak resident memory of 256 MB with the state
! over the chain; 186 s and 183 MB with checkpoint = 10 (100 block
! solves again); 315 s and 175 MB with snapshots = 5 (285 again, the
! binomial count 380 with the forward sweep); 175 MB is the peak of
! five instants, one block's solve, so that nothing in proportion to
! the chain is left; J and the reactions equal in every digit.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program taylor_green_study

  use util_precision             , only : dp
  use iso_fortran_env            , only : int64
  use util_verbosity             , only : set_verbosity
  use operation_manifold         , only : continuous_manifold, discrete_manifold, interval, region, instants, mesh
  use operation_field            , only : continuous_field, discrete_field, integral, sin, cos, exp, &
       &                                  operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual         , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family           , only : family, dirk, bdf, adams, chain
  use operation_finite_difference, only : finite_difference
  use operation_finite_volume    , only : finite_volume

  implicit none

  real(dp), parameter :: nu = 0.01_dp

  real(dp) :: T_final = 1.0_dp
  integer  :: num_instants = 10, order = 2, dirk_order = 2, bdf_order = 2, adams_order = 2, from2 = 5, from3 = 8
  integer  :: length = 0, with_adjoint = 0, checkpoint = 0, snapshots = 0, stages = 0
  character(len=256) :: mesh_file = 'box.msh'
  character(len=16)  :: method = 'difference'
  type(family), allocatable :: schemes(:)
  integer     , allocatable :: from(:)
  integer :: k, num_adams

  type(continuous_manifold) :: omega, domega, tau
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, y, u_exact, v_exact, p_exact, exact
  type(continuous_field)    :: u, v, p, lambda, mu, momentum_x, momentum_y, pressure, lambda_r, energy
  type(discrete_field)      :: estimate, solution, error, functional
  real(dp) :: J(1)
  type(continuous_residual) :: r, g, gauge, L
  type(discrete_residual)   :: L_h
  integer(int64) :: t0, t1, rate

  call arguments()

  ! the chain: dirk, bdf, then adams once to the end or in blocks of
  ! the length given
  if (length < 0) then
     allocate(schemes(num_instants), from(num_instants))
     do k = 1, num_instants
        if (stages > 0) then
           schemes(k) = dirk(dirk_order, stages=stages)
        else
           schemes(k) = dirk(dirk_order)
        end if
        from(k) = k
     end do
  else
     num_adams = 1
     if (length > 0) num_adams = max(1, (num_instants - from3) / length + 1)
     allocate(schemes(2 + num_adams), from(2 + num_adams))
     schemes(1) = dirk(dirk_order); from(1) = 1
     schemes(2) = bdf(bdf_order);   from(2) = from2
     do k = 1, num_adams
        schemes(2 + k) = adams(adams_order)
        from(2 + k)    = from3 + (k - 1) * length
     end do
  end if

  print '(a,f6.2,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,a,a,*(i0,1x))', 'study  T = ', T_final, '  instants = ', num_instants, &
       & '  ', trim(method), ' of order ', order, '  chain dirk(', dirk_order, ') bdf(', bdf_order, &
       & ') adams(', adams_order, ')  mesh ', trim(mesh_file), '  blocks from ', from

  omega   = continuous_manifold(time=interval(0.0_dp, T_final), space=region('box.geo'))
  domega  = omega % boundary(time=0.0_dp)
  tau     = omega % time()
  omega_h = omega % discretize(time=instants(num_instants), space=mesh(trim(mesh_file)))

  t = omega % coordinate('t')
  x = omega % coordinate('x')
  y = omega % coordinate('y')

  u_exact =  sin(x)*cos(y)*exp(-2.0_dp*nu*t)
  v_exact = -cos(x)*sin(y)*exp(-2.0_dp*nu*t)
  p_exact = (cos(2.0_dp*x) + cos(2.0_dp*y))*exp(-4.0_dp*nu*t)/4.0_dp
  exact   = continuous_field(omega, [u_exact, v_exact, p_exact])

  u      = omega % unknown('u', [t, x, y])
  v      = omega % unknown('v', [t, x, y])
  p      = omega % unknown('p', [t, x, y])

  momentum_x = u % derivative([t]) + u*u % derivative([x]) + v*u % derivative([y]) + p % derivative([x]) &
       &     - nu*(u % derivative([x, x]) + u % derivative([y, y]))
  momentum_y = v % derivative([t]) + u*v % derivative([x]) + v*v % derivative([y]) + p % derivative([y]) &
       &     - nu*(v % derivative([x, x]) + v % derivative([y, y]))
  pressure   = p % derivative([x, x]) + p % derivative([y, y]) &
       &     + u % derivative([x])**2 + 2.0_dp*u % derivative([y])*v % derivative([x]) + v % derivative([y])**2

  r      = continuous_residual(omega,  [momentum_x, momentum_y, pressure])
  g      = continuous_residual(domega, [u - u_exact, v - v_exact])
  gauge  = continuous_residual(tau,    [integral(p, over=omega % space())])
  lambda = g % multiplier('lambda', [x, y])
  mu     = gauge % multiplier('mu', [t])
  L      = r + lambda*g + mu*gauge
  if (with_adjoint == 1) then
     lambda_r = r % multiplier('adjoint', [t, x, y])
     energy   = integral((u*u + v*v)/2.0_dp, over=omega)
     L        = energy + lambda_r*r + lambda*g + mu*gauge
  end if

  select case (trim(method))
  case ('difference')
     L_h = L % discretize(omega_h, time=chain(schemes, from=from), space=finite_difference(order=order))
  case ('volume')
     L_h = L % discretize(omega_h, time=chain(schemes, from=from), space=finite_volume(order=order))
  case default
     error stop 'taylor_green_study: the spatial method is difference or volume'
  end select

  if (this_image() == 1) call set_verbosity(1)
  estimate = exact % discretize(omega_h)
  call system_clock(t0, rate)
  if (checkpoint > 0) then
     call L_h % minimize(estimate, solution, checkpoint=checkpoint)
  else if (snapshots > 0) then
     call L_h % minimize(estimate, solution, snapshots=snapshots)
  else
     call L_h % minimize(estimate, solution)
  end if
  call system_clock(t1)

  error = solution % fields(['u', 'v', 'p']) - estimate
  print '(a, es12.4, a, f8.2, a)', 'error against the exact solution ', error % norm(), &
       & '   minimize ', real(t1 - t0, dp) / real(rate, dp), ' s'
  if (with_adjoint == 1) then
     functional = energy % at(solution)
     call functional % values(J)
     error = solution % fields(['lambda'])
     print '(a, es24.16, a, es24.16)', 'J = ', J(1), '   the norm of the reactions of the initial data ', error % norm()
  end if

contains

  ! the command line: each argument present replaces its default
  subroutine arguments()
    character(len=256) :: item
    integer :: n
    n = command_argument_count()
    if (n >= 1) then; call get_command_argument(1, item); read(item, *) T_final;      end if
    if (n >= 2) then; call get_command_argument(2, item); read(item, *) num_instants; end if
    if (n >= 3) then; call get_command_argument(3, item); read(item, *) order;        end if
    if (n >= 4) then; call get_command_argument(4, item); read(item, *) dirk_order;   end if
    if (n >= 5) then; call get_command_argument(5, item); read(item, *) bdf_order;    end if
    if (n >= 6) then; call get_command_argument(6, item); read(item, *) adams_order;  end if
    if (n >= 7) then; call get_command_argument(7, item); read(item, *) from2;        end if
    if (n >= 8) then; call get_command_argument(8, item); read(item, *) from3;        end if
    if (n >= 9) then; call get_command_argument(9, mesh_file);                        end if
    if (n >= 10) then; call get_command_argument(10, item); read(item, *) length;     end if
    if (n >= 11) then; call get_command_argument(11, method);                        end if
    if (n >= 12) then; call get_command_argument(12, item); read(item, *) with_adjoint; end if
    if (n >= 13) then; call get_command_argument(13, item); read(item, *) checkpoint;   end if
    if (n >= 14) then; call get_command_argument(14, item); read(item, *) snapshots;    end if
    if (n >= 15) then; call get_command_argument(15, item); read(item, *) stages;       end if
  end subroutine arguments

end program taylor_green_study
