!=====================================================================!
! THE VAN DER POL OSCILLATOR OF THE GRAPH TIME INTEGRATORS PAPER
! (Boopathy 2026, sections IV and V), stated in this framework:
!
!      u_tt - nu (1 - u^2) u_t + u = 0,   u(0) = 2,  u_t(0) = 0,
!      F(nu) = integral over [0, T] of nu u^2 / 2,   nu_0 = 1,
!
! and g_k = d^k F / d nu^k, k = 0 .. 4, by the expansion along nu and
! by the adjoint and its expansion, the integral by the tableau's
! weights at the stages. The paper's eq. 18 writes the integrand with
! u_t; its Table 4 is reproduced by u, the potential energy the text
! names, so u it is here.
!
! Table 4, IMID-2 at h = 0.02, is reproduced to every digit printed:
! 1.65208997, 1.79779439, 0.19985514, -0.22354341, 0.22194265. The
! three-window chain of Table 7 gives dF/dnu = 2.50287158092 against
! the paper's 2.502871579137, eight significant digits, the
! difference 1.8e-9 not resolved.
!
!      van_der_pol_study [mode] [order] [stages] [h] [degree] [checkpoint] [snapshots] [verbosity]
!
! mode 'table' (default): one family over [0, 1] at the step h, the
! paper's Table 4 at h = 0.02 - IMID-2 = dirk(2), DIRK-3 = dirk(3),
! DIRK-4 = dirk(4), SDIRK-4 = dirk(4, stages=5); mode 'chain': the
! same family over [0, 1] at the step h as a chain of one block per
! instant, the form the wavefront over the blocks and the orders
! pipelines on degree + 1 images: at h = 0.0005, SDIRK-4, one image
! takes 1.24 s at degree 0 and 0.55 s more per degree, 2.90 s at
! degree 3; the wavefronts on 2, 3, 4 images take 1.38, 1.51, 1.65 s
! at degrees 1, 2, 3, the same digits. Mode 'mixed' is the chain
! dirk, bdf, adams of the order given, one block per instant, whose
! multistep rules weigh history instants. With checkpoint = k or
! snapshots = s the state is stored at snapshots and the blocks are
! solved again in the reverse sweep: at h = 0.0005, degree 3, the
! coarray build takes 2.96 s with the state over the chain; 4.11 s
! with k = 45 (2001 block solves again), 7.79 s with s = 10 (8191
! again, the binomial count 10186 with the forward sweep) and 6.39 s
! with s = 20, the digits equal to fourteen places; two images, one
! solving again while the other forms the multipliers, take 3.08,
! 6.59 and 5.23 s;
! mode 'windows': the
! three-window chain of the paper's Table 7 over [0, 1.5], SDIRK-4 at
! h = 0.025, DIRK-4 at h = 0.02, SDIRK-4 at h = 0.0125, the instants
! listed, and dF/dnu against the paper's 2.502871579137.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol_study

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, instants, parameter, expansion
  use operation_field   , only : continuous_field, discrete_field, integral, &
       &                         operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family  , only : family, dirk, bdf, adams, chain

  implicit none

  real(dp), parameter :: nu_design = 1.0_dp, u_initial = 2.0_dp, v_initial = 0.0_dp
  integer :: degree = 4, checkpoint = 0, snapshots = 0, verbose = 0

  character(len=16) :: mode = 'table'
  integer  :: order = 2, stages = 0
  real(dp) :: h = 0.02_dp, T_final
  real(dp), allocatable :: times(:)
  type(family), allocatable :: schemes(:)
  integer, allocatable :: from(:)

  type(continuous_manifold) :: omega, domega
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, u, lambda_r, lambda, kappa, oscillator, u0, v0, energy, rest
  type(discrete_field)      :: estimate, solution, kappa_h, sensitivity
  type(continuous_residual) :: r, g, d, L
  type(discrete_residual)   :: L_h
  real(dp), allocatable :: values(:), reverse(:)
  character(len=32) :: item
  integer :: k, n

  call arguments()
  if (this_image() == 1) call set_verbosity(verbose)
  allocate(values(1 + degree), reverse(0:degree))

  select case (trim(mode))
  case ('table')
     T_final = 1.0_dp
     n = nint(T_final / h) + 1
     times = [(real(k - 1, dp) * h, k = 1, n)]
     allocate(schemes(1), from(1))
     schemes(1) = tableau(order, stages)
     from(1) = 1
  case ('chain')
     T_final = 1.0_dp
     n = nint(T_final / h) + 1
     times = [(real(k - 1, dp) * h, k = 1, n)]
     allocate(schemes(n), from(n))
     do k = 1, n
        schemes(k) = tableau(order, stages)
        from(k) = k
     end do
  case ('mixed')
     ! dirk over the first three instants, bdf over the next four,
     ! adams to the end, each of the order given, one block per instant
     T_final = 1.0_dp
     n = nint(T_final / h) + 1
     times = [(real(k - 1, dp) * h, k = 1, n)]
     allocate(schemes(n), from(n))
     do k = 1, n
        if (k <= 3) then
           schemes(k) = dirk(order)
        else if (k <= 7) then
           schemes(k) = bdf(order)
        else
           schemes(k) = adams(order)
        end if
        from(k) = k
     end do
  case ('windows')
     T_final = 1.5_dp
     times = [(0.025_dp * real(k, dp), k = 0, 19), (0.5_dp + 0.02_dp * real(k, dp), k = 0, 24), &
          &   (1.0_dp + 0.0125_dp * real(k, dp), k = 0, 40)]
     n = size(times)
     allocate(schemes(3), from(3))
     schemes(1) = dirk(4, stages=5); from(1) = 1
     schemes(2) = dirk(4);           from(2) = 21
     schemes(3) = dirk(4, stages=5); from(3) = 46
  case default
     error stop 'van_der_pol_study: the mode is table, chain, mixed or windows'
  end select
  if (this_image() == 1) then
     print '(a,a,a,i0,a,i0,a,i0)', 'mode ', trim(mode), '  instants ', n, '  dirk order ', order, '  stages ', stages
  end if

  omega   = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_design))
  domega  = omega % boundary(time=0.0_dp)
  omega_h = omega % discretize(time=instants(times), design=expansion(order=degree))
  t  = omega % coordinate('t')
  nu = omega % coordinate('nu')
  u  = omega % unknown('u', [t, nu])

  oscillator = u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u
  u0 = u - u_initial
  v0 = u % derivative([t]) - v_initial
  energy   = integral(nu*u*u/2.0_dp, over=omega % time())
  r        = continuous_residual(omega,  [oscillator])
  g        = continuous_residual(domega, [u0, v0])
  d        = continuous_residual(omega % design(), [nu - nu_design])
  lambda_r = r % multiplier('adjoint', [t, nu])
  lambda   = g % multiplier('lambda', [nu])
  kappa    = d % multiplier('kappa')
  L        = energy + lambda_r*r + lambda*g + kappa*d
  L_h      = L % discretize(omega_h, time=chain(schemes, from=from))

  rest     = continuous_field(omega, [u_initial])
  estimate = rest % discretize(omega_h)
  ! checkpoint > 0: the stages retained over segments of that many
  ! blocks alone, formed again in the reverse sweep
  if (checkpoint > 0) then
     call L_h % minimize(estimate, solution, checkpoint=checkpoint)
  else if (snapshots > 0) then
     call L_h % minimize(estimate, solution, snapshots=snapshots)
  else
     call L_h % minimize(estimate, solution)
  end if

  sensitivity = energy % at(solution)
  call sensitivity % values(values)
  kappa_h = solution % fields(['kappa'])
  reverse(0) = -kappa_h % value(1, 1)
  do k = 1, degree - 1
     kappa_h    = kappa_h % derivative([nu])
     reverse(k) = -kappa_h % value(1, 1)
  end do

  ! the table, printed by image 1: every image has the same solution
  if (this_image() == 1) then
     print '(a)', '  k   g_k by the expansion       g_k by the adjoint expansion'
     print '(i3, es24.16)', 0, values(1)
     do k = 1, degree
        if (k < degree) then
           print '(i3, 2es24.16)', k, values(1 + k), reverse(k - 1)
        else
           print '(i3, es24.16, a)', k, values(1 + k), '   (the adjoint expansion of order ' // trim(item_of(degree - 1)) // ')'
        end if
     end do
     if (trim(mode) == 'table') then
        print '(a)', 'the paper''s Table 4 at h = 0.02, IMID-2: 1.65208997, 1.79779439, 0.19985514, -0.22354341, 0.22194265'
     else if (trim(mode) == 'windows') then
        print '(a)', 'the paper''s eq. 20: dF/dnu = 2.502871579137'
     end if
  end if

contains

  function tableau(order, stages) result(scheme)
    integer, intent(in) :: order, stages
    type(family) :: scheme
    if (stages > 0) then
       scheme = dirk(order, stages=stages)
    else
       scheme = dirk(order)
    end if
  end function tableau

  function item_of(i) result(text)
    integer, intent(in) :: i
    character(len=8) :: text
    write(text, '(i0)') i
  end function item_of

  subroutine arguments()
    integer :: count
    count = command_argument_count()
    if (count >= 1) call get_command_argument(1, mode)
    if (count >= 2) then; call get_command_argument(2, item); read(item, *) order;  end if
    if (count >= 3) then; call get_command_argument(3, item); read(item, *) stages; end if
    if (count >= 4) then; call get_command_argument(4, item); read(item, *) h;      end if
    if (count >= 5) then; call get_command_argument(5, item); read(item, *) degree; end if
    if (count >= 6) then; call get_command_argument(6, item); read(item, *) checkpoint; end if
    if (count >= 7) then; call get_command_argument(7, item); read(item, *) snapshots; end if
    if (count >= 8) then; call get_command_argument(8, item); read(item, *) verbose; end if
  end subroutine arguments

end program van_der_pol_study
