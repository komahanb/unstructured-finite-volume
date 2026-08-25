! The three families' coefficients read off their couplings.
!
! bdf and adams on a uniform grid beside the tabulated values, and on
! a non-uniform grid beside the variable-step row operation_step
! carries for bdf of order two. dirk on a stage coupling beside its
! tableau.
program family_coefficients

  use util_precision  , only : dp
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_action      , only : variation
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : dirk_family, crouzeix_two_stage

  implicit none

  integer, parameter :: order = 2
  integer :: k

  call bdf_on('uniform',     [(0.5_dp, k = 1, 2 * order + 1)])
  call bdf_on('non-uniform', [0.0_dp, 0.3_dp, 0.2_dp, 0.4_dp, 0.25_dp])
  call adams_on('uniform',     3, [(0.5_dp, k = 1, 3)])
  call adams_on('non-uniform', 3, [0.0_dp, 0.3_dp, 0.2_dp])
  call dirk_on(crouzeix_two_stage())
  call bdf_step_sensitivity([0.0_dp, 0.3_dp, 0.2_dp, 0.4_dp, 0.25_dp])

contains

  !-------------------------------------------------------------------!
  ! Coefficients of a family on a coupling into one vertex: the given
  ! sources into the given head, each edge labelled with its source
  ! degree and the degree its constraint determines.
  !-------------------------------------------------------------------!

  subroutine coefficients(scheme, num_vertices, tails, head, source_degree, determines, dt, c)

    class(family), intent(in)  :: scheme
    integer      , intent(in)  :: num_vertices, tails(:), head, source_degree(:), determines(:)
    real(dp)     , intent(in)  :: dt(:)
    real(dp), allocatable, intent(out) :: c(:)

    type(stored_directed_graph) :: coupling
    type(stored_field) :: steps, degrees, conditions
    class(field), allocatable :: out

    coupling = stored_directed_graph(num_vertices, tails=tails, heads=[(head, k = 1, size(tails))])

    steps      = stored_field('dt', coupling % vertex_set(), num_vertices)
    degrees    = stored_field('source degree', coupling % edge_set(), size(tails))
    conditions = stored_field('determines', coupling % edge_set(), size(tails))
    call steps      % set_real_vector(dt)
    call degrees    % set_integer_vector(source_degree)
    call conditions % set_integer_vector(determines)

    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(c)

  end subroutine coefficients

  subroutine bdf_on(label, dt)

    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: dt(:)

    integer, parameter :: last = 2 * order + 1
    real(dp), allocatable :: c(:)
    real(dp) :: h0, h1

    call coefficients(bdf_family(order), last, &
         & [(last - k, k = 0, order), (last - k, k = 0, 2 * order)], last, &
         & [(0, k = 0, order), (0, k = 0, 2 * order)], &
         & [(1, k = 0, order), (2, k = 0, 2 * order)], dt, c)

    write(*,'(a)') ' '
    write(*,'(a)') ' bdf 2 on a ' // label // ' grid'
    write(*,'(a,3f10.5)') '   velocity     alpha_0..2      ', c(1:order + 1)
    write(*,'(a,5f10.5)') '   acceleration beta_0..4       ', c(order + 2:)

    if (label == 'uniform') then
       write(*,'(a,3f10.5)') '   tabulated    alpha           ', [1.5_dp, -2.0_dp, 0.5_dp]
       write(*,'(a,5f10.5)') '   convolution  beta            ', [2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp]
    else
       h0 = dt(last)
       h1 = dt(last - 1)
       write(*,'(a,3f10.5)') '   set_bdf row  alpha           ', &
            & [(2.0_dp * h0 + h1) / (h0 + h1), -(h0 + h1) / h1, h0 * h0 / (h1 * (h0 + h1))]
    end if

  end subroutine bdf_on

  subroutine adams_on(label, p, dt)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: p
    real(dp)        , intent(in) :: dt(:)

    real(dp), allocatable :: c(:)

    call coefficients(adams_family(p), p, [(p - k, k = 0, p - 1)], p, &
         & [(2, k = 0, p - 1)], [(1, k = 0, p - 1)], dt, c)

    write(*,'(a)') ' '
    write(*,'(a)') ' adams-moulton 3 on a ' // label // ' grid'
    write(*,'(a,3f10.5)') '   quadrature   alpha_0..2      ', c
    if (label == 'uniform') then
       write(*,'(a,3f10.5)') '   tabulated    alpha           ', [5.0_dp, 8.0_dp, -1.0_dp] / 12.0_dp
    end if

  end subroutine adams_on

  !-------------------------------------------------------------------!
  ! Every stage-to-stage edge of a two-stage tableau, and the two
  ! recovery edges into the arriving instant.
  !-------------------------------------------------------------------!

  subroutine dirk_on(scheme)

    type(dirk_family), intent(in) :: scheme

    real(dp), allocatable :: c(:)
    integer :: s

    s = scheme % num_stages()

    call coefficients(scheme, 2 + s, [2, 2, 3], 3, [2, 2, 2], [1, 1, 1], &
         & [(0.5_dp, k = 1, 2 + s)], c)
    write(*,'(a)') ' '
    write(*,'(a)') ' crouzeix two-stage, stage 2 from stages 1, 1, 2'
    write(*,'(a,3f10.5)') '   a_21, a_21, a_22                ', c

    call coefficients(scheme, 2 + s, [2, 3], 2 + s, [2, 2], [2, 2], &
         & [(0.5_dp, k = 1, 2 + s)], c)
    write(*,'(a,2f10.5)') '   b_1, b_2 into the arriving instant', c
    write(*,'(a,3f10.5)') '   tableau gamma, 1 - 2 gamma, b   ', &
         & (3.0_dp + sqrt(3.0_dp)) / 6.0_dp, 1.0_dp - (3.0_dp + sqrt(3.0_dp)) / 3.0_dp, 0.5_dp

  end subroutine dirk_on

  !-------------------------------------------------------------------!
  ! The partial of bdf 2's velocity coefficients in the last step,
  ! from the family's own partial_action, beside a central finite
  ! difference of apply.
  !-------------------------------------------------------------------!

  subroutine bdf_step_sensitivity(dt)

    real(dp), intent(in) :: dt(:)

    integer , parameter :: last = 2 * order + 1
    real(dp), parameter :: delta = 1.0e-6_dp

    type(stored_directed_graph) :: coupling
    type(stored_field) :: steps, degrees, conditions, direction
    type(bdf_family) :: scheme
    class(field), allocatable :: out
    real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)

    scheme   = bdf_family(order)
    coupling = stored_directed_graph(last, tails=[(last - k, k = 0, order)], &
         & heads=[(last, k = 0, order)])

    steps      = stored_field('dt', coupling % vertex_set(), last)
    degrees    = stored_field('source degree', coupling % edge_set(), order + 1)
    conditions = stored_field('determines', coupling % edge_set(), order + 1)
    direction  = stored_field('v', coupling % vertex_set(), last)
    call degrees    % set_integer_vector([(0, k = 0, order)])
    call conditions % set_integer_vector([(1, k = 0, order)])

    allocate(v(last), source=0.0_dp)
    v(last) = 1.0_dp
    call direction % set_real_vector(v)

    call steps % set_real_vector(dt)
    call scheme % partial_action(coupling, [steps, degrees, conditions], &
         & [variation(scheme % argument(1), direction)], out)
    call out % real_vector(exact)

    call steps % set_real_vector(dt + delta * v)
    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(plus)
    call steps % set_real_vector(dt - delta * v)
    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(minus)

    write(*,'(a)') ' '
    write(*,'(a)') ' bdf 2, partial of alpha_0..2 in the last step, non-uniform grid'
    write(*,'(a,3f12.6)') '   partial_action, degree one    ', exact
    write(*,'(a,3f12.6)') '   central difference            ', (plus - minus) / (2.0_dp * delta)

    call bdf_second_partials(scheme, coupling, steps, degrees, conditions, dt, v)

  end subroutine bdf_step_sensitivity

  !-------------------------------------------------------------------!
  ! Degree two: the second partial in the last step beside a second
  ! central difference of apply, and the mixed partial in the last two
  ! steps beside a central difference of the degree-one partial.
  !-------------------------------------------------------------------!

  subroutine bdf_second_partials(scheme, coupling, steps, degrees, conditions, dt, v)

    type(bdf_family)           , intent(in)    :: scheme
    type(stored_directed_graph), intent(in)    :: coupling
    type(stored_field)         , intent(inout) :: steps
    type(stored_field)         , intent(in)    :: degrees, conditions
    real(dp)                   , intent(in)    :: dt(:), v(:)

    real(dp), parameter :: delta = 1.0e-4_dp
    type(stored_field) :: along_v, along_w
    class(field), allocatable :: out
    real(dp), allocatable :: plus(:), at(:), minus(:), w(:)
    real(dp), allocatable :: second(:), mixed_partial(:)

    along_v = stored_field('v', coupling % vertex_set(), size(dt))
    along_w = stored_field('w', coupling % vertex_set(), size(dt))
    call along_v % set_real_vector(v)
    w = 0.0_dp * v
    w(size(dt) - 1) = 1.0_dp
    call along_w % set_real_vector(w)

    call steps % set_real_vector(dt)
    call scheme % partial_action(coupling, [steps, degrees, conditions], &
         & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_v)], out)
    call out % real_vector(second)
    call scheme % partial_action(coupling, [steps, degrees, conditions], &
         & [variation(scheme % argument(1), along_v), variation(scheme % argument(1), along_w)], out)
    call out % real_vector(mixed_partial)

    call steps % set_real_vector(dt + delta * v)
    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(plus)
    call steps % set_real_vector(dt)
    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(at)
    call steps % set_real_vector(dt - delta * v)
    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(minus)

    write(*,'(a)') ' '
    write(*,'(a)') ' bdf 2, second partials of alpha_0..2, non-uniform grid'
    write(*,'(a,3f12.6)') '   partial_action, (last, last)  ', second
    write(*,'(a,3f12.6)') '   second central difference     ', (plus - 2.0_dp * at + minus) / delta**2

    call steps % set_real_vector(dt + delta * w)
    call scheme % partial_action(coupling, [steps, degrees, conditions], &
         & [variation(scheme % argument(1), along_v)], out)
    call out % real_vector(plus)
    call steps % set_real_vector(dt - delta * w)
    call scheme % partial_action(coupling, [steps, degrees, conditions], &
         & [variation(scheme % argument(1), along_v)], out)
    call out % real_vector(minus)

    write(*,'(a,3f12.6)') '   partial_action, (last, before)', mixed_partial
    write(*,'(a,3f12.6)') '   difference of degree one      ', (plus - minus) / (2.0_dp * delta)

  end subroutine bdf_second_partials

end program family_coefficients
