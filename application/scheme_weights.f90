! The two fields a coupling carries and their product: the step
! scaling tau, the family's dimensionless alpha, and the weight
! tau * alpha that a constraint row holds for each source.
!
! The weights are then applied to sampled polynomials. A row that
! determines degree d and reads sources of degree d' must satisfy
!
!      sum over edges of  weight * q^(d')(t_tail)  =  q^(d)(t_head)
!
! exactly for every polynomial the scheme reproduces. The residual is
! printed across a span of polynomial degrees so that the degree at
! which exactness ends is shown rather than asserted.
!
! Both families are exercised because they scale in opposite
! directions: a backward difference reads values into a derivative
! and its exponents are negative, while an Adams quadrature reads a
! derivative into a lower one and its exponents are zero and
! positive. One family alone would not catch the exponent written
! the wrong way round.
program scheme_weights

  use util_precision  , only : dp
  use view_directed_stored   , only : stored_directed_graph
  use field_calculus         , only : field
  use field_stored           , only : stored_field
  use operation_action       , only : variation
  use operation_family       , only : family
  use operation_family_bdf   , only : bdf_family
  use operation_family_adams , only : adams_family
  use operation_coupling     , only : coupling_inputs
  use operation_weight       , only : scheme_weight

  implicit none

  integer :: k

  call bdf_rows(2, 'uniform',     [0.0_dp, (0.5_dp, k = 2, 5)])
  call bdf_rows(2, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp])
  call adams_row(3, 'uniform',     [0.0_dp, 0.5_dp, 0.5_dp])
  call adams_row(3, 'non-uniform', [0.0_dp, 0.30_dp, 0.20_dp])
  call weight_partials(2, [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp])

contains

  !-------------------------------------------------------------------!
  ! The instants a step field reaches: the first carries no step.
  !-------------------------------------------------------------------!

  pure function instants(dt) result(t)

    real(dp), intent(in) :: dt(:)
    real(dp) :: t(size(dt))

    integer :: i

    t(1) = 0.0_dp
    do i = 2, size(dt)
       t(i) = t(i - 1) + dt(i)
    end do

  end function instants

  !-------------------------------------------------------------------!
  ! The d-th derivative of t raised to m, which is zero once d passes
  ! m.
  !-------------------------------------------------------------------!

  pure real(dp) function power_derivative(m, d, t) result(q)

    integer , intent(in) :: m, d
    real(dp), intent(in) :: t

    integer :: i

    if (d > m) then
       q = 0.0_dp
       return
    end if

    q = 1.0_dp
    do i = 0, d - 1
       q = q * real(m - i, dp)
    end do
    q = q * t**(m - d)

  end function power_derivative

  !-------------------------------------------------------------------!
  ! The three fields on one row: tau, alpha and their product.
  !-------------------------------------------------------------------!

  subroutine row_fields(scheme, nv, tails, head, source_degree, determines, dt, &
       & tau, alpha, w)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: nv, tails(:), head, source_degree(:), determines(:)
    real(dp)     , intent(in) :: dt(:)
    real(dp), allocatable, intent(out) :: tau(:), alpha(:), w(:)

    type(stored_directed_graph) :: coupling
    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: out
    type(scheme_weight) :: weights
    integer :: e

    call coupling_inputs(nv, tails, [(head, e = 1, size(tails))], dt, source_degree, determines, &
         & coupling, inputs)

    tau = [(step_power(dt(head), source_degree(e) - determines(e)), e = 1, size(tails))]

    call scheme % apply(coupling, inputs, out)
    call out % real_vector(alpha)

    weights = scheme_weight(scheme)
    call weights % apply(coupling, inputs, out)
    call out % real_vector(w)

  end subroutine row_fields

  !-------------------------------------------------------------------!
  ! The step raised to an integer exponent of either sign, as
  ! repeated products, matching the exact arithmetic's power.
  !-------------------------------------------------------------------!

  pure real(dp) function step_power(h, n) result(p)

    real(dp), intent(in) :: h
    integer , intent(in) :: n

    integer :: i

    p = 1.0_dp
    do i = 1, abs(n)
       p = p * h
    end do
    if (n < 0) p = 1.0_dp / p

  end function step_power

  !-------------------------------------------------------------------!
  ! The residual of one row on the polynomial t raised to m.
  !-------------------------------------------------------------------!

  pure real(dp) function row_residual(w, tails, head, source_degree, determines, t, m) &
       & result(r)

    real(dp), intent(in) :: w(:), t(:)
    integer , intent(in) :: tails(:), head, source_degree(:), determines(:), m

    integer :: e

    r = -power_derivative(m, determines(1), t(head))
    do e = 1, size(w)
       r = r + w(e) * power_derivative(m, source_degree(e), t(tails(e)))
    end do

  end function row_residual

  !-------------------------------------------------------------------!
  ! One row printed, then its residual across polynomial degrees.
  !-------------------------------------------------------------------!

  subroutine one_row(title, scheme, nv, tails, head, source_degree, determines, dt, top)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme
    integer         , intent(in) :: nv, tails(:), head, source_degree(:), determines(:), top
    real(dp)        , intent(in) :: dt(:)

    real(dp), allocatable :: tau(:), alpha(:), w(:)
    real(dp) :: t(nv), residual(0:top)
    integer :: m

    call row_fields(scheme, nv, tails, head, source_degree, determines, dt, tau, alpha, w)
    t = instants(dt)

    write(*,'(a)') ' '
    write(*,'(a)')        ' ' // title
    write(*,'(a,9f11.5)') '   tau                        ', tau
    write(*,'(a,9f11.5)') '   alpha                      ', alpha
    write(*,'(a,9f11.5)') '   weight                     ', w

    do m = 0, top
       residual(m) = row_residual(w, tails, head, source_degree, determines, t, m)
    end do

    write(*,'(a,9i11)')     '   on t**m, m =              ', [(m, m = 0, top)]
    write(*,'(a,9es11.2)')  '   residual                  ', residual

  end subroutine one_row

  subroutine bdf_rows(p, label, dt)

    integer         , intent(in) :: p
    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: dt(:)

    integer :: last

    last = 2 * p + 1

    call one_row('bdf ' // digit(p) // ' velocity row, ' // label // ' grid', &
         & bdf_family(p), last, [(last - k, k = 0, p)], last, &
         & [(0, k = 0, p)], [(1, k = 0, p)], dt, p + 2)

    call one_row('bdf ' // digit(p) // ' acceleration row, ' // label // ' grid', &
         & bdf_family(p), last, [(last - k, k = 0, 2 * p)], last, &
         & [(0, k = 0, 2 * p)], [(2, k = 0, 2 * p)], dt, p + 2)

  end subroutine bdf_rows

  !-------------------------------------------------------------------!
  ! The Adams velocity row: the velocity one instant back carried
  ! unchanged, and the accelerations quadratured over the last step.
  !-------------------------------------------------------------------!

  subroutine adams_row(p, label, dt)

    integer         , intent(in) :: p
    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: dt(:)

    call one_row('adams-moulton ' // digit(p) // ' velocity row, ' // label // ' grid', &
         & adams_family(p), p, [p - 1, (p - k, k = 0, p - 1)], p, &
         & [1, (2, k = 0, p - 1)], [(1, k = 0, p)], dt, p + 2)

  end subroutine adams_row

  !-------------------------------------------------------------------!
  ! The partial of the weights in the last step. The product rule
  ! must carry it through both factors: the coefficient's own partial
  ! times the scaling, plus the coefficient times the scaling's. A
  ! partial that came from one factor alone would not match.
  !-------------------------------------------------------------------!

  subroutine weight_partials(p, dt)

    integer , intent(in) :: p
    real(dp), intent(in) :: dt(:)

    real(dp), parameter :: delta = 1.0e-6_dp

    type(stored_directed_graph) :: coupling
    type(stored_field), allocatable :: inputs(:)
    type(stored_field) :: direction
    type(scheme_weight) :: weights
    class(field), allocatable :: out
    real(dp), allocatable :: exact(:), plus(:), minus(:), v(:)
    integer , allocatable :: tails(:)
    integer :: last, e, j

    last  = 2 * p + 1
    tails = [(last - j, j = 0, p)]

    call coupling_inputs(last, tails, [(last, e = 1, size(tails))], dt, &
         & [(0, e = 1, size(tails))], [(1, e = 1, size(tails))], coupling, inputs)

    direction = stored_field('v', coupling % vertex_set(), last)
    allocate(v(last), source=0.0_dp)
    v(last) = 1.0_dp
    call direction % set_real_vector(v)

    weights = scheme_weight(bdf_family(p))

    call weights % partial_action(coupling, inputs, [variation(weights % argument(1), direction)], out)
    call out % real_vector(exact)

    call inputs(1) % set_real_vector(dt + delta * v)
    call weights % apply(coupling, inputs, out)
    call out % real_vector(plus)
    call inputs(1) % set_real_vector(dt - delta * v)
    call weights % apply(coupling, inputs, out)
    call out % real_vector(minus)

    write(*,'(a)') ' '
    write(*,'(a)') ' bdf 2 velocity weights, partial in the last step'
    write(*,'(a,3f13.5)') '   partial_action             ', exact
    write(*,'(a,3f13.5)') '   central difference         ', (plus - minus) / (2.0_dp * delta)

  end subroutine weight_partials

  pure function digit(n) result(c)

    integer, intent(in) :: n
    character(len=1) :: c

    write(c,'(i1)') n

  end function digit

end program scheme_weights
