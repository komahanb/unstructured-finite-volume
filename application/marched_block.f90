! One block solved.
!
! The rows a scheme makes, the row the physics makes, and the rows of
! the instants carried in are added into one statement, and newton
! drives it to zero. The tangent it linearizes with is exact: the
! stencil is its own jacobian and the physics differentiates its own
! rule, so nothing anywhere in the chain is differenced.
!
! The check is exact rather than comparative. Van der Pol at a design
! of zero is
!
!      q" - 0 (1 - q^2) q' + q  =  q" + q  =  0 ,
!
! the harmonic oscillator, whose solution through q(0) = 1, q'(0) = 0
! is the cosine. So the marched value is compared against cos(t)
! directly, and the steps are then halved: the error must fall by two
! raised to the scheme's order, which is what the printed ratio is.
!
! The last part marches the true oscillator, at a design of one,
! where the statement is nonlinear and newton has to iterate.
program marched_block

  use util_precision  , only : dp
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use physics_vanderpol     , only : van_der_pol
  use gti_block             , only : block_residual
  use gti_march             , only : partition, block_of, solved

  implicit none

  integer , parameter :: max_state_degree = 2
  integer , parameter :: degrees = max_state_degree + 1
  real(dp), parameter :: duration = 2.0_dp

  call order_of('bdf 2',           bdf_family(2),   2.0_dp)
  call order_of('bdf 3',           bdf_family(3),   3.0_dp)
  call order_of('adams-moulton 3', adams_family(3), 3.0_dp)
  call nonlinear()

contains

  pure integer function unknown(instant, degree)

    integer, intent(in) :: instant, degree

    unknown = (instant - 1) * degrees + degree + 1

  end function unknown



  !-------------------------------------------------------------------!
  ! March a block: the steps, the rows, the carried instants, and
  ! newton over all of them.
  !-------------------------------------------------------------------!

  subroutine march(scheme, n, design_value, q, t, achieved)

    class(family), intent(in)  :: scheme
    integer      , intent(in)  :: n
    real(dp)     , intent(in)  :: design_value
    real(dp), allocatable, intent(out) :: q(:), t(:)
    real(dp)     , intent(out) :: achieved

    type(block_residual) :: rows
    real(dp), allocatable :: dt(:), held(:)
    integer :: h, k, d

    call partition(duration, n, dt, t)
    h = scheme % history_depth(degrees - 1)

    held = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, h)]
    rows = block_of(scheme, van_der_pol(max_state_degree), degrees, n, dt, held)

    call solved(rows, design_value, q, achieved)

  end subroutine march

  !-------------------------------------------------------------------!
  ! The d-th derivative of the cosine.
  !-------------------------------------------------------------------!

  pure real(dp) function exact(d, t) result(q)

    integer , intent(in) :: d
    real(dp), intent(in) :: t

    select case (mod(d, 4))
    case (0)
       q =  cos(t)
    case (1)
       q = -sin(t)
    case (2)
       q = -cos(t)
    case default
       q =  sin(t)
    end select

  end function exact

  pure real(dp) function worst(q, t) result(e)

    real(dp), intent(in) :: q(:), t(:)

    integer :: k

    e = 0.0_dp
    do k = 1, size(t)
       e = max(e, abs(q(unknown(k, 0)) - exact(0, t(k))))
    end do

  end function worst

  !-------------------------------------------------------------------!
  ! The error at two resolutions, and the ratio between them.
  !-------------------------------------------------------------------!

  subroutine order_of(title, scheme, expected)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme
    real(dp)        , intent(in) :: expected

    real(dp), allocatable :: q(:), t(:)
    real(dp) :: e(4), achieved
    integer :: level, steps

    do level = 1, 4
       steps = 10 * 2 ** (level - 1)
       call march(scheme, steps + 1, 0.0_dp, q, t, achieved)
       e(level) = worst(q, t)
    end do

    write(*,'(a)')          ' '
    write(*,'(a)')          ' ' // title // ' on the harmonic oscillator'
    write(*,'(a,4i11)')     '   steps                      ', [(10 * 2 ** (level - 1), level = 1, 4)]
    write(*,'(a,4es11.3)')  '   worst error                ', e
    write(*,'(a,33x,3f11.3)') '   ratio                    ', e(1:3) / e(2:4)
    write(*,'(a,f11.3)')    '   two to the scheme order    ', 2.0_dp ** expected
    write(*,'(a,es11.3)')   '   residual newton achieved   ', achieved

  end subroutine order_of

  !-------------------------------------------------------------------!
  ! The true oscillator, where the statement is nonlinear.
  !-------------------------------------------------------------------!

  subroutine nonlinear()

    real(dp), allocatable :: q(:), t(:)
    real(dp) :: achieved

    call march(bdf_family(2), 41, 1.0_dp, q, t, achieved)

    write(*,'(a)')        ' '
    write(*,'(a)')        ' van der pol at a design of one, bdf 2, 40 steps'
    write(*,'(a,es11.3)') '   residual newton achieved   ', achieved
    write(*,'(a,3f11.5)') '   the last instant, q q'' q"  ', &
         & q(unknown(size(t), 0)), q(unknown(size(t), 1)), q(unknown(size(t), 2))

  end subroutine nonlinear

end program marched_block
