! A block marched by stages.
!
! The physics is enforced at the stages and the instants are
! recovered from them, so a stage block's unknowns are not one set
! per instant: each step carries its stages and then the instant it
! arrives at, and the block's first instant carries nothing but
! itself.
!
! Checked the same way a multistep block was: van der Pol at a design
! of zero is the harmonic oscillator, whose solution is the cosine,
! so the marched value is compared against it and the steps halved.
program marched_stages

  use util_precision  , only : dp
  use operation_family     , only : family
  use operation_family_dirk, only : dirk_family, implicit_midpoint, &
       & crouzeix_two_stage, crouzeix_three_stage
  use physics_vanderpol    , only : van_der_pol
  use gti_block            , only : block_residual
  use gti_march            , only : partition, solved
  use gti_stage            , only : stage_block_of, instant_at

  implicit none

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  real(dp), parameter :: duration = 2.0_dp

  call order_of('implicit midpoint',  implicit_midpoint())
  call order_of('crouzeix two-stage', crouzeix_two_stage())
  call order_of('crouzeix three-stage', crouzeix_three_stage())

contains

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

  !-------------------------------------------------------------------!
  ! One stage block marched, and the worst error in the value at the
  ! instants it recovered.
  !-------------------------------------------------------------------!

  real(dp) function worst_error(scheme, n) result(e)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: n

    type(block_residual) :: rows
    real(dp), allocatable :: dt(:), t(:), q(:)
    real(dp) :: achieved
    integer :: k, d, s

    call partition(duration, n, dt, t)
    s = scheme % num_stages()

    rows = stage_block_of(scheme, van_der_pol(state_degree), degrees, n, dt, &
         & [(exact(d, t(1)), d = 0, degrees - 1)])
    call solved(rows, 0.0_dp, q, achieved)

    e = 0.0_dp
    do k = 1, n
       e = max(e, abs(q(instant_at(k, s, degrees) + 1) - exact(0, t(k))))
    end do

  end function worst_error

  subroutine order_of(title, scheme)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme

    real(dp) :: e(3)
    integer :: level

    do level = 1, 3
       e(level) = worst_error(scheme, 5 * 2 ** (level - 1) + 1)
    end do

    write(*,'(a)')            ' '
    write(*,'(a)')            ' ' // title // ' on the harmonic oscillator'
    write(*,'(a,i0,a)')       '   stages                     ', scheme % num_stages(), ''
    write(*,'(a,3i11)')       '   steps                      ', [(5 * 2 ** (level - 1), level = 1, 3)]
    write(*,'(a,3es11.3)')    '   worst error                ', e
    write(*,'(a,22x,2f11.3)') '   ratio                    ', e(1:2) / e(2:3)

  end subroutine order_of

end program marched_stages
