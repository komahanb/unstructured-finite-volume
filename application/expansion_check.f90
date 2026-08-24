! The functional and its derivatives in the design, to fourth order,
! each checked against a difference of the one below it.
!
! The expansion is exact given the orders beneath it, so the m-th
! derivative and a central difference of the (m-1)-th must agree to
! the accuracy of the difference and no better. Running the whole
! expansion at three designs gives every check at once, and each one
! is independent: the m-th derivative is computed from the statement,
! never from the (m-1)-th.
program expansion_check

  use iso_fortran_env      , only : dp => REAL64
  use operation_family     , only : family
  use operation_family_bdf , only : bdf_family
  use operation_family_adams, only : adams_family
  use physics_vanderpol    , only : van_der_pol, van_der_pol_energy
  use gti_march            , only : partition
  use gti_taylor           , only : block_expansion

  implicit none

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  integer , parameter :: instants = 21
  integer , parameter :: max_order = 4
  real(dp), parameter :: duration = 2.0_dp
  real(dp), parameter :: design = 1.0_dp
  real(dp), parameter :: delta = 1.0e-4_dp

  call expansion_of('bdf 2', bdf_family(2))
  call expansion_of('adams-moulton 3', adams_family(3))

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

  subroutine expanded(scheme, design_value, f)

    class(family), intent(in) :: scheme
    real(dp)     , intent(in) :: design_value
    real(dp), allocatable, intent(out) :: f(:)

    real(dp), allocatable :: dt(:), t(:), held(:), q(:)
    real(dp) :: achieved
    integer :: k, d

    call partition(duration, instants, dt, t)
    held = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, scheme % history_depth())]

    call block_expansion(scheme, van_der_pol(state_degree), &
         & van_der_pol_energy(state_degree), degrees, instants, dt, held, &
         & design_value, max_order, q, f, achieved)

  end subroutine expanded

  subroutine expansion_of(title, scheme)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme

    real(dp), allocatable :: f(:), plus(:), minus(:)
    real(dp) :: differenced(max_order)
    integer :: m

    call expanded(scheme, design, f)
    call expanded(scheme, design + delta, plus)
    call expanded(scheme, design - delta, minus)

    do m = 1, max_order
       differenced(m) = (plus(m - 1) - minus(m - 1)) / (2.0_dp * delta)
    end do

    write(*,'(a)')          ' '
    write(*,'(a)')          ' ' // title // ', van der pol at a design of one'
    write(*,'(a,5i15)')     '   derivative order    ', [(m, m = 0, max_order)]
    write(*,'(a,5f15.8)')   '   from the expansion  ', f
    write(*,'(a,15x,4f15.8)') '   differenced       ', differenced
    write(*,'(a,15x,4es15.2)') '   apart             ', abs(f(1:max_order) - differenced)

  end subroutine expansion_of

end program expansion_check
