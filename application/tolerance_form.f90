! The analytical form of the weight a block carries, identified from
! the family's own coefficients rather than measured from its matrix.
!
! A family's coefficients are dimensionless; the weight an edge
! finally carries is the coefficient times dt raised to the difference
! of the two degrees it joins. For a difference family every source is
! a value and a row determining the d-th derivative therefore carries
! dt^-d.
!
! The d-th row is the velocity row composed d times, so its
! coefficients are bounded by the d-th power of the velocity row's,
! and this reports whether that bound is attained. If it is, then
!
!      ||A||_inf  =  1  +  (sum_j |a_j|)^d / dt^d
!
! with a the velocity coefficients and the one being the unit a row
! carries on the column it determines. That is a closed form for the
! norm, hence for the floor eps ||A|| ||q|| a march can reach, and
! hence for a tolerance that need not be chosen.
!
! Both sides are computed here: the left from the assembled jacobian,
! the right from the family through the same coupling the scheme is
! built on.
program tolerance_form

  use iso_fortran_env       , only : dp => REAL64
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_grid        , only : uniform_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use gti_march             , only : partition
  use gti_chain             , only : chain_block, march_chain, &
       & chain_system, chain_systems

  implicit none

  write(*,'(a)') ' '
  write(*,'(a)') '  the velocity row, and the d-th row it composes to'
  write(*,'(a)') '  scheme    d    sum|a|    sum|c(d)|   (sum|a|)^d   attained'
  call composed(1, 2)
  call composed(2, 2)
  call composed(3, 2)
  call composed(4, 2)
  call composed(2, 3)
  call composed(3, 3)

  write(*,'(a)') ' '
  write(*,'(a)') '  ||A||_inf, from the family against the assembled jacobian'
  write(*,'(a)') '  scheme    d   instants       dt      from family     assembled   agreement'
  call against(1, 3, 21)
  call against(2, 3, 21)
  call against(3, 3, 21)
  call against(2, 3, 61)
  call against(2, 4, 61)
  call against(3, 3, 41)

contains

  !-------------------------------------------------------------------!
  ! The absolute sum of a family's coefficients on the row determining
  ! the given degree, from the coupling that row reads.
  !-------------------------------------------------------------------!

  real(dp) function row_sum(scheme, order, determines) result(total)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: order, determines

    type(stored_directed_graph) :: coupling
    type(stored_field) :: steps, degrees, conditions
    class(field), allocatable :: out
    real(dp), allocatable :: c(:)
    integer :: reach, last, k

    reach = determines * order
    last  = reach + 1

    coupling = stored_directed_graph(last, &
         & tails=[(last - k, k = 0, reach)], heads=[(last, k = 0, reach)])

    steps      = stored_field('dt', coupling % vertex_set(), last)
    degrees    = stored_field('source degree', coupling % edge_set(), reach + 1)
    conditions = stored_field('determines', coupling % edge_set(), reach + 1)

    call steps      % set_real_vector([(1.0_dp, k = 1, last)])
    call degrees    % set_integer_vector([(0, k = 0, reach)])
    call conditions % set_integer_vector([(determines, k = 0, reach)])

    call scheme % apply(coupling, [steps, degrees, conditions], out)
    call out % real_vector(c)

    total = sum(abs(c))

  end function row_sum

  subroutine composed(order, determines)

    integer, intent(in) :: order, determines

    real(dp) :: velocity, derived, powered
    character(len=8) :: named

    velocity = row_sum(bdf_family(order), order, 1)
    derived  = row_sum(bdf_family(order), order, determines)
    powered  = velocity ** determines

    write(named,'(a,i0)') 'bdf ', order
    write(*,'(a,a,i5,3f12.4,a)') '  ', named, determines, velocity, derived, powered, &
         & merge('   yes', '    no', abs(derived - powered) <= 1.0e-10_dp * powered)

  end subroutine composed

  !-------------------------------------------------------------------!
  ! The closed form beside the assembled jacobian's infinity norm.
  !-------------------------------------------------------------------!

  subroutine against(order, degrees, instants)

    integer, intent(in) :: order, degrees, instants

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    type(chain_system), allocatable :: systems(:)
    type(bdf_family) :: scheme
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:)
    real(dp) :: achieved, duration, design, predicted, assembled
    character(len=8) :: named
    integer :: k, d, top

    duration = 3.0_dp
    design   = 1.0_dp
    scheme   = bdf_family(order)
    top      = degrees - 1

    allocate(schemes(1))
    allocate(schemes(1) % scheme, source=scheme)
    added = [instants]

    call partition(duration, instants, dt, t)
    held = [((0.0_dp, d = 0, degrees - 1), k = 1, &
         &   scheme % history_depth(degrees - 1))]

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)
    call chain_systems(chain, van_der_pol_energy(degrees - 1), degrees, dt, &
         & design, systems)

    predicted = 1.0_dp + row_sum(scheme, order, top) / dt(size(dt)) ** top
    assembled = largest_row(systems(1) % a)

    write(named,'(a,i0)') 'bdf ', order
    write(*,'(a,a,i5,i10,f10.5,2es15.5,f11.3,a)') '  ', named, top, instants, &
         & dt(size(dt)), predicted, assembled, &
         & 100.0_dp * (1.0_dp - abs(predicted - assembled) / assembled), ' %'

  end subroutine against

  real(dp) function largest_row(a) result(most)

    real(dp), intent(in) :: a(:,:)

    real(dp) :: total
    integer  :: i, j

    most = 0.0_dp
    do i = 1, size(a, 1)
       total = 0.0_dp
       do j = 1, size(a, 2)
          total = total + abs(a(i, j))
       end do
       most = max(most, total)
    end do

  end function largest_row

end program tolerance_form
