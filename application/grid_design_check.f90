! THE GRID AS A DESIGN. The steps are the weights of a designed grid,
! normalised to the duration, so a functional of the march is a
! function of those weights and of the physics' parameter. Its first
! derivatives come by the two routes over the whole table, one row
! per functional and one column per design, and are checked three
! ways:
!
!   the routes     tangent against adjoint over the table
!   homogeneity    the steps do not change when every weight is
!                  scaled, so sum_j p_j df/dp_j = 0 exactly
!   differences    a central difference of the march in a weight,
!                  and in the parameter, against the route: the
!                  march is known to a relative tolerance tau, so the
!                  step is tau^(1/3) and the agreement tau^(2/3)
!
! and then the derivatives of every order from two to the order asked
! for, by the reverse route, one table per order with one column per
! multiset of designs, checked three ways at every order: the entries
! of one multiset - one per distinct design of the multiset against the rest -
! agree in theory and are not made to; the entry of the parameter
! alone is the expansion's coefficient of that order; and a central
! difference of the table one order below, in a weight or in the
! parameter, gives the entries holding that design to tau^(2/3).
!
! The families read three instants back, so a startup block marches
! the first three from the given state, on steps split four ways;
! it is part of the chain, and its own dependence on the weights and
! the parameter is carried through the junction like any other.
!
!      ./grid_design_check [tolerance] [order]
program grid_design_check

  use util_precision        , only : dp
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : designed_grid
  use operation_minimization, only : relative, by_rate
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use gti_expansion         , only : family_holder, expansion
  use gti_march             , only : set_stopping, consistent_state, imbalance
  use gti_chain             , only : chain_block, march_chain, chain_expansion, chain_system, &
       & chain_systems, functional_holder, one_functional, &
       & chain_derivative, asymmetry, multiset_count, multiset_rank, multiset_of
  use gti_sweeps            , only : route_of, forward_route, reverse_route

  implicit none

  integer , parameter :: state_degree = 2, degrees = state_degree + 1
  integer , parameter :: instants = 21, checked(3) = [1, 7, 20]
  real(dp), parameter :: duration = 4.0_dp, design = 0.8_dp

  type(family_holder)     :: schemes(2)
  type(functional_holder) :: functionals(2)
  type(chain_block) , allocatable :: chain(:)
  type(expansion), allocatable, target :: tower
  type(chain_system), allocatable :: systems(:)
  real(dp), allocatable :: p(:), dt(:), t(:), v(:,:), f(:,:), tangent(:,:), adjoint(:,:)
  real(dp), allocatable :: plus(:,:), minus(:,:), q0(:), table(:,:), entries(:,:,:)
  real(dp), allocatable :: below(:,:), above(:,:), by_class(:)
  real(dp) :: tau, delta, achieved, worst
  character(len=32) :: argument
  integer :: k, j, i, route, order, max_order, nd

  tau       = 1.0e-12_dp
  max_order = 3
  call get_command_argument(1, argument)
  if (len_trim(argument) > 0) read(argument, *) tau
  call get_command_argument(2, argument)
  if (len_trim(argument) > 0) read(argument, *) max_order
  call set_stopping(tau, relative, by_rate, 100)
  delta = tau ** (1.0_dp / 3.0_dp)

  ! two blocks of three-history families over weights that are not
  ! uniform, so that no step is like another
  allocate(schemes(1) % scheme, source=bdf_family(3))
  allocate(schemes(2) % scheme, source=adams_family(3))
  functionals(1) = one_functional(van_der_pol_energy(state_degree))
  functionals(2) = one_functional(van_der_pol_dissipation(state_degree))
  p  = [(1.0_dp + 0.5_dp * sin(real(k, dp)), k = 1, instants - 1)]
  q0 = consistent_state(van_der_pol(state_degree), degrees, [1.0_dp, 0.0_dp], design)

  call marched(p, design, f)
  call tower % step_partials(v)
  call chain_systems(chain, tower, functionals, degrees, systems)
  call chain_derivative(chain, tower, systems, functionals, degrees, 1, forward_route, tangent)
  call chain_derivative(chain, tower, systems, functionals, degrees, 1, reverse_route, adjoint)
  route   = route_of(size(tangent, 2), size(tangent, 1), 1)

  write(*,'(a,es9.2,a,es9.2,a,es9.2)') ' relative tolerance', tau, '   difference step', delta, &
       & '   expected agreement tau^(2/3)', tau ** (2.0_dp / 3.0_dp)
  write(*,'(a,i0,a,i0,a,a)') ' designs ', size(tangent, 2), '   functionals ', size(tangent, 1), &
       & '   the gate chooses the ', trim(merge('forward', 'reverse', route == forward_route))
  write(*,'(a,es10.2)') ' tangent against adjoint over the table, relative  ', &
       & maxval(abs(tangent - adjoint)) / maxval(abs(tangent))
  write(*,'(a,es10.2,a,es10.2)') ' physics column against the expansion, relative   ', &
       & abs(adjoint(1, 1) - f(1, 1)) / abs(f(1, 1)), '  ', abs(adjoint(2, 1) - f(1, 2)) / abs(f(1, 2))
  do i = 1, 2
     write(*,'(a,i0,a,es10.2,a,es10.2,a,es10.2)') ' homogeneity, functional ', i, &
          & ':  p . df/dp / |p||df/dp|  ', &
          & dot_product(p, adjoint(i, 2:)) / (norm2(p) * norm2(adjoint(i, 2:))), &
          & '   f ', f(0, i), '   |df/dp| ', norm2(adjoint(i, 2:))
  end do

  ! central differences in three weights and in the parameter
  worst = 0.0_dp
  do k = 1, size(checked)
     j = checked(k)
     call marched(p + delta * unit(j), design, plus)
     call marched(p - delta * unit(j), design, minus)
     do i = 1, 2
        worst = max(worst, abs((plus(0, i) - minus(0, i)) / (2.0_dp * delta) - adjoint(i, 1 + j)) &
             & / max(1.0_dp, abs(adjoint(i, 1 + j))))
     end do
  end do
  write(*,'(a,es10.2)') ' differenced in three weights against the route, worst   ', worst
  call marched(p, design + delta, plus)
  call marched(p, design - delta, minus)
  write(*,'(a,es10.2)') ' differenced in the parameter against the route, worst   ', &
       & maxval(abs((plus(0, :) - minus(0, :)) / (2.0_dp * delta) - adjoint(:, 1)) / &
       &        max(1.0_dp, abs(adjoint(:, 1))))

  nd = size(tangent, 2)

  ! every order above one by the reverse route
  do order = 2, max_order
     call marched(p, design, f, order)
     call chain_systems(chain, tower, functionals, degrees, systems)
     call chain_derivative(chain, tower, systems, functionals, degrees, order, reverse_route, &
          & table, entries=entries)
     write(*,'(a)') ' '
     write(*,'(a,i0,a,i0,a,i0,a,i0)') ' derivatives of order ', order, ' by the reverse route: ', &
          & size(table, 1), ' tables of ', size(table, 2), ' multisets over ', nd
     do i = 1, 2
        write(*,'(a,i0,a,es10.2,a,es10.2)') ' functional ', i, &
             & ':  departure among the entries of a multiset, relative ', &
             & asymmetry(entries(i:i, :, :), nd, order), &
             & '   parameter entry against the expansion ', &
             & abs(table(i, 1) - f(order, i)) / abs(f(order, i))
     end do
     ! central differences of the table one order below, in three
     ! weights and in the parameter, against the entries holding that
     ! design, the departure by the count of parameter indices among
     ! the order's
     allocate(by_class(0:order), source=0.0_dp)
     do k = 1, size(checked)
        j = checked(k)
        call differenced_table(p + delta * unit(j), design, order - 1, above)
        call differenced_table(p - delta * unit(j), design, order - 1, below)
        call classed((above - below) / (2.0_dp * delta), 1 + j)
     end do
     call differenced_table(p, design + delta, order - 1, above)
     call differenced_table(p, design - delta, order - 1, below)
     call classed((above - below) / (2.0_dp * delta), 1)
     write(*,'(a,i0,a,*(es10.2))') ' differenced tables of order ', order - 1, &
          & ' against the entries, worst by parameter count from ', &
          & by_class(order:0:-1) / max(1.0_dp, maxval(abs(table)))
     deallocate(by_class)
  end do

contains

  subroutine marched(weights, nu, f, order)

    real(dp), intent(in) :: weights(:), nu
    real(dp), allocatable, intent(out) :: f(:,:)
    integer , intent(in), optional :: order

    type(imbalance) :: left
    integer :: m

    m = 1
    if (present(order)) m = order
    call march_chain(schemes, [11, 10], van_der_pol(state_degree), degrees, &
         & designed_grid(duration), nu, q0, chain, tower, dt, t, achieved, grid_design=weights, &
         & left=left, startup=4)
    if (.not. left % converged) error stop 'grid_design_check: the march converged'
    call chain_expansion(chain, tower, functionals, degrees, m, f)

  end subroutine marched

  ! the table of one order at other weights or parameter, by the
  ! reverse route
  subroutine differenced_table(weights, nu, order, t)

    real(dp), intent(in) :: weights(:), nu
    integer , intent(in) :: order
    real(dp), allocatable, intent(out) :: t(:,:)

    real(dp), allocatable :: f(:,:)

    call marched(weights, nu, f)
    call chain_systems(chain, tower, functionals, degrees, systems)
    call chain_derivative(chain, tower, systems, functionals, degrees, order, reverse_route, t)

  end subroutine differenced_table

  ! the departure of a differenced table in design l from the entries
  ! holding l, by the count of the parameter among the order's designs
  subroutine classed(e, l)

    real(dp), intent(in) :: e(:,:)
    integer , intent(in) :: l

    integer, allocatable :: s(:), with(:)
    integer :: rank, c

    do rank = 1, size(e, 2)
       s    = multiset_of(rank, order - 1, nd)
       with = [s(1:count(s < l)), l, s(count(s < l) + 1:)]
       c    = count(with == 1)
       by_class(c) = max(by_class(c), maxval(abs(e(:, rank) - table(:, multiset_rank(with, nd)))))
    end do

  end subroutine classed


  pure function unit(j) result(e)

    integer, intent(in) :: j
    real(dp) :: e(instants - 1)

    e    = 0.0_dp
    e(j) = 1.0_dp

  end function unit

end program grid_design_check
