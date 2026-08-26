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
! and then the second derivatives by the reverse route at order two,
! one hessian per functional, checked three ways: the table is
! symmetric in theory and is not made so; its parameter entry is the
! expansion's second order; and a central difference of the first
! derivatives, themselves known to tau, gives the rows to tau^(2/3).
!
! The families read three instants back, so a startup block marches
! the first three from the given state, on steps split four ways;
! it is part of the chain, and its own dependence on the weights and
! the parameter is carried through the junction like any other.
!
!      ./grid_design_check [tolerance]
program grid_design_check

  use util_precision        , only : dp
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : designed_grid
  use operation_minimization, only : relative, by_rate
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy, van_der_pol_dissipation
  use gti_expansion         , only : family_holder
  use gti_march             , only : set_stopping, consistent_state, step_partials, imbalance
  use gti_chain             , only : chain_block, march_chain, chain_expansion, chain_system, &
       & chain_systems, chain_by_tangent, chain_by_adjoint, functional_holder, one_functional, &
       & chain_hessian
  use gti_sweeps            , only : route_of, forward_route

  implicit none

  integer , parameter :: state_degree = 2, degrees = state_degree + 1
  integer , parameter :: instants = 21, checked(3) = [1, 7, 20]
  real(dp), parameter :: duration = 4.0_dp, design = 0.8_dp

  type(family_holder)     :: schemes(2)
  type(functional_holder) :: functionals(2)
  type(chain_block) , allocatable :: chain(:)
  type(chain_system), allocatable :: systems(:)
  real(dp), allocatable :: p(:), dt(:), t(:), v(:,:), f(:,:), tangent(:,:), adjoint(:,:)
  real(dp), allocatable :: plus(:,:), minus(:,:), q0(:), hessian(:,:,:), dplus(:,:), dminus(:,:)
  real(dp) :: tau, delta, achieved, worst
  character(len=32) :: argument
  integer :: k, j, i, route

  tau = 1.0e-12_dp
  call get_command_argument(1, argument)
  if (len_trim(argument) > 0) read(argument, *) tau
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
  call step_partials(designed_grid(duration), instants, p, v)
  call chain_systems(chain, functionals, degrees, design, systems, step_partials=v)
  tangent = chain_by_tangent(chain, systems, degrees, design)
  adjoint = chain_by_adjoint(chain, systems, degrees, design)
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

  ! the second derivatives: the hessian of every functional, by the
  ! reverse route at order two
  call marched(p, design, f, 2)
  call chain_hessian(chain, systems, functionals, degrees, design, hessian, &
       & step_partials=v, steps=designed_grid(duration), grid_design=p)
  write(*,'(a)') ' '
  write(*,'(a,i0,a,i0,a,i0)') ' second derivatives by the reverse route: hessians ', &
       & size(hessian, 1), ' of ', size(hessian, 2), ' x ', size(hessian, 3)
  do i = 1, 2
     write(*,'(a,i0,a,es10.2,a,es10.2)') ' functional ', i, &
          & ':  symmetry, relative ', maxval(abs(hessian(i, :, :) - transpose(hessian(i, :, :)))) &
          & / maxval(abs(hessian(i, :, :))), &
          & '   parameter entry against the expansion ', abs(hessian(i, 1, 1) - f(2, i)) / abs(f(2, i))
  end do
  worst = 0.0_dp
  do k = 1, size(checked)
     j = checked(k)
     call differenced(p + delta * unit(j), design, dplus)
     call differenced(p - delta * unit(j), design, dminus)
     worst = max(worst, maxval(abs((dplus - dminus) / (2.0_dp * delta) - hessian(:, 1 + j, :)) &
          & / max(1.0_dp, maxval(abs(hessian)))))
  end do
  call differenced(p, design + delta, dplus)
  call differenced(p, design - delta, dminus)
  worst = max(worst, maxval(abs((dplus - dminus) / (2.0_dp * delta) - hessian(:, 1, :)) &
       & / max(1.0_dp, maxval(abs(hessian)))))
  write(*,'(a,es10.2)') ' differenced first derivatives against the hessian rows, worst ', worst

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
         & designed_grid(duration), nu, q0, chain, dt, t, achieved, grid_design=weights, &
         & left=left, startup=4)
    if (.not. left % converged) error stop 'grid_design_check: the march converged'
    call chain_expansion(chain, van_der_pol(state_degree), functionals, degrees, nu, m, f)

  end subroutine marched

  ! the first derivatives by the adjoint at other weights or parameter
  subroutine differenced(weights, nu, df)

    real(dp), intent(in) :: weights(:), nu
    real(dp), allocatable, intent(out) :: df(:,:)

    real(dp), allocatable :: f(:,:), vw(:,:)

    call marched(weights, nu, f)
    call step_partials(designed_grid(duration), instants, weights, vw)
    call chain_systems(chain, functionals, degrees, nu, systems, step_partials=vw)
    df = chain_by_adjoint(chain, systems, degrees, nu)

  end subroutine differenced

  pure function unit(j) result(e)

    integer, intent(in) :: j
    real(dp) :: e(instants - 1)

    e    = 0.0_dp
    e(j) = 1.0_dp

  end function unit

end program grid_design_check
