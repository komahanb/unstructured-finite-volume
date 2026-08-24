!=====================================================================!
! The functional and every derivative of it in the design, to any
! order.
!
! The statement R(q(x), x) = 0 determines the trajectory from the
! design, so the trajectory has a Taylor expansion in the design and
! its coefficients come out one order at a time. Differentiating the
! statement m times and collecting,
!
!      J q^(m)  +  r_m  =  0 ,
!
! where r_m is the m-th coefficient of the statement computed with
! q^(m) held at zero and every lower coefficient already known. So
! each order is one solve against the same jacobian the march already
! formed, and the whole expansion costs one factorisation and m right
! hand sides.
!
! Only the governing rows contribute to r_m at all. The scheme's rows
! are linear in the state, so their m-th coefficient is the jacobian
! times q^(m), which is held at zero. The carried rows hold given
! numbers, so theirs is zero as well - and a carried instant has a
! governing row like any other, so that row has to be struck out
! after the physics is placed, or the instants a block was given
! would appear to move with the design when they do not. Every term
! that remains is the physics, and the physics is nodal, so r_m is
! computed one instant at a time.
!
!             WHERE A COEFFICIENT IS TAKEN
!
! At the points the rule is stated on, which are given rather than
! assumed: a multistep block's physics is evaluated at its instants
! and a stage block's at its stages, while the functional is
! evaluated at the instants either way, so that one number means the
! same thing whichever family produced it.
!
!             HOW A COEFFICIENT IS TAKEN
!
! A quantity is handed to the physics as derivative_terms seeded
! symmetrically: the subsets of size k hold the k-th derivative of
! that quantity in the design. The arithmetic then composes them the
! way Leibniz does, and the coefficient of the full subset is the
! m-th derivative of whatever the rule built. The design itself is
! seeded with its value and a one, since it is the parameter.
!
! Nothing is differenced and no order is approximated: each
! coefficient is exact given the ones below it.
!
!             WHAT IS REFUSED
!
! A negative order, and a state series whose extent does not match
! the order asked for.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_taylor

  use iso_fortran_env      , only : dp => REAL64
  use util_derivative_terms, only : derivative_terms, mixed_partial
  use view_directed_stored , only : stored_directed_graph
  use field_stored         , only : stored_field
  use operation_family     , only : family
  use physics_integrand    , only : nodal_integrand
  use gti_block            , only : block_residual
  use gti_march            , only : block_of, solved, unknowns_graph
  use gti_sweeps           , only : jacobian_of, dense_solve

  implicit none

  private
  public :: nodal_coefficient, block_expansion

contains

  !===================================================================!
  ! The m-th coefficient of a nodal rule at every instant, from the
  ! state's coefficients up to that order and the design.
  !===================================================================!

  subroutine nodal_coefficient(integrand, degrees, at, series, design, order, values)

    class(nodal_integrand), intent(in) :: integrand
    integer               , intent(in) :: degrees, at(:), order
    real(dp)              , intent(in) :: series(0:, :)
    real(dp)              , intent(in) :: design
    real(dp), allocatable , intent(out) :: values(:)

    type(derivative_terms) :: q(0:degrees - 1), knob
    integer :: p, d, m

    if (order < 0) then
       error stop 'gti_taylor: an order is not negative'
    end if
    if (ubound(series, 1) < order) then
       error stop 'gti_taylor: the series carries every order asked for'
    end if

    knob = derivative_terms(design, order)
    if (order >= 1) call knob % set_symmetric(1, 1.0_dp)

    allocate(values(size(at)))

    do p = 1, size(at)
       do d = 0, degrees - 1
          q(d) = derivative_terms(0.0_dp, order)
          do m = 0, order
             call q(d) % set_symmetric(m, series(m, at(p) + d + 1))
          end do
       end do
       values(p) = mixed_partial(integrand % at_instant(q, knob))
    end do

  end subroutine nodal_coefficient

  !===================================================================!
  ! The governing coefficient placed on the rows it belongs to.
  !===================================================================!

  pure subroutine placed(governing, at, primary, carried, r)

    real(dp), intent(in)  :: governing(:)
    integer , intent(in)  :: at(:), primary, carried
    real(dp), intent(out) :: r(:)

    integer :: p

    r = 0.0_dp

    do p = 1, size(at)
       r(at(p) + primary + 1) = governing(p)
    end do

    r(1:carried) = 0.0_dp

  end subroutine placed

  !===================================================================!
  ! One block marched, then expanded: the functional and every
  ! derivative of it in the design, to the order asked for.
  !===================================================================!

  subroutine block_expansion(rows, physics, integrand, degrees, primary, &
       & instants_at, dt, design, max_order, q, f, achieved)

    type(block_residual)  , intent(in) :: rows
    class(nodal_integrand), intent(in) :: physics, integrand
    integer               , intent(in) :: degrees, primary, instants_at(:), max_order
    real(dp)              , intent(in) :: dt(:), design
    real(dp), allocatable , intent(out) :: q(:), f(:)
    real(dp)              , intent(out) :: achieved

    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, knobs
    real(dp), allocatable :: a(:,:), series(:,:)
    integer :: unknown_count

    unknown_count = rows % num_unknowns()
    call solved(rows, design, q, achieved)

    unknowns = stored_directed_graph(unknown_count, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), unknown_count)
    knobs    = stored_field('design', unknowns % vertex_set(), rows % num_points())
    call state % set_real_vector(q)
    call knobs % set_real_vector(spread(design, 1, rows % num_points()))

    call jacobian_of(rows, unknowns, [state, knobs], unknown_count, &
         & unknowns % vertex_set(), a)

    call state_series(physics, a, degrees, rows % points_at(), primary, &
         & rows % num_carried(), design, max_order, q, series)
    call functional_series(integrand, degrees, instants_at, design, max_order, &
         & series, dt, f)

  end subroutine block_expansion

  !===================================================================!
  ! The trajectory's coefficients, one order at a time. Each is one
  ! solve, against the jacobian the march already formed, for a right
  ! side the orders beneath it determine.
  !===================================================================!

  subroutine state_series(physics, a, degrees, at, primary, carried, design, &
       & max_order, q, series)

    class(nodal_integrand), intent(in) :: physics
    real(dp)              , intent(in) :: a(:,:), design, q(:)
    integer               , intent(in) :: degrees, at(:), primary, carried, max_order
    real(dp), allocatable , intent(out) :: series(:,:)

    real(dp), allocatable :: frozen(:,:), r(:), w(:), coefficient(:)
    integer :: m

    allocate(series(0:max_order, size(q)), source=0.0_dp)
    allocate(r(size(q)))
    series(0, :) = q

    do m = 1, max_order
       frozen = series
       frozen(m, :) = 0.0_dp
       call nodal_coefficient(physics, degrees, at, frozen, design, m, coefficient)
       call placed(coefficient, at, primary, carried, r)
       call dense_solve(a, -r, .false., w)
       series(m, :) = w
    end do

  end subroutine state_series

  !===================================================================!
  ! The functional and its derivatives, read off the trajectory's
  ! coefficients once they are all known.
  !===================================================================!

  subroutine functional_series(integrand, degrees, at, design, max_order, series, dt, f)

    class(nodal_integrand), intent(in) :: integrand
    integer               , intent(in) :: degrees, at(:), max_order
    real(dp)              , intent(in) :: design, series(0:, :), dt(:)
    real(dp), allocatable , intent(out) :: f(:)

    real(dp), allocatable :: coefficient(:)
    integer :: m

    allocate(f(0:max_order))

    do m = 0, max_order
       call nodal_coefficient(integrand, degrees, at, series, design, m, coefficient)
       f(m) = sum(dt * coefficient)
    end do

  end subroutine functional_series

end module gti_taylor
