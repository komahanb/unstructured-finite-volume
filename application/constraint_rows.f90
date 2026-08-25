! The constraints a coupling carries: the derived rows a scheme
! makes, and the governing row the physics makes.
!
! One block of five instants of a degree-two problem under a backward
! difference of order two. Each instant holds three components, so
! the block is square in fifteen unknowns and each instant owns three
! rows:
!
!    (k-1)*3 + 1   the governing constraint, van der Pol
!    (k-1)*3 + 2   the derived row that determines the velocity
!    (k-1)*3 + 3   the derived row that determines the acceleration
!
! The row the governing constraint takes is the family's primary
! degree, which for a backward difference is the value.
!
! The derived rows are checked on a state sampled from a polynomial
! the scheme reproduces: they must be zero wherever they exist. The
! governing row is not zero there, because a polynomial does not
! solve van der Pol, and printing the whole residual shows which row
! is which.
!
! The physics partials are then extracted as a row of numbers, one
! call per degree with a direction that is one at that degree, and
! compared both against the closed form and against a central
! difference. Degree three is included because its partial in the
! first derivative is exactly zero, which a partial indexed one place
! out would not be.
program constraint_rows

  use util_precision  , only : dp
  use view_directed_stored       , only : stored_directed_graph
  use field_calculus             , only : field
  use field_stored               , only : stored_field
  use operation_action           , only : variation
  use operation_stencil          , only : stencil
  use operation_scheme_stencil   , only : derived_constraints
  use operation_family_bdf       , only : bdf_family
  use operation_weight           , only : scheme_weight
  use physics_vanderpol          , only : van_der_pol

  implicit none

  integer , parameter :: order = 2
  integer , parameter :: num_instants = 5
  integer , parameter :: num_degrees = 3
  integer , parameter :: num_unknowns = num_instants * num_degrees
  real(dp), parameter :: nu = 1.0_dp

  call derived_rows()
  call physics_partials(2, [0.7_dp, -0.4_dp], [0.3_dp, 0.9_dp], [1.1_dp, -0.5_dp])
  call physics_partials(3, [0.7_dp, -0.4_dp], [0.3_dp, 0.9_dp], [1.1_dp, -0.5_dp])

contains

  pure integer function unknown(instant, degree)

    integer, intent(in) :: instant, degree

    unknown = (instant - 1) * num_degrees + degree + 1

  end function unknown

  !-------------------------------------------------------------------!
  ! The derived rows, on a state sampled from t squared: a backward
  ! difference of order two reproduces it in both rows.
  !-------------------------------------------------------------------!

  subroutine derived_rows()

    real(dp), parameter :: dt(num_instants) = [0.0_dp, 0.30_dp, 0.20_dp, 0.40_dp, 0.25_dp]

    type(stencil) :: rows
    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, direction
    class(field), allocatable :: out
    type(van_der_pol) :: physics
    real(dp), allocatable :: weight(:), residual(:), acted(:), governing(:)
    real(dp) :: q(num_unknowns), t(num_instants)
    integer , allocatable :: tails(:), heads(:), determines(:), source_degree(:)
    integer :: j, k

    call scheme_reach(tails, heads, source_degree, determines)
    call edge_weights(dt, tails, heads, source_degree, determines, weight)

    rows = derived_constraints( &
         & [(unknown(heads(j), determines(j)), j = 1, size(heads))], &
         & [(unknown(tails(j), source_degree(j)), j = 1, size(tails))], &
         & weight, num_unknowns, 'derived constraints')

    t(1) = 0.0_dp
    do k = 2, num_instants
       t(k) = t(k - 1) + dt(k)
    end do

    do k = 1, num_instants
       q(unknown(k, 0)) = t(k) * t(k)
       q(unknown(k, 1)) = 2.0_dp * t(k)
       q(unknown(k, 2)) = 2.0_dp
    end do

    unknowns = stored_directed_graph(num_unknowns, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), num_unknowns)
    call state % set_real_vector(q)

    call rows % apply(unknowns, [state], out)
    call out % real_vector(residual)

    direction = stored_field('v', unknowns % vertex_set(), num_unknowns)
    call direction % set_real_vector(q)
    call rows % partial_action(unknowns, [state], &
         & [variation(rows % argument(1), direction)], out)
    call out % real_vector(acted)

    physics = van_der_pol(2)
    call governing_rows(physics, dt, q, governing)
    call show_block(governing, residual, maxval(abs(residual - acted)))

  end subroutine derived_rows

  !-------------------------------------------------------------------!
  ! Every velocity row that fits, then every acceleration row that
  ! fits.
  !-------------------------------------------------------------------!

  subroutine scheme_reach(tails, heads, source_degree, determines)

    integer, allocatable, intent(out) :: tails(:), heads(:)
    integer, allocatable, intent(out) :: source_degree(:), determines(:)

    integer :: j, k

    tails = [((k - j, j = 0, order), k = order + 1, num_instants), &
         &   ((k - j, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    heads = [((k, j = 0, order), k = order + 1, num_instants), &
         &   ((k, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    determines = [((1, j = 0, order), k = order + 1, num_instants), &
         &        ((2, j = 0, 2 * order), k = 2 * order + 1, num_instants)]
    allocate(source_degree(size(tails)), source=0)

  end subroutine scheme_reach

  subroutine show_block(governing, residual, jacobian_gap)

    real(dp), intent(in) :: governing(:), residual(:), jacobian_gap

    integer :: k

    write(*,'(a)') ' the block on a state sampled from t squared'
    write(*,'(a)') '   instant   governing    velocity   acceleration'

    do k = 1, num_instants
       write(*,'(i10,3es13.2)') k, governing(k), &
            & residual(unknown(k, 1)), residual(unknown(k, 2))
    end do

    write(*,'(a)')        ' '
    write(*,'(a)')        ' the derived rows are linear, so the stencil is its own jacobian:'
    write(*,'(a,es11.2)') '   largest difference between apply and partial action ', &
         & jacobian_gap

  end subroutine show_block

  !-------------------------------------------------------------------!
  ! The governing row at every instant.
  !-------------------------------------------------------------------!

  subroutine governing_rows(physics, dt, q, r)

    type(van_der_pol), intent(in) :: physics
    real(dp)         , intent(in) :: dt(:), q(:)
    real(dp), allocatable, intent(out) :: r(:)

    type(stored_directed_graph) :: instants
    type(stored_field) :: state, design
    class(field), allocatable :: out

    associate (u1 => dt); end associate

    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', instants % vertex_set(), num_unknowns)
    design   = stored_field('nu', instants % vertex_set(), num_instants)
    call state  % set_real_vector(q)
    call design % set_real_vector(spread(nu, 1, num_instants))

    call physics % apply(instants, [state, design], out)
    call out % real_vector(r)

  end subroutine governing_rows

  !-------------------------------------------------------------------!
  ! The whole partial row of the governing constraint, one call per
  ! degree, beside the closed form and a central difference.
  !-------------------------------------------------------------------!

  subroutine physics_partials(degree, q0, q_top, design)

    integer , intent(in) :: degree
    real(dp), intent(in) :: q0(:), q_top(:), design(:)

    real(dp), parameter :: delta = 1.0e-6_dp
    integer , parameter :: instants = 2

    type(van_der_pol) :: physics
    type(stored_directed_graph) :: graph_of
    type(stored_field) :: state, nu_field, direction
    class(field), allocatable :: out
    real(dp), allocatable :: exact(:)
    real(dp) :: q(instants * (degree + 1)), v(instants * (degree + 1))
    real(dp) :: closed(0:degree), taken(0:degree), differenced(0:degree)
    integer :: d, k, nd

    nd      = degree + 1
    physics = van_der_pol(degree)

    call sample(degree, q0, q_top, design, instants, q)

    graph_of = stored_directed_graph(instants, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', graph_of % vertex_set(), size(q))
    nu_field = stored_field('nu', graph_of % vertex_set(), instants)
    direction = stored_field('v', graph_of % vertex_set(), size(q))
    call state    % set_real_vector(q)
    call nu_field % set_real_vector(spread(nu, 1, instants))

    do d = 0, degree
       v = 0.0_dp
       do k = 1, instants
          v((k - 1) * nd + d + 1) = 1.0_dp
       end do

       call direction % set_real_vector(v)
       call physics % partial_action(graph_of, [state, nu_field], &
            & [variation(physics % argument(1), direction)], out)
       call out % real_vector(exact)
       taken(d) = exact(1)

       differenced(d) = state_difference(physics, graph_of, state, nu_field, q, v)
       closed(d)      = closed_form(degree, d, q(1:nd))
    end do

    call show_row(degree, taken, closed, differenced)
    call design_partial(physics, graph_of, state, nu_field, q, instants, nd)

  end subroutine physics_partials

  subroutine show_row(degree, taken, closed, differenced)

    integer , intent(in) :: degree
    real(dp), intent(in) :: taken(0:), closed(0:), differenced(0:)

    integer :: d

    write(*,'(a)')         ' '
    write(*,'(a,i0)')      ' van der pol at degree ', degree
    write(*,'(a,9i12)')    '   partial in degree      ', [(d, d = 0, degree)]
    write(*,'(a,9f12.6)')  '   from partial_action    ', taken
    write(*,'(a,9f12.6)')  '   closed form            ', closed
    write(*,'(a,9f12.6)')  '   central difference     ', differenced

  end subroutine show_row

  !-------------------------------------------------------------------!
  ! One component of every instant set from the arguments given, the
  ! rest left at zero.
  !-------------------------------------------------------------------!

  subroutine sample(degree, q0, q_top, design, instants, q)

    integer , intent(in)  :: degree, instants
    real(dp), intent(in)  :: q0(:), q_top(:), design(:)
    real(dp), intent(out) :: q(:)

    integer :: k, nd

    nd = degree + 1
    q  = 0.0_dp

    do k = 1, instants
       q((k - 1) * nd + 1)  = q0(k)
       q((k - 1) * nd + nd) = q_top(k)
       if (degree >= 2) q((k - 1) * nd + nd - 1) = design(k)
    end do

  end subroutine sample

  !-------------------------------------------------------------------!
  ! A central difference of the residual along one state direction,
  ! the state left as it was found.
  !-------------------------------------------------------------------!

  function state_difference(physics, graph_of, state, nu_field, q, v) result(d)

    type(van_der_pol)          , intent(in)    :: physics
    type(stored_directed_graph), intent(in)    :: graph_of
    type(stored_field)         , intent(inout) :: state
    type(stored_field)         , intent(in)    :: nu_field
    real(dp)                   , intent(in)    :: q(:), v(:)
    real(dp) :: d

    real(dp), parameter :: delta = 1.0e-6_dp
    class(field), allocatable :: out
    real(dp), allocatable :: plus(:), minus(:)

    call state % set_real_vector(q + delta * v)
    call physics % apply(graph_of, [state, nu_field], out)
    call out % real_vector(plus)

    call state % set_real_vector(q - delta * v)
    call physics % apply(graph_of, [state, nu_field], out)
    call out % real_vector(minus)

    call state % set_real_vector(q)
    d = (plus(1) - minus(1)) / (2.0_dp * delta)

  end function state_difference

  !-------------------------------------------------------------------!
  ! The partial of the governing constraint in the design.
  !-------------------------------------------------------------------!

  subroutine design_partial(physics, graph_of, state, nu_field, q, instants, nd)

    type(van_der_pol)          , intent(in)    :: physics
    type(stored_directed_graph), intent(in)    :: graph_of
    type(stored_field)         , intent(inout) :: state, nu_field
    real(dp)                   , intent(in)    :: q(:)
    integer                    , intent(in)    :: instants, nd

    real(dp), parameter :: delta = 1.0e-6_dp
    type(stored_field) :: direction
    class(field), allocatable :: out
    real(dp), allocatable :: exact(:), plus(:), minus(:)
    real(dp) :: w(instants), q0, q_below

    associate (u1 => state); end associate

    w    = 1.0_dp
    q0      = q(1)
    q_below = q(nd - 1)

    direction = stored_field('w', graph_of % vertex_set(), instants)
    call direction % set_real_vector(w)

    call physics % partial_action(graph_of, [state, nu_field], &
         & [variation(physics % argument(2), direction)], out)
    call out % real_vector(exact)

    call nu_field % set_real_vector(spread(nu, 1, instants) + delta * w)
    call physics % apply(graph_of, [state, nu_field], out)
    call out % real_vector(plus)
    call nu_field % set_real_vector(spread(nu, 1, instants) - delta * w)
    call physics % apply(graph_of, [state, nu_field], out)
    call out % real_vector(minus)
    call nu_field % set_real_vector(spread(nu, 1, instants))

    write(*,'(a,f12.6)') '   partial in the design  ', exact(1)
    write(*,'(a,f12.6)') '   closed form            ', -(1.0_dp - q0 * q0) * q_below
    write(*,'(a,f12.6)') '   central difference     ', (plus(1) - minus(1)) / (2.0_dp * delta)

  end subroutine design_partial

  pure real(dp) function closed_form(degree, d, q) result(c)

    integer , intent(in) :: degree, d
    real(dp), intent(in) :: q(0:)

    if (d == degree) then
       c = 1.0_dp
    else if (d == degree - 1) then
       c = -nu * (1.0_dp - q(0) * q(0))
    else if (d == 0) then
       c = 2.0_dp * nu * q(0) * q(degree - 1) + 1.0_dp
    else
       c = 0.0_dp
    end if

  end function closed_form

  subroutine edge_weights(dt, tails, heads, source_degree, determines, w)

    real(dp), intent(in) :: dt(:)
    integer , intent(in) :: tails(:), heads(:), source_degree(:), determines(:)
    real(dp), allocatable, intent(out) :: w(:)

    type(stored_directed_graph) :: edges
    type(stored_field) :: steps, degrees, conditions
    type(scheme_weight) :: weights
    class(field), allocatable :: out

    edges  = stored_directed_graph(num_instants, tails=tails, heads=heads)
    steps  = stored_field('dt', edges % vertex_set(), num_instants)
    degrees    = stored_field('source degree', edges % edge_set(), size(tails))
    conditions = stored_field('determines', edges % edge_set(), size(tails))
    call steps      % set_real_vector(dt)
    call degrees    % set_integer_vector(source_degree)
    call conditions % set_integer_vector(determines)

    weights = scheme_weight(bdf_family(order))
    call weights % apply(edges, [steps, degrees, conditions], out)
    call out % real_vector(w)

  end subroutine edge_weights

end program constraint_rows
