! THE SPATIO-TEMPORAL EXAMPLE: the van der pol field
!
!      q^(N) - nu (1 - q^2) q^(N-1) + q - kappa laplacian(q) = 0
!
! over a two-dimensional domain and a duration, marched by the same
! block, the same newton and the same expansion as one node's
! equation, with the spatial level attached beneath the block.
!
! Three things are shown, each against something known:
!
!   the ode      at kappa = 0 with a constant field every node is one
!                node's equation, so the field's functional over the
!                area is the node's, and so is every derivative
!   the mode     at nu = 0 on a rectangle the field q = cos(pi x / a)
!                cos(pi y / b) cos(omega t) is exact, with omega^2 =
!                1 + kappa pi^2 (1/a^2 + 1/b^2), so the error at the
!                last instant is measured against it, and against the
!                semi-discrete solution that isolates the time error
!   the routes   the tangent and the adjoint agree on the field, both
!                against one factorisation
!
! and a run may write every instant as a vtu file for paraview.
!
!      ./spatio_temporal --config=field [--setting=value ...]
program spatio_temporal

  use util_precision        , only : dp
  use iso_fortran_env       , only : int64
  use gti_driver            , only : settings, steps_of, clock
  use gti_configuration     , only : configuration, read_configuration, override, show, &
       & worded, lists, refuse_unknown
  use gti_space             , only : room, spatial_mesh, spatial_operator, written_paraview, &
       & geometry_of, cartesian
  use operation_stencil     , only : stencil
  use field_calculus        , only : field
  use gti_field             , only : field_measure, field_startup, field_aggregates, spatial_rows
  use gti_march             , only : partitioned, set_stopping, block_of, unknowns_graph, &
       & unknown, consistent_states, frozen_inputs, &
       & set_sweep, swept, imbalance, by_tangent, by_adjoint, fresh_stamp
  use gti_taylor            , only : block_expansion
  use gti_sweeps            , only : design_partial, functional_gradient, &
       & set_linear_solver, set_aggregates, set_assembly, set_storage, set_multigrid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : uniform_grid, random_grid
  use operation_minimization, only : relative, absolute, by_count, by_rate
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use gti_block             , only : block_residual

  implicit none

  type(configuration) :: cfg
  type(room)          :: space
  real(dp), allocatable :: dt(:), t(:), q0(:), measure(:)
  real(dp) :: a, b, kappa, began
  integer  :: nd, n1, n2
  real(dp) :: x, y

  call settings('field', cfg)
  call show(cfg)

  call refuse_unknown(cfg % families, ['bdf  ', 'adams'], 'families')
  call refuse_unknown(cfg % spatial_grid, ['uniform', 'random '], 'spatial_grid')
  call refuse_unknown(cfg % initial_field, ['constant', 'mode    ', 'bump    '], 'initial_field')
  call refuse_unknown(cfg % export, ['none    ', 'paraview'], 'export')
  call refuse_unknown(cfg % check, ['none    ', 'ode     ', 'mode    ', 'operator'], 'check')

  call set_linear_solver(cfg % linear_solver)
  call set_assembly(cfg % assembly)
  call set_storage(cfg % storage)
  call set_multigrid(cfg % multigrid)
  call set_sweep(cfg % sweep)
  call set_stopping(cfg % tolerance, &
       & merge(relative, absolute, trim(cfg % tolerance_criterion) == 'relative'), &
       & merge(by_rate, by_count, trim(cfg % iteration_criterion) == 'by_rate'), &
       & cfg % max_iterations)

  nd    = cfg % state_degree + 1
  kappa = cfg % diffusion
  call pair_of(cfg % spatial_extent, a, b, 'extents')
  call pair_of(cfg % spatial_counts, x, y, 'counts')
  n1 = nint(x)
  n2 = nint(y)
  if (real(n1, dp) /= x .or. real(n2, dp) /= y) error stop 'spatio_temporal: a count is whole'

  began = clock()
  space = spatial_mesh(geometry_of(cfg % spatial_geometry), a, b, n1, n2, &
       & trim(cfg % spatial_grid) == 'random', cfg % seed)

  write(*,'(a,i0,a,i0,a,f12.6,a,i0,a,f9.3,a)') '   spatial mesh: cells ', space % num_cells, &
       & '   faces ', space % num_faces, '   area ', sum(space % volume), &
       & '   form degree ', cfg % spatial_order, '   built in ', clock() - began, ' s'

  call steps_of(cfg, dt, t)
  q0      = initial_field(cfg, space, nd, kappa)
  measure = field_measure(dt, space)

  if (trim(cfg % check) == 'operator') then
     call against_the_laplacian(cfg, space, kappa)
  else
     call table(cfg, space, dt, t, q0, measure, nd, kappa)
  end if

contains

  !-------------------------------------------------------------------!
  ! The operator alone against the laplacian of the mode, cell by
  ! cell: the balance over the area against kappa times minus pi^2
  ! (1/a^2 + 1/b^2) times the mode, which has no normal derivative at
  ! any wall. The error is reported over the cells that touch no
  ! wall, those that touch one, and those that touch two, so a wall's
  ! treatment is told apart from the interior's.
  !-------------------------------------------------------------------!

  subroutine against_the_laplacian(cfg, space, kappa)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    real(dp)           , intent(in) :: kappa

    real(dp), allocatable :: shape(:), balanced(:), exact(:)
    real(dp) :: pi, err(0:2), norm(0:2)
    integer  :: i, walls, count(0:2)

    if (space % geometry /= cartesian) then
       error stop 'spatio_temporal: the laplacian check is the rectangle''s'
    end if

    pi    = acos(-1.0_dp)
    shape = mode_shape(space, a, b)
    exact = -kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2) * shape

    began = clock()
    call balance_of(space, kappa, cfg % spatial_order, shape, balanced)
    write(*,'(a,f9.3,a)') '   the operator, built and applied in ', clock() - began, ' s'

    err   = 0.0_dp
    norm  = 0.0_dp
    count = 0
    do i = 1, space % num_cells
       walls = 0
       if (space % cell_ij(1, i) == 1 .or. space % cell_ij(1, i) == n2) walls = walls + 1
       if (space % cell_ij(2, i) == 1 .or. space % cell_ij(2, i) == n1) walls = walls + 1
       err(walls)   = err(walls)   + (balanced(i) / space % volume(i) - exact(i)) ** 2
       norm(walls)  = norm(walls)  + exact(i) ** 2
       count(walls) = count(walls) + 1
    end do

    write(*,'(a,i0,a,i0,a)') '   the operator against kappa laplacian of the mode, ', &
         & space % num_cells, ' cells, form degree ', cfg % spatial_order, ':'
    write(*,'(a,3(a,es10.3))') '   relative rms error', &
         & '   interior ', sqrt(err(0) / max(norm(0), tiny(1.0_dp))), &
         & '   one wall ', sqrt(err(1) / max(norm(1), tiny(1.0_dp))), &
         & '   corner ',   sqrt(err(2) / max(norm(2), tiny(1.0_dp)))

  end subroutine against_the_laplacian

  !-------------------------------------------------------------------!
  ! The rectangle's mode at every cell centre, cos(pi x / a) cos(pi y
  ! / b): the shape every check on the rectangle reads.
  !-------------------------------------------------------------------!

  pure function mode_shape(space, a, b) result(shape)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: a, b
    real(dp), allocatable :: shape(:)

    real(dp) :: pi
    integer  :: i

    pi = acos(-1.0_dp)
    shape = [(cos(pi * space % centre(1, i) / a) * cos(pi * space % centre(2, i) / b), &
         &    i = 1, space % num_cells)]

  end function mode_shape

  !-------------------------------------------------------------------!
  ! Two numbers from a setting, one per coordinate; a count is a
  ! whole one.
  !-------------------------------------------------------------------!

  subroutine pair_of(text, x, y, subject)

    character(len=*), intent(in)  :: text, subject
    real(dp)        , intent(out) :: x, y

    character(len=32), allocatable :: w(:)

    w = worded(text)
    if (size(w) /= 2) error stop 'spatio_temporal: two ' // subject // ', one per coordinate'
    read(w(1), *) x
    read(w(2), *) y

  end subroutine pair_of

  !-------------------------------------------------------------------!
  ! The field at the first instant: the components below the highest
  ! at every node, constant from the configuration or the rectangle's
  ! mode, and the highest solved from the physics with the level
  ! below attached.
  !-------------------------------------------------------------------!

  function initial_field(cfg, space, nd, kappa) result(q)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    integer            , intent(in) :: nd
    real(dp)           , intent(in) :: kappa
    real(dp), allocatable :: q(:)

    character(len=32), allocatable :: given(:)
    real(dp), allocatable :: lower(:,:)
    real(dp) :: pi
    integer  :: i, d

    allocate(lower(nd - 1, space % num_cells), source=0.0_dp)
    pi = acos(-1.0_dp)

    select case (trim(cfg % initial_field))
    case ('constant')
       given = worded(cfg % initial_state)
       if (size(given) > nd - 1) then
          error stop 'spatio_temporal: the initial state is given below the highest derivative'
       end if
       do d = 1, size(given)
          read(given(d), *) lower(d, 1)
       end do
       do i = 2, space % num_cells
          lower(:, i) = lower(:, 1)
       end do
    case ('mode')
       if (space % geometry /= cartesian) then
          error stop 'spatio_temporal: the mode is the rectangle''s'
       end if
       lower(1, :) = mode_shape(space, a, b)
    case ('bump')
       ! one plus half the rectangle's mode, on any geometry: a field
       ! that is not uniform, so the level below has something to do
       lower(1, :) = 1.0_dp + 0.5_dp * mode_shape(space, a, b)
    end select

    call set_aggregates(field_aggregates(space, 1, nd))
    q = consistent_states(van_der_pol(nd - 1), nd, lower, cfg % design, &
         & spatial_rows(space, kappa, cfg % spatial_order, nd, nd - 1, 1))

  end function initial_field

  !-------------------------------------------------------------------!
  ! One row per family and order. Stage families keep their instants
  ! between stages and are not yet laid out over a field.
  !-------------------------------------------------------------------!

  subroutine table(cfg, space, dt, t, q0, measure, nd, kappa)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    real(dp)           , intent(in) :: dt(:), t(:), q0(:), measure(:), kappa
    integer            , intent(in) :: nd

    class(family), allocatable :: scheme
    character(len=8), allocatable :: names(:)
    character(len=16) :: label
    integer :: i, order

    names = [character(len=8) :: 'bdf', 'adams']

    write(*,'(a)') ' '
    write(*,'(a)') '  scheme        f' // repeat(' ', 20) // 'derivatives in the design ...' // &
         & repeat(' ', 4) // 'tangent - adjoint     imbalance'

    do i = 1, size(names)
       if (.not. lists(cfg % families, trim(names(i)))) cycle
       do order = 1, cfg % max_discretization_order
          if (names(i) == 'bdf') then
             allocate(scheme, source=bdf_family(order))
          else
             allocate(scheme, source=adams_family(order))
          end if
          write(label,'(a,i0)') trim(names(i)), order
          call one_row(cfg, space, dt, t, q0, measure, nd, kappa, scheme, trim(label))
          deallocate(scheme)
       end do
    end do

  end subroutine table

  subroutine one_row(cfg, space, dt, t, q0, measure, nd, kappa, scheme, label)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    real(dp)           , intent(in) :: dt(:), t(:), q0(:), measure(:), kappa
    integer            , intent(in) :: nd
    class(family)      , intent(in) :: scheme
    character(len=*)   , intent(in) :: label

    type(block_residual) :: rows
    type(imbalance) :: left
    real(dp), allocatable :: held(:), q(:), f(:), marched(:)
    real(dp) :: achieved, tangent, adjoint, unused
    integer  :: h, n, m
    character(len=:), allocatable :: line
    character(len=20) :: cell

    n = cfg % instants
    h = scheme % history_depth(nd - 1)
    if (h >= n) return

    call set_aggregates(field_aggregates(space, 1 + (h - 1) * max(cfg % startup_refinement, 1), nd))
    held = field_startup(adams_family(2), van_der_pol(nd - 1), nd, h, &
         & max(cfg % startup_refinement, 1), dt, space, kappa, cfg % spatial_order, &
         & cfg % design, q0)

    rows = block_of(scheme, van_der_pol(nd - 1), nd, n, dt, held, nodes=space % num_cells, &
         & spatial=spatial_rows(space, kappa, cfg % spatial_order, nd, scheme % primary_degree(nd - 1), n))

    call set_aggregates(field_aggregates(space, n, nd))

    ! the state by the sweep chosen; the expansion from that state
    call swept(rows, cfg % design, space % num_cells, marched, achieved, left)
    call block_expansion(rows, van_der_pol(nd - 1), van_der_pol_energy(nd - 1), nd, &
         & scheme % primary_degree(nd - 1), rows % points_at(), measure, cfg % design, &
         & cfg % max_derivative_degree, q, f, unused, given=marched, nodes=space % num_cells)

    ! the two routes are compared where a derivative was asked for
    tangent = 0.0_dp
    adjoint = 0.0_dp
    if (cfg % max_derivative_degree >= 1) then
       call both_routes(rows, space, nd, n, measure, cfg % design, q, tangent, adjoint)
    end if

    line = '  ' // label // repeat(' ', max(1, 10 - len(label)))
    do m = lbound(f, 1), ubound(f, 1)
       write(cell,'(es20.11)') f(m)
       line = line // cell
    end do
    write(cell,'(es20.2)') tangent - adjoint
    line = line // cell
    write(cell,'(es14.2)') achieved
    line = line // cell
    if (.not. left % converged) line = line // '   unconverged'
    write(*,'(a)') line
    if (.not. left % converged) then
       write(*,'(a,es10.3,a,es10.3,a)') '      imbalance ', left % norm, ' against ', &
            & left % began, ' where the sweep began'
    end if

    if (trim(cfg % check) == 'ode')  call against_the_ode(cfg, space, dt, nd, scheme, held, f)
    if (trim(cfg % check) == 'mode') call against_the_mode(cfg, space, t, nd, kappa, q, n)
    if (trim(cfg % export) == 'paraview') call exported(cfg, space, nd, n, q, label)

  end subroutine one_row

  !-------------------------------------------------------------------!
  ! The gradient in the design by both routes, each one substitution
  ! against the same factorisation.
  !-------------------------------------------------------------------!

  subroutine both_routes(rows, space, nd, n, measure, design, q, tangent, adjoint)

    type(block_residual), intent(in)  :: rows
    type(room)          , intent(in)  :: space
    integer             , intent(in)  :: nd, n
    real(dp)            , intent(in)  :: measure(:), design, q(:)
    real(dp)            , intent(out) :: tangent, adjoint

    type(stored_directed_graph) :: unknowns, points
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: g(:), rate(:)
    integer :: num_points, mark

    num_points = n * space % num_cells
    points     = stored_directed_graph(num_points, tails=[integer ::], heads=[integer ::])
    call frozen_inputs(q, design, num_points, unknowns, inputs)

    call functional_gradient(van_der_pol_energy(nd - 1), points, inputs, &
         & measure, num_points, nd, unknowns % vertex_set(), g)
    call design_partial(rows, unknowns, inputs, num_points, unknowns % vertex_set(), rate)

    mark    = fresh_stamp()
    tangent = by_tangent(rows, unknowns, inputs, g, rate, 0.0_dp, space % num_cells, mark)
    adjoint = by_adjoint(rows, unknowns, inputs, g, rate, 0.0_dp, space % num_cells, mark)

  end subroutine both_routes

  !-------------------------------------------------------------------!
  ! At kappa = 0 with a constant field every node is one node's
  ! equation: the field's functional over the area is the node's,
  ! order by order. The node's march is the ordinary block from the
  ! same history, read at the first node.
  !-------------------------------------------------------------------!

  subroutine against_the_ode(cfg, space, dt, nd, scheme, held, f_field)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    real(dp)           , intent(in) :: dt(:), held(:), f_field(:)
    integer            , intent(in) :: nd
    class(family)      , intent(in) :: scheme

    type(block_residual) :: rows
    real(dp), allocatable :: held_node(:), q(:), f(:)
    real(dp) :: achieved, area
    integer  :: h, k, d, nodes
    character(len=:), allocatable :: line
    character(len=20) :: cell

    nodes = space % num_cells
    h     = scheme % history_depth(nd - 1)
    area  = sum(space % volume)

    held_node = [((held(unknown(k, d, nd, 1, nodes)), d = 0, nd - 1), k = 1, h)]

    rows = block_of(scheme, van_der_pol(nd - 1), nd, cfg % instants, dt, held_node)
    call block_expansion(rows, van_der_pol(nd - 1), van_der_pol_energy(nd - 1), nd, &
         & scheme % primary_degree(nd - 1), rows % points_at(), dt, cfg % design, &
         & cfg % max_derivative_degree, q, f, achieved)

    line = '      field / area over the node, less one:'
    do d = 0, ubound(f, 1) - lbound(f, 1)
       write(cell,'(es14.2)') f_field(lbound(f_field, 1) + d) / area / f(lbound(f, 1) + d) - 1.0_dp
       line = line // cell
    end do
    write(*,'(a)') line

  end subroutine against_the_ode

  !-------------------------------------------------------------------!
  ! At nu = 0 on a rectangle, the last instant against the exact mode
  ! and against the semi-discrete mode, whose frequency carries the
  ! discrete laplacian's eigenvalue on a uniform grid of spacings
  ! a / n1 and b / n2. The first error holds space and time, the
  ! second time alone.
  !-------------------------------------------------------------------!

  subroutine against_the_mode(cfg, space, t, nd, kappa, q, n)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    real(dp)           , intent(in) :: t(:), kappa, q(:)
    integer            , intent(in) :: nd, n

    real(dp) :: pi, omega, omega_h, exact, semi, e_exact, e_semi, area, mode
    real(dp), allocatable :: shape(:), balanced(:)
    integer  :: i, nodes

    if (space % geometry /= cartesian .or. cfg % design /= 0.0_dp) return

    nodes = space % num_cells
    pi    = acos(-1.0_dp)
    omega = sqrt(1.0_dp + kappa * pi ** 2 * (1.0_dp / a ** 2 + 1.0_dp / b ** 2))

    ! the semi-discrete frequency from the operator as built, by the
    ! rayleigh quotient of the mode: minus the mode against its own
    ! balance, over the mode against itself by area
    shape = mode_shape(space, a, b)
    call balance_of(space, kappa, cfg % spatial_order, shape, balanced)
    omega_h = sqrt(1.0_dp - dot_product(shape, balanced) / &
         & dot_product(shape, space % volume * shape))

    e_exact = 0.0_dp
    e_semi  = 0.0_dp
    area    = sum(space % volume)

    do i = 1, nodes
       mode  = shape(i)
       exact = mode * cos(omega   * t(n))
       semi  = mode * cos(omega_h * t(n))
       e_exact = e_exact + space % volume(i) * (q(unknown(n, 0, nd, i, nodes)) - exact) ** 2
       e_semi  = e_semi  + space % volume(i) * (q(unknown(n, 0, nd, i, nodes)) - semi) ** 2
    end do

    write(*,'(a,es12.3,a,es12.3,a,f10.6,a,f10.6)') &
         & '      error at the last instant, against the mode ', sqrt(e_exact / area), &
         & '   semi-discrete ', sqrt(e_semi / area), '   omega ', omega, '   omega_h ', omega_h

  end subroutine against_the_mode

  !-------------------------------------------------------------------!
  ! The operator applied to a field over the cells: the integrated
  ! flux balance of that field.
  !-------------------------------------------------------------------!

  subroutine balance_of(space, kappa, degree, values, balanced)

    type(room), intent(in) :: space
    real(dp)  , intent(in) :: kappa, values(:)
    integer   , intent(in) :: degree
    real(dp), allocatable, intent(out) :: balanced(:)

    type(stencil) :: op
    type(stored_field) :: given
    class(field), allocatable :: out

    op    = spatial_operator(space, kappa, degree)
    given = stored_field('values', op % pattern % vertex_set(), size(values))
    call given % set_real_vector(values)
    call op % apply(op % pattern, [given], out)
    call out % real_vector(balanced)

  end subroutine balance_of

  !-------------------------------------------------------------------!
  ! Every instant as one vtk file, numbered, so paraview reads the
  ! series as time.
  !-------------------------------------------------------------------!

  subroutine exported(cfg, space, nd, n, q, label)

    type(configuration), intent(in) :: cfg
    type(room)         , intent(in) :: space
    integer            , intent(in) :: nd, n
    real(dp)           , intent(in) :: q(:)
    character(len=*)   , intent(in) :: label

    character(len=8), allocatable :: names(:)
    real(dp), allocatable :: values(:,:)
    character(len=256) :: path
    integer :: k, i, d, nodes

    nodes = space % num_cells
    allocate(names(nd), values(nodes, nd))
    do d = 0, nd - 1
       write(names(d + 1),'(a,i0)') 'q', d
    end do

    do k = 1, n
       do i = 1, nodes
          do d = 0, nd - 1
             values(i, d + 1) = q(unknown(k, d, nd, i, nodes))
          end do
       end do
       write(path,'(a,a,a,a,i4.4,a)') trim(cfg % export_path), '_', label, '_', k, '.vtu'
       call written_paraview(space, trim(path), names, values)
    end do

    write(*,'(a,i0,a,a,a)') '      written ', n, ' files ', trim(cfg % export_path) // '_' // label, '_*.vtu'

  end subroutine exported

end program spatio_temporal
