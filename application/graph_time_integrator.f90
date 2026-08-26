! The graph time integrator.
!
! It assembles what src holds and prints one table: a row for every
! scheme a configuration asks for, and a column for the functional
! and each of its derivatives in the design.
!
! Every row is one march and one expansion. The march solves the
! block whole, over a partition of the duration that is not uniform -
! the steps are drawn from a seed and scaled so that they sum to the
! duration exactly - and the expansion then reads the functional and
! its derivatives off the same jacobian, one solve per order.
!
!      ./graph_time_integrator --config=homogeneous
!      ./graph_time_integrator --config=homogeneous --max-derivative-degree=5
!
! A setting on the command line overrides the one in the file, and a
! setting that is not a setting stops the run rather than being
! passed over.
!
!             WHERE EVERY ROW STARTS
!
! A scheme cannot take its first step until it has instants behind it
! to look back at, and how many differs by family: a backward
! difference of order four wants eight, an Adams quadrature of order
! one wants a single one, and a stage family wants a single one
! whatever its order. Whatever fills those instants is not solved by
! the scheme; it is handed to it.
!
! If it were filled from a formula the rows would not be comparable.
! The widest scheme would hold a third of the horizon at numbers that
! are not a trajectory, would only integrate what remained, and would
! begin that from a state the equation would never have produced. Its
! functional would be mostly the formula and the narrowest scheme's
! mostly a solution, and the two would have no reason to agree.
!
! So the instants are integrated rather than invented. A stage family
! needs one instant and therefore no filler at all, so one is marched
! first over a refined grid across the startup, and every row takes
! its own reach from what that produced. That march is itself a chain
! of short blocks rather than one long one, because a block is solved
! whole and a long one costs far more than the several it could have
! been - which is the same junction the table's own rows use. Every row then begins from
! the same trajectory, holding one instant of it or eight is equally
! sound, and what separates the rows is how well each integrates,
! which is what the table is for.
!
! How many instants a row is handed and how many it works out is
! printed beside it, because they differ sharply: a backward
! difference of order four on an equation of degree four looks back
! over sixteen, so on a horizon of twenty-one it integrates five. Its
! functional is then mostly what it was given, and a reader who did
! not know that would take it for a peer of a row that integrated
! twenty.
!
! That is what automatic_order_conservation asks for. Turned off, the
! startup would have to be filled some other way, and there is no
! other way here that keeps the rows comparable, so the run says so
! and stops rather than printing a table that cannot be read across.
program graph_time_integrator

  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_grid        , only : uniform_grid, random_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use operation_grid        , only : grid
  use gti_march             , only : partitioned, set_stopping, consistent_state, imbalance, set_sweep, &
       & weight_of, precision_needed
  use util_precision        , only : precision_named
  use iso_fortran_env       , only : real128
  use gti_expansion         , only : family_holder
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & expansion_substitutions, &
       & started => startup_trajectory
  use gti_sweeps            , only : set_linear_solver, set_assembly, set_storage, set_multigrid
  use operation_minimization, only : relative, absolute, by_count, by_rate
  use gti_driver            , only : settings, chosen_grid, steps_of, family_named
  use gti_configuration     , only : configuration, read_configuration, override, show, &
       & lists, refuse_unknown, worded
  use util_tally            , only : tally_open, tally_close, tally_order, &
       & tally_enter, tally_leave, tally_amount, tally_event_of, &
       & tally_num_levels, tally_level_name, tally_event_name, &
       & at_expansion, at_horizon, wall_time

  implicit none


  type(configuration) :: cfg

  call settings('homogeneous', cfg)
  call show(cfg)
  call set_linear_solver(cfg % linear_solver)
  call set_assembly(cfg % assembly)
  call set_storage(cfg % storage)
  call set_multigrid(cfg % multigrid)
  call set_sweep(cfg % sweep)
  call table(cfg)

contains

  !-------------------------------------------------------------------!
  ! The one instant a stage family needs, and it is consistent with
  ! the equation rather than merely plausible: the value and every
  ! derivative below the highest are chosen, and the highest is what
  ! the governing constraint then requires.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The state at the first instant: the components given, and the
  ! highest solved from the physics so that the state is consistent
  ! with it. What a hand would have set here is exactly what comes
  ! back, for the one physics and one value a hand had in mind, and
  ! for every other the hand would have been wrong.
  !-------------------------------------------------------------------!

  function at_rest(cfg) result(q)

    type(configuration), intent(in) :: cfg
    real(dp), allocatable :: q(:)

    character(len=32), allocatable :: given(:)
    real(dp), allocatable :: lower(:)
    integer :: nd, i

    nd    = cfg % state_degree + 1
    given = worded(cfg % initial_state)

    if (size(given) > nd - 1) then
       write(*,'(a,i0,a)') ' the initial state holds the ', nd - 1, &
            & ' components below the highest, which the physics gives.'
       error stop 'graph_time_integrator: the initial state is given below the highest derivative'
    end if

    allocate(lower(nd - 1), source=0.0_dp)
    do i = 1, size(given)
       read(given(i), *) lower(i)
    end do

    q = consistent_state(van_der_pol(cfg % state_degree), nd, lower, cfg % design)

  end function at_rest

  !-------------------------------------------------------------------!
  ! How many instants the widest row that fits looks back over. A row
  ! that reaches past the horizon is not built, so it does not decide
  ! how long a startup the others need; zero means none of them fit.
  !-------------------------------------------------------------------!

  integer function widest_reach(cfg) result(widest)

    type(configuration), intent(in) :: cfg

    character(len=8) :: every(3)
    class(family), allocatable :: scheme
    logical :: staged, ok
    integer :: i, order, reach

    every  = ['bdf     ', 'adams   ', 'dirk    ']
    widest = 0

    ! how far back a family looks is the family's own answer
    do i = 1, 3
       if (index(cfg % families, trim(every(i))) == 0) cycle
       do order = 1, cfg % max_discretization_order
          call chosen(trim(every(i)), order, scheme, staged, ok)
          if (.not. ok) cycle
          reach = scheme % history_depth(cfg % state_degree)
          if (reach < cfg % instants) widest = max(widest, reach)
       end do
    end do

  end function widest_reach

  !-------------------------------------------------------------------!
  ! The instants every row starts from, integrated rather than
  ! invented: a stage family over the startup, on a grid refined
  ! within each of its steps, sampled back at the coarse instants.
  !-------------------------------------------------------------------!
  ! The first widest instants, marched by the chain's startup from the
  ! state at rest.
  !-------------------------------------------------------------------!

  subroutine startup_trajectory(cfg, dt, widest, held)

    type(configuration), intent(in)  :: cfg
    real(dp)           , intent(in)  :: dt(:)
    integer            , intent(in)  :: widest
    real(dp), allocatable, intent(out) :: held(:)

    call started(van_der_pol(cfg % state_degree), cfg % state_degree + 1, widest, &
         & cfg % startup_refinement, dt, cfg % design, at_rest(cfg), held)

  end subroutine startup_trajectory

  !-------------------------------------------------------------------!
  ! One family, by name and order. A stage family says that it is
  ! one, since its block is laid out differently.
  !-------------------------------------------------------------------!

  subroutine chosen(name, order, scheme, staged, ok)

    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: staged, ok

    call family_named(name, order, scheme, ok)
    staged = .false.
    if (ok) staged = scheme % num_stages() > 1

  end subroutine chosen


  function labelled(names, orders) result(text)

    character(len=*), intent(in) :: names(:)
    integer         , intent(in) :: orders(:)
    character(len=:), allocatable :: text

    character(len=2) :: digit
    integer :: b

    text = ''
    do b = 1, size(names)
       write(digit,'(i0)') orders(b)
       if (b > 1) text = text // '-'
       text = text // trim(names(b)) // trim(digit)
    end do

  end function labelled

  !-------------------------------------------------------------------!
  ! The heading, then one line per row.
  !-------------------------------------------------------------------!

  subroutine shown_initial(cfg)

    type(configuration), intent(in) :: cfg

    real(dp), allocatable :: q(:)
    character(len=:), allocatable :: line
    character(len=18) :: cell
    integer :: i

    q = at_rest(cfg)
    line = '   initial state, consistent'
    do i = 1, size(q)
       write(cell,'(es18.10)') q(i)
       line = line // cell
    end do
    write(*,'(a)') line

  end subroutine shown_initial

  subroutine heading(cfg)

    type(configuration), intent(in) :: cfg

    character(len=:), allocatable :: line, name
    integer :: m

    line = '  scheme' // repeat(' ', 14) // 'solved' // repeat(' ', 10)

    ! Each name sits over its own column, right against the digits.
    do m = 0, cfg % max_derivative_degree
       name = order_named(m)
       line = line // repeat(' ', 20 - len(name)) // name // ' '
    end do

    write(*,'(a)') ' '
    write(*,'(a)') line

  end subroutine heading

  subroutine show_row(label, solved, f, left, columns)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: solved
    real(dp)        , intent(in) :: f(0:)
    type(imbalance) , intent(in) :: left
    integer         , intent(in) :: columns

    character(len=21) :: cell
    character(len=6)  :: counted
    character(len=:), allocatable :: line
    integer :: m

    line = '  ' // label // repeat(' ', max(2, 20 - len(label)))
    write(counted,'(i6)') solved
    line = line // counted // repeat(' ', 10)

    ! A column the row holds no expansion for is left empty rather
    ! than filled, there being no number to state under it.
    do m = 0, columns
       if (m <= ubound(f, 1)) then
          write(cell,'(es20.11)') f(m)
       else
          write(cell,'(a20)') '-'
       end if
       line = line // cell
    end do

    if (.not. left % converged) then
       if (left % diverging) then
          line = line // '   diverging'
       else
          line = line // '   unconverged'
       end if
    end if

    write(*,'(a)') line

    if (.not. left % converged) call shown_aspect(left)

  end subroutine show_row

  !-------------------------------------------------------------------!
  ! What the march left, by aspect, beneath the row that did not
  ! converge: how the norm splits by degree, where the largest entry
  ! sits, and which state the norm is steepest in.
  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!
  ! The precision a row's target needs, from the weight its family
  ! carries at its smallest step and the state it reached: shown
  ! beneath every row under accounting, and beneath any row whose
  ! need exceeds this build.
  !-------------------------------------------------------------------!

  subroutine shown_precision(scheme, nd, dt, chain, left, cfg)

    class(family)      , intent(in) :: scheme
    integer            , intent(in) :: nd
    real(dp)           , intent(in) :: dt(:)
    type(chain_block)  , intent(in) :: chain(:)
    type(imbalance)    , intent(in) :: left
    type(configuration), intent(in) :: cfg

    real(dp) :: weight, state_size
    real(real128) :: needed
    character(len=:), allocatable :: least
    integer :: b

    weight     = weight_of(scheme, nd, minval(dt(2:)))
    state_size = 0.0_dp
    do b = 1, size(chain)
       state_size = max(state_size, maxval(abs(chain(b) % state)))
    end do

    call precision_needed(weight, state_size, left % began, needed, least)

    if (.not. cfg % accounting .and. least == precision_named()) return
    if (.not. cfg % accounting .and. least == 'single') return

    write(*,'(a,es9.2,a,es9.2,a,es9.2,a,a,a,a)') '      precision  ||A|| ', weight, &
         & '  ||q|| ', state_size, '  spacing needed ', real(needed, dp), &
         & '  least kind ', least, '  this build ', precision_named()

  end subroutine shown_precision

  subroutine shown_aspect(left)

    type(imbalance), intent(in) :: left

    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: d

    write(*,'(a,es10.3,a,es10.3,a)') '      imbalance ', left % norm, &
         & ' against ', left % began, ' where the march began'

    line = '      by degree '
    do d = 0, ubound(left % by_degree, 1)
       write(cell,'(es14.3)') left % by_degree(d)
       line = line // cell
    end do
    write(*,'(a)') line

    write(*,'(a,i0,a,i0)') '      largest entry at slot ', left % worst_slot, &
         & ' degree ', left % worst_degree
    write(*,'(a,i0,a,i0,a,es10.3)') '      steepest in the state at slot ', &
         & left % steepest_slot, ' degree ', left % steepest_degree, &
         & ', d||r||/dq = ', left % steepest

  end subroutine shown_aspect

  !-------------------------------------------------------------------!
  ! One row: a chain of blocks, marched and then expanded. A chain of
  ! one is a homogeneous row and takes the same path.
  !-------------------------------------------------------------------!

  subroutine one_row(cfg, startup, names, orders, printed)

    type(configuration), intent(in)    :: cfg
    real(dp)           , intent(in)    :: startup(:)
    character(len=*)   , intent(in)    :: names(:)
    integer            , intent(in)    :: orders(:)
    integer            , intent(inout) :: printed

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    integer , allocatable :: added(:)
    real(dp), allocatable :: dt(:), t(:), held(:), f(:)
    type(imbalance) :: left
    real(dp) :: achieved
    integer :: b, nd, given, reported
    logical :: ok

    nd = cfg % state_degree + 1
    allocate(schemes(size(names)), added(size(names)))
    call assembled(cfg, names, orders, schemes, added, ok)
    if (.not. ok) return

    call steps_of(cfg, dt, t)
    given = schemes(1) % scheme % history_depth(nd - 1)
    if (given * nd > size(startup)) return
    held  = startup(1:given * nd)

    call tally_enter(at_expansion)
    call tally_order(0)

    call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
         & chosen_grid(cfg), cfg % design, held, chain, dt, t, achieved, left=left)

    ! Every derivative is taken at the state the march reached, so a
    ! row that did not converge has none to take and only its value is
    ! expanded.
    if (.not. left % converged) then
       reported = 0
    else
       reported = cfg % max_derivative_degree
    end if

    call chain_expansion(chain, van_der_pol(cfg % state_degree), &
         & van_der_pol_energy(cfg % state_degree), nd, dt, cfg % design, &
         & reported, f)

    call tally_leave()

    call show_row(labelled(names, orders), cfg % instants - given, f, left, &
         & cfg % max_derivative_degree)
    call shown_precision(schemes(1) % scheme, nd, dt, chain, left, cfg)
    printed = printed + 1

    associate (u1 => b); end associate

  end subroutine one_row

  !-------------------------------------------------------------------!
  ! The families a row names, and the instants split among them. A
  ! row whose blocks would add no more instants than their families
  ! reach back over is not built.
  !-------------------------------------------------------------------!

  subroutine assembled(cfg, names, orders, schemes, added, ok)

    type(configuration), intent(in)  :: cfg
    character(len=*)   , intent(in)  :: names(:)
    integer            , intent(in)  :: orders(:)
    type(family_holder), intent(inout) :: schemes(:)
    integer            , intent(inout) :: added(:)
    logical            , intent(out)   :: ok

    class(family), allocatable :: scheme
    logical :: staged, exists
    integer :: b, blocks, share

    blocks = size(names)
    ok = .true.

    share = cfg % instants / blocks
    added = share
    added(1) = cfg % instants - share * (blocks - 1)

    do b = 1, blocks
       call chosen(names(b), orders(b), scheme, staged, exists)
       if (.not. exists) then
          ok = .false.
          cycle
       end if
       allocate(schemes(b) % scheme, source=scheme)
       deallocate(scheme)
       if (added(b) <= schemes(b) % scheme % history_depth(cfg % state_degree)) ok = .false.
    end do

  end subroutine assembled


  !-------------------------------------------------------------------!
  ! Every row the configuration asks for.
  !-------------------------------------------------------------------!

  subroutine table(cfg)

    type(configuration), intent(in) :: cfg

    real(dp), allocatable :: startup(:), dt(:), t(:)
    integer :: widest, printed

    ! Before anything is measured against the horizon, since a word
    ! this program has nothing for reaches back over nothing and
    ! would be reported as a horizon too narrow to hold it.
    !
    ! physics is refused rather than dispatched on because one
    ! integrand is built. Were it neither, the run would state a
    ! physics in its heading and integrate a different one.
    call refuse_unknown(cfg % physics, ['vanderpol'], 'physics')
    call refuse_unknown(cfg % tolerance_criterion, ['relative', 'absolute'], &
         & 'tolerance_criterion')
    call refuse_unknown(cfg % iteration_criterion, ['by_rate ', 'by_count'], &
         & 'iteration_criterion')

    call set_stopping(cfg % tolerance, &
         & merge(relative, absolute, trim(cfg % tolerance_criterion) == 'relative'), &
         & merge(by_rate, by_count, trim(cfg % iteration_criterion) == 'by_rate'), &
         & cfg % max_iterations)
    if (cfg % accounting) then
       call refuse_unknown(cfg % measurements, &
            & ['wall_time     ', 'primal_loops  ', 'tangent_loops ', &
            &  'adjoint_loops ', 'newton_solves ', 'linear_solves ', &
            &  'factorisations'], 'measurements')
    end if
    call refuse_unknown(cfg % families, ['bdf     ', 'adams   ', 'dirk    '], 'families')
    call refuse_unknown(cfg % combinations, &
         & ['homogeneous', 'pairs      ', 'triples    '], 'combinations')

    widest = widest_reach(cfg)

    if (.not. cfg % automatic_order_conservation) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' the widest row here looks back over instants: ', widest
       write(*,'(a)')    ' filling them by any other means leaves the rows solving different'
       write(*,'(a)')    ' problems from different starting states, and no table read across'
       write(*,'(a)')    ' such rows means anything.'
       error stop 'graph_time_integrator: order conservation is the only startup built'
    end if

    if (widest == 0) then
       write(*,'(a)')    ' '
       write(*,'(a,i0)') ' every family and order asked for looks further back than the'
       write(*,'(a,i0)') ' horizon holds, which is instants: ', cfg % instants
       error stop 'graph_time_integrator: no row fits in this horizon'
    end if

    call steps_of(cfg, dt, t)
    call startup_trajectory(cfg, dt, widest, startup)

    call shown_initial(cfg)
    write(*,'(a,a)') '   precision of this build  ', precision_named()
    call heading(cfg)

    if (cfg % accounting) call tally_open(cfg % max_derivative_degree)

    printed = 0
    if (asked(cfg, 'homogeneous')) call tuple_rows(cfg, startup, 1, printed)
    if (asked(cfg, 'pairs'))       call tuple_rows(cfg, startup, 2, printed)
    if (asked(cfg, 'triples'))     call tuple_rows(cfg, startup, 3, printed)

    if (cfg % accounting) then
       call tally_close()
       call accounted(cfg)
    end if

    if (printed == 0) then
       write(*,'(a)') ' '
       write(*,'(a)') ' no row was built. A family has no scheme at every order - a stage'
       write(*,'(a)') ' family has none below order two - and a row whose blocks would add'
       write(*,'(a)') ' no more instants than they look back over is not built either.'
    end if

  end subroutine table

  pure logical function asked(cfg, what) result(yes)

    type(configuration), intent(in) :: cfg
    character(len=*)   , intent(in) :: what

    yes = lists(cfg % combinations, what)

  end function asked


  !-------------------------------------------------------------------!
  ! The names a configuration lists, in the order it lists them.
  !-------------------------------------------------------------------!

  function listed(cfg) result(list)

    type(configuration), intent(in) :: cfg
    character(len=8), allocatable :: list(:)

    character(len=8) :: every(3)
    integer :: i, n

    every = ['bdf     ', 'adams   ', 'dirk    ']
    n = 0
    do i = 1, 3
       if (lists(cfg % families, trim(every(i)))) n = n + 1
    end do

    allocate(list(n))
    n = 0
    do i = 1, 3
       if (lists(cfg % families, trim(every(i)))) then
          n = n + 1
          list(n) = every(i)
       end if
    end do

  end function listed

  !-------------------------------------------------------------------!
  ! Every ordered tuple of this many distinct families, at every
  ! order - one order for the whole tuple, or, when mixed orders are
  ! asked for, every tuple of orders. One arity serves the
  ! homogeneous rows, the pairs and the triples alike.
  !-------------------------------------------------------------------!

  subroutine tuple_rows(cfg, startup, arity, printed)

    type(configuration), intent(in)    :: cfg
    real(dp)           , intent(in)    :: startup(:)
    integer            , intent(in)    :: arity
    integer            , intent(inout) :: printed

    character(len=8), allocatable :: names(:)
    integer :: which(arity), orders(arity)
    integer :: m, code, k, r, order

    names = listed(cfg)
    m     = size(names)

    do code = 0, m ** arity - 1
       ! the tuple of families, the last position varying fastest
       r = code
       do k = arity, 1, -1
          which(k) = mod(r, m) + 1
          r        = r / m
       end do
       if (any([(any(which(1:k-1) == which(k)), k = 2, arity)])) cycle

       if (cfg % mixed_orders .and. arity >= 2) then
          code_of_orders: block
            integer :: oc
            do oc = 0, cfg % max_discretization_order ** arity - 1
               r = oc
               do k = arity, 1, -1
                  orders(k) = mod(r, cfg % max_discretization_order) + 1
                  r         = r / cfg % max_discretization_order
               end do
               call one_row(cfg, startup, names(which), orders, printed)
            end do
          end block code_of_orders
       else
          do order = 1, cfg % max_discretization_order
             orders = order
             call one_row(cfg, startup, names(which), orders, printed)
          end do
       end if
    end do

  end subroutine tuple_rows

  !-------------------------------------------------------------------!
  ! What the run spent, one table per measurement asked for: the
  ! amount at each level of the hierarchy against the derivative order
  ! it was spent on, and then the same amounts as ratios of one order
  ! to another.
  !
  ! The ratio is what a higher order costs against a lower one, so the
  ! entry at row i and column j is the amount at order i over the
  ! amount at order j. A column whose order spent nothing leaves its
  ! ratio empty rather than dividing by it.
  !-------------------------------------------------------------------!

  subroutine accounted(cfg)

    type(configuration), intent(in) :: cfg

    character(len=32), allocatable :: wanted(:)
    integer :: i, event

    wanted = worded(cfg % measurements)

    do i = 1, size(wanted)
       event = tally_event_of(trim(wanted(i)))
       call one_measurement(cfg, event)
    end do

    call route_note(cfg)
    call cliff_note(cfg)

  end subroutine accounted

  !-------------------------------------------------------------------!
  ! What the route's cost model said each order would cost in
  ! substitutions, beside what was counted. One block per row is what
  ! the homogeneous table builds, so the model is read at one block.
  !-------------------------------------------------------------------!

  subroutine route_note(cfg)

    type(configuration), intent(in) :: cfg

    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: m
    real(dp) :: counted, rows

    rows = over_levels(0, 5)
    if (rows <= 0.0_dp) return

    write(*,'(a)') ' '
    write(*,'(a)') '   tangent substitutions per row, the model against the count'
    line = '   model            '
    do m = 1, cfg % max_derivative_degree
       write(cell,'(i14)') expansion_substitutions(1, m)
       line = line // cell
    end do
    write(*,'(a)') line

    line = '   counted          '
    do m = 1, cfg % max_derivative_degree
       counted = over_levels(m, 3)
       write(cell,'(f14.2)') counted / rows
       line = line // cell
    end do
    write(*,'(a)') line

  end subroutine route_note

  subroutine one_measurement(cfg, event)

    type(configuration), intent(in) :: cfg
    integer            , intent(in) :: event

    real(dp), allocatable :: whole(:)
    character(len=:), allocatable :: line
    integer :: level, m, top

    top = cfg % max_derivative_degree
    allocate(whole(0:top), source=0.0_dp)

    ! Levels nest, so a level's time already holds the time of the
    ! levels opened inside it and a sum over levels would count the
    ! same seconds again. The expansion is opened once for a whole
    ! row and closed after every order has been taken, so its time
    ! belongs to no single order and is filed where the row began.
    ! The horizon is opened once per order, which is what a time
    ! against an order means, so it is the one the ratios are taken
    ! from. A count is filed at one level only and does sum.
    do m = 0, top
       if (event == wall_time) then
          whole(m) = tally_amount(at_horizon, m, event)
       else
          whole(m) = over_levels(m, event)
       end if
    end do

    write(*,'(a)') ' '
    write(*,'(a)') ' accounting: ' // tally_event_name(event)
    if (event == wall_time) then
       write(*,'(a)') '   seconds. A level holds the levels opened inside it, and the'
       write(*,'(a)') '   expansion spans every order, so the ratios are the horizon.'
    end if
    write(*,'(a)') ' '

    line = '   at each level    '
    do m = 0, top
       line = line // right(order_named(m))
    end do
    write(*,'(a)') line

    do level = 1, tally_num_levels()
       line = '   ' // tally_level_name(level) // &
            & repeat(' ', max(1, 18 - len(tally_level_name(level))))
       do m = 0, top
          line = line // right(amount_text(tally_amount(level, m, event), event))
       end do
       write(*,'(a)') line
    end do

    if (event == wall_time) then
       line = '   per order        '
    else
       line = '   whole run        '
    end if
    do m = 0, top
       line = line // right(amount_text(whole(m), event))
    end do
    write(*,'(a)') line

    call ratio_matrix(whole, top)

  end subroutine one_measurement

  !-------------------------------------------------------------------!
  ! Row over column, the whole run. An order that spent nothing is no
  ! denominator, and its column is left empty.
  !-------------------------------------------------------------------!

  subroutine ratio_matrix(whole, top)

    real(dp), intent(in) :: whole(0:)
    integer , intent(in) :: top

    character(len=:), allocatable :: line
    character(len=14) :: cell
    integer :: i, j

    write(*,'(a)') ' '
    line = '   row over column  '
    do j = 0, top
       line = line // right(order_named(j))
    end do
    write(*,'(a)') line

    do i = 0, top
       line = '   ' // order_named(i) // repeat(' ', max(1, 18 - len(order_named(i))))
       do j = 0, top
          if (whole(j) > 0.0_dp) then
             write(cell,'(f14.2)') whole(i) / whole(j)
          else
             write(cell,'(a14)') '-'
          end if
          line = line // cell
       end do
       write(*,'(a)') line
    end do

  end subroutine ratio_matrix

  !-------------------------------------------------------------------!
  ! Which side of the iteration cap the run sits on. Past it every
  ! order spends the whole budget instead of converging, and a ratio
  ! measured there reports the cap and not the order.
  !-------------------------------------------------------------------!

  subroutine cliff_note(cfg)

    type(configuration), intent(in) :: cfg

    real(dp) :: loops, solves
    integer  :: m

    loops  = 0.0_dp
    solves = 0.0_dp

    do m = 0, cfg % max_derivative_degree
       loops  = loops  + over_levels(m, 2)
       solves = solves + over_levels(m, 5)
    end do

    if (solves <= 0.0_dp) return

    write(*,'(a)') ' '
    write(*,'(a,f8.1)') '   primal loops per newton solve      ', loops / solves
    if (loops / solves >= 39.0_dp) then
       write(*,'(a)') '   at the iteration budget: the march is not converging, so'
       write(*,'(a)') '   these ratios report the budget and not the derivative order.'
    end if

  end subroutine cliff_note

  function amount_text(spent, event) result(text)

    real(dp), intent(in) :: spent
    integer , intent(in) :: event
    character(len=:), allocatable :: text

    character(len=14) :: cell

    if (event == wall_time) then
       write(cell,'(f14.4)') spent
    else
       write(cell,'(i14)') nint(spent)
    end if
    text = trim(adjustl(cell))

  end function amount_text

  !-------------------------------------------------------------------!
  ! An amount summed over every level of the hierarchy, for one order
  ! and one event.
  !-------------------------------------------------------------------!

  real(dp) function over_levels(m, event) result(total)

    integer, intent(in) :: m, event

    integer :: level

    total = 0.0_dp
    do level = 1, tally_num_levels()
       total = total + tally_amount(level, m, event)
    end do

  end function over_levels

  function order_named(m) result(named)

    integer, intent(in) :: m
    character(len=:), allocatable :: named

    character(len=2) :: digit

    write(digit,'(i0)') m
    if (m == 0) then
       named = 'f'
    else if (m == 1) then
       named = 'dfdx'
    else
       named = 'd' // trim(digit) // 'fdx' // trim(digit)
    end if

  end function order_named

  function right(text) result(cell)

    character(len=*), intent(in) :: text
    character(len=14) :: cell

    write(cell,'(a14)') text

  end function right

end program graph_time_integrator
