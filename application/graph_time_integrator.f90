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

  use iso_fortran_env       , only : dp => REAL64
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : implicit_midpoint, crouzeix_two_stage, &
       & crouzeix_three_stage
  use operation_grid        , only : uniform_grid, random_grid, designed_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use operation_grid        , only : grid
  use gti_march             , only : partitioned
  use gti_expansion         , only : family_holder
  use gti_chain             , only : chain_block, march_chain, chain_expansion, &
       & instant_components
  use gti_sweeps            , only : set_krylov_above
  use gti_configuration     , only : configuration, read_configuration, override, show

  implicit none

  ! The residual at or below which a march has converged. A row above
  ! it is reported unconverged and carries no derivative, the state a
  ! gradient would be taken at not having been reached.
  real(dp), parameter :: converged_residual = 1.0e-8_dp

  type(configuration) :: cfg

  call settings(cfg)
  call show(cfg)
  call set_krylov_above(cfg % krylov_above)
  call table(cfg)

contains

  !-------------------------------------------------------------------!
  ! The configuration named on the command line, then every other
  ! argument applied over it.
  !-------------------------------------------------------------------!

  subroutine settings(cfg)

    type(configuration), intent(out) :: cfg

    character(len=256) :: argument
    integer :: i

    call read_configuration(named(), cfg)

    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(argument, '--config=') == 1) cycle
       call override(cfg, argument)
    end do

  end subroutine settings

  function named() result(name)

    character(len=:), allocatable :: name

    character(len=256) :: argument
    integer :: i

    name = 'homogeneous'

    do i = 1, command_argument_count()
       call get_command_argument(i, argument)
       if (index(argument, '--config=') == 1) name = trim(argument(10:))
    end do

  end function named

  !-------------------------------------------------------------------!
  ! The one instant a stage family needs, and it is consistent with
  ! the equation rather than merely plausible: the value and every
  ! derivative below the highest are chosen, and the highest is what
  ! the governing constraint then requires.
  !-------------------------------------------------------------------!

  pure function at_rest(cfg) result(q)

    type(configuration), intent(in) :: cfg
    real(dp), allocatable :: q(:)

    integer :: nd

    nd = cfg % state_degree + 1
    allocate(q(nd), source=0.0_dp)

    q(1)  =  1.0_dp
    q(nd) = -1.0_dp

  end function at_rest

  !-------------------------------------------------------------------!
  ! How far back a family of one order looks, which for a backward
  ! difference grows with the degree of the equation as well as with
  ! the order, and for the others does not.
  !-------------------------------------------------------------------!

  pure integer function reach_of(cfg, name, order) result(reach)

    type(configuration), intent(in) :: cfg
    character(len=*)   , intent(in) :: name
    integer            , intent(in) :: order

    select case (name)
    case ('bdf')
       reach = cfg % state_degree * order
    case ('adams')
       reach = max(order - 1, 1)
    case default
       reach = 1
    end select

  end function reach_of

  !-------------------------------------------------------------------!
  ! How many instants the widest row that fits looks back over. A row
  ! that reaches past the horizon is not built, so it does not decide
  ! how long a startup the others need; zero means none of them fit.
  !-------------------------------------------------------------------!

  pure integer function widest_reach(cfg) result(widest)

    type(configuration), intent(in) :: cfg

    character(len=8) :: every(3)
    integer :: i, order, reach

    every  = ['bdf     ', 'adams   ', 'dirk    ']
    widest = 0

    do i = 1, 3
       if (index(cfg % families, trim(every(i))) == 0) cycle
       do order = 1, cfg % max_discretization_order
          reach = reach_of(cfg, trim(every(i)), order)
          if (reach < cfg % instants) widest = max(widest, reach)
       end do
    end do

  end function widest_reach

  !-------------------------------------------------------------------!
  ! The instants every row starts from, integrated rather than
  ! invented: a stage family over the startup, on a grid refined
  ! within each of its steps, sampled back at the coarse instants.
  !-------------------------------------------------------------------!

  subroutine startup_trajectory(cfg, dt, widest, held)

    type(configuration), intent(in)  :: cfg
    real(dp)           , intent(in)  :: dt(:)
    integer            , intent(in)  :: widest
    real(dp), allocatable, intent(out) :: held(:)

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    integer , allocatable :: added(:)
    real(dp), allocatable :: fine_dt(:), fine_t(:), substeps(:)
    real(dp) :: achieved, span
    integer :: nd, k, r, b

    nd = cfg % state_degree + 1

    if (widest == 1) then
       held = at_rest(cfg)
       return
    end if

    r        = max(cfg % startup_refinement, 1)
    substeps = [(dt(1 + (k - 1) / r + 1) / real(r, dp), k = 1, (widest - 1) * r)]
    span     = sum(substeps)

    added = in_pieces((widest - 1) * r + 1)
    allocate(schemes(size(added)))
    do b = 1, size(added)
       allocate(schemes(b) % scheme, source=crouzeix_three_stage())
    end do

    call march_chain(schemes, added, van_der_pol(cfg % state_degree), &
         & nd, designed_grid(span), cfg % design, at_rest(cfg), chain, &
         & fine_dt, fine_t, achieved, grid_design = substeps)

    call sampled(chain, nd, widest, r, held)

  end subroutine startup_trajectory

  !-------------------------------------------------------------------!
  ! A count of instants split into blocks short enough that each is
  ! cheap to solve. Every block after the first shares one instant
  ! with the one before it, which is what a stage family looks back
  ! over, so the pieces add to the whole.
  !-------------------------------------------------------------------!

  pure function in_pieces(instants) result(added)

    integer, intent(in) :: instants
    integer, allocatable :: added(:)

    integer, parameter :: piece = 4
    integer :: blocks

    if (instants <= piece + 1) then
       added = [instants]
       return
    end if

    blocks = (instants - 1) / piece
    allocate(added(blocks))

    added    = piece
    added(1) = instants - piece * (blocks - 1)

  end function in_pieces

  !-------------------------------------------------------------------!
  ! The refined trajectory read back at the coarse instants, which
  ! are every r-th instant of it.
  !-------------------------------------------------------------------!

  subroutine sampled(chain, nd, widest, r, held)

    type(chain_block), intent(in)  :: chain(:)
    integer          , intent(in)  :: nd, widest, r
    real(dp), allocatable, intent(out) :: held(:)

    integer :: k

    allocate(held(widest * nd))

    do k = 1, widest
       held((k - 1) * nd + 1:k * nd) = &
            & instant_components(chain, 1 + (k - 1) * r, nd)
    end do

  end subroutine sampled

  !-------------------------------------------------------------------!
  ! One family, by name and order. A stage family says that it is
  ! one, since its block is laid out differently.
  !-------------------------------------------------------------------!

  subroutine chosen(name, order, scheme, staged, ok)

    character(len=*), intent(in)  :: name
    integer         , intent(in)  :: order
    class(family), allocatable, intent(out) :: scheme
    logical         , intent(out) :: staged, ok

    staged = .false.
    ok     = .true.

    select case (name)
    case ('bdf')
       allocate(scheme, source=bdf_family(order))
    case ('adams')
       allocate(scheme, source=adams_family(order))
    case ('dirk')
       staged = .true.
       select case (order)
       case (2)
          allocate(scheme, source=implicit_midpoint())
       case (3)
          allocate(scheme, source=crouzeix_two_stage())
       case (4)
          allocate(scheme, source=crouzeix_three_stage())
       case default
          ok = .false.
       end select
    case default
       ok = .false.
    end select

  end subroutine chosen


  subroutine steps_of(cfg, dt, t)

    type(configuration), intent(in) :: cfg
    real(dp), allocatable, intent(out) :: dt(:), t(:)

    select case (trim(cfg % grid))
    case ('uniform')
       call partitioned(uniform_grid(cfg % time_duration), cfg % instants, dt, t)
    case ('random')
       call partitioned(random_grid(cfg % time_duration, cfg % seed), cfg % instants, dt, t)
    case default
       error stop 'graph_time_integrator: a grid is uniform or random'
    end select

  end subroutine steps_of

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

  subroutine heading(cfg)

    type(configuration), intent(in) :: cfg

    character(len=:), allocatable :: line
    character(len=2) :: digit
    integer :: m

    line = '  scheme' // repeat(' ', 14) // 'solved' // repeat(' ', 8) // 'f'

    do m = 1, cfg % max_derivative_degree
       write(digit,'(i0)') m
       if (m == 1) then
          line = line // repeat(' ', 12) // 'dfdx'
       else
          line = line // repeat(' ', 9) // 'd' // trim(digit) // 'fdx' // trim(digit)
       end if
    end do

    write(*,'(a)') ' '
    write(*,'(a)') line

  end subroutine heading

  subroutine show_row(label, solved, f, achieved, columns)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: solved
    real(dp)        , intent(in) :: f(0:), achieved
    integer         , intent(in) :: columns

    character(len=16) :: cell
    character(len=:), allocatable :: line
    integer :: m

    line = '  ' // label // repeat(' ', max(2, 20 - len(label)))
    write(cell,'(i6)') solved
    line = line // cell

    ! A column the row holds no expansion for is left empty rather
    ! than filled, there being no number to state under it.
    do m = 0, columns
       if (m <= ubound(f, 1)) then
          write(cell,'(es15.6)') f(m)
       else
          write(cell,'(a15)') '-'
       end if
       line = line // cell
    end do

    if (achieved > converged_residual) line = line // '   unconverged'

    write(*,'(a)') line

  end subroutine show_row

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

    call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
         & chosen_grid(cfg), cfg % design, held, chain, dt, t, achieved)

    ! Every derivative is taken at the state the march reached, so a
    ! row that did not converge has none to take and only its value is
    ! expanded.
    if (achieved > converged_residual) then
       reported = 0
    else
       reported = cfg % max_derivative_degree
    end if

    call chain_expansion(chain, van_der_pol(cfg % state_degree), &
         & van_der_pol_energy(cfg % state_degree), nd, dt, cfg % design, &
         & reported, f)

    call show_row(labelled(names, orders), cfg % instants - given, f, achieved, &
         & cfg % max_derivative_degree)
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

  function chosen_grid(cfg) result(steps)

    type(configuration), intent(in) :: cfg
    class(grid), allocatable :: steps

    select case (trim(cfg % grid))
    case ('uniform')
       allocate(steps, source=uniform_grid(cfg % time_duration))
    case ('random')
       allocate(steps, source=random_grid(cfg % time_duration, cfg % seed))
    case default
       error stop 'graph_time_integrator: a grid is uniform or random'
    end select

  end function chosen_grid

  !-------------------------------------------------------------------!
  ! Every row the configuration asks for.
  !-------------------------------------------------------------------!

  subroutine table(cfg)

    type(configuration), intent(in) :: cfg

    real(dp), allocatable :: startup(:), dt(:), t(:)
    integer :: widest, printed

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

    call heading(cfg)

    printed = 0
    if (asked(cfg, 'homogeneous')) call homogeneous_rows(cfg, startup, printed)
    if (asked(cfg, 'pairs'))       call pair_rows(cfg, startup, printed)
    if (asked(cfg, 'triples'))     call triple_rows(cfg, startup, printed)

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

    yes = index(cfg % combinations, what) > 0

  end function asked

  !-------------------------------------------------------------------!
  ! The names a configuration lists, in the order it lists them.
  !-------------------------------------------------------------------!

  function listed(cfg) result(names)

    type(configuration), intent(in) :: cfg
    character(len=8), allocatable :: names(:)

    character(len=8) :: every(3)
    integer :: i, n

    every = ['bdf     ', 'adams   ', 'dirk    ']
    n = 0
    do i = 1, 3
       if (index(cfg % families, trim(every(i))) > 0) n = n + 1
    end do

    allocate(names(n))
    n = 0
    do i = 1, 3
       if (index(cfg % families, trim(every(i))) > 0) then
          n = n + 1
          names(n) = every(i)
       end if
    end do

  end function listed

  subroutine homogeneous_rows(cfg, startup, printed)

    type(configuration), intent(in) :: cfg
    real(dp)           , intent(in) :: startup(:)
    integer            , intent(inout) :: printed

    character(len=8), allocatable :: names(:)
    integer :: i, order

    names = listed(cfg)

    do i = 1, size(names)
       do order = 1, cfg % max_discretization_order
          call one_row(cfg, startup, [names(i)], [order], printed)
       end do
    end do

  end subroutine homogeneous_rows

  !-------------------------------------------------------------------!
  ! Every ordered pair of distinct families, at one order or, when
  ! mixed orders are asked for, at every pair of orders.
  !-------------------------------------------------------------------!

  subroutine pair_rows(cfg, startup, printed)

    type(configuration), intent(in) :: cfg
    real(dp)           , intent(in) :: startup(:)
    integer            , intent(inout) :: printed

    character(len=8), allocatable :: names(:)
    integer :: i, j, p, q

    names = listed(cfg)

    do i = 1, size(names)
       do j = 1, size(names)
          if (i == j) cycle
          do p = 1, cfg % max_discretization_order
             if (cfg % mixed_orders) then
                do q = 1, cfg % max_discretization_order
                   call one_row(cfg, startup, [names(i), names(j)], [p, q], printed)
                end do
             else
                call one_row(cfg, startup, [names(i), names(j)], [p, p], printed)
             end if
          end do
       end do
    end do

  end subroutine pair_rows

  !-------------------------------------------------------------------!
  ! Every permutation of the families listed, at one order.
  !-------------------------------------------------------------------!

  subroutine triple_rows(cfg, startup, printed)

    type(configuration), intent(in) :: cfg
    real(dp)           , intent(in) :: startup(:)
    integer            , intent(inout) :: printed

    character(len=8), allocatable :: names(:)
    integer :: i, j, k, order

    names = listed(cfg)

    do i = 1, size(names)
       do j = 1, size(names)
          if (j == i) cycle
          do k = 1, size(names)
             if (k == i .or. k == j) cycle
             do order = 1, cfg % max_discretization_order
                call one_row(cfg, startup, [names(i), names(j), names(k)], [order, order, order], printed)
             end do
          end do
       end do
    end do

  end subroutine triple_rows

end program graph_time_integrator
