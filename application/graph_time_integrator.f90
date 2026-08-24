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
program graph_time_integrator

  use iso_fortran_env       , only : dp => REAL64
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : implicit_midpoint, crouzeix_two_stage, &
       & crouzeix_three_stage
  use operation_grid        , only : uniform_grid, random_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use operation_grid        , only : grid
  use gti_march             , only : partitioned
  use gti_expansion         , only : family_holder
  use gti_chain             , only : chain_block, march_chain, chain_expansion
  use gti_configuration     , only : configuration, read_configuration, override, show

  implicit none

  type(configuration) :: cfg

  call settings(cfg)
  call show(cfg)
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
  ! The d-th derivative of the cosine, which is what every row starts
  ! from, so that the rows are comparable.
  !-------------------------------------------------------------------!

  pure real(dp) function initial(d, t) result(q)

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

  end function initial

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

    line = '  scheme' // repeat(' ', 20) // 'f'

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

  subroutine show_row(label, f, achieved)

    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: f(0:), achieved

    character(len=16) :: cell
    character(len=:), allocatable :: line
    integer :: m

    line = '  ' // label // repeat(' ', max(2, 20 - len(label)))

    do m = 0, ubound(f, 1)
       write(cell,'(es15.6)') f(m)
       line = line // cell
    end do

    if (achieved > 1.0e-8_dp) line = line // '   unconverged'

    write(*,'(a)') line

  end subroutine show_row

  !-------------------------------------------------------------------!
  ! One row: a chain of blocks, marched and then expanded. A chain of
  ! one is a homogeneous row and takes the same path.
  !-------------------------------------------------------------------!

  subroutine one_row(cfg, names, orders)

    type(configuration), intent(in) :: cfg
    character(len=*)   , intent(in) :: names(:)
    integer            , intent(in) :: orders(:)

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    integer , allocatable :: added(:)
    real(dp), allocatable :: dt(:), t(:), held(:), f(:)
    real(dp) :: achieved
    integer :: b, nd, k, d, given
    logical :: ok

    nd = cfg % state_degree + 1
    allocate(schemes(size(names)), added(size(names)))
    call assembled(cfg, names, orders, schemes, added, ok)
    if (.not. ok) return

    call steps_of(cfg, dt, t)
    given = schemes(1) % scheme % history_depth()
    held  = [((initial(d, t(k)), d = 0, nd - 1), k = 1, given)]

    call march_chain(schemes, added, van_der_pol(cfg % state_degree), nd, &
         & chosen_grid(cfg), cfg % design, held, chain, dt, t, achieved)

    call chain_expansion(chain, van_der_pol(cfg % state_degree), &
         & van_der_pol_energy(cfg % state_degree), nd, dt, cfg % design, &
         & cfg % max_derivative_degree, f)

    call show_row(labelled(names, orders), f, achieved)

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
       if (added(b) <= schemes(b) % scheme % history_depth()) ok = .false.
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

    call heading(cfg)

    if (asked(cfg, 'homogeneous')) call homogeneous_rows(cfg)
    if (asked(cfg, 'pairs'))       call pair_rows(cfg)
    if (asked(cfg, 'triples'))     call triple_rows(cfg)

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

  subroutine homogeneous_rows(cfg)

    type(configuration), intent(in) :: cfg

    character(len=8), allocatable :: names(:)
    integer :: i, order

    names = listed(cfg)

    do i = 1, size(names)
       do order = 1, cfg % max_discretization_order
          call one_row(cfg, [names(i)], [order])
       end do
    end do

  end subroutine homogeneous_rows

  !-------------------------------------------------------------------!
  ! Every ordered pair of distinct families, at one order or, when
  ! mixed orders are asked for, at every pair of orders.
  !-------------------------------------------------------------------!

  subroutine pair_rows(cfg)

    type(configuration), intent(in) :: cfg

    character(len=8), allocatable :: names(:)
    integer :: i, j, p, q

    names = listed(cfg)

    do i = 1, size(names)
       do j = 1, size(names)
          if (i == j) cycle
          do p = 1, cfg % max_discretization_order
             if (cfg % mixed_orders) then
                do q = 1, cfg % max_discretization_order
                   call one_row(cfg, [names(i), names(j)], [p, q])
                end do
             else
                call one_row(cfg, [names(i), names(j)], [p, p])
             end if
          end do
       end do
    end do

  end subroutine pair_rows

  !-------------------------------------------------------------------!
  ! Every permutation of the families listed, at one order.
  !-------------------------------------------------------------------!

  subroutine triple_rows(cfg)

    type(configuration), intent(in) :: cfg

    character(len=8), allocatable :: names(:)
    integer :: i, j, k, order

    names = listed(cfg)

    do i = 1, size(names)
       do j = 1, size(names)
          if (j == i) cycle
          do k = 1, size(names)
             if (k == i .or. k == j) cycle
             do order = 1, cfg % max_discretization_order
                call one_row(cfg, [names(i), names(j), names(k)], [order, order, order])
             end do
          end do
       end do
    end do

  end subroutine triple_rows

end program graph_time_integrator
