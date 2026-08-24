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
  use gti_block             , only : block_residual
  use gti_march             , only : partitioned, block_of
  use gti_stage             , only : stage_block_of, instant_at
  use gti_taylor            , only : block_expansion
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

  !-------------------------------------------------------------------!
  ! One row: march the block and expand the functional over it.
  !-------------------------------------------------------------------!

  subroutine one_row(cfg, name, order)

    type(configuration), intent(in) :: cfg
    character(len=*)   , intent(in) :: name
    integer            , intent(in) :: order

    class(family), allocatable :: scheme
    type(block_residual) :: rows
    real(dp), allocatable :: dt(:), t(:), q(:), f(:)
    integer , allocatable :: at(:)
    real(dp) :: achieved
    integer :: nd
    logical :: staged, ok

    call chosen(name, order, scheme, staged, ok)
    if (.not. ok) return

    nd = cfg % state_degree + 1
    call steps_of(cfg, dt, t)
    if (scheme % history_depth() >= cfg % instants) return

    call block_and_instants(cfg, scheme, staged, nd, dt, t, rows, at)

    call block_expansion(rows, van_der_pol(cfg % state_degree), &
         & van_der_pol_energy(cfg % state_degree), nd, &
         & scheme % primary_degree(nd - 1), at, dt, cfg % design, &
         & cfg % max_derivative_degree, q, f, achieved)

    call show_row(labelled(name, order), f, achieved)

  end subroutine one_row

  !-------------------------------------------------------------------!
  ! The block a family makes, and where its instants sit among the
  ! unknowns: one set per instant for a multistep family, and between
  ! the stages for a stage family.
  !-------------------------------------------------------------------!

  subroutine block_and_instants(cfg, scheme, staged, nd, dt, t, rows, at)

    type(configuration) , intent(in)  :: cfg
    class(family)       , intent(in)  :: scheme
    logical             , intent(in)  :: staged
    integer             , intent(in)  :: nd
    real(dp)            , intent(in)  :: dt(:), t(:)
    type(block_residual), intent(out) :: rows
    integer, allocatable, intent(out) :: at(:)

    real(dp), allocatable :: held(:)
    integer :: k, d

    held = [((initial(d, t(k)), d = 0, nd - 1), k = 1, scheme % history_depth())]

    if (staged) then
       rows = stage_block_of(scheme, van_der_pol(cfg % state_degree), nd, &
            & cfg % instants, dt, held)
       at   = [(instant_at(k, scheme % num_stages(), nd), k = 1, cfg % instants)]
    else
       rows = block_of(scheme, van_der_pol(cfg % state_degree), nd, &
            & cfg % instants, dt, held)
       at   = [((k - 1) * nd, k = 1, cfg % instants)]
    end if

  end subroutine block_and_instants

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

  function labelled(name, order) result(text)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: order
    character(len=:), allocatable :: text

    character(len=2) :: digit

    write(digit,'(i0)') order
    text = trim(name) // trim(digit)

  end function labelled

  !-------------------------------------------------------------------!
  ! The heading, then one line per row.
  !-------------------------------------------------------------------!

  subroutine heading(cfg)

    type(configuration), intent(in) :: cfg

    character(len=:), allocatable :: line
    character(len=2) :: digit
    integer :: m

    line = '  scheme  ' // repeat(' ', 6) // 'f'

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

    line = '  ' // label // repeat(' ', max(1, 8 - len(label)))

    do m = 0, ubound(f, 1)
       write(cell,'(es15.6)') f(m)
       line = line // cell
    end do

    if (achieved > 1.0e-8_dp) line = line // '   unconverged'

    write(*,'(a)') line

  end subroutine show_row

  !-------------------------------------------------------------------!
  ! Every row the configuration asks for.
  !-------------------------------------------------------------------!

  subroutine table(cfg)

    type(configuration), intent(in) :: cfg

    integer :: order

    if (trim(cfg % combinations) /= 'homogeneous') then
       write(*,'(a)') ' '
       write(*,'(a)') ' only homogeneous rows are built: a horizon of two families needs'
       write(*,'(a)') ' the junction to map one block layout onto another, which is not here.'
       error stop 'graph_time_integrator: the combinations asked for are not built'
    end if

    call heading(cfg)

    if (index(cfg % families, 'bdf') > 0) then
       do order = 1, cfg % max_discretization_order
          call one_row(cfg, 'bdf', order)
       end do
    end if

    if (index(cfg % families, 'adams') > 0) then
       do order = 1, cfg % max_discretization_order
          call one_row(cfg, 'adams', order)
       end do
    end if

    if (index(cfg % families, 'dirk') > 0) then
       do order = 2, min(4, cfg % max_discretization_order)
          call one_row(cfg, 'dirk', order)
       end do
    end if

  end subroutine table

end program graph_time_integrator
