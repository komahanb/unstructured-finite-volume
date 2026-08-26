!=====================================================================!
! The sensitivity of a functional to a design, by the tangent and
! by the adjoint.
!
! A block's statement R(q, x) = 0 determines q from x, and the
! functional is a sum over the instants,
!
!      f  =  sum over k of  dt_k F(q_k, x) .
!
! Differentiating the statement gives J dq/dx = -dR/dx, so
!
!      tangent    solve J w = -dR/dx  once, then  df/dx = f_x + g.w
!      adjoint    solve J^T l = g     once, then  df/dx = f_x - l.dR/dx
!
! with J the jacobian in the state, g the gradient of f in the state
! and f_x its own partial in the design. The two read the same
! three objects and must agree to round-off; that they do is what
! this module exists to make checkable.
!
! Each of the three comes from a partial action that is exact: the
! scheme's rows are linear and are their own jacobian, and the
! physics and the integrand differentiate their own rules. Nothing
! here is differenced.
!
!             HOW THE JACOBIAN IS FORMED
!
! For the adjoint, column by column: one partial action per unknown,
! and then solved densely. That does not scale, and what would is
! the stencil's compiled transpose rather than a dense one.
!
! For the tangent no matrix is formed at all past a large block. The
! statement's partial action is already a matvec, so freezing it at
! the trajectory gives a linear operation a krylov solver can be
! driven with directly. Below the same threshold gti_march uses, a
! dense factorisation is the faster of the two and is taken.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_sweeps

  use gti_configuration, only : refuse_unknown
  use util_precision  , only : dp
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use graph_fractal         , only : graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_dense_direct, only : dense_direct
  use operation_multigrid   , only : multigrid
  use operation_gauss_seidel, only : gauss_seidel
  use operation_gmres       , only : gmres
  use operation_minimization, only : minimizer

  implicit none

  private
  public :: functional_of, functional_gradient
  public :: design_partial, jacobian_of
  public :: forward_route, reverse_route, route_of, route_substitutions

  integer, parameter :: forward_route = 1
  integer, parameter :: reverse_route = 2
  public :: set_linear_solver, set_assembly, set_storage, set_multigrid
  public :: set_aggregates, aggregates_given, aggregates_of, assembly_present, multigrid_on
  public :: take_inner, keep_inner, forget_inner

  !===================================================================!
  ! HOW A LINEAR SYSTEM IS SOLVED: four specifications, each its own.
  !
  !      linear_solver   direct | iterative    factorise, or iterate
  !      assembly        matrix | free         a matrix is formed, or
  !                                            only a matvec is attached
  !      storage         dense | sparse        the matrix formed is a
  !                                            square, or a stencil
  !      multigrid       yes | no              a two-grid over aggregates:
  !                                            block gauss-seidel over a
  !                                            point's components smooths,
  !                                            the solver named serves
  !                                            the coarse level
  !
  ! The corners that mean nothing are refused by name: direct on a
  ! free assembly, sparse direct (not built), dense iterative, and
  ! multigrid on a free assembly, its coarse operator being read
  ! through the aggregates from a stencil.
  !===================================================================!

  character(len=16), save :: chosen_solver   = 'direct'
  character(len=16), save :: chosen_assembly = 'matrix'
  character(len=16), save :: chosen_storage  = 'dense'
  logical          , save :: chosen_multigrid = .false.
  integer, allocatable, save :: chosen_aggregates(:)

  ! The inner minimizer kept between solves, so that a direct one
  ! keeps its factors across the statements stamped alike.
  class(minimizer), allocatable, save :: kept_inner

contains

  subroutine set_linear_solver(name)

    character(len=*), intent(in) :: name

    call refuse_unknown(name, ['direct   ', 'iterative'], 'linear_solver')
    chosen_solver = name
    call forget_inner()

  end subroutine set_linear_solver

  subroutine set_assembly(name)

    character(len=*), intent(in) :: name

    call refuse_unknown(name, ['matrix', 'free  '], 'assembly')
    chosen_assembly = name
    call forget_inner()

  end subroutine set_assembly

  subroutine set_storage(name)

    character(len=*), intent(in) :: name

    call refuse_unknown(name, ['dense ', 'sparse'], 'storage')
    chosen_storage = name
    call forget_inner()

  end subroutine set_storage

  subroutine set_multigrid(on)

    logical, intent(in) :: on

    chosen_multigrid = on
    call forget_inner()

  end subroutine set_multigrid

  pure logical function assembly_present() result(yes)

    yes = trim(chosen_assembly) == 'matrix'

  end function assembly_present

  pure logical function multigrid_on() result(yes)

    yes = chosen_multigrid

  end function multigrid_on

  !===================================================================!
  ! The aggregates multigrid coarsens by: one block per unknown. A
  ! caller with no field clears them.
  !===================================================================!

  subroutine set_aggregates(aggregates)

    integer, intent(in), optional :: aggregates(:)

    if (allocated(chosen_aggregates)) deallocate(chosen_aggregates)
    if (present(aggregates)) chosen_aggregates = aggregates

  end subroutine set_aggregates

  subroutine aggregates_of(aggregates)

    integer, allocatable, intent(out) :: aggregates(:)

    if (allocated(chosen_aggregates)) aggregates = chosen_aggregates

  end subroutine aggregates_of

  pure logical function aggregates_given() result(yes)

    yes = allocated(chosen_aggregates)

  end function aggregates_given

  !===================================================================!
  ! The minimizer the specifications name, built for a system of the
  ! given count. A singular pivot in a direct one is reported, the
  ! matrix being a tangent at an intermediate iterate.
  !===================================================================!

  function inner_minimizer(count, width) result(inner)

    integer, intent(in) :: count, width
    class(minimizer), allocatable :: inner

    class(minimizer), allocatable :: named
    type(gmres)        :: krylov
    type(dense_direct) :: factorisation
    type(multigrid)    :: levels
    type(gauss_seidel) :: sweeps

    if (trim(chosen_assembly) == 'free' .and. trim(chosen_solver) == 'direct') then
       error stop 'gti_sweeps: a free assembly has no matrix to factorise; its solver iterates'
    end if
    if (trim(chosen_solver) == 'direct' .and. trim(chosen_storage) == 'sparse') then
       error stop 'gti_sweeps: a sparse direct solve is not built'
    end if
    if (trim(chosen_solver) == 'iterative' .and. trim(chosen_assembly) == 'matrix' &
         & .and. trim(chosen_storage) == 'dense') then
       error stop 'gti_sweeps: an iterative solve reads the sparse stencil; dense storage is for factorising'
    end if
    if (chosen_multigrid .and. trim(chosen_assembly) == 'free') then
       error stop 'gti_sweeps: multigrid coarsens a stencil, which a free assembly has not'
    end if

    select case (trim(chosen_solver))
    case ('direct')
       factorisation = dense_direct()
       factorisation % singular_reported = .true.
       allocate(named, source=factorisation)
    case ('iterative')
       krylov = gmres()
       krylov % restart        = min(count, 60)
       krylov % tolerance      = 1.0e-13_dp
       krylov % max_iterations = 4
       allocate(named, source=krylov)
    end select

    if (.not. chosen_multigrid) then
       call move_alloc(named, inner)
       return
    end if

    if (.not. allocated(chosen_aggregates)) then
       error stop 'gti_sweeps: multigrid coarsens by aggregates, and none were given'
    end if
    if (size(chosen_aggregates) /= count) then
       error stop 'gti_sweeps: one aggregate per unknown'
    end if

    ! gauss-seidel smooths, a point's components at a time - the
    ! coupling within a point being what no point smoother damps -
    ! and the solver named serves the coarse level
    sweeps % max_iterations = 2
    sweeps % block_width    = width
    allocate(levels % smoother, source=sweeps)
    call move_alloc(named, levels % coarse)
    levels % block_width    = width
    levels % aggregates     = chosen_aggregates
    levels % tolerance      = 1.0e-13_dp
    levels % max_iterations = 200
    allocate(inner, source=levels)

  end function inner_minimizer

  !===================================================================!
  ! The inner minimizer taken for a solve and kept after it, so that
  ! what it holds - a direct solver's factors - outlives one solve.
  !===================================================================!

  subroutine take_inner(inner, count, width)

    class(minimizer), allocatable, intent(out) :: inner
    integer                      , intent(in)  :: count, width

    ! multigrid is built afresh for every statement, its aggregates
    ! being the statement's; a kept one would carry another's
    if (allocated(kept_inner) .and. .not. chosen_multigrid) then
       call move_alloc(kept_inner, inner)
    else
       call forget_inner()
       allocate(inner, source=inner_minimizer(count, width))
    end if

  end subroutine take_inner

  subroutine keep_inner(inner)

    class(minimizer), allocatable, intent(inout) :: inner

    if (allocated(kept_inner)) deallocate(kept_inner)
    call move_alloc(inner, kept_inner)

  end subroutine keep_inner

  subroutine forget_inner()

    if (allocated(kept_inner)) deallocate(kept_inner)

  end subroutine forget_inner

  subroutine applied(action, on, inputs, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp), allocatable, intent(out) :: y(:)

    class(field), allocatable :: out

    call action % apply(on, inputs, out)
    call out % real_vector(y)

  end subroutine applied

  !===================================================================!
  ! The partial action of a statement along one direction in one of
  ! its arguments.
  !===================================================================!

  subroutine varied(action, on, inputs, which, domain, v, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: which
    type(graph)          , intent(in) :: domain
    real(dp)             , intent(in) :: v(:)
    real(dp), allocatable, intent(out) :: y(:)

    type(stored_field) :: direction
    class(field), allocatable :: out

    direction = stored_field('direction', domain, size(v))
    call direction % set_real_vector(v)

    call action % partial_action(on, inputs, &
         & [variation(action % argument(which), direction)], out)
    call out % real_vector(y)

  end subroutine varied

  !===================================================================!
  ! The functional: the integrand at every instant, weighted by the
  ! step that ends there. The first instant carries no step and so
  ! contributes nothing.
  !===================================================================!

  real(dp) function functional_of(integrand, instants, inputs, dt) result(f)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)

    real(dp), allocatable :: values(:)

    call applied(integrand, instants, inputs, values)
    f = sum(dt * values)

  end function functional_of

  !===================================================================!
  ! Its gradient in the state. The integrand at one instant reads
  ! only that instant's components, so one partial action per degree
  ! gives the whole gradient rather than one per unknown.
  !===================================================================!

  subroutine functional_gradient(integrand, instants, inputs, dt, n, degrees, &
       & state_domain, g)

    class(operation)     , intent(in) :: integrand
    class(directed_graph), intent(in) :: instants
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: dt(:)
    integer              , intent(in) :: n, degrees
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: g(:)

    real(dp), allocatable :: v(:), rate(:)
    integer :: d, k

    allocate(g(n * degrees), source=0.0_dp)
    allocate(v(n * degrees))

    do d = 0, degrees - 1
       v = 0.0_dp
       do k = 1, n
          v((k - 1) * degrees + d + 1) = 1.0_dp
       end do

       call varied(integrand, instants, inputs, 1, state_domain, v, rate)

       do k = 1, n
          g((k - 1) * degrees + d + 1) = dt(k) * rate(k)
       end do
    end do

  end subroutine functional_gradient

  !===================================================================!
  ! The statement's partial in the design, along the direction that
  ! varies every instant's design value together, which is a single
  ! design number for the whole block.
  !===================================================================!

  subroutine design_partial(rows, unknowns, inputs, n, design_domain, d)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: n
    type(graph)          , intent(in) :: design_domain
    real(dp), allocatable, intent(out) :: d(:)

    call varied(rows, unknowns, inputs, 2, design_domain, spread(1.0_dp, 1, n), d)

  end subroutine design_partial

  !===================================================================!
  ! The jacobian in the state, one column per unknown.
  !===================================================================!

  subroutine jacobian_of(rows, unknowns, inputs, num_unknowns, state_domain, a)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: num_unknowns
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: a(:,:)

    real(dp), allocatable :: v(:), column(:), w(:)
    integer , allocatable :: r(:), c(:)
    logical :: available
    integer :: j, e

    allocate(a(num_unknowns, num_unknowns), source=0.0_dp)

    ! A statement that compiles its tangent hands over its triples,
    ! and the square is filled from them; any other is probed one
    ! column at a time, a partial action each.
    call rows % compiled_tangent(unknowns, inputs, 1, r, c, w, available)
    if (available) then
       do e = 1, size(r)
          a(r(e), c(e)) = a(r(e), c(e)) + w(e)
       end do
       return
    end if

    allocate(v(num_unknowns), source=0.0_dp)
    do j = 1, num_unknowns
       v    = 0.0_dp
       v(j) = 1.0_dp
       call varied(rows, unknowns, inputs, 1, state_domain, v, column)
       a(:, j) = column
    end do

  end subroutine jacobian_of

  !===================================================================!
  ! THE GATE. Which route computes the m-th derivative of n_f
  ! functionals in n_d designs, by the count of substitutions each
  ! costs against one kept factorisation:
  !
  !      forward, m times            C(n_d + m - 1, m)
  !      forward m-1 times over      (1 + n_f) C(n_d + m - 2, m - 1)
  !      one reverse
  !
  ! The ratio of the first to the second is (n_d + m - 1) / (m (1 + n_f)),
  ! so the forward route is the cheaper exactly where n_d <= m n_f. At
  ! m = 1 that is the rule for a gradient, n_d <= n_f. At one design
  ! and one functional the forward route costs m substitutions and the
  ! other 2 m, at every order.
  !
  ! A count below one, or an order below one, stops the program.
  !===================================================================!

  pure integer function route_of(num_designs, num_functionals, order) result(route)

    integer, intent(in) :: num_designs, num_functionals, order

    if (num_designs < 1 .or. num_functionals < 1 .or. order < 1) then
       error stop 'gti_sweeps: a route is chosen for at least one design, one functional and order one'
    end if

    if (num_designs <= order * num_functionals) then
       route = forward_route
    else
       route = reverse_route
    end if

  end function route_of

  !===================================================================!
  ! The substitutions a route costs at one order, per block. What the
  ! accounting layer counts is compared against this.
  !===================================================================!

  pure integer function route_substitutions(route, num_designs, num_functionals, order) &
       & result(count)

    integer, intent(in) :: route, num_designs, num_functionals, order

    select case (route)
    case (forward_route)
       count = choose(num_designs + order - 1, order)
    case (reverse_route)
       count = (1 + num_functionals) * choose(num_designs + order - 2, order - 1)
    case default
       error stop 'gti_sweeps: a route is forward or reverse'
    end select

  end function route_substitutions

  pure integer function choose(n, k) result(c)

    integer, intent(in) :: n, k

    integer :: i

    c = 1
    do i = 1, k
       c = c * (n - k + i) / i
    end do

  end function choose

end module gti_sweeps
