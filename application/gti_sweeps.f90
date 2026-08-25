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

  use util_precision  , only : dp
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use graph_fractal         , only : graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_stencil     , only : stencil
  use operation_dense_direct, only : dense_direct
  use operation_multigrid   , only : multigrid
  use operation_gauss_seidel, only : gauss_seidel
  use util_tally, only : tally_record, tangent_loops, adjoint_loops
  use util_factorisation, only : dense_factorisation
  use operation_gmres       , only : gmres
  use operation_minimization, only : minimizer
  use operation_linearization, only : linearization, tangent_of

  implicit none

  private
  public :: functional_of, functional_gradient
  public :: design_partial, jacobian_of
  public :: by_tangent, by_adjoint, tangent_solve
  public :: forward_route, reverse_route, route_of, route_substitutions

  integer, parameter :: forward_route = 1
  integer, parameter :: reverse_route = 2
  public :: linear_solver_named, set_linear_solver, set_aggregates, inner_minimizer
  public :: set_assembly, assembly_present

  !===================================================================!
  ! WHICH LINEAR SOLVER the tangent systems go to: named, not chosen
  ! by a count. dense factorises; gmres iterates without a matrix;
  ! multigrid smooths and detours to the aggregates it was given,
  ! which a field supplies from its mesh. A name this module has
  ! nothing for stops the program.
  !===================================================================!

  character(len=16), save :: chosen_solver = 'dense'
  integer, allocatable, save :: chosen_aggregates(:)

  !===================================================================!
  ! WHETHER A MATRIX IS ASSEMBLED. present: the statement's compiled
  ! tangent, a sparse stencil, from which a dense factorisation is
  ! formed where the solver is dense. free: no matrix anywhere - the
  ! linearization's matvec is what the inner minimizer sees, so it
  ! must iterate. A dense solver on a free assembly is refused.
  !===================================================================!

  character(len=16), save :: chosen_assembly = 'present'

contains

  subroutine set_assembly(name)

    character(len=*), intent(in) :: name

    select case (trim(name))
    case ('present', 'free')
       chosen_assembly = name
    case default
       write(*,'(a)') ' assembly names ' // trim(name) // ', which this program has nothing for.'
       error stop 'gti_sweeps: an assembly is present or free'
    end select

  end subroutine set_assembly

  pure logical function assembly_present() result(yes)

    yes = trim(chosen_assembly) == 'present'

  end function assembly_present

  pure function linear_solver_named() result(name)

    character(len=:), allocatable :: name

    name = trim(chosen_solver)

  end function linear_solver_named

  subroutine set_linear_solver(name)

    character(len=*), intent(in) :: name

    select case (trim(name))
    case ('dense', 'gmres', 'multigrid')
       chosen_solver = name
    case default
       write(*,'(a)') ' linear_solver names ' // trim(name) // &
            & ', which this program has nothing for.'
       error stop 'gti_sweeps: a linear solver is dense, gmres or multigrid'
    end select

  end subroutine set_linear_solver

  !===================================================================!
  ! The aggregates multigrid coarsens by: one block per unknown. A
  ! caller with no field clears them.
  !===================================================================!

  subroutine set_aggregates(aggregates)

    integer, intent(in), optional :: aggregates(:)

    if (allocated(chosen_aggregates)) deallocate(chosen_aggregates)
    if (present(aggregates)) chosen_aggregates = aggregates

  end subroutine set_aggregates

  !===================================================================!
  ! The minimizer named, built for a system of the given count. A
  ! singular pivot in the dense one is reported, the matrix being a
  ! tangent at an intermediate iterate. multigrid without aggregates
  ! of the right count stops the program.
  !===================================================================!

  function inner_minimizer(count) result(inner)

    integer, intent(in) :: count
    class(minimizer), allocatable :: inner

    type(gmres)        :: krylov
    type(dense_direct) :: factorisation
    type(multigrid)    :: levels
    type(gauss_seidel) :: sweeps

    if (.not. assembly_present() .and. trim(chosen_solver) /= 'gmres') then
       error stop 'gti_sweeps: a free assembly has no matrix to factorise; its solver iterates'
    end if

    select case (trim(chosen_solver))
    case ('dense')
       factorisation = dense_direct()
       factorisation % singular_reported = .true.
       allocate(inner, source=factorisation)
    case ('gmres')
       krylov = gmres()
       krylov % restart        = min(count, 60)
       krylov % tolerance      = 1.0e-13_dp
       krylov % max_iterations = 4
       allocate(inner, source=krylov)
    case ('multigrid')
       if (.not. allocated(chosen_aggregates)) then
          error stop 'gti_sweeps: multigrid coarsens by aggregates, and none were given'
       end if
       if (size(chosen_aggregates) /= count) then
          error stop 'gti_sweeps: one aggregate per unknown'
       end if
       sweeps % max_iterations = 2
       sweeps % tolerance      = 1.0e-13_dp
       allocate(levels % smoother, source=sweeps)
       factorisation = dense_direct()
       allocate(levels % coarse, source=factorisation)
       levels % aggregates     = chosen_aggregates
       levels % tolerance      = 1.0e-13_dp
       levels % max_iterations = 200
       allocate(inner, source=levels)
    end select

  end function inner_minimizer

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
  ! One solve against the statement's tangent in the state, frozen at
  ! the inputs given. No matrix is formed past a small block: the
  ! statement's own partial action is the matvec.
  !===================================================================!

  subroutine tangent_solve(rows, on, inputs, b, x)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: on
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp)             , intent(in) :: b(:)
    real(dp), allocatable, intent(out) :: x(:)

    type(linearization) :: jacobian
    type(stencil) :: compiled
    class(minimizer), allocatable :: solver
    integer , allocatable :: r(:), c(:)
    real(dp), allocatable :: w(:)
    logical :: available
    real(dp) :: achieved

    allocate(solver, source=inner_minimizer(size(b)))

    ! the statement's own compiled tangent where it offers one, the
    ! linearization otherwise
    available = .false.
    if (assembly_present()) call rows % compiled_tangent(on, inputs, 1, r, c, w, available)
    if (available) then
       compiled = stencil(r, c, w, spread(0.0_dp, 1, size(b)), 'compiled tangent')
       call solver % attach(compiled, compiled % pattern, on % vertex_set(), size(b), &
            & coupling = compiled % pattern)
    else
       jacobian = tangent_of(rows, rows % argument(1))
       call jacobian % freeze(inputs)
       call solver % attach(jacobian, on, on % vertex_set(), size(b))
    end if

    allocate(x(size(b)), source=0.0_dp)
    call solver % solve(b, x, achieved)

  end subroutine tangent_solve

  !===================================================================!
  ! The tangent: one solve in the state, then the gradient read
  ! along what it gives.
  !===================================================================!

  real(dp) function by_tangent(factor, g, design_rate, explicit) result(df)

    type(dense_factorisation), intent(in) :: factor
    real(dp)                 , intent(in) :: g(:), design_rate(:), explicit

    real(dp), allocatable :: w(:)

    call tally_record(tangent_loops)
    call factor % substitute(-design_rate, w, transposed=.false.)
    df = explicit + dot_product(g, w)

  end function by_tangent

  !===================================================================!
  ! The adjoint: one substitution against the transpose, then the
  ! statement's design partial read along what it gives.
  !===================================================================!

  real(dp) function by_adjoint(factor, g, design_rate, explicit) result(df)

    type(dense_factorisation), intent(in) :: factor
    real(dp)                 , intent(in) :: g(:), design_rate(:), explicit

    real(dp), allocatable :: lambda(:)

    call tally_record(adjoint_loops)
    call factor % substitute(g, lambda, transposed=.true.)
    df = explicit - dot_product(lambda, design_rate)

  end function by_adjoint

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
