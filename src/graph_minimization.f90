!=====================================================================!
! LEVEL 2 OF THE STRATIFICATION . THE MINIMIZATION
!
! The first level with a goal: drive a residual toward zero. It holds
! TWO citizens and no abstraction of its own, because its
! abstractions are the grammar's - and the grammar has exactly two
! nouns, so there are exactly two things a minimizer can vary:
!
!      fit_minimizer    varies the VALUES     - an operation
!      form_minimizer   varies the STRUCTURE  - a transform
!
! and no third is possible, because nothing else exists to vary.
!
!=====================================================================!
!
!                        ONE ACT, WRITTEN ONCE
!
! Every iterative method ever named is the same step: choose a
! correction space, choose a constraint, move.
!
!      x  <-  x + V y        subject to     W' R(x + V y) = 0
!
! V spans the directions a correction may use; W says which parts of
! the residual must vanish. Nothing else is going on. What separates
! jacobi from gmres is not the act but WHERE THEY STAND in a product
! of four finite choices:
!
!      M         how a direction is made from a residual
!      memory    how many past directions are kept
!      step      how far to go: a fixed length, or the length the
!                constraint picks
!      tangent   whether the operator is re-read at the new state
!
! and every famous name is one point in that product:
!
!   method        M              memory  step        tangent
!   -----------------------------------------------------------------
!   jacobi        diagonal       0       fixed w=1   no
!   gauss-seidel  colour         0       fixed w     no      (sor: w)
!   conj gradient any            full    galerkin    no
!   gmres         any            k       min resid   no
!   multigrid     a cycle        0       fixed 1     no
!   substitution  causal block   0       fixed 1     no
!   newton        any            any     any         yes
!
! So the seven are ONE type with four arguments, and the names are
! CONSTRUCTORS - the module is the catalogue and the compiler is its
! index. This is the admission law's absorption rule applied where it
! bites hardest, and it is the same shape the engine already uses for
! its eight reduction rules, four walks, three partition rules and
! two coarsening rules: a new method costs a case, never a class.
!
! Multigrid is not a kind of solver; it is a choice of M. So is
! substitution, which is the M you get when the pattern is
! triangular. Newton is not a kind of solver; it is the last column.
!
!=====================================================================!
!
!                     THE INHERITED VOCABULARY
!
! The solver words are thin delegations to engine seats, so the
! engine never learns a solver word and the solver never says apply
! or measure:
!
!      matvec ········· the operation applied, minus its constant
!      inner_product ·· a sum reduction with the second field as
!                       the measure
!      norm ··········· the norm reduction
!      sweep_order ···· the colouring walk
!      diagonal ······· matvec probed by colour: applied to one
!                       colour's indicator, the answer at a member
!                       IS its diagonal entry, because no two
!                       neighbours share a colour
!      constant ······· the affine part of the attached operation
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_minimization

  use iso_fortran_env       , only : dp => REAL64
  use graph_grammar         , only : graph, graph_field, graph_operation
  use graph_grammar         , only : graph_transform
  use graph_calculus        , only : GRAPH_SIDE_VERTEX, graph_functional
  use class_graph_support   , only : support
  use class_graph_field     , only : field
  use class_graph_stencil   , only : stencil_operator
  use class_graph_linearization, only : difference_linearization
  use class_graph_reduction , only : reduction, REDUCE_SUM, REDUCE_NORM
  use class_graph_walk      , only : walk, WALK_COLOURING
  use graph_forms           , only : form

  implicit none

  private
  public :: fit_minimizer, form_minimizer

  ! the named points of the product space
  public :: jacobi, gauss_seidel, conjugate_gradient, gmres
  public :: multigrid, newton
  public :: pruner

  ! the axes themselves, for a caller standing somewhere unnamed
  public :: PRECONDITION_NONE, PRECONDITION_DIAGONAL
  public :: PRECONDITION_COLOUR, PRECONDITION_CYCLE
  public :: MEMORY_NONE, MEMORY_FULL
  public :: STEP_FIXED, STEP_GALERKIN, STEP_MINIMUM_RESIDUAL
  public :: FORM_PRUNE

  !-------------------------------------------------------------------!
  ! M: how a direction is made from a residual.
  !-------------------------------------------------------------------!

  integer, parameter :: PRECONDITION_NONE     = 1   ! take it as it is
  integer, parameter :: PRECONDITION_DIAGONAL = 2   ! divide by the diagonal
  integer, parameter :: PRECONDITION_COLOUR   = 3   ! sweep, colour by colour
  integer, parameter :: PRECONDITION_CYCLE    = 4   ! smooth, coarsen, correct

  !-------------------------------------------------------------------!
  ! Memory: how many past directions are kept. Zero is a relaxation,
  ! which remembers nothing; full is a krylov space, which remembers
  ! everything; k between them restarts.
  !-------------------------------------------------------------------!

  integer, parameter :: MEMORY_NONE = 0
  integer, parameter :: MEMORY_FULL = -1

  !-------------------------------------------------------------------!
  ! The step: a fixed length, or the one a constraint picks.
  !-------------------------------------------------------------------!

  integer, parameter :: STEP_FIXED            = 1
  integer, parameter :: STEP_GALERKIN         = 2   ! W = V
  integer, parameter :: STEP_MINIMUM_RESIDUAL = 3   ! W = A V

  !-------------------------------------------------------------------!
  ! What a form minimizer does to a form.
  !-------------------------------------------------------------------!

  integer, parameter :: FORM_PRUNE = 1

  !===================================================================!
  ! THE VALUE SIDE. One type, standing wherever its four arguments
  ! put it.
  !===================================================================!

  type, extends(graph_operation) :: fit_minimizer

     !----------------------------------------------------------------!
     ! The attached statement and the seats its answers live on.
     !----------------------------------------------------------------!

     class(graph_operation), allocatable :: action
     class(graph)          , allocatable :: on

     type(support) :: cells
     type(support) :: numbers

     integer :: ncomp = 1

     real(dp), allocatable :: affine(:)

     integer  :: max_iterations = 1000
     real(dp) :: tolerance      = 1.0d-10

     !----------------------------------------------------------------!
     ! Where this one stands in the product space.
     !----------------------------------------------------------------!

     integer  :: precondition = PRECONDITION_NONE
     real(dp) :: relaxation   = 1.0_dp
     integer  :: memory       = MEMORY_NONE
     integer  :: step         = STEP_FIXED

     !----------------------------------------------------------------!
     ! What the cycle needs: which block each cell joins, and the two
     ! minimizers it governs. A minimizer holding minimizers is
     ! composition, not a second family.
     !----------------------------------------------------------------!

     integer, allocatable :: aggregates(:)
     integer              :: nblocks = 0

     type(fit_minimizer), allocatable :: smoother
     type(fit_minimizer), allocatable :: coarse

     !----------------------------------------------------------------!
     ! The last column: the minimizer that answers the linear
     ! question a re-read tangent poses.
     !----------------------------------------------------------------!

     type(fit_minimizer), allocatable :: inner

     character(len=:), allocatable :: label

   contains

     procedure :: attach
     procedure :: matvec
     procedure :: inner_product
     procedure :: norm
     procedure :: sweep_order
     procedure :: diagonal
     procedure :: constant
     procedure :: setup

     ! The operation face: a solver IS an operation - the one that
     ! answers the attached statement. apply solves from zero, so a
     ! solver composes wherever operations go; a preconditioner is
     ! exactly this face of an inner solver.
     procedure :: name   => solver_name
     procedure :: domain => solver_domain
     procedure :: apply  => solver_apply

     procedure :: solve

  end type fit_minimizer

  !===================================================================!
  ! THE STRUCTURE SIDE. It answers with a graph, because a form is
  ! one, so its lineage is the transform's while its role keeps it
  ! here: dependency sets the floor, role sets the home.
  !===================================================================!

  type, extends(graph_transform) :: form_minimizer

     integer  :: rule      = FORM_PRUNE
     real(dp) :: threshold = 1.0d-12

   contains

     procedure :: defined_on_graph => form_defined_on_graph
     procedure :: adapt

  end type form_minimizer

contains

  !===================================================================!
  ! THE CATALOGUE. Each name is a point, and the point is the whole
  ! difference between the methods.
  !===================================================================!

  pure type(fit_minimizer) function jacobi(omega) result(this)

    real(dp), intent(in), optional :: omega

    this % precondition = PRECONDITION_DIAGONAL
    this % memory       = MEMORY_NONE
    this % step         = STEP_FIXED
    if (present(omega)) this % relaxation = omega
    this % label        = 'jacobi'

  end function jacobi

  pure type(fit_minimizer) function gauss_seidel(omega) result(this)

    real(dp), intent(in), optional :: omega

    this % precondition = PRECONDITION_COLOUR
    this % memory       = MEMORY_NONE
    this % step         = STEP_FIXED
    if (present(omega)) this % relaxation = omega
    this % label        = 'gauss-seidel'

  end function gauss_seidel

  pure type(fit_minimizer) function conjugate_gradient() result(this)

    this % precondition = PRECONDITION_NONE
    this % memory       = MEMORY_FULL
    this % step         = STEP_GALERKIN
    this % label        = 'conjugate gradient'

  end function conjugate_gradient

  pure type(fit_minimizer) function gmres(restart) result(this)

    integer, intent(in), optional :: restart

    this % precondition = PRECONDITION_NONE
    this % memory       = 30
    this % step         = STEP_MINIMUM_RESIDUAL
    if (present(restart)) this % memory = max(restart, 1)
    this % label        = 'gmres'

  end function gmres

  type(fit_minimizer) function multigrid(smoother, coarse) result(this)

    type(fit_minimizer), intent(in) :: smoother
    type(fit_minimizer), intent(in) :: coarse

    this % precondition = PRECONDITION_CYCLE
    this % memory       = MEMORY_NONE
    this % step         = STEP_FIXED
    allocate(this % smoother, source=smoother)
    allocate(this % coarse  , source=coarse)
    this % label        = 'two-grid multigrid'

  end function multigrid

  !===================================================================!
  ! The last column. The statement is re-read at every state through
  ! a level-1 tangent, and the linear question it poses is handed to
  ! the minimizer this one governs. Nonlinear is not another family:
  ! it is one argument.
  !===================================================================!

  type(fit_minimizer) function newton(inner) result(this)

    type(fit_minimizer), intent(in) :: inner

    this % precondition = PRECONDITION_NONE
    this % memory       = MEMORY_NONE
    this % step         = STEP_FIXED
    allocate(this % inner, source=inner)
    this % label        = 'newton'

  end function newton

  pure type(form_minimizer) function pruner(threshold) result(this)

    real(dp), intent(in), optional :: threshold

    this % rule = FORM_PRUNE
    if (present(threshold)) this % threshold = threshold

  end function pruner

  !===================================================================!
  ! Take the operation and the graph it reads. The affine part is
  ! measured here, once: the operation applied to nothing is what
  ! its boundary values and sources say by themselves.
  !===================================================================!

  subroutine attach(this, action, on, ncomp)

    class(fit_minimizer)  , intent(inout) :: this
    class(graph_operation), intent(in)    :: action
    class(graph)          , intent(in)    :: on
    integer, intent(in), optional         :: ncomp

    real(dp), allocatable :: zero(:)
    integer :: v, nv

    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % on)) deallocate(this % on)
    allocate(this % on, source=on)

    this % ncomp = 1
    if (present(ncomp)) this % ncomp = max(ncomp, 1)

    nv = on % num_vertices()
    this % cells   = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    this % numbers = support(GRAPH_SIDE_VERTEX, &
         & [(v, v = 1, nv * this % ncomp)])

    allocate(zero(nv * this % ncomp))
    zero = 0.0_dp
    call raw_apply(this, zero, this % affine)

  end subroutine attach

  !===================================================================!
  ! The cycle's setup: read the fine stencil through the aggregates
  ! and hand the governed pair their statements. The smoother sweeps
  ! the fine one; the coarse minimizer answers the block one.
  !===================================================================!

  subroutine setup(this, aggregates)

    class(fit_minimizer), intent(inout) :: this
    integer, intent(in) :: aggregates(:)

    type(stencil_operator) :: block_statement
    type(support) :: blocks
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:), zeros(:)
    integer :: e, ne, b

    this % aggregates = aggregates
    this % nblocks    = maxval(aggregates)

    ! The Galerkin road: every dependency lands between blocks.
    select type (fine => this % action)

    type is (stencil_operator)

       ne = fine % pattern % num_edges()
       allocate(rows(ne), columns(ne), weights(ne))

       call fine % weights % get_real_vector(weights)
       do e = 1, ne
          rows(e)    = aggregates(fine % pattern % edge_head(e))
          columns(e) = aggregates(fine % pattern % edge_tail(e))
       end do

       allocate(zeros(this % nblocks))
       zeros = 0.0_dp

       block_statement = stencil_operator(rows, columns, weights, &
            & zeros, label='block statement')

    class default

       error stop 'the cycle: attach a compiled (stencil) operator'

    end select

    call this % smoother % attach(this % action, this % on)

    blocks = support(GRAPH_SIDE_VERTEX, [(b, b = 1, this % nblocks)])
    call this % coarse % attach(block_statement, blocks)

  end subroutine setup

  !===================================================================!
  ! The operation applied as it stands, affine part and all.
  !===================================================================!

  subroutine raw_apply(this, x, y)

    class(fit_minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    type(field) :: state
    class(graph_field), allocatable :: answer

    state = field('state', this % cells, ncomp=this % ncomp)
    call state % set_real_vector(x)

    call this % action % apply(this % on, [state], answer)
    call answer % get_real_vector(y)

  end subroutine raw_apply

  !===================================================================!
  ! The solver's words, each a delegation.
  !===================================================================!

  subroutine matvec(this, x, y)

    class(fit_minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    call raw_apply(this, x, y)
    y = y - this % affine

  end subroutine matvec

  real(dp) function inner_product(this, u, v) result(prod)

    class(fit_minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:), v(:)

    type(reduction) :: total
    type(field) :: uf, vf
    class(graph_functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = field('u', this % numbers)
    call uf % set_real_vector(u)
    vf = field('v', this % numbers)
    call vf % set_real_vector(v)

    total = reduction(REDUCE_SUM)
    call total % reduce(uf, answer, measure=vf)

    call answer % get_real_vector(got)
    prod = got(1)

  end function inner_product

  real(dp) function norm(this, u) result(length)

    class(fit_minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:)

    type(reduction) :: measure_of
    type(field) :: uf
    class(graph_functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = field('u', this % numbers)
    call uf % set_real_vector(u)

    measure_of = reduction(REDUCE_NORM)
    call measure_of % reduce(uf, answer)

    call answer % get_real_vector(got)
    length = got(1)

  end function norm

  subroutine sweep_order(this, colours)

    class(fit_minimizer), intent(in)  :: this
    integer, allocatable, intent(out) :: colours(:)

    type(walk) :: colouring
    class(graph_field), allocatable :: answer

    colouring = walk(WALK_COLOURING)
    call colouring % apply(this % on, output=answer)
    call answer % get_integer_vector(colours)

  end subroutine sweep_order

  !===================================================================!
  ! The diagonal, probed by colour. Applied to the indicator of one
  ! colour class, the answer at a member is that member's diagonal
  ! entry, because none of its neighbours is in the class. One
  ! matvec per colour, a handful in all.
  !===================================================================!

  subroutine diagonal(this, d)

    class(fit_minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: d(:)

    integer , allocatable :: colours(:)
    real(dp), allocatable :: indicator(:), y(:)
    integer :: nv, col, v

    if (this % ncomp > 1) then
       ! The probe reads one answer per cell, and a wide entry has
       ! several. A block probe is the honest generalization and no
       ! citizen has asked for one yet.
       error stop 'diagonal: the coloured probe answers one number per cell'
    end if

    nv = size(this % affine)
    allocate(d(nv), indicator(nv))
    d = 0.0_dp

    call this % sweep_order(colours)

    do col = 1, maxval(colours)

       indicator = 0.0_dp
       do v = 1, nv
          if (colours(v) == col) indicator(v) = 1.0_dp
       end do

       call this % matvec(indicator, y)

       do v = 1, nv
          if (colours(v) == col) d(v) = y(v)
       end do

    end do

  end subroutine diagonal

  !===================================================================!
  ! The affine part, for the caller assembling an equation: the
  ! statement action(q) = 0 reads matvec(q) = -constant.
  !===================================================================!

  subroutine constant(this, g)

    class(fit_minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: g(:)

    g = this % affine

  end subroutine constant

  !===================================================================!
  ! THE SOLVE. One dispatch on where this minimizer stands: the
  ! tangent column first, then the step rule, which is the only
  ! thing that changes the recurrence.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(fit_minimizer), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    if (allocated(this % inner)) then
       call solve_by_relinearizing(this, rhs, x, achieved)
       return
    end if

    select case (this % step)
    case (STEP_GALERKIN)
       call solve_by_conjugacy(this, rhs, x, achieved)
    case (STEP_MINIMUM_RESIDUAL)
       call solve_by_arnoldi(this, rhs, x, achieved)
    case default
       call solve_by_relaxing(this, rhs, x, achieved)
    end select

  end subroutine solve

  !===================================================================!
  ! MEMORY NONE, FIXED STEP: move a fixed distance along the
  ! preconditioned residual, remembering nothing. Which
  ! preconditioner is the whole difference between jacobi, a swept
  ! gauss-seidel, and a multigrid cycle.
  !===================================================================!

  subroutine solve_by_relaxing(this, rhs, x, achieved)

    class(fit_minimizer), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: d(:), y(:), r(:), rc(:), ec(:)
    integer , allocatable :: colours(:)
    real(dp) :: goal, elsewhere
    integer :: it, col, v

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    if (this % precondition == PRECONDITION_DIAGONAL .or. &
         & this % precondition == PRECONDITION_COLOUR) then
       call this % diagonal(d)
       ! A zero diagonal cannot correct its cell; leave that cell
       ! alone rather than divide by nothing.
       do v = 1, size(d)
          if (abs(d(v)) < tiny(1.0_dp)) d(v) = huge(1.0_dp)
       end do
    end if

    if (this % precondition == PRECONDITION_COLOUR) then
       call this % sweep_order(colours)
    end if

    if (this % precondition == PRECONDITION_CYCLE) then
       allocate(rc(this % nblocks), ec(this % nblocks))
    end if

    do it = 1, this % max_iterations

       select case (this % precondition)

       case (PRECONDITION_COLOUR)

          ! Each colour sees the corrections already made, because no
          ! two cells of one colour share a face.
          do col = 1, maxval(colours)
             call this % matvec(x, y)
             r = rhs - y
             do v = 1, size(x)
                if (colours(v) == col) then
                   x(v) = x(v) + this % relaxation * r(v) / d(v)
                end if
             end do
          end do

          call this % matvec(x, y)
          achieved = this % norm(rhs - y)
          if (achieved < goal) return

       case (PRECONDITION_CYCLE)

          ! Smooth, send the remainder down, bring the answer back,
          ! smooth again.
          call this % smoother % solve(rhs, x, elsewhere)

          call this % matvec(x, y)
          r = rhs - y

          achieved = this % norm(r)
          if (achieved < goal) return

          rc = 0.0_dp
          do v = 1, size(r)
             rc(this % aggregates(v)) = rc(this % aggregates(v)) + r(v)
          end do

          ec = 0.0_dp
          call this % coarse % solve(rc, ec, elsewhere)

          do v = 1, size(x)
             x(v) = x(v) + ec(this % aggregates(v))
          end do

          call this % smoother % solve(rhs, x, elsewhere)

       case default

          ! All cells at once, from the state as it stands.
          call this % matvec(x, y)
          r = rhs - y

          achieved = this % norm(r)
          if (achieved < goal) return

          if (this % precondition == PRECONDITION_DIAGONAL) then
             x = x + this % relaxation * r / d
          else
             x = x + this % relaxation * r
          end if

       end select

    end do

    call this % matvec(x, y)
    achieved = this % norm(rhs - y)

  end subroutine solve_by_relaxing

  !===================================================================!
  ! FULL MEMORY, GALERKIN STEP. When the constraint is the search
  ! space itself and the operator is symmetric, the whole history
  ! collapses into a three-term recurrence - the projection is still
  ! over every direction searched, but only one need be held:
  !
  !      alpha = (r, r) / (p, A p)         walk this far
  !      x <- x + alpha p                  along this direction
  !      r <- r - alpha A p                what remains
  !      p <- r + beta p                   the next direction, kept
  !                                        conjugate to all before it
  !===================================================================!

  subroutine solve_by_conjugacy(this, rhs, x, achieved)

    class(fit_minimizer), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: r(:), p(:), ap(:), y(:)
    real(dp) :: rr, rr_next, alpha, beta, pap, goal
    integer :: it

    call this % matvec(x, y)
    r = rhs - y
    p = r

    rr   = this % inner_product(r, r)
    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do it = 1, this % max_iterations

       achieved = this % norm(r)
       if (achieved < goal) return

       call this % matvec(p, ap)
       pap = this % inner_product(p, ap)
       if (abs(pap) < tiny(1.0_dp)) return

       alpha = rr / pap
       x = x + alpha * p
       r = r - alpha * ap

       rr_next = this % inner_product(r, r)
       beta    = rr_next / rr
       rr      = rr_next

       p = r + beta * p

    end do

    achieved = this % norm(r)

  end subroutine solve_by_conjugacy

  !===================================================================!
  ! BOUNDED MEMORY, MINIMUM RESIDUAL. The constraint is A V, so the
  ! coefficients are the ones that make the residual smallest over
  ! the whole space held; with no symmetry to lean on, the space is
  ! built explicitly and the small problem solved by rotations. The
  ! memory is the restart: how much is remembered before beginning
  ! again from where it stands.
  !===================================================================!

  subroutine solve_by_arnoldi(this, rhs, x, achieved)

    class(fit_minimizer), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: basis(:,:), h(:,:), cs(:), sn(:), s(:)
    real(dp), allocatable :: r(:), w(:), y(:)
    real(dp) :: goal, beta, hik, radius, subdiag
    integer :: n, m, outer, i, j, k

    n = size(x)
    m = max(this % memory, 1)

    allocate(basis(n, m + 1), h(m + 1, m), cs(m), sn(m), s(m + 1))

    goal = this % tolerance * (1.0_dp + this % norm(rhs))

    do outer = 1, this % max_iterations

       call this % matvec(x, r)
       r = rhs - r
       beta = this % norm(r)

       achieved = beta
       if (beta < goal) return

       basis = 0.0_dp
       h  = 0.0_dp
       s  = 0.0_dp
       s(1) = beta
       basis(:, 1) = r / beta

       do j = 1, m

          ! One more direction from the operator, kept orthonormal.
          ! The raw subdiagonal is kept aside: the rotations will
          ! overwrite its slot, and both the next basis vector and
          ! the breakdown test need the true number.
          call this % matvec(basis(:, j), w)
          do i = 1, j
             h(i, j) = this % inner_product(w, basis(:, i))
             w = w - h(i, j) * basis(:, i)
          end do
          subdiag     = this % norm(w)
          h(j + 1, j) = subdiag

          if (j < m .and. subdiag > tiny(1.0_dp)) then
             basis(:, j + 1) = w / subdiag
          end if

          ! The rotations that came before, then the new one.
          do i = 1, j - 1
             hik        = cs(i) * h(i, j) + sn(i) * h(i + 1, j)
             h(i + 1, j) = -sn(i) * h(i, j) + cs(i) * h(i + 1, j)
             h(i, j)     = hik
          end do

          radius = sqrt(h(j, j)**2 + h(j + 1, j)**2)
          if (radius < tiny(1.0_dp)) then
             k = j
             exit
          end if
          cs(j) = h(j, j) / radius
          sn(j) = h(j + 1, j) / radius
          h(j, j)     = radius
          h(j + 1, j) = 0.0_dp

          s(j + 1) = -sn(j) * s(j)
          s(j)     =  cs(j) * s(j)

          k = j
          achieved = abs(s(j + 1))
          if (achieved < goal) exit
          if (subdiag < tiny(1.0_dp)) exit

       end do

       ! The best answer this basis can express.
       allocate(y(k))
       do i = k, 1, -1
          y(i) = s(i)
          do j = i + 1, k
             y(i) = y(i) - h(i, j) * y(j)
          end do
          y(i) = y(i) / h(i, i)
       end do
       do i = 1, k
          x = x + y(i) * basis(:, i)
       end do
       deallocate(y)

       if (achieved < goal) then
          call this % matvec(x, r)
          achieved = this % norm(rhs - r)
          return
       end if

    end do

  end subroutine solve_by_arnoldi

  !===================================================================!
  ! THE TANGENT COLUMN. Drive the statement toward rhs from where we
  ! stand by re-reading it at every state: freeze the tangent, hand
  ! the linear question to the governed minimizer, step. This owns no
  ! derivative mathematics - the tangent is a level-1 citizen, and an
  ! exactly-linearizing statement takes the same seat with no change
  ! here.
  !===================================================================!

  subroutine solve_by_relinearizing(this, rhs, x, achieved)

    class(fit_minimizer), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    type(difference_linearization) :: jacobian
    real(dp), allocatable :: residual(:), g(:), y(:), dq(:)
    real(dp) :: elsewhere
    integer :: it

    allocate(dq(size(x)))

    call this % constant(g)

    jacobian = difference_linearization(this % action)

    do it = 1, this % max_iterations

       ! Where we stand: the full statement, whatever its linearity.
       call this % matvec(x, y)
       residual = y + g - rhs

       achieved = this % norm(residual)
       if (achieved < this % tolerance) return

       ! The linear question at this point, answered by the governed
       ! minimizer.
       call jacobian % freeze(x, base=y + g)

       call this % inner % attach(jacobian, this % on, ncomp = this % ncomp)
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, elsewhere)

       x = x + dq

    end do

    call this % matvec(x, y)
    achieved = this % norm(y + g - rhs)

  end subroutine solve_by_relinearizing

  !===================================================================!
  ! The operation face.
  !===================================================================!

  pure function solver_name(this) result(name)

    class(fit_minimizer), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = 'fit minimizer'
    end if

  end function solver_name

  subroutine solver_domain(this, input_graph, domain)

    class(fit_minimizer), intent(in)       :: this
    class(graph), intent(in)               :: input_graph
    class(graph), allocatable, intent(out) :: domain

    associate (u1 => this); end associate

    call input_graph % all_vertices(domain)

  end subroutine solver_domain

  subroutine solver_apply(this, input_graph, input_data, output)

    class(fit_minimizer), intent(in)               :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(fit_minimizer) :: worker
    type(field) :: out
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    allocate(x(size(this % affine)))
    x = 0.0_dp

    if (present(input_data)) then
       call input_data(1) % get_real_vector(rhs)
       worker = this
       call worker % solve(rhs, x, achieved)
    end if

    out = field('solution', this % cells, ncomp=this % ncomp)
    call out % set_real_vector(x)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine solver_apply

  !===================================================================!
  ! THE STRUCTURE SIDE. A form arrives, a smaller form leaves. The
  ! rule says which members go; the restriction itself is the form's
  ! own act, so nothing here reaches into anyone.
  !===================================================================!

  pure logical function form_defined_on_graph(this, input_graph)

    class(form_minimizer), intent(in) :: this
    class(graph)         , intent(in) :: input_graph

    associate (u1 => this); end associate

    form_defined_on_graph = input_graph % num_vertices() > 0

  end function form_defined_on_graph

  !===================================================================!
  ! Strike what the points cannot see. A basis member whose column of
  ! values vanishes over the whole constellation carries no
  ! information about it, and a fit that keeps it is solving for a
  ! coefficient nothing determines. The members are read about the
  ! constellation's own middle, so what is invisible does not depend
  ! on where the points happen to sit. The constant member always
  ! survives: every point sees it.
  !===================================================================!

  subroutine adapt(this, shape, positions, adapted)

    class(form_minimizer), intent(in)      :: this
    class(form)          , intent(in)      :: shape
    real(dp)             , intent(in)      :: positions(:)
    class(form), allocatable, intent(out)  :: adapted

    real(dp), allocatable :: phi(:), seen(:)
    real(dp) :: middle(3)
    integer :: nc, npts, j, m

    nc   = shape % size_of()
    npts = size(positions) / 3

    middle = 0.0_dp
    do j = 1, npts
       middle = middle + positions(3 * j - 2 : 3 * j)
    end do
    middle = middle / real(max(npts, 1), dp)

    allocate(phi(nc), seen(nc))
    seen = 0.0_dp

    do j = 1, npts
       call shape % values(positions(3 * j - 2 : 3 * j), middle, phi)
       do m = 1, nc
          seen(m) = seen(m) + phi(m) * phi(m)
       end do
    end do

    allocate(adapted, source=shape)
    call adapted % restrict(pack([(m, m = 1, nc)], seen > this % threshold))

  end subroutine adapt

end module graph_minimization
