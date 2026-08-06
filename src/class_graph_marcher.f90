!=====================================================================!
! The marcher: time as a graph, walked both ways.
!
! LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension
! here - it is one more graph: instants are vertices, steps are
! edges, and a step size is a number riding an edge, exactly as a
! spacing rides a mesh face,
!
!      (t0) --h--> (t1) --h--> (t2) --h--> ... --h--> (tn)
!
! THE FORWARD WALK. At every edge the attached statement is read
! and the state moves against it. Three rules, absorbed:
!
!      MARCH_FORWARD    q <- q - h * action(q)          explicit
!      MARCH_BACKWARD   q - qold + h * action(q) = 0    implicit,
!                                                       one governed
!                                                       solve per edge
!      MARCH_BDF2       (3q - 4qold + qolder)/2
!                              + h * action(q) = 0      second order,
!                                                       started by one
!                                                       backward step
!
! The implicit rules GOVERN: the marcher holds a minimizer and hands
! it one step operator per edge - the level-1 citizen from the time
! shelf, time's own discretization stencil - the same family
! discipline as newton over its linear question. The statement
! returns MINUS the velocity, matching the house convention that a
! balance measures what a cell has left over; z -> z^2 + c at h = 1
! is the forward walk on S = z - z^2 - c, an identity.
!
! THE REVERSE WALK IS THE ADJOINT - OF THE EXPLICIT RULE. The same
! chain traversed head to tail carries sensitivities backward:
! handed the TRANSPOSED statement, lambda <- lambda - h *
! transposed(lambda) per edge, in reverse edge order. The pairing
! <lambda, q> is then invariant across every step - the duality the
! suite holds to machine precision. For statements that change along
! the walk the reversed order is not a courtesy; it is the
! derivative's chain rule. Under an implicit rule this walk refuses
! rather than lies: the adjoint of a solved step is the transpose of
! that solve, and nothing here holds the trajectory it would need.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_marcher

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph, graph_field, graph_operation
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use class_graph        , only : stored_graph
  use class_graph_step   , only : step_operator, bdf
  use graph_minimization , only : fit_minimizer

  implicit none

  private
  public :: marcher
  public :: MARCH_FORWARD, MARCH_BACKWARD, MARCH_BDF2

  integer, parameter :: MARCH_FORWARD  = 1
  integer, parameter :: MARCH_BACKWARD = 2
  integer, parameter :: MARCH_BDF2     = 3

  type :: marcher

     integer  :: rule = MARCH_FORWARD
     real(dp) :: step = 1.0_dp

     type(fit_minimizer), allocatable :: inner

   contains

     procedure :: instants
     procedure :: march
     procedure :: march_adjoint

  end type marcher

contains

  !===================================================================!
  ! The time graph itself: one vertex per instant, one edge per
  ! step. Callers who want to see time as structure hold this.
  !===================================================================!

  subroutine instants(this, nsteps, chain)

    class(marcher), intent(in) :: this
    integer, intent(in) :: nsteps
    type(stored_graph), intent(out) :: chain

    integer :: n

    associate (u1 => this); end associate

    chain = stored_graph(nsteps + 1, &
         & tails=[(n, n = 1, nsteps)], heads=[(n + 1, n = 1, nsteps)])

  end subroutine instants

  !===================================================================!
  ! Walk the chain forward, by the rule.
  !===================================================================!

  subroutine march(this, action, on, q, nsteps)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: on
    real(dp), intent(inout)            :: q(:)
    integer, intent(in)                :: nsteps

    type(stored_graph) :: chain
    type(step_operator) :: statement
    real(dp), allocatable :: s(:), qold(:), qolder(:), zeros(:)
    real(dp) :: answered
    integer :: e

    call this % instants(nsteps, chain)

    if (this % rule == MARCH_FORWARD) then

       do e = 1, chain % num_edges()
          call read_statement(action, on, q, s)
          q = q - this % step * s
       end do
       return

    end if

    ! The implicit rules: one governed solve per edge; bdf2 starts
    ! with a single backward step, as it must. The step comes off
    ! the time shelf; the per-edge state and the bdf2 start are
    ! written onto it as the walk proceeds.
    statement = bdf(1, action, this % step)
    allocate(zeros(size(q)))
    zeros = 0.0_dp

    qold = q

    do e = 1, chain % num_edges()

       if (this % rule == MARCH_BDF2 .and. e > 1) then
          statement % a0 = 1.5_dp
          statement % a1 = -2.0_dp
          statement % a2 = 0.5_dp
          statement % reach = 2
          statement % qolder = qolder
       else
          statement % a0 = 1.0_dp
          statement % a1 = -1.0_dp
          statement % a2 = 0.0_dp
          statement % reach = 1
       end if

       statement % hs   = this % step
       statement % qold = qold

       ! The state's width travels with it: a cell carrying several
       ! numbers is measured whole, not by its first stripe.
       call this % inner % attach(statement, on, &
            & ncomp = size(q) / max(on % num_vertices(), 1))
       call this % inner % solve(zeros, q, answered)

       qolder = qold
       qold   = q

    end do

  end subroutine march

  !===================================================================!
  ! Walk the chain in reverse, carrying the pairing: handed the
  ! transposed statement, the sensitivities travel head to tail and
  ! <lambda, q> stays put.
  !
  ! ONLY THE EXPLICIT RULE, AND IT SAYS SO. The adjoint of a walk is
  ! the adjoint of the walk that was actually taken: if the forward
  ! rule solves an implicit step at every edge, the reverse walk owes
  ! the transpose of that solve, and doing plain euler backwards
  ! instead returns a lambda that pairs with nothing. Worse, an
  ! implicit reverse walk needs the forward trajectory to linearize
  ! against, and this signature is handed no trajectory - only a
  ! statement and a lambda. So the refusal below is the honest state
  ! of the art, not a stub: the trajectory arrives when this citizen
  ! becomes backward substitution on the chain, and the reverse walk
  ! becomes the same verb as the forward one.
  !===================================================================!

  subroutine march_adjoint(this, transposed, on, lambda, nsteps)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: transposed
    class(graph), intent(in)           :: on
    real(dp), intent(inout)            :: lambda(:)
    integer, intent(in)                :: nsteps

    type(stored_graph) :: chain
    real(dp), allocatable :: s(:)
    integer :: e

    if (this % rule /= MARCH_FORWARD) then
       error stop 'march_adjoint: the reverse walk answers the explicit rule only'
    end if

    call this % instants(nsteps, chain)

    do e = chain % num_edges(), 1, -1
       call read_statement(transposed, on, lambda, s)
       lambda = lambda - this % step * s
    end do

  end subroutine march_adjoint

  !===================================================================!
  ! One read of a statement at the standing state.
  !===================================================================!

  subroutine read_statement(action, on, q, s)

    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: on
    real(dp), intent(in)               :: q(:)
    real(dp), allocatable, intent(out) :: s(:)

    type(support) :: cells
    type(field)   :: state
    class(graph_field), allocatable :: answer
    integer :: nv, ncomp, v

    nv    = on % num_vertices()
    ncomp = size(q) / max(nv, 1)

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    state = field('state', cells, ncomp=ncomp)
    call state % set_real_vector(q)

    call action % apply(on, [state], answer)
    call answer % get_real_vector(s)

  end subroutine read_statement

end module class_graph_marcher
