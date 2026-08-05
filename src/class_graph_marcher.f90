!=====================================================================!
! The marcher: time as a graph, walked forward.
!
! LEVEL 2 OF THE STRATIFICATION. Time is not a special dimension
! here - it is one more graph: instants are vertices, steps are
! edges, and a step size is a number riding an edge, exactly as a
! spacing rides a mesh face,
!
!      (t0) --h--> (t1) --h--> (t2) --h--> ... --h--> (tn)
!
! A march walks that chain once, forward. At every edge the attached
! statement is read where the state stands and the state moves
! against it,
!
!      q  <-  q - h * action(q)
!
! which is the explicit euler step; the statement returns MINUS the
! velocity, matching the house convention that a balance measures
! what a cell has left over. The map z -> z^2 + c is this walk with
! h = 1 on the law S = z - z^2 - c: an identity, not an
! approximation.
!
! THE ADJOINT IS THE REVERSE WALK. The same chain traversed
! head-to-tail carries sensitivities backward through the steps;
! that walk, and the implicit members that need a minimizer per
! step, are this citizen's family still to arrive.
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

  implicit none

  private
  public :: marcher

  type :: marcher

     real(dp) :: step = 1.0_dp

   contains

     procedure :: instants
     procedure :: march

  end type marcher

contains

  !===================================================================!
  ! The time graph itself: one vertex per instant, one edge per
  ! step. Callers who want to see time as structure - or walk it in
  ! reverse - hold this.
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
  ! Walk the chain forward: one statement read and one move per
  ! edge, the state carried in place.
  !===================================================================!

  subroutine march(this, action, on, q, nsteps)

    class(marcher), intent(inout)      :: this
    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: on
    real(dp), intent(inout)            :: q(:)
    integer, intent(in)                :: nsteps

    type(stored_graph) :: chain
    type(support) :: cells
    type(field)   :: state
    class(graph_field), allocatable :: answer
    real(dp), allocatable :: s(:)
    integer :: nv, ncomp, e, v

    nv    = on % num_vertices()
    ncomp = size(q) / max(nv, 1)

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])

    call this % instants(nsteps, chain)

    do e = 1, chain % num_edges()

       state = field('state', cells, ncomp=ncomp)
       call state % set_real_vector(q)

       call action % apply(on, [state], answer)
       call answer % get_real_vector(s)

       q = q - this % step * s

    end do

  end subroutine march

end module class_graph_marcher
