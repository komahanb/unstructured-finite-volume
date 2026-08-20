!=====================================================================!
! VISUALIZATION TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question: CAN THE TWO AXES AND THE EXECUTION
! CONTEXT COEXIST IN ONE COMPOSITION WITHOUT BEING CONFLATED?
!
! One object at a time, all four alive together:
!
!      INDEPENDENT-AXIS STENCIL     step % dependencies()
!                                   BDF2's fan-in, 1->3 2->3 3->3
!
!      DEPENDENT-AXIS STENCIL       stencil % dependencies()
!                                   which unknown feeds which
!
!      EXECUTION CONTEXT            the graph handed to apply,
!                                   DELIBERATELY DIFFERENT from both
!
!      MINIMIZER                    handed the dependent-axis stencil
!                                   EXPLICITLY, and nothing else
!
! The constitution is that these are three distinct structures that
! compose, not one structure wearing three names.
!
!                       WHAT IS DELIBERATELY ABSENT
!
! No Kronecker product is formed. No product-space type is introduced.
! The two axes are shown to remain distinct and composable, and the
! tower stops there - building the product would be choosing an
! abstraction this tower has not earned.
!
! The step is never applied and never solved. It is constructed, asked
! for its own axis, and set beside the others.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_8

  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use iso_fortran_env      , only : dp => REAL64
  use visualization_assert , only : report, verdict
  use graph_fractal        , only : set_graph => graph
  use view_directed  , only : directed_graph
  use view_directed_stored          , only : directed_stored_graph
  use class_graph_stencil  , only : stencil_operator
  use class_graph_step     , only : step_operator, bdf, backward_euler
  use class_graph_jacobi   , only : jacobi
  use production_pattern_renderer_fixture, only : pattern_picture
  use production_pattern_renderer_fixture, only : production_has
  use structural_renderer_fixture        , only : picture

  implicit none

  integer , parameter :: N = 3
  real(dp), parameter :: TOL = 1.0e-12_dp
  real(dp), parameter :: H_PROBE = 0.5_dp

  real(dp), parameter :: X_PROBE(N)   = [1.0_dp, 2.0_dp, 3.0_dp]
  real(dp), parameter :: A_TIMES_X(N) = [6.0_dp, 14.0_dp, 20.0_dp]
  real(dp), parameter :: TRUE_DIAG(N) = [4.0_dp, 5.0_dp, 6.0_dp]

  type(stencil_operator)    :: a
  type(step_operator)       :: clock, euler
  class(directed_graph), allocatable :: dependent, independent, independent_be
  type(directed_stored_graph)        :: context
  type(jacobi)              :: solver
  type(set_map)     :: sets
  type(label_map)     :: labels

  integer , allocatable :: colours(:)
  real(dp), allocatable :: mv(:), d(:)

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  ! ---- the dependent axis: a 3x3 state coupling.
  a = stencil_operator(rows     = [1, 2, 3, 1, 2, 2, 3], &
       &               columns  = [1, 2, 3, 2, 1, 3, 2], &
       &               weights  = [4.0_dp, 5.0_dp, 6.0_dp, &
       &                           1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
       &               constant = [0.0_dp, 0.0_dp, 0.0_dp], &
       &               label    = 'A')
  call a % dependencies(dependent)
  if (.not. sets % describes(dependent % vertex_set())) &
       & call sets % bind(dependent % vertex_set(), &
       &      counted_set_representation(dependent % num_vertices()))
  if (.not. sets % describes(dependent % edge_set())) &
       & call sets % bind(dependent % edge_set(), &
       &      counted_set_representation(dependent % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(dependent % vertex_set())) &
       & call labels % bind(dependent % vertex_set(), 'vertices')
  if (.not. labels % labelled(dependent % edge_set())) &
       & call labels % bind(dependent % edge_set(), 'edges')

  ! ---- the independent axis: BDF2 around that same action.
  clock = bdf(2, a, H_PROBE)
  call clock % dependencies(independent)
  if (.not. sets % describes(independent % vertex_set())) &
       & call sets % bind(independent % vertex_set(), &
       &      counted_set_representation(independent % num_vertices()))
  if (.not. sets % describes(independent % edge_set())) &
       & call sets % bind(independent % edge_set(), &
       &      counted_set_representation(independent % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(independent % vertex_set())) &
       & call labels % bind(independent % vertex_set(), 'vertices')
  if (.not. labels % labelled(independent % edge_set())) &
       & call labels % bind(independent % edge_set(), 'edges')

  euler = backward_euler(a, H_PROBE)
  call euler % dependencies(independent_be)
  if (.not. sets % describes(independent_be % vertex_set())) &
       & call sets % bind(independent_be % vertex_set(), &
       &      counted_set_representation(independent_be % num_vertices()))
  if (.not. sets % describes(independent_be % edge_set())) &
       & call sets % bind(independent_be % edge_set(), &
       &      counted_set_representation(independent_be % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(independent_be % vertex_set())) &
       & call labels % bind(independent_be % vertex_set(), 'vertices')
  if (.not. labels % labelled(independent_be % edge_set())) &
       & call labels % bind(independent_be % edge_set(), 'edges')

  ! ---- the execution context, deliberately neither of them.
  context = directed_stored_graph(N, tails=[integer ::], heads=[integer ::])
  if (.not. sets % describes(context % vertex_set())) &
       & call sets % bind(context % vertex_set(), &
       &      counted_set_representation(context % num_vertices()))
  if (.not. sets % describes(context % edge_set())) &
       & call sets % bind(context % edge_set(), &
       &      counted_set_representation(context % num_edges()))
  if (.not. labels % labelled(context % vertex_set())) &
       & call labels % bind(context % vertex_set(), 'vertices')
  if (.not. labels % labelled(context % edge_set())) &
       & call labels % bind(context % edge_set(), 'edges')

  ! ---- the solver, handed the DEPENDENT axis and nothing else.
  call solver % attach(a, context, context % vertex_set(), &
       &               context % num_vertices(), coupling = dependent)

  call solver % matvec(X_PROBE, mv)
  call solver % sweep_order(colours)
  call solver % diagonal(d)

  call say_the_constitution()

  call check_the_two_axes_are_distinct(nfail)
  call check_the_independent_axis(nfail)
  call check_the_context_is_a_third_thing(nfail)
  call check_the_solver_took_the_dependent_axis(nfail)

  call verdict(nfail, "level 8")

contains

  subroutine say_the_constitution()

    type(picture) :: pic

    write(*,'(1x,a)') "---------------------------------------------"

    write(*,'(4x,a)') "INDEPENDENT-VARIABLE STENCIL   step % dependencies()"
    pic = pattern_picture(independent, '', sets, labels); call say_grid(pic)

    write(*,'(4x,a)') "DEPENDENT-VARIABLE STENCIL     stencil % dependencies()"
    pic = pattern_picture(dependent, '', sets, labels); call say_grid(pic)

    write(*,'(4x,a)') "EXECUTION CONTEXT              handed to apply"
    pic = pattern_picture(context, '', sets, labels); call say_grid(pic)

    write(*,'(4x,a)') "SOLVER"
    write(*,'(8x,a)') "coupling ............ dependent-variable stencil"
    write(*,'(8x,a)') "does NOT use ........ the BDF stencil"
    write(*,'(8x,a)') "does NOT infer ...... coupling from the context"
    call say_ints ("colours ............. ", colours)
    call say_reals("matvec([1,2,3]) ..... ", mv)
    call say_reals("diagonal ............ ", d)

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_constitution

  subroutine say_grid(pic)

    type(picture), intent(in) :: pic

    integer :: k

    do k = 2, pic % rows()
       write(*,'(4x,a)') pic % at(k)
    end do
    write(*,*)

  end subroutine say_grid

  subroutine say_ints(label, v)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: v(:)

    integer :: k

    write(*,'(8x,a)', advance='no') label
    do k = 1, size(v)
       write(*,'(1x,i0)', advance='no') v(k)
    end do
    write(*,*)

  end subroutine say_ints

  subroutine say_reals(label, v)

    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: v(:)

    integer :: k

    write(*,'(8x,a)', advance='no') label
    do k = 1, size(v)
       write(*,'(1x,f0.1)', advance='no') v(k)
    end do
    write(*,*)

  end subroutine say_reals

  !===================================================================!
  ! Two stencils, three vertices each, and not the same structure.
  !===================================================================!

  subroutine check_the_two_axes_are_distinct(nfail)

    integer, intent(inout) :: nfail

    call report(dependent % num_vertices() .eq. N .and. &
         &      independent % num_vertices() .eq. N, &
         & "both axes carry three members here - the coincidence " // &
         & "that makes the distinction worth checking", nfail)

    call report(.not. same_pattern(dependent, independent), &
         & "AND THEY ARE NOT THE SAME STRUCTURE: one says which " // &
         & "unknown feeds which, the other which instants the " // &
         & "residual reads", nfail)

    call report(production_has(dependent, 2, 1) .and. &
         &      .not. production_has(independent, 2, 1), &
         & "the dependent axis couples unknown 2 to unknown 1; the " // &
         & "independent axis has no such arrow", nfail)

    call report(production_has(independent, 1, 3) .and. &
         &      .not. production_has(dependent, 1, 3), &
         & "the independent axis reads instant 1 into instant 3; the " // &
         & "dependent axis leaves 1 and 3 uncoupled", nfail)

  end subroutine check_the_two_axes_are_distinct

  !===================================================================!
  ! The independent axis is a stencil, and it is the scheme's own.
  !===================================================================!

  subroutine check_the_independent_axis(nfail)

    integer, intent(inout) :: nfail

    call report(independent % num_edges() .eq. 3 .and. &
         &      production_has(independent, 1, 3) .and. &
         &      production_has(independent, 2, 3) .and. &
         &      production_has(independent, 3, 3), &
         & "BDF2 fans in on the instant being solved for: 1->3, " // &
         & "2->3, 3->3", nfail)

    call report(independent_be % num_vertices() .eq. 2 .and. &
         &      independent_be % num_edges() .eq. 2 .and. &
         &      production_has(independent_be, 1, 2) .and. &
         &      production_has(independent_be, 2, 2), &
         & "and backward euler answers the same shape one instant " // &
         & "narrower: 1->2, 2->2", nfail)

    call report(.not. production_has(independent, 1, 2) .and. &
         &      production_has(independent, 3, 3), &
         & "A STENCIL IS NOT A CHRONOLOGY - succession would hold " // &
         & "1->2 and would carry no self-arrow", nfail)

  end subroutine check_the_independent_axis

  !===================================================================!
  ! And the context is a third object, matching neither.
  !===================================================================!

  subroutine check_the_context_is_a_third_thing(nfail)

    integer, intent(inout) :: nfail

    call report(context % num_vertices() .eq. N .and. &
         &      context % num_edges() .eq. 0, &
         & "the execution context has the right dimension and no " // &
         & "structure at all", nfail)

    call report(.not. same_pattern(context, dependent) .and. &
         &      .not. same_pattern(context, independent), &
         & "THREE DISTINCT STRUCTURES ALIVE AT ONCE: context, " // &
         & "dependent stencil, independent stencil", nfail)

  end subroutine check_the_context_is_a_third_thing

  !===================================================================!
  ! THE CONSTITUTION. The solver took the dependent axis, was not
  ! confused by the independent one, and did not infer anything from
  ! the context it runs on.
  !===================================================================!

  subroutine check_the_solver_took_the_dependent_axis(nfail)

    integer, intent(inout) :: nfail

    call report(alike(mv, A_TIMES_X), &
         & "matvec([1,2,3]) = [6,14,20] - the action executed over " // &
         & "the context, as it should", nfail)

    call report(properly_coloured(colours, dependent), &
         & "the colouring is a valid colouring OF THE DEPENDENT " // &
         & "STENCIL", nfail)

    call report(maxval(colours) .gt. 1, &
         & "and it is not the flat colouring the empty context would " // &
         & "have produced", nfail)

    call report(.not. properly_coloured(colours, independent), &
         & "nor is it a colouring of the INDEPENDENT stencil - the " // &
         & "solver did not take the BDF fan-in for its coupling", nfail)

    call report(alike(d, TRUE_DIAG), &
         & "diagonal = [4,5,6] - THE AXES REMAIN DISTINCT AND THE " // &
         & "ANSWER IS THE OPERATOR'S OWN", nfail)

  end subroutine check_the_solver_took_the_dependent_axis

  !===================================================================!
  ! Helpers.
  !===================================================================!

  logical function alike(got, want)

    real(dp), intent(in) :: got(:), want(:)

    alike = (size(got) .eq. size(want))
    if (alike) alike = all(abs(got - want) .lt. TOL)

  end function alike

  logical function same_pattern(p, q)

    class(directed_graph), intent(in) :: p, q

    integer :: i, j

    same_pattern = (p % num_vertices() .eq. q % num_vertices())
    if (.not. same_pattern) return

    do j = 1, p % num_vertices()
       do i = 1, p % num_vertices()
          same_pattern = same_pattern .and. &
               & (production_has(p, j, i) .eqv. production_has(q, j, i))
       end do
    end do

  end function same_pattern

  logical function properly_coloured(colours, coupling)

    integer     , intent(in) :: colours(:)
    class(directed_graph), intent(in) :: coupling

    integer :: i, j

    properly_coloured = .true.
    do j = 1, coupling % num_vertices()
       do i = 1, coupling % num_vertices()
          if (i .eq. j) cycle
          if (production_has(coupling, j, i)) then
             properly_coloured = properly_coloured .and. &
                  & (colours(i) .ne. colours(j))
          end if
       end do
    end do

  end function properly_coloured

end program visualization_level_8
