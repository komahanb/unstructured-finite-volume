!=====================================================================!
! VISUALIZATION TOWER . LEVEL 9 . STATEMENT
!
! The level answers one question: WHAT COMPLETE STATEMENT DOES THIS
! TOWER MAKE?
!
! One output, five things in it, and none of them the same thing:
!
!      INDEPENDENT-VARIABLE STENCIL   which instants the residual
!                                     reads
!      DEPENDENT-VARIABLE STENCIL     which unknown feeds which
!      DEPENDENT-VARIABLE VALUES      the coefficients on that
!                                     structure, zero included
!      EXECUTION CONTEXT              where the action runs, and it
!                                     may differ from both
!      SOLVER STRUCTURE               what the minimizer coloured,
!                                     and the diagonal it read
!
! and one executable claim:
!
!      THE CORRECT STRUCTURE IS DETERMINED BY THE MATHEMATICAL AXIS,
!      NOT BY WHICHEVER GRAPH HAPPENS TO BE NEARBY.
!
! Every representation below comes from the tower's own test-local
! renderers. No production visualization API exists, no visualize()
! was added to any root, and no representation type was introduced.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_9

  use iso_fortran_env      , only : dp => REAL64
  use visualization_assert , only : report, verdict
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use relation_finitary       , only : relation
  use relation_binary, only : csr_relation
  use view_directed  , only : directed_graph
  use view_directed_stored          , only : directed_stored_graph
  use field_stored    , only : field
  use operation_stencil  , only : stencil_operator
  use operation_step     , only : step_operator, bdf
  use operation_jacobi   , only : jacobi
  use visualization_carriers_fixture , only : structural_carriers
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_algebra_fixture  , only : derive_dependency
  use visualization_values_fixture   , only : coefficients_of_a2
  use structural_renderer_fixture    , only : picture, sparsity_picture
  use valued_renderer_fixture        , only : valued_sparsity_picture
  use production_pattern_renderer_fixture, only : pattern_picture
  use production_pattern_renderer_fixture, only : production_has
  use production_pattern_renderer_fixture, only : signature_of_relation

  implicit none

  integer , parameter :: N = 3
  real(dp), parameter :: TOL = 1.0e-12_dp

  real(dp), parameter :: X_PROBE(N)   = [1.0_dp, 2.0_dp, 3.0_dp]
  real(dp), parameter :: A_TIMES_X(N) = [6.0_dp, 14.0_dp, 20.0_dp]
  real(dp), parameter :: TRUE_DIAG(N) = [4.0_dp, 5.0_dp, 6.0_dp]
  real(dp), parameter :: A_TIMES_ONE(N) = [5.0_dp, 7.0_dp, 7.0_dp]

  ! ---- the relational half, as Gate A left it
  type(set_graph)  :: x0, x1, x2, x3, e1, e2, e3
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(csr_relation) :: t2, h2, d2
  type(field)        :: w2

  ! ---- the production half
  type(stencil_operator)    :: a
  type(step_operator)       :: clock
  class(directed_graph), allocatable :: dependent, independent
  type(directed_stored_graph)        :: context
  type(jacobi)              :: solver

  integer , allocatable :: colours(:)
  real(dp), allocatable :: mv(:), d(:)

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 9 . statement"
  write(*,'(1x,a)') "============================================="

  call structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)
  call occurrences_of_a2(e2, x1, x2, t2, h2, sets)
  d2 = derive_dependency('D2', t2, h2, sets)
  w2 = coefficients_of_a2(e2, sets)

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

  clock = bdf(2, a, 0.5_dp)
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

  call solver % attach(a, context, context % vertex_set(), &
       &               context % num_vertices(), coupling = dependent)
  call solver % matvec(X_PROBE, mv)
  call solver % sweep_order(colours)
  call solver % diagonal(d)

  call say_the_statement()

  call check_the_statement(nfail)

  call verdict(nfail, "level 9")

contains

  subroutine say_the_statement()

    type(picture) :: pic

    write(*,'(1x,a)') "---------------------------------------------"

    write(*,'(4x,a)') "INDEPENDENT-VARIABLE STENCIL"
    pic = pattern_picture(independent, '', sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "DEPENDENT-VARIABLE STENCIL"
    pic = pattern_picture(dependent, '', sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "DEPENDENT-VARIABLE VALUES"
    write(*,'(4x,a)') "signature: " // signature_of_relation(d2, labels)
    pic = valued_sparsity_picture(d2, t2, h2, e2, w2, sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "EXECUTION CONTEXT"
    pic = pattern_picture(context, '', sets, labels); call say_grid(pic, 2)

    write(*,'(4x,a)') "SOLVER STRUCTURE"
    write(*,'(8x,a)') "coupling ............ dependent-variable stencil"
    call say_ints ("colours ............. ", colours)
    write(*,*)

    write(*,'(4x,a)') "CONTEXT-INVARIANT DIAGONAL"
    call say_reals("diagonal ............ ", d)

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_statement

  subroutine say_grid(pic, from)

    type(picture), intent(in) :: pic
    integer      , intent(in) :: from

    integer :: k

    do k = from, pic % rows()
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
  ! THE STATEMENT.
  !===================================================================!

  subroutine check_the_statement(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: from, to
    type(picture)                  :: valued

    ! ---- structure, from relations, with no numbers in it
    from = d2 % domain(1)
    to   = d2 % domain(2)
    call report(d2 % num_tuples() .eq. 4 .and. .not. from % same_as(to), &
         & "RELATION -> STRUCTURE: D2 : X1 -> X2, derived from " // &
         & "occurrences, two declared carriers", nfail)

    ! ---- values, inhabiting that structure without defining it
    valued = valued_sparsity_picture(d2, t2, h2, e2, w2, sets, labels)
    call report(valued % at(4) .eq. 'v          .   -2    .' .and. &
         &      d2 % has([2, 2]), &
         & "FIELD -> VALUES: the coefficients sit on the occurrences, " // &
         & "and the structure is untouched by them", nfail)

    ! ---- dependencies(), read on each concrete type's own axis
    call report(production_has(dependent, 2, 1) .and. &
         &      .not. production_has(dependent, 1, 3), &
         & "DEPENDENCIES -> THE STENCIL ON THE DEPENDENT AXIS: which " // &
         & "unknown feeds which", nfail)

    call report(production_has(independent, 1, 3) .and. &
         &      production_has(independent, 3, 3) .and. &
         &      .not. production_has(independent, 1, 2), &
         & "DEPENDENCIES -> THE STENCIL ON THE INDEPENDENT AXIS: " // &
         & "which instants the residual reads - a fan-in, not a " // &
         & "chronology", nfail)

    ! ---- the three structures, all alive, none of them each other
    call report(.not. same_pattern(dependent, independent) .and. &
         &      .not. same_pattern(dependent, context) .and. &
         &      .not. same_pattern(independent, context), &
         & "CONTEXT /= INDEPENDENT STENCIL /= DEPENDENT STENCIL - " // &
         & "three vertices each, three different structures", nfail)

    ! ---- the minimizer, given its coupling explicitly
    call report(alike(mv, A_TIMES_X) .and. properly_coloured(colours, dependent), &
         & "MINIMIZER -> EXPLICIT DEPENDENT-VARIABLE COUPLING: it ran " // &
         & "over the context and coloured the stencil", nfail)

    ! ---- THE STATEMENT
    call report(alike(d, TRUE_DIAG) .and. .not. alike(d, A_TIMES_ONE), &
         & "THE CORRECT STRUCTURE IS DETERMINED BY THE MATHEMATICAL " // &
         & "AXIS, NOT BY WHICHEVER GRAPH HAPPENS TO BE NEARBY: " // &
         & "diagonal = [4,5,6] over a context that knows nothing", &
         & nfail)

  end subroutine check_the_statement

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

end program visualization_level_9
