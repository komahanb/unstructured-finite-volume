!=====================================================================!
! VISUALIZATION TOWER . LEVEL 7 . MINIMIZATION
!
! The level answers one question: WHAT STRUCTURE DOES MINIMIZATION
! ACTUALLY CONSUME - and is it the structure of the operator being
! minimized?
!
!                        LEVEL 7 IS INHABITED
!
! It is not N/A. Minimization already consumes graph structure:
!
!      minimizer % sweep_order()   colours a graph
!      minimizer % diagonal()      probes one indicator per colour
!      jacobi % solve()            consumes diagonal()
!      gauss_seidel % solve()      consumes diagonal() and sweep_order()
!
! So the question is not WHETHER structure is consumed but WHOSE.
!
!                          THE THREE STRUCTURES
!
!      P_A     the operator's own numerical coupling - what
!              stencil % dependencies() answers, and what
!              stencil % apply() actually computes with
!
!      H       the host graph handed through the legacy operation and
!              minimizer interface
!
!      C(H)    the colouring sweep_order() computes - and it is
!              computed from H, never from P_A
!
!                    THE RED, AND WHAT RESOLVED IT
!
! When this level was first written, sweep_order() coloured the
! ATTACHED HOST. The experiment below then reported
!
!      diagonal on H_match   [4, 5, 6]
!      diagonal on H_empty   [5, 7, 7]      <-- A times the all-ones
!
! because an edgeless host colours everything 1, so the coloured probe
! fired one all-ones indicator and read a matvec where it expected a
! diagonal. The diagonal of an unchanged matrix had moved because an
! irrelevant context graph changed.
!
! THE RESOLUTION WAS NOT TO ASK THE ACTION. `dependencies()` is
! axis-relative - a stencil answers on the dependent variable, a step
! on the independent one - so a minimizer reaching in for it would
! have coloured a BDF2 step's unknowns by a fan-in over instants.
! Instead the minimizer gained a SEAT:
!
!      on         the execution context, handed to apply
!      coupling   the dependent-variable stencil, coloured by
!                 sweep_order, supplied EXPLICITLY at attach and
!                 never defaulted from `on`
!
! The caller knows which object owns the dependent axis. Here it is
! the stencil, and the stencil says so through its own
! dependencies(). Below, both attachments are given coupling = P_A
! and differ only in `on` - so the diagonal must not move, and does
! not.
!
!                            THE EXPERIMENT
!
! ONE stencil, TWO hosts, both of the correct numerical dimension:
!
!      A = [ 4 1 0 ]        H_match  = P_A itself, the most
!          [ 1 5 1 ]                   favourable host possible
!          [ 0 1 6 ]        H_empty  = three vertices, no edges
!
! The control comes first and is load-bearing: the stencil does its
! arithmetic entirely from its OWN stored pattern, so
!
!      A_match x  =  A_empty x  =  [6, 14, 20]      for x = [1,2,3]
!
! and the host is irrelevant to the numerical map. Therefore any later
! difference caused by host topology alone belongs to the MINIMIZER'S
! STRUCTURAL INTERPRETATION and to nothing else.
!
!                     THE ACCEPTANCE CRITERION
!
!      diagonal(A)  =  [4, 5, 6]
!
! under BOTH hosts. diagonal() means the diagonal of the attached
! numerical action. Changing an irrelevant context graph must not
! change the diagonal of an unchanged matrix.
!
! The hostile host was never weakened to reach it. H_empty is still
! three vertices and no edges; what changed is that the solver is no
! longer asked to take its structure from it.
!
!                     NOTHING IS APPLIED DIRECTLY
!
! The tower never calls action % apply. Numerical probing at this
! level is earned only through the minimizer's own vocabulary -
! matvec, sweep_order, diagonal - which is what makes the measurement
! about the minimizer rather than about the stencil.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_7

  use map_set_representation, only : counted_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use iso_fortran_env      , only : dp => REAL64
  use visualization_assert , only : report, verdict
  use graph_fractal        , only : set_graph => graph
  use graph_directed_view  , only : directed_graph
  use class_graph          , only : directed_stored_graph
  use class_graph_stencil  , only : stencil_operator
  use class_graph_jacobi   , only : jacobi
  use production_pattern_renderer_fixture, only : pattern_picture
  use production_pattern_renderer_fixture, only : production_has
  use structural_renderer_fixture        , only : picture

  implicit none

  !-------------------------------------------------------------------!
  ! The specimen, and the oracles. Worked out on paper from
  !
  !      A = [ 4 1 0 ]
  !          [ 1 5 1 ]
  !          [ 0 1 6 ]
  !
  ! and never obtained from the machinery under test.
  !-------------------------------------------------------------------!

  integer , parameter :: N        = 3
  real(dp), parameter :: TOL      = 1.0e-12_dp

  real(dp), parameter :: X_PROBE(N)  = [1.0_dp, 2.0_dp, 3.0_dp]
  real(dp), parameter :: A_TIMES_X(N) = [6.0_dp, 14.0_dp, 20.0_dp]
  real(dp), parameter :: TRUE_DIAG(N) = [4.0_dp, 5.0_dp, 6.0_dp]

  ! What the coloured probe would report if it read ONE colour for
  ! everything: A times the all-ones vector. Present as a diagnostic
  ! oracle only - never as an acceptance criterion.
  real(dp), parameter :: A_TIMES_ONE(N) = [5.0_dp, 7.0_dp, 7.0_dp]

  type(stencil_operator)    :: a
  class(directed_graph), allocatable :: pa
  type(directed_stored_graph)        :: h_empty
  type(jacobi)              :: on_match, on_empty
  type(set_map)     :: sets
  type(label_map)     :: labels

  integer , allocatable :: col_match(:), col_empty(:)
  real(dp), allocatable :: mv_match(:), mv_empty(:)
  real(dp), allocatable :: d_match(:), d_empty(:)

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 7 . minimization"
  write(*,'(1x,a)') "============================================="

  ! ---- ONE production stencil, built through the actual API.
  !      edge k runs column -> row with weight k.
  a = stencil_operator(rows     = [1, 2, 3, 1, 2, 2, 3], &
       &               columns  = [1, 2, 3, 2, 1, 3, 2], &
       &               weights  = [4.0_dp, 5.0_dp, 6.0_dp, &
       &                           1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
       &               constant = [0.0_dp, 0.0_dp, 0.0_dp], &
       &               label    = 'A')

  ! ---- the operator's OWN dependency structure, from production.
  call a % dependencies(pa)
  if (.not. sets % describes(pa % vertex_set())) &
       & call sets % bind(pa % vertex_set(), &
       &      counted_set_representation(pa % num_vertices()))
  if (.not. sets % describes(pa % edge_set())) &
       & call sets % bind(pa % edge_set(), &
       &      counted_set_representation(pa % num_edges()))
  ! the production graph names its own two domains
  if (.not. labels % labelled(pa % vertex_set())) &
       & call labels % bind(pa % vertex_set(), 'vertices')
  if (.not. labels % labelled(pa % edge_set())) &
       & call labels % bind(pa % edge_set(), 'edges')

  ! ---- two hosts, both of the right numerical dimension.
  !      H_match IS P_A - the most favourable host there is.
  !      H_empty has the same three vertices and no edges at all.
  h_empty = directed_stored_graph(N, tails=[integer ::], heads=[integer ::])
  if (.not. sets % describes(h_empty % vertex_set())) &
       & call sets % bind(h_empty % vertex_set(), &
       &      counted_set_representation(h_empty % num_vertices()))
  if (.not. sets % describes(h_empty % edge_set())) &
       & call sets % bind(h_empty % edge_set(), &
       &      counted_set_representation(h_empty % num_edges()))
  if (.not. labels % labelled(h_empty % vertex_set())) &
       & call labels % bind(h_empty % vertex_set(), 'vertices')
  if (.not. labels % labelled(h_empty % edge_set())) &
       & call labels % bind(h_empty % edge_set(), 'edges')

  call attach_to(on_match, pa)
  call attach_to(on_empty, h_empty)

  call on_match % matvec(X_PROBE, mv_match)
  call on_empty % matvec(X_PROBE, mv_empty)

  call on_match % sweep_order(col_match)
  call on_empty % sweep_order(col_empty)

  call on_match % diagonal(d_match)
  call on_empty % diagonal(d_empty)

  call say_the_experiment()

  call check_the_operator_owns_its_structure(nfail)
  call check_the_hosts_differ(nfail)
  call check_the_numerical_map_is_host_independent(nfail)
  call check_the_colouring_is_host_dependent(nfail)
  call check_the_diagonal(nfail)

  call verdict(nfail, "level 7")

contains

  !===================================================================!
  ! Attach a jacobi to the one stencil over a given host, using that
  ! host's own vertex carrier as the unknown domain.
  !===================================================================!

  subroutine attach_to(solver, on)

    type(jacobi), intent(out) :: solver
    class(directed_graph), intent(in)  :: on

    type(set_graph) :: unknowns

    ! `on` is a production host whose carriers this scope must
    ! describe before it can size anything on them.
    unknowns = on % vertex_set()
    if (.not. sets % describes(unknowns)) &
         & call sets % bind(unknowns, &
         &      counted_set_representation(on % num_vertices()))

    ! Two seats, and the experiment varies only the first. `on` is
    ! where the action executes; `coupling` is which unknowns feed
    ! which, and the caller knows the stencil owns that axis.
    call solver % attach(a, on, unknowns, sets % size_of(unknowns), coupling = pa)

  end subroutine attach_to

  !===================================================================!
  ! THE STRUCTURAL PROVENANCE ERROR, MADE VISIBLE. A tower that hid
  ! this in one scalar assertion would have wasted seven levels.
  !===================================================================!

  subroutine say_the_experiment()

    type(picture) :: pic

    write(*,'(1x,a)') "---------------------------------------------"

    write(*,'(4x,a)') "OPERATOR A"
    write(*,'(8x,a)') "[ 4 1 0 ]"
    write(*,'(8x,a)') "[ 1 5 1 ]"
    write(*,'(8x,a)') "[ 0 1 6 ]"
    write(*,'(4x,a)') "true diagonal ....... 4  5  6"
    write(*,*)

    write(*,'(4x,a)') "OPERATOR DEPENDENCY P_A   (stencil % dependencies())"
    pic = pattern_picture(pa, '', sets, labels); call say_grid(pic)

    write(*,'(4x,a)') "HOST H_match  = P_A"
    call say_edges(pa)
    call say_ints ("solver colours .....", col_match)
    call say_reals("matvec([1,2,3]) ....", mv_match)
    call say_reals("diagonal ...........", d_match)
    write(*,*)

    write(*,'(4x,a)') "HOST H_empty  = 3 vertices, no edges"
    call say_edges(h_empty)
    call say_ints ("solver colours .....", col_empty)
    call say_reals("matvec([1,2,3]) ....", mv_empty)
    call say_reals("diagonal ...........", d_empty)
    write(*,*)

    write(*,'(4x,a)') "THE COMPARISON"
    call say_flag("same numerical operator under both hosts", &
         &        alike(mv_match, mv_empty))
    call say_flag("same matvec as the hand oracle [6,14,20]", &
         &        alike(mv_match, A_TIMES_X) .and. alike(mv_empty, A_TIMES_X))
    call say_flag("same solver colouring", same_ints(col_match, col_empty))
    call say_flag("diagonal host-dependent", &
         &        .not. alike(d_match, d_empty))
    call say_flag("diagonal correct on H_match", alike(d_match, TRUE_DIAG))
    call say_flag("diagonal correct on H_empty", alike(d_empty, TRUE_DIAG))
    call say_flag("the old RED [5,7,7] is gone", &
         &        .not. alike(d_empty, A_TIMES_ONE))

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_experiment

  subroutine say_grid(pic)

    type(picture), intent(in) :: pic

    integer :: k

    do k = 2, pic % rows()
       write(*,'(4x,a)') pic % at(k)
    end do
    write(*,*)

  end subroutine say_grid

  subroutine say_edges(g)

    class(directed_graph), intent(in) :: g

    integer :: e

    if (g % num_edges() .eq. 0) then
       write(*,'(8x,a)') "(no edges)"
       return
    end if
    do e = 1, g % num_edges()
       write(*,'(8x,i0,a,i0)') g % edge_tail(e), " -> ", g % edge_head(e)
    end do

  end subroutine say_edges

  subroutine say_ints(label, v)

    character(len=*), intent(in) :: label
    integer         , intent(in) :: v(:)

    integer :: k

    write(*,'(4x,a)', advance='no') label
    do k = 1, size(v)
       write(*,'(1x,i0)', advance='no') v(k)
    end do
    write(*,*)

  end subroutine say_ints

  subroutine say_reals(label, v)

    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: v(:)

    integer :: k

    write(*,'(4x,a)', advance='no') label
    do k = 1, size(v)
       write(*,'(1x,f0.1)', advance='no') v(k)
    end do
    write(*,*)

  end subroutine say_reals

  subroutine say_flag(label, yes)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: yes

    character(len=52) :: dots

    dots = repeat('.', max(1, 52 - len(label)))
    if (yes) then
       write(*,'(8x,a,1x,a,1x,a)') label, trim(dots), "YES"
    else
       write(*,'(8x,a,1x,a,1x,a)') label, trim(dots), "NO"
    end if

  end subroutine say_flag

  !===================================================================!
  ! P_A is the operator's own coupling, and it is the same object
  ! whichever host the operator is later attached to - because
  ! dependencies() reads the stencil's stored pattern and nothing
  ! else.
  !===================================================================!

  subroutine check_the_operator_owns_its_structure(nfail)

    integer, intent(inout) :: nfail

    class(directed_graph), allocatable :: again

    call report(pa % num_vertices() .eq. N .and. pa % num_edges() .eq. 7, &
         & "P_A has three vertices and seven arrows - one per stencil " // &
         & "coefficient, diagonal included", nfail)

    call report(production_has(pa, 1, 1) .and. production_has(pa, 2, 2) .and. &
         &      production_has(pa, 3, 3) .and. production_has(pa, 2, 1) .and. &
         &      production_has(pa, 1, 2) .and. production_has(pa, 3, 2) .and. &
         &      production_has(pa, 2, 3) .and. &
         &      .not. production_has(pa, 3, 1) .and. &
         &      .not. production_has(pa, 1, 3), &
         & "and it is A's coupling exactly: 1 and 3 are not coupled, " // &
         & "everything else is", nfail)

    call a % dependencies(again)
    if (.not. sets % describes(again % vertex_set())) &
         & call sets % bind(again % vertex_set(), &
         &      counted_set_representation(again % num_vertices()))
    if (.not. sets % describes(again % edge_set())) &
         & call sets % bind(again % edge_set(), &
         &      counted_set_representation(again % num_edges()))
    ! the production graph names its own two domains
    if (.not. labels % labelled(again % vertex_set())) &
         & call labels % bind(again % vertex_set(), 'vertices')
    if (.not. labels % labelled(again % edge_set())) &
         & call labels % bind(again % edge_set(), 'edges')
    call report(same_pattern(pa, again), &
         & "P_A DOES NOT DEPEND ON A HOST - dependencies() reads the " // &
         & "stencil's own stored pattern, and answers the same twice", &
         & nfail)

  end subroutine check_the_operator_owns_its_structure

  !===================================================================!
  ! The two hosts are genuinely different, and both are legal.
  !===================================================================!

  subroutine check_the_hosts_differ(nfail)

    integer, intent(inout) :: nfail

    call report(pa % num_vertices() .eq. N .and. &
         &      h_empty % num_vertices() .eq. N, &
         & "both hosts carry the correct numerical dimension - three " // &
         & "vertices each", nfail)

    call report(pa % num_edges() .eq. 7 .and. h_empty % num_edges() .eq. 0, &
         & "and their topology differs completely: H_match IS P_A, " // &
         & "H_empty has no edges at all", nfail)

    call report(.not. same_pattern(pa, h_empty), &
         & "H_match /= H_empty - the only thing the experiment " // &
         & "changes", nfail)

  end subroutine check_the_hosts_differ

  !===================================================================!
  ! THE LOAD-BEARING CONTROL. The stencil computes from its own
  ! pattern, so the host cannot touch the numerical map. Measured
  ! through the minimizer's matvec, never by calling apply.
  !===================================================================!

  subroutine check_the_numerical_map_is_host_independent(nfail)

    integer, intent(inout) :: nfail

    call report(alike(mv_match, A_TIMES_X), &
         & "on H_match, matvec([1,2,3]) = [6,14,20] - the hand oracle", &
         & nfail)

    call report(alike(mv_empty, A_TIMES_X), &
         & "on H_empty, matvec([1,2,3]) = [6,14,20] - THE SAME", nfail)

    call report(alike(mv_match, mv_empty), &
         & "SO THE HOST IS IRRELEVANT TO THIS STENCIL'S NUMERICAL " // &
         & "MAP: the stencil owns its own sparse pattern", nfail)

  end subroutine check_the_numerical_map_is_host_independent

  !===================================================================!
  ! And yet the minimizer's structural interpretation changes. This
  ! is what makes Level 7 inhabited: minimization consumes graph
  ! structure, and it consumes the HOST'S.
  !===================================================================!

  subroutine check_the_colouring_is_host_dependent(nfail)

    integer, intent(inout) :: nfail

    call report(size(col_match) .eq. N .and. size(col_empty) .eq. N, &
         & "sweep_order answers one colour per unknown on both hosts", &
         & nfail)

    call report(properly_coloured(col_match, pa), &
         & "on H_match no two coupled unknowns share a colour - the " // &
         & "property, not the literal", nfail)

    call report(properly_coloured(col_empty, pa), &
         & "and on H_empty TOO - though that host has no edges at " // &
         & "all, so it could never have produced this colouring", nfail)

    call report(same_ints(col_match, col_empty), &
         & "SAME OPERATOR, SAME COUPLING => SAME SWEEP STRUCTURE, " // &
         & "under two completely different execution contexts", nfail)

    call report(maxval(col_empty) .gt. 1, &
         & "THE EMPTY HOST DID NOT FLATTEN THE COLOURING. Colouring " // &
         & "H_empty itself would give every unknown one colour; " // &
         & "colouring the COUPLING does not", nfail)

  end subroutine check_the_colouring_is_host_dependent

  !===================================================================!
  ! THE ACCEPTANCE CRITERION.
  !
  !      diagonal() means the diagonal of the attached numerical
  !      action.
  !
  ! An unchanged matrix has an unchanged diagonal. The oracle is
  ! [4,5,6] under every host, and this level does not negotiate it.
  !===================================================================!

  subroutine check_the_diagonal(nfail)

    integer, intent(inout) :: nfail

    call report(alike(d_match, TRUE_DIAG), &
         & "diagonal(A) on H_match = [4,5,6]", nfail)

    call report(alike(d_empty, TRUE_DIAG), &
         & "diagonal(A) on H_empty = [4,5,6] - THE DIAGONAL OF AN " // &
         & "UNCHANGED MATRIX DOES NOT DEPEND ON AN IRRELEVANT " // &
         & "CONTEXT GRAPH", nfail)

    call report(alike(d_match, d_empty), &
         & "and the two agree with each other", nfail)

    ! The RED this level captured, named so it cannot come back
    ! silently: [5,7,7] is A times the all-ones vector, which is what
    ! a probe reports when one colour covers everything.
    call report(.not. alike(d_empty, A_TIMES_ONE), &
         & "and neither is [5,7,7] - EXECUTION CONTEXT DOES NOT " // &
         & "DETERMINE SOLVER COUPLING", nfail)

  end subroutine check_the_diagonal

  !===================================================================!
  ! Helpers.
  !===================================================================!

  logical function alike(got, want)

    real(dp), intent(in) :: got(:), want(:)

    alike = (size(got) .eq. size(want))
    if (alike) alike = all(abs(got - want) .lt. TOL)

  end function alike

  logical function same_ints(a1, a2)

    integer, intent(in) :: a1(:), a2(:)

    same_ints = (size(a1) .eq. size(a2))
    if (same_ints) same_ints = all(a1 .eq. a2)

  end function same_ints

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

  !-------------------------------------------------------------------!
  ! Is this colouring a valid colouring of THAT coupling: no two
  ! distinct members joined by an off-diagonal arrow share a colour.
  ! A self-arrow is not a coupling between two unknowns.
  !-------------------------------------------------------------------!

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

end program visualization_level_7
