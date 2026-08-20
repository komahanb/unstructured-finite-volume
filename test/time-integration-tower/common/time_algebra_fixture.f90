!=====================================================================!
! THE TIME ALGEBRA FIXTURE - earned at LEVEL 2.
!
! Two derivations, and neither is written down as data. Both follow
! from the Level-1 primitives by transpose and composition alone,
! reading the repository's convention
!
!      compose_binary(P_AB, P_BC) = P_BC o P_AB.
!
! ONE-STEP REACH. An instant reaches another when some step leaves
! the first and enters the second:
!
!      A1 = Head o Tail^T : T -> T        t -> e -> head(e)
!
! written in code as compose_binary(Tail^T, Head). The chain
! t0 -> t1 -> t2 -> t3 -> t4 comes back without anyone writing it.
!
! TWO-STEP REACH. Compose the one-step reach with itself:
!
!      A2 = A1 o A1 : T -> T
!
! written as compose_binary(A1, A1), giving t0 -> t2, t1 -> t3,
! t2 -> t4.
!
!                    WHAT A2 IS NOT
!
! A2 is NOT BDF2. It is not a scheme, not a dependency stencil, and
! not a discretization of anything. It says exactly one thing:
!
!      an instant two steps later is structurally REACHABLE.
!
! A scheme such as BDF2 may LATER consume present, one-step history
! and two-step history - but that consumption is an interpretation
! laid on this structure by a level that has not been built. Naming
! these functions for reach rather than for schemes is not fussiness;
! it is where the conflation would first get in.
!
!      temporal reach  !=  temporal discretization scheme
!
! No union and no transitive closure appear here. Neither has a
! caller that earns it, and writing A1 U A2 merely because the
! notation is available would be inventing algebra for a reader
! rather than for a client.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_algebra_fixture

  use graph_fractal         , only : set_graph => graph
  use map_set         , only : set_map
  use relation_finitary        , only : relation
  use relation_algebra, only : compose_binary
  use relation_binary , only : binary_relation, csr_relation, &
       &                             transposed_view, transpose_of

  implicit none

  private
  public :: derive_one_step_reach, derive_two_step_reach

contains

  !===================================================================!
  ! A1 = Head o Tail^T : T -> T, through the steps.
  !===================================================================!

  type(csr_relation) function derive_one_step_reach(tail, head, sets) &
       & result(a1)

    class(binary_relation), intent(in), target :: tail, head
    type(set_map)         , intent(in)         :: sets

    type(transposed_view) :: tail_t

    tail_t = transpose_of(tail)                 ! T -> E
    a1     = compose_binary(tail_t, head, sets) ! T -> E -> T

  end function derive_one_step_reach

  !===================================================================!
  ! A2 = A1 o A1 : T -> T, the two-step REACH - and no more than
  ! that.
  !===================================================================!

  type(csr_relation) function derive_two_step_reach(a1, sets) result(a2)

    class(binary_relation), intent(in), target :: a1
    type(set_map)         , intent(in)         :: sets

    a2 = compose_binary(a1, a1, sets)    ! T -> T -> T

  end function derive_two_step_reach

end module time_algebra_fixture
