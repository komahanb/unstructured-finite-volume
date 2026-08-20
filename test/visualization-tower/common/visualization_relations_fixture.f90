!=====================================================================!
! THE SPECIMEN'S PRIMITIVE INCIDENCE - earned at LEVEL 1.
!
! For each operator A_k, two typed relations and no more:
!
!      T_k  <=  E_k x X_(k-1)      where the occurrence reads FROM
!      H_k  <=  E_k x X_k          where the occurrence writes TO
!
! Twelve occurrences in all, each with exactly one tail and exactly
! one head, and each end living in its own correctly typed carrier.
!
!      A1        A2        A3
!      --------  --------  --------
!      e11 a p   e21 p u   e31 u m
!      e12 b p   e22 q u   e32 v n
!      e13 b q   e23 q v   e33 w n
!      e14 c q   e24 r w
!      e15 d r
!
!                       WHAT IS *NOT* WRITTEN HERE
!
! D1, D2 and D3 are absent, and their absence is the level's whole
! discipline. The dependency of an operator is NOT a primitive of
! this specimen: it is what follows algebraically from the
! occurrences, and Level 2 derives it. Writing D1 out by hand beside
! T1 and H1 would make the derivation a comparison of two hand-typed
! tables, which proves nothing about relation algebra.
!
! Coefficients are absent too, and will be at Gate A's end. T and H
! say WHICH members an occurrence touches. Nobody has yet said by how
! much, and nothing above needs to know.
!
!                    THIS IS INCIDENCE, NOT A JACOBIAN
!
! T and H are the tail and head of a dependency occurrence. They are
! not a matrix, not a sparsity pattern, and not the derivative of
! anything: a Jacobian would need an A to differentiate, and at this
! level no A has been given numbers to differentiate. The resemblance
! to a mesh's edge incidence is exact and deliberate - it is the same
! mathematics, read over a different specimen.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_relations_fixture

  use graph_fractal        , only : set_graph => graph
  use map_set        , only : set_map
  use relation_binary, only : csr_relation
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use visualization_assert , only : E1_1, E1_2, E1_3, E1_4, E1_5
  use visualization_assert , only : E2_1, E2_2, E2_3, E2_4
  use visualization_assert , only : E3_1, E3_2, E3_3

  implicit none

  private
  public :: occurrences_of_a1, occurrences_of_a2, occurrences_of_a3

contains

  !===================================================================!
  ! A1's five occurrences, split into their two ends.
  !
  !      e11 : a -> p        e14 : c -> q
  !      e12 : b -> p        e15 : d -> r
  !      e13 : b -> q
  !===================================================================!

  subroutine occurrences_of_a1(e1, x0, x1, tail, head, sets)

    type(set_graph) , intent(in)  :: e1, x0, x1
    type(csr_relation), intent(out) :: tail, head
    type(set_map)  , intent(in) :: sets

    tail = csr_relation('T1', e1, x0, reshape( &
         & [E1_1, X0_A, &
         &  E1_2, X0_B, &
         &  E1_3, X0_B, &
         &  E1_4, X0_C, &
         &  E1_5, X0_D], [2, 5]), sets)

    head = csr_relation('H1', e1, x1, reshape( &
         & [E1_1, X1_P, &
         &  E1_2, X1_P, &
         &  E1_3, X1_Q, &
         &  E1_4, X1_Q, &
         &  E1_5, X1_R], [2, 5]), sets)

  end subroutine occurrences_of_a1

  !===================================================================!
  ! A2's four occurrences.
  !
  !      e21 : p -> u        e23 : q -> v
  !      e22 : q -> u        e24 : r -> w
  !===================================================================!

  subroutine occurrences_of_a2(e2, x1, x2, tail, head, sets)

    type(set_graph) , intent(in)  :: e2, x1, x2
    type(csr_relation), intent(out) :: tail, head
    type(set_map)  , intent(in) :: sets

    tail = csr_relation('T2', e2, x1, reshape( &
         & [E2_1, X1_P, &
         &  E2_2, X1_Q, &
         &  E2_3, X1_Q, &
         &  E2_4, X1_R], [2, 4]), sets)

    head = csr_relation('H2', e2, x2, reshape( &
         & [E2_1, X2_U, &
         &  E2_2, X2_U, &
         &  E2_3, X2_V, &
         &  E2_4, X2_W], [2, 4]), sets)

  end subroutine occurrences_of_a2

  !===================================================================!
  ! A3's three occurrences.
  !
  !      e31 : u -> m        e33 : w -> n
  !      e32 : v -> n
  !===================================================================!

  subroutine occurrences_of_a3(e3, x2, x3, tail, head, sets)

    type(set_graph) , intent(in)  :: e3, x2, x3
    type(csr_relation), intent(out) :: tail, head
    type(set_map)  , intent(in) :: sets

    tail = csr_relation('T3', e3, x2, reshape( &
         & [E3_1, X2_U, &
         &  E3_2, X2_V, &
         &  E3_3, X2_W], [2, 3]), sets)

    head = csr_relation('H3', e3, x3, reshape( &
         & [E3_1, X3_M, &
         &  E3_2, X3_N, &
         &  E3_3, X3_N], [2, 3]), sets)

  end subroutine occurrences_of_a3

end module visualization_relations_fixture
