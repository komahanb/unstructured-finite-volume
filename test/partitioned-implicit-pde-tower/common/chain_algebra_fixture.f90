!=====================================================================!
! THE CHAIN ALGEBRA FIXTURE - earned at LEVEL 2.
!
! Three derivations, and none is written down as data. All follow
! from the Level-1 primitives by transpose and composition alone,
! reading the repository's convention
!
!      compose_binary(P_AB, P_BC, sets) = P_BC o P_AB.
!
! ADJACENCY.  A vertex follows another when an edge leaves the first
! and enters the second:
!
!      A = Head o Tail^T : V -> V         v -> e -> head(e)
!
! written in code as compose_binary(Tail^T, Head, sets).
!
! CANDIDATE EDGE-OWNERSHIP POLICIES.  An edge has two vertices, and
! each of them has an owner. Composing through either one gives a
! perfectly good map from edges to parts:
!
!      TailOwner = Own^T o Tail : E -> K   e -> tail(e) -> owner
!      HeadOwner = Own^T o Head : E -> K   e -> head(e) -> owner
!
! written in code as compose_binary(Tail, Own^T, sets) and
! compose_binary(Head, Own^T, sets).
!
! Both are total. Both are single-valued. Both therefore satisfy the
! reconstruction law - one edge, one owner - and on this chain they
! agree everywhere except at the crossing edge e3 = 3->4, where
! TailOwner says part1 and HeadOwner says part2.
!
! That is the point of deriving BOTH. The algebra tells you the
! consequences of a policy once you have named one; it does not name
! one for you. Choosing the tail as the ownership anchor is a
! semantic decision that lives outside these two lines of
! composition, and Level 4 discovers - empirically, by reading
! edge_owner_part - which of the two production made.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module chain_algebra_fixture

  use relation_finitary        , only : relation
  use graph_fractal         , only : set_graph => graph
  use map_set         , only : set_map
  use relation_algebra, only : compose_binary
  use relation_binary , only : binary_relation, csr_relation, &
       &                             transposed_relation, transpose_of

  implicit none

  private
  public :: derive_adjacency, derive_tail_owner, derive_head_owner

contains

  !===================================================================!
  ! A = Head o Tail^T : V -> V, through the edges.
  !===================================================================!

  type(csr_relation) function derive_adjacency(tail, head, sets) result(adj)

    class(binary_relation), intent(in), target :: tail, head
    type(set_map)         , intent(in)         :: sets

    type(transposed_relation) :: tail_t

    tail_t = transpose_of(tail)          ! V -> E
    adj    = compose_binary(tail_t, head, sets) ! V -> E -> V

  end function derive_adjacency

  !===================================================================!
  ! TailOwner = Own^T o Tail : E -> K, anchored at the vertex the
  ! edge LEAVES. A candidate policy, not the policy.
  !===================================================================!

  type(csr_relation) function derive_tail_owner(tail, own, sets) result(eo)

    class(binary_relation), intent(in), target :: tail, own
    type(set_map)         , intent(in)         :: sets

    type(transposed_relation) :: own_t

    own_t = transpose_of(own)          ! V -> K
    eo    = compose_binary(tail, own_t, sets) ! E -> V -> K

  end function derive_tail_owner

  !===================================================================!
  ! HeadOwner = Own^T o Head : E -> K, anchored at the vertex the
  ! edge ENTERS. Equally derivable, equally total, equally
  ! single-valued - and it disagrees with TailOwner exactly where a
  ! cut falls.
  !===================================================================!

  type(csr_relation) function derive_head_owner(head, own, sets) result(eo)

    class(binary_relation), intent(in), target :: head, own
    type(set_map)         , intent(in)         :: sets

    type(transposed_relation) :: own_t

    own_t = transpose_of(own)          ! V -> K
    eo    = compose_binary(head, own_t, sets) ! E -> V -> K

  end function derive_head_owner

end module chain_algebra_fixture
