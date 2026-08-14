!=====================================================================!
! THE CHAIN ALGEBRA FIXTURE - earned at LEVEL 2.
!
! Two derivations, and neither is written down as data. Both follow
! from the Level-1 primitives by transpose and composition alone,
! reading the repository's convention
!
!      compose_binary(R_AB, R_BC) = R_BC o R_AB.
!
! ADJACENCY.  A vertex follows another when an edge leaves the first
! and enters the second:
!
!      A = Head o Tail^T : V -> V         v -> e -> head(e)
!
! written in code as compose_binary(Tail^T, Head).
!
! EDGE OWNERSHIP.  An edge belongs to whichever part owns the vertex
! it LEAVES:
!
!      EdgeOwner = Own^T o Tail : E -> K   e -> tail(e) -> owner
!
! written in code as compose_binary(Tail, Own^T).
!
! That second derivation matters out of proportion to its size. The
! tail-ownership rule was discovered OPERATIONALLY by the earlier
! gate work - by imposing an assembly law on a probe field and then
! reading what the partitioner had done. Here it is stated as
! relational mathematics BEFORE any partitioner exists, so that
! Level 4 can check production against a law rather than against a
! previous observation of production.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module chain_algebra_fixture

  use graph_relation        , only : relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : binary_relation, csr_relation, &
       &                             transposed_view, transpose_of

  implicit none

  private
  public :: derive_adjacency, derive_edge_owner

contains

  !===================================================================!
  ! A = Head o Tail^T : V -> V, through the edges.
  !===================================================================!

  type(csr_relation) function derive_adjacency(tail, head) result(adj)

    class(binary_relation), intent(in), target :: tail, head

    type(transposed_view) :: tail_t

    tail_t = transpose_of(tail)          ! V -> E
    adj    = compose_binary(tail_t, head) ! V -> E -> V

  end function derive_adjacency

  !===================================================================!
  ! EdgeOwner = Own^T o Tail : E -> K, through the tail vertex.
  !===================================================================!

  type(csr_relation) function derive_edge_owner(tail, own) result(eo)

    class(binary_relation), intent(in), target :: tail, own

    type(transposed_view) :: own_t

    own_t = transpose_of(own)          ! V -> K
    eo    = compose_binary(tail, own_t) ! E -> V -> K

  end function derive_edge_owner

end module chain_algebra_fixture
