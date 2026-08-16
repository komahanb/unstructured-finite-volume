!=====================================================================!
! THE CHAIN RELATIONS FIXTURE - earned at LEVEL 1, and deliberately
! independent of any graph object.
!
! It holds three primitive facts and nothing else:
!
!      Tail <= E x V     e_i -> i
!      Head <= E x V     e_i -> i+1
!      Own  <= K x V     part1 -> 1,2,3   part2 -> 4,5,6
!
! The carriers V, E and K are NOT declared here. They are Level 0's
! property, they live in chain_carriers_fixture, and every
! constructor below receives them as arguments. This file cannot
! name a set into existence; it can only state facts over sets
! somebody else has already declared - which is precisely what
! makes it a Level-1 file and not a Level-0 one.
!
! These are STRUCTURAL truths only. There is no overlap here, no
! borrowed member, no local numbering, and no ordinary graph: those
! are interpretations that Levels 3 and 4 will impose on top.
!
! Note the hazard this file is careful about: all three carriers
! enumerate their members from one, so their raw ids COLLIDE - the
! integer 1 is a vertex, an edge and a part all at once. Nothing
! here may be read positionally across carriers, and orientation is
! carried by the SIGNATURE, never by the numbers.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module chain_relations_fixture

  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_binary_relation, only : csr_relation

  implicit none

  private
  public :: tail_relation, head_relation, own_relation

contains

  !===================================================================!
  ! Tail <= E x V : which vertex each edge leaves.
  !===================================================================!

  type(csr_relation) function tail_relation(e, v, sets) result(tail)

    type(set_graph), intent(in) :: e, v
    type(set_map)  , intent(in) :: sets

    integer :: table(2, 5), i

    do i = 1, 5
       table(:, i) = [i, i]
    end do
    tail = csr_relation('tail', e, v, table, sets)

  end function tail_relation

  !===================================================================!
  ! Head <= E x V : which vertex each edge enters.
  !===================================================================!

  type(csr_relation) function head_relation(e, v, sets) result(head)

    type(set_graph), intent(in) :: e, v
    type(set_map)  , intent(in) :: sets

    integer :: table(2, 5), i

    do i = 1, 5
       table(:, i) = [i, i + 1]
    end do
    head = csr_relation('head', e, v, table, sets)

  end function head_relation

  !===================================================================!
  ! Own <= K x V : the INTENDED ownership, stated before any
  ! partitioner exists to realize it.
  !===================================================================!

  type(csr_relation) function own_relation(k, v, sets) result(own)

    type(set_graph), intent(in) :: k, v
    type(set_map)  , intent(in) :: sets

    integer :: table(2, 6)

    table(:, 1) = [1, 1]
    table(:, 2) = [1, 2]
    table(:, 3) = [1, 3]
    table(:, 4) = [2, 4]
    table(:, 5) = [2, 5]
    table(:, 6) = [2, 6]
    own = csr_relation('owns', k, v, table, sets)

  end function own_relation

end module chain_relations_fixture
