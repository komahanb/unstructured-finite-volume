!=====================================================================!
! THE CHAIN SETS FIXTURE - earned at LEVEL 0, and the lowest
! rung of the fixture ladder.
!
! Before the chain is a relation it is three sets, and nothing but
! three sets:
!
!      V = { 1 2 3 4 5 6 }        global vertices
!      E = { e1 e2 e3 e4 e5 }     global edges
!      K = { part1 part2 }        partition labels
!
! No edge knows a vertex here, no part owns anything, and the words
! borrowed and overlap have no meaning yet. Those arrive at Levels 1
! and 4 respectively, on top of what this file declares.
!
! The hazard is worth naming at the very bottom, because everything
! above depends on it being handled here first: all three sets
! enumerate their members from one, so the integer 1 is a vertex, an
! edge AND a part. No count and no numeral separates these sets -
! only identity does.
!
! This file imports graph_set and nothing else. It is the only
! place in the tower where the three sets are constructed, and
! every level from 0 upward takes them from here rather than
! declaring its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module chain_sets_fixture

  use graph_set, only : index_set

  implicit none

  private
  public :: chain_sets

contains

  !===================================================================!
  ! The three sets, declared and nothing more.
  !===================================================================!

  subroutine chain_sets(v, e, k)

    type(index_set), intent(out) :: v, e, k

    v = index_set('global vertices', 6)
    e = index_set('global edges'   , 5)
    k = index_set('partition labels', 2)

  end subroutine chain_sets

end module chain_sets_fixture
