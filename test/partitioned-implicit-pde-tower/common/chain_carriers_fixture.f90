!=====================================================================!
! THE CHAIN CARRIERS FIXTURE - earned at LEVEL 0, and the lowest
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
! above depends on it being handled here first: all three carriers
! enumerate their members from one, so the integer 1 is a vertex, an
! edge AND a part. No count and no numeral separates these sets -
! only identity does.
!
! This file imports the set modules and nothing else. It is the only
! place in the tower where the three carriers are constructed, and
! every level from 0 upward takes them from here rather than
! declaring its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module chain_carriers_fixture

  use fractal_graph           , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map           , only : set_map

  implicit none

  private
  public :: chain_carriers

contains

  !===================================================================!
  ! The three carriers, declared and nothing more.
  !===================================================================!

  subroutine chain_carriers(sets, v, e, k)

    type(set_map)  , intent(inout) :: sets
    type(set_graph), intent(out)   :: v, e, k

    call v % declare()      ! global vertices
    call e % declare()      ! global edges
    call k % declare()      ! partition labels

    call sets % bind(v, counted_set_representation(6))
    call sets % bind(e, counted_set_representation(5))
    call sets % bind(k, counted_set_representation(2))

  end subroutine chain_carriers

end module chain_carriers_fixture
