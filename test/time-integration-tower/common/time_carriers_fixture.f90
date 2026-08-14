!=====================================================================!
! THE TIME CARRIERS FIXTURE - earned at LEVEL 0, and the lowest rung
! of the fixture ladder.
!
! Before time has a direction it is sets, and there are three of
! them - because time integration has TWO conceptually independent
! axes and one connective tissue between them:
!
!      Q = { x  y }                 state coordinates
!      T = { t0 t1 t2 t3 t4 }       time instants
!      E = { e1 e2 e3 e4 }          time steps
!
! Q IS PRESENT FROM THE BOTTOM, deliberately, and it is present
! before any value exists to live on it. The temptation this tower
! exists to refuse is the quiet collapse
!
!      state coordinate  ~  time instant  ~  graph vertex
!
! which no line of this file permits: Q, T and E are three declared
! identities, and nothing here relates any of them to any other.
!
! The hazard is the usual one, met by a second independent
! specimen: all three carriers enumerate from one, so the integer 1
! is a member of Q, of T and of E - meaning x, t0 and e1
! respectively. No count and no numeral separates them; only
! identity does.
!
! This file imports graph_carrier and nothing else. It is the only
! place in the tower where the three carriers are constructed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_carriers_fixture

  use graph_carrier, only : counted_set

  implicit none

  private
  public :: time_carriers

contains

  !===================================================================!
  ! The three carriers, declared and nothing more. No direction, no
  ! incidence, no value, no step size, no graph.
  !===================================================================!

  subroutine time_carriers(q, t, e)

    type(counted_set), intent(out) :: q, t, e

    q = counted_set('state coordinates', 2)
    t = counted_set('time instants'    , 5)
    e = counted_set('time steps'       , 4)

  end subroutine time_carriers

end module time_carriers_fixture
