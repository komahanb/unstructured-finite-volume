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
! This file is the only place in the tower where the three carriers are
! constructed, and so the only place granted a representation. They
! ESCAPE it, so the map that says what they contain escapes with them:
! a domain nobody described is a domain nobody can ask anything of.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_carriers_fixture

  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation
  use map_set           , only : set_map

  implicit none

  private
  public :: time_carriers

contains

  !===================================================================!
  ! The three carriers, declared and nothing more. No direction, no
  ! incidence, no value, no step size, no graph.
  !===================================================================!

  subroutine time_carriers(sets, q, t, e)

    type(set_map)  , intent(inout) :: sets
    type(graph), intent(out)   :: q, t, e

    call q % declare()      ! state coordinates
    call t % declare()      ! time instants
    call e % declare()      ! time steps

    call sets % bind(q, counted_set_representation(2))
    call sets % bind(t, counted_set_representation(5))
    call sets % bind(e, counted_set_representation(4))

  end subroutine time_carriers

end module time_carriers_fixture
