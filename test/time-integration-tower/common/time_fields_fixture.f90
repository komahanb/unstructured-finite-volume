!=====================================================================!
! THE TIME FIELDS FIXTURE - earned at LEVEL 5.
!
! Three fields on three different domains, and the point of the
! file is that they are three DIFFERENT domains:
!
!      q0   : Q -> reals  [2, 0]                 the state
!      time : T -> reals  [0, 1/2, 1, 3/2, 2]    the instants' coords
!      h    : E -> reals  [1/2, 1/2, 1/2, 1/2]   the step sizes
!
! THE STATE FIELD NEEDS NO GRAPH. A field is a function over one
! member set - that is the field ontology as production already
! states it - so q0 lives on Q directly. No graph is constructed
! here, and none is needed: manufacturing a two-vertex graph so
! that q0 could live on its vertices would be inventing the very
! conflation this tower exists to refuse, one level before the
! level that tests it.
!
!                    VALUES ARE NOT STRUCTURE
!
! The instants t0..t4 are MEMBERS of T. The reals 0, 1/2, 1, 3/2, 2
! are VALUES of a field on T. Those are different things, and the
! tower keeps four objects apart where a looser client would keep
! one:
!
!      T                  the carrier of instants
!      Tail/Head/A1/A2    the temporal structure over it
!      time : T -> reals  numerical coordinates
!      h    : E -> reals  numerical step sizes
!
! The consistency between the last two and the structure is a
! THEOREM this fixture's consumer proves, not an assumption it
! makes:
!
!      time(head(e)) - time(tail(e)) = h(e)      for every step e
!
! which is field calculus over already-earned structure, and is not
! a time scheme. Nothing here knows about Euler, BDF, a residual or
! a solver.
!
! The carriers arrive as arguments, as they have since Level 1.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_fields_fixture

  use iso_fortran_env  , only : dp => REAL64
  use graph_fractal    , only : graph
  use field_stored, only : stored_field
  use time_assert      , only : NQ, NT, NE, H_STEP, TIME_COORD, Q0

  implicit none

  private
  public :: state_field, instant_coordinates, step_sizes

contains

  !===================================================================!
  ! q0 : Q -> reals, and NOT a graph in sight.
  !===================================================================!

  type(stored_field) function state_field(q) result(f)

    type(graph), intent(in) :: q

    f = stored_field('state', q, NQ, num_components=1)
    call f % set_real_vector(Q0)

  end function state_field

  !===================================================================!
  ! time : T -> reals, the numerical coordinate of each instant.
  !===================================================================!

  type(stored_field) function instant_coordinates(t) result(f)

    type(graph), intent(in) :: t

    f = stored_field('instant coordinate', t, NT, num_components=1)
    call f % set_real_vector(TIME_COORD)

  end function instant_coordinates

  !===================================================================!
  ! h : E -> reals, one step size per step. A FIELD ON THE STEPS, not a
  ! scalar hidden in a scheme: the tower will use a uniform h, and
  ! the uniformity is a property of these values rather than of the
  ! type that holds them.
  !===================================================================!

  type(stored_field) function step_sizes(e) result(f)

    type(graph), intent(in) :: e

    real(dp) :: values(NE)

    values = H_STEP
    f = stored_field('step size', e, NE, num_components=1)
    call f % set_real_vector(values)

  end function step_sizes

end module time_fields_fixture
