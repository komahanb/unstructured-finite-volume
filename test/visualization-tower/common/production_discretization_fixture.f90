!=====================================================================!
! THE PRODUCTION WITNESSES - earned at LEVEL 6, and the first place in
! this tower where production discretization machinery is named at
! all.
!
! Three real stencil_operators and the step_operators that wrap them.
! Every one is built with the ACTUAL production constructor; nothing
! here is a stand-in, a mock, or a re-implementation.
!
!                    WHY THESE THREE, AND NOT OTHERS
!
!   d2_coordinate_stencil    the same Boolean occupancy as the
!                            tower's D2 : X1 -> X2, expressed in
!                            production's own coordinates, and
!                            carrying Level 5's w2 = [1, 5, -2, 2]
!                            as its weights. The primary probe.
!
!   d1_coordinate_stencil    the RECTANGULAR witness. D1 : X0 -> X1
!                            runs 4 -> 3, and production's
!                            constructor takes ONE vertex count. What
!                            it can and cannot say about that is the
!                            level's second question.
!
!   diagonal_stencil         a second wrapped action whose state
!                            sparsity differs from the first, used
!                            once to find out whether a step's
!                            dependency answer depends on what it
!                            wraps.
!
!                        NOTHING IS EVER APPLIED
!
! Level 6 is structural introspection. These objects are CONSTRUCTED
! and INTERROGATED, and apply() is never called on any of them - not
! on a stencil, not on a step. The weights and constants exist only
! because the production constructor requires a complete object, and
! the step size exists only because bdf requires one.
!
!                    THE COORDINATES ARE PRODUCTION'S
!
! stencil_operator(rows, columns, weights, constant) builds
!
!      stored_graph(nv, tails = columns, heads = rows)
!
! with nv = size(constant). So an edge runs COLUMN -> ROW, and the
! whole object stands on a single vertex count. Both facts are read
! off the production source rather than assumed, and both are what
! Level 6 measures.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module production_discretization_fixture

  use iso_fortran_env    , only : dp => REAL64
  use operation_action, only : graph_operation
  use class_graph_stencil, only : stencil_operator
  use class_graph_step   , only : step_operator, bdf

  implicit none

  private
  public :: d2_coordinate_stencil, d1_coordinate_stencil
  public :: diagonal_stencil, bdf2_around
  public :: H_PROBE, BDF2_ORDER

  !-------------------------------------------------------------------!
  ! A harmless nonzero step size, and the order. Neither is ever used
  ! arithmetically here - bdf simply refuses to exist without them.
  !-------------------------------------------------------------------!

  real(dp), parameter :: H_PROBE    = 0.5_dp
  integer , parameter :: BDF2_ORDER = 2

contains

  !===================================================================!
  ! D2's occupancy in production coordinates.
  !
  !      D2 = { p->u, q->u, q->v, r->w }
  !
  ! becomes, by declaration position,
  !
  !      (column, row) = (1,1), (2,1), (2,2), (3,3)
  !
  ! and the weights are Level 5's w2, unchanged. The constant vector's
  ! LENGTH is what sets the vertex count - three, here.
  !===================================================================!

  type(stencil_operator) function d2_coordinate_stencil() result(s)

    s = stencil_operator(rows     = [1, 1, 2, 3], &
         &               columns  = [1, 2, 2, 3], &
         &               weights  = [1.0_dp, 5.0_dp, -2.0_dp, 2.0_dp], &
         &               constant = [0.0_dp, 0.0_dp, 0.0_dp], &
         &               label    = 'D2 in production coordinates')

  end function d2_coordinate_stencil

  !===================================================================!
  ! D1's occupancy in production coordinates - the rectangular one.
  !
  !      D1 = { a->p, b->p, b->q, c->q, d->r }
  !
  !      (column, row) = (1,1), (2,1), (2,2), (3,2), (4,3)
  !
  ! The columns run to 4 and the rows only to 3, and the constructor
  ! takes ONE count for both. Four is the only count that can hold the
  ! columns, so four is what it is given - and what that costs is the
  ! level's second finding, not an assumption made here.
  !
  ! The weights are Level 5's w1, zero and all.
  !===================================================================!

  type(stencil_operator) function d1_coordinate_stencil() result(s)

    s = stencil_operator(rows     = [1, 1, 2, 2, 3], &
         &               columns  = [1, 2, 2, 3, 4], &
         &               weights  = [2.0_dp, -1.0_dp, 0.0_dp, 3.0_dp, 4.0_dp], &
         &               constant = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         &               label    = 'D1 in production coordinates')

  end function d1_coordinate_stencil

  !===================================================================!
  ! A second action with a DIFFERENT state sparsity - the identity
  ! pattern, three vertices, three self-edges and nothing else.
  !
  ! Its only job is to be wrapped, so that two steps differing in
  ! everything except their order can be asked the same question.
  !===================================================================!

  type(stencil_operator) function diagonal_stencil() result(s)

    s = stencil_operator(rows     = [1, 2, 3], &
         &               columns  = [1, 2, 3], &
         &               weights  = [1.0_dp, 1.0_dp, 1.0_dp], &
         &               constant = [0.0_dp, 0.0_dp, 0.0_dp], &
         &               label    = 'a diagonal action')

  end function diagonal_stencil

  !===================================================================!
  ! BDF2 around whatever action it is handed. The production shelf,
  ! called as production calls it.
  !===================================================================!

  type(step_operator) function bdf2_around(action) result(clock)

    class(graph_operation), intent(in) :: action

    clock = bdf(BDF2_ORDER, action, H_PROBE)

  end function bdf2_around

end module production_discretization_fixture
