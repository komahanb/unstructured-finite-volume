!=====================================================================!
! THE PRODUCTION WITNESSES - earned at LEVEL 6, and the first place in
! this tower where production discretization machinery is named at
! all.
!
! Two real stencil operators and the BDF family beside them. Every
! one is built with the ACTUAL production constructor; nothing here is
! a stand-in, a mock, or a re-implementation.
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
!                        NOTHING IS EVER APPLIED
!
! Level 6 is structural introspection. These objects are CONSTRUCTED
! and INTERROGATED, and apply() is never called on any of them - not
! on a stencil, not on a family. The weights and constants exist only
! because the production constructor requires a complete object. A
! family is dimensionless: it takes no step size and no action, and
! its connectivity over a block of instants is read from the family
! alone (src commit b9944c3).
!
!                    THE COORDINATES ARE PRODUCTION'S
!
! stencil(rows, columns, weights, constant) builds
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
  use view_directed    , only : directed_graph
  use operation_stencil, only : stencil
  use operation_family , only : family, bdf_family

  implicit none

  private
  public :: d2_coordinate_stencil, d1_coordinate_stencil
  public :: dependency_pattern, temporal_pattern, bdf2_family
  public :: BDF2_ORDER, FIRST_ORDER

  integer , parameter :: BDF2_ORDER = 2

  ! the equation degree the temporal pattern is read at: a first-order
  ! state, so the block spans the state and its one derivative
  integer , parameter :: FIRST_ORDER = 1

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

  type(stencil) function d2_coordinate_stencil() result(s)

    s = stencil(rows     = [1, 1, 2, 3], &
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

  type(stencil) function d1_coordinate_stencil() result(s)

    s = stencil(rows     = [1, 1, 2, 2, 3], &
         &               columns  = [1, 2, 2, 3, 4], &
         &               weights  = [2.0_dp, -1.0_dp, 0.0_dp, 3.0_dp, 4.0_dp], &
         &               constant = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         &               label    = 'D1 in production coordinates')

  end function d1_coordinate_stencil

  !===================================================================!
  ! The stencil on the dependent axis: the stored pattern, one edge per
  ! coefficient, column -> row. Read from the stencil's own component,
  ! which is what the retired dependencies() accessor returned
  ! (src commit 9fd1732).
  !===================================================================!

  subroutine dependency_pattern(s, pattern)

    type(stencil), intent(in) :: s
    class(directed_graph), allocatable, intent(out) :: pattern

    allocate(pattern, source=s % pattern)

  end subroutine dependency_pattern

  !===================================================================!
  ! BDF2, the production family of that order.
  !===================================================================!

  type(family) function bdf2_family() result(scheme)

    scheme = bdf_family(BDF2_ORDER)

  end function bdf2_family

  !===================================================================!
  ! The stencil on the independent axis: the family's connectivity
  ! over the shortest block its row pattern fits in - history depth
  ! plus one instants, the state and its derivative as the degrees.
  ! For BDF-k the one row it places reads the state at offsets 0..k
  ! into the newest instant: a fan-in, and no succession arrow.
  !===================================================================!

  subroutine temporal_pattern(scheme, pattern)

    type(family), intent(in) :: scheme
    class(directed_graph), allocatable, intent(out) :: pattern

    integer :: num_instants, num_degrees

    num_instants = scheme % history_depth(FIRST_ORDER) + 1
    num_degrees  = FIRST_ORDER + 1

    allocate(pattern, source=scheme % block_connectivity(num_degrees, num_instants))

  end subroutine temporal_pattern

end module production_discretization_fixture
