!=====================================================================!
! DERIVATIVE ACTION TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: WHAT PRIMAL VALUES EXIST, AND
! WHICH SLOTS ARE COMPUTED RATHER THAN SEEDED. The base point
! enters - and each number has exactly one home:
!
!      X = { y, x }     independent inputs     [3, 2]
!      C = { u, z }     computed later         NO FIELD AT ALL
!
! The declaration order X = {y, x} is deliberately nonsemantic:
! storage follows domain enumeration - y's value stands first - and
! every read walks local_index. C carries no field, no zero, no
! NaN, no sentinel: u and z are not zero, they are UNCOMPUTED, and
! not constructing a field is a stronger truth than constructing an
! empty one.
!
! And the derivative negative truths are this rung's own: NO
! tangent field, NO cotangent field, NO seed, NO derivative value.
! Primal state uses the ordinary field on an ordinary domain -
! derivative vectors have not yet earned a separate type, and Gate
! A never asks them to. The operation laws are still unspoken:
! knowing x = 2 is not knowing what product does with it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program derivative_level_5

  use iso_fortran_env  , only : dp => REAL64
  use derivative_assert, only : report, verdict
  use derivative_assert, only : SLOT_X, SLOT_Y, SLOT_U, SLOT_Z
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use field_stored, only : stored_field

  implicit none

  type(graph) :: v
  type(graph)  :: x_dom, c
  type(stored_field)       :: qx
  integer           :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 5 . fields"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(4))

  ! The two roles, declared in deliberately nonsemantic order.
  call x_dom % declare()
  call sets       % bind(x_dom, listed_set_representation([SLOT_Y, SLOT_X]))
  call inclusions % include_in(x_dom, v)
  call c % declare()
  call sets       % bind(c, listed_set_representation([SLOT_U, SLOT_Z]))
  call inclusions % include_in(c, v)

  call check_partition(nfail)
  call check_base_point(nfail)

  ! C deliberately receives NO field here: no q_C, no zeros, no
  ! sentinels. And no tangent, cotangent, or seed exists anywhere:
  ! Gate A has primal base data and no derivative arithmetic.

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! Extensions, embeddings, and the one-home law: every member of V
  ! belongs to exactly one of X, C - disjointness and coverage
  ! proved together, composed locally from membership.
  !===================================================================!

  subroutine check_partition(nfail)

    integer, intent(inout) :: nfail

    integer :: i, m, homes
    logical :: ok

    call report(sets % num_members_of(x_dom) .eq. 2 .and. sets % has(x_dom, SLOT_Y) .and. &
         &      sets % has(x_dom, SLOT_X) .and. .not. sets % has(x_dom, SLOT_U), &
         & "X = { y, x }, the independent inputs", nfail)
    call report(sets % num_members_of(c) .eq. 2 .and. sets % has(c, SLOT_U) .and. &
         &      sets % has(c, SLOT_Z) .and. .not. sets % has(c, SLOT_X), &
         & "C = { u, z }, where answers will one day live", nfail)

    call report(declared_subobject(x_dom, v, inclusions) .and. &
         &      declared_subobject(c, v, inclusions), &
         & "both stand embedded in the value slots", nfail)

    ok = .true.
    do i = 1, sets % num_members_of(v)
       m = sets % member_of(v, i)
       homes = count([sets % has(x_dom, m), sets % has(c, m)])
       ok = ok .and. (homes .eq. 1)
    end do
    call report(ok, &
         & "every slot has exactly one home: disjoint, and covering", &
         & nfail)

  end subroutine check_partition

  !===================================================================!
  ! The base point: [3, 2] on X = { y, x } - storage follows the
  ! DECLARATION, y's value first, and every read walks local_index.
  !===================================================================!

  subroutine check_base_point(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: dom
    real(dp), allocatable          :: val(:)

    qx = stored_field('base point', x_dom, sets % num_members_of(x_dom))
    call qx % set_real_vector([3.0_dp, 2.0_dp])

    dom = qx % domain()
    call report(dom % same_as(x_dom), &
         & "the base field's domain is X, by identity", nfail)
    call report(qx % num_entries() .eq. 2 .and. &
         &      qx % num_components() .eq. 1, &
         & "two entries, one component", nfail)

    call qx % real_vector(val)
    call report(abs(val(sets % index_in(x_dom, SLOT_Y)) - 3.0_dp) &
         &      < 1.0d-14, &
         & "y = 3, read through the domain map", nfail)
    call report(abs(val(sets % index_in(x_dom, SLOT_X)) - 2.0_dp) &
         &      < 1.0d-14, &
         & "x = 2", nfail)
    call report(abs(val(1) - 3.0_dp) < 1.0d-14 .and. &
         &      abs(val(2) - 2.0_dp) < 1.0d-14, &
         & "storage follows declaration: y's value stands first", nfail)

  end subroutine check_base_point

end program derivative_level_5
