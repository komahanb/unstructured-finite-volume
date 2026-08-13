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
  use graph_carrier    , only : counted_set, subset_set, member_set
  use class_graph_field, only : field

  implicit none

  type(counted_set) :: v
  type(subset_set)  :: x_dom, c
  type(field)       :: qx
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "derivative action tower . level 5 . fields"
  write(*,'(1x,a)') "============================================="

  v = counted_set('value-slots', 4)

  ! The two roles, declared in deliberately nonsemantic order.
  x_dom = subset_set('independent', v, [SLOT_Y, SLOT_X])
  c     = subset_set('computed'   , v, [SLOT_U, SLOT_Z])

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

    call report(x_dom % size() .eq. 2 .and. x_dom % has(SLOT_Y) .and. &
         &      x_dom % has(SLOT_X) .and. .not. x_dom % has(SLOT_U), &
         & "X = { y, x }, the independent inputs", nfail)
    call report(c % size() .eq. 2 .and. c % has(SLOT_U) .and. &
         &      c % has(SLOT_Z) .and. .not. c % has(SLOT_X), &
         & "C = { u, z }, where answers will one day live", nfail)

    call report(x_dom % is_subobject_of(v) .and. &
         &      c % is_subobject_of(v), &
         & "both stand embedded in the value slots", nfail)

    ok = .true.
    do i = 1, v % size()
       m = v % member(i)
       homes = count([x_dom % has(m), c % has(m)])
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

    class(member_set), allocatable :: dom
    real(dp), allocatable          :: val(:)

    qx = field('base point', x_dom)
    call qx % set_real_vector([3.0_dp, 2.0_dp])

    call qx % domain(dom)
    call report(dom % same_as(x_dom), &
         & "the base field's domain is X, by identity", nfail)
    call report(qx % num_entries() .eq. 2 .and. &
         &      qx % num_components() .eq. 1, &
         & "two entries, one component", nfail)

    call qx % get_real_vector(val)
    call report(abs(val(x_dom % local_index(SLOT_Y)) - 3.0_dp) &
         &      < 1.0d-14, &
         & "y = 3, read through the domain map", nfail)
    call report(abs(val(x_dom % local_index(SLOT_X)) - 2.0_dp) &
         &      < 1.0d-14, &
         & "x = 2", nfail)
    call report(abs(val(1) - 3.0_dp) < 1.0d-14 .and. &
         &      abs(val(2) - 2.0_dp) < 1.0d-14, &
         & "storage follows declaration: y's value stands first", nfail)

  end subroutine check_base_point

end program derivative_level_5
