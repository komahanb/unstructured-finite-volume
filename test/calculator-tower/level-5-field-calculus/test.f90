!=====================================================================!
! CALCULATOR TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: WHAT VALUES LIVE ON A DOMAIN. The
! calculator's first numbers arrive - and ONLY the known ones:
!
!      K = { a, b, d }  c-->  X        q(a)=2, q(b)=3, q(d)=4
!      U = { c, e }     c-->  X        recorded, never fabricated
!
! No value is invented for c or e: the known restriction q|_K is
! instantiated, and U only records where computed values will live
! when later levels produce them. K is declared in the NONNUMERIC
! order { d, a, b }, so storage order is domain enumeration and a
! raw-member indexing habit fails loudly here. The import list is
! the negative truth: a field needs a DOMAIN - calculator_assert,
! graph_set, graph_field_calculus, class_graph_field - and no
! graph container anywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_5

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use graph_set    , only : index_set, subset, set
  use graph_field_calculus, only : graph_field
  use class_graph_field, only : field

  implicit none

  type(index_set) :: x
  type(subset)  :: k, u, e
  type(field)       :: qk, qe
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 5 . field calculus"
  write(*,'(1x,a)') "============================================="

  x = index_set('value-slots', 5)

  ! The known support, declared in NONNUMERIC order on purpose.
  k = subset('known'  , x, [SLOT_D, SLOT_A, SLOT_B])
  u = subset('unknown', x, [SLOT_C, SLOT_E])

  call check_supports(nfail)
  call check_known_field(nfail)
  call check_empty_field(nfail)

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! K and U are subobjects of X with exactly the architect-owned
  ! extensions, and their enumerations round-trip.
  !===================================================================!

  subroutine check_supports(nfail)

    integer, intent(inout) :: nfail

    integer :: i
    logical :: ok

    call report(k % size() .eq. 3 .and. u % size() .eq. 2, &
         & "|K| = 3 and |U| = 2", nfail)
    call report(k % has(SLOT_A) .and. k % has(SLOT_B) .and. &
         &      k % has(SLOT_D), &
         & "K holds a, b and d", nfail)
    call report(.not. k % has(SLOT_C) .and. .not. k % has(SLOT_E), &
         & "and neither c nor e", nfail)
    call report(u % has(SLOT_C) .and. u % has(SLOT_E), &
         & "U holds c and e, where answers will one day live", nfail)

    ok = .true.
    do i = 1, k % size()
       ok = ok .and. (k % member(k % local_index(k % member(i))) &
            &         .eq. k % member(i))
    end do
    do i = 1, u % size()
       ok = ok .and. (u % member(u % local_index(u % member(i))) &
            &         .eq. u % member(i))
    end do
    call report(ok, "member and local_index invert on both supports", nfail)

    call report(k % is_subobject_of(x) .and. u % is_subobject_of(x), &
         & "both stand embedded in the value slots", nfail)

  end subroutine check_supports

  !===================================================================!
  ! q|_K with values [4, 2, 3] in the DECLARED order {d, a, b} -
  ! read back through the domain map, never by raw member value.
  !===================================================================!

  subroutine check_known_field(nfail)

    integer, intent(inout) :: nfail

    class(set), allocatable :: dom
    real(dp), allocatable          :: v(:)

    qk = field('q', k)
    call qk % set_real_vector([4.0_dp, 2.0_dp, 3.0_dp])

    call qk % domain(dom)
    call report(dom % equals(k), &
         & "the field's domain is K, by identity", nfail)
    call report(qk % num_entries() .eq. k % size(), &
         & "and its entries count K's members", nfail)

    call qk % get_real_vector(v)
    call report(abs(v(k % local_index(SLOT_A)) - 2.0_dp) < 1.0d-14, &
         & "q(a) = 2, read through the domain map", nfail)
    call report(abs(v(k % local_index(SLOT_B)) - 3.0_dp) < 1.0d-14, &
         & "q(b) = 3", nfail)
    call report(abs(v(k % local_index(SLOT_D)) - 4.0_dp) < 1.0d-14, &
         & "q(d) = 4", nfail)
    call report(abs(v(1) - 4.0_dp) < 1.0d-14, &
         & "storage follows declaration: d's value stands first", nfail)

  end subroutine check_known_field

  !===================================================================!
  ! The empty subset carries a lawful zero-entry field.
  !===================================================================!

  subroutine check_empty_field(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: v(:)

    e  = subset('nowhere', x, [integer ::])
    qe = field('q', e)
    call qe % set_real_vector([real(dp) ::])

    call report(qe % num_entries() .eq. 0, &
         & "a field on the empty subset has zero entries", nfail)
    call qe % get_real_vector(v)
    call report(size(v) .eq. 0, &
         & "and a zero-length vector, no fake member anywhere", nfail)

  end subroutine check_empty_field

end program calculator_level_5
