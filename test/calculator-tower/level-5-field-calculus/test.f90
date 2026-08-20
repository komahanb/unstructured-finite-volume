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
! the set modules, field_calculus, field_stored - and no
! graph container anywhere.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_5

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use field_calculus, only : field
  use field_stored, only : stored_field

  implicit none

  type(graph) :: x
  type(graph)  :: k, u, e
  type(stored_field)       :: qk, qe
  integer           :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 5 . field calculus"
  write(*,'(1x,a)') "============================================="

  call x % declare()
  call sets % bind(x, counted_set_representation(5))

  ! The known support, declared in NONNUMERIC order on purpose.
  call k % declare()
  call sets       % bind(k, listed_set_representation([SLOT_D, SLOT_A, SLOT_B]))
  call inclusions % include_in(k, x)
  call u % declare()
  call sets       % bind(u, listed_set_representation([SLOT_C, SLOT_E]))
  call inclusions % include_in(u, x)

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

    call report(sets % num_members_of(k) .eq. 3 .and. sets % num_members_of(u) .eq. 2, &
         & "|K| = 3 and |U| = 2", nfail)
    call report(sets % has(k, SLOT_A) .and. sets % has(k, SLOT_B) .and. &
         &      sets % has(k, SLOT_D), &
         & "K holds a, b and d", nfail)
    call report(.not. sets % has(k, SLOT_C) .and. .not. sets % has(k, SLOT_E), &
         & "and neither c nor e", nfail)
    call report(sets % has(u, SLOT_C) .and. sets % has(u, SLOT_E), &
         & "U holds c and e, where answers will one day live", nfail)

    ok = .true.
    do i = 1, sets % num_members_of(k)
       ok = ok .and. (sets % member_of(k, sets % index_in(k, sets % member_of(k, i))) &
            &         .eq. sets % member_of(k, i))
    end do
    do i = 1, sets % num_members_of(u)
       ok = ok .and. (sets % member_of(u, sets % index_in(u, sets % member_of(u, i))) &
            &         .eq. sets % member_of(u, i))
    end do
    call report(ok, "member and local_index invert on both supports", nfail)

    call report(declared_subobject(k, x, inclusions) .and. declared_subobject(u, x, inclusions), &
         & "both stand embedded in the value slots", nfail)

  end subroutine check_supports

  !===================================================================!
  ! q|_K with values [4, 2, 3] in the DECLARED order {d, a, b} -
  ! read back through the domain map, never by raw member value.
  !===================================================================!

  subroutine check_known_field(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: dom
    real(dp), allocatable          :: v(:)

    qk = stored_field('q', k, sets % num_members_of(k))
    call qk % set_real_vector([4.0_dp, 2.0_dp, 3.0_dp])

    dom = qk % domain()
    call report(dom % same_as(k), &
         & "the field's domain is K, by identity", nfail)
    call report(qk % num_entries() .eq. sets % num_members_of(k), &
         & "and its entries count K's members", nfail)

    call qk % real_vector(v)
    call report(abs(v(sets % index_in(k, SLOT_A)) - 2.0_dp) < 1.0d-14, &
         & "q(a) = 2, read through the domain map", nfail)
    call report(abs(v(sets % index_in(k, SLOT_B)) - 3.0_dp) < 1.0d-14, &
         & "q(b) = 3", nfail)
    call report(abs(v(sets % index_in(k, SLOT_D)) - 4.0_dp) < 1.0d-14, &
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

    call e % declare()
    call sets       % bind(e, listed_set_representation([integer ::]))
    call inclusions % include_in(e, x)
    qe = stored_field('q', e, sets % num_members_of(e))
    call qe % set_real_vector([real(dp) ::])

    call report(qe % num_entries() .eq. 0, &
         & "a field on the empty subset has zero entries", nfail)
    call qe % real_vector(v)
    call report(size(v) .eq. 0, &
         & "and a zero-length vector, no fake member anywhere", nfail)

  end subroutine check_empty_field

end program calculator_level_5
