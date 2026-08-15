!=====================================================================!
! VISUALIZATION TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: WHAT CHANGES WHEN NUMERICAL VALUES
! ARE ATTACHED TO DEPENDENCY OCCURRENCES.
!
! The answer is: NOTHING STRUCTURAL CHANGES.
!
!      w1 : E1 -> reals = [ 2, -1,  0,  3,  4 ]
!      w2 : E2 -> reals = [ 1,  5, -2,  2 ]
!      w3 : E3 -> reals = [ 3, -1,  4 ]
!
! and after all three exist, D1, D2, D3, D2:1 and D3:1 hold exactly
! the tuples Gate A derived, and the Level-4 pictures are the same
! strings character for character.
!
!                        THE LOAD-BEARING ASSERTION
!
!      e13 : b -> q        w1(e13) = 0        b -> q IS IN D1
!
! One occurrence, one zero, and the dependency survives it. That is
! the whole level in one line, and everything else here exists to keep
! it from being an accident.
!
!                          ZERO IS NOT ABSENCE
!
! In the coefficient picture:
!
!      #   structurally present            (the Level-4 view)
!      0   present, carrying zero
!      .   ABSENT - there is no dependency, so there is nothing to
!          carry a value
!
! (b,q) prints 0 and (a,q) prints '.', and they are different facts.
! A representation that wrote 0 for both would have thrown away the
! distinction this level was built to establish.
!
!                    VALUE MAY NEVER ANSWER A STRUCTURAL QUESTION
!
! The forbidden inference is
!
!      if (coefficient /= 0) then present
!
! and it is forbidden in both directions. The last check below runs a
! SECOND field on the same E1 - w1_alt = [9, 8, 7, 6, 5], sharing no
! value with the first and containing no zero - and finds the
! structural picture identical while the coefficient picture is not.
! So D1 does not determine w1, and w1 does not determine D1.
!
!                       WHAT IS DELIBERATELY ABSENT
!
! No coefficient is computed for D2:1 or D3:1. Composing values would
! mean choosing an algebra for numerical composition - sums over
! intermediate members - which is operator mathematics and belongs to
! no level yet built. No w^T exists either; transposing a coefficient
! is a numerical act, and Gate A's structural reverse asked for none.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_5

  use iso_fortran_env      , only : dp => REAL64
  use visualization_assert , only : report, verdict
  use visualization_assert , only : ND1, ND2, ND3, ND21, ND31
  use visualization_assert , only : NX0, NE2
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use visualization_assert , only : E1_1, E1_2, E1_3, E1_4, E1_5
  use visualization_assert , only : E2_1, E2_4, E3_1, E3_2, E3_3
  use graph_carrier        , only : counted_set, member_set
  use graph_relation       , only : relation
  use graph_binary_relation, only : csr_relation
  use graph_field_calculus , only : GRAPH_FIELD_REAL
  use class_graph_field    , only : field
  use visualization_carriers_fixture , only : structural_carriers, label_for
  use visualization_relations_fixture, only : occurrences_of_a1
  use visualization_relations_fixture, only : occurrences_of_a2
  use visualization_relations_fixture, only : occurrences_of_a3
  use visualization_algebra_fixture  , only : derive_dependency
  use visualization_algebra_fixture  , only : derive_composition
  use visualization_algebra_fixture  , only : materialized_transpose
  use visualization_algebra_fixture  , only : same_extension
  use structural_renderer_fixture    , only : picture, sparsity_picture
  use structural_renderer_fixture    , only : glyph_at
  use visualization_values_fixture   , only : coefficients_of_a1
  use visualization_values_fixture   , only : coefficients_of_a2
  use visualization_values_fixture   , only : coefficients_of_a3
  use visualization_values_fixture   , only : alternate_coefficients_of_a1
  use visualization_values_fixture   , only : COEFF_A1, COEFF_A2, COEFF_A3
  use visualization_values_fixture   , only : COEFF_A1_ALT, TOL
  use valued_renderer_fixture        , only : valued_sparsity_picture
  use valued_renderer_fixture        , only : coefficients_fit, value_at
  use valued_renderer_fixture        , only : occurrence_joining
  use valued_renderer_fixture        , only : occurrences_joining
  use valued_renderer_fixture        , only : ABSENT, value_token

  implicit none

  type(counted_set)          :: x0, x1, x2, x3, e1, e2, e3
  type(csr_relation)         :: t1, h1, t2, h2, t3, h3
  type(csr_relation), target :: d1, d2, d3, d21, d31
  type(field)                :: w1, w2, w3, w1_alt, decoy
  type(picture)              :: before_d1, before_d2, before_d3
  integer                    :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 5 . field calculus"
  write(*,'(1x,a)') "============================================="

  ! ---- the persistent structure, exactly as Gate A left it.
  call structural_carriers(x0, x1, x2, x3, e1, e2, e3)
  call occurrences_of_a1(e1, x0, x1, t1, h1)
  call occurrences_of_a2(e2, x1, x2, t2, h2)
  call occurrences_of_a3(e3, x2, x3, t3, h3)

  d1 = derive_dependency('D1', t1, h1)
  d2 = derive_dependency('D2', t2, h2)
  d3 = derive_dependency('D3', t3, h3)

  ! ---- the structural pictures, taken BEFORE any number exists.
  before_d1 = sparsity_picture(d1)
  before_d2 = sparsity_picture(d2)
  before_d3 = sparsity_picture(d3)

  ! ---- and now the numbers arrive.
  w1     = coefficients_of_a1(e1)
  w2     = coefficients_of_a2(e2)
  w3     = coefficients_of_a3(e3)
  w1_alt = alternate_coefficients_of_a1(e1)

  ! A perfectly valid field that is simply not about E2. |X0| = |E2|.
  decoy = field('four numbers on X0', x0, ncomp=1)
  call decoy % set_real_vector([7.0_dp, 7.0_dp, 7.0_dp, 7.0_dp])

  ! ---- the compositions, derived AFTER the fields exist.
  d21 = derive_composition('D2:1', d1,  d2)
  d31 = derive_composition('D3:1', d21, d3)

  call say_the_structure_and_the_values()

  call check_the_fields_exist(nfail)
  call check_the_exact_values(nfail)
  call check_the_domains_by_identity(nfail)
  call check_the_occurrence_seat_is_unique(nfail)
  call check_the_zero_witness(nfail)
  call check_zero_is_not_absence(nfail)
  call check_the_two_views_agree_on_presence(nfail)
  call check_the_structural_pictures_are_unchanged(nfail)
  call check_the_derivation_is_unchanged(nfail)
  call check_value_and_structure_are_independent(nfail)

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! What the level produced, printed for a person to look at: the
  ! structural view and the coefficient view of the same relation, on
  ! the same page, row for row.
  !===================================================================!

  subroutine say_the_structure_and_the_values()

    type(picture) :: pic

    write(*,'(1x,a)') "---------------------------------------------"

    write(*,'(4x,a)') "D1 STRUCTURE"
    pic = sparsity_picture(d1); call say_grid(pic)

    write(*,'(4x,a)') "D1 VALUES"
    pic = valued_sparsity_picture(d1, t1, h1, e1, w1); call say_grid(pic)

    write(*,'(4x,a)') "D2 VALUES"
    pic = valued_sparsity_picture(d2, t2, h2, e2, w2); call say_grid(pic)

    write(*,'(4x,a)') "D3 VALUES"
    pic = valued_sparsity_picture(d3, t3, h3, e3, w3); call say_grid(pic)

    write(*,'(4x,a)') "ZERO IS NOT ABSENCE"
    call say_one_cell(X0_B, X1_Q)
    call say_one_cell(X0_A, X1_Q)
    write(*,*)

    write(*,'(4x,a)') "THE SAME STRUCTURE, DIFFERENT NUMBERS"
    pic = valued_sparsity_picture(d1, t1, h1, e1, w1_alt); call say_grid(pic)

    write(*,'(1x,a)') "---------------------------------------------"

  end subroutine say_the_structure_and_the_values

  !-------------------------------------------------------------------!
  ! The grid alone. A picture carries its own name line, and the
  ! headings above already say which relation this is.
  !-------------------------------------------------------------------!

  subroutine say_grid(pic)

    type(picture), intent(in) :: pic

    integer :: k

    do k = 2, pic % rows()
       write(*,'(4x,a)') pic % at(k)
    end do
    write(*,*)

  end subroutine say_grid

  !-------------------------------------------------------------------!
  ! One cell of D1, in the reader's own names, said in full.
  !-------------------------------------------------------------------!

  subroutine say_one_cell(from, to)

    integer, intent(in) :: from, to

    integer :: e

    write(*,'(6x,a)') "(" // label_for(x0, from) // "," // &
         &            label_for(x1, to) // ") in X0 x X1:"

    if (d1 % has([from, to])) then
       e = occurrence_joining(t1, h1, e1, from, to)
       write(*,'(10x,a)') "structural = PRESENT"
       write(*,'(10x,a)') "occurrence = " // label_for(e1, e)
       write(*,'(10x,a)') "value      = " // value_token(value_at(w1, e1, e))
    else
       write(*,'(10x,a)') "structural = ABSENT"
       write(*,'(10x,a)') "value      = " // ABSENT // &
            &             "   (nothing there to carry one)"
    end if

  end subroutine say_one_cell

  !===================================================================!
  ! Three real scalar fields, on the three occurrence carriers.
  !===================================================================!

  subroutine check_the_fields_exist(nfail)

    integer, intent(inout) :: nfail

    call report(w1 % num_entries() .eq. e1 % size() .and. &
         &      w2 % num_entries() .eq. e2 % size() .and. &
         &      w3 % num_entries() .eq. e3 % size(), &
         & "w1, w2, w3 have one entry per occurrence: 5, 4, 3 - the " // &
         & "field takes its count from the domain it lives on", nfail)

    call report(w1 % num_components() .eq. 1 .and. &
         &      w2 % num_components() .eq. 1 .and. &
         &      w3 % num_components() .eq. 1, &
         & "each carries one number per occurrence", nfail)

    call report(w1 % value_kind() .eq. GRAPH_FIELD_REAL .and. &
         &      w2 % value_kind() .eq. GRAPH_FIELD_REAL .and. &
         &      w3 % value_kind() .eq. GRAPH_FIELD_REAL, &
         & "and all three are REAL fields - THE FIRST NUMBERS IN " // &
         & "THIS TOWER", nfail)

  end subroutine check_the_fields_exist

  !===================================================================!
  ! All twelve coefficients, recovered BY OCCURRENCE MEMBER - through
  ! the domain's own local_index, never by assuming a member is a
  ! position.
  !===================================================================!

  subroutine check_the_exact_values(nfail)

    integer, intent(inout) :: nfail

    call report(near(value_at(w1, e1, E1_1), COEFF_A1(1)) .and. &
         &      near(value_at(w1, e1, E1_2), COEFF_A1(2)) .and. &
         &      near(value_at(w1, e1, E1_3), COEFF_A1(3)) .and. &
         &      near(value_at(w1, e1, E1_4), COEFF_A1(4)) .and. &
         &      near(value_at(w1, e1, E1_5), COEFF_A1(5)), &
         & "w1 = [2, -1, 0, 3, 4] on e11..e15", nfail)

    call report(near(value_at(w2, e2, E2_1), COEFF_A2(1)) .and. &
         &      near(value_at(w2, e2, E2_4), COEFF_A2(4)) .and. &
         &      all_values_recovered(w2, e2, COEFF_A2), &
         & "w2 = [1, 5, -2, 2] on e21..e24", nfail)

    call report(all_values_recovered(w3, e3, COEFF_A3), &
         & "w3 = [3, -1, 4] on e31..e33", nfail)

    call report(all_values_recovered(w1, e1, COEFF_A1) .and. &
         &      all_values_recovered(w1_alt, e1, COEFF_A1_ALT), &
         & "twelve coefficients recovered by occurrence member, and " // &
         & "five more for the probe", nfail)

  end subroutine check_the_exact_values

  !===================================================================!
  ! THE COEFFICIENTS LIVE ON E, AND NOWHERE ELSE.
  !
  ! Checked by identity. |X0| = |E2| = 4 in this specimen, so a field
  ! of exactly the right length exists on the wrong carrier - and is
  ! refused, because same_as is not a size comparison.
  !===================================================================!

  subroutine check_the_domains_by_identity(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: on

    call w1 % domain(on)
    call report(on % same_as(e1) .and. .not. on % same_as(x0) .and. &
         &      .not. on % same_as(x1), &
         & "domain(w1) IS E1 - not X0 which it reads, and not X1 " // &
         & "which it writes", nfail)

    call w2 % domain(on)
    call report(on % same_as(e2), "domain(w2) is E2", nfail)
    call w3 % domain(on)
    call report(on % same_as(e3), "domain(w3) is E3", nfail)

    call report(coefficients_fit(w1, e1) .and. coefficients_fit(w2, e2) .and. &
         &      coefficients_fit(w3, e3), &
         & "and each field fits the occurrences it was made for", nfail)

    ! The hostile case the specimen supplies for free.
    call report(decoy % num_entries() .eq. NX0 .and. NX0 .eq. NE2, &
         & "a four-valued field on X0 exists and is perfectly valid - " // &
         & "|X0| = |E2| = 4", nfail)

    call report(.not. coefficients_fit(decoy, e2), &
         & "AND IT IS REFUSED AS E2's COEFFICIENTS - by domain " // &
         & "IDENTITY, not by cardinality, which would have accepted it", &
         & nfail)

    call report(.not. coefficients_fit(w1, e2) .and. &
         &      .not. coefficients_fit(w2, e1) .and. &
         &      .not. coefficients_fit(decoy, e1), &
         & "and no field is mistaken for another operator's " // &
         & "coefficients", nfail)

  end subroutine check_the_domains_by_identity

  !===================================================================!
  ! A coefficient picture is only well defined where the seat is
  ! unique. For the three direct dependencies it is: exactly one
  ! occurrence joins each tuple, and none joins a tuple that is
  ! absent.
  !===================================================================!

  subroutine check_the_occurrence_seat_is_unique(nfail)

    integer, intent(inout) :: nfail

    call report(seats_are_unique(d1, t1, h1, e1) .and. &
         &      seats_are_unique(d2, t2, h2, e2) .and. &
         &      seats_are_unique(d3, t3, h3, e3), &
         & "every tuple of D1, D2, D3 is joined by EXACTLY ONE " // &
         & "occurrence, and every absent pair by none", nfail)

    call report(occurrence_joining(t1, h1, e1, X0_B, X1_Q) .eq. E1_3 .and. &
         &      occurrence_joining(t1, h1, e1, X0_A, X1_Q) .eq. 0, &
         & "b->q is seated at e13; a->q is seated nowhere, because " // &
         & "a->q is not a dependency", nfail)

  end subroutine check_the_occurrence_seat_is_unique

  !===================================================================!
  ! THE LOAD-BEARING ASSERTION.
  !===================================================================!

  subroutine check_the_zero_witness(nfail)

    integer, intent(inout) :: nfail

    call report(t1 % has([E1_3, X0_B]) .and. h1 % has([E1_3, X1_Q]), &
         & "e13 reads b and writes q - T1(e13) = b, H1(e13) = q", nfail)

    call report(abs(value_at(w1, e1, E1_3)) .lt. TOL, &
         & "w1(e13) = 0 EXACTLY", nfail)

    call report(d1 % has([X0_B, X1_Q]), &
         & "AND b->q IS STILL IN D1 - a zero coefficient does not " // &
         & "remove a dependency", nfail)

    call report(d1 % num_tuples() .eq. ND1, &
         & "D1 still holds five tuples, not four: the zero cost it " // &
         & "nothing", nfail)

    call report(glyph_at(d1, X0_B, X1_Q) .eq. '#', &
         & "and the STRUCTURAL picture still writes # at (b,q), " // &
         & "because it never asked what the value was", nfail)

  end subroutine check_the_zero_witness

  !===================================================================!
  ! ZERO IS NOT ABSENCE, in the representation.
  !===================================================================!

  subroutine check_zero_is_not_absence(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: pic

    pic = valued_sparsity_picture(d1, t1, h1, e1, w1)

    call report(pic % at(1) .eq. 'D1 VALUES' .and. &
         &      pic % at(2) .eq. '           a    b    c    d' .and. &
         &      pic % at(3) .eq. 'p          2   -1    .    .' .and. &
         &      pic % at(4) .eq. 'q          .    0    3    .' .and. &
         &      pic % at(5) .eq. 'r          .    .    .    4', &
         & "D1's coefficient picture: 0 stands at (b,q) and '.' at " // &
         & "(a,q) - ZERO IS NOT ABSENCE", nfail)

    call report(index(pic % at(4), '0') .gt. 0 .and. &
         &      index(pic % at(4), ABSENT) .gt. 0 .and. &
         &      index(pic % at(4), '0') .ne. index(pic % at(4), ABSENT), &
         & "the q row carries BOTH a zero and a dot, in different " // &
         & "columns: they are different facts", nfail)

    pic = valued_sparsity_picture(d2, t2, h2, e2, w2)
    call report(pic % at(2) .eq. '           p    q    r' .and. &
         &      pic % at(3) .eq. 'u          1    5    .' .and. &
         &      pic % at(4) .eq. 'v          .   -2    .' .and. &
         &      pic % at(5) .eq. 'w          .    .    2', &
         & "D2's coefficient picture carries [1, 5, -2, 2] at its " // &
         & "four seats", nfail)

    pic = valued_sparsity_picture(d3, t3, h3, e3, w3)
    call report(pic % at(2) .eq. '           u    v    w' .and. &
         &      pic % at(3) .eq. 'm          3    .    .' .and. &
         &      pic % at(4) .eq. 'n          .   -1    4', &
         & "D3's coefficient picture carries [3, -1, 4]", nfail)

  end subroutine check_zero_is_not_absence

  !===================================================================!
  ! THE TWO-VIEW INVARIANT. For every cell of every direct
  ! dependency, the structural view and the coefficient view agree on
  ! PRESENCE - and presence is the relation's answer in both.
  !
  ! Positive, negative and zero coefficients all read as present. Only
  ! a missing tuple reads as absent.
  !===================================================================!

  subroutine check_the_two_views_agree_on_presence(nfail)

    integer, intent(inout) :: nfail

    call report(views_agree(d1, t1, h1, e1, w1) .and. &
         &      views_agree(d2, t2, h2, e2, w2) .and. &
         &      views_agree(d3, t3, h3, e3, w3), &
         & "on every cell of D1, D2 and D3 the two views agree on " // &
         & "presence - and BOTH take it from relation % has", nfail)

    call report(views_agree(d1, t1, h1, e1, w1_alt), &
         & "and they still agree when the coefficients are replaced " // &
         & "wholesale", nfail)

  end subroutine check_the_two_views_agree_on_presence

  !===================================================================!
  ! The Level-4 pictures, character for character, before and after
  ! the numbers existed - and against Gate A's own literals.
  !===================================================================!

  subroutine check_the_structural_pictures_are_unchanged(nfail)

    integer, intent(inout) :: nfail

    type(picture) :: after

    after = sparsity_picture(d1)
    call report(same_picture(before_d1, after) .and. &
         &      after % at(2) .eq. '        a b c d' .and. &
         &      after % at(3) .eq. 'p       # # . .' .and. &
         &      after % at(4) .eq. 'q       . # # .' .and. &
         &      after % at(5) .eq. 'r       . . . #', &
         & "D1's STRUCTURAL picture is the same string it was before " // &
         & "w1 existed, and the same one Gate A sealed", nfail)

    after = sparsity_picture(d2)
    call report(same_picture(before_d2, after) .and. &
         &      after % at(3) .eq. 'u       # # .', &
         & "D2's structural picture is unchanged", nfail)

    after = sparsity_picture(d3)
    call report(same_picture(before_d3, after) .and. &
         &      after % at(4) .eq. 'n       . # #', &
         & "D3's structural picture is unchanged", nfail)

  end subroutine check_the_structural_pictures_are_unchanged

  !===================================================================!
  ! And the algebra itself. Every relation Gate A derived is derived
  ! again with three fields alive in the program, and holds exactly
  ! the same tuples. VALUES HAVE NO PATH BACK INTO RELATION ALGEBRA.
  !===================================================================!

  subroutine check_the_derivation_is_unchanged(nfail)

    integer, intent(inout) :: nfail

    type(csr_relation), target :: again1, again2, again3, again21, again31
    type(csr_relation)         :: d31t

    again1 = derive_dependency('D1', t1, h1)
    again2 = derive_dependency('D2', t2, h2)
    again3 = derive_dependency('D3', t3, h3)
    again21 = derive_composition('D2:1', again1,  again2)
    again31 = derive_composition('D3:1', again21, again3)

    call report(same_extension(again1, d1) .and. same_extension(again2, d2) .and. &
         &      same_extension(again3, d3), &
         & "D1, D2, D3 re-derived with the fields alive: same " // &
         & "extensions", nfail)

    call report(again1 % num_tuples() .eq. ND1 .and. &
         &      again2 % num_tuples() .eq. ND2 .and. &
         &      again3 % num_tuples() .eq. ND3 .and. &
         &      again21 % num_tuples() .eq. ND21 .and. &
         &      again31 % num_tuples() .eq. ND31, &
         & "5, 4, 3, 6, 6 - and D2:1 still collapses its two " // &
         & "witnesses of b->u to one tuple", nfail)

    call report(same_extension(again21, d21) .and. same_extension(again31, d31), &
         & "D2:1 and D3:1 unchanged - and NO COEFFICIENT WAS " // &
         & "COMPOSED for either, which would be operator mathematics", &
         & nfail)

    ! The reverse stays value-free, exactly as Gate A left it.
    d31t = materialized_transpose('D3:1^T', again31)
    call report(d31t % num_tuples() .eq. ND31 .and. &
         &      .not. same_extension(again31, d31t), &
         & "and the structural reverse is untouched: no w^T was " // &
         & "built, and none was needed", nfail)

  end subroutine check_the_derivation_is_unchanged

  !===================================================================!
  ! THE INDEPENDENCE PROBE. A second field on the SAME E1, sharing no
  ! value with the first and containing no zero.
  !===================================================================!

  subroutine check_value_and_structure_are_independent(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: on
    type(picture)                  :: with_w1, with_alt, structural

    call w1_alt % domain(on)
    call report(on % same_as(e1) .and. coefficients_fit(w1_alt, e1), &
         & "w1_alt lives on the SAME E1 - the probe changes the " // &
         & "numbers and nothing else", nfail)

    structural = sparsity_picture(d1)
    call report(same_picture(structural, before_d1), &
         & "the structural picture of D1 is what it always was, with " // &
         & "either field in the program", nfail)

    with_w1  = valued_sparsity_picture(d1, t1, h1, e1, w1)
    with_alt = valued_sparsity_picture(d1, t1, h1, e1, w1_alt)

    call report(.not. same_picture(with_w1, with_alt), &
         & "AND THE COEFFICIENT PICTURES DIFFER - so D1 does not " // &
         & "determine w1", nfail)

    call report(with_alt % at(3) .eq. 'p         9   8   .   .' .and. &
         &      with_alt % at(4) .eq. 'q         .   7   6   .' .and. &
         &      with_alt % at(5) .eq. 'r         .   .   .   5', &
         & "the probe writes 7 where w1 wrote 0, and '.' in exactly " // &
         & "the same places - SO w1 DOES NOT DETERMINE D1", nfail)

    call report(near(value_at(w1, e1, E1_3), 0.0_dp) .and. &
         &      near(value_at(w1_alt, e1, E1_3), 7.0_dp) .and. &
         &      d1 % has([X0_B, X1_Q]), &
         & "e13 carries 0 in one field and 7 in the other, and b->q " // &
         & "is a dependency in both", nfail)

  end subroutine check_value_and_structure_are_independent

  !===================================================================!
  ! Helpers.
  !===================================================================!

  logical function near(got, want)

    real(dp), intent(in) :: got, want

    near = abs(got - want) .lt. TOL

  end function near

  logical function all_values_recovered(w, occurrences, want)

    class(field)     , intent(in) :: w
    class(member_set), intent(in) :: occurrences
    real(dp)         , intent(in) :: want(:)

    integer :: k

    all_values_recovered = (occurrences % size() .eq. size(want))
    do k = 1, occurrences % size()
       all_values_recovered = all_values_recovered .and. &
            & near(value_at(w, occurrences, occurrences % member(k)), want(k))
    end do

  end function all_values_recovered

  !-------------------------------------------------------------------!
  ! Exactly one occurrence per present tuple, none per absent pair -
  ! over every cell of the grid, not a sample.
  !-------------------------------------------------------------------!

  logical function seats_are_unique(d, tail, head, occurrences)

    class(relation)  , intent(in) :: d, tail, head
    class(member_set), intent(in) :: occurrences

    class(member_set), allocatable :: cols, rows
    integer                        :: i, j, want

    cols = d % domain(1)
    rows = d % domain(2)

    seats_are_unique = .true.
    do j = 1, cols % size()
       do i = 1, rows % size()
          if (d % has([cols % member(j), rows % member(i)])) then
             want = 1
          else
             want = 0
          end if
          seats_are_unique = seats_are_unique .and. &
               & (occurrences_joining(tail, head, occurrences, &
               &                      cols % member(j), rows % member(i)) &
               &  .eq. want)
       end do
    end do

  end function seats_are_unique

  !-------------------------------------------------------------------!
  ! THE TWO VIEWS, CELL BY CELL, with a third witness.
  !
  ! Both pictures are read by splitting each row into whitespace-
  ! separated tokens - the row label, then one token per column. No
  ! column position or field width is assumed, so the two layouts may
  ! differ freely and the comparison still holds.
  !
  !      in_structure   the token is '#'
  !      in_values      the token is not the absent glyph
  !      in_relation    the relation says so, asked directly
  !
  ! All three must agree, on every cell.
  !-------------------------------------------------------------------!

  logical function views_agree(d, tail, head, occurrences, w)

    class(relation)  , intent(in) :: d, tail, head
    class(member_set), intent(in) :: occurrences
    class(field)     , intent(in) :: w

    class(member_set), allocatable :: cols, rows
    type(picture)                  :: structural, valued
    character(len=32)              :: sbit(16), vbit(16)
    integer :: i, j, ns, nv
    logical :: in_structure, in_values, in_relation

    cols = d % domain(1)
    rows = d % domain(2)

    structural = sparsity_picture(d)
    valued     = valued_sparsity_picture(d, tail, head, occurrences, w)

    views_agree = .true.
    do i = 1, rows % size()

       call split(structural % at(2 + i), sbit, ns)
       call split(valued     % at(2 + i), vbit, nv)

       ! One row label plus one token per column, in both.
       views_agree = views_agree .and. (ns .eq. cols % size() + 1) &
            &                    .and. (nv .eq. cols % size() + 1)
       if (.not. views_agree) return

       do j = 1, cols % size()
          in_structure = (trim(sbit(j + 1)) .eq. '#')
          in_values    = (trim(vbit(j + 1)) .ne. ABSENT)
          in_relation  = d % has([cols % member(j), rows % member(i)])

          views_agree = views_agree .and. (in_structure .eqv. in_relation) &
               &                    .and. (in_values    .eqv. in_relation)
       end do
    end do

  end function views_agree

  !-------------------------------------------------------------------!
  ! Split a rendered line on blanks.
  !-------------------------------------------------------------------!

  subroutine split(line, bits, n)

    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: bits(:)
    integer         , intent(out) :: n

    integer :: k, start

    n     = 0
    start = 0
    do k = 1, len(line) + 1
       if (k .le. len(line)) then
          if (line(k:k) .ne. ' ') then
             if (start .eq. 0) start = k
             cycle
          end if
       end if
       if (start .ne. 0) then
          n = n + 1
          if (n .gt. size(bits)) error stop 'level 5: that line has more columns than the splitter expects'
          bits(n) = line(start : k - 1)
          start   = 0
       end if
    end do

  end subroutine split

  logical function same_picture(one, other)

    type(picture), intent(in) :: one, other

    integer :: k

    same_picture = (one % rows() .eq. other % rows())
    if (.not. same_picture) return
    do k = 1, one % rows()
       same_picture = same_picture .and. (one % at(k) .eq. other % at(k))
    end do

  end function same_picture

end program visualization_level_5
