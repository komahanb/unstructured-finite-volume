!=====================================================================!
! VISUALIZATION TOWER . LEVEL 0 . SETS
!
! The level answers one question: WHAT STRUCTURAL SETS EXIST BEFORE
! DEPENDENCY EXISTS.
!
! Seven of them, and no more:
!
!      X0 = { a  b  c  d }         where the chain starts
!      X1 = { p  q  r }
!      X2 = { u  v  w }
!      X3 = { m  n }               where it ends
!
!      E1 = { e11 e12 e13 e14 e15 }
!      E2 = { e21 e22 e23 e24 }    the DEPENDENCY OCCURRENCES
!      E3 = { e31 e32 e33 }
!
! Level 0 knows nothing of relations, operators, graphs, fields,
! visualization, matrices, Jacobians or transpose. It knows that
! seven domains have been declared, how many members each holds, and
! that no two of them are the same domain.
!
!                  THE OCCURRENCE SET IS A CHOICE
!
! E1 is not scaffolding. Declaring the five dependency occurrences of
! A1 as MEMBERS OF A DOMAIN - rather than as pairs, or as nonzeros,
! or as anything else - gives each occurrence a first-class identity
! before any coefficient exists to hang on it. Level 1 will relate
! those occurrences to their two ends; Level 2 will forget them again
! on the way to D1. That E1 can be forgotten later is precisely why
! it has to be somewhere now.
!
!                      IDENTITY IS NOT CARDINALITY
!
! The specimen is built so that counting cannot be mistaken for
! knowing:
!
!      |X1| = |X2| = 3     and X1 is not X2
!      |X0| = |E2| = 4     and a state set is not an occurrence
!                          set
!
! and so that the integers cannot be mistaken for the members:
!
!      1 belongs to all seven sets, and means seven things.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program visualization_level_0

  use visualization_assert , only : report, verdict
  use visualization_assert , only : NX0, NX1, NX2, NX3, NE1, NE2, NE3
  use visualization_assert , only : X0_A, X0_B, X0_C, X0_D
  use visualization_assert , only : X1_P, X1_Q, X1_R
  use visualization_assert , only : X2_U, X2_V, X2_W
  use visualization_assert , only : X3_M, X3_N
  use visualization_assert , only : E1_1, E1_5, E2_1, E2_4, E3_1, E3_3
  use graph_set        , only : index_set
  use visualization_sets_fixture, only : structural_sets, label_for

  implicit none

  type(index_set) :: x0, x1, x2, x3, e1, e2, e3
  type(index_set) :: roll(7)
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 0 . sets"
  write(*,'(1x,a)') "============================================="

  call structural_sets(x0, x1, x2, x3, e1, e2, e3)
  roll = [x0, x1, x2, x3, e1, e2, e3]

  call check_the_seven_exist(nfail)
  call check_the_seven_are_distinct(nfail)
  call check_size_is_not_identity(nfail)
  call check_an_integer_is_not_a_member(nfail)
  call check_enumeration_laws(nfail)
  call check_the_reader_s_names(nfail)

  call verdict(nfail, "level 0")

contains

  !===================================================================!
  ! Seven declared domains, with the cardinalities the specimen
  ! names. The counts are compared against visualization_assert's
  ! constants, which the fixture never reads - two statements held
  ! against each other, not one read back to itself.
  !===================================================================!

  subroutine check_the_seven_exist(nfail)

    integer, intent(inout) :: nfail

    call report(x0 % size() .eq. NX0 .and. x1 % size() .eq. NX1 .and. &
         &      x2 % size() .eq. NX2 .and. x3 % size() .eq. NX3, &
         & "|X0|=4  |X1|=3  |X2|=3  |X3|=2 - the four state sets " // &
         & "of the operator chain", nfail)

    call report(e1 % size() .eq. NE1 .and. e2 % size() .eq. NE2 .and. &
         &      e3 % size() .eq. NE3, &
         & "|E1|=5  |E2|=4  |E3|=3 - the dependency OCCURRENCES, " // &
         & "first-class before any coefficient exists", nfail)

    call report(x0 % name() .eq. 'X0' .and. e3 % name() .eq. 'E3', &
         & "each set carries the name it was declared with - " // &
         & "metadata for the reader, and no part of the mathematics", nfail)

  end subroutine check_the_seven_exist

  !===================================================================!
  ! All twenty-one pairs, both ways round, and each set equals
  ! itself. Identity is declared once and answers for life.
  !===================================================================!

  subroutine check_the_seven_are_distinct(nfail)

    integer, intent(inout) :: nfail

    integer :: i, j
    logical :: self, apart

    self  = .true.
    apart = .true.
    do i = 1, 7
       self = self .and. roll(i) % equals(roll(i))
       do j = 1, 7
          if (i .eq. j) cycle
          apart = apart .and. (.not. roll(i) % equals(roll(j)))
       end do
    end do

    call report(self, &
         & "every set is itself - identity signed once, at " // &
         & "declaration", nfail)

    call report(apart, &
         & "and no two of the seven are the same domain: all 42 " // &
         & "ordered pairs stand apart", nfail)

  end subroutine check_the_seven_are_distinct

  !===================================================================!
  ! The two coincidences the specimen plants on purpose.
  !===================================================================!

  subroutine check_size_is_not_identity(nfail)

    integer, intent(inout) :: nfail

    call report(x1 % size() .eq. x2 % size() .and. &
         &      .not. x1 % equals(x2), &
         & "|X1| = |X2| = 3, AND X1 IS NOT X2 - counting the same " // &
         & "is not being the same", nfail)

    call report(x0 % size() .eq. e2 % size() .and. &
         &      .not. x0 % equals(e2), &
         & "|X0| = |E2| = 4, and a state set is still not an " // &
         & "occurrence set", nfail)

  end subroutine check_size_is_not_identity

  !===================================================================!
  ! One raw integer, seven meanings. a, p, u, m and e11 are all
  ! written 1, and the tower must never confuse them; that is what
  ! the named constants in visualization_assert are for.
  !===================================================================!

  subroutine check_an_integer_is_not_a_member(nfail)

    integer, intent(inout) :: nfail

    integer :: i
    logical :: everywhere

    everywhere = .true.
    do i = 1, 7
       everywhere = everywhere .and. roll(i) % has(1)
    end do

    call report(everywhere, &
         & "the raw integer 1 is a member of all seven sets", nfail)

    call report(X0_A .eq. 1 .and. X1_P .eq. 1 .and. X2_U .eq. 1 .and. &
         &      X3_M .eq. 1 .and. E1_1 .eq. 1 .and. E2_1 .eq. 1 .and. &
         &      E3_1 .eq. 1, &
         & "a, p, u, m, e11, e21 and e31 are all written 1 - SEVEN " // &
         & "OBJECTS WEARING ONE INTEGER", nfail)

    call report(.not. x0 % has(NX0 + 1) .and. .not. x3 % has(NX3 + 1) .and. &
         &      .not. e3 % has(NE3 + 1), &
         & "and each set refuses the member just past its last: " // &
         & "membership is the set's own answer", nfail)

    call report(x0 % has(X0_D) .and. x1 % has(X1_R) .and. &
         &      x2 % has(X2_W) .and. x3 % has(X3_N) .and. &
         &      e1 % has(E1_5) .and. e2 % has(E2_4) .and. e3 % has(E3_3), &
         & "every named member of the specimen is held by the set " // &
         & "that names it", nfail)

  end subroutine check_an_integer_is_not_a_member

  !===================================================================!
  ! The two enumeration laws, on every set and every member. The
  ! renderer three levels above walks sets by position and asks
  ! for members; these laws are what make that walk well-defined.
  !===================================================================!

  subroutine check_enumeration_laws(nfail)

    integer, intent(inout) :: nfail

    integer :: i, k
    logical :: forward, backward

    forward  = .true.
    backward = .true.
    do i = 1, 7
       do k = 1, roll(i) % size()
          forward  = forward .and. &
               &     (roll(i) % local_index(roll(i) % member(k)) .eq. k)
          backward = backward .and. &
               &     (roll(i) % member(roll(i) % local_index(k)) .eq. k)
       end do
    end do

    call report(forward, &
         & "local_index(member(i)) = i on all seven sets - " // &
         & "enumeration is injective", nfail)

    call report(backward, &
         & "member(local_index(a)) = a on all seven sets - and so " // &
         & "position and member determine each other", nfail)

    call report(x0 % local_index(NX0 + 1) .eq. 0, &
         & "an outsider stands nowhere: local_index answers zero " // &
         & "rather than guessing", nfail)

  end subroutine check_enumeration_laws

  !===================================================================!
  ! The reader's names, which no law reads. They are checked here so
  ! that Level 4 can print a picture without inventing anything - and
  ! checked to FALL BACK, so the renderer can be pointed at a set
  ! this fixture has never heard of.
  !===================================================================!

  subroutine check_the_reader_s_names(nfail)

    integer, intent(inout) :: nfail

    type(index_set) :: nameless

    call report(label_for(x0, X0_A) .eq. 'a' .and. &
         &      label_for(x0, X0_D) .eq. 'd' .and. &
         &      label_for(x1, X1_Q) .eq. 'q' .and. &
         &      label_for(x2, X2_W) .eq. 'w' .and. &
         &      label_for(x3, X3_N) .eq. 'n', &
         & "the state members are called a b c d / p q r / u v w / m n", &
         & nfail)

    call report(label_for(e1, E1_5) .eq. 'e15' .and. &
         &      label_for(e2, E2_4) .eq. 'e24' .and. &
         &      label_for(e3, E3_3) .eq. 'e33', &
         & "and the occurrences are called e11..e15 / e21..e24 / e31..e33", &
         & nfail)

    nameless = index_set('somewhere else', 3)
    call report(label_for(nameless, 2) .eq. '2' .and. &
         &      label_for(x0, NX0 + 1) .eq. '5', &
         & "a set the fixture does not know gets its member " // &
         & "printed as the integer it is - nothing is invented", nfail)

  end subroutine check_the_reader_s_names

end program visualization_level_0
