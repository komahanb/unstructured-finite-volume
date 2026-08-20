!=====================================================================!
! VISUALIZATION TOWER . LEVEL 0 . CARRIERS
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
!                  THE OCCURRENCE CARRIER IS A CHOICE
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
!      |X0| = |E2| = 4     and a state carrier is not an occurrence
!                          carrier
!
! and so that the integers cannot be mistaken for the members:
!
!      1 belongs to all seven carriers, and means seven things.
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
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use visualization_carriers_fixture, only : structural_carriers, label_for

  implicit none

  type(graph) :: x0, x1, x2, x3, e1, e2, e3
  type(graph) :: roll(7)
  type(set_map)     :: sets
  type(label_map)     :: labels
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "visualization tower . level 0 . carriers"
  write(*,'(1x,a)') "============================================="

  call structural_carriers(x0, x1, x2, x3, e1, e2, e3, sets, labels)
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

    call report(sets % num_members_of(x0) .eq. NX0 .and. sets % num_members_of(x1) .eq. NX1 .and. &
         &      sets % num_members_of(x2) .eq. NX2 .and. sets % num_members_of(x3) .eq. NX3, &
         & "|X0|=4  |X1|=3  |X2|=3  |X3|=2 - the four state carriers " // &
         & "of the operator chain", nfail)

    call report(sets % num_members_of(e1) .eq. NE1 .and. sets % num_members_of(e2) .eq. NE2 .and. &
         &      sets % num_members_of(e3) .eq. NE3, &
         & "|E1|=5  |E2|=4  |E3|=3 - the dependency OCCURRENCES, " // &
         & "first-class before any coefficient exists", nfail)

    call report(labels % label_of(x0) .eq. 'X0' .and. labels % label_of(e3) .eq. 'E3', &
         & "each carrier carries the name it was declared with - " // &
         & "metadata for the reader, and no part of the mathematics", nfail)

  end subroutine check_the_seven_exist

  !===================================================================!
  ! All twenty-one pairs, both ways round, and each carrier same_as
  ! itself. Identity is declared once and answers for life.
  !===================================================================!

  subroutine check_the_seven_are_distinct(nfail)

    integer, intent(inout) :: nfail

    integer :: i, j
    logical :: self, apart

    self  = .true.
    apart = .true.
    do i = 1, 7
       self = self .and. roll(i) % same_as(roll(i))
       do j = 1, 7
          if (i .eq. j) cycle
          apart = apart .and. (.not. roll(i) % same_as(roll(j)))
       end do
    end do

    call report(self, &
         & "every carrier is itself - identity signed once, at " // &
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

    call report(sets % num_members_of(x1) .eq. sets % num_members_of(x2) .and. &
         &      .not. x1 % same_as(x2), &
         & "|X1| = |X2| = 3, AND X1 IS NOT X2 - counting the same " // &
         & "is not being the same", nfail)

    call report(sets % num_members_of(x0) .eq. sets % num_members_of(e2) .and. &
         &      .not. x0 % same_as(e2), &
         & "|X0| = |E2| = 4, and a state carrier is still not an " // &
         & "occurrence carrier", nfail)

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
       everywhere = everywhere .and. sets % has(roll(i), 1)
    end do

    call report(everywhere, &
         & "the raw integer 1 is a member of all seven carriers", nfail)

    call report(X0_A .eq. 1 .and. X1_P .eq. 1 .and. X2_U .eq. 1 .and. &
         &      X3_M .eq. 1 .and. E1_1 .eq. 1 .and. E2_1 .eq. 1 .and. &
         &      E3_1 .eq. 1, &
         & "a, p, u, m, e11, e21 and e31 are all written 1 - SEVEN " // &
         & "OBJECTS WEARING ONE INTEGER", nfail)

    call report(.not. sets % has(x0, NX0 + 1) .and. .not. sets % has(x3, NX3 + 1) .and. &
         &      .not. sets % has(e3, NE3 + 1), &
         & "and each carrier refuses the member just past its last: " // &
         & "membership is the carrier's own answer", nfail)

    call report(sets % has(x0, X0_D) .and. sets % has(x1, X1_R) .and. &
         &      sets % has(x2, X2_W) .and. sets % has(x3, X3_N) .and. &
         &      sets % has(e1, E1_5) .and. sets % has(e2, E2_4) .and. sets % has(e3, E3_3), &
         & "every named member of the specimen is held by the carrier " // &
         & "that names it", nfail)

  end subroutine check_an_integer_is_not_a_member

  !===================================================================!
  ! The two enumeration laws, on every carrier and every member. The
  ! renderer three levels above walks carriers by position and asks
  ! for members; these laws are what make that walk well-defined.
  !===================================================================!

  subroutine check_enumeration_laws(nfail)

    integer, intent(inout) :: nfail

    integer :: i, k
    logical :: forward, backward

    forward  = .true.
    backward = .true.
    do i = 1, 7
       do k = 1, sets % num_members_of(roll(i))
          forward  = forward .and. &
               &     (sets % index_in(roll(i), sets % member_of(roll(i), k)) .eq. k)
          backward = backward .and. &
               &     (sets % member_of(roll(i), sets % index_in(roll(i), k)) .eq. k)
       end do
    end do

    call report(forward, &
         & "local_index(member(i)) = i on all seven carriers - " // &
         & "enumeration is injective", nfail)

    call report(backward, &
         & "member(local_index(a)) = a on all seven carriers - and so " // &
         & "position and member determine each other", nfail)

    call report(sets % index_in(x0, NX0 + 1) .eq. 0, &
         & "an outsider stands nowhere: local_index answers zero " // &
         & "rather than guessing", nfail)

  end subroutine check_enumeration_laws

  !===================================================================!
  ! The reader's names, which no law reads. They are checked here so
  ! that Level 4 can print a picture without inventing anything - and
  ! checked to FALL BACK, so the renderer can be pointed at a carrier
  ! this fixture has never heard of.
  !===================================================================!

  subroutine check_the_reader_s_names(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: nameless

    call report(label_for(x0, X0_A, labels) .eq. 'a' .and. &
         &      label_for(x0, X0_D, labels) .eq. 'd' .and. &
         &      label_for(x1, X1_Q, labels) .eq. 'q' .and. &
         &      label_for(x2, X2_W, labels) .eq. 'w' .and. &
         &      label_for(x3, X3_N, labels) .eq. 'n', &
         & "the state members are called a b c d / p q r / u v w / m n", &
         & nfail)

    call report(label_for(e1, E1_5, labels) .eq. 'e15' .and. &
         &      label_for(e2, E2_4, labels) .eq. 'e24' .and. &
         &      label_for(e3, E3_3, labels) .eq. 'e33', &
         & "and the occurrences are called e11..e15 / e21..e24 / e31..e33", &
         & nfail)

    call nameless % declare()
    call sets % bind(nameless, counted_set_representation(3))
    call report(label_for(nameless, 2, labels) .eq. '2' .and. &
         &      label_for(x0, NX0 + 1, labels) .eq. '5', &
         & "a carrier the fixture does not know gets its member " // &
         & "printed as the integer it is - nothing is invented", nfail)

  end subroutine check_the_reader_s_names

end program visualization_level_0
