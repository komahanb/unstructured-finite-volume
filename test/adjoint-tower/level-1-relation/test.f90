!=====================================================================!
! ADJOINT TOWER . LEVEL 1 . THE DEPENDENCY SOURCE
!
! The level answers one question: WHO MAY PARTICIPATE IN WHAT. One
! relation carries it,
!
!      R_dep  <=  T x V
!
! and for this specimen it is DENSE: each of r1, r2 and f may draw
! on each of p, u and v - nine facts, handed in ten times with one
! duplicate, because a relation is a set.
!
! That density is worth stating plainly. The derivative action
! tower's message came from SPARSITY - which slots could not reach
! the response. Here every target touches every variable, and the
! whole structural content of the tower is the ROLE PARTITION and
! the ORIENTATION built on it. A dense relation is not a degenerate
! one: it says every derivative block is structurally full, and the
! four blocks are told apart by their domains alone.
!
! And the relation knows nothing numerical. It does not know 2, 1,
! 3, 4, -4, -11, or that the response adds twice the parameter:
! those are law data, and they wait for constitution at Level 8.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_1

  use adjoint_assert, only : report, verdict
  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_fractal        , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use relation_finitary, only : stored_relation

  implicit none

  type(graph)     :: v, t
  type(stored_relation) :: dep
  integer               :: table(2, 10)
  integer               :: nfail
  type(set_map)     :: sets

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 1 . dependency"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call t % declare()
  call sets % bind(t, counted_set_representation(3))

  ! Nine facts - and the first handed twice.
  table(:,  1) = [TGT_R1, VAR_P]
  table(:,  2) = [TGT_R1, VAR_U]
  table(:,  3) = [TGT_R1, VAR_V]
  table(:,  4) = [TGT_R2, VAR_P]
  table(:,  5) = [TGT_R2, VAR_U]
  table(:,  6) = [TGT_R2, VAR_V]
  table(:,  7) = [TGT_F , VAR_P]
  table(:,  8) = [TGT_F , VAR_U]
  table(:,  9) = [TGT_F , VAR_V]
  table(:, 10) = [TGT_R1, VAR_P]

  dep = stored_relation('dependency', [t, v], table, sets)

  call check_signature(nfail)
  call check_complete_extension(nfail)
  call check_absences(nfail)

  call verdict(nfail, "level 1")

contains

  !===================================================================!
  ! Arity two, and the ordered signature answers the two DECLARED
  ! parents by identity - equal sizes buy nothing, and T and V both
  ! hold three members.
  !===================================================================!

  subroutine check_signature(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: d

    call report(dep % arity() .eq. 2, &
         & "the dependency is binary", nfail)

    d = dep % domain(1)
    call report(d % same_as(t) .and. .not. d % same_as(v), &
         & "slot one is the targets, by identity - not merely a " // &
         & "set of size three", nfail)
    d = dep % domain(2)
    call report(d % same_as(v) .and. .not. d % same_as(t), &
         & "slot two is the variables", nfail)

  end subroutine check_signature

  !===================================================================!
  ! The complete extension: nine held from ten handed, and every
  ! target-variable pair present. Nine present in a nine-element
  ! set: no tenth exists.
  !===================================================================!

  subroutine check_complete_extension(nfail)

    integer, intent(inout) :: nfail

    integer :: i, j
    logical :: ok

    call report(dep % num_tuples() .eq. 9, &
         & "ten handed, nine held: a relation is a set", nfail)

    ok = .true.
    do i = 1, sets % num_members_of(t)
       do j = 1, sets % num_members_of(v)
          ok = ok .and. dep % has([sets % member_of(t, i), sets % member_of(v, j)])
       end do
    end do
    call report(ok, &
         & "every target may draw on every variable: the specimen " // &
         & "is structurally dense", nfail)

    call report(dep % has([TGT_F, VAR_P]), &
         & "the response depends on the parameter DIRECTLY - the " // &
         & "fact that makes f_p live", nfail)
    call report(dep % has([TGT_R1, VAR_P]) .and. &
         &      dep % has([TGT_R2, VAR_P]), &
         & "and both residual rows depend on it too", nfail)

  end subroutine check_complete_extension

  !===================================================================!
  ! What a dense relation can and cannot say. It has NO in-range
  ! absence: every pair of members is a fact, so membership alone
  ! can never report a missing dependency here.
  !
  ! It can say even less than that. The two carriers enumerate their
  ! members from one, so their raw ids OVERLAP - TGT_R1 and VAR_P
  ! are both the integer 1 - and a tuple written in the wrong
  ! orientation, [VAR_U, TGT_R1] = [2, 1], is therefore a perfectly
  ! good member: it reads as (r2, p). Membership cannot probe
  ! orientation when the id spaces collide.
  !
  ! Orientation is carried by the SIGNATURE, checked by identity in
  ! check_signature above - slot one is T, slot two is V, and no
  ! integer comparison establishes that. This is the tower's lesson
  ! arriving early: identity distinguishes what values cannot.
  !===================================================================!

  subroutine check_absences(nfail)

    integer, intent(inout) :: nfail

    integer :: i, j
    logical :: any_absent

    any_absent = .false.
    do i = 1, sets % num_members_of(t)
       do j = 1, sets % num_members_of(v)
          any_absent = any_absent .or. &
               & .not. dep % has([sets % member_of(t, i), sets % member_of(v, j)])
       end do
    end do
    call report(.not. any_absent, &
         & "a dense relation has no in-range absence to report", &
         & nfail)

    call report(dep % has([VAR_U, TGT_R1]), &
         & "and ids collide across carriers, so even a reversed " // &
         & "pair reads as a member: orientation is the signature's " // &
         & "business, never membership's", nfail)

    call report(.not. dep % has([TGT_R1, VAR_U, VAR_V]), &
         & "a tuple of the wrong length belongs to nothing", nfail)

  end subroutine check_absences

end program adjoint_level_1
