!=====================================================================!
! THE TIME RELATIONS FIXTURE - earned at LEVEL 1.
!
! It holds the two primitive facts that give time its DIRECTION, and
! nothing else:
!
!      Tail <= E x T     e_i -> t_(i-1)     which instant a step leaves
!      Head <= E x T     e_i -> t_i         which instant it enters
!
! Direction is not the order of a loop index and not the sign of a
! subtraction. It is these two relations, and a step knows which end
! is which because the SIGNATURE says so, never because 2 is bigger
! than 1.
!
! The carriers Q, T and E are NOT declared here. They are Level 0's
! property, they live in time_carriers_fixture, and both
! constructors below receive them as arguments. This file cannot
! name a set into existence; it can only state facts over sets
! somebody else has already declared - which is what makes it a
! Level-1 file and not a Level-0 one.
!
! Q PARTICIPATES IN NOTHING HERE, and that absence is the level's
! sharpest content. Time has structure; state coordinates exist; and
! not one fact in this file says a state coordinate IS a time
! instant.
!
! A numbering hazard worth naming, because it will read as familiar
! and is not: over this specimen
!
!      Tail = { [1,1] [2,2] [3,3] [4,4] }      signature  E x T
!
! is tuple-for-tuple what a chain graph's tail map would be over a
! different pair of carriers entirely. The integers carry no
! meaning; the signature carries all of it.
!
! Nothing here knows about q(t), a history, a derivative, a scheme
! coefficient, Euler, BDF or a marcher.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_relations_fixture

  use graph_carrier        , only : counted_set
  use graph_binary_relation, only : csr_relation

  implicit none

  private
  public :: tail_relation, head_relation

contains

  !===================================================================!
  ! Tail <= E x T : the instant each step LEAVES.
  !
  !      e1 -> t0    e2 -> t1    e3 -> t2    e4 -> t3
  !
  ! and in the carriers' own numbering, where t0 is member one,
  ! that is the pair [i, i].
  !===================================================================!

  type(csr_relation) function tail_relation(e, t) result(tail)

    type(counted_set), intent(in) :: e, t

    integer :: table(2, 4), i

    do i = 1, 4
       table(:, i) = [i, i]
    end do
    tail = csr_relation('tail', e, t, table)

  end function tail_relation

  !===================================================================!
  ! Head <= E x T : the instant each step ENTERS.
  !
  !      e1 -> t1    e2 -> t2    e3 -> t3    e4 -> t4
  !
  ! the pair [i, i+1] in the carriers' numbering.
  !===================================================================!

  type(csr_relation) function head_relation(e, t) result(head)

    type(counted_set), intent(in) :: e, t

    integer :: table(2, 4), i

    do i = 1, 4
       table(:, i) = [i, i + 1]
    end do
    head = csr_relation('head', e, t, table)

  end function head_relation

end module time_relations_fixture
