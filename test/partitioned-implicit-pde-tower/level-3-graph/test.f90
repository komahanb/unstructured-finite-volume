!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 3 . GRAPH
!
! The level answers one question: DOES THE ORDINARY GRAPH FAITHFULLY
! REALIZE THE STRUCTURE THE LOWER LEVELS DESCRIBED. This is the rung
! where
!
!      relation structure  ->  ordinary graph realization
!
! and nothing else happens here: no partition, no field, no
! operator, no solver.
!
! The comparison is EXTENSIONAL and signature-aware, not by carrier
! identity. G builds its own vertex and edge carriers, and there is
! no reason on earth they should be the same objects as the
! independent oracle's V and E - only that they hold the same
! members and relate them the same way. Demanding identity here
! would be demanding that two people who agree must be the same
! person.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_3

  use partitioned_pde_assert , only : report, verdict
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation
  use graph_set_map        , only : set_map
  use graph_binary_relation  , only : csr_relation
  use class_graph            , only : stored_graph
  use chain_carriers_fixture , only : chain_carriers
  use chain_relations_fixture, only : tail_relation, head_relation

  implicit none

  type(set_graph)  :: v, e, k
  type(set_map)  :: sets
  type(csr_relation) :: tail, head
  type(stored_graph) :: g
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 3 . graph"
  write(*,'(1x,a)') "============================================="

  call chain_carriers(sets, v, e, k)
  tail = tail_relation(e, v, sets)
  head = head_relation(e, v, sets)

  g = stored_graph(6, tails=[1,2,3,4,5], heads=[2,3,4,5,6])
  call sets % bind(g % vertex_set(), &
       & counted_set_representation(g % num_vertices()))
  call sets % bind(g % edge_set(), &
       & counted_set_representation(g % num_edges()))

  call check_counts(nfail)
  call check_carriers_extensionally(nfail)
  call check_incidence_against_the_oracle(nfail)
  call check_it_is_not_a_part(nfail)

  call verdict(nfail, "level 3")

contains

  subroutine check_counts(nfail)

    integer, intent(inout) :: nfail

    call report(g % num_vertices() .eq. sets % size_of(v) .and. &
         &      g % num_edges() .eq. sets % size_of(e), &
         & "G realizes six vertices and five edges", nfail)

  end subroutine check_counts

  !===================================================================!
  ! Extensional agreement: the same members, member for member. NOT
  ! carrier identity - G's carriers are its own.
  !===================================================================!

  subroutine check_carriers_extensionally(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: gv, ge
    integer           :: i
    logical           :: ok

    gv = g % vertex_set()
    ge = g % edge_set()

    ok = sets % size_of(gv) .eq. sets % size_of(v)
    do i = 1, sets % size_of(v)
       ok = ok .and. sets % has_in(gv, sets % member_of(v, i))
    end do
    call report(ok, &
         & "G's vertex carrier holds exactly V's members", nfail)

    ok = sets % size_of(ge) .eq. sets % size_of(e)
    do i = 1, sets % size_of(e)
       ok = ok .and. sets % has_in(ge, sets % member_of(e, i))
    end do
    call report(ok, &
         & "and its edge carrier exactly E's", nfail)

    call report(.not. gv % same_as(v), &
         & "yet G's carrier is NOT the oracle's V: agreement is " // &
         & "extensional, and two independent declarations are two " // &
         & "different domains", nfail)

  end subroutine check_carriers_extensionally

  !===================================================================!
  ! Every incidence G reports agrees with the Level-1 relations, in
  ! both directions of the check.
  !===================================================================!

  subroutine check_incidence_against_the_oracle(nfail)

    integer, intent(inout) :: nfail

    integer :: ge
    logical :: ok

    ok = .true.
    do ge = 1, g % num_edges()
       ok = ok .and. tail % has([ge, g % edge_tail(ge)])
       ok = ok .and. g % edge_has_head(ge)
       ok = ok .and. head % has([ge, g % edge_head(ge)])
    end do
    call report(ok, &
         & "every edge G reports leaves and enters where Tail and " // &
         & "Head say it should", nfail)

    ! ...and nothing the oracle forbids is present.
    ok = .true.
    do ge = 1, g % num_edges()
       ok = ok .and. .not. tail % has([ge, g % edge_head(ge)])
    end do
    call report(ok, &
         & "and no edge leaves the vertex it enters: the realization " // &
         & "adds nothing the relations did not license", nfail)

  end subroutine check_incidence_against_the_oracle

  subroutine check_it_is_not_a_part(nfail)

    integer, intent(inout) :: nfail

    call report(.not. g % has_part_relation(), &
         & "G is a whole graph, not a part of one", nfail)

  end subroutine check_it_is_not_a_part

end program partitioned_pde_level_3
