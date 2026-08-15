!=====================================================================!
! ADJOINT TOWER . LEVEL 3 . RELATIONAL OWNERSHIP
!
! The level answers one question: CAN THE WHOLE STRUCTURAL MODEL
! STAND AS ONE OWNED STRUCTURE. It can:
!
!      G = ( { V, T, P, Q, Y, Z }, { R_dep, J_Q, J_P, F_Q, F_P } )
!
! six sets - two parents and their four role subobjects, seated
! side by side as the ordinary member sets they are - and five
! relations, the source and the four blocks DERIVED once more by
! the certified Level-2 road and then admitted. The container
! stores structure; it derives nothing and infers nothing.
!
! The rung's real work is signature closure over a MIXED ownership:
! R_dep's slots resolve to the parents, while J_Q's resolve to the
! subobjects, and every one of them must be a set the graph
! holds. A subobject is not a lesser citizen here - it is a domain
! with its own identity, and relations over it close exactly as
! relations over its parent do.
!
!      the complete adjoint problem has one structural ownership
!      environment before it has an adjoint solution
!
! There is no adjoint_graph. There is no second structure for the
! transpose. And still no number: nothing here knows what any of
! these blocks will one day multiply.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_3

  use adjoint_assert, only : report, verdict
  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_set , only : index_set, subset, set, unrelated_graph
  use graph_relation, only : stored_relation, relation
  use graph_relation_algebra, only : compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of
  use graph_structure, only : related_graph, declared_set, declared_relation

  implicit none

  type(index_set)              :: v, t
  type(subset)               :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_z, inc_q, inc_p
  type(transposed_view)          :: inc_q_t, inc_p_t
  type(csr_relation)             :: jq, jp, fq, fp
  type(related_graph), target :: g
  integer                        :: table(2, 9)
  integer                        :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 3 . ownership"
  write(*,'(1x,a)') "============================================="

  v = index_set('variables', 3)
  t = index_set('targets'  , 3)

  p_dom = subset('parameter', v, [VAR_P])
  q_dom = subset('state'    , v, [VAR_U, VAR_V])
  y_dom = subset('residual' , t, [TGT_R1, TGT_R2])
  z_dom = subset('response' , t, [TGT_F])

  table(:, 1) = [TGT_R1, VAR_P]
  table(:, 2) = [TGT_R1, VAR_U]
  table(:, 3) = [TGT_R1, VAR_V]
  table(:, 4) = [TGT_R2, VAR_P]
  table(:, 5) = [TGT_R2, VAR_U]
  table(:, 6) = [TGT_R2, VAR_V]
  table(:, 7) = [TGT_F , VAR_P]
  table(:, 8) = [TGT_F , VAR_U]
  table(:, 9) = [TGT_F , VAR_V]
  dep = stored_relation('dependency', [t, v], table)

  ! The certified Level-2 road, walked once more.
  inc_y = inclusion_of(y_dom)
  inc_z = inclusion_of(z_dom)
  inc_q = inclusion_of(q_dom)
  inc_p = inclusion_of(p_dom)
  inc_q_t = transpose_of(inc_q)
  inc_p_t = transpose_of(inc_p)

  jq = compose_binary(compose_binary(inc_y, dep), inc_q_t)
  jp = compose_binary(compose_binary(inc_y, dep), inc_p_t)
  fq = compose_binary(compose_binary(inc_z, dep), inc_q_t)
  fp = compose_binary(compose_binary(inc_z, dep), inc_p_t)

  g = related_graph('adjoint specimen', &
       & [declared_set(v), declared_set(t), declared_set(p_dom), declared_set(q_dom), &
       &  declared_set(y_dom), declared_set(z_dom)], &
       & [declared_relation(dep), declared_relation(jq), declared_relation(jp), &
       &  declared_relation(fq), declared_relation(fp)])

  call check_ownership(nfail)
  call check_roles_are_citizens(nfail)
  call check_blocks_survive(nfail)
  call check_signature_closure(nfail)
  call check_graph_identity(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Six sets and five relations, owned by identity.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(g % num_sets() .eq. 6 .and. &
         &      g % num_relations() .eq. 5, &
         & "the graph owns six member sets and five relations", nfail)

    call report(g % holds_set(v) .and. g % holds_set(t), &
         & "the two parents are its own, by identity", nfail)

    call report(graph_holds_relation(g, dep), &
         & "the dependency source is owned", nfail)
    call report(graph_holds_relation(g, jq) .and. &
         &      graph_holds_relation(g, jp) .and. &
         &      graph_holds_relation(g, fq) .and. &
         &      graph_holds_relation(g, fp), &
         & "and so are all four derived blocks", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! A subobject is an ordinary citizen: the four roles are seated
  ! beside their parents, each by its own identity.
  !===================================================================!

  subroutine check_roles_are_citizens(nfail)

    integer, intent(inout) :: nfail

    call report(g % holds_set(p_dom) .and. g % holds_set(q_dom) .and. &
         &      g % holds_set(y_dom) .and. g % holds_set(z_dom), &
         & "all four role subdomains are seated as sets", nfail)

    call report(g % holds_set(q_dom) .and. g % holds_set(y_dom) .and. &
         &      .not. q_dom % equals(y_dom), &
         & "Q and Y sit side by side, still not one another", nfail)

  end subroutine check_roles_are_citizens

  !===================================================================!
  ! The blocks survive ownership unchanged: J_Q keeps its four
  ! facts and its Y x Q orientation inside the container.
  !===================================================================!

  subroutine check_blocks_survive(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    class(set), allocatable :: dom
    integer                        :: k

    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % equals(jq)) then
          call report(rp % num_tuples() .eq. 4 .and. &
               &      rp % has([TGT_R1, VAR_U]) .and. &
               &      rp % has([TGT_R2, VAR_V]), &
               & "the owned J_Q still holds its four derived facts", &
               & nfail)
          dom = rp % domain(1)
          call report(dom % equals(y_dom), &
               & "its first slot is still the residual rows", nfail)
          dom = rp % domain(2)
          call report(dom % equals(q_dom), &
               & "and its second still the state", nfail)
       end if
       if (rp % equals(dep)) then
          call report(rp % num_tuples() .eq. 9, &
               & "and the dense source survives with all nine", nfail)
       end if
    end do

  end subroutine check_blocks_survive

  !===================================================================!
  ! Signature closure over MIXED ownership: R_dep resolves to the
  ! parents, the four blocks to the subobjects, and every slot of
  ! every owned relation lands on a set the graph holds.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    class(set), allocatable :: dom
    integer                        :: k, s
    logical                        :: ok

    ok = .true.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. g % holds_set(dom)
       end do
    end do
    call report(ok, &
         & "every slot of every owned relation resolves to an owned " // &
         & "set - parents and subobjects alike", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! G is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(unrelated_graph) :: g2

    call report(g % equals(g), &
         & "the specimen graph is itself", nfail)

    g2 = unrelated_graph('adjoint specimen again', &
         & [declared_set(v), declared_set(t)])
    call report(.not. g % equals(g2), &
         & "and no identically stocked twin is it", nfail)

  end subroutine check_graph_identity

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! num_relations and relation_at, no convenience API.
  !===================================================================!

  logical function graph_holds_relation(g, r)

    type(related_graph), target, intent(in) :: g
    class(relation)               , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, g % num_relations()
       rp => g % relation_at(k)
       if (rp % equals(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program adjoint_level_3
