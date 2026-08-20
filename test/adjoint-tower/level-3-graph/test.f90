!=====================================================================!
! ADJOINT TOWER . LEVEL 3 . RELATIONAL OWNERSHIP
!
! The level answers one question: CAN THE WHOLE STRUCTURAL MODEL
! STAND AS ONE OWNED STRUCTURE. It can:
!
!      GAMMA = ( { V, T, P, Q, Y, Z }, { R_dep, J_Q, J_P, F_Q, F_P } )
!
! six carriers - two parents and their four role subobjects, seated
! side by side as the ordinary member sets they are - and five
! relations, the source and the four blocks DERIVED once more by
! the certified Level-2 road and then admitted. The container
! stores structure; it derives nothing and infers nothing.
!
! The rung's real work is signature closure over a MIXED ownership:
! R_dep's slots resolve to the parents, while J_Q's resolve to the
! subobjects, and every one of them must be a carrier the graph
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
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_label      , only : label_map
  use map_inclusion  , only : inclusion_map, declared_subobject
  use relation_finitary, only : stored_relation, relation
  use relation_algebra, only : compose_binary
  use relation_binary , only : csr_relation, transposed_relation, &
       &                             transpose_of, inclusion_of
  use graph_fractal        , only : graph, known_branch, null_branch
  use view_relational, only : relational_binding, &
       & num_member_sets, member_set_at, num_relations, relation_at, &
       & holds_set

  implicit none

  type(set_graph)              :: v, t
  type(set_graph)               :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)          :: dep
  type(csr_relation), target     :: inc_y, inc_z, inc_q, inc_p
  type(transposed_relation)          :: inc_q_t, inc_p_t
  type(csr_relation)             :: jq, jp, fq, fp
  type(graph)             , target :: g
  type(graph)             , target :: scell(6), selem(6)
  type(graph)             , target :: rcell(5), relem(5)
  type(relational_binding)         :: bnd
  integer                          :: k
  integer                        :: table(2, 9)
  integer                        :: nfail
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 3 . ownership"
  write(*,'(1x,a)') "============================================="

  call v % declare()
  call sets % bind(v, counted_set_representation(3))
  call t % declare()
  call sets % bind(t, counted_set_representation(3))

  call p_dom % declare()
  call sets       % bind(p_dom, listed_set_representation([VAR_P]))
  call inclusions % include_in(p_dom, v)
  call q_dom % declare()
  call sets       % bind(q_dom, listed_set_representation([VAR_U, VAR_V]))
  call inclusions % include_in(q_dom, v)
  call y_dom % declare()
  call sets       % bind(y_dom, listed_set_representation([TGT_R1, TGT_R2]))
  call inclusions % include_in(y_dom, t)
  call z_dom % declare()
  call sets       % bind(z_dom, listed_set_representation([TGT_F]))
  call inclusions % include_in(z_dom, t)

  table(:, 1) = [TGT_R1, VAR_P]
  table(:, 2) = [TGT_R1, VAR_U]
  table(:, 3) = [TGT_R1, VAR_V]
  table(:, 4) = [TGT_R2, VAR_P]
  table(:, 5) = [TGT_R2, VAR_U]
  table(:, 6) = [TGT_R2, VAR_V]
  table(:, 7) = [TGT_F , VAR_P]
  table(:, 8) = [TGT_F , VAR_U]
  table(:, 9) = [TGT_F , VAR_V]
  dep = stored_relation('dependency', [t, v], table, sets)

  ! The certified Level-2 road, walked once more.
  inc_y = inclusion_of(y_dom, t, sets, labels)
  inc_z = inclusion_of(z_dom, t, sets, labels)
  inc_q = inclusion_of(q_dom, v, sets, labels)
  inc_p = inclusion_of(p_dom, v, sets, labels)
  inc_q_t = transpose_of(inc_q)
  inc_p_t = transpose_of(inc_p)

  jq = compose_binary(compose_binary(inc_y, dep, sets), inc_q_t, sets)
  jp = compose_binary(compose_binary(inc_y, dep, sets), inc_p_t, sets)
  fq = compose_binary(compose_binary(inc_z, dep, sets), inc_q_t, sets)
  fp = compose_binary(compose_binary(inc_z, dep, sets), inc_p_t, sets)

  ! 'adjoint specimen': (S, P) as one sequence on each branch.
  call g % declare()
  do k = 1, 6
     call scell(k) % declare()
     call selem(k) % declare()
  end do
  do k = 1, 5
     call rcell(k) % declare()
     call relem(k) % declare()
  end do

  call bnd % bind_set(selem(1), v)
  call bnd % bind_set(selem(2), t)
  call bnd % bind_set(selem(3), p_dom)
  call bnd % bind_set(selem(4), q_dom)
  call bnd % bind_set(selem(5), y_dom)
  call bnd % bind_set(selem(6), z_dom)
  call bnd % bind_relation(relem(1), dep)
  call bnd % bind_relation(relem(2), jq)
  call bnd % bind_relation(relem(3), jp)
  call bnd % bind_relation(relem(4), fq)
  call bnd % bind_relation(relem(5), fp)

  do k = 1, 6
     scell(k) % branch(1) = known_branch(selem(k))
     if (k .lt. 6) scell(k) % branch(2) = &
          & known_branch(scell(k + 1))
  end do
  do k = 1, 5
     rcell(k) % branch(1) = known_branch(relem(k))
     if (k .lt. 5) rcell(k) % branch(2) = &
          & known_branch(rcell(k + 1))
  end do

  g % branch(1) = known_branch(scell(1))
  g % branch(2) = known_branch(rcell(1))

  call check_ownership(nfail)
  call check_roles_are_citizens(nfail)
  call check_blocks_survive(nfail)
  call check_signature_closure(nfail)
  call check_graph_identity(nfail)

  call verdict(nfail, "level 3")

contains

  !===================================================================!
  ! Six carriers and five relations, owned by identity.
  !===================================================================!

  subroutine check_ownership(nfail)

    integer, intent(inout) :: nfail

    call report(num_member_sets(g) .eq. 6 .and. &
         &      num_relations(g) .eq. 5, &
         & "the graph owns six member sets and five relations", nfail)

    call report(holds_set(g, bnd, v) .and. holds_set(g, bnd, t), &
         & "the two parents are its own, by identity", nfail)

    call report(graph_holds_relation(g, bnd, dep), &
         & "the dependency source is owned", nfail)
    call report(graph_holds_relation(g, bnd, jq) .and. &
         &      graph_holds_relation(g, bnd, jp) .and. &
         &      graph_holds_relation(g, bnd, fq) .and. &
         &      graph_holds_relation(g, bnd, fp), &
         & "and so are all four derived blocks", nfail)

  end subroutine check_ownership

  !===================================================================!
  ! A subobject is an ordinary citizen: the four roles are seated
  ! beside their parents, each by its own identity.
  !===================================================================!

  subroutine check_roles_are_citizens(nfail)

    integer, intent(inout) :: nfail

    call report(holds_set(g, bnd, p_dom) .and. holds_set(g, bnd, q_dom) .and. &
         &      holds_set(g, bnd, y_dom) .and. holds_set(g, bnd, z_dom), &
         & "all four role subdomains are seated as carriers", nfail)

    call report(holds_set(g, bnd, q_dom) .and. holds_set(g, bnd, y_dom) .and. &
         &      .not. q_dom % same_as(y_dom), &
         & "Q and Y sit side by side, still not one another", nfail)

  end subroutine check_roles_are_citizens

  !===================================================================!
  ! The blocks survive ownership unchanged: J_Q keeps its four
  ! facts and its Y x Q orientation inside the container.
  !===================================================================!

  subroutine check_blocks_survive(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    type(set_graph) :: dom
    integer                        :: k

    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       if (rp % same_as(jq)) then
          call report(rp % num_tuples() .eq. 4 .and. &
               &      rp % has([TGT_R1, VAR_U]) .and. &
               &      rp % has([TGT_R2, VAR_V]), &
               & "the owned J_Q still holds its four derived facts", &
               & nfail)
          dom = rp % domain(1)
          call report(dom % same_as(y_dom), &
               & "its first slot is still the residual rows", nfail)
          dom = rp % domain(2)
          call report(dom % same_as(q_dom), &
               & "and its second still the state", nfail)
       end if
       if (rp % same_as(dep)) then
          call report(rp % num_tuples() .eq. 9, &
               & "and the dense source survives with all nine", nfail)
       end if
    end do

  end subroutine check_blocks_survive

  !===================================================================!
  ! Signature closure over MIXED ownership: R_dep resolves to the
  ! parents, the four blocks to the subobjects, and every slot of
  ! every owned relation lands on a carrier the graph holds.
  !===================================================================!

  subroutine check_signature_closure(nfail)

    integer, intent(inout) :: nfail

    class(relation), pointer       :: rp
    type(set_graph) :: dom
    integer                        :: k, s
    logical                        :: ok

    ok = .true.
    do k = 1, num_relations(g)
       rp => relation_at(g, bnd, k)
       do s = 1, rp % arity()
          dom = rp % domain(s)
          ok = ok .and. holds_set(g, bnd, dom)
       end do
    end do
    call report(ok, &
         & "every slot of every owned relation resolves to an owned " // &
         & "carrier - parents and subobjects alike", nfail)

  end subroutine check_signature_closure

  !===================================================================!
  ! G is a declared graph: itself, and not an identically built
  ! twin - extension equality never collapses identity.
  !===================================================================!

  subroutine check_graph_identity(nfail)

    integer, intent(inout) :: nfail

    type(graph)             , target :: g2
    type(graph)             , target :: scell2(2), selem2(2)
    type(relational_binding)         :: bnd2
    integer                          :: k2

    call report(g % same_as(g), &
         & "the specimen graph is itself", nfail)

    ! 'adjoint specimen again': (S, P) as one sequence on each branch.
    call g2 % declare()
    do k2 = 1, 2
       call scell2(k2) % declare()
       call selem2(k2) % declare()
    end do

    call bnd2 % bind_set(selem2(1), v)
    call bnd2 % bind_set(selem2(2), t)

    do k2 = 1, 2
       scell2(k2) % branch(1) = known_branch(selem2(k2))
       if (k2 .lt. 2) scell2(k2) % branch(2) = &
            & known_branch(scell2(k2 + 1))
    end do

    g2 % branch(1) = known_branch(scell2(1))
    g2 % branch(2) = null_branch()
    call report(.not. g % same_as(g2), &
         & "and no identically stocked twin is it", nfail)

  end subroutine check_graph_identity

  !===================================================================!
  ! Is this declared relation among the graph's own - composed from
  ! num_relations and relation_at, no convenience API.
  !===================================================================!

  logical function graph_holds_relation(g, b, r)

    type(graph)             , intent(in) :: g
    type(relational_binding), intent(in) :: b
    class(relation)         , intent(in) :: r

    class(relation), pointer :: rp
    integer                  :: k

    graph_holds_relation = .false.
    do k = 1, num_relations(g)
       rp => relation_at(g, b, k)
       if (rp % same_as(r)) graph_holds_relation = .true.
    end do

  end function graph_holds_relation

end program adjoint_level_3
