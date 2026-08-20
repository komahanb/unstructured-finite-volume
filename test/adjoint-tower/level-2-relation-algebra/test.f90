!=====================================================================!
! ADJOINT TOWER . LEVEL 2 . THE DERIVED SUPPORTS
!
! The level answers one question: WHERE DO THE FOUR DERIVATIVE
! BLOCKS COME FROM. Not from four hand-kept relations - from ONE
! dependency source and the roles' own relational faces. Each
! subobject S c--> A states its membership as an inclusion
! I_S <= S x A, and then the algebra does the rest:
!
!      J_Q = I_Q^T o R_dep o I_Y   <=  Y x Q      Y -> T -> V -> Q
!      J_P = I_P^T o R_dep o I_Y   <=  Y x P      Y -> T -> V -> P
!      F_Q = I_Q^T o R_dep o I_Z   <=  Z x Q      Z -> T -> V -> Q
!      F_P = I_P^T o R_dep o I_Z   <=  Z x P      Z -> T -> V -> P
!
! written as relational composition, read right to left. In code the
! same road runs left to right, because this repository reads
! compose_binary(P_AB, P_BC, sets) = P_BC o P_AB.
!
! The inclusion nearest R_dep on the right selects which targets the
! block answers ON; the transposed inclusion on the left - borrowed
! rather than rebuilt - selects which variables it answers FOR.
! Four blocks, one truth, no duplication.
!
! Their extensions are exact and their signatures are checked BY
! IDENTITY - which is the whole point, because J_Q <= Y x Q and
! F_Q <= Z x Q share a second slot while differing in the first,
! and J_Q and J_P share a first slot while differing in the second.
! Sizes could never tell them apart. Still no number anywhere: this
! is support, not derivative.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_2

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

  implicit none

  type(set_graph)          :: v, t
  type(set_graph)           :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)      :: dep
  type(csr_relation), target :: inc_y, inc_z, inc_q, inc_p
  type(transposed_relation)      :: inc_q_t, inc_p_t
  type(csr_relation)         :: on_y, on_z
  type(csr_relation)         :: jq, jp, fq, fp
  integer                    :: table(2, 9)
  integer                    :: nfail
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 2 . derived supports"
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

  ! The roles' own relational faces - membership, stated as algebra.
  inc_y = inclusion_of(y_dom, t, sets, labels)
  inc_z = inclusion_of(z_dom, t, sets, labels)
  inc_q = inclusion_of(q_dom, v, sets, labels)
  inc_p = inclusion_of(p_dom, v, sets, labels)

  ! Borrowed, never rebuilt: the right-hand selectors act reversed.
  inc_q_t = transpose_of(inc_q)
  inc_p_t = transpose_of(inc_p)

  ! Restrict the targets, then the variables.
  on_y = compose_binary(inc_y, dep, sets)      ! Y -> T -> V
  on_z = compose_binary(inc_z, dep, sets)      ! Z -> T -> V

  jq = compose_binary(on_y, inc_q_t, sets)     ! Y -> V -> Q
  jp = compose_binary(on_y, inc_p_t, sets)     ! Y -> V -> P
  fq = compose_binary(on_z, inc_q_t, sets)     ! Z -> V -> Q
  fp = compose_binary(on_z, inc_p_t, sets)     ! Z -> V -> P

  call check_inclusions(nfail)
  call check_restriction_to_targets(nfail)
  call check_state_block(nfail)
  call check_parameter_block(nfail)
  call check_response_blocks(nfail)
  call check_blocks_are_distinct(nfail)

  call verdict(nfail, "level 2")

contains

  !===================================================================!
  ! An inclusion says what membership already said, in the algebra's
  ! own language: (s, s) for every member, and nothing else.
  !===================================================================!

  subroutine check_inclusions(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = inc_q % domain(1)
    call report(dom % same_as(q_dom), &
         & "I_Q runs from the state domain", nfail)
    dom = inc_q % domain(2)
    call report(dom % same_as(v), &
         & "into the variables it is embedded in", nfail)
    call report(inc_q % num_tuples() .eq. 2 .and. &
         &      inc_q % has([VAR_U, VAR_U]) .and. &
         &      inc_q % has([VAR_V, VAR_V]) .and. &
         &      .not. inc_q % has([VAR_U, VAR_V]), &
         & "and holds each member to itself, nothing more", nfail)

    call report(inc_y % num_tuples() .eq. 2 .and. &
         &      inc_z % num_tuples() .eq. 1 .and. &
         &      inc_p % num_tuples() .eq. 1, &
         & "the other three inclusions count their own members", nfail)

    ! The transposed selector is a borrowed view of the same truth.
    call report(inc_q_t % num_tuples() .eq. 2 .and. &
         &      inc_q_t % has([VAR_U, VAR_U]) .and. &
         &      .not. inc_q_t % materialized(), &
         & "I_Q^T is a borrowed view, not a second inclusion", nfail)

  end subroutine check_inclusions

  !===================================================================!
  ! The first composition selects which targets a block answers on.
  !===================================================================!

  subroutine check_restriction_to_targets(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = on_y % domain(1)
    call report(dom % same_as(y_dom), &
         & "the residual restriction runs from Y", nfail)
    dom = on_y % domain(2)
    call report(dom % same_as(v), &
         & "into the whole variable carrier, still", nfail)
    call report(on_y % num_tuples() .eq. 6 .and. &
         &      .not. on_y % has([TGT_F, VAR_P]), &
         & "six facts survive: the response's are gone", nfail)

    call report(on_z % num_tuples() .eq. 3 .and. &
         &      on_z % has([TGT_F, VAR_P]) .and. &
         &      .not. on_z % has([TGT_R1, VAR_P]), &
         & "and the response restriction keeps exactly its three", &
         & nfail)

  end subroutine check_restriction_to_targets

  !===================================================================!
  ! J_Q <= Y x Q - the support of the state linearization Rq : Q -> Y.
  ! Note the orientation: the relation's FIRST slot is the operator's
  ! CODOMAIN.
  !===================================================================!

  subroutine check_state_block(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = jq % domain(1)
    call report(dom % same_as(y_dom), &
         & "J_Q runs from the residual rows", nfail)
    dom = jq % domain(2)
    call report(dom % same_as(q_dom), &
         & "into the state slots", nfail)

    call report(jq % num_tuples() .eq. 4 .and. &
         &      jq % has([TGT_R1, VAR_U]) .and. &
         &      jq % has([TGT_R1, VAR_V]) .and. &
         &      jq % has([TGT_R2, VAR_U]) .and. &
         &      jq % has([TGT_R2, VAR_V]), &
         & "J_Q = { (r1,u), (r1,v), (r2,u), (r2,v) } - exactly", nfail)
    call report(.not. jq % has([TGT_R1, VAR_P]), &
         & "and the parameter is not in it: that block is J_P", nfail)

  end subroutine check_state_block

  !===================================================================!
  ! J_P <= Y x P - the support of the parameter action Rp : P -> Y.
  !===================================================================!

  subroutine check_parameter_block(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = jp % domain(1)
    call report(dom % same_as(y_dom), &
         & "J_P runs from the residual rows", nfail)
    dom = jp % domain(2)
    call report(dom % same_as(p_dom), &
         & "into the parameter", nfail)

    call report(jp % num_tuples() .eq. 2 .and. &
         &      jp % has([TGT_R1, VAR_P]) .and. &
         &      jp % has([TGT_R2, VAR_P]), &
         & "J_P = { (r1,p), (r2,p) } - exactly", nfail)

  end subroutine check_parameter_block

  !===================================================================!
  ! F_Q <= Z x Q and F_P <= Z x P - the response's two dependencies,
  ! the indirect one through the state and the DIRECT one on the
  ! parameter. The second is what keeps f_p alive in the total
  ! derivative law.
  !===================================================================!

  subroutine check_response_blocks(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dom

    dom = fq % domain(1)
    call report(dom % same_as(z_dom), &
         & "F_Q runs from the response", nfail)
    dom = fq % domain(2)
    call report(dom % same_as(q_dom), &
         & "into the state slots", nfail)
    call report(fq % num_tuples() .eq. 2 .and. &
         &      fq % has([TGT_F, VAR_U]) .and. &
         &      fq % has([TGT_F, VAR_V]), &
         & "F_Q = { (f,u), (f,v) } - exactly", nfail)

    dom = fp % domain(2)
    call report(dom % same_as(p_dom), &
         & "F_P runs into the parameter", nfail)
    call report(fp % num_tuples() .eq. 1 .and. &
         &      fp % has([TGT_F, VAR_P]), &
         & "F_P = { (f,p) }: the response's DIRECT parameter " // &
         & "dependence stands on its own", nfail)

  end subroutine check_response_blocks

  !===================================================================!
  ! Four blocks, told apart by domain identity alone: J_Q and F_Q
  ! share their second slot, J_Q and J_P share their first, and no
  ! comparison of sizes could separate any of them.
  !===================================================================!

  subroutine check_blocks_are_distinct(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: a, b

    a = jq % domain(1)
    b = fq % domain(1)
    call report(.not. a % same_as(b), &
         & "J_Q and F_Q answer on different targets", nfail)

    a = jq % domain(2)
    b = fq % domain(2)
    call report(a % same_as(b), &
         & "though both answer for the state", nfail)

    a = jq % domain(2)
    b = jp % domain(2)
    call report(.not. a % same_as(b), &
         & "J_Q and J_P answer for different variables", nfail)

    a = jq % domain(1)
    b = jp % domain(1)
    call report(a % same_as(b), &
         & "though both answer on the residual rows", nfail)

  end subroutine check_blocks_are_distinct

end program adjoint_level_2
