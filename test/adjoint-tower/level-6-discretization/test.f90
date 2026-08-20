!=====================================================================!
! ADJOINT TOWER . LEVEL 6 . OPERATOR SUPPORT AND ORIENTATION
!
! The level answers one question: WHAT SHAPE WILL THE FOUR
! OPERATORS HAVE, AND WHICH WAY DO THEY POINT. The four derived
! blocks are read as the discrete supports of
!
!      Rq : Q -> Y        Rp : P -> Y
!      fq : Q -> Z        fp : P -> Z
!
! and the adjoint support is the TRANSPOSE VIEW of the same J_Q:
!
!      J_Q^T <= Q x Y = { (u,r1), (u,r2), (v,r1), (v,r2) }
!
! borrowed, not materialized, and never a second relation.
!
! Two laws are pinned here and they are the tower's spine.
!
! FIRST, orientation is about IDENTITY, not indices. Transposing a
! 2x2 array swaps two integers; transposing this support swaps the
! actual DOMAINS - J_Q answers on Y for Q, J_Q^T answers on Q for
! Y - and since |Q| = |Y| = 2, no dimension check could ever notice
! the difference. Only same_as can.
!
! SECOND, a support relation is written (row, column), so its first
! slot is the OPERATOR'S CODOMAIN: J_Q <= Y x Q supports Rq: Q -> Y.
! Reading a support as if it were the operator's direction is the
! easiest way to get an adjoint backwards, and it is exactly the
! mistake a non-symmetric A would then hide from a size check.
!
! Still no coefficient: this is where the operators will stand, not
! what they will multiply.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_6

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
  use relation_binary , only : csr_relation, transposed_view, &
       &                             transpose_of, inclusion_of

  implicit none

  type(set_graph)          :: v, t
  type(set_graph)           :: p_dom, q_dom, y_dom, z_dom
  type(stored_relation)      :: dep
  type(csr_relation), target :: inc_y, inc_z, inc_q, inc_p, jq
  type(transposed_view)      :: inc_q_t, inc_p_t, jq_t
  type(csr_relation)         :: jp, fq, fp
  integer                    :: table(2, 9)
  integer                    :: nfail
  type(set_map)     :: sets
  type(label_map)     :: labels
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 6 . operator support"
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

  inc_y   = inclusion_of(y_dom, t, sets, labels)
  inc_z   = inclusion_of(z_dom, t, sets, labels)
  inc_q   = inclusion_of(q_dom, v, sets, labels)
  inc_p   = inclusion_of(p_dom, v, sets, labels)
  inc_q_t = transpose_of(inc_q)
  inc_p_t = transpose_of(inc_p)

  jq = compose_binary(compose_binary(inc_y, dep, sets), inc_q_t, sets)
  jp = compose_binary(compose_binary(inc_y, dep, sets), inc_p_t, sets)
  fq = compose_binary(compose_binary(inc_z, dep, sets), inc_q_t, sets)
  fp = compose_binary(compose_binary(inc_z, dep, sets), inc_p_t, sets)

  ! The adjoint support: the SAME J_Q, read the other way.
  jq_t = transpose_of(jq)

  call check_operator_supports(nfail)
  call check_transpose_extension(nfail)
  call check_orientation_is_identity(nfail)
  call check_one_stored_truth(nfail)

  call verdict(nfail, "level 6")

contains

  !===================================================================!
  ! Where the four operators will stand. A support is written
  ! (row, column): the FIRST slot is the operator's codomain.
  !===================================================================!

  subroutine check_operator_supports(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: row, col

    row = jq % domain(1)
    col = jq % domain(2)
    call report(row % same_as(y_dom) .and. col % same_as(q_dom), &
         & "J_Q <= Y x Q supports Rq : Q -> Y", nfail)

    row = jp % domain(1)
    col = jp % domain(2)
    call report(row % same_as(y_dom) .and. col % same_as(p_dom), &
         & "J_P <= Y x P supports Rp : P -> Y", nfail)

    row = fq % domain(1)
    col = fq % domain(2)
    call report(row % same_as(z_dom) .and. col % same_as(q_dom), &
         & "F_Q <= Z x Q supports fq : Q -> Z", nfail)

    row = fp % domain(1)
    col = fp % domain(2)
    call report(row % same_as(z_dom) .and. col % same_as(p_dom), &
         & "F_P <= Z x P supports fp : P -> Z", nfail)

  end subroutine check_operator_supports

  !===================================================================!
  ! The adjoint support, exact - every pair of J_Q with its ends
  ! swapped, and nothing else.
  !===================================================================!

  subroutine check_transpose_extension(nfail)

    integer, intent(inout) :: nfail

    integer, allocatable :: tab(:,:)
    integer              :: k
    logical              :: ok

    call report(jq_t % num_tuples() .eq. 4 .and. &
         &      jq_t % has([VAR_U, TGT_R1]) .and. &
         &      jq_t % has([VAR_U, TGT_R2]) .and. &
         &      jq_t % has([VAR_V, TGT_R1]) .and. &
         &      jq_t % has([VAR_V, TGT_R2]), &
         & "J_Q^T = { (u,r1), (u,r2), (v,r1), (v,r2) }", nfail)

    call jq % tuples(tab)
    ok = .true.
    do k = 1, size(tab, 2)
       ok = ok .and. jq_t % has([tab(2, k), tab(1, k)])
    end do
    call report(ok, &
         & "and it holds every fact of J_Q, ends swapped: one truth, " // &
         & "read backwards", nfail)

  end subroutine check_transpose_extension

  !===================================================================!
  ! THE law of the rung: the transpose swaps DOMAIN IDENTITIES.
  ! With |Q| = |Y| = 2 no dimension check could tell J_Q from
  ! J_Q^T - and the specimen's A is non-symmetric, so a solver that
  ! confused them would answer a different, plausible-looking
  ! number. Only same_as separates them.
  !===================================================================!

  subroutine check_orientation_is_identity(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: row, col

    row = jq_t % domain(1)
    col = jq_t % domain(2)

    call report(row % same_as(q_dom) .and. col % same_as(y_dom), &
         & "J_Q^T <= Q x Y supports Rq^T : Y -> Q - the domains " // &
         & "themselves are exchanged", nfail)

    call report(sets % size_of(row) .eq. 2 .and. sets % size_of(col) .eq. 2 .and. &
         &      .not. row % same_as(col), &
         & "both slots hold two members and are still different " // &
         & "domains: no size check could catch a reversed adjoint", &
         & nfail)

    call report(.not. row % same_as(y_dom) .and. &
         &      .not. col % same_as(q_dom), &
         & "and neither slot is what the primal support had there", &
         & nfail)

  end subroutine check_orientation_is_identity

  !===================================================================!
  ! One stored truth: the primal support is materialized, the
  ! adjoint support borrows it. No second reverse relation exists
  ! to drift from the first.
  !===================================================================!

  subroutine check_one_stored_truth(nfail)

    integer, intent(inout) :: nfail

    call report(jq % materialized() .and. .not. jq_t % materialized(), &
         & "J_Q is stored; J_Q^T merely borrows it", nfail)

  end subroutine check_one_stored_truth

end program adjoint_level_6
