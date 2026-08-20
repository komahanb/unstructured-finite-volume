!=====================================================================!
! ADJOINT TOWER . LEVEL 0 . DOMAINS
!
! The level answers one question: WHAT ROLES EXIST, AND ARE THEY
! GENUINELY DIFFERENT. Two parent carriers are declared, and four
! roles are carved from them as subobjects:
!
!      V = { p  u  v }              T = { r1  r2  f }
!      P = { p }     c--> V         Y = { r1, r2 } c--> T
!      Q = { u, v }  c--> V         Z = { f }      c--> T
!
! and nothing else is true yet. The load-bearing truth of the rung
! is a NEGATIVE one:
!
!      |Q| = |Y| = 2      and yet      Q is NOT Y
!
! The state domain and the residual domain of an implicit problem
! have the same cardinality and no relationship whatever. Every
! square system in this repository will meet that temptation; this
! is where the tower refuses it, before any equation exists to be
! confused. Roles are DOMAINS - there is no parameter_field, no
! state_field, no residual_field, no adjoint_field - and the
! imports of this file ARE the negative truth: adjoint_assert,
! the set modules, and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_0

  use adjoint_assert, only : report, verdict
  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_fractal        , only : set_graph => graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set        , only : set_map
  use map_inclusion  , only : inclusion_map, declared_subobject

  implicit none

  type(set_graph) :: v, t
  type(set_graph)  :: p_dom, q_dom, y_dom, z_dom
  integer           :: nfail
  type(set_map)     :: sets
  type(inclusion_map)     :: inclusions

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 0 . domains"
  write(*,'(1x,a)') "============================================="

  ! The two parents, and the four roles carved from them.
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

  call check_cardinalities(nfail)
  call check_enumeration_round_trips(nfail)
  call check_membership_boundary(nfail)
  call check_embeddings(nfail)
  call check_role_identities(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(sets % size_of(v) .eq. 3 .and. sets % size_of(t) .eq. 3, &
         & "V and T each count three members", nfail)
    call report(sets % size_of(p_dom) .eq. 1 .and. sets % size_of(q_dom) .eq. 2, &
         & "P holds the parameter, Q the two state slots", nfail)
    call report(sets % size_of(y_dom) .eq. 2 .and. sets % size_of(z_dom) .eq. 1, &
         & "Y holds the two residual rows, Z the response", nfail)

  end subroutine check_cardinalities

  !===================================================================!
  ! The two enumeration laws on every member of every carrier and
  ! every role: member(local_index(m)) = m and
  ! local_index(member(i)) = i. Note the deliberate offsets - v is
  ! member 3 of V but entry 2 of Q - so a raw member id is never a
  ! storage position.
  !===================================================================!

  subroutine check_enumeration_round_trips(nfail)

    integer, intent(inout) :: nfail

    call report(round_trips(v) .and. round_trips(t), &
         & "member and local_index invert on V and T", nfail)
    call report(round_trips(p_dom) .and. round_trips(q_dom) .and. &
         &      round_trips(y_dom) .and. round_trips(z_dom), &
         & "and on every role subdomain", nfail)

    call report(sets % index_in(q_dom, VAR_V) .eq. 2 .and. &
         &      sets % index_in(v, VAR_V) .eq. 3, &
         & "v is entry 2 of Q but member 3 of V: an id is not an " // &
         & "index", nfail)
    call report(sets % index_in(z_dom, TGT_F) .eq. 1 .and. &
         &      sets % index_in(t, TGT_F) .eq. 3, &
         & "and f is entry 1 of Z but member 3 of T", nfail)

  end subroutine check_enumeration_round_trips

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(sets % has_in(p_dom, VAR_P) .and. .not. sets % has_in(p_dom, VAR_U), &
         & "P holds p alone", nfail)
    call report(sets % has_in(q_dom, VAR_U) .and. sets % has_in(q_dom, VAR_V) .and. &
         &      .not. sets % has_in(q_dom, VAR_P), &
         & "Q holds u and v, and not the parameter", nfail)
    call report(sets % has_in(y_dom, TGT_R1) .and. sets % has_in(y_dom, TGT_R2) .and. &
         &      .not. sets % has_in(y_dom, TGT_F), &
         & "Y holds the residual rows, and not the response", nfail)
    call report(sets % has_in(z_dom, TGT_F) .and. .not. sets % has_in(z_dom, TGT_R1), &
         & "Z holds the response alone", nfail)

    call report(.not. sets % has_in(v, 4) .and. .not. sets % has_in(v, 0) .and. &
         &      .not. sets % has_in(t, 4), &
         & "an outsider is rejected by both parents", nfail)

  end subroutine check_membership_boundary

  subroutine check_embeddings(nfail)

    integer, intent(inout) :: nfail

    call report(declared_subobject(p_dom, v, inclusions) .and. &
         &      declared_subobject(q_dom, v, inclusions), &
         & "P and Q stand embedded in the variables", nfail)
    call report(declared_subobject(y_dom, t, inclusions) .and. &
         &      declared_subobject(z_dom, t, inclusions), &
         & "Y and Z stand embedded in the targets", nfail)
    call report(.not. declared_subobject(q_dom, t, inclusions) .and. &
         &      .not. declared_subobject(y_dom, v, inclusions), &
         & "and neither is embedded in the other's parent", nfail)

  end subroutine check_embeddings

  !===================================================================!
  ! THE rung's truth: the roles are distinct by identity, and
  ! equal cardinality buys nothing. Q and Y are the pair that will
  ! tempt every square-system shortcut above this level.
  !===================================================================!

  subroutine check_role_identities(nfail)

    integer, intent(inout) :: nfail

    call report(.not. v % same_as(t), &
         & "V is not T", nfail)
    call report(.not. p_dom % same_as(q_dom), &
         & "the parameter domain is not the state domain", nfail)
    call report(.not. y_dom % same_as(z_dom), &
         & "the residual domain is not the response domain", nfail)

    call report(sets % size_of(q_dom) .eq. sets % size_of(y_dom) .and. &
         &      .not. q_dom % same_as(y_dom), &
         & "|Q| = |Y| = 2, and Q is STILL not Y: equal dimensions " // &
         & "are not equal domains", nfail)

    call report(q_dom % same_as(q_dom) .and. y_dom % same_as(y_dom), &
         & "and each role is itself", nfail)

  end subroutine check_role_identities

  !===================================================================!
  ! Both enumeration laws, on every member of any carrier - parents
  ! and subobjects alike, through the one abstract contract.
  !===================================================================!

  logical function round_trips(s)

    type(set_graph), intent(in) :: s

    integer :: i, m

    round_trips = .true.
    do i = 1, sets % size_of(s)
       m = sets % member_of(s, i)
       round_trips = round_trips .and. &
            & (sets % member_of(s, sets % index_in(s, m)) .eq. m)
       round_trips = round_trips .and. (sets % index_in(s, m) .eq. i)
    end do

  end function round_trips

end program adjoint_level_0
