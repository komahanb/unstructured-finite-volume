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
! graph_carrier, and nothing above.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program adjoint_level_0

  use adjoint_assert, only : report, verdict
  use adjoint_assert, only : VAR_P, VAR_U, VAR_V
  use adjoint_assert, only : TGT_R1, TGT_R2, TGT_F
  use graph_carrier , only : counted_set, subset_set, member_set

  implicit none

  type(counted_set) :: v, t
  type(subset_set)  :: p_dom, q_dom, y_dom, z_dom
  integer           :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "adjoint tower . level 0 . domains"
  write(*,'(1x,a)') "============================================="

  ! The two parents, and the four roles carved from them.
  v = counted_set('variables', 3)
  t = counted_set('targets'  , 3)

  p_dom = subset_set('parameter', v, [VAR_P])
  q_dom = subset_set('state'    , v, [VAR_U, VAR_V])
  y_dom = subset_set('residual' , t, [TGT_R1, TGT_R2])
  z_dom = subset_set('response' , t, [TGT_F])

  call check_cardinalities(nfail)
  call check_enumeration_round_trips(nfail)
  call check_membership_boundary(nfail)
  call check_embeddings(nfail)
  call check_role_identities(nfail)

  call verdict(nfail, "level 0")

contains

  subroutine check_cardinalities(nfail)

    integer, intent(inout) :: nfail

    call report(v % size() .eq. 3 .and. t % size() .eq. 3, &
         & "V and T each count three members", nfail)
    call report(p_dom % size() .eq. 1 .and. q_dom % size() .eq. 2, &
         & "P holds the parameter, Q the two state slots", nfail)
    call report(y_dom % size() .eq. 2 .and. z_dom % size() .eq. 1, &
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

    call report(q_dom % local_index(VAR_V) .eq. 2 .and. &
         &      v % local_index(VAR_V) .eq. 3, &
         & "v is entry 2 of Q but member 3 of V: an id is not an " // &
         & "index", nfail)
    call report(z_dom % local_index(TGT_F) .eq. 1 .and. &
         &      t % local_index(TGT_F) .eq. 3, &
         & "and f is entry 1 of Z but member 3 of T", nfail)

  end subroutine check_enumeration_round_trips

  subroutine check_membership_boundary(nfail)

    integer, intent(inout) :: nfail

    call report(p_dom % has(VAR_P) .and. .not. p_dom % has(VAR_U), &
         & "P holds p alone", nfail)
    call report(q_dom % has(VAR_U) .and. q_dom % has(VAR_V) .and. &
         &      .not. q_dom % has(VAR_P), &
         & "Q holds u and v, and not the parameter", nfail)
    call report(y_dom % has(TGT_R1) .and. y_dom % has(TGT_R2) .and. &
         &      .not. y_dom % has(TGT_F), &
         & "Y holds the residual rows, and not the response", nfail)
    call report(z_dom % has(TGT_F) .and. .not. z_dom % has(TGT_R1), &
         & "Z holds the response alone", nfail)

    call report(.not. v % has(4) .and. .not. v % has(0) .and. &
         &      .not. t % has(4), &
         & "an outsider is rejected by both parents", nfail)

  end subroutine check_membership_boundary

  subroutine check_embeddings(nfail)

    integer, intent(inout) :: nfail

    call report(p_dom % is_subobject_of(v) .and. &
         &      q_dom % is_subobject_of(v), &
         & "P and Q stand embedded in the variables", nfail)
    call report(y_dom % is_subobject_of(t) .and. &
         &      z_dom % is_subobject_of(t), &
         & "Y and Z stand embedded in the targets", nfail)
    call report(.not. q_dom % is_subobject_of(t) .and. &
         &      .not. y_dom % is_subobject_of(v), &
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

    call report(q_dom % size() .eq. y_dom % size() .and. &
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

    class(member_set), intent(in) :: s

    integer :: i, m

    round_trips = .true.
    do i = 1, s % size()
       m = s % member(i)
       round_trips = round_trips .and. &
            & (s % member(s % local_index(m)) .eq. m)
       round_trips = round_trips .and. (s % local_index(m) .eq. i)
    end do

  end function round_trips

end program adjoint_level_0
