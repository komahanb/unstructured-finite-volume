!=====================================================================!
! MUTATION AND IDENTITY
!
! Evidence for the identity semantics, gathered against the shipped
! kernel. Nothing here modifies graph_fractal.f90.
!
!     T . the nine branch-state transitions
!     I . identity under lawful branch mutation
!     R . references to graphs without identity
!     X . structure that legitimately changes, in both directions
!     S . compiled representations as snapshots
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program mutation

  use token_identity, only : token
  use graph_fractal , only : graph, branch, &
       & BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN, &
       & null_branch, unknown_branch, known_branch
  use graph_views   , only : dp, attribute_map, csr, &
       & relation_view, residual, compile_csr

  implicit none

  integer :: failures = 0

  !===================================================================!
  ! T . All nine transitions of one branch, on one graph. Each is
  ! applied in place; after each, the branch invariant must hold and
  ! the graph's identity must be unchanged.
  !
  ! The kernel states three branch states. It states no transition
  ! law, and none of the nine is refused here.
  !===================================================================!

  transition_block: block

    type(graph), target :: g, ref
    integer             :: from, to
    logical             :: ok

    call ref % declare()

    ok = .true.
    do from = BRANCH_NULL, BRANCH_KNOWN
       do to = BRANCH_NULL, BRANCH_KNOWN
          ok = ok .and. transition(from, to, ref)
       end do
    end do
    call check('T  all nine transitions apply, invariant and identity intact', ok)

    ! And a transition may be applied repeatedly without accumulating
    ! state: the branch holds its current value, not a history.
    call g % declare()
    g % branch(1) = unknown_branch()
    g % branch(1) = known_branch(ref)
    g % branch(1) = unknown_branch()
    call check('T  KNOWN -> UNKNOWN is applied, so knowledge is not monotone', &
         & g % branch(1) % status() .eq. BRANCH_UNKNOWN .and. &
         & .not. associated(g % branch(1) % known()))

  end block transition_block

  !===================================================================!
  ! I . Identity is independent of branch state, before and after.
  !===================================================================!

  identity_block: block

    type(graph), target :: g, h
    type(token)         :: id0, id1

    call g % declare(); call h % declare()

    id0 = g % id()
    g % branch(1) = known_branch(h)
    g % branch(2) = unknown_branch()
    id1 = g % id()

    call check('I  id(g) is unchanged by lawful branch mutation', &
         & id0 % matches(id1))
    call check('I  g remains same_as itself across the mutation', g % same_as(g))
    call check('I  and distinct from h, whose branches are identical to none', &
         & .not. g % same_as(h))

  end block identity_block

  !===================================================================!
  ! R . A graph without identity is not same_as itself, so equality is
  ! not reflexive on it. A KNOWN branch to such a graph would place a
  ! non-reflexive element in the reachable set, and every map keyed on
  ! identity would fail to find it. known_branch refuses it instead.
  !===================================================================!

  reference_block: block

    type(graph), target :: declared, undeclared

    call declared % declare()

    call check('R  a declared graph is same_as itself', declared % same_as(declared))
    call check('R  an undeclared graph is not same_as itself', &
         & .not. undeclared % same_as(undeclared))
    call check('R  so identity cannot order a reference to it', &
         & .not. undeclared % same_as(declared) .and. &
         & .not. declared % same_as(undeclared))

  end block reference_block

  !===================================================================!
  ! X . Structure that legitimately changes, in both directions.
  !
  ! Three cells in a strip. Face f4 is a boundary: its far member is
  ! NULL. A fourth cell is attached, and f4 becomes interior; then it
  ! is detached, and f4 becomes a boundary again.
  !
  !     NULL -> KNOWN     attach
  !     KNOWN -> NULL     detach
  !
  ! f4 keeps one identity throughout. It is the same face, in three
  ! successive states. A monotone-knowledge law would forbid the
  ! second transition, and the second transition runs.
  !===================================================================!

  adaptive_block: block

    type(graph), target   :: c(4), f(4), node(4)
    type(attribute_map)   :: m
    type(token)           :: f4
    integer , allocatable :: left(:), right(:)
    real(dp), allocatable :: res(:)
    integer               :: i

    do i = 1, 4
       call c(i) % declare()
       call m % bind(c(i), index = i)
       call f(i) % declare()
       call node(i) % declare()
    end do

    f(1) % branch(1) = null_branch()      ; f(1) % branch(2) = known_branch(c(1))
    f(2) % branch(1) = known_branch(c(1)) ; f(2) % branch(2) = known_branch(c(2))
    f(3) % branch(1) = known_branch(c(2)) ; f(3) % branch(2) = known_branch(c(3))
    f(4) % branch(1) = known_branch(c(3)) ; f(4) % branch(2) = null_branch()

    do i = 1, 4
       node(i) % branch(1) = known_branch(f(i))
       if (i .lt. 4) then
          node(i) % branch(2) = known_branch(node(i + 1))
       else
          node(i) % branch(2) = null_branch()
       end if
    end do

    f4 = f(4) % id()
    call relation_view(node(1), m, left, right)
    call check('X  f4 is a boundary: (3,0)', left(4) .eq. 3 .and. right(4) .eq. 0)

    f(4) % branch(2) = known_branch(c(4))                    ! attach
    call relation_view(node(1), m, left, right)
    call check('X  after attach f4 is interior: (3,4)', &
         & left(4) .eq. 3 .and. right(4) .eq. 4)
    call check('X  and f4 is the same face: identity unchanged', &
         & f4 % matches(f(4) % id()))

    call residual(node(1), m, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], res)
    call check('X  the operator view follows the new structure: R(const) = 0', &
         & all(res .eq. 0.0_dp) .and. size(res) .eq. 4)

    f(4) % branch(2) = null_branch()                         ! detach
    call relation_view(node(1), m, left, right)
    call check('X  after detach f4 is a boundary again: (3,0)', &
         & left(4) .eq. 3 .and. right(4) .eq. 0)
    call check('X  identity still unchanged across both transitions', &
         & f4 % matches(f(4) % id()))

  end block adaptive_block

  !===================================================================!
  ! S . Compiled representations are snapshots.
  !
  !     G(t0) -> compile -> C0
  !     G mutates
  !     G(t1) -> compile -> C1
  !
  ! C0 stays internally valid and does not follow G. C1 reflects the
  ! new structure. Neither carries a graph version: each is the value
  ! that compilation returned.
  !===================================================================!

  snapshot_block: block

    type(graph), target :: c(4), f(4), node(4)
    type(attribute_map) :: m
    type(csr)           :: c0, c1
    integer             :: i

    do i = 1, 4
       call c(i) % declare()
       call m % bind(c(i), index = i)
       call f(i) % declare()
       call node(i) % declare()
    end do

    f(1) % branch(1) = null_branch()      ; f(1) % branch(2) = known_branch(c(1))
    f(2) % branch(1) = known_branch(c(1)) ; f(2) % branch(2) = known_branch(c(2))
    f(3) % branch(1) = known_branch(c(2)) ; f(3) % branch(2) = known_branch(c(3))
    f(4) % branch(1) = known_branch(c(3)) ; f(4) % branch(2) = null_branch()

    do i = 1, 4
       node(i) % branch(1) = known_branch(f(i))
       if (i .lt. 4) then
          node(i) % branch(2) = known_branch(node(i + 1))
       else
          node(i) % branch(2) = null_branch()
       end if
    end do

    c0 = compile_csr(node(1), m, 4)
    call check('S  C0 = [1,2,4,5,5] over four cells, the fourth isolated', &
         & all(c0 % rowptr .eq. [1, 2, 4, 5, 5]))

    f(4) % branch(2) = known_branch(c(4))
    c1 = compile_csr(node(1), m, 4)

    call check('S  C1 reflects the new structure', &
         & all(c1 % rowptr .eq. [1, 2, 4, 6, 7]) .and. &
         & any(c1 % colidx .eq. 4))
    call check('S  C0 is unchanged and internally valid', &
         & all(c0 % rowptr .eq. [1, 2, 4, 5, 5]) .and. &
         & size(c0 % colidx) .eq. c0 % rowptr(5) - 1)
    call check('S  C0 and C1 are independent values, neither carrying a version', &
         & size(c0 % colidx) .ne. size(c1 % colidx))

  end block snapshot_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'mutation: a proposition failed'
  end if

contains

  !===================================================================!
  ! One transition, applied in place: put the branch in state a, then
  ! in state b, and report whether the invariant and the identity both
  ! survive.
  !===================================================================!

  logical function transition(a, b, ref) result(good)

    integer            , intent(in) :: a, b
    type(graph), target, intent(in) :: ref

    type(graph), target :: t
    type(token)         :: before, after

    call t % declare()
    before = t % id()

    t % branch(1) = branch_with(a, ref)
    t % branch(1) = branch_with(b, ref)
    after = t % id()

    good = (t % branch(1) % status() .eq. b)                          .and. &
         & ((t % branch(1) % status() .eq. BRANCH_KNOWN)               .eqv. &
         &  associated(t % branch(1) % known()))                      .and. &
         & before % matches(after)

  end function transition

  type(branch) function branch_with(s, ref) result(b)

    integer            , intent(in) :: s
    type(graph), target, intent(in) :: ref

    select case (s)
    case (BRANCH_NULL)
       b = null_branch()
    case (BRANCH_UNKNOWN)
       b = unknown_branch()
    case default
       b = known_branch(ref)
    end select

  end function branch_with

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program mutation
