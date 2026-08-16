!=====================================================================!
! THE BINDING LIFETIME SUITE
!
! relational_binding lends pointers into the objects it owns, and
! graph_profile keeps such a pointer for the life of its view. This
! suite holds the storage law that makes that safe:
!
!     a borrowed pointer stays valid for the life of the binding,
!     across any number of later bindings.
!
! The law is a property of the STORAGE, not of graph. A graph mutates
! freely under stable identity; an object that lends pointers may
! impose stricter conditions than the graph does.
!
! WHAT WAS MEASURED AND REJECTED. With the object held in the row as an
! allocatable component, growth relocates the row array and a borrowed
! pointer then reads freed storage. The measurement, before the fix:
!
!     A before growth : same_as = T, name = R1
!     B after  growth : same_as = F              <- silently wrong
!       name through the same pointer            <- SIGSEGV
!
! Silently wrong first, fatal second. Every case below would have been
! a corrupted graph_profile had the cutover gone first.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program lifetime

  use fractal_graph        , only : graph, null_branch, known_branch
  use graph_carrier        , only : member_set, counted_set
  use graph_relation       , only : relation, stored_relation, slot
  use graph_binary_relation, only : binary_relation, csr_relation
  use graph_relational_view, only : relational_binding, relation_at, &
       & member_set_at

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "binding lifetime suite"

  !===================================================================!
  ! A . A pointer borrowed after all binding is complete.
  ! B . then bind_relation grows the relation storage.
  ! C . then bind_set grows the independent set storage.
  !===================================================================!

  growth_block: block

    type(graph), target      :: e1, e2, e3
    type(relational_binding) :: b
    type(counted_set)        :: s
    type(stored_relation)    :: r1, r2
    class(relation), pointer :: p, q

    call e1 % declare(); call e2 % declare(); call e3 % declare()
    s  = counted_set('E', 4)
    r1 = stored_relation('R1', [slot(s)], reshape([1, 2], [1, 2]))
    r2 = stored_relation('R2', [slot(s)], reshape([3, 4], [1, 2]))

    call b % bind_relation(e1, r1)
    p => b % relation_for(e1)

    call check('A  a borrowed pointer denotes the bound object', &
         & p % same_as(r1) .and. p % name() .eq. 'R1')

    call b % bind_relation(e2, r2)                  ! relation storage grows
    call check('B  it still denotes it after the relation storage grows', &
         & p % same_as(r1) .and. p % name() .eq. 'R1')
    q => b % relation_for(e1)
    call check('B  and a fresh borrow is the same storage', associated(p, q))

    call b % bind_set(e3, s)                        ! set storage grows
    call check('C  and after the independent set storage grows', &
         & p % same_as(r1) .and. p % name() .eq. 'R1')

  end block growth_block

  !===================================================================!
  ! D . The binding is intrinsically assigned to another binding. The
  ! copy owns copies, so the original's borrowers are untouched and
  ! nothing is freed twice.
  !===================================================================!

  assignment_block: block

    type(graph), target      :: e1
    type(relational_binding) :: b, c
    type(counted_set)        :: s
    type(stored_relation)    :: r1
    class(relation), pointer :: p, pc

    call e1 % declare()
    s  = counted_set('E', 4)
    r1 = stored_relation('R1', [slot(s)], reshape([1, 2], [1, 2]))

    call b % bind_relation(e1, r1)
    p => b % relation_for(e1)

    c = b
    call check('D  the pointer into the original survives the assignment', &
         & p % same_as(r1) .and. p % name() .eq. 'R1')

    pc => c % relation_for(e1)
    call check('D  the copy denotes an equal object in its own storage', &
         & pc % same_as(r1) .and. .not. associated(pc, p))

  end block assignment_block

  !===================================================================!
  ! E . The binding variable is overwritten by a second assignment.
  ! The overwritten binding's own storage is released, and the storage
  ! borrowed from the surviving one is not.
  !===================================================================!

  overwrite_block: block

    type(graph), target      :: e1, e2
    type(relational_binding) :: b, c, d
    type(counted_set)        :: s
    type(stored_relation)    :: r1, r2
    class(relation), pointer :: p

    call e1 % declare(); call e2 % declare()
    s  = counted_set('E', 4)
    r1 = stored_relation('R1', [slot(s)], reshape([1, 2], [1, 2]))
    r2 = stored_relation('R2', [slot(s)], reshape([3, 4], [1, 2]))

    call b % bind_relation(e1, r1)
    call d % bind_relation(e2, r2)

    c = b
    p => b % relation_for(e1)
    c = d                                            ! c is overwritten

    call check('E  the pointer into b is unaffected by overwriting c', &
         & p % same_as(r1) .and. p % name() .eq. 'R1')

  end block overwrite_block

  !===================================================================!
  ! F . The graph_profile pattern: borrow, narrow the dynamic type,
  ! keep the narrowed pointer, then grow the binding.
  !===================================================================!

  profile_pattern_block: block

    type(graph), target             :: g, rcell(2), relem(2), scell, selem
    type(relational_binding)        :: b
    type(counted_set)               :: e, v
    type(csr_relation)              :: t
    type(stored_relation)           :: extra
    class(relation), pointer        :: p
    class(binary_relation), pointer :: tails
    integer                         :: i

    call g % declare(); call scell % declare(); call selem % declare()
    do i = 1, 2
       call rcell(i) % declare(); call relem(i) % declare()
    end do

    e = counted_set('E', 2)
    v = counted_set('V', 3)
    t     = csr_relation('T', e, v, reshape([1, 1, 2, 2], [2, 2]))
    extra = stored_relation('X', [slot(e)], reshape([1, 2], [1, 2]))

    call b % bind_set(selem, e)
    call b % bind_relation(relem(1), t)

    scell % branch(1) = known_branch(selem); scell % branch(2) = null_branch()
    rcell(1) % branch(1) = known_branch(relem(1))
    rcell(1) % branch(2) = null_branch()
    g % branch(1) = known_branch(scell)
    g % branch(2) = known_branch(rcell(1))

    p => relation_at(g, b, 1)
    select type (p)
    class is (binary_relation)
       tails => p
    class default
       call check('F  the bound relation is binary', .false.)
    end select

    call check('F  the narrowed pointer answers before growth', &
         & tails % same_as(t) .and. tails % arity() .eq. 2)

    call b % bind_relation(relem(2), extra)          ! growth after narrowing
    call check('F  and after growth: this is the graph_profile pattern', &
         & tails % same_as(t) .and. tails % arity() .eq. 2)

  end block profile_pattern_block

  !===================================================================!
  ! G . TARGET. A row holds a pointer to the object, so the pointer
  ! returned does not point into the binding: it survives the return
  ! whether or not the actual argument has TARGET. Both are exercised;
  ! neither declares the binding TARGET.
  !===================================================================!

  target_block: block

    type(graph), target        :: g, scell, selem, rcell, relem
    type(relational_binding)   :: b                  ! no TARGET
    type(counted_set)          :: e
    type(stored_relation)      :: p
    class(relation), pointer   :: rp
    class(member_set), pointer :: sp

    call g % declare(); call scell % declare(); call selem % declare()
    call rcell % declare(); call relem % declare()

    e = counted_set('E', 2)
    p = stored_relation('P', [slot(e)], reshape([1, 2], [1, 2]))

    call b % bind_set(selem, e)
    call b % bind_relation(relem, p)

    scell % branch(1) = known_branch(selem); scell % branch(2) = null_branch()
    rcell % branch(1) = known_branch(relem); rcell % branch(2) = null_branch()
    g % branch(1) = known_branch(scell)
    g % branch(2) = known_branch(rcell)

    rp => relation_at(g, b, 1)                       ! through the view
    sp => member_set_at(g, b, 1)
    call check('G  relation_at survives return with no TARGET on the binding', &
         & rp % same_as(p))
    call check('G  member_set_at likewise', sp % same_as(e))

    rp => b % relation_for(relem)                    ! directly
    call check('G  relation_for likewise', rp % same_as(p))

  end block target_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'lifetime: a proposition failed'
  end if

contains

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

end program lifetime
