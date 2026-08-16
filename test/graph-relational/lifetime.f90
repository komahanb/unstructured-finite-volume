!=====================================================================!
! THE BINDING LIFETIME SUITE
!
! relational_binding lends pointers into the objects it owns, and
! graph_profile keeps such a pointer for the life of its view. This
! suite holds the storage law that makes that safe:
!
!     relational_binding is not assignable; bind_* preserves every
!     outstanding object pointer until the binding is destroyed.
!
! The law is a property of the STORAGE, not of graph. A graph mutates
! freely under stable identity; an object that lends pointers may
! impose stricter conditions than the graph does.
!
! Extension is here; refusal of replacement is in refusal.f90.
!
! WHAT WAS MEASURED AND REJECTED, twice.
!
! First, the object held in the row as an allocatable component. Growth
! relocates the row array, and a borrowed pointer then reads freed
! storage:
!
!     A before growth : same_as = T, name = R1
!     B after  growth : same_as = F              <- silently wrong
!       name through the same pointer            <- SIGSEGV
!
! Second, a deep copy on assignment. It finalizes its left-hand side
! first, so overwriting a lender freed what the lender had lent:
!
!     b = d, then read p borrowed from b
!       native   : associated, and answers R2    <- silently wrong
!       valgrind : invalid read of a freed block, and answers R1
!
! Two allocators, two different wrong answers, and the first is one no
! tool reports. No caller assigned a binding, so the copy that made
! this reachable is gone and assignment refuses.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program lifetime

  use fractal_graph        , only : graph, null_branch, known_branch
  use fractal_graph        , only : set_graph => graph
  use graph_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use graph_set_map        , only : set_map
  use graph_label_map      , only : label_map
  use graph_relation       , only : relation, stored_relation
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
    type(set_graph)        :: s
    type(stored_relation)    :: r1, r2
    class(relation), pointer :: p, q
    type(set_map)     :: sets

    call e1 % declare(); call e2 % declare(); call e3 % declare()
    call s % declare()
    call sets % bind(s, counted_set_representation(4))
    r1 = stored_relation('R1', [s], reshape([1, 2], [1, 2]), sets)
    r2 = stored_relation('R2', [s], reshape([3, 4], [1, 2]), sets)

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
  ! D . Two bindings coexist, each lending. Each owns its own objects,
  ! so growth in one is invisible to the other's borrowers. This is the
  ! whole of what replacement was ever wanted for.
  !===================================================================!

  independent_block: block

    type(graph), target      :: e1, e2, e3
    type(relational_binding) :: b, d
    type(set_graph)        :: s
    type(stored_relation)    :: r1, r2, r3
    class(relation), pointer :: p, q
    type(set_map)     :: sets

    call e1 % declare(); call e2 % declare(); call e3 % declare()
    call s % declare()
    call sets % bind(s, counted_set_representation(4))
    r1 = stored_relation('R1', [s], reshape([1, 2], [1, 2]), sets)
    r2 = stored_relation('R2', [s], reshape([3, 4], [1, 2]), sets)
    r3 = stored_relation('R3', [s], reshape([1, 4], [1, 2]), sets)

    call b % bind_relation(e1, r1)
    call d % bind_relation(e2, r2)
    p => b % relation_for(e1)
    q => d % relation_for(e2)

    call check('D  two bindings lend distinct storage', &
         & .not. associated(p, q) .and. p % name() .eq. 'R1' &
         & .and. q % name() .eq. 'R2')

    call d % bind_relation(e3, r3)                   ! only d grows
    call check('D  growth in one binding is invisible to the other', &
         & p % name() .eq. 'R1' .and. q % name() .eq. 'R2')

  end block independent_block

  !===================================================================!
  ! E . The graph_profile pattern: borrow, narrow the dynamic type,
  ! keep the narrowed pointer, then grow the binding.
  !===================================================================!

  profile_pattern_block: block

    type(graph), target             :: g, rcell(2), relem(2), scell, selem
    type(relational_binding)        :: b
    type(set_graph)               :: e, v
    type(csr_relation)              :: t
    type(stored_relation)           :: extra
    class(relation), pointer        :: p
    class(binary_relation), pointer :: tails
    integer                         :: i
    type(set_map)     :: sets

    call g % declare(); call scell % declare(); call selem % declare()
    do i = 1, 2
       call rcell(i) % declare(); call relem(i) % declare()
    end do

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    call v % declare()
    call sets % bind(v, counted_set_representation(3))
    t     = csr_relation('T', e, v, reshape([1, 1, 2, 2], [2, 2]), sets)
    extra = stored_relation('X', [e], reshape([1, 2], [1, 2]), sets)

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
       call check('E  the bound relation is binary', .false.)
    end select

    call check('E  the narrowed pointer answers before growth', &
         & tails % same_as(t) .and. tails % arity() .eq. 2)

    call b % bind_relation(relem(2), extra)          ! growth after narrowing
    call check('E  and after growth: this is the graph_profile pattern', &
         & tails % same_as(t) .and. tails % arity() .eq. 2)

  end block profile_pattern_block

  !===================================================================!
  ! F . TARGET. A row holds a pointer to the object, so the pointer
  ! returned does not point into the binding: it survives the return
  ! whether or not the actual argument has TARGET. Both are exercised;
  ! neither declares the binding TARGET.
  !===================================================================!

  target_block: block

    type(graph), target        :: g, scell, selem, rcell, relem
    type(relational_binding)   :: b                  ! no TARGET
    type(set_graph)          :: e
    type(stored_relation)      :: p
    class(relation), pointer   :: rp
    type(set_graph), pointer :: sp
    type(set_map)     :: sets

    call g % declare(); call scell % declare(); call selem % declare()
    call rcell % declare(); call relem % declare()

    call e % declare()
    call sets % bind(e, counted_set_representation(2))
    p = stored_relation('P', [e], reshape([1, 2], [1, 2]), sets)

    call b % bind_set(selem, e)
    call b % bind_relation(relem, p)

    scell % branch(1) = known_branch(selem); scell % branch(2) = null_branch()
    rcell % branch(1) = known_branch(relem); rcell % branch(2) = null_branch()
    g % branch(1) = known_branch(scell)
    g % branch(2) = known_branch(rcell)

    rp => relation_at(g, b, 1)                       ! through the view
    sp => member_set_at(g, b, 1)
    call check('F  relation_at survives return with no TARGET on the binding', &
         & rp % same_as(p))
    call check('F  member_set_at likewise', sp % same_as(e))

    rp => b % relation_for(relem)                    ! directly
    call check('F  relation_for likewise', rp % same_as(p))

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
