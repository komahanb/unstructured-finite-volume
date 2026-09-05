!=====================================================================!
! KERNEL AND VIEW TESTS
!
! A-H are the required constructions, INV the invariant closure, I and
! J the evidence for the realization-model comparison. One graph type
! throughout; no subtype of graph is declared anywhere below.
!
! Navigation past one function reference is not a data-ref in Fortran,
! so a traversal deeper than one level binds a pointer or an ASSOCIATE
! name. Both forms appear below.
!
! Lifetime discipline: every graph is declared TARGET in the scope
! that owns it, and no procedure returns a pointer to a graph it
! constructed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test

  use graph_fractal, only : graph, branch, &
       & BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN, &
       & null_branch, unknown_branch, known_branch
  use graph_views  , only : dp, attribute_map, csr, &
       & evaluate, num_atoms, relation_view, residual, compile_csr

  implicit none

  integer :: failures = 0

  !===================================================================!
  ! A . Atoms. Identity is independent of branch state.
  !===================================================================!

  atom_block: block

    type(graph), target :: a, b

    call a % declare()
    call b % declare()

    call check('A  (NULL,NULL) holds by default initialization', &
         & a % branch(1) % status() .eq. BRANCH_NULL .and. &
         & a % branch(2) % status() .eq. BRANCH_NULL)
    call check('A  a same_as a', a % same_as(a))
    call check('A  a /= b although both are (NULL,NULL)', .not. a % same_as(b))
    call check('A  NULL implies .not. associated(known())', &
         & .not. associated(a % branch(1) % known()))

  end block atom_block

  !===================================================================!
  ! B . Sharing. Both branches of a reference one b.
  !===================================================================!

  sharing_block: block

    type(graph), target  :: a, b
    type(graph), pointer :: x, y

    call a % declare()
    call b % declare()

    a % branch(1) = known_branch(b)
    a % branch(2) = known_branch(b)

    call check('B  associated(a%branch(1)%known(), a%branch(2)%known())', &
         & associated(a % branch(1) % known(), a % branch(2) % known()))

    x => a % branch(1) % known()
    y => a % branch(2) % known()
    call check('B  a%branch(1)%known() same_as a%branch(2)%known() same_as b', &
         & x % same_as(b) .and. y % same_as(b))

    b % branch(1) = unknown_branch()
    call check('B  a change to b is observed on both branches: no copy', &
         & x % branch(1) % status() .eq. BRANCH_UNKNOWN .and. &
         & y % branch(1) % status() .eq. BRANCH_UNKNOWN)

  end block sharing_block

  !===================================================================!
  ! C . Cycle. a -> b -> a.
  !===================================================================!

  cycle_block: block

    type(graph), target  :: a, b
    type(graph), pointer :: x, y

    call a % declare()
    call b % declare()

    a % branch(1) = known_branch(b)
    b % branch(1) = known_branch(a)

    x => a % branch(1) % known()
    y => x % branch(1) % known()

    call check('C  a%branch(1)%known() same_as b', x % same_as(b))
    call check('C  b%branch(1)%known() same_as a', y % same_as(a))
    call check('C  two traversals return to a, by identity', y % same_as(a))
    call check('C  and by reference: closure without copy', associated(y, a))

  end block cycle_block

  !===================================================================!
  ! INV . Invariant closure.
  !
  !     status() == BRANCH_KNOWN  iff  associated(known())
  !
  ! Compile-time protection is asserted by run.sh against the negative
  ! fixtures in fortran-recursion/: no external scope can assign
  ! status_, rebind known_, or invoke the structure constructor.
  !===================================================================!

  invariant_block: block

    type(graph), target  :: ref, g(3)
    type(graph), pointer :: x
    integer              :: i, s
    logical              :: iff

    call ref % declare()
    do i = 1, 3
       call g(i) % declare()
    end do

    g(1) % branch = [null_branch()      , null_branch()      ]
    g(2) % branch = [unknown_branch()   , unknown_branch()   ]
    g(3) % branch = [known_branch(ref)  , known_branch(ref)  ]

    call check('INV KNOWN implies associated(known())', &
         & associated(g(3) % branch(1) % known()) .and. &
         & associated(g(3) % branch(2) % known()))
    call check('INV NULL implies .not. associated(known())', &
         & .not. associated(g(1) % branch(1) % known()) .and. &
         & .not. associated(g(1) % branch(2) % known()))
    call check('INV UNKNOWN implies .not. associated(known())', &
         & .not. associated(g(2) % branch(1) % known()) .and. &
         & .not. associated(g(2) % branch(2) % known()))

    iff = .true.
    do i = 1, 3
       do s = 1, 2
          iff = iff .and. ((g(i) % branch(s) % status() .eq. BRANCH_KNOWN) &
               &     .eqv. associated(g(i) % branch(s) % known()))
       end do
    end do
    call check('INV status()==KNOWN .eqv. associated(known()), every branch', iff)

    ! Navigation remains concise at depth through ASSOCIATE.
    associate (r => g(3) % branch(1) % known())
      call check('INV g%branch(1)%known() then %branch(2)%status() reads directly', &
           & r % same_as(ref) .and. r % branch(2) % status() .eq. BRANCH_NULL)
    end associate

    ! A whole branch value may be replaced; the iff survives it.
    g(3) % branch(1) = null_branch()
    call check('INV whole-branch replacement preserves the iff', &
         & g(3) % branch(1) % status() .eq. BRANCH_NULL .and. &
         & .not. associated(g(3) % branch(1) % known()))

  end block invariant_block

  !===================================================================!
  ! D . Expression view. (2 + 3) * 4 = 20, with 2, 3, 4, + and *
  ! bound in the attribute map and absent from the kernel.
  !===================================================================!

  expression_block: block

    type(graph), target :: two, three, four, plus, times
    type(attribute_map) :: m

    call two % declare(); call three % declare(); call four % declare()
    call plus % declare(); call times % declare()

    plus  % branch(1) = known_branch(two)
    plus  % branch(2) = known_branch(three)
    times % branch(1) = known_branch(plus)
    times % branch(2) = known_branch(four)

    call m % bind(two  , number = 2.0_dp)
    call m % bind(three, number = 3.0_dp)
    call m % bind(four , number = 4.0_dp)
    call m % bind(plus , symbol = '+')
    call m % bind(times, symbol = '*')

    call check('D  the operands are (NULL,NULL): structure does not encode 2 or 3', &
         & two % branch(1) % status() .eq. BRANCH_NULL .and. &
         & three % branch(1) % status() .eq. BRANCH_NULL)
    call check('D  evaluate(times) = 20', evaluate(times, m) .eq. 20.0_dp)
    call check('D  evaluate(plus) = 5', evaluate(plus, m) .eq. 5.0_dp)

  end block expression_block

  !===================================================================!
  ! E, F, G . Three views of one graph.
  !
  ! Three cells, four faces, the outer two with a NULL member:
  !
  !        |   c1   |   c2   |   c3   |
  !        f1       f2       f3       f4
  !
  ! The sequence of faces is constructed once, then read by the
  ! relation view (E), the operator view (F), and compile_csr, and it
  ! is finite-volume connectivity throughout (G).
  !===================================================================!

  connectivity_block: block

    type(graph), target   :: c(3), f(4), node(4)
    type(attribute_map)   :: m
    type(csr)             :: k
    integer , allocatable :: left(:), right(:)
    real(dp), allocatable :: res(:)
    integer               :: i
    logical               :: faithful

    do i = 1, 3
       call c(i) % declare()
       call m % bind(c(i), index = i)
    end do
    do i = 1, 4
       call f(i) % declare()
       call node(i) % declare()
    end do

    f(1) % branch(1) = null_branch()          ! boundary member
    f(1) % branch(2) = known_branch(c(1))
    f(2) % branch(1) = known_branch(c(1))
    f(2) % branch(2) = known_branch(c(2))
    f(3) % branch(1) = known_branch(c(2))
    f(3) % branch(2) = known_branch(c(3))
    f(4) % branch(1) = known_branch(c(3))
    f(4) % branch(2) = null_branch()          ! boundary member

    do i = 1, 4
       node(i) % branch(1) = known_branch(f(i))
       if (i .lt. 4) then
          node(i) % branch(2) = known_branch(node(i + 1))
       else
          node(i) % branch(2) = null_branch()
       end if
    end do

    ! E . relation view
    call relation_view(node(1), m, left, right)
    call check('E  |relation| = 4', size(left) .eq. 4)
    call check('E  interior tuples are (1,2) and (2,3)', &
         & left(2) .eq. 1 .and. right(2) .eq. 2 .and. &
         & left(3) .eq. 2 .and. right(3) .eq. 3)

    ! F . operator view
    call residual(node(1), m, [1.0_dp, 1.0_dp, 1.0_dp], res)
    call check('F  R(Q) = 0 for Q constant', all(res .eq. 0.0_dp))
    call residual(node(1), m, [1.0_dp, 2.0_dp, 4.0_dp], res)
    call check('F  R(Q) /= 0 otherwise, and sum(R) = 0', &
         & any(res .ne. 0.0_dp) .and. sum(res) .eq. 0.0_dp)

    ! G . connectivity, boundary by NULL
    call check('G  a boundary face has a NULL member', &
         & f(1) % branch(1) % status() .eq. BRANCH_NULL .and. &
         & f(4) % branch(2) % status() .eq. BRANCH_NULL)
    call check('G  only a NULL member maps to index 0', &
         & left(1) .eq. 0 .and. right(4) .eq. 0)

    ! Compiled representation of the same sequence.
    k = compile_csr(node(1), m, 3)
    call check('11 rowptr = [1,2,4,5]', all(k % rowptr .eq. [1, 2, 4, 5]))
    call check('11 colidx = [2,1,3,2]', all(k % colidx .eq. [2, 1, 3, 2]))

    faithful = .true.
    do i = 1, size(left)
       if (left(i) .eq. 0 .or. right(i) .eq. 0) cycle
       faithful = faithful .and. &
            & any(k % colidx(k % rowptr(left(i)):k % rowptr(left(i)+1)-1) .eq. right(i))
       faithful = faithful .and. &
            & any(k % colidx(k % rowptr(right(i)):k % rowptr(right(i)+1)-1) .eq. left(i))
    end do
    call check('11 every tuple of the relation view appears in the compiled rows', &
         & faithful)

  end block connectivity_block

  !===================================================================!
  ! H . The nine branch-state combinations, each instantiated.
  !===================================================================!

  states_block: block

    type(graph), target :: ref
    type(graph), target :: nn, nu, nk, un, uu, uk, kn, ku, kk

    call ref % declare()
    call nn % declare(); call nu % declare(); call nk % declare()
    call un % declare(); call uu % declare(); call uk % declare()
    call kn % declare(); call ku % declare(); call kk % declare()

    nn % branch(1) = null_branch()      ; nn % branch(2) = null_branch()
    nu % branch(1) = null_branch()      ; nu % branch(2) = unknown_branch()
    nk % branch(1) = null_branch()      ; nk % branch(2) = known_branch(ref)
    un % branch(1) = unknown_branch()   ; un % branch(2) = null_branch()
    uu % branch(1) = unknown_branch()   ; uu % branch(2) = unknown_branch()
    uk % branch(1) = unknown_branch()   ; uk % branch(2) = known_branch(ref)
    kn % branch(1) = known_branch(ref)  ; kn % branch(2) = null_branch()
    ku % branch(1) = known_branch(ref)  ; ku % branch(2) = unknown_branch()
    kk % branch(1) = known_branch(ref)  ; kk % branch(2) = known_branch(ref)

    call state('H  (N,N)', nn, BRANCH_NULL   , BRANCH_NULL   )
    call state('H  (N,U)', nu, BRANCH_NULL   , BRANCH_UNKNOWN)
    call state('H  (N,K)', nk, BRANCH_NULL   , BRANCH_KNOWN  )
    call state('H  (U,N)', un, BRANCH_UNKNOWN, BRANCH_NULL   )
    call state('H  (U,U)', uu, BRANCH_UNKNOWN, BRANCH_UNKNOWN)
    call state('H  (U,K)', uk, BRANCH_UNKNOWN, BRANCH_KNOWN  )
    call state('H  (K,N)', kn, BRANCH_KNOWN  , BRANCH_NULL   )
    call state('H  (K,U)', ku, BRANCH_KNOWN  , BRANCH_UNKNOWN)
    call state('H  (K,K)', kk, BRANCH_KNOWN  , BRANCH_KNOWN  )

    call check('H  NULL /= UNKNOWN by status, both disassociated', &
         & nu % branch(1) % status() .ne. nu % branch(2) % status() .and. &
         & .not. associated(nu % branch(1) % known()) .and. &
         & .not. associated(nu % branch(2) % known()))

  end block states_block

  !===================================================================!
  ! I . Mutable realization. UNKNOWN -> KNOWN after publication leaves
  ! identity fixed while the value of a view changes, so any result
  ! cached against identity becomes false.
  !===================================================================!

  mutable_block: block

    type(graph), target :: t, u, v
    integer             :: cached, recomputed

    call t % declare(); call u % declare(); call v % declare()

    t % branch(1) = known_branch(u)
    t % branch(2) = unknown_branch()

    cached = num_atoms(t)                   ! published, then evaluated
    t % branch(2) = known_branch(v)         ! UNKNOWN -> KNOWN, after
    recomputed = num_atoms(t)

    call check('I  num_atoms changed: 1 -> 2', &
         & cached .eq. 1 .and. recomputed .eq. 2)
    call check('I  identity did not change, so cached /= recomputed is undetectable', &
         & t % same_as(t) .and. cached .ne. recomputed)

  end block mutable_block

  !===================================================================!
  ! J . Persistent realization. A new identity per change leaves every
  ! predecessor referencing the previous graph, so the change must be
  ! propagated upward; around a cycle the propagation has no fixpoint.
  !===================================================================!

  persistent_block: block

    type(graph), target  :: a, b, c, b1, a1
    type(graph), target  :: p, q
    type(graph), pointer :: x, y

    call a % declare(); call b % declare(); call c % declare()

    b % branch(1) = unknown_branch()
    a % branch(1) = known_branch(b)

    call b1 % declare()                     ! b with the branch realized
    b1 % branch(1) = known_branch(c)

    x => a % branch(1) % known()
    call check('J  b1 /= b', .not. b1 % same_as(b))
    call check('J  a%branch(1)%known() same_as b, still UNKNOWN', &
         & x % same_as(b) .and. x % branch(1) % status() .eq. BRANCH_UNKNOWN)

    call a1 % declare()                     ! so a must be rebuilt as a1
    a1 % branch(1) = known_branch(b1)

    x => a1 % branch(1) % known()
    call check('J  a1 /= a and a1%branch(1)%known() same_as b1', &
         & .not. a1 % same_as(a) .and. x % same_as(b1))

    call p % declare(); call q % declare()  ! propagation around a cycle
    p % branch(1) = known_branch(q)
    q % branch(1) = known_branch(p)
    q % branch(2) = unknown_branch()

    x => q % branch(1) % known()
    y => x % branch(1) % known()
    call check('J  two traversals from q return to q: q1 requires q1', &
         & y % same_as(q))

  end block persistent_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'test: a proposition failed'
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

  subroutine state(label, g, s1, s2)

    character(len=*), intent(in) :: label
    type(graph)     , intent(in) :: g
    integer         , intent(in) :: s1, s2

    logical :: ok

    ok = g % branch(1) % status() .eq. s1 .and. g % branch(2) % status() .eq. s2
    ok = ok .and. (associated(g % branch(1) % known()) .eqv. (s1 .eq. BRANCH_KNOWN))
    ok = ok .and. (associated(g % branch(2) % known()) .eqv. (s2 .eq. BRANCH_KNOWN))
    call check(label, ok)

  end subroutine state

end program test
