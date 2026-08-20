!=====================================================================!
! THE EPISTEMIC VIEW SUITE
!
! One graph ontology, read as the pair (Q, R):
!
!     branch(1) = Q     branch(2) = R
!
! The type under test is graph_fractal's graph. graph_state and its
! computational_graph are retired: the seats were the branches all
! along, the four states were four of nine, and the borrowed host was
! a privilege no kernel should hold.
!
! What the retired suite proved, and this one still proves:
!
!     realized is not solved
!     UNKNOWN is not empty, and UNKNOWN is not NULL
!     identity is independent of branch state
!
! What is new: NULL is a third primitive state, all nine combinations
! are instantiated, and the reading names only the four it defines.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test

  use token_identity     , only : token
  use graph_fractal      , only : graph, GRAPH_NULL, GRAPH_UNKNOWN, GRAPH_KNOWN, &
       & null_branch, unknown_branch, known_branch
  use view_epistemic, only : has_data, has_operator, &
       & epistemic_defined, epistemic_name, data_of, residual_of

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "epistemic view suite (AGENTS.md, The graph ontology)"

  !===================================================================!
  ! The four named readings, each built from branch values.
  !===================================================================!

  named_block: block

    type(graph), target  :: q, r
    type(graph), target  :: void, dat, op, realized
    type(graph), pointer :: pq, pr

    call q % declare(); call r % declare()
    call void % declare(); call dat % declare()
    call op % declare(); call realized % declare()

    void     % branch = [unknown_branch()  , unknown_branch()]
    dat      % branch = [known_branch(q)   , unknown_branch()]
    op       % branch = [unknown_branch()  , known_branch(r) ]
    realized % branch = [known_branch(q)   , known_branch(r) ]

    call check('void     = (UNKNOWN, UNKNOWN)', epistemic_name(void)     .eq. 'void')
    call check('data     = (KNOWN  , UNKNOWN)', epistemic_name(dat)      .eq. 'data')
    call check('operator = (UNKNOWN, KNOWN  )', epistemic_name(op)       .eq. 'operator')
    call check('realized = (KNOWN  , KNOWN  )', epistemic_name(realized) .eq. 'realized')

    call check('exactly one name holds: has_data and has_operator decide it', &
         & (.not. has_data(void))     .and. (.not. has_operator(void))     .and. &
         & (      has_data(dat))      .and. (.not. has_operator(dat))      .and. &
         & (.not. has_data(op))       .and. (      has_operator(op))       .and. &
         & (      has_data(realized)) .and. (      has_operator(realized)))

    pq => data_of(realized)
    pr => residual_of(realized)
    call check('Q and R are graphs, reached by reference', &
         & pq % same_as(q) .and. pr % same_as(r))

  end block named_block

  !===================================================================!
  ! REALIZED IS NOT SOLVED. Occupancy of both branches says nothing
  ! about R(Q) = 0. Two realized graphs, one consistent pair and one
  ! deliberately inconsistent, are indistinguishable to this reading.
  !===================================================================!

  realized_block: block

    type(graph), target :: q, r, consistent, inconsistent

    call q % declare(); call r % declare()
    call consistent % declare(); call inconsistent % declare()

    consistent   % branch = [known_branch(q), known_branch(r)]
    inconsistent % branch = [known_branch(r), known_branch(q)]   ! Q and R swapped

    call check('both are realized; the reading cannot tell them apart', &
         & epistemic_name(consistent)   .eq. 'realized' .and. &
         & epistemic_name(inconsistent) .eq. 'realized')
    call check('and they are different graphs', &
         & .not. consistent % same_as(inconsistent))

  end block realized_block

  !===================================================================!
  ! UNKNOWN IS NOT EMPTY. A Q that references a graph with no members
  ! is KNOWN; only an unrealized branch is UNKNOWN.
  !===================================================================!

  bottom_block: block

    type(graph), target :: empty_q, holds_empty, holds_nothing

    call empty_q % declare()                     ! (NULL, NULL): an atom
    call holds_empty % declare(); call holds_nothing % declare()

    holds_empty   % branch = [known_branch(empty_q), unknown_branch()]
    holds_nothing % branch = [unknown_branch()     , unknown_branch()]

    call check('a Q with no members is still KNOWN', has_data(holds_empty))
    call check('an unrealized Q is UNKNOWN, and that is a different absence', &
         & .not. has_data(holds_nothing))
    call check('the two readings differ: data /= void', &
         & epistemic_name(holds_empty) .ne. epistemic_name(holds_nothing))

  end block bottom_block

  !===================================================================!
  ! IDENTITY IS INDEPENDENT OF BRANCH STATE.
  !===================================================================!

  identity_block: block

    type(graph), target :: g, q, r
    type(token)         :: before, after

    call g % declare(); call q % declare(); call r % declare()

    g % branch = [unknown_branch(), unknown_branch()]
    before = g % id()
    g % branch = [known_branch(q), known_branch(r)]
    after = g % id()

    call check('void becomes realized under one identity', &
         & before % matches(after) .and. epistemic_name(g) .eq. 'realized')
    call check('two void graphs are still two graphs', &
         & .not. q % same_as(r))

  end block identity_block

  !===================================================================!
  ! ALL NINE COMBINATIONS. Four are named; the five involving NULL are
  ! outside the reading's domain and are reported as such, not forced
  ! into a name.
  !===================================================================!

  nine_block: block

    type(graph), target :: ref
    type(graph), target :: nn, nu, nk, un, uu, uk, kn, ku, kk

    call ref % declare()
    call nn % declare(); call nu % declare(); call nk % declare()
    call un % declare(); call uu % declare(); call uk % declare()
    call kn % declare(); call ku % declare(); call kk % declare()

    nn % branch = [null_branch()    , null_branch()    ]
    nu % branch = [null_branch()    , unknown_branch() ]
    nk % branch = [null_branch()    , known_branch(ref)]
    un % branch = [unknown_branch() , null_branch()    ]
    uu % branch = [unknown_branch() , unknown_branch() ]
    uk % branch = [unknown_branch() , known_branch(ref)]
    kn % branch = [known_branch(ref), null_branch()    ]
    ku % branch = [known_branch(ref), unknown_branch() ]
    kk % branch = [known_branch(ref), known_branch(ref)]

    call state('(N,N)', nn, GRAPH_NULL   , GRAPH_NULL   , .false.)
    call state('(N,U)', nu, GRAPH_NULL   , GRAPH_UNKNOWN, .false.)
    call state('(N,K)', nk, GRAPH_NULL   , GRAPH_KNOWN  , .false.)
    call state('(U,N)', un, GRAPH_UNKNOWN, GRAPH_NULL   , .false.)
    call state('(U,U)', uu, GRAPH_UNKNOWN, GRAPH_UNKNOWN, .true. )
    call state('(U,K)', uk, GRAPH_UNKNOWN, GRAPH_KNOWN  , .true. )
    call state('(K,N)', kn, GRAPH_KNOWN  , GRAPH_NULL   , .false.)
    call state('(K,U)', ku, GRAPH_KNOWN  , GRAPH_UNKNOWN, .true. )
    call state('(K,K)', kk, GRAPH_KNOWN  , GRAPH_KNOWN  , .true. )

    call check('four of nine combinations carry an epistemic name', &
         & count([epistemic_defined(nn), epistemic_defined(nu), epistemic_defined(nk), &
         &        epistemic_defined(un), epistemic_defined(uu), epistemic_defined(uk), &
         &        epistemic_defined(kn), epistemic_defined(ku), epistemic_defined(kk)]) .eq. 4)

  end block nine_block

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

  subroutine state(label, g, s1, s2, defined)

    character(len=*), intent(in) :: label
    type(graph)     , intent(in) :: g
    integer         , intent(in) :: s1, s2
    logical         , intent(in) :: defined

    logical :: ok

    ok = g % branch(1) % status() .eq. s1 .and. g % branch(2) % status() .eq. s2
    ok = ok .and. (epistemic_defined(g) .eqv. defined)
    call check(label, ok)

  end subroutine state

end program test
