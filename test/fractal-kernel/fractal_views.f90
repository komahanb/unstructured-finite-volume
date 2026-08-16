!=====================================================================!
! THE VIEWS . INTERPRETATION ABOVE THE KERNEL               (spike)
!
! The kernel holds generators; everything here is generated from
! them, and nothing here adds storage. One encoding law, stated
! once, that every reader below agrees on:
!
!    SEQUENCE   a cell's branch(1) is the element, KNOWN always;
!               branch(2) is the continuation - KNOWN the next
!               cell, NULL the end. The EMPTY sequence is a NULL
!               branch at whoever holds it: there is no empty
!               cell, so no second spelling of emptiness exists.
!               An UNKNOWN tail is a sequence still being learned
!               - it has no extent yet, and the walkers refuse it.
!
!    PAIR       branch(1) and branch(2), both KNOWN - except where
!               a reading permits absence, as a boundary face's
!               missing far side: NULL, and no imaginary member.
!
! One shape, many readings. The same list of pairs is read below
! as a binary relation, as directed adjacency, and as an ordinary
! neighbourhood; the same evaluator reads a calculator expression
! and a residual's right-hand side. Where two old concepts differ
! only by interpretation, they collapse onto one representation -
! that collapse is the acceptance criterion of this spike.
!
! Three laws the red team earned:
!
!    FINALITY    every status a walker's answer depends on is
!                final at the moment of answering - NULL and
!                KNOWN never change, and no walker answers
!                through UNKNOWN; it refuses. So no answer is
!                ever invalidated by a later realize.
!
!    OCCURRENCE  the pair list is an occurrence (edge-list)
!                reading: counts count occurrences, and rel_has
!                quotients them. Both are lawful; neither is the
!                other.
!
!    SHADOWING   an environment's newest binding wins - env_bind
!                conses in front, env_value answers the first
!                match. A solver that must not overwrite asks
!                env_knows first; a compressed numeral is never
!                a member.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module fractal_views

  use iso_fortran_env, only : dp => REAL64
  use fractal_graph  , only : graph, graph_arena, known_branch, &
       &                      null_branch, BRANCH_NULL, BRANCH_KNOWN, &
       &                      BRANCH_UNKNOWN

  implicit none

  private
  public :: is_atom, pair_of, seq_of, seq_len, seq_item
  public :: rel_has, rel_image_count
  public :: tail_of, head_of, out_degree
  public :: neighbour_count
  public :: env_of, env_bind, env_value, env_knows
  public :: evaluate
  public :: block_value
  public :: branch_name

contains

  !===================================================================!
  ! Shapes: the atom test, the pair, the sequence and its walkers.
  !===================================================================!

  logical function is_atom(g)

    type(graph), intent(in) :: g

    is_atom = g % status(1) .eq. BRANCH_NULL .and. &
         &    g % status(2) .eq. BRANCH_NULL

  end function is_atom

  type(graph) function pair_of(a, x, y) result(g)

    type(graph_arena), target, intent(inout) :: a
    type(graph)              , intent(in)    :: x, y

    g = a % node(known_branch(x), known_branch(y))

  end function pair_of

  type(graph) function seq_of(a, items) result(g)

    type(graph_arena), target, intent(inout) :: a
    type(graph)              , intent(in)    :: items(:)

    integer :: k

    if (size(items) .eq. 0) then
       error stop 'fractal view: the empty sequence is a null branch at its holder, never a cell'
    end if
    g = a % node(known_branch(items(size(items))), null_branch())
    do k = size(items) - 1, 1, -1
       g = a % node(known_branch(items(k)), known_branch(g))
    end do

  end function seq_of

  integer function seq_len(g) result(n)

    type(graph), intent(in) :: g

    type(graph) :: cell
    integer     :: steps

    n = 0
    cell = g
    do steps = 0, g % universe_size()
       if (cell % status(1) .ne. BRANCH_KNOWN) then
          error stop 'fractal view: a sequence cell holds its element as knowledge'
       end if
       n = n + 1
       select case (cell % status(2))
       case (BRANCH_NULL)
          return
       case (BRANCH_UNKNOWN)
          error stop 'fractal view: an unknown tail has no extent yet'
       case (BRANCH_KNOWN)
          cell = cell % known(2)
       end select
    end do
    error stop 'fractal view: this spine never ends - a cycle has no extent'

  end function seq_len

  type(graph) function seq_item(g, k) result(e)

    type(graph), intent(in) :: g
    integer    , intent(in) :: k

    type(graph) :: cell
    integer     :: step

    if (k .lt. 1) error stop 'fractal view: a sequence counts from one'
    cell = g
    do step = 2, k
       select case (cell % status(2))
       case (BRANCH_NULL)
          error stop 'fractal view: this sequence ends before that item'
       case (BRANCH_UNKNOWN)
          error stop 'fractal view: an unknown tail has no extent yet'
       end select
       cell = cell % known(2)
    end do
    e = cell % known(1)

  end function seq_item

  !===================================================================!
  ! Three readings of one list of pairs. The relation asks whether
  ! members stand together; the directed reading asks tail and
  ! head; the ordinary reading forgets order. A pair whose second
  ! branch is NULL is a boundary occurrence: it relates nobody,
  ! heads nowhere, and neighbours nothing - no imaginary member.
  !===================================================================!

  logical function rel_has(rel, x, y)

    type(graph), intent(in) :: rel, x, y

    type(graph) :: p, s, t
    integer     :: k

    rel_has = .false.
    do k = 1, seq_len(rel)
       p = seq_item(rel, k)
       s = p % known(1)
       if (.not. s % same_as(x)) cycle
       select case (p % status(2))
       case (BRANCH_NULL)
          cycle                          ! a boundary occurrence relates nobody
       case (BRANCH_UNKNOWN)
          error stop 'fractal view: an unlearned far side has no answer yet'
       end select
       t = p % known(2)
       if (t % same_as(y)) then
          rel_has = .true.
          return
       end if
    end do

  end function rel_has

  integer function rel_image_count(rel, x) result(n)

    type(graph), intent(in) :: rel, x

    type(graph) :: p, s
    integer     :: k

    n = 0
    do k = 1, seq_len(rel)
       p = seq_item(rel, k)
       s = p % known(1)
       if (.not. s % same_as(x)) cycle
       select case (p % status(2))
       case (BRANCH_NULL)
          cycle
       case (BRANCH_UNKNOWN)
          error stop 'fractal view: an unlearned far side has no answer yet'
       end select
       n = n + 1
    end do

  end function rel_image_count

  type(graph) function tail_of(rel, k) result(v)

    type(graph), intent(in) :: rel
    integer    , intent(in) :: k

    type(graph) :: p

    p = seq_item(rel, k)
    v = p % known(1)

  end function tail_of

  type(graph) function head_of(rel, k) result(v)

    type(graph), intent(in) :: rel
    integer    , intent(in) :: k

    type(graph) :: p

    p = seq_item(rel, k)
    select case (p % status(2))
    case (BRANCH_NULL)
       error stop 'fractal view: a boundary occurrence heads nowhere'
    case (BRANCH_UNKNOWN)
       error stop 'fractal view: an unlearned far side has no answer yet'
    end select
    v = p % known(2)

  end function head_of

  integer function out_degree(rel, v) result(n)

    type(graph), intent(in) :: rel, v

    n = rel_image_count(rel, v)

  end function out_degree

  integer function neighbour_count(rel, v) result(n)

    type(graph), intent(in) :: rel, v

    type(graph) :: p, s, t
    integer     :: k

    n = 0
    do k = 1, seq_len(rel)
       p = seq_item(rel, k)
       select case (p % status(2))
       case (BRANCH_NULL)
          cycle
       case (BRANCH_UNKNOWN)
          error stop 'fractal view: an unlearned far side has no answer yet'
       end select
       s = p % known(1)
       t = p % known(2)
       if (s % same_as(v) .or. t % same_as(v)) n = n + 1
    end do

  end function neighbour_count

  !===================================================================!
  ! An environment is itself a graph: a sequence of (member, value
  ! atom) pairs. Binding is consing; lookup is identity search.
  !===================================================================!

  type(graph) function env_of(a, member, val) result(env)

    type(graph_arena), target, intent(inout) :: a
    type(graph)              , intent(in)    :: member
    real(dp)                 , intent(in)    :: val

    type(graph) :: carrier, binding

    if (member % carries_value()) then
       error stop 'fractal view: a compressed numeral is never a member'
    end if
    carrier = a % value_atom(val)
    binding = pair_of(a, member, carrier)
    env = a % node(known_branch(binding), null_branch())

  end function env_of

  type(graph) function env_bind(a, env, member, val) result(bound)

    type(graph_arena), target, intent(inout) :: a
    type(graph)              , intent(in)    :: env, member
    real(dp)                 , intent(in)    :: val

    type(graph) :: carrier, binding

    if (member % carries_value()) then
       error stop 'fractal view: a compressed numeral is never a member'
    end if
    carrier = a % value_atom(val)
    binding = pair_of(a, member, carrier)
    bound = a % node(known_branch(binding), known_branch(env))

  end function env_bind

  real(dp) function env_value(env, member) result(v)

    type(graph), intent(in) :: env, member

    type(graph) :: p, s, t
    integer     :: k

    do k = 1, seq_len(env)
       p = seq_item(env, k)
       s = p % known(1)
       if (s % same_as(member)) then
          t = p % known(2)
          v = t % value()
          return
       end if
    end do
    error stop 'fractal view: this environment binds no value to that member'

  end function env_value

  logical function env_knows(env, member)

    type(graph), intent(in) :: env, member

    type(graph) :: p, s
    integer     :: k

    env_knows = .false.
    do k = 1, seq_len(env)
       p = seq_item(env, k)
       s = p % known(1)
       if (s % same_as(member)) then
          env_knows = .true.
          return
       end if
    end do

  end function env_knows

  !===================================================================!
  ! ONE evaluator. A value atom answers its value; a bare atom is a
  ! member looked up in the environment; a pair is an application:
  ! branch(1) the operation atom, branch(2) the argument sequence.
  ! Meaning attaches to the operation's IDENTITY through this
  ! reading - the calculator's law table, fractally. The same
  ! function evaluates (2+3)x4 and a residual's right-hand side.
  !===================================================================!

  function evaluate(expr, plus, times, env) result(v)

    type(graph), intent(in)           :: expr, plus, times
    type(graph), intent(in), optional :: env
    real(dp)                          :: v

    v = eval_bounded(expr, plus, times, expr % universe_size(), env)

  end function evaluate

  recursive function eval_bounded(expr, plus, times, budget, env) &
       & result(v)

    type(graph), intent(in)           :: expr, plus, times
    integer    , intent(in)           :: budget
    type(graph), intent(in), optional :: env
    real(dp)                          :: v

    type(graph) :: op, args
    integer     :: k, n

    ! any terminating path visits at most the population once; a
    ! deeper path has repeated a node, and a cycle has no value.
    if (budget .lt. 0) then
       error stop 'fractal view: this expression never bottoms out - a cycle has no value'
    end if

    if (expr % carries_value()) then
       v = expr % value()
       return
    end if

    if (is_atom(expr)) then
       if (.not. present(env)) then
          error stop 'fractal view: a bare atom means nothing without an environment'
       end if
       v = env_value(env, expr)
       return
    end if

    op = expr % known(1)

    ! the empty argument sequence is the NULL branch at its holder,
    ! and the identities answer it: a sum over nobody is 0, a
    ! product over nobody is 1.
    if (expr % status(2) .eq. BRANCH_NULL) then
       if (op % same_as(plus)) then
          v = 0.0_dp
       else if (op % same_as(times)) then
          v = 1.0_dp
       else
          error stop 'fractal view: this reading binds no law to that operation'
       end if
       return
    end if

    args = expr % known(2)
    n    = seq_len(args)

    if (op % same_as(plus)) then
       v = 0.0_dp
       do k = 1, n
          v = v + eval_bounded(seq_item(args, k), plus, times, &
               &               budget - 1, env)
       end do
    else if (op % same_as(times)) then
       v = 1.0_dp
       do k = 1, n
          v = v * eval_bounded(seq_item(args, k), plus, times, &
               &               budget - 1, env)
       end do
    else
       error stop 'fractal view: this reading binds no law to that operation'
    end if

  end function eval_bounded

  !===================================================================!
  ! The storage law, demonstrated: value atoms minted consecutively
  ! form a block, and a holder of the first may read the k-th in
  ! one step. Semantics above, contiguity below - CSR's whole idea,
  ! kept as a promise a view holds, never as a second ontology.
  !===================================================================!

  real(dp) function block_value(a, first, k) result(v)

    type(graph_arena), target, intent(in) :: a
    type(graph)              , intent(in) :: first
    integer                  , intent(in) :: k

    type(graph) :: g

    if (k .lt. 1) error stop 'fractal view: a block counts from one'
    g = a % citizen(first % seat() + k - 1)
    v = g % value()

  end function block_value

  !===================================================================!
  ! The canonical names, for diagnostics and for nobody's synonyms.
  ! A view, not a generator: whoever accepts raw integers guards
  ! them.
  !===================================================================!

  function branch_name(s) result(said)

    integer, intent(in)           :: s
    character(len=:), allocatable :: said

    select case (s)
    case (BRANCH_NULL);    said = 'null'
    case (BRANCH_UNKNOWN); said = 'unknown'
    case (BRANCH_KNOWN);   said = 'known'
    case default
       error stop 'fractal view: there is no fourth status'
    end select

  end function branch_name

end module fractal_views
