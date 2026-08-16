!=====================================================================!
! The fractal kernel suite: one primitive, four inhabitations.
!
! The laws first - identity between two atoms nothing else could
! tell apart, sharing without duplication, a cycle tied through
! the one lawful transition, nine status pairs as data, growth
! that invalidates nobody, one prescribed encoding for tuples,
! and the two absences kept apart at the primitive itself.
!
! Then the four inhabitations of ONE kernel, unchanged between
! them: the bare atom; the calculator (2+3)x4 = 20; a strip of
! finite-volume cells with boundary faces as NULL and a field
! read both semantically and through the contiguity promise; and
! the computational pair G = (Q, R) whose Q branch is realized by
! solving - UNKNOWN to KNOWN, once, with realized still not
! meaning solved.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_fractal_kernel

  use iso_fortran_env, only : dp => REAL64, int64
  use fractal_graph  , only : graph_arena, graph, branch_spec, &
       &                      BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN, &
       &                      null_branch, unknown_branch, known_branch
  use fractal_views  , only : is_atom, pair_of, seq_of, seq_len, &
       &                      seq_item, rel_has, rel_image_count, &
       &                      tail_of, head_of, out_degree, &
       &                      neighbour_count, env_of, env_bind, &
       &                      env_value, env_knows, evaluate, block_value

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "fractal kernel suite (feasibility spike)"
  write(*,'(1x,a)') "============================================="

  call check_atoms_and_identity(nfail)
  call check_two_absences(nfail)
  call check_nine_states(nfail)
  call check_sharing(nfail)
  call check_cycle(nfail)
  call check_growth_keeps_handles(nfail)
  call check_encoding_law(nfail)
  call check_values(nfail)
  call check_calculator(nfail)
  call check_finite_volume(nfail)
  call check_computational_pair(nfail)
  call bench()

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all fractal checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " fractal check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Inhabitation 1: G = (NULL, NULL) is a graph. Two of them are
  ! two atoms, and only identity can say so.
  !===================================================================!

  subroutine check_atoms_and_identity(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: x, y, never

    a = graph_arena()
    x = a % node(null_branch(), null_branch())
    y = a % node(null_branch(), null_branch())

    call report(is_atom(x) .and. is_atom(y), &
         & "a (null, null) node is already a graph - the atom", nfail)
    call report(x % same_as(x) .and. .not. x % same_as(y), &
         & "two atoms of one shape are two atoms, by seat alone", nfail)
    call report(.not. never % same_as(never), &
         & "an unhatched handle equals nothing, itself included", nfail)

  end subroutine check_atoms_and_identity

  !===================================================================!
  ! NULL is not UNKNOWN: absence and bottom are different absences,
  ! and the primitive itself keeps them apart.
  !===================================================================!

  subroutine check_two_absences(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: g

    a = graph_arena()
    g = a % node(unknown_branch(), null_branch())

    call report(g % status(1) .eq. BRANCH_UNKNOWN .and. &
         &      g % status(2) .eq. BRANCH_NULL .and. &
         &      g % status(1) .ne. g % status(2), &
         & "bottom and absence are two statuses, never one storage trick", &
         & nfail)

  end subroutine check_two_absences

  !===================================================================!
  ! Nine status pairs, one type: 3 x 3 combinations are values of
  ! data, and no subclass anywhere answers for any of them.
  !===================================================================!

  subroutine check_nine_states(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: t, g
    type(branch_spec)         :: specs(3)
    integer                   :: i, j, seen
    logical                   :: ok

    a = graph_arena()
    t = a % node(null_branch(), null_branch())
    specs(1) = null_branch()
    specs(2) = unknown_branch()
    specs(3) = known_branch(t)

    seen = 0
    ok   = .true.
    do i = 1, 3
       do j = 1, 3
          g = a % node(specs(i), specs(j))
          ok = ok .and. g % status(1) .eq. i - 1 &
               & .and.  g % status(2) .eq. j - 1
          seen = seen + 1
       end do
    end do

    call report(ok .and. seen .eq. 9, &
         & "nine branch-state combinations, zero subclasses", nfail)

  end subroutine check_nine_states

  !===================================================================!
  ! Sharing: one realization, two branches, no duplication - the
  ! diamond A = (B, B) holds the SAME B twice.
  !===================================================================!

  subroutine check_sharing(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: b, d, l, r

    a = graph_arena()
    b = a % node(null_branch(), null_branch())
    d = a % node(known_branch(b), known_branch(b))

    l = d % known(1)
    r = d % known(2)
    call report(l % same_as(r) .and. l % same_as(b), &
         & "both branches hold the SAME realization, not a twin", nfail)

  end subroutine check_sharing

  !===================================================================!
  ! Cycles: the knot is tied by the one lawful transition. A is
  ! minted ignorant, B points at A, and realizing A's branch at B
  ! closes the loop - a graph, not a tree.
  !===================================================================!

  subroutine check_cycle(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: ga, gb, walk

    a = graph_arena()
    ga = a % node(unknown_branch(), null_branch())
    gb = a % node(known_branch(ga), null_branch())
    call ga % realize(1, gb)

    walk = ga % known(1)          ! A -> B
    walk = walk % known(1)        ! B -> A
    call report(walk % same_as(ga) .and. &
         &      ga % status(1) .eq. BRANCH_KNOWN, &
         & "two hops return home: the representation is a graph", nfail)

  end subroutine check_cycle

  !===================================================================!
  ! Growth: a seat is an index, never an address. Reallocation
  ! behind the arena invalidates no handle anyone holds.
  !===================================================================!

  subroutine check_growth_keeps_handles(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: first, late, again
    integer                   :: k

    a = graph_arena()
    first = a % node(unknown_branch(), null_branch())
    do k = 1, 300
       late = a % node(null_branch(), null_branch())
    end do

    again = a % citizen(1)
    call report(first % same_as(again) .and. &
         &      first % status(1) .eq. BRANCH_UNKNOWN, &
         & "three hundred births later, the first handle still answers", &
         & nfail)

    call first % realize(1, late)
    call report(first % status(1) .eq. BRANCH_KNOWN, &
         & "and can still learn: the transition survives the growth", &
         & nfail)

  end subroutine check_growth_keeps_handles

  !===================================================================!
  ! The encoding law: (a, b, c) has one prescribed spelling - the
  ! right-nested spine, elements on branch one, NULL closing the
  ! last cell. A one-item sequence is a cell, never the bare atom.
  !===================================================================!

  subroutine check_encoding_law(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: x, y, z, t, cell, e
    type(graph)               :: items(3)
    logical                   :: ok

    a = graph_arena()
    x = a % node(null_branch(), null_branch())
    y = a % node(null_branch(), null_branch())
    z = a % node(null_branch(), null_branch())

    items(1) = x
    items(2) = y
    items(3) = z
    t = seq_of(a, items)

    e = seq_item(t, 1)
    ok = e % same_as(x)
    e = seq_item(t, 2)
    ok = ok .and. e % same_as(y)
    e = seq_item(t, 3)
    ok = ok .and. e % same_as(z)
    call report(seq_len(t) .eq. 3 .and. ok, &
         & "(a, b, c) reads back as a, b, c - one spelling, one law", &
         & nfail)

    cell = t % known(2)
    cell = cell % known(2)
    call report(cell % status(2) .eq. BRANCH_NULL, &
         & "the last cell closes with NULL - no terminator atom exists", &
         & nfail)

    t = seq_of(a, items(1:1))
    call report(seq_len(t) .eq. 1 .and. .not. is_atom(t), &
         & "a one-item sequence is a cell, never the bare atom", nfail)

  end subroutine check_encoding_law

  !===================================================================!
  ! Values: an atom may carry one number - the compressed numeral.
  ! And the two absences persist at the container: a NULL branch is
  ! an empty collection; an UNKNOWN branch is a collection not yet
  ! learned. Bottom is not empty, here as everywhere.
  !===================================================================!

  subroutine check_values(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: v, bare, holder_empty, holder_ignorant

    a = graph_arena()
    v    = a % value_atom(2.0_dp)
    bare = a % node(null_branch(), null_branch())

    call report(v % carries_value() .and. v % value() .eq. 2.0_dp &
         & .and. .not. bare % carries_value() .and. is_atom(v), &
         & "a value atom is still an atom; a bare atom carries nothing", &
         & nfail)

    holder_empty    = a % node(null_branch(),    null_branch())
    holder_ignorant = a % node(unknown_branch(), null_branch())
    call report(holder_empty % status(1)    .eq. BRANCH_NULL .and. &
         &      holder_ignorant % status(1) .eq. BRANCH_UNKNOWN, &
         & "an empty collection is NULL; an unlearned one is UNKNOWN", &
         & nfail)

  end subroutine check_values

  !===================================================================!
  ! Inhabitation 2: the calculator. Operations are ATOMS - meaning
  ! attaches to identity through the reading, exactly as the
  ! constitution has always said - and (2+3)x4 answers 20.
  !===================================================================!

  subroutine check_calculator(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: plus, times, two, three, four
    type(graph)               :: sum_args(2), prod_args(2)
    type(graph)               :: sum_seq, prod_seq, sum_expr, expr
    real(dp)                  :: answer

    a = graph_arena()
    plus  = a % node(null_branch(), null_branch())
    times = a % node(null_branch(), null_branch())
    two   = a % value_atom(2.0_dp)
    three = a % value_atom(3.0_dp)
    four  = a % value_atom(4.0_dp)

    sum_args(1) = two
    sum_args(2) = three
    sum_seq  = seq_of(a, sum_args)
    sum_expr = pair_of(a, plus, sum_seq)

    prod_args(1) = sum_expr
    prod_args(2) = four
    prod_seq = seq_of(a, prod_args)
    expr = pair_of(a, times, prod_seq)

    answer = evaluate(expr, plus, times)
    call report(answer .eq. 20.0_dp, &
         & "(2 + 3) x 4 = 20, read fractally off four atoms", nfail)
    write(*,'(1x,a,i0)') "FRACTAL_CALCULATOR_RESULT = ", nint(answer)

    sum_expr = a % node(known_branch(plus), null_branch())
    expr     = a % node(known_branch(times), null_branch())
    call report(evaluate(sum_expr, plus, times) .eq. 0.0_dp .and. &
         &      evaluate(expr, plus, times) .eq. 1.0_dp, &
         & "a sum over nobody is 0, a product over nobody is 1", nfail)

  end subroutine check_calculator

  !===================================================================!
  ! Inhabitation 3: a strip of finite-volume cells. Faces are
  ! pairs; a boundary face's far side is NULL - the absence law,
  ! now primitive. One list of pairs is read three ways: relation,
  ! directed adjacency, ordinary neighbourhood. The field rides as
  ! contiguous value atoms, read semantically and through the
  ! contiguity promise, and interior differences telescope.
  !===================================================================!

  subroutine check_finite_volume(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: c(4), q(4), faces(5), t, h
    type(graph)               :: incidence, field, cells_seq, vals_seq
    real(dp)                  :: qv(4), flux, semantic, blocked
    integer                   :: k
    logical                   :: ok

    a = graph_arena()
    do k = 1, 4
       c(k) = a % node(null_branch(), null_branch())
    end do

    faces(1) = a % node(known_branch(c(1)), null_branch())
    faces(2) = pair_of(a, c(1), c(2))
    faces(3) = pair_of(a, c(2), c(3))
    faces(4) = pair_of(a, c(3), c(4))
    faces(5) = a % node(known_branch(c(4)), null_branch())
    incidence = seq_of(a, faces)

    call report(rel_has(incidence, c(1), c(2)) .and. &
         &      .not. rel_has(incidence, c(2), c(1)), &
         & "the relation reading: (c1, c2) stands, (c2, c1) does not", &
         & nfail)

    t = tail_of(incidence, 2)
    h = head_of(incidence, 2)
    call report(t % same_as(c(1)) .and. h % same_as(c(2)) .and. &
         &      out_degree(incidence, c(1)) .eq. 1, &
         & "the directed reading: tails, heads, degrees - same pairs", &
         & nfail)

    call report(neighbour_count(incidence, c(2)) .eq. 2 .and. &
         &      neighbour_count(incidence, c(1)) .eq. 1, &
         & "the ordinary reading forgets order; boundaries add nobody", &
         & nfail)

    qv = [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp]
    do k = 1, 4
       q(k) = a % value_atom(qv(k))
    end do
    cells_seq = seq_of(a, c)
    vals_seq  = seq_of(a, q)
    field = pair_of(a, cells_seq, vals_seq)

    ok = .true.
    do k = 1, 4
       t = field % known(2)
       t = seq_item(t, k)
       semantic = t % value()
       blocked  = block_value(a, q(1), k)
       ok = ok .and. semantic .eq. qv(k) .and. blocked .eq. qv(k)
    end do
    call report(ok, &
         & "the field answers alike down the semantic and block roads", &
         & nfail)

    ! the member's local index is seat arithmetic on the cell block;
    ! its value is the same index into the value block - local_index,
    ! compiled down to two subtractions.
    flux = 0.0_dp
    do k = 2, 4
       t = tail_of(incidence, k)
       h = head_of(incidence, k)
       flux = flux + block_value(a, q(1), h % seat() - c(1) % seat() + 1) &
            &      - block_value(a, q(1), t % seat() - c(1) % seat() + 1)
    end do
    call report(flux .eq. qv(4) - qv(1), &
         & "interior differences telescope: conservation, fractally", &
         & nfail)

  end subroutine check_finite_volume

  !===================================================================!
  ! Inhabitation 4: the computational pair. G = (Q, R) is a NODE:
  ! branch one the data, branch two the residual, statuses saying
  ! exactly what COMPUTATIONAL-GRAPH.md's four names say. Solving
  ! is the one lawful transition, and realized is still not solved.
  !===================================================================!

  subroutine check_computational_pair(nfail)

    integer, intent(inout) :: nfail

    type(graph_arena), target :: a
    type(graph)               :: qa, qb, qc, qd, qe, plus, times
    type(graph)               :: args(2), eqs(2), pairG, wrong
    type(graph)               :: args_seq, expr
    type(graph)               :: residual, env, eq, target_m, rhs
    real(dp)                  :: v, defect
    integer                   :: k

    a = graph_arena()
    qa = a % node(null_branch(), null_branch())
    qb = a % node(null_branch(), null_branch())
    qc = a % node(null_branch(), null_branch())
    qd = a % node(null_branch(), null_branch())
    qe = a % node(null_branch(), null_branch())
    plus  = a % node(null_branch(), null_branch())
    times = a % node(null_branch(), null_branch())

    args(1) = qa
    args(2) = qb
    args_seq = seq_of(a, args)
    expr     = pair_of(a, plus, args_seq)
    eqs(1)   = pair_of(a, qc, expr)
    args(1) = qc
    args(2) = qd
    args_seq = seq_of(a, args)
    expr     = pair_of(a, times, args_seq)
    eqs(2)   = pair_of(a, qe, expr)
    residual = seq_of(a, eqs)

    pairG = a % node(unknown_branch(), known_branch(residual))
    call report(pairG % status(1) .eq. BRANCH_UNKNOWN .and. &
         &      pairG % status(2) .eq. BRANCH_KNOWN, &
         & "the pair reads 'operator graph': R known, Q bottom", nfail)

    env = env_of(a, qa, 2.0_dp)
    env = env_bind(a, env, qb, 3.0_dp)
    env = env_bind(a, env, qd, 4.0_dp)

    do k = 1, seq_len(residual)
       eq = seq_item(residual, k)
       target_m = eq % known(1)
       rhs      = eq % known(2)
       if (env_knows(env, target_m)) then
          ! the solver's own law: bind each member at most once, and
          ! boundary data never - shadowing is lawful, overwriting
          ! history is not.
          error stop 'solve: a solver binds each member at most once'
       end if
       v = evaluate(rhs, plus, times, env)
       env = env_bind(a, env, target_m, v)
    end do

    call pairG % realize(1, env)
    call report(pairG % status(1) .eq. BRANCH_KNOWN .and. &
         &      pairG % status(2) .eq. BRANCH_KNOWN, &
         & "solve is the transition: bottom became knowledge, once", &
         & nfail)

    defect = 0.0_dp
    do k = 1, seq_len(residual)
       eq = seq_item(residual, k)
       target_m = eq % known(1)
       rhs      = eq % known(2)
       defect = defect + abs(env_value(env, target_m) - &
            &                evaluate(rhs, plus, times, env))
    end do
    call report(defect .eq. 0.0_dp .and. &
         &      env_value(env, qc) .eq. 5.0_dp .and. &
         &      env_value(env, qe) .eq. 20.0_dp, &
         & "R(Q) = 0: q(c) = 5, q(e) = 20 - the old oracle, new bones", &
         & nfail)
    write(*,'(1x,a,i0)') "FRACTAL_REALIZED_RESULT = ", &
         & nint(env_value(env, qe))

    wrong = env_of(a, qc, 99.0_dp)
    wrong = a % node(known_branch(wrong), known_branch(residual))
    call report(wrong % status(1) .eq. BRANCH_KNOWN .and. &
         &      wrong % status(2) .eq. BRANCH_KNOWN, &
         & "an inconsistent pair still reads 'realized' - not 'solved'", &
         & nfail)

  end subroutine check_computational_pair

  !===================================================================!
  ! The honest bill: the semantic road pays per hop; the block road
  ! pays an index. Informational, not asserted - flags are the
  ! suite's own, unoptimized.
  !===================================================================!

  subroutine bench()

    integer, parameter :: n = 100000

    type(graph_arena), target :: a
    type(graph)               :: v, cell, spine, first
    real(dp), allocatable     :: raw(:)
    real(dp)                  :: s_raw, s_block, s_spine
    integer(int64)            :: t0, t1, rate
    integer                   :: k
    real(dp)                  :: ns_raw, ns_block, ns_spine

    a = graph_arena()
    first = a % value_atom(1.0_dp)
    do k = 2, n
       v = a % value_atom(real(k, dp))
    end do

    spine = a % node(known_branch(first), null_branch())
    do k = n, 1, -1
       v = a % citizen(k)
       spine = a % node(known_branch(v), known_branch(spine))
    end do

    allocate(raw(n))
    do k = 1, n
       raw(k) = real(k, dp)
    end do

    call system_clock(t0, rate)
    s_raw = 0.0_dp
    do k = 1, n
       s_raw = s_raw + raw(k)
    end do
    call system_clock(t1)
    ns_raw = 1.0e9_dp * real(t1 - t0, dp) / real(rate, dp) / real(n, dp)

    call system_clock(t0)
    s_block = 0.0_dp
    do k = 1, n
       s_block = s_block + block_value(a, first, k)
    end do
    call system_clock(t1)
    ns_block = 1.0e9_dp * real(t1 - t0, dp) / real(rate, dp) / real(n, dp)

    call system_clock(t0)
    s_spine = 0.0_dp
    cell = spine
    do
       v = cell % known(1)
       if (v % carries_value()) s_spine = s_spine + v % value()
       if (cell % status(2) .ne. BRANCH_KNOWN) exit
       cell = cell % known(2)
    end do
    call system_clock(t1)
    ns_spine = 1.0e9_dp * real(t1 - t0, dp) / real(rate, dp) &
         &   / real(n + 1, dp)

    write(*,'(1x,a)') "---------------------------------------------"
    write(*,'(1x,a,f8.2,a)') "raw array walk ........ ", ns_raw,   " ns/element"
    write(*,'(1x,a,f8.2,a)') "block (storage law) ... ", ns_block, " ns/element"
    write(*,'(1x,a,f8.2,a)') "spine (semantic road) . ", ns_spine, " ns/element"
    write(*,'(1x,a,f20.1)') "checksums agree at ", s_raw
    if (abs(s_raw - s_block) .gt. 0.0_dp) then
       write(*,'(1x,a)') "FAIL : block checksum drifted"
       error stop
    end if

  end subroutine bench

end program test_fractal_kernel
