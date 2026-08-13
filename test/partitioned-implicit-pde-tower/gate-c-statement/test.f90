!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . GATE C . THE STATEMENT
!
! The sealing gate. It asks whether
!
!      partition + overlap refresh + local topology actions
!                + owned assembly
!
! compose into ONE global matrix-free operation that ordinary
! production GMRES can solve with:
!
!      A_partitioned(v) = A_global(v)   for arbitrary v
!
!      GMRES(A_part, b) = GMRES(A_global, b) = q*
!
! Four probes settle the first law - the exact state, a mixed-sign
! vector, and the two INTERFACE basis vectors e3 and e4, which are
! the ones that force information across the cut. Then the same
! composite instance is applied five times in a row, ending where
! it began: y1 = y5. That sequence is what convicts a stale halo.
! The structure is cut ONCE at construction; the overlap values are
! rebuilt on EVERY call, from the state handed in at that moment,
! because the composite owns no numerical state to reuse.
!
! Equivalence is proved twice over: directly between the fixtures,
! and again at solver % matvec - the exact seat GMRES consumes -
! because that is where the reverse review's question actually
! lives.
!
! And the composite is decomposition-context-bound: the six-vertex
! star that Gate B used to convict the host is refused here, in
! domain() before attach can even finish. Same cardinality, wrong
! decomposition.
!
! What this is NOT: distributed. One process, one global trial
! vector, sequential parts, in-process assembly, global-array inner
! products and norms. The honest name is a PARTITIONED MATRIX-FREE
! SOLVE - a serial semantic model of a partitioned operator.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_gate_c

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, verdict
  use partitioned_pde_assert, only : NV, Q_EXACT, B_EXACT
  use graph_carrier    , only : counted_set, member_set
  use graph_grammar    , only : graph, graph_field
  use class_graph      , only : stored_graph
  use class_graph_field, only : field
  use class_graph_gmres, only : gmres
  use shifted_laplacian_fixture, only : shifted_laplacian
  use partitioned_shifted_laplacian_fixture, only : &
       & partitioned_shifted_laplacian

  implicit none

  ! The four probes, in global vertex order.
  real(dp), parameter :: MIXED(NV) = &
       & [3.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, -9.0_dp]
  real(dp), parameter :: E3(NV) = &
       & [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: E4(NV) = &
       & [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]

  type(stored_graph)                  :: g
  type(shifted_laplacian)             :: direct
  type(partitioned_shifted_laplacian) :: composite
  type(gmres)                         :: solver_global, solver_part
  real(dp)                            :: q_part(NV), q_global(NV)
  integer                             :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . gate C . the statement"
  write(*,'(1x,a)') "============================================="

  g = stored_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])

  ! STRUCTURAL PARTITION - once, here, and never again.
  composite = partitioned_shifted_laplacian(g)

  call check_composite_context(nfail)
  call check_action_probes(nfail)
  call check_no_stale_overlap(nfail)
  call check_solver_matvec_equivalence(nfail)
  call check_the_two_solves(nfail)

  ! The tower's answer: the partitioned solution, as it stands.
  write(*,'(1x,a)', advance='no') "PARTITIONED_PDE_RESULT ="
  call write_field(q_part)

  call verdict(nfail, "gate C")

contains

  !===================================================================!
  ! The composite answers on the decomposition it was built from -
  ! V(G) to V(G) - and on nothing else.
  !===================================================================!

  subroutine check_composite_context(nfail)

    integer, intent(inout) :: nfail

    class(member_set), allocatable :: dom
    type(counted_set)              :: vs

    vs = g % vertex_set()
    call composite % domain(g, dom)
    call report(dom % same_as(vs), &
         & "the composite maps V(G) -> V(G), by identity", nfail)

  end subroutine check_composite_context

  !===================================================================!
  ! The extensional law, on four probes. The two interface basis
  ! vectors are the ones that force information across the cut: e3
  ! is owned by part 1 and borrowed by part 2, e4 the reverse.
  !===================================================================!

  subroutine check_action_probes(nfail)

    integer, intent(inout) :: nfail

    call one_probe(Q_EXACT, "q*", nfail)
    call one_probe(MIXED  , "the mixed-sign probe", nfail)
    call one_probe(E3     , "the interface basis e3", nfail)
    call one_probe(E4     , "the interface basis e4", nfail)

  end subroutine check_action_probes

  subroutine one_probe(v, tag, nfail)

    real(dp)        , intent(in)    :: v(:)
    character(len=*), intent(in)    :: tag
    integer         , intent(inout) :: nfail

    real(dp) :: direct_answer(NV), part_answer(NV)

    direct_answer = act(direct_on(v))
    part_answer   = act(part_on(v))

    call report(maxval(abs(direct_answer - part_answer)) < 1.0d-12, &
         & "A_partitioned = A_global on " // tag // ", by member", &
         & nfail)

  end subroutine one_probe

  !===================================================================!
  ! Five applications of ONE composite instance, ending where they
  ! began. A halo cached from a previous matvec would survive into
  ! the next and break y1 = y5 - and would break the interleaved
  ! probes long before that.
  !===================================================================!

  subroutine check_no_stale_overlap(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: y1(NV), y2(NV), y3(NV), y4(NV), y5(NV)
    logical  :: ok

    y1 = act(part_on(Q_EXACT))
    y2 = act(part_on(MIXED))
    y3 = act(part_on(E3))
    y4 = act(part_on(E4))
    y5 = act(part_on(Q_EXACT))

    call report(maxval(abs(y1 - y5)) < 1.0d-14, &
         & "the same composite, applied five times, returns to its " // &
         & "first answer: no overlap survives a call", nfail)

    ok = maxval(abs(y1 - act(direct_on(Q_EXACT)))) < 1.0d-12 .and. &
       & maxval(abs(y2 - act(direct_on(MIXED  )))) < 1.0d-12 .and. &
       & maxval(abs(y3 - act(direct_on(E3     )))) < 1.0d-12 .and. &
       & maxval(abs(y4 - act(direct_on(E4     )))) < 1.0d-12
    call report(ok, &
         & "and every one of the five agreed with the global action " // &
         & "at the time it was asked", nfail)

    call report(maxval(abs(y1 - y2)) > 1.0_dp, &
         & "the interleaved probes really did differ - the sequence " // &
         & "is not a repeated no-op", nfail)

  end subroutine check_no_stale_overlap

  !===================================================================!
  ! The same equivalence at the seat GMRES actually consumes.
  !===================================================================!

  subroutine check_solver_matvec_equivalence(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)     :: vs
    real(dp), allocatable :: gv(:)

    vs = g % vertex_set()

    call solver_global % attach(direct, g, vs)
    call solver_part % attach(composite, g, vs)
    solver_global % tolerance      = 1.0d-12
    solver_global % max_iterations = 200
    solver_part % tolerance        = 1.0d-12
    solver_part % max_iterations   = 200

    call solver_global % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "the global action's affine constant is zero", nfail)
    call solver_part % constant(gv)
    call report(maxval(abs(gv)) < 1.0d-12, &
         & "and so is the partitioned action's: both are linear", nfail)

    call one_matvec(Q_EXACT, "q*", nfail)
    call one_matvec(MIXED  , "the mixed-sign probe", nfail)
    call one_matvec(E3     , "the interface basis e3", nfail)
    call one_matvec(E4     , "the interface basis e4", nfail)

  end subroutine check_solver_matvec_equivalence

  subroutine one_matvec(v, tag, nfail)

    real(dp)        , intent(in)    :: v(:)
    character(len=*), intent(in)    :: tag
    integer         , intent(inout) :: nfail

    real(dp), allocatable :: yg(:), yp(:)

    call solver_global % matvec(v, yg)
    call solver_part % matvec(v, yp)

    call report(maxval(abs(yg - yp)) < 1.0d-12, &
         & "solver matvec agrees on " // tag // ": the composite is " // &
         & "the operator GMRES sees", nfail)

  end subroutine one_matvec

  !===================================================================!
  ! The statement: solve A q = b twice - once through the global
  ! action, once through the partitioned one - and meet on q*.
  !===================================================================!

  subroutine check_the_two_solves(nfail)

    integer, intent(inout) :: nfail

    type(field)                     :: rhs
    type(counted_set)               :: vs
    class(graph_field), allocatable :: sol
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: v(:)

    vs = g % vertex_set()
    rhs = field('b', vs)
    call rhs % set_real_vector(B_EXACT)

    ! The partitioned solve - the tower's own road.
    call solver_part % apply(g, [rhs], sol)
    call sol % domain(dom)
    call sol % get_real_vector(v)
    q_part = v(1:NV)

    call report(dom % same_as(vs), &
         & "the partitioned solution is a field on V(G)", nfail)
    call report(by_member(q_part, vs, Q_EXACT), &
         & "and it is q* = [1,2,4,7,11,16], by global member", nfail)

    ! The global baseline, run independently.
    call solver_global % apply(g, [rhs], sol)
    call sol % get_real_vector(v)
    q_global = v(1:NV)

    call report(by_member(q_global, vs, Q_EXACT), &
         & "the global solve independently reaches q* as well", nfail)
    call report(maxval(abs(q_part - q_global)) < 1.0d-9, &
         & "q_partitioned = q_global = q*: the decomposition changed " // &
         & "the road, not the answer", nfail)

  end subroutine check_the_two_solves

  !===================================================================!
  ! Helpers. act() runs an operation and returns its values; the
  ! two *_on builders make the state field each road wants.
  !===================================================================!

  function direct_on(v) result(answer)

    real(dp), intent(in)            :: v(:)
    class(graph_field), allocatable :: answer

    type(field) :: q

    q = field('probe', g % vertex_set())
    call q % set_real_vector(v)
    call direct % apply(g, [q], answer)

  end function direct_on

  function part_on(v) result(answer)

    real(dp), intent(in)            :: v(:)
    class(graph_field), allocatable :: answer

    type(field) :: q

    q = field('probe', g % vertex_set())
    call q % set_real_vector(v)
    call composite % apply(g, [q], answer)

  end function part_on

  function act(answer) result(v)

    class(graph_field), intent(in) :: answer
    real(dp)                       :: v(NV)

    real(dp), allocatable :: got(:)

    call answer % get_real_vector(got)
    v = got(1:NV)

  end function act

  logical function by_member(v, dom, expect)

    real(dp)         , intent(in) :: v(:)
    class(member_set), intent(in) :: dom
    real(dp)         , intent(in) :: expect(:)

    integer :: i, m

    by_member = .true.
    do i = 1, dom % size()
       m = dom % member(i)
       by_member = by_member .and. &
            & (abs(v(dom % local_index(m)) - expect(m)) < 1.0d-9)
    end do

  end function by_member

  ! The solution field, serialized as it stands - one real per
  ! global vertex, in global enumeration order, unrounded.
  subroutine write_field(v)

    real(dp), intent(in) :: v(:)

    integer :: i

    do i = 1, size(v)
       write(*,'(es24.16)', advance='no') v(i)
    end do
    write(*,'(a)') ""

  end subroutine write_field

end program partitioned_pde_gate_c
