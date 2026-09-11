!=====================================================================!
! PARTITIONED IMPLICIT PDE TOWER . LEVEL 8 . CONSTITUTION
!
! The level answers one question: CAN THE SAME GLOBAL OPERATOR BE
! CONSTITUTED BY A PARTITIONED REALIZATION INSTEAD OF DIRECT GLOBAL
! TRAVERSAL?
!
!      SAME MATHEMATICS, DIFFERENT COMPUTATIONAL CONSTITUTION.
!
! Level 6 established the local discrete law. This level COMPOSES
! that law into a complete realization of the global operation:
!
!      G -> {G1,G2}     structural partition, ONCE
!      q -> {q1,q2}     overlap refresh, EVERY APPLY
!      A_G1, A_G2       local constituted actions
!      owned assembly   EVERY APPLY
!                       -> A q
!
! The composite owns STRUCTURE and no mutable numerical state: no
! cached q, no cached halo, no cached residual, no previous matvec.
! Nothing in it can go stale because there is nothing to go stale -
! and that is proved, not asserted, by applying ONE instance five
! times with different probes and requiring it to return to its
! first answer.
!
! It is also DECOMPOSITION-CONTEXT-BOUND, unlike the Level-6 action
! which is merely graph-parameterized: its G1 and G2 were cut from
! one particular G, so Level 7's six-vertex star is refused
! outright.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program partitioned_pde_level_8

  use iso_fortran_env  , only : dp => REAL64
  use partitioned_pde_assert, only : report, assert_all
  use partitioned_pde_assert, only : NV, Q_EXACT
  use graph_fractal        , only : graph
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use view_directed_stored      , only : stored_directed_graph
  use field_stored, only : stored_field
  use shifted_laplacian_fixture, only : shifted_laplacian
  use partitioned_shifted_laplacian_fixture, only : &
       & partitioned_shifted_laplacian

  implicit none

  real(dp), parameter :: MIXED(NV) = &
       & [3.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, -9.0_dp]
  real(dp), parameter :: E3(NV) = &
       & [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: E4(NV) = &
       & [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]

  type(stored_directed_graph)                  :: g
  type(shifted_laplacian)             :: direct
  type(partitioned_shifted_laplacian) :: composite
  integer                             :: nfail
  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "partitioned pde tower . level 8 . constitution"
  write(*,'(1x,a)') "============================================="

  g = stored_directed_graph(NV, tails=[1,2,3,4,5], heads=[2,3,4,5,6])

  ! STRUCTURAL PARTITION - once, here, and never again.
  composite = partitioned_shifted_laplacian(g)

  call check_composite_context(nfail)
  call check_action_values(nfail)
  call check_no_stale_overlap(nfail)

  call assert_all(nfail, "level 8")

contains
  !===================================================================!
  ! The composite answers on the decomposition it was built from -
  ! V(G) to V(G) - and on nothing else.
  !===================================================================!

  subroutine check_composite_context(nfail)

    integer, intent(inout) :: nfail

    type(graph) :: dom
    type(graph)              :: vs
    integer         :: n_dom

    vs = g % vertex_set()
    call composite % domain(g, dom, n_dom)
    call report(dom % same_as(vs), &
         & "the composite maps V(G) -> V(G), by identity", nfail)

  end subroutine check_composite_context
  !===================================================================!
  ! The extensional law, on four probes. The two interface basis
  ! vectors are the ones that force information across the cut: e3
  ! is owned by part 1 and borrowed by part 2, e4 the reverse.
  !===================================================================!

  subroutine check_action_values(nfail)

    integer, intent(inout) :: nfail

    call compare_actions(Q_EXACT, "q*", nfail)
    call compare_actions(MIXED  , "the mixed-sign probe", nfail)
    call compare_actions(E3     , "the interface basis e3", nfail)
    call compare_actions(E4     , "the interface basis e4", nfail)

  end subroutine check_action_values
  subroutine compare_actions(v, tag, nfail)

    real(dp)        , intent(in)    :: v(:)
    character(len=*), intent(in)    :: tag
    integer         , intent(inout) :: nfail

    real(dp) :: direct_image(NV), partitioned_image(NV)

    direct_image = act(direct_on(v))
    partitioned_image   = act(part_on(v))

    call report(maxval(abs(direct_image - partitioned_image)) < 1.0d-12, &
         & "A_partitioned = A_global on " // tag // ", by member", &
         & nfail)

  end subroutine compare_actions
  !===================================================================!
  ! Five applications of ONE composite instance, ending where they
  ! began. A halo cached from a previous matvec would survive into
  ! the next and break y1 = y5 - and would break the interleaved
  ! probes long before that.
  !===================================================================!

  subroutine check_no_stale_overlap(nfail)

    integer, intent(inout) :: nfail

    real(dp) :: y1(NV), y2(NV), y3(NV), y4(NV), y5(NV)
    logical  :: satisfied

    y1 = act(part_on(Q_EXACT))
    y2 = act(part_on(MIXED))
    y3 = act(part_on(E3))
    y4 = act(part_on(E4))
    y5 = act(part_on(Q_EXACT))

    call report(maxval(abs(y1 - y5)) < 1.0d-14, &
         & "the same composite, applied five times, returns to its " // &
         & "first answer: no overlap survives a call", nfail)

    satisfied = maxval(abs(y1 - act(direct_on(Q_EXACT)))) < 1.0d-12 .and. &
       & maxval(abs(y2 - act(direct_on(MIXED  )))) < 1.0d-12 .and. &
       & maxval(abs(y3 - act(direct_on(E3     )))) < 1.0d-12 .and. &
       & maxval(abs(y4 - act(direct_on(E4     )))) < 1.0d-12
    call report(satisfied, &
         & "and every one of the five agreed with the global action " // &
         & "at the time it was asked", nfail)

    call report(maxval(abs(y1 - y2)) > 1.0_dp, &
         & "the interleaved probes really did differ - the sequence " // &
         & "is not a repeated no-op", nfail)

  end subroutine check_no_stale_overlap
  !===================================================================!
  ! Helpers. act() runs an operation and returns its values; the
  ! two *_on builders make the state field each road wants.
  !===================================================================!

  function direct_on(v) result(image)

    real(dp), intent(in)            :: v(:)
    class(field), allocatable :: image

    type(stored_field) :: q

    q = stored_field('probe', g % vertex_set(), g % num_vertices())
    call q % set_real_vector(v)
    direct = shifted_laplacian()
    call direct % apply(g, direct % bind([q]), image)

  end function direct_on
  function part_on(v) result(image)

    real(dp), intent(in)            :: v(:)
    class(field), allocatable :: image

    type(stored_field) :: q

    q = stored_field('probe', g % vertex_set(), g % num_vertices())
    call q % set_real_vector(v)
    call composite % apply(g, composite % bind([q]), image)

  end function part_on
  function act(image) result(v)

    class(field), intent(in) :: image
    real(dp)                       :: v(NV)

    real(dp), allocatable :: field_values(:)

    call image % real_vector(field_values)
    v = field_values(1:NV)

  end function act
end program partitioned_pde_level_8
