!=====================================================================!
! CALCULATOR TOWER . LEVEL 6 . DISCRETIZATION
!
! The level answers one question: HOW STRUCTURE BECOMES AN EQUATION
! LAYOUT. One genuinely new fact enters - the LOCATION relation,
!
!      L <= Y x X        (r_c, c)   (r_e, e)
!
! where Y = { r_c, r_e } are residual rows, born at this level and
! nowhere below. EVERYTHING ELSE IS ALGEBRA: which operation
! produces each computed slot, which operation each row belongs to,
! which slots each operation touches - restricted, projected,
! composed - until the structural Jacobian
!
!      J <= Y x X        r_c -> {a,b,c}      r_e -> {c,d,e}
!
! stands DERIVED, its six pairs never written down. The reverse
! structure is the transpose VIEW of the same J - one dependency
! description for tangent and adjoint, never two. No arithmetic
! anywhere: no 2, no 3, no sum, no product - participation is
! structure, and meaning waits for Level 8.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program calculator_level_6

  use calculator_assert, only : report, verdict
  use calculator_assert, only : SLOT_A, SLOT_B, SLOT_C, SLOT_D, SLOT_E
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : counted_set, subset_set, member_set
  use graph_relation   , only : stored_relation, relation
  use graph_relation_algebra, only : restrict_slot, project_slots, &
       &                             compose_binary
  use graph_binary_relation , only : csr_relation, transposed_view, &
       &                             transpose_of

  implicit none

  ! Residual rows exist only from this level upward.
  integer, parameter :: ROW_C = 1
  integer, parameter :: ROW_E = 2

  type(counted_set)          :: x, o, p, y
  type(subset_set)           :: p_out
  type(stored_relation)      :: flow, backwards, located
  type(stored_relation)      :: q_prod, a_part
  type(csr_relation), target :: j, j2
  type(transposed_view)      :: jt
  integer                    :: table(3, 6)
  integer                    :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "calculator tower . level 6 . discretization"
  write(*,'(1x,a)') "============================================="

  x = counted_set('value-slots'  , 5)
  o = counted_set('operations'   , 2)
  p = counted_set('ports'        , 3)
  y = counted_set('residual-rows', 2)

  table(:, 1) = [OP_PLUS , SLOT_A, PORT_IN1]
  table(:, 2) = [OP_PLUS , SLOT_B, PORT_IN2]
  table(:, 3) = [OP_PLUS , SLOT_C, PORT_OUT]
  table(:, 4) = [OP_TIMES, SLOT_C, PORT_IN1]
  table(:, 5) = [OP_TIMES, SLOT_D, PORT_IN2]
  table(:, 6) = [OP_TIMES, SLOT_E, PORT_OUT]
  flow = stored_relation('flow', [o, x, p], table)

  ! THE one new level-6 fact: where each residual row is located.
  located = stored_relation('located', [y, x], &
       & reshape([ROW_C, SLOT_C,  ROW_E, SLOT_E], [2, 2]))

  p_out = subset_set('output-port', p, [PORT_OUT])

  j = derive_jacobian(flow)

  call check_jacobian(j, "derived", nfail)
  call check_row_supports(nfail)
  call check_order_invariance(nfail)
  call check_reverse_is_the_view(nfail)

  call verdict(nfail, "level 6")

contains

  !===================================================================!
  ! The whole derivation, algebra end to end: Q who-produces,
  ! B row-to-operation, A operation-participation, J = A o B.
  !===================================================================!

  function derive_jacobian(t_flow) result(jac)

    type(stored_relation), intent(in) :: t_flow
    type(csr_relation)                :: jac

    type(stored_relation)        :: q_, a_
    class(relation), allocatable :: b_

    q_  = project_slots(restrict_slot(t_flow, 3, p_out), [2, 1])
    b_  = compose_binary(located, q_)
    a_  = project_slots(t_flow, [1, 2])
    jac = compose_binary(b_, a_)

    ! keep the intermediates visible for their own checks
    q_prod = q_
    a_part = a_

  end function derive_jacobian

  !===================================================================!
  ! J holds exactly the six architect-owned pairs.
  !===================================================================!

  subroutine check_jacobian(jac, what, nfail)

    type(csr_relation), intent(in) :: jac
    character(len=*)  , intent(in) :: what
    integer           , intent(inout) :: nfail

    class(member_set), allocatable :: dom

    dom = jac % domain(1)
    call report(dom % same_as(y), &
         & "J runs from the residual rows (" // what // ")", nfail)
    dom = jac % domain(2)
    call report(dom % same_as(x), &
         & "into the value slots", nfail)

    call report(jac % num_tuples() .eq. 6, &
         & "|J| = 6: no pair beyond the six", nfail)
    call report(jac % has([ROW_C, SLOT_A]) .and. &
         &      jac % has([ROW_C, SLOT_B]) .and. &
         &      jac % has([ROW_C, SLOT_C]), &
         & "r_c depends on a, b and c", nfail)
    call report(jac % has([ROW_E, SLOT_C]) .and. &
         &      jac % has([ROW_E, SLOT_D]) .and. &
         &      jac % has([ROW_E, SLOT_E]), &
         & "r_e depends on c, d and e", nfail)

  end subroutine check_jacobian

  !===================================================================!
  ! Row supports, composed locally from membership by scanning X in
  ! its own order - the generation rule at work: no production
  ! support_of_row() exists, none is needed.
  !===================================================================!

  subroutine check_row_supports(nfail)

    integer, intent(inout) :: nfail

    type(subset_set)     :: src, sre
    integer              :: i, n1, n2
    integer              :: kept1(x % size()), kept2(x % size())

    n1 = 0
    n2 = 0
    do i = 1, x % size()
       if (j % has([ROW_C, x % member(i)])) then
          n1 = n1 + 1
          kept1(n1) = x % member(i)
       end if
       if (j % has([ROW_E, x % member(i)])) then
          n2 = n2 + 1
          kept2(n2) = x % member(i)
       end if
    end do
    src = subset_set('support of r_c', x, kept1(1:n1))
    sre = subset_set('support of r_e', x, kept2(1:n2))

    call report(src % size() .eq. 3 .and. src % has(SLOT_A) .and. &
         &      src % has(SLOT_B) .and. src % has(SLOT_C), &
         & "support(r_c) = { a, b, c }, composed, not stored", nfail)
    call report(sre % size() .eq. 3 .and. sre % has(SLOT_C) .and. &
         &      sre % has(SLOT_D) .and. sre % has(SLOT_E), &
         & "support(r_e) = { c, d, e }", nfail)
    call report(src % is_subobject_of(x) .and. sre % is_subobject_of(x), &
         & "and both stand embedded in the value slots", nfail)

  end subroutine check_row_supports

  !===================================================================!
  ! The Jacobian cannot remember how the flow was written: the same
  ! six facts handed backwards derive the same J, as sets.
  !===================================================================!

  subroutine check_order_invariance(nfail)

    integer, intent(inout) :: nfail

    integer              :: rev(3, 6), k
    integer, allocatable :: jt_(:,:)
    logical              :: ok

    do k = 1, 6
       rev(:, k) = table(:, 7 - k)
    end do
    backwards = stored_relation('flow backwards', [o, x, p], rev)

    j2 = derive_jacobian(backwards)

    call report(j2 % num_tuples() .eq. j % num_tuples(), &
         & "|J1| = |J2|", nfail)
    call j % tuples(jt_)
    ok = .true.
    do k = 1, size(jt_, 2)
       ok = ok .and. j2 % has(jt_(:, k))
    end do
    call report(ok, &
         & "every pair of J1 stands in J2: equal as sets", nfail)

  end subroutine check_order_invariance

  !===================================================================!
  ! One dependency description: the reverse structure is the
  ! transpose VIEW of J - the six inverted pairs, no second stored
  ! pattern anywhere.
  !===================================================================!

  subroutine check_reverse_is_the_view(nfail)

    integer, intent(inout) :: nfail

    jt = transpose_of(j)

    call report(jt % has([SLOT_A, ROW_C]) .and. &
         &      jt % has([SLOT_B, ROW_C]) .and. &
         &      jt % has([SLOT_C, ROW_C]) .and. &
         &      jt % has([SLOT_C, ROW_E]) .and. &
         &      jt % has([SLOT_D, ROW_E]) .and. &
         &      jt % has([SLOT_E, ROW_E]), &
         & "the reverse pattern is the transpose view, pair for pair", &
         & nfail)
    call report(jt % num_tuples() .eq. 6, &
         & "and holds exactly the six, no more", nfail)

  end subroutine check_reverse_is_the_view

end program calculator_level_6
