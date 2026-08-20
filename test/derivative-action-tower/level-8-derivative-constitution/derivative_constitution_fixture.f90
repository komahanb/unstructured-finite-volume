!=====================================================================!
! THE DERIVATIVE CONSTITUTION FIXTURE - test-local, deliberately
! outside src/: one derivative specimen has not earned a production
! derivative API. This module is the ONE place in the tower where
! the symbols mean something:
!
!      product(a, b) = a * b          sum(a, b) = a + b
!
! and where each operation owes exactly ONE more thing: its local
! linearization at given primal inputs,
!
!      L_product(a, b) = [b, a]       L_sum(a, b) = [1, 1]
!
! a row of port-relative coefficients - coeff(1) belongs to in1,
! coeff(2) to in2 - knowing the operation symbol, the primal input
! VALUES, and the port positions, and NOTHING about slots: no
! SLOT_X, no SLOT_Y, no graph identity anywhere. Structure chooses
! the operands; the law sees a symbol and two numbers.
!
! The Gate-B load-bearing design: ONE local linearization supports
! BOTH actions. The tangent applies L_o forward,
!
!      dot(out) = coeff(1)*dot(in1) + coeff(2)*dot(in2)
!
! and the reverse applies its transpose backward,
!
!      bar(in1) += coeff(1)*bar(out)
!      bar(in2) += coeff(2)*bar(out)
!
! There is NO separate reverse derivative formula in this file - the
! transpose action reads the same coefficients, and the += is
! architecturally essential: numerical action occurs per
! operation/input-port incidence, an operation is visited once, and
! contributions landing on the same value slot ACCUMULATE. No path
! is ever enumerated. The reverse accumulator starts at zero
! because zero is addition's identity - the lawful start of an
! accumulator, not a sentinel for "uncomputed" (the primal
! workspace keeps its separate availability flags).
!
! Everything is parameterized by the abstract contracts -
! class(relation), set identities - handed in by the caller: no
! derivative singleton, no stored-type demands, no graph argument.
! What the evaluators consume is T_flow incidences, the DERIVED
! order, primal values, and the law table: the structural J-pattern
! is never an input here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module derivative_constitution_fixture

  use iso_fortran_env, only : dp => REAL64
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_fractal        , only : set_graph => graph
  use map_set        , only : set_map
  use relation_finitary   , only : relation

  implicit none

  private
  public :: apply_law, local_linearization, slot_for_port
  public :: primal_execution, tangent_action, reverse_action

contains

  !===================================================================!
  ! The primal law table: the exact line where meaning enters. An
  ! unbound symbol refuses.
  !===================================================================!

  real(dp) function apply_law(op, a, b) result(c)

    integer , intent(in) :: op
    real(dp), intent(in) :: a, b

    select case (op)
    case (OP_PRODUCT)
       c = a * b
    case (OP_SUM)
       c = a + b
    case default
       error stop 'derivative constitution: no law binds this operation symbol'
    end select

  end function apply_law

  !===================================================================!
  ! The one local linear shadow: port-relative coefficients at the
  ! given primal inputs. Tangent and reverse both read THIS row;
  ! no second derivative formula exists anywhere.
  !===================================================================!

  function local_linearization(op, a, b) result(coeff)

    integer , intent(in) :: op
    real(dp), intent(in) :: a, b
    real(dp)             :: coeff(2)

    select case (op)
    case (OP_PRODUCT)
       coeff = [b, a]
    case (OP_SUM)
       coeff = [1.0_dp, 1.0_dp]
    case default
       error stop 'derivative constitution: no law binds this operation symbol'
    end select

  end function local_linearization

  !===================================================================!
  ! The unique slot standing on a port of an operation - discovered
  ! by scanning the value slots against T_flow, refused at zero or
  ! many. The facts live in the relation, never in the evaluators.
  !===================================================================!

  integer function slot_for_port(flow, slots, sets, op, port) result(found)

    class(relation)  , intent(in) :: flow
    type(set_graph), intent(in) :: slots
    type(set_map)  , intent(in) :: sets
    integer          , intent(in) :: op, port

    integer :: j, hits

    hits  = 0
    found = 0
    do j = 1, sets % size_of(slots)
       if (flow % has([op, sets % member_of(slots, j), port])) then
          hits  = hits + 1
          found = sets % member_of(slots, j)
       end if
    end do
    if (hits .ne. 1) then
       error stop 'derivative constitution: a port names exactly one slot'
    end if

  end function slot_for_port

  !===================================================================!
  ! The primal execution: seed the independent values by DOMAIN
  ! enumeration, execute the caller-DERIVED order into the computed
  ! slots. Availability is a separate flag - no number, zero or
  ! otherwise, ever stands for "not yet". The finished workspace is
  ! the base point at which local linearizations are evaluated.
  !===================================================================!

  subroutine primal_execution(flow, slots, sets, indep, indep_values, &
       & computed, order, values, available)

    class(relation)  , intent(in)  :: flow
    type(set_graph), intent(in)  :: slots      ! V
    type(set_map)  , intent(in) :: sets
    type(set_graph), intent(in)  :: indep      ! X
    real(dp)         , intent(in)  :: indep_values(:)
    type(set_graph), intent(in)  :: computed   ! C
    integer          , intent(in)  :: order(:)   ! derived by caller
    real(dp)         , intent(out) :: values(:)
    logical          , intent(out) :: available(:)

    integer :: i, m, op, in1, in2, out

    if (size(indep_values) .ne. sets % size_of(indep) .or. &
         & size(values) .ne. sets % size_of(slots) .or. &
         & size(available) .ne. sets % size_of(slots)) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    available = .false.
    values    = 0.0_dp

    do i = 1, sets % size_of(indep)
       m = sets % member_of(indep, i)
       values(sets % index_in(slots, m)) = &
            & indep_values(sets % index_in(indep, m))
       available(sets % index_in(slots, m)) = .true.
    end do

    do i = 1, size(order)
       op  = order(i)
       in1 = slot_for_port(flow, slots, sets, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, sets, op, PORT_IN2)
       out = slot_for_port(flow, slots, sets, op, PORT_OUT)

       if (.not. available(sets % index_in(slots, in1)) .or. &
            & .not. available(sets % index_in(slots, in2))) then
          error stop 'derivative constitution: an operation was scheduled before its primal inputs exist'
       end if
       if (.not. sets % has_in(computed, out)) then
          error stop 'derivative constitution: an operation must produce a computed slot'
       end if

       values(sets % index_in(slots, out)) = apply_law(op, &
            & values(sets % index_in(slots, in1)), &
            & values(sets % index_in(slots, in2)))
       available(sets % index_in(slots, out)) = .true.
    end do

  end subroutine primal_execution

  !===================================================================!
  ! The tangent action: Jv by forward traversal of the SAME derived
  ! order. At each operation the local linearization is evaluated at
  ! the stored primal inputs and applied to the incoming tangents.
  ! A tangent may not be read before it is established.
  !===================================================================!

  subroutine tangent_action(flow, slots, sets, indep, seed_values, &
       & computed, order, primal, dot, dot_available)

    class(relation)  , intent(in)  :: flow
    type(set_graph), intent(in)  :: slots
    type(set_map)  , intent(in) :: sets
    type(set_graph), intent(in)  :: indep
    real(dp)         , intent(in)  :: seed_values(:)   ! on X
    type(set_graph), intent(in)  :: computed
    integer          , intent(in)  :: order(:)
    real(dp)         , intent(in)  :: primal(:)        ! full V workspace
    real(dp)         , intent(out) :: dot(:)
    logical          , intent(out) :: dot_available(:)

    real(dp) :: coeff(2)
    integer  :: i, m, op, in1, in2, out

    if (size(seed_values) .ne. sets % size_of(indep) .or. &
         & size(primal) .ne. sets % size_of(slots) .or. &
         & size(dot) .ne. sets % size_of(slots) .or. &
         & size(dot_available) .ne. sets % size_of(slots)) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    dot_available = .false.
    dot           = 0.0_dp

    do i = 1, sets % size_of(indep)
       m = sets % member_of(indep, i)
       dot(sets % index_in(slots, m)) = &
            & seed_values(sets % index_in(indep, m))
       dot_available(sets % index_in(slots, m)) = .true.
    end do

    do i = 1, size(order)
       op  = order(i)
       in1 = slot_for_port(flow, slots, sets, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, sets, op, PORT_IN2)
       out = slot_for_port(flow, slots, sets, op, PORT_OUT)

       if (.not. dot_available(sets % index_in(slots, in1)) .or. &
            & .not. dot_available(sets % index_in(slots, in2))) then
          error stop 'derivative constitution: a tangent was demanded before its input tangents exist'
       end if
       if (.not. sets % has_in(computed, out)) then
          error stop 'derivative constitution: an operation must produce a computed slot'
       end if

       coeff = local_linearization(op, &
            & primal(sets % index_in(slots, in1)), &
            & primal(sets % index_in(slots, in2)))

       dot(sets % index_in(slots, out)) = &
            & coeff(1) * dot(sets % index_in(slots, in1)) + &
            & coeff(2) * dot(sets % index_in(slots, in2))
       dot_available(sets % index_in(slots, out)) = .true.
    end do

  end subroutine tangent_action

  !===================================================================!
  ! The reverse action: J^T zbar by traversing the SAME derived
  ! order backwards, applying the SAME local linearization with its
  ! ends swapped. The += per input-port incidence is the whole
  ! architecture: contributions landing on one slot accumulate, and
  ! no path is counted. The optional hits array records - for every
  ! slot, generically - how many INPUT-PORT INCIDENCE ACCUMULATION
  ! EVENTS landed on it. An event is a traversal fact, not a
  ! numerical one: it is counted whenever
  !
  !      bar(in_i) += L_i * bar(out)
  !
  ! executes, even where that added value happens to be zero.
  ! Incidence multiplicity is therefore NOT numerical nonzero
  ! multiplicity - the coefficients and the seed decide the latter.
  !
  ! The computed domain is deliberately NOT a dummy here: reverse
  ! traversal reads the primal workspace and the flow's incidences
  ! and needs no schema check of its own. Signature symmetry with
  ! the forward action was never a mathematical necessity.
  !===================================================================!

  subroutine reverse_action(flow, slots, sets, indep, order, &
       & primal, response, seed_values, result_values, hits)

    class(relation)  , intent(in)  :: flow
    type(set_graph), intent(in)  :: slots
    type(set_map)  , intent(in) :: sets
    type(set_graph), intent(in)  :: indep
    integer          , intent(in)  :: order(:)
    real(dp)         , intent(in)  :: primal(:)        ! full V workspace
    type(set_graph), intent(in)  :: response         ! Z
    real(dp)         , intent(in)  :: seed_values(:)   ! on Z
    real(dp)         , intent(out) :: result_values(:) ! on X
    integer, allocatable, intent(out), optional :: hits(:)

    real(dp), allocatable :: bar(:)
    integer, allocatable  :: nhits(:)
    real(dp) :: coeff(2)
    integer  :: i, m, op, in1, in2, out

    allocate(bar(sets % size_of(slots)), nhits(sets % size_of(slots)))
    if (size(seed_values) .ne. sets % size_of(response) .or. &
         & size(primal) .ne. sets % size_of(slots) .or. &
         & size(result_values) .ne. sets % size_of(indep)) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    ! Zero is addition's identity: the lawful start of an
    ! accumulator, not a sentinel.
    bar   = 0.0_dp
    nhits = 0

    do i = 1, sets % size_of(response)
       m = sets % member_of(response, i)
       bar(sets % index_in(slots, m)) = &
            & seed_values(sets % index_in(response, m))
    end do

    do i = size(order), 1, -1
       op  = order(i)
       in1 = slot_for_port(flow, slots, sets, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, sets, op, PORT_IN2)
       out = slot_for_port(flow, slots, sets, op, PORT_OUT)

       coeff = local_linearization(op, &
            & primal(sets % index_in(slots, in1)), &
            & primal(sets % index_in(slots, in2)))

       bar(sets % index_in(slots, in1)) = &
            & bar(sets % index_in(slots, in1)) + &
            & coeff(1) * bar(sets % index_in(slots, out))
       nhits(sets % index_in(slots, in1)) = &
            & nhits(sets % index_in(slots, in1)) + 1

       bar(sets % index_in(slots, in2)) = &
            & bar(sets % index_in(slots, in2)) + &
            & coeff(2) * bar(sets % index_in(slots, out))
       nhits(sets % index_in(slots, in2)) = &
            & nhits(sets % index_in(slots, in2)) + 1
    end do

    do i = 1, sets % size_of(indep)
       m = sets % member_of(indep, i)
       result_values(sets % index_in(indep, m)) = &
            & bar(sets % index_in(slots, m))
    end do

    if (present(hits)) hits = nhits

  end subroutine reverse_action

end module derivative_constitution_fixture
