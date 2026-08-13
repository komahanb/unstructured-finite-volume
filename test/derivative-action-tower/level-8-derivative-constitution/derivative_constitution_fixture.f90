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
! class(relation), class(member_set) - handed in by the caller: no
! derivative singleton, no stored-type demands, no graph argument.
! What the evaluators consume is R_flow incidences, the DERIVED
! order, primal values, and the law table: the structural J-pattern
! is never an input here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module derivative_constitution_fixture

  use iso_fortran_env, only : dp => REAL64
  use derivative_assert, only : OP_PRODUCT, OP_SUM
  use derivative_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_carrier    , only : member_set
  use graph_relation   , only : relation

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
  ! by scanning the value slots against R_flow, refused at zero or
  ! many. The facts live in the relation, never in the evaluators.
  !===================================================================!

  integer function slot_for_port(flow, slots, op, port) result(found)

    class(relation)  , intent(in) :: flow
    class(member_set), intent(in) :: slots
    integer          , intent(in) :: op, port

    integer :: j, hits

    hits  = 0
    found = 0
    do j = 1, slots % size()
       if (flow % has([op, slots % member(j), port])) then
          hits  = hits + 1
          found = slots % member(j)
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

  subroutine primal_execution(flow, slots, indep, indep_values, &
       & computed, order, values, available)

    class(relation)  , intent(in)  :: flow
    class(member_set), intent(in)  :: slots      ! V
    class(member_set), intent(in)  :: indep      ! X
    real(dp)         , intent(in)  :: indep_values(:)
    class(member_set), intent(in)  :: computed   ! C
    integer          , intent(in)  :: order(:)   ! derived by caller
    real(dp)         , intent(out) :: values(:)
    logical          , intent(out) :: available(:)

    integer :: i, m, op, in1, in2, out

    if (size(indep_values) .ne. indep % size() .or. &
         & size(values) .ne. slots % size() .or. &
         & size(available) .ne. slots % size()) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    available = .false.
    values    = 0.0_dp

    do i = 1, indep % size()
       m = indep % member(i)
       values(slots % local_index(m)) = &
            & indep_values(indep % local_index(m))
       available(slots % local_index(m)) = .true.
    end do

    do i = 1, size(order)
       op  = order(i)
       in1 = slot_for_port(flow, slots, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, op, PORT_IN2)
       out = slot_for_port(flow, slots, op, PORT_OUT)

       if (.not. available(slots % local_index(in1)) .or. &
            & .not. available(slots % local_index(in2))) then
          error stop 'derivative constitution: an operation was scheduled before its primal inputs exist'
       end if
       if (.not. computed % has(out)) then
          error stop 'derivative constitution: an operation must produce a computed slot'
       end if

       values(slots % local_index(out)) = apply_law(op, &
            & values(slots % local_index(in1)), &
            & values(slots % local_index(in2)))
       available(slots % local_index(out)) = .true.
    end do

  end subroutine primal_execution

  !===================================================================!
  ! The tangent action: Jv by forward traversal of the SAME derived
  ! order. At each operation the local linearization is evaluated at
  ! the stored primal inputs and applied to the incoming tangents.
  ! A tangent may not be read before it is established.
  !===================================================================!

  subroutine tangent_action(flow, slots, indep, seed_values, &
       & computed, order, primal, dot, dot_available)

    class(relation)  , intent(in)  :: flow
    class(member_set), intent(in)  :: slots
    class(member_set), intent(in)  :: indep
    real(dp)         , intent(in)  :: seed_values(:)   ! on X
    class(member_set), intent(in)  :: computed
    integer          , intent(in)  :: order(:)
    real(dp)         , intent(in)  :: primal(:)        ! full V workspace
    real(dp)         , intent(out) :: dot(:)
    logical          , intent(out) :: dot_available(:)

    real(dp) :: coeff(2)
    integer  :: i, m, op, in1, in2, out

    if (size(seed_values) .ne. indep % size() .or. &
         & size(primal) .ne. slots % size() .or. &
         & size(dot) .ne. slots % size() .or. &
         & size(dot_available) .ne. slots % size()) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    dot_available = .false.
    dot           = 0.0_dp

    do i = 1, indep % size()
       m = indep % member(i)
       dot(slots % local_index(m)) = &
            & seed_values(indep % local_index(m))
       dot_available(slots % local_index(m)) = .true.
    end do

    do i = 1, size(order)
       op  = order(i)
       in1 = slot_for_port(flow, slots, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, op, PORT_IN2)
       out = slot_for_port(flow, slots, op, PORT_OUT)

       if (.not. dot_available(slots % local_index(in1)) .or. &
            & .not. dot_available(slots % local_index(in2))) then
          error stop 'derivative constitution: a tangent was demanded before its input tangents exist'
       end if
       if (.not. computed % has(out)) then
          error stop 'derivative constitution: an operation must produce a computed slot'
       end if

       coeff = local_linearization(op, &
            & primal(slots % local_index(in1)), &
            & primal(slots % local_index(in2)))

       dot(slots % local_index(out)) = &
            & coeff(1) * dot(slots % local_index(in1)) + &
            & coeff(2) * dot(slots % local_index(in2))
       dot_available(slots % local_index(out)) = .true.
    end do

  end subroutine tangent_action

  !===================================================================!
  ! The reverse action: J^T zbar by traversing the SAME derived
  ! order backwards, applying the SAME local linearization with its
  ! ends swapped. The += per input-port incidence is the whole
  ! architecture: contributions landing on one slot accumulate, and
  ! no path is counted. The optional hits array records - for every
  ! slot, generically - how many incidences accumulated into it.
  !===================================================================!

  subroutine reverse_action(flow, slots, indep, computed, order, &
       & primal, response, seed_values, result_values, hits)

    class(relation)  , intent(in)  :: flow
    class(member_set), intent(in)  :: slots
    class(member_set), intent(in)  :: indep
    class(member_set), intent(in)  :: computed
    integer          , intent(in)  :: order(:)
    real(dp)         , intent(in)  :: primal(:)        ! full V workspace
    class(member_set), intent(in)  :: response         ! Z
    real(dp)         , intent(in)  :: seed_values(:)   ! on Z
    real(dp)         , intent(out) :: result_values(:) ! on X
    integer, allocatable, intent(out), optional :: hits(:)

    real(dp) :: bar(slots % size())
    integer  :: nhits(slots % size())
    real(dp) :: coeff(2)
    integer  :: i, m, op, in1, in2, out

    if (size(seed_values) .ne. response % size() .or. &
         & size(primal) .ne. slots % size() .or. &
         & size(result_values) .ne. indep % size()) then
       error stop 'derivative constitution: every vector is sized by its domain'
    end if

    ! Zero is addition's identity: the lawful start of an
    ! accumulator, not a sentinel.
    bar   = 0.0_dp
    nhits = 0

    do i = 1, response % size()
       m = response % member(i)
       bar(slots % local_index(m)) = &
            & seed_values(response % local_index(m))
    end do

    do i = size(order), 1, -1
       op  = order(i)
       in1 = slot_for_port(flow, slots, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, op, PORT_IN2)
       out = slot_for_port(flow, slots, op, PORT_OUT)

       coeff = local_linearization(op, &
            & primal(slots % local_index(in1)), &
            & primal(slots % local_index(in2)))

       bar(slots % local_index(in1)) = &
            & bar(slots % local_index(in1)) + &
            & coeff(1) * bar(slots % local_index(out))
       nhits(slots % local_index(in1)) = &
            & nhits(slots % local_index(in1)) + 1

       bar(slots % local_index(in2)) = &
            & bar(slots % local_index(in2)) + &
            & coeff(2) * bar(slots % local_index(out))
       nhits(slots % local_index(in2)) = &
            & nhits(slots % local_index(in2)) + 1
    end do

    do i = 1, indep % size()
       m = indep % member(i)
       result_values(indep % local_index(m)) = &
            & bar(slots % local_index(m))
    end do

    if (present(hits)) hits = nhits

  end subroutine reverse_action

end module derivative_constitution_fixture
