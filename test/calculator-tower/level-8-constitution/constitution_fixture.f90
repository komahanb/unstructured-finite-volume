!=====================================================================!
! THE ARITHMETIC CONSTITUTION FIXTURE - test-local, deliberately
! outside src/: one calculator has not earned a universal
! constitution interface. This module is the ONE place in the tower
! where the symbols mean something:
!
!      OP_PLUS  -> x + y          OP_TIMES -> x * y
!
! Everything else here is structural delegation, parameterized by
! explicit abstract contracts - class(relation), class(set) -
! handed in by the caller: no global calculator singleton, no
! stored-type demands. The producer of a slot, the slot on a port,
! the located slot of a row are DISCOVERED by uniqueness scans and
! refused otherwise; no "if not plus then times" lives anywhere.
! Level 8 proves this machinery; Level 9 composes it with the
! solver.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module arithmetic_constitution_fixture

  use iso_fortran_env  , only : dp => REAL64
  use calculator_assert, only : OP_PLUS, OP_TIMES
  use calculator_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_set    , only : set
  use graph_relation   , only : relation

  implicit none

  private
  public :: apply_law, generated_residual, constitution_support

contains

  !===================================================================!
  ! The law table: the exact line where meaning enters the tower.
  !===================================================================!

  real(dp) function apply_law(op, x, y) result(z)

    integer , intent(in) :: op
    real(dp), intent(in) :: x, y

    select case (op)
    case (OP_PLUS)
       z = x + y
    case (OP_TIMES)
       z = x * y
    case default
       error stop 'constitution: no law binds this operation symbol'
    end select

  end function apply_law

  !===================================================================!
  ! The unique operation whose OUT port names this slot - found by
  ! the structural law alone, refused at zero or many.
  !===================================================================!

  integer function producer_of(flow, ops, x_out) result(op)

    class(relation)  , intent(in) :: flow
    class(set), intent(in) :: ops
    integer          , intent(in) :: x_out

    integer :: j, n

    n  = 0
    op = 0
    do j = 1, ops % size()
       if (flow % has([ops % member(j), x_out, PORT_OUT])) then
          n  = n + 1
          op = ops % member(j)
       end if
    end do
    if (n /= 1) error stop 'constitution: one producer per slot - or refusal'

  end function producer_of

  !===================================================================!
  ! The unique slot the flow holds for (operation, port).
  !===================================================================!

  integer function slot_for_port(flow, slots, op, port) result(slot)

    class(relation)  , intent(in) :: flow
    class(set), intent(in) :: slots
    integer          , intent(in) :: op, port

    integer :: i, n

    n    = 0
    slot = 0
    do i = 1, slots % size()
       if (flow % has([op, slots % member(i), port])) then
          n    = n + 1
          slot = slots % member(i)
       end if
    end do
    if (n /= 1) error stop 'constitution: one port, one slot - or refusal'

  end function slot_for_port

  !===================================================================!
  ! The slot where a residual row is located.
  !===================================================================!

  integer function located_slot(located, slots, row) result(x_out)

    class(relation)  , intent(in) :: located
    class(set), intent(in) :: slots
    integer          , intent(in) :: row

    integer :: j, n

    n     = 0
    x_out = 0
    do j = 1, slots % size()
       if (located % has([row, slots % member(j)])) then
          n     = n + 1
          x_out = slots % member(j)
       end if
    end do
    if (n /= 1) error stop 'constitution: one location per residual row - or refusal'

  end function located_slot

  !===================================================================!
  ! The generated residual: structure chooses every slot, the law
  ! supplies every number. Values are answered by the two declared
  ! domains - known K with its vector, unknown U with the state -
  ! and a slot with no home is refused.
  !===================================================================!

  subroutine generated_residual(flow, located, slots, ops, rows, &
       &                        known, known_values, unknown, ustate, r)

    class(relation)  , intent(in)      :: flow, located
    class(set), intent(in)      :: slots, ops, rows
    class(set), intent(in)      :: known, unknown
    real(dp), intent(in)               :: known_values(:), ustate(:)
    real(dp), allocatable, intent(out) :: r(:)

    integer  :: i, row, x_out, op, in1, in2
    real(dp) :: produced

    allocate(r(rows % size()))
    do i = 1, rows % size()
       row   = rows % member(i)
       x_out = located_slot(located, slots, row)
       op    = producer_of(flow, ops, x_out)
       in1   = slot_for_port(flow, slots, op, PORT_IN1)
       in2   = slot_for_port(flow, slots, op, PORT_IN2)

       produced = apply_law(op, q_of(in1), q_of(in2))
       r(rows % local_index(row)) = q_of(x_out) - produced
    end do

  contains

    real(dp) function q_of(slot) result(v)
      integer, intent(in) :: slot
      if (known % has(slot)) then
         v = known_values(known % local_index(slot))
      else if (unknown % has(slot)) then
         v = ustate(unknown % local_index(slot))
      else
         error stop 'constitution: a slot with no home holds no value'
      end if
    end function q_of

  end subroutine generated_residual

  !===================================================================!
  ! The slots one row's constituted equation touches - out, in1,
  ! in2 - scanned canonically against the slot domain, by the SAME
  ! structural interpretation the evaluator uses.
  !===================================================================!

  subroutine constitution_support(flow, located, slots, ops, row, members)

    class(relation)  , intent(in)      :: flow, located
    class(set), intent(in)      :: slots, ops
    integer          , intent(in)      :: row
    integer, allocatable, intent(out)  :: members(:)

    integer :: x_out, op, in1, in2, i, n, m
    integer :: keep(slots % size())

    x_out = located_slot(located, slots, row)
    op    = producer_of(flow, ops, x_out)
    in1   = slot_for_port(flow, slots, op, PORT_IN1)
    in2   = slot_for_port(flow, slots, op, PORT_IN2)

    n = 0
    do i = 1, slots % size()
       m = slots % member(i)
       if (m == x_out .or. m == in1 .or. m == in2) then
          n = n + 1
          keep(n) = m
       end if
    end do
    members = keep(1:n)

  end subroutine constitution_support

end module arithmetic_constitution_fixture
