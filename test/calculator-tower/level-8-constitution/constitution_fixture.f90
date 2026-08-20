!=====================================================================!
! THE ARITHMETIC CONSTITUTION FIXTURE - test-local, deliberately
! outside src/: one calculator has not earned a universal
! constitution interface. This module is the ONE place in the tower
! where the symbols mean something:
!
!      OP_PLUS  -> x + y          OP_TIMES -> x * y
!
! Everything else here is structural delegation, parameterized by
! explicit abstract contracts - class(relation), and set identities -
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
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use relation_finitary   , only : relation

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

  integer function producer_of(flow, ops, sets, x_out) result(op)

    class(relation)  , intent(in) :: flow
    type(graph), intent(in) :: ops
    type(set_map)  , intent(in) :: sets
    integer        , intent(in) :: x_out

    integer :: j, n

    n  = 0
    op = 0
    do j = 1, sets % num_members_of(ops)
       if (flow % has([sets % member_of(ops, j), x_out, PORT_OUT])) then
          n  = n + 1
          op = sets % member_of(ops, j)
       end if
    end do
    if (n /= 1) error stop 'constitution: one producer per slot - or refusal'

  end function producer_of

  !===================================================================!
  ! The unique slot the flow holds for (operation, port).
  !===================================================================!

  integer function slot_for_port(flow, slots, sets, op, port) result(slot)

    class(relation)  , intent(in) :: flow
    type(graph), intent(in) :: slots
    type(set_map)  , intent(in) :: sets
    integer        , intent(in) :: op, port

    integer :: i, n

    n    = 0
    slot = 0
    do i = 1, sets % num_members_of(slots)
       if (flow % has([op, sets % member_of(slots, i), port])) then
          n    = n + 1
          slot = sets % member_of(slots, i)
       end if
    end do
    if (n /= 1) error stop 'constitution: one port, one slot - or refusal'

  end function slot_for_port

  !===================================================================!
  ! The slot where a residual row is located.
  !===================================================================!

  integer function located_slot(located, slots, sets, row) result(x_out)

    class(relation)  , intent(in) :: located
    type(graph), intent(in) :: slots
    type(set_map)  , intent(in) :: sets
    integer        , intent(in) :: row

    integer :: j, n

    n     = 0
    x_out = 0
    do j = 1, sets % num_members_of(slots)
       if (located % has([row, sets % member_of(slots, j)])) then
          n     = n + 1
          x_out = sets % member_of(slots, j)
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

  subroutine generated_residual(flow, located, slots, ops, rows, sets, &
       &                        known, known_values, unknown, ustate, r)

    class(relation)  , intent(in)      :: flow, located
    type(graph), intent(in)      :: slots, ops, rows
    type(set_map)  , intent(in)      :: sets
    type(graph), intent(in)      :: known, unknown
    real(dp), intent(in)               :: known_values(:), ustate(:)
    real(dp), allocatable, intent(out) :: r(:)

    integer  :: i, row, x_out, op, in1, in2
    real(dp) :: produced

    allocate(r(sets % num_members_of(rows)))
    do i = 1, sets % num_members_of(rows)
       row   = sets % member_of(rows, i)
       x_out = located_slot(located, slots, sets, row)
       op    = producer_of(flow, ops, sets, x_out)
       in1   = slot_for_port(flow, slots, sets, op, PORT_IN1)
       in2   = slot_for_port(flow, slots, sets, op, PORT_IN2)

       produced = apply_law(op, q_of(in1), q_of(in2))
       r(sets % index_in(rows, row)) = q_of(x_out) - produced
    end do

  contains

    real(dp) function q_of(slot) result(v)
      integer, intent(in) :: slot
      ! sets, known and unknown come from the host by association.
      if (sets % has(known, slot)) then
         v = known_values(sets % index_in(known, slot))
      else if (sets % has(unknown, slot)) then
         v = ustate(sets % index_in(unknown, slot))
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

  subroutine constitution_support(flow, located, slots, ops, sets, row, members)

    class(relation)  , intent(in)      :: flow, located
    type(graph)  , intent(in)      :: slots, ops
    type(set_map)    , intent(in)      :: sets
    integer          , intent(in)      :: row
    integer, allocatable, intent(out)  :: members(:)

    integer              :: x_out, op, in1, in2, i, n, m
    integer, allocatable :: keep(:)

    allocate(keep(sets % num_members_of(slots)))

    x_out = located_slot(located, slots, sets, row)
    op    = producer_of(flow, ops, sets, x_out)
    in1   = slot_for_port(flow, slots, sets, op, PORT_IN1)
    in2   = slot_for_port(flow, slots, sets, op, PORT_IN2)

    n = 0
    do i = 1, sets % num_members_of(slots)
       m = sets % member_of(slots, i)
       if (m == x_out .or. m == in1 .or. m == in2) then
          n = n + 1
          keep(n) = m
       end if
    end do
    members = keep(1:n)

  end subroutine constitution_support

end module arithmetic_constitution_fixture
