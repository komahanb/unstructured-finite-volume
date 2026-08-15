!=====================================================================!
! THE LEARNING CONSTITUTION FIXTURE - test-local, deliberately
! outside src/: one learning example has not earned a model, layer,
! activation or loss in production. This module is the ONE place in
! the learning tower where the symbols mean something:
!
!      OP_PREDICT -> a * b          OP_ERROR -> a - b
!
! and its evaluation semantics are the learning tower's OWN - not
! the calculator's. There, unknown outputs are state and a residual
! reads q(out) - law(inputs). Here the outputs are COMPUTED: the
! law executes INTO the out slot, availability spreads along the
! derived order, and the residual is simply the VALUE at the slot
! located by L. There is no independent q(e); before evaluation the
! computed domain owns no numbers, and no zero, NaN or sentinel
! ever stands for "not yet" - absence is a separate flag.
!
! Everything else is structural delegation, parameterized by the
! abstract contracts - class(relation), class(set) - handed
! in by the caller: no learning singleton, no stored-type demands.
! The slot on a port and the located slot of a row are DISCOVERED
! by uniqueness scans and refused otherwise; no slot name is wired
! to any operation anywhere. Structure chooses the operands; the
! law sees a symbol and two numbers. No derivative, no chain rule,
! no gradient: forward laws only.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module learning_constitution_fixture

  use iso_fortran_env, only : dp => REAL64
  use learning_assert, only : OP_PREDICT, OP_ERROR
  use learning_assert, only : PORT_IN1, PORT_IN2, PORT_OUT
  use graph_set  , only : set
  use graph_relation , only : relation

  implicit none

  private
  public :: apply_law, slot_for_port, located_slot, generated_residual

contains

  !===================================================================!
  ! The law table: the exact line where meaning enters the tower.
  ! An unbound symbol refuses - there is no "if not predict then
  ! error" anywhere.
  !===================================================================!

  real(dp) function apply_law(op, a, b) result(c)

    integer , intent(in) :: op
    real(dp), intent(in) :: a, b

    select case (op)
    case (OP_PREDICT)
       c = a * b
    case (OP_ERROR)
       c = a - b
    case default
       error stop 'constitution: no law binds this operation symbol'
    end select

  end function apply_law

  !===================================================================!
  ! The unique slot standing on a port of an operation - found by
  ! scanning the value slots against R_flow, refused at zero or
  ! many. The facts live in the relation, never in the evaluator.
  !===================================================================!

  integer function slot_for_port(flow, slots, op, port) result(found)

    class(relation)  , intent(in) :: flow
    class(set), intent(in) :: slots
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
       error stop 'constitution: a port names exactly one slot'
    end if

  end function slot_for_port

  !===================================================================!
  ! The unique home of a residual row, read FROM the location
  ! relation by scanning the value slots - exactly one, or refusal.
  !===================================================================!

  integer function located_slot(located, slots, row) result(home)

    class(relation)  , intent(in) :: located
    class(set), intent(in) :: slots
    integer          , intent(in) :: row

    integer :: j, hits

    hits = 0
    home = 0
    do j = 1, slots % size()
       if (located % has([row, slots % member(j)])) then
          hits = hits + 1
          home = slots % member(j)
       end if
    end do
    if (hits .ne. 1) then
       error stop 'constitution: a residual row lives at exactly one value slot'
    end if

  end function located_slot

  !===================================================================!
  ! The constituted evaluation. Seed the observed and trainable
  ! values by DOMAIN enumeration, execute the caller-DERIVED order -
  ! each operation's slots discovered from R_flow, each input
  ! demanded available, each output required to land in the
  ! computed domain - then read every row's residual as the VALUE
  ! at the home L locates. The caller may ask for the final
  ! workspace (trace) and for the trainable members actually read
  ! as inputs (touched) - the constituted support, to be held
  ! against Level 6's structural one.
  !===================================================================!

  subroutine generated_residual(flow, located, slots, rows, &
       & observed, observed_values, trainable, trainable_values, &
       & computed, order, residual, touched, trace)

    class(relation)  , intent(in)  :: flow, located
    class(set), intent(in)  :: slots        ! V
    class(set), intent(in)  :: rows         ! Y
    class(set), intent(in)  :: observed     ! K
    real(dp)         , intent(in)  :: observed_values(:)
    class(set), intent(in)  :: trainable    ! Theta
    real(dp)         , intent(in)  :: trainable_values(:)
    class(set), intent(in)  :: computed     ! U
    integer          , intent(in)  :: order(:)     ! derived by caller
    real(dp)         , intent(out) :: residual(:)
    integer , allocatable, intent(out), optional :: touched(:)
    real(dp), allocatable, intent(out), optional :: trace(:)

    real(dp) :: values(slots % size())
    logical  :: available(slots % size())
    integer  :: reads(trainable % size())
    integer  :: nreads, i, m, op, in1, in2, out

    if (size(observed_values) .ne. observed % size() .or. &
         & size(trainable_values) .ne. trainable % size() .or. &
         & size(residual) .ne. rows % size()) then
       error stop 'constitution: every vector is sized by its domain'
    end if

    ! Before execution nothing has a value - and no number, zero or
    ! otherwise, is fabricated to say so.
    available = .false.
    values    = 0.0_dp
    nreads    = 0

    ! Seed the observed state through K's own enumeration.
    do i = 1, observed % size()
       m = observed % member(i)
       values(slots % local_index(m)) = &
            & observed_values(observed % local_index(m))
       available(slots % local_index(m)) = .true.
    end do

    ! Seed the trainable state through Theta's own enumeration -
    ! an INPUT to evaluation, never modified here.
    do i = 1, trainable % size()
       m = trainable % member(i)
       values(slots % local_index(m)) = &
            & trainable_values(trainable % local_index(m))
       available(slots % local_index(m)) = .true.
    end do

    ! Execute the derived order: discover, demand, compute, admit.
    do i = 1, size(order)
       op  = order(i)
       in1 = slot_for_port(flow, slots, op, PORT_IN1)
       in2 = slot_for_port(flow, slots, op, PORT_IN2)
       out = slot_for_port(flow, slots, op, PORT_OUT)

       if (.not. available(slots % local_index(in1)) .or. &
            & .not. available(slots % local_index(in2))) then
          error stop 'constitution: an operation was scheduled before its inputs exist'
       end if

       ! The learning schema law: operations produce the slots the
       ! partition classified as computed.
       if (.not. computed % has(out)) then
          error stop 'constitution: an operation must produce a computed slot'
       end if

       call note_read(in1)
       call note_read(in2)

       values(slots % local_index(out)) = apply_law(op, &
            & values(slots % local_index(in1)), &
            & values(slots % local_index(in2)))
       available(slots % local_index(out)) = .true.
    end do

    ! The learning residual rule, complete: the value at the home.
    do i = 1, rows % size()
       m = located_slot(located, slots, rows % member(i))
       if (.not. available(slots % local_index(m))) then
          error stop 'constitution: a residual home was never computed'
       end if
       residual(rows % local_index(rows % member(i))) = &
            & values(slots % local_index(m))
    end do

    if (present(touched)) touched = reads(1:nreads)
    if (present(trace))   trace   = values

  contains

    subroutine note_read(slot)
      integer, intent(in) :: slot
      if (trainable % has(slot)) then
         if (.not. any(reads(1:nreads) .eq. slot)) then
            nreads        = nreads + 1
            reads(nreads) = slot
         end if
      end if
    end subroutine note_read

  end subroutine generated_residual

end module learning_constitution_fixture
