!=====================================================================!
! The chain rule, to any order.
!
! One differentiable statement, one standing state, and per-slot
! derivative paths; this seat assembles the total derivative of
! the composition S(x(s)) at any degree n, through the statement's
! exact partial actions alone. Each symbolic term of the
! composition is one PATTERN: a nondecreasing positive tuple
! d = [d1, ..., dm] with d1 + ... + dm = n, carrying the
! multinomial count
!
!      c(d) = n! / ( prod_i d_i!  *  prod_j multiplicity_j! )
!
! where multiplicity_j counts repeated equal entries of d. The
! patterns are GENERATED, never tabulated - increasing slot count
! first, lexicographic inside one slot count - and every pattern
! slot expands over all channels providing that derivative seat,
! in ordered tuples, which is exactly the factor a symmetric mixed
! partial carries. The coefficient multiplies OUTSIDE the
! statement's calculus; no derivative tensor is ever stored.
!
! A CHANNEL IS ONE INPUT SLOT'S PATH, by construction: it names
! the input argument it perturbs once, and its seats carry only
! direction fields. The old worry that one channel might mix two
! arguments' derivatives is not refused here - it is
! unrepresentable. Two channels naming the same slot would double
! count, and are refused. An unoccupied seat is an absent
! derivative, read as zero: the terms it would have fed are simply
! not assembled.
!
! Degree 0 is the statement's own value. The degree is data: this
! seat refuses only a negative degree, a multinomial count beyond
! the integer range it stores, and an order past the statement's
! declared calculus.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_chain_rule

  use iso_fortran_env     , only : dp => REAL64, int64
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph       , only : set_graph => graph
  use graph_calculus      , only : differentiable_operation
  use class_graph_field   , only : field

  implicit none

  private
  public :: chain_rule
  public :: chain_channel
  public :: chain_seat

  !===================================================================!
  ! One derivative seat: occupied, it carries x^(k) as a direction
  ! field; unoccupied, it is an absent derivative - zero.
  !===================================================================!

  type :: chain_seat

     logical     :: occupied = .false.
     type(field) :: direction

  end type chain_seat

  !===================================================================!
  ! One input slot's path: which argument of the statement this
  ! channel perturbs, and its derivative seats - seat k holds
  ! x^(k).
  !===================================================================!

  type :: chain_channel

     integer :: slot = 0
     type(chain_seat), allocatable :: derivative(:)

   contains

     procedure :: has_degree => channel_has_degree

  end type chain_channel

  !===================================================================!
  ! The stateless assembler verb.
  !===================================================================!

  type :: chain_rule

   contains

     procedure :: assemble

  end type chain_rule

  !===================================================================!
  ! One symbolic term of the composition: a multinomial count and
  ! the derivative order each slot carries. Private: patterns are
  ! how this seat thinks, not what it promises.
  !===================================================================!

  type :: chain_term_pattern

     integer(int64) :: coefficient = 1_int64
     integer, allocatable :: slot_degree(:)

  end type chain_term_pattern

contains

  !===================================================================!
  ! Does the channel carry x^(order)? True when the seat exists
  ! and is occupied.
  !===================================================================!

  pure function channel_has_degree(this, order) result(has)

    class(chain_channel), intent(in) :: this
    integer             , intent(in) :: order
    logical :: has

    has = .false.
    if (.not. allocated(this % derivative)) return
    if (order < 1 .or. order > size(this % derivative)) return

    has = this % derivative(order) % occupied

  end function channel_has_degree

  !===================================================================!
  ! The one public verb: assemble the total derivative of the
  ! given degree. The degree and the channels are proven lawful,
  ! and only then does any statement call run.
  !===================================================================!

  subroutine assemble(this, statement, input_graph, input_data, degree, &
       & channels, output)

    class(chain_rule)              , intent(in)    :: this
    class(differentiable_operation), intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(field)                    , intent(in)    :: input_data(:)
    integer                        , intent(in)    :: degree
    type(chain_channel)            , intent(in)    :: channels(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(chain_term_pattern), allocatable :: patterns(:)
    real(dp), allocatable :: running(:)
    logical :: started
    integer :: p, ncomp

    associate(unread => this)
    end associate

    if (degree < 0) then
       error stop 'chain_rule: degree is supported'
    end if

    call require_lawful_channels(channels, size(input_data))

    !----------------------------------------------------------------!
    ! Degree 0 is the statement's own value.
    !----------------------------------------------------------------!

    if (degree == 0) then
       call statement % apply(input_graph, input_data, output)
       return
    end if

    call get_patterns(degree, patterns)

    started = .false.
    ncomp   = 1

    do p = 1, size(patterns)
       call assemble_pattern(statement, input_graph, input_data, &
            & patterns(p), channels, running, started, ncomp)
    end do

    call land_output(statement, input_graph, input_data, running, &
         & started, ncomp, output)

  end subroutine assemble

  !===================================================================!
  ! The channel laws: every channel names an input slot the
  ! statement actually takes, and no two channels name the same
  ! slot - a duplicated slot would double count chain-rule terms.
  ! Mixing arguments inside one channel is unrepresentable here.
  !===================================================================!

  pure subroutine require_lawful_channels(channels, nslots)

    type(chain_channel), intent(in) :: channels(:)
    integer            , intent(in) :: nslots

    integer :: i, j

    do i = 1, size(channels)
       if (channels(i) % slot < 1 .or. channels(i) % slot > nslots) then
          error stop 'chain_rule: a channel names an input slot'
       end if
    end do

    do i = 1, size(channels)
       do j = i + 1, size(channels)
          if (channels(i) % slot == channels(j) % slot) then
             error stop 'chain_rule: duplicate slot channel is refused'
          end if
       end do
    end do

  end subroutine require_lawful_channels

  !===================================================================!
  ! The patterns of the composition at one degree, GENERATED:
  ! every nondecreasing positive tuple summing to the degree, in
  ! the stable order - increasing slot count first, lexicographic
  ! inside one slot count - each carrying its multinomial count.
  !===================================================================!

  subroutine get_patterns(degree, patterns)

    integer                              , intent(in)  :: degree
    type(chain_term_pattern), allocatable, intent(out) :: patterns(:)

    integer, allocatable :: tuple(:)
    integer :: slots

    allocate(patterns(0))

    do slots = 1, degree
       allocate(tuple(slots))
       call generate_patterns(degree, tuple, 1, 1, patterns)
       deallocate(tuple)
    end do

  end subroutine get_patterns

  recursive subroutine generate_patterns(remaining, tuple, position, &
       & minimum, patterns)

    integer                              , intent(in)    :: remaining
    integer                              , intent(inout) :: tuple(:)
    integer                              , intent(in)    :: position
    integer                              , intent(in)    :: minimum
    type(chain_term_pattern), allocatable, intent(inout) :: patterns(:)

    integer :: slots_left, entry_degree

    slots_left = size(tuple) - position + 1

    if (slots_left == 1) then
       if (remaining >= minimum) then
          tuple(position) = remaining
          call append_pattern(patterns, tuple)
       end if
       return
    end if

    do entry_degree = minimum, remaining / slots_left
       tuple(position) = entry_degree
       call generate_patterns(remaining - entry_degree, tuple, &
            & position + 1, entry_degree, patterns)
    end do

  end subroutine generate_patterns

  subroutine append_pattern(patterns, tuple)

    type(chain_term_pattern), allocatable, intent(inout) :: patterns(:)
    integer                              , intent(in)    :: tuple(:)

    type(chain_term_pattern), allocatable :: grown(:)
    integer :: n

    n = size(patterns)
    allocate(grown(n + 1))
    grown(1:n) = patterns
    grown(n + 1) % slot_degree = tuple
    grown(n + 1) % coefficient = pattern_coefficient(sum(tuple), tuple)
    call move_alloc(grown, patterns)

  end subroutine append_pattern

  !===================================================================!
  ! The multinomial count of one pattern, as one exact integer
  ! division; the tuple's nondecreasing order makes every
  ! multiplicity a contiguous run.
  !===================================================================!

  pure function pattern_coefficient(total, slot_degree) result(coefficient)

    integer, intent(in) :: total
    integer, intent(in) :: slot_degree(:)
    integer(int64)      :: coefficient

    integer(int64) :: denominator
    integer :: j, run

    denominator = 1_int64
    do j = 1, size(slot_degree)
       denominator = denominator * factorial_int64(slot_degree(j))
    end do

    run = 1
    do j = 2, size(slot_degree)
       if (slot_degree(j) == slot_degree(j - 1)) then
          run = run + 1
       else
          denominator = denominator * factorial_int64(run)
          run = 1
       end if
    end do
    denominator = denominator * factorial_int64(run)

    coefficient = factorial_int64(total) / denominator

  end function pattern_coefficient

  pure function factorial_int64(n) result(factorial)

    integer, intent(in) :: n
    integer(int64)      :: factorial

    integer :: k

    if (n > 20) then
       error stop 'chain_rule: pattern coefficient is representable'
    end if

    factorial = 1_int64
    do k = 2, n
       factorial = factorial * int(k, int64)
    end do

  end function factorial_int64

  !===================================================================!
  ! Expand one pattern over the channels: every pattern slot
  ! ranges over every channel providing that derivative seat, in
  ! ordered tuples, and each admitted tuple becomes one partial
  ! action, scaled by the pattern's count.
  !===================================================================!

  subroutine assemble_pattern(statement, input_graph, input_data, pattern, &
       & channels, running, started, ncomp)

    class(differentiable_operation), intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(field)                    , intent(in)    :: input_data(:)
    type(chain_term_pattern)       , intent(in)    :: pattern
    type(chain_channel)            , intent(in)    :: channels(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: ncomp

    integer :: chosen(size(pattern % slot_degree))
    integer :: k, j, nchannels
    logical :: admitted

    k         = size(pattern % slot_degree)
    nchannels = size(channels)
    if (nchannels == 0) return

    chosen = 1

    do

       admitted = .true.
       do j = 1, k
          if (.not. channels(chosen(j)) % &
               & has_degree(pattern % slot_degree(j))) then
             admitted = .false.
             exit
          end if
       end do

       if (admitted) then
          call emit_term(statement, input_graph, input_data, pattern, &
               & channels, chosen, running, started, ncomp)
       end if

       ! the odometer: advance the last slot, carrying leftwards
       j = k
       do
          chosen(j) = chosen(j) + 1
          if (chosen(j) <= nchannels) exit
          chosen(j) = 1
          j = j - 1
          if (j == 0) return
       end do

    end do

  end subroutine assemble_pattern

  !===================================================================!
  ! One admitted tuple, one partial action: the slot list and the
  ! direction family are written from the chosen seats themselves,
  ! the statement's calculus depth is proven before the call, and
  ! the count scales the accumulated term outside its calculus.
  !===================================================================!

  subroutine emit_term(statement, input_graph, input_data, pattern, &
       & channels, chosen, running, started, ncomp)

    class(differentiable_operation), intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(field)                    , intent(in)    :: input_data(:)
    type(chain_term_pattern)       , intent(in)    :: pattern
    type(chain_channel)            , intent(in)    :: channels(:)
    integer                        , intent(in)    :: chosen(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: ncomp

    class(graph_field), allocatable :: answer
    type(field) :: directions(size(chosen))
    integer     :: slots(size(chosen))
    real(dp), allocatable :: term(:)
    integer :: j, k

    k = size(chosen)

    if (k > statement % max_degree()) then
       error stop 'chain_rule: the statement supports the requested order'
    end if

    do j = 1, k
       slots(j)      = channels(chosen(j)) % slot
       directions(j) = channels(chosen(j)) % &
            & derivative(pattern % slot_degree(j)) % direction
    end do

    call statement % partial_action(input_graph, input_data, slots, &
         & directions, answer)

    call answer % get_real_vector(term)

    if (started) then
       if (size(term) /= size(running)) then
          error stop 'chain_rule: accumulated terms share one shape'
       end if
       running = running + real(pattern % coefficient, dp) * term
    else
       running = real(pattern % coefficient, dp) * term
       ncomp   = answer % num_components()
       started = .true.
    end if

  end subroutine emit_term

  !===================================================================!
  ! Land the total on the statement's own domain. When no term was
  ! admitted - every seat the degree asked for was absent - the
  ! total is the codomain's zero, its shape learned from the
  ! statement's own value.
  !===================================================================!

  subroutine land_output(statement, input_graph, input_data, running, &
       & started, ncomp, output)

    class(differentiable_operation), intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(field)                    , intent(in)    :: input_data(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(in)    :: started
    integer                        , intent(in)    :: ncomp
    class(graph_field), allocatable, intent(inout) :: output

    class(graph_field), allocatable :: probe
    type(field)     :: total
    type(set_graph) :: on
    integer         :: n_on, width

    call statement % domain(input_graph, on, n_on)

    width = ncomp
    if (.not. started) then
       call statement % apply(input_graph, input_data, probe)
       call probe % get_real_vector(running)
       running = 0.0_dp
       width   = probe % num_components()
    end if

    total = field('total derivative', on, size(running) / width, ncomp=width)
    call total % set_real_vector(running)

    if (allocated(output)) deallocate(output)
    allocate(output, source=total)

  end subroutine land_output

end module class_graph_chain_rule
