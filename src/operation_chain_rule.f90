!=====================================================================!
! The chain rule to any order: for a composition S(x_1(s), ...,
! x_m(s)) of one differentiable statement with per-argument paths
! x_j(s), assemble the total derivative d^n/ds^n from the
! statement's exact partial actions.
!
! Each term is indexed by an integer partition of n: a
! nondecreasing positive tuple d = [d_1, ..., d_k] with
! d_1 + ... + d_k = n, carrying the multinomial count
!
!      c(d) = n! / ( prod_i d_i!  *  prod_j multiplicity_j! )
!
! where multiplicity_j counts repeated equal entries of d.
! Partitions are generated, not tabulated, in a fixed order:
! increasing entry count first, lexicographic within one entry
! count. Every partition entry ranges over all argument paths
! providing that derivative order, in ordered tuples - the factor
! a symmetric mixed partial requires. The count multiplies the
! result outside the statement's calculus; no derivative tensor is
! stored.
!
! An argument_path names one input slot and carries the derivative
! sequence x^(1), ..., x^(k) of that argument, so one path cannot
! mix two arguments' derivatives. Two paths naming the same slot
! would double count and stop the program. An unoccupied
! derivative is read as zero: the terms it would feed are not
! assembled.
!
! Degree 0 is the statement's value. Inputs that stop the program:
! a negative degree, a multinomial count outside int64, a path
! naming a slot the statement does not take, a duplicated slot,
! and a partition needing a partial past the statement's
! max_degree.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_chain_rule

  use iso_fortran_env     , only : dp => REAL64, int64
  use view_directed , only : directed_graph
  use field_calculus, only : field
  use graph_fractal       , only : graph
  use operation_action    , only : operation
  use field_stored   , only : stored_field

  implicit none

  private
  public :: chain_rule
  public :: argument_path
  public :: path_derivative

  !===================================================================!
  ! One derivative of a path: occupied, it carries x^(k) as a
  ! direction field; unoccupied, it is read as zero.
  !===================================================================!

  type :: path_derivative

     logical     :: occupied = .false.
     type(stored_field) :: direction

  end type path_derivative

  !===================================================================!
  ! One argument's path: the input slot it perturbs, and its
  ! derivative sequence - derivative(k) holds x^(k).
  !===================================================================!

  type :: argument_path

     integer :: slot = 0
     type(path_derivative), allocatable :: derivative(:)

   contains

     procedure :: has_degree => path_has_degree

  end type argument_path

  !===================================================================!
  ! The stateless assembler.
  !===================================================================!

  type :: chain_rule

   contains

     procedure :: assemble

  end type chain_rule

  !===================================================================!
  ! One integer partition of the degree, with its multinomial
  ! count. Private to this module.
  !===================================================================!

  type :: derivative_partition

     integer(int64) :: coefficient = 1_int64
     integer, allocatable :: slot_degree(:)

  end type derivative_partition

contains

  !===================================================================!
  ! True when derivative(order) exists and is occupied.
  !===================================================================!

  pure function path_has_degree(this, order) result(has)

    class(argument_path), intent(in) :: this
    integer             , intent(in) :: order
    logical :: has

    has = .false.
    if (.not. allocated(this % derivative)) return
    if (order < 1 .or. order > size(this % derivative)) return

    has = this % derivative(order) % occupied

  end function path_has_degree

  !===================================================================!
  ! Assemble the total derivative of the given degree. The degree
  ! and the paths are checked before any statement call runs.
  !===================================================================!

  subroutine assemble(this, statement, input_graph, input_data, degree, &
       & paths, output)

    class(chain_rule)              , intent(in)    :: this
    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)                    , intent(in)    :: input_data(:)
    integer                        , intent(in)    :: degree
    type(argument_path)            , intent(in)    :: paths(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_partition), allocatable :: partitions(:)
    real(dp), allocatable :: running(:)
    logical :: started
    integer :: p, num_components

    associate(unread => this)
    end associate

    if (degree < 0) then
       error stop 'chain_rule: degree is supported'
    end if

    call require_valid_paths(paths, size(input_data))

    ! degree 0 is the statement's own value
    if (degree == 0) then
       call statement % apply(input_graph, input_data, output)
       return
    end if

    call enumerate_partitions(degree, partitions)

    started = .false.
    num_components   = 1

    do p = 1, size(partitions)
       call assemble_partition(statement, input_graph, input_data, &
            & partitions(p), paths, running, started, num_components)
    end do

    call write_output(statement, input_graph, input_data, running, &
         & started, num_components, output)

  end subroutine assemble

  !===================================================================!
  ! Check the paths: each must name an input slot the statement
  ! takes, and no two may name the same slot, because a duplicated
  ! slot would double count chain-rule terms. Both violations stop
  ! the program.
  !===================================================================!

  pure subroutine require_valid_paths(paths, nslots)

    type(argument_path), intent(in) :: paths(:)
    integer            , intent(in) :: nslots

    integer :: i, j

    do i = 1, size(paths)
       if (paths(i) % slot < 1 .or. paths(i) % slot > nslots) then
          error stop 'chain_rule: a path names an input slot'
       end if
    end do

    do i = 1, size(paths)
       do j = i + 1, size(paths)
          if (paths(i) % slot == paths(j) % slot) then
             error stop 'chain_rule: duplicate slot path is refused'
          end if
       end do
    end do

  end subroutine require_valid_paths

  !===================================================================!
  ! Generate the partitions of one degree: every nondecreasing
  ! positive tuple summing to the degree, ordered by entry count
  ! then lexicographically, each with its multinomial count.
  !===================================================================!

  subroutine enumerate_partitions(degree, partitions)

    integer                                , intent(in)  :: degree
    type(derivative_partition), allocatable, intent(out) :: partitions(:)

    integer, allocatable :: tuple(:)
    integer :: entries

    allocate(partitions(0))

    do entries = 1, degree
       allocate(tuple(entries))
       call generate_partitions(degree, tuple, 1, 1, partitions)
       deallocate(tuple)
    end do

  end subroutine enumerate_partitions

  recursive subroutine generate_partitions(remaining, tuple, position, &
       & minimum, partitions)

    integer                                , intent(in)    :: remaining
    integer                                , intent(inout) :: tuple(:)
    integer                                , intent(in)    :: position
    integer                                , intent(in)    :: minimum
    type(derivative_partition), allocatable, intent(inout) :: partitions(:)

    integer :: entries_left, entry_degree

    entries_left = size(tuple) - position + 1

    if (entries_left == 1) then
       if (remaining >= minimum) then
          tuple(position) = remaining
          call append_partition(partitions, tuple)
       end if
       return
    end if

    do entry_degree = minimum, remaining / entries_left
       tuple(position) = entry_degree
       call generate_partitions(remaining - entry_degree, tuple, &
            & position + 1, entry_degree, partitions)
    end do

  end subroutine generate_partitions

  subroutine append_partition(partitions, tuple)

    type(derivative_partition), allocatable, intent(inout) :: partitions(:)
    integer                                , intent(in)    :: tuple(:)

    type(derivative_partition), allocatable :: grown(:)
    integer :: n

    n = size(partitions)
    allocate(grown(n + 1))
    grown(1:n) = partitions
    grown(n + 1) % slot_degree = tuple
    grown(n + 1) % coefficient = partition_coefficient(sum(tuple), tuple)
    call move_alloc(grown, partitions)

  end subroutine append_partition

  !===================================================================!
  ! The multinomial count of one partition, as one exact integer
  ! division; the tuple's nondecreasing order makes every
  ! multiplicity a contiguous run.
  !===================================================================!

  pure function partition_coefficient(total, slot_degree) result(coefficient)

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

  end function partition_coefficient

  pure function factorial_int64(n) result(factorial)

    integer, intent(in) :: n
    integer(int64)      :: factorial

    integer :: k

    ! 21! overflows int64; stop rather than wrap
    if (n > 20) then
       error stop 'chain_rule: partition coefficient is representable'
    end if

    factorial = 1_int64
    do k = 2, n
       factorial = factorial * int(k, int64)
    end do

  end function factorial_int64

  !===================================================================!
  ! Expand one partition over the paths: every partition entry
  ! ranges over every path providing that derivative order, in
  ! ordered tuples; each admitted tuple becomes one partial action
  ! scaled by the partition's count.
  !===================================================================!

  subroutine assemble_partition(statement, input_graph, input_data, &
       & partition, paths, running, started, num_components)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)                    , intent(in)    :: input_data(:)
    type(derivative_partition)     , intent(in)    :: partition
    type(argument_path)            , intent(in)    :: paths(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: num_components

    integer :: chosen(size(partition % slot_degree))
    integer :: k, j, npaths
    logical :: admitted

    k      = size(partition % slot_degree)
    npaths = size(paths)
    if (npaths == 0) return

    chosen = 1

    do

       admitted = .true.
       do j = 1, k
          if (.not. paths(chosen(j)) % &
               & has_degree(partition % slot_degree(j))) then
             admitted = .false.
             exit
          end if
       end do

       if (admitted) then
          call emit_term(statement, input_graph, input_data, partition, &
               & paths, chosen, running, started, num_components)
       end if

       ! the odometer: advance the last entry, carrying leftwards
       j = k
       do
          chosen(j) = chosen(j) + 1
          if (chosen(j) <= npaths) exit
          chosen(j) = 1
          j = j - 1
          if (j == 0) return
       end do

    end do

  end subroutine assemble_partition

  !===================================================================!
  ! One admitted tuple, one partial action. The order (the tuple
  ! length) must not exceed the statement's max_degree; violation
  ! stops the program, because the statement declared no partials
  ! past that order. The multinomial count scales the term after
  ! the statement call, and every accumulated term must have one
  ! shape.
  !===================================================================!

  subroutine emit_term(statement, input_graph, input_data, partition, &
       & paths, chosen, running, started, num_components)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)                    , intent(in)    :: input_data(:)
    type(derivative_partition)     , intent(in)    :: partition
    type(argument_path)            , intent(in)    :: paths(:)
    integer                        , intent(in)    :: chosen(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: num_components

    class(field), allocatable :: output
    type(stored_field) :: directions(size(chosen))
    integer     :: slots(size(chosen))
    real(dp), allocatable :: term(:)
    integer :: j, k

    k = size(chosen)

    if (k > statement % max_degree()) then
       error stop 'chain_rule: the statement supports the requested order'
    end if

    do j = 1, k
       slots(j)      = paths(chosen(j)) % slot
       directions(j) = paths(chosen(j)) % &
            & derivative(partition % slot_degree(j)) % direction
    end do

    call statement % partial_action(input_graph, input_data, slots, &
         & directions, output)

    call output % real_vector(term)

    if (started) then
       if (size(term) /= size(running)) then
          error stop 'chain_rule: accumulated terms share one shape'
       end if
       running = running + real(partition % coefficient, dp) * term
    else
       running = real(partition % coefficient, dp) * term
       num_components   = output % num_components()
       started = .true.
    end if

  end subroutine emit_term

  !===================================================================!
  ! Write the total onto the statement's domain. When no term was
  ! admitted (every required derivative was unoccupied) the total
  ! is zero, with its shape taken from the statement's own value.
  !===================================================================!

  subroutine write_output(statement, input_graph, input_data, running, &
       & started, num_components, output)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)                    , intent(in)    :: input_data(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(in)    :: started
    integer                        , intent(in)    :: num_components
    class(field), allocatable, intent(inout) :: output

    class(field), allocatable :: value
    type(stored_field)     :: total
    type(graph) :: on
    integer         :: n_on, width

    call statement % domain(input_graph, on, n_on)

    width = num_components
    if (.not. started) then
       call statement % apply(input_graph, input_data, value)
       call value % real_vector(running)
       running = 0.0_dp
       width   = value % num_components()
    end if

    total = stored_field('total derivative', on, size(running) / width, num_components=width)
    call total % set_real_vector(running)

    if (allocated(output)) deallocate(output)
    allocate(output, source=total)

  end subroutine write_output

end module operation_chain_rule
