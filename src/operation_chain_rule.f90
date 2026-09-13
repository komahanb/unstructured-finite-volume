!=====================================================================!
! The chain rule to any order: for a composition S(x_1(s), ...,
! x_m(s)) of one differentiable statement with per-argument paths
! x_j(s), assemble the total derivative d^n/ds^n from the
! statement's exact partial actions.
!
! Each term is indexed by an integer partition of n: a
! nondecreasing positive tuple d = [d_1, ..., d_k] with
! d_1 + ... + d_k = n, with the multinomial count
!
!      c(d) = n! / ( prod_i d_i!  *  prod_j multiplicity_j! )
!
! where multiplicity_j counts repeated equal entries of d.
! The partitions of every degree up to the assembler's highest are
! enumerated once, when it is constructed, in a fixed order:
! increasing entry count first, lexicographic within one entry
! count. They are the index structure of the chain rule; an
! algorithm supplies only the paths, so Newton's corrections and a
! march's tangents each read one table. Every partition entry
! ranges over all argument paths providing that derivative order,
! in ordered tuples - the factor
! a symmetric mixed partial requires. The count multiplies the
! result outside the statement's calculus; no derivative tensor is
! stored.
!
! An argument_path names one input slot and stores the derivative
! sequence x^(1), ..., x^(k) of that argument, so one path cannot
! mix two arguments' derivatives. Two paths naming the same slot
! would double count and stop the program. An unoccupied
! derivative is read as zero: the terms it would enter are not
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

  use iso_fortran_env, only : int64
  use util_precision  , only : dp
  use view_directed , only : directed_graph
  use field_calculus, only : field
  use graph_fractal       , only : graph
  use operation_action    , only : operation, argument, variation
  use operation_action    , only : binding, emit_real
  use field_stored   , only : stored_field

  implicit none

  private
  public :: total_derivative, derivative_of
  public :: argument_path
  public :: path_derivative

  !===================================================================!
  ! One derivative of a path: occupied, it stores x^(k) as a
  ! direction field; unoccupied, it is read as zero.
  !===================================================================!

  type :: path_derivative

     logical     :: occupied = .false.
     type(stored_field) :: direction

  end type path_derivative

  !===================================================================!
  ! One argument's path: the argument of the statement it perturbs,
  ! obtained from the statement itself, and its derivative sequence
  ! - derivative(k) stores x^(k).
  !===================================================================!

  type :: argument_path

     type(argument) :: wrt
     type(path_derivative), allocatable :: derivative(:)

   contains

     procedure :: has_degree => path_has_degree

  end type argument_path

  interface argument_path
     module procedure create_argument_path
  end interface argument_path

  !===================================================================!
  ! One integer partition of the degree, with its multinomial
  ! count. Private to this module.
  !===================================================================!

  type :: derivative_partition

     integer(int64) :: coefficient = 1_int64
     integer, allocatable :: path_degree(:)

  end type derivative_partition

  !===================================================================!
  ! The partitions of one degree: the terms of d^n/ds^n.
  !===================================================================!

  type :: degree_partitions

     type(derivative_partition), allocatable :: term(:)

  end type degree_partitions

  !===================================================================!
  ! The assembler, constructed for a highest degree. An assembler
  ! that was not constructed stops the program at assemble.
  !===================================================================!

  !===================================================================!
  ! THE CHAIN RULE IS AN OPERATION. Composed over a statement, a
  ! degree and the paths its arguments follow, it is a map from the
  ! same inputs the statement reads to the total derivative of that
  ! order - which is what an operation is. `assemble` remains, for a
  ! caller that supplies the statement and paths at the call rather
  ! than storing them; `apply` is that call with them stored.
  !===================================================================!

  type, extends(operation) :: total_derivative

     private

     type(degree_partitions), allocatable :: of_degree(:)

     ! what the composition is over, when it is stored rather than passed
     class(operation)   , allocatable :: statement
     type(argument_path), allocatable :: along(:)
     integer                          :: order = -1

   contains

     procedure :: name  => total_derivative_name
     procedure :: apply => total_derivative_apply
     procedure :: assemble

  end type total_derivative

  interface total_derivative
     module procedure create_total_derivative
  end interface total_derivative

contains

  !===================================================================!
  ! Create the path for one input slot with derivative slots
  ! 1..degree. Filling those slots is an algorithm decision: Newton
  ! fills the correction sequence, a march fills state histories, and
  ! an empty slot contributes zero to the chain-rule assembly.
  !===================================================================!

  !===================================================================!
  ! The name of a composed rule: the statement's, and the order of
  ! the derivative taken of it.
  !===================================================================!

  pure function total_derivative_name(this) result(name)

    class(total_derivative), intent(in) :: this
    character(len=:), allocatable :: name
    character(len=12) :: digits

    write(digits,'(i0)') this % order
    if (allocated(this % statement)) then
       name = 'derivative ' // trim(digits) // ' of ' // this % statement % name()
    else
       name = 'chain rule'
    end if

  end function total_derivative_name

  !===================================================================!
  ! THE DERIVATIVE OF A STATEMENT, of the given order, along the
  ! given directions. The result is an operation: applying it to the
  ! statement's own inputs returns that derivative.
  !===================================================================!

  function derivative_of(statement, order, along) result(this)

    class(operation)   , intent(in) :: statement
    integer            , intent(in) :: order
    type(argument_path), intent(in) :: along(:)
    type(total_derivative) :: this
    character(len=250) :: message

    if (order < 0) then
       write(message,'(a,i0)') 'total_derivative: the order of a derivative must not be &
            &negative; order = ', order
       error stop trim(message)
    end if

    this = total_derivative(order)
    allocate(this % statement, source=statement)
    this % order = order
    this % along = along
    call this % declare_arguments(statement % num_arguments(), statement % contracts())

  end function derivative_of

  !===================================================================!
  ! The composed rule applied: the total derivative of the stored
  ! order, from the stored statement's partial actions along the
  ! stored paths.
  !===================================================================!

  subroutine total_derivative_apply(this, input_graph, inputs, output)

    class(total_derivative)        , intent(in)    :: this
    class(directed_graph)    , intent(in)    :: input_graph
    type(binding)             , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    if (.not. allocated(this % statement)) then
       error stop 'total_derivative: total_derivative_apply was called on a value not &
            &constructed by derivative_of() - this % statement is not allocated'
    end if
    if (.not. present(inputs)) then
       error stop 'total_derivative: total_derivative_apply was called without the &
            &statement''s inputs'
    end if

    call this % assemble(this % statement, input_graph, &
         & this % statement % bind(inputs), this % order, this % along, output)

  end subroutine total_derivative_apply

  function create_argument_path(wrt, degree) result(path)

    type(argument), intent(in) :: wrt
    integer       , intent(in) :: degree
    type(argument_path)        :: path
    character(len=250) :: message

    if (degree < 0) then
       write(message,'(a,i0)') 'argument_path: the derivative degree must be nonnegative; &
            &degree = ', degree
       error stop trim(message)
    end if

    path % wrt = wrt
    allocate(path % derivative(degree))

  end function create_argument_path

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
  ! Enumerate the partitions of every degree up to the highest. A
  ! negative highest degree stops the program.
  !===================================================================!

  function create_total_derivative(max_degree) result(this)

    integer, intent(in) :: max_degree
    type(total_derivative) :: this

    integer :: degree
    character(len=250) :: message

    if (max_degree < 0) then
       write(message,'(a,i0)') 'total_derivative: the highest degree must be nonnegative; &
            &max_degree = ', max_degree
       error stop trim(message)
    end if

    allocate(this % of_degree(max_degree))
    do degree = 1, max_degree
       call enumerate_partitions(degree, this % of_degree(degree) % term)
    end do

  end function create_total_derivative

  !===================================================================!
  ! Assemble the total derivative of the given degree. The degree
  ! and the paths are checked before any statement call runs.
  !===================================================================!

  subroutine assemble(this, statement, input_graph, inputs, degree, &
       & paths, output)

    class(total_derivative)              , intent(in)    :: this
    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(binding)                         , intent(in)    :: inputs(:)
    integer                        , intent(in)    :: degree
    type(argument_path)            , intent(in)    :: paths(:)
    class(field), allocatable, intent(inout) :: output

    real(dp), allocatable :: derivative_sum(:)
    logical :: sum_initialized
    integer :: p, num_components
    character(len=250) :: message

    if (degree < 0) then
       write(message,'(a,i0)') 'total_derivative: assemble does not support a negative degree; &
            &degree = ', degree
       error stop trim(message)
    end if
    if (.not. allocated(this % of_degree)) then
       error stop 'total_derivative: assemble was called before the assembler was constructed &
            &for a highest degree - this % of_degree is not allocated'
    end if

    call require_valid_paths(paths, statement)

    ! degree 0 is the statement's own value
    if (degree == 0) then
       call statement % apply(input_graph, inputs, output)
       return
    end if

    if (degree > size(this % of_degree)) then
       write(message,'(a,i0,a,i0)') 'total_derivative: degree must be no larger than the &
            &constructed highest degree; degree = ', degree, ', size(of_degree) = ', &
            & size(this % of_degree)
       error stop trim(message)
    end if

    sum_initialized = .false.
    num_components   = 1

    associate (partitions => this % of_degree(degree) % term)
      do p = 1, size(partitions)
         call assemble_partition(statement, input_graph, inputs, &
              & partitions(p), paths, derivative_sum, sum_initialized, num_components)
      end do
    end associate

    call write_output(statement, input_graph, inputs, derivative_sum, &
         & sum_initialized, num_components, output)

  end subroutine assemble

  !===================================================================!
  ! Check the paths: each must name an argument of the statement -
  ! one of its own declared positions, not another operation's - and
  ! no two may name the same argument, because a duplicated argument
  ! would double count chain-rule terms. Both violations stop the
  ! program.
  !===================================================================!

  pure subroutine require_valid_paths(paths, statement)

    type(argument_path), intent(in) :: paths(:)
    class(operation)   , intent(in) :: statement

    integer :: i, j

    do i = 1, size(paths)
       if (.not. statement % owns(paths(i) % wrt)) then
          error stop 'total_derivative: paths(i) % wrt is not owned by the statement - a path &
               &must name one of its declared arguments'
       end if
    end do

    do i = 1, size(paths)
       do j = i + 1, size(paths)
          if (paths(i) % wrt % matches(paths(j) % wrt)) then
             error stop 'total_derivative: paths(i) % wrt matches paths(j) % wrt for i /= j - &
                  &a duplicate argument path is rejected'
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

    integer :: entries_remaining, entry_degree

    entries_remaining = size(tuple) - position + 1

    if (entries_remaining == 1) then
       if (remaining >= minimum) then
          tuple(position) = remaining
          call append_partition(partitions, tuple)
       end if
       return
    end if

    do entry_degree = minimum, remaining / entries_remaining
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
    grown(n + 1) % path_degree = tuple
    grown(n + 1) % coefficient = partition_coefficient(sum(tuple), tuple)
    call move_alloc(grown, partitions)

  end subroutine append_partition

  !===================================================================!
  ! The multinomial count of one partition, as one exact integer
  ! division; the tuple's nondecreasing order makes every
  ! equal degrees adjacent for multiplicity counting.
  !===================================================================!

  pure function partition_coefficient(total, path_degree) result(coefficient)

    integer, intent(in) :: total
    integer, intent(in) :: path_degree(:)
    integer(int64)      :: coefficient

    integer(int64) :: denominator
    integer :: j, multiplicity

    denominator = 1_int64
    do j = 1, size(path_degree)
       denominator = denominator * factorial_int64(path_degree(j))
    end do

    multiplicity = 1
    do j = 2, size(path_degree)
       if (path_degree(j) == path_degree(j - 1)) then
          multiplicity = multiplicity + 1
       else
          denominator = denominator * factorial_int64(multiplicity)
          multiplicity = 1
       end if
    end do
    denominator = denominator * factorial_int64(multiplicity)

    coefficient = factorial_int64(total) / denominator

  end function partition_coefficient

  pure function factorial_int64(n) result(factorial)

    integer, intent(in) :: n
    integer(int64)      :: factorial

    integer :: k

    ! 21! overflows int64; stop rather than wrap
    if (n > 20) then
       error stop 'total_derivative: n exceeds 20 - 21! overflows int64, so the partition &
            &coefficient is not representable'
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

  subroutine assemble_partition(statement, input_graph, inputs, &
       & partition, paths, derivative_sum, sum_initialized, num_components)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(binding)                         , intent(in)    :: inputs(:)
    type(derivative_partition)     , intent(in)    :: partition
    type(argument_path)            , intent(in)    :: paths(:)
    real(dp), allocatable          , intent(inout) :: derivative_sum(:)
    logical                        , intent(inout) :: sum_initialized
    integer                        , intent(inout) :: num_components

    integer :: path_indices(size(partition % path_degree))
    integer :: k, j, npaths
    logical :: admitted

    k      = size(partition % path_degree)
    npaths = size(paths)
    if (npaths == 0) return

    path_indices = 1

    do

       admitted = .true.
       do j = 1, k
          if (.not. paths(path_indices(j)) % &
               & has_degree(partition % path_degree(j))) then
             admitted = .false.
             exit
          end if
       end do

       if (admitted) then
          call emit_term(statement, input_graph, inputs, partition, &
               & paths, path_indices, derivative_sum, sum_initialized, num_components)
       end if

       ! mixed-radix increment: advance the last entry, with overflow into
       ! the left
       j = k
       do
          path_indices(j) = path_indices(j) + 1
          if (path_indices(j) <= npaths) exit
          path_indices(j) = 1
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

  subroutine emit_term(statement, input_graph, inputs, partition, &
       & paths, path_indices, derivative_sum, sum_initialized, num_components)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(binding)                         , intent(in)    :: inputs(:)
    type(derivative_partition)     , intent(in)    :: partition
    type(argument_path)            , intent(in)    :: paths(:)
    integer                        , intent(in)    :: path_indices(:)
    real(dp), allocatable          , intent(inout) :: derivative_sum(:)
    logical                        , intent(inout) :: sum_initialized
    integer                        , intent(inout) :: num_components

    class(field), allocatable :: output
    type(variation) :: variations(size(path_indices))
    real(dp), allocatable :: term(:)
    integer :: j, k
    character(len=250) :: message

    k = size(path_indices)

    if (k > statement % max_degree()) then
       write(message,'(a,i0,a,i0)') 'total_derivative: the statement does not support the &
            &requested order; k = ', k, ', max_degree() = ', statement % max_degree()
       error stop trim(message)
    end if

    ! one factor per chosen path: its argument, and the derivative
    ! of the order this partition entry specifies as the direction
    do j = 1, k
       variations(j) = variation(paths(path_indices(j)) % wrt, &
            & paths(path_indices(j)) % derivative(partition % path_degree(j)) % direction)
    end do

    call statement % partial_action(input_graph, inputs, variations, output)

    call output % real_vector(term)

    if (sum_initialized) then
       if (size(term) /= size(derivative_sum)) then
          write(message,'(a,i0,a,i0)') 'total_derivative: accumulated terms must share one &
               &shape; size(term) = ', size(term), ', size(derivative_sum) = ', size(derivative_sum)
          error stop trim(message)
       end if
       derivative_sum = derivative_sum + real(partition % coefficient, dp) * term
    else
       derivative_sum = real(partition % coefficient, dp) * term
       num_components   = output % num_components()
       sum_initialized = .true.
    end if

  end subroutine emit_term

  !===================================================================!
  ! Write the total onto the statement's domain. When no term was
  ! admitted (every required derivative was unoccupied) the total
  ! is zero, with its shape taken from the statement's own value.
  !===================================================================!

  subroutine write_output(statement, input_graph, inputs, derivative_sum, &
       & sum_initialized, num_components, output)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(binding)                         , intent(in)    :: inputs(:)
    real(dp), allocatable          , intent(inout) :: derivative_sum(:)
    logical                        , intent(in)    :: sum_initialized
    integer                        , intent(in)    :: num_components
    class(field), allocatable, intent(inout) :: output

    class(field), allocatable :: value
    type(graph) :: domain
    integer         :: n_domain, width

    call statement % domain(input_graph, domain, n_domain)

    width = num_components
    if (.not. sum_initialized) then
       call statement % apply(input_graph, inputs, value)
       call value % real_vector(derivative_sum)
       derivative_sum = 0.0_dp
       width   = value % num_components()
    end if

    call emit_real('total derivative', domain, size(derivative_sum) / width, derivative_sum, output, width)

  end subroutine write_output

end module operation_chain_rule
