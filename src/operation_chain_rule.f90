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
! PARTIAL CONNECTIVITY. The partition list is also a small index
! graph. A term node is one partition, for example [1, 2] in degree
! three. Its ports read path-derivative degrees one and two; the term's
! coefficient is the multiplicity with which that slot pattern appears
! in the total derivative. The graph carries only this connectivity:
! no statement, no input tuple, no direction fields. Those values arrive
! later through argument_path, and assemble evaluates the operation's
! partial_action on the indexed slots.
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
  use field_stored   , only : stored_field

  implicit none

  private
  public :: chain_rule
  public :: argument_path
  public :: path_derivative
  public :: partial_connectivity_element
  public :: partial_connectivity
  public :: partial_connectivity_graph

  !===================================================================!
  ! One derivative of a path: occupied, it carries x^(k) as a
  ! direction field; unoccupied, it is read as zero.
  !===================================================================!

  type :: path_derivative

     logical     :: occupied = .false.
     type(stored_field) :: direction

  end type path_derivative

  !===================================================================!
  ! One argument's path: the argument of the statement it perturbs,
  ! obtained from the statement itself, and its derivative sequence
  ! - derivative(k) holds x^(k).
  !===================================================================!

  type :: argument_path

     type(argument) :: wrt
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
  ! One connectivity element in the derivative-degree graph.
  !
  ! In a mesh, the element connectivity [b, c, d] says that one
  ! triangular element reads nodes b, c and d. Here the connectivity
  ! [1, 2] says that one chain-rule term reads a first path derivative
  ! and a second path derivative. The coefficient is the multiplicity
  ! of that slot pattern.
  !===================================================================!

  type :: partial_connectivity_element

     private

     integer(int64) :: coefficient_value = 1_int64
     integer, allocatable :: slot_degree(:)

   contains

     procedure :: num_slots     => element_num_slots
     procedure :: source_degree => element_source_degree
     procedure :: coefficient   => element_coefficient

  end type partial_connectivity_element

  !===================================================================!
  ! The index connectivity of one total derivative degree.
  !
  ! For degree n, every partition is a term node:
  !
  !      [n]        reads one path derivative of degree n
  !      [1, n-1]   reads two path derivatives, degrees 1 and n-1
  !      [1, 1, 1]  reads three first-derivative slots
  !
  ! The partition entries are the labelled ports of that node. When a
  ! port says degree k, assemble ranges it over every argument_path that
  ! has x^(k) occupied. The coefficient is the exact multiplicity of
  ! the symmetric slot pattern. This is the reusable part: Newton,
  ! marching, or another driver may supply different paths and fields
  ! while sharing the same degree graph.
  !===================================================================!

  type :: partial_connectivity_graph

     private

     integer :: degree = -1
     type(partial_connectivity_element), allocatable :: elements(:)

   contains

     procedure :: order         => connectivity_order
     procedure :: num_terms     => connectivity_num_terms
     procedure :: num_elements  => connectivity_num_terms
     procedure :: num_edges     => connectivity_num_edges
     procedure :: element       => connectivity_element
     procedure :: term_size     => connectivity_term_size
     procedure :: source_degree => connectivity_source_degree
     procedure :: coefficient   => connectivity_coefficient

  end type partial_connectivity_graph

  interface partial_connectivity
     module procedure create_partial_connectivity
  end interface partial_connectivity

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
  ! and the paths are checked before any statement call runs. If a
  ! partial_connectivity_graph is supplied, its term/slot graph selects
  ! the partial_action calls; otherwise the same graph is generated here
  ! for the call.
  !===================================================================!

  subroutine assemble(this, statement, input_graph, input_data, degree, &
       & paths, output, connectivity)

    class(chain_rule)              , intent(in)    :: this
    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)                    , intent(in)    :: input_data(:)
    integer                        , intent(in)    :: degree
    type(argument_path)            , intent(in)    :: paths(:)
    class(field), allocatable, intent(inout) :: output
    type(partial_connectivity_graph), intent(in), optional :: connectivity

    type(partial_connectivity_element), allocatable :: partitions(:)
    real(dp), allocatable :: running(:)
    logical :: started
    integer :: num_components

    associate(unread => this)
    end associate

    if (degree < 0) then
       error stop 'chain_rule: degree is supported'
    end if

    call require_valid_paths(paths, statement)

    ! degree 0 is the statement's own value
    if (degree == 0) then
       call statement % apply(input_graph, input_data, output)
       return
    end if

    started = .false.
    num_components   = 1

    if (present(connectivity)) then
       call require_connectivity(connectivity, degree)
       call assemble_partitions(statement, input_graph, input_data, &
            & connectivity % elements, paths, running, started, &
            & num_components)
    else
       call enumerate_partitions(degree, partitions)
       call assemble_partitions(statement, input_graph, input_data, &
            & partitions, paths, running, started, num_components)
    end if

    call write_output(statement, input_graph, input_data, running, &
         & started, num_components, output)

  end subroutine assemble

  !===================================================================!
  ! Construct the index graph for one derivative degree. This freezes
  ! the combinatorics of D^degree S(x(s)): which term nodes exist,
  ! which derivative degree each slot reads, and the multiplicity of
  ! each term. It does not freeze S, x, or any direction value.
  !===================================================================!

  function create_partial_connectivity(degree) result(this)

    integer, intent(in) :: degree
    type(partial_connectivity_graph) :: this

    if (degree < 0) then
       error stop 'chain_rule: degree is supported'
    end if

    this % degree = degree

    if (degree == 0) then
       allocate(this % elements(0))
    else
       call enumerate_partitions(degree, this % elements)
    end if

  end function create_partial_connectivity

  pure integer function element_num_slots(this) result(num_slots)

    class(partial_connectivity_element), intent(in) :: this

    if (allocated(this % slot_degree)) then
       num_slots = size(this % slot_degree)
    else
       num_slots = 0
    end if

  end function element_num_slots

  pure integer function element_source_degree(this, slot) result(degree)

    class(partial_connectivity_element), intent(in) :: this
    integer                            , intent(in) :: slot

    call require_element_slot(this, slot)

    degree = this % slot_degree(slot)

  end function element_source_degree

  pure integer(int64) function element_coefficient(this) result(coefficient)

    class(partial_connectivity_element), intent(in) :: this

    coefficient = this % coefficient_value

  end function element_coefficient

  pure subroutine require_element_slot(element, slot)

    class(partial_connectivity_element), intent(in) :: element
    integer                            , intent(in) :: slot

    if (.not. allocated(element % slot_degree)) then
       error stop 'chain_rule: the connectivity element is constructed'
    end if
    if (slot < 1 .or. slot > size(element % slot_degree)) then
       error stop 'chain_rule: the connectivity element slot exists'
    end if

  end subroutine require_element_slot

  pure integer function connectivity_order(this) result(order)

    class(partial_connectivity_graph), intent(in) :: this

    order = this % degree

  end function connectivity_order

  pure integer function connectivity_num_terms(this) result(num_terms)

    class(partial_connectivity_graph), intent(in) :: this

    if (allocated(this % elements)) then
       num_terms = size(this % elements)
    else
       num_terms = 0
    end if

  end function connectivity_num_terms

  pure integer function connectivity_num_edges(this) result(num_edges)

    class(partial_connectivity_graph), intent(in) :: this

    integer :: p

    num_edges = 0
    if (.not. allocated(this % elements)) return

    do p = 1, size(this % elements)
       num_edges = num_edges + size(this % elements(p) % slot_degree)
    end do

  end function connectivity_num_edges

  pure function connectivity_element(this, term) result(element)

    class(partial_connectivity_graph), intent(in) :: this
    integer                          , intent(in) :: term
    type(partial_connectivity_element) :: element

    call require_term(this, term)

    element = this % elements(term)

  end function connectivity_element

  pure integer function connectivity_term_size(this, term) result(term_size)

    class(partial_connectivity_graph), intent(in) :: this
    integer                          , intent(in) :: term

    call require_term(this, term)

    term_size = size(this % elements(term) % slot_degree)

  end function connectivity_term_size

  pure integer function connectivity_source_degree(this, term, slot) result(degree)

    class(partial_connectivity_graph), intent(in) :: this
    integer                          , intent(in) :: term, slot

    call require_slot(this, term, slot)

    degree = this % elements(term) % slot_degree(slot)

  end function connectivity_source_degree

  pure integer(int64) function connectivity_coefficient(this, term) result(coefficient)

    class(partial_connectivity_graph), intent(in) :: this
    integer                          , intent(in) :: term

    call require_term(this, term)

    coefficient = this % elements(term) % coefficient_value

  end function connectivity_coefficient

  pure subroutine require_connectivity(connectivity, degree)

    type(partial_connectivity_graph), intent(in) :: connectivity
    integer                         , intent(in) :: degree

    if (connectivity % degree /= degree) then
       error stop 'chain_rule: the connectivity degree matches the requested degree'
    end if
    if (.not. allocated(connectivity % elements)) then
       error stop 'chain_rule: the connectivity graph is constructed'
    end if

  end subroutine require_connectivity

  pure subroutine require_term(connectivity, term)

    class(partial_connectivity_graph), intent(in) :: connectivity
    integer                          , intent(in) :: term

    if (.not. allocated(connectivity % elements)) then
       error stop 'chain_rule: the connectivity graph is constructed'
    end if
    if (term < 1 .or. term > size(connectivity % elements)) then
       error stop 'chain_rule: the connectivity term exists'
    end if

  end subroutine require_term

  pure subroutine require_slot(connectivity, term, slot)

    class(partial_connectivity_graph), intent(in) :: connectivity
    integer                          , intent(in) :: term, slot

    call require_term(connectivity, term)

    if (slot < 1 .or. slot > size(connectivity % elements(term) % slot_degree)) then
       error stop 'chain_rule: the connectivity slot exists'
    end if

  end subroutine require_slot

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
          error stop 'chain_rule: a path names an argument of the statement'
       end if
    end do

    do i = 1, size(paths)
       do j = i + 1, size(paths)
          if (paths(i) % wrt % matches(paths(j) % wrt)) then
             error stop 'chain_rule: duplicate argument path is refused'
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
    type(partial_connectivity_element), allocatable, intent(out) :: partitions(:)

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
    type(partial_connectivity_element), allocatable, intent(inout) :: partitions(:)

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

    type(partial_connectivity_element), allocatable, intent(inout) :: partitions(:)
    integer                                , intent(in)    :: tuple(:)

    type(partial_connectivity_element), allocatable :: grown(:)
    integer :: n

    n = size(partitions)
    allocate(grown(n + 1))
    grown(1:n) = partitions
    grown(n + 1) % slot_degree = tuple
    grown(n + 1) % coefficient_value = partition_coefficient(sum(tuple), tuple)
    call move_alloc(grown, partitions)

  end subroutine append_partition

  !===================================================================!
  ! The multinomial count of one partition, as one exact integer
  ! division; the tuple's nondecreasing order makes every
  ! multiplicity a contiguous run.
  !===================================================================!

  pure function partition_coefficient(total, path_degree) result(coefficient)

    integer, intent(in) :: total
    integer, intent(in) :: path_degree(:)
    integer(int64)      :: coefficient

    integer(int64) :: denominator
    integer :: j, run

    denominator = 1_int64
    do j = 1, size(path_degree)
       denominator = denominator * factorial_int64(path_degree(j))
    end do

    run = 1
    do j = 2, size(path_degree)
       if (path_degree(j) == path_degree(j - 1)) then
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
    type(partial_connectivity_element)     , intent(in)    :: partition
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

  subroutine assemble_partitions(statement, input_graph, input_data, &
       & partitions, paths, running, started, num_components)

    class(operation)               , intent(in)    :: statement
    class(directed_graph)          , intent(in)    :: input_graph
    type(stored_field)             , intent(in)    :: input_data(:)
    type(partial_connectivity_element)     , intent(in)    :: partitions(:)
    type(argument_path)            , intent(in)    :: paths(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: num_components

    integer :: p

    do p = 1, size(partitions)
       call assemble_partition(statement, input_graph, input_data, &
            & partitions(p), paths, running, started, num_components)
    end do

  end subroutine assemble_partitions

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
    type(partial_connectivity_element)     , intent(in)    :: partition
    type(argument_path)            , intent(in)    :: paths(:)
    integer                        , intent(in)    :: chosen(:)
    real(dp), allocatable          , intent(inout) :: running(:)
    logical                        , intent(inout) :: started
    integer                        , intent(inout) :: num_components

    class(field), allocatable :: output
    type(variation) :: variations(size(chosen))
    real(dp), allocatable :: term(:)
    integer :: j, k

    k = size(chosen)

    if (k > statement % max_degree()) then
       error stop 'chain_rule: the statement supports the requested order'
    end if

    ! one factor per chosen path: its argument, and the derivative
    ! of the order this partition entry asks for as the direction
    do j = 1, k
       variations(j) = variation(paths(chosen(j)) % wrt, &
            & paths(chosen(j)) % derivative(partition % slot_degree(j)) % direction)
    end do

    call statement % partial_action(input_graph, input_data, variations, output)

    call output % real_vector(term)

    if (started) then
       if (size(term) /= size(running)) then
          error stop 'chain_rule: accumulated terms share one shape'
       end if
       running = running + real(partition % coefficient_value, dp) * term
    else
       running = real(partition % coefficient_value, dp) * term
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
