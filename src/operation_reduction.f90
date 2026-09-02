!=====================================================================!
! Concrete graph reductions.
!
! A reduction maps a field on a support to one number. There is a
! single concrete type here storing a rule, rather than one class per
! rule, for the same reason there is one concrete field per side: a
! caller can then store reductions in a plain array, and adding a new
! rule requires a case rather than a class.
!
!=====================================================================!
!
!                          WHY FOUR STEPS
!
! Summation appears to need one procedure. It needs four, because the
! graph may be partitioned:
!
!   part 1  [2 2 2]  --accumulate-->  (sum 6, count 3) --+
!                                                        +--> combine
!   part 2  [5 9]    --accumulate-->  (sum 14, count 2) -+        |
!                                                                 v
!                                                    (sum 20, count 5)
!                                                                 |
!                                                             finalize
!                                                                 |
!                                                                 v
!                                                            J = 4.0
!
! Part one averages to 2 and part two to 7. The mean of those is 4.5,
! which is incorrect; the result is 4, because 20/5 = 4. The sum and
! the count must be combined together and the division must occur
! once, at the end. Finalizing each part separately makes a parallel
! run differ from a serial one without any error report.
!
! Minimum and maximum may be finalized per part; average and norm may
! not. The four steps add no operations for the first group and are
! required by the second, so every rule uses them.
!
!=====================================================================!
!
!                        THE MEASURE
!
! Pass a measure and a bare sum becomes an integral. Weight each cell
! by its volume, or each face by its area, and the result no longer
! depends on the mesh resolution:
!
!      sum        J = sum q_i
!      integral   J = sum q_i V_i          <- measure is the volume
!      average    J = sum q_i V_i / sum V_i
!      norm       J = ( sum |q_i|^p V_i )^(1/p)
!
! The measure stores one value per entry. A field several components
! wide weights every component of an entry by that entry's measure.
!
! The measure position is also the inner product's second field: a sum
! reduced with measure v returns the sum of q times v, the product
! <q, v>.
!
!=====================================================================!
!
!                        WHICH KINDS
!
! Summing works on real and on complex values, and the complex case is
! required: a complex-step objective is a weighted sum, and its
! derivative is the imaginary part. Without it the reason the
! functional stores complex values is removed.
!
! Ordering rules - minimum, maximum, norm - are real only, because
! the complex numbers are not ordered.
!
! All and any work on logical fields. They let a predicate such as
! "this graph is acyclic" evaluate to true or false rather than to a
! one or a zero.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_reduction

  use util_precision  , only : dp
  use view_directed   , only : directed_graph
  use field_calculus  , only : field
  use graph_fractal      , only : graph
  use field_calculus  , only : FIELD_REAL, FIELD_COMPLEX
  use field_calculus  , only : FIELD_LOGICAL
  use operation_action  , only : operation, contract
  use operation_action, only : binding, bound_value
  use operation_action, only : emit
  use field_calculus  , only : functional
  use view_directed   , only : SIDE_VERTEX
  use field_stored    , only : stored_field
  use field_functional, only : stored_functional

  implicit none

  private
  public :: reduction
  public :: REDUCE_SUM, REDUCE_AVERAGE, REDUCE_MINIMUM, REDUCE_MAXIMUM
  public :: REDUCE_NORM, REDUCE_COUNT, REDUCE_ALL, REDUCE_ANY
  public :: broadcast
  public :: BROADCAST_COPY, BROADCAST_SHARE

  integer, parameter :: REDUCE_SUM     = 1
  integer, parameter :: REDUCE_AVERAGE = 2
  integer, parameter :: REDUCE_MINIMUM = 3
  integer, parameter :: REDUCE_MAXIMUM = 4
  integer, parameter :: REDUCE_NORM    = 5
  integer, parameter :: REDUCE_COUNT   = 6
  integer, parameter :: REDUCE_ALL     = 7
  integer, parameter :: REDUCE_ANY     = 8

  integer, parameter :: BROADCAST_COPY  = 1   ! transpose of a sum
  integer, parameter :: BROADCAST_SHARE = 2   ! transpose of an average

  !===================================================================!
  ! REDUCTION. Many values become one: field -> functional, storing
  ! the rule it follows. Four staged steps so a partitioned run
  ! combines partial states before finishing once - a mean of
  ! part-means is not the mean - plus the one-call form. The measure
  ! argument turns a bare sum into an integral and is the inner
  ! product's second field.
  !===================================================================!

  type, extends(operation) :: reduction

     integer  :: rule  = REDUCE_SUM
     real(dp) :: power = 2.0_dp     ! which norm, when the rule is a norm

     !----------------------------------------------------------------!
     ! The one-entry domain the result is stored on, declared once at
     ! construction. Declaring it per call to domain() would make two
     ! calls two domains, and a functional built on the first would
     ! fail same_as against the second.
     !----------------------------------------------------------------!

     type(graph), private :: result_domain

   contains

     !----------------------------------------------------------------!
     ! Initialize, accumulate values, combine two parts, finalize once.
     !----------------------------------------------------------------!

     procedure :: initialize
     procedure :: accumulate
     procedure :: combine
     procedure :: finalize

     !----------------------------------------------------------------!
     ! All four steps in one call, for a caller with the whole field.
     !----------------------------------------------------------------!

     procedure :: reduce

     !----------------------------------------------------------------!
     ! The operation interface: the field reduced over its own domain,
     ! the measure passed as the second input field.
     !----------------------------------------------------------------!

     procedure :: name   => reduction_name
     procedure :: domain => reduction_domain
     procedure :: apply  => reduction_apply

  end type reduction

  !===================================================================!
  ! Constructor. Specify the rule; the power is read only for a norm.
  !===================================================================!

  interface reduction
     module procedure create
  end interface reduction

  interface broadcast
     module procedure create_broadcast
  end interface broadcast

  !===================================================================!
  ! The reduction's transpose: one value fills a field. Copy is the
  ! transpose of a sum, share the transpose of an average, and the
  ! identity reduce(broadcast(J)) = J fixes them. The constructor
  ! broadcast(BROADCAST_SHARE) builds one and declares its one
  ! argument.
  !===================================================================!

  type, extends(operation) :: broadcast

     integer :: rule = BROADCAST_COPY

   contains

     procedure :: broadcast => broadcast_functional

     !----------------------------------------------------------------!
     ! The operation interface: the functional is the one input
     ! field, and the filled field is returned on the graph's own
     ! vertices.
     !----------------------------------------------------------------!

     procedure :: name   => broadcast_name
     procedure :: apply  => broadcast_apply

  end type broadcast

contains

  !===================================================================!
  ! Build a reduction that follows one rule.
  !===================================================================!

  type(reduction) function create(rule, power) result(this)

    integer , intent(in)           :: rule
    real(dp), intent(in), optional :: power

    this % rule = rule

    if (present(power)) this % power = power

    call this % result_domain % declare()

    ! two readable positions, the values and the measure; a call may
    ! pass one or none
    select case (rule)
    case (REDUCE_ALL, REDUCE_ANY)
       call this % declare_arguments(2, [contract([FIELD_LOGICAL], 1), &
            & contract([FIELD_LOGICAL], 1)])
    case (REDUCE_SUM, REDUCE_AVERAGE)
       call this % declare_arguments(2, [contract([FIELD_REAL, FIELD_COMPLEX], 1), &
            & contract([FIELD_REAL, FIELD_COMPLEX], 1)])
    case default
       call this % declare_arguments(2, [contract(FIELD_REAL, 1), contract(FIELD_REAL, 1)])
    end select

  end function create

  !===================================================================!
  ! Build a broadcast that follows one rule: one argument, the
  ! functional it distributes.
  !===================================================================!

  type(broadcast) function create_broadcast(rule) result(this)

    integer, intent(in) :: rule

    this % rule = rule

    call this % declare_arguments(1, [contract([FIELD_REAL, FIELD_COMPLEX], 1)])

  end function create_broadcast

  !===================================================================!
  ! The initial state of a reduction. Zero for a sum, huge() for a
  ! minimum, true for an "all".
  !
  ! Summing starts real. If the first field is complex, accumulate
  ! promotes the state - a complex zero and a real zero are the same
  ! number, so either start is exact.
  !===================================================================!

  pure subroutine initialize(this, state)

    class(reduction), intent(in)                            :: this
    class(functional), allocatable, intent(inout)     :: state

    if (allocated(state)) deallocate(state)
    allocate(stored_functional :: state)

    select type (state)
    type is (stored_functional)

       state % tally  = 0.0_dp
       state % weight = 0.0_dp

       select case (this % rule)
       case (REDUCE_MINIMUM)
          call state % set_real_value(huge(1.0_dp))
       case (REDUCE_MAXIMUM)
          call state % set_real_value(-huge(1.0_dp))
       case (REDUCE_ALL)
          call state % set_logical_value(.true.)
       case (REDUCE_ANY)
          call state % set_logical_value(.false.)
       case default
          call state % set_real_value(0.0_dp)
       end select

    end select

  end subroutine initialize

  !===================================================================!
  ! Add one part's values into the running state.
  !===================================================================!

  pure subroutine accumulate(this, values, state, measure)

    class(reduction)       , intent(in)    :: this
    class(field)     , intent(in)    :: values
    class(functional), intent(inout) :: state
    class(field)     , intent(in), optional :: measure

    real(dp)   , allocatable :: v(:), m(:)
    complex(dp), allocatable :: cv(:)
    logical    , allocatable :: lv(:)

    real(dp)    :: acc, wi
    complex(dp) :: cacc
    logical     :: lacc, weighted
    integer     :: i, c, k, num_components, nentry

    num_components  = values % num_components()
    nentry = values % num_entries()

    call weights_of(measure, nentry, m, weighted)

    select case (this % rule)

    case (REDUCE_ALL, REDUCE_ANY)

       call values % logical_vector(lv)
       call state % logical_value(lacc)
       if (this % rule == REDUCE_ALL) then
          lacc = lacc .and. all(lv)
       else
          lacc = lacc .or. any(lv)
       end if
       call state % set_logical_value(lacc)

    case (REDUCE_COUNT)

       select type (state)
       type is (stored_functional)
          state % tally = state % tally + real(nentry, dp)
       end select

    case default

       ! Complex values follow the complex branch; all others the
       ! real branch. Only summing is defined for complex values,
       ! since the complex numbers are not ordered.
       if (values % value_kind() == FIELD_COMPLEX) then

          call values % complex_vector(cv)
          call state % complex_value(cacc)
          do i = 1, nentry
             wi = 1.0_dp
             if (weighted) wi = m(i)
             do c = 1, num_components
                k = (i - 1) * num_components + c
                if (k <= size(cv)) cacc = cacc + cv(k) * wi
             end do
          end do
          call state % set_complex_value(cacc)

       else

          call values % real_vector(v)

          select case (this % rule)

          case (REDUCE_SUM)
             call state % real_value(acc)
             do i = 1, nentry
                wi = 1.0_dp
                if (weighted) wi = m(i)
                do c = 1, num_components
                   k = (i - 1) * num_components + c
                   if (k <= size(v)) acc = acc + v(k) * wi
                end do
             end do
             call state % set_real_value(acc)

          case (REDUCE_AVERAGE, REDUCE_NORM)
             ! Both store a running total and a running weight, and
             ! both divide or take the root only at finalize.
             select type (state)
             type is (stored_functional)
                do i = 1, nentry
                   wi = 1.0_dp
                   if (weighted) wi = m(i)
                   do c = 1, num_components
                      k = (i - 1) * num_components + c
                      if (k <= size(v)) then
                         if (this % rule == REDUCE_AVERAGE) then
                            state % tally = state % tally + v(k) * wi
                         else
                            state % tally = state % tally + abs(v(k))**this % power * wi
                         end if
                         state % weight = state % weight + wi
                      end if
                   end do
                end do
             end select

          case (REDUCE_MINIMUM)
             call state % real_value(acc)
             do k = 1, size(v)
                acc = min(acc, v(k))
             end do
             call state % set_real_value(acc)

          case (REDUCE_MAXIMUM)
             call state % real_value(acc)
             do k = 1, size(v)
                acc = max(acc, v(k))
             end do
             call state % set_real_value(acc)

          end select

       end if

    end select

  end subroutine accumulate

  !===================================================================!
  ! One weight per entry: the measure if there is one, otherwise one.
  !===================================================================!

  pure subroutine weights_of(measure, nentry, m, weighted)

    class(field), intent(in), optional :: measure
    integer           , intent(in)           :: nentry
    real(dp), allocatable, intent(out)       :: m(:)
    logical              , intent(out)       :: weighted

    real(dp), allocatable :: raw(:)

    ! an absent measure weights every entry by one; that is a property
    ! of the reduction, so no vector of ones is allocated and read
    weighted = .false.
    if (present(measure)) then
       call measure % real_vector(raw)
       if (size(raw) >= nentry) then
          allocate(m(max(nentry, 1)))
          m = 1.0_dp
          m(1:nentry) = raw(1:nentry)
          weighted    = .true.
       end if
    end if
    if (.not. weighted) allocate(m(0))

  end subroutine weights_of

  !===================================================================!
  ! Combine two partial results. The result must not depend on the
  ! order of the parts; otherwise a parallel run would depend on
  ! which image finishes first.
  !===================================================================!

  pure subroutine combine(this, left, right, combined)

    class(reduction)       , intent(in)    :: this
    class(functional), intent(in)    :: left
    class(functional), intent(in)    :: right
    class(functional), allocatable, intent(inout) :: combined

    real(dp)    :: a, b
    complex(dp) :: ca, cb
    logical     :: la, lb

    if (allocated(combined)) deallocate(combined)
    allocate(stored_functional :: combined)

    select case (this % rule)

    case (REDUCE_ALL, REDUCE_ANY)

       call left % logical_value(la)
       call right % logical_value(lb)
       if (this % rule == REDUCE_ALL) then
          call combined % set_logical_value(la .and. lb)
       else
          call combined % set_logical_value(la .or. lb)
       end if

    case (REDUCE_AVERAGE, REDUCE_NORM, REDUCE_COUNT)

       ! The running totals add; the division and the root are
       ! deferred to finalize.
       select type (combined)
       type is (stored_functional)
          select type (left)
          type is (stored_functional)
             select type (right)
             type is (stored_functional)
                combined % tally  = left % tally  + right % tally
                combined % weight = left % weight + right % weight
             end select
          end select
       end select

    case (REDUCE_MINIMUM, REDUCE_MAXIMUM)

       call left % real_value(a)
       call right % real_value(b)
       if (this % rule == REDUCE_MINIMUM) then
          call combined % set_real_value(min(a, b))
       else
          call combined % set_real_value(max(a, b))
       end if

    case default

       ! Summing, in the complex or the real branch according to the
       ! parts' value kinds.
       if (left % value_kind() == FIELD_COMPLEX .or. &
            & right % value_kind() == FIELD_COMPLEX) then
          call left % complex_value(ca)
          call right % complex_value(cb)
          call combined % set_complex_value(ca + cb)
       else
          call left % real_value(a)
          call right % real_value(b)
          call combined % set_real_value(a + b)
       end if

    end select

  end subroutine combine

  !===================================================================!
  ! Finalize, once, after every part has been combined. Here an
  ! average divides and a norm takes its root; doing either earlier
  ! is the error the four steps prevent.
  !===================================================================!

  pure subroutine finalize(this, state, scalar)

    class(reduction)       , intent(in) :: this
    class(functional), intent(in) :: state
    class(functional), allocatable, intent(inout) :: scalar

    if (allocated(scalar)) deallocate(scalar)
    allocate(scalar, source=state)

    select case (this % rule)

    case (REDUCE_AVERAGE)

       select type (state)
       type is (stored_functional)
          if (state % weight > 0.0_dp) then
             call scalar % set_real_value(state % tally / state % weight)
          else
             call scalar % set_real_value(0.0_dp)
          end if
       end select

    case (REDUCE_NORM)

       select type (state)
       type is (stored_functional)
          if (state % tally > 0.0_dp) then
             call scalar % set_real_value(state % tally**(1.0_dp / this % power))
          else
             call scalar % set_real_value(0.0_dp)
          end if
       end select

    case (REDUCE_COUNT)

       select type (state)
       type is (stored_functional)
          call scalar % set_integer_vector([nint(state % tally)])
       end select

    end select

  end subroutine finalize

  !===================================================================!
  ! All four steps for a caller with the whole graph.
  !
  ! Not pure. This is the one procedure where a reduction distributed
  ! across images may communicate with the other images, and a
  ! distributed reduction would sum here before finalizing.
  !===================================================================!

  subroutine reduce(this, values, scalar, measure)

    class(reduction)    , intent(in) :: this
    class(field)  , intent(in) :: values
    class(functional), allocatable, intent(inout) :: scalar
    class(field)  , intent(in), optional :: measure

    class(functional), allocatable :: state

    call this % initialize(state)
    call this % accumulate(values, state, measure)
    call this % finalize(state, scalar)

  end subroutine reduce

  !===================================================================!
  ! The reduction's operation interface. The field is reduced over
  ! its own domain; a second input field is the measure; the
  ! functional is returned as the output, since a functional IS a
  ! field.
  !===================================================================!

  pure function reduction_name(this) result(name)

    class(reduction), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'reduction'

  end function reduction_name

  subroutine reduction_domain(this, input_graph, domain, num_entries)

    class(reduction), intent(in)           :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    associate (u1 => input_graph); end associate

    ! The result's domain is this reduction's own one-entry domain.
    domain   = this % result_domain
    num_entries = 1

  end subroutine reduction_domain

  subroutine reduction_apply(this, input_graph, inputs, output)

    class(reduction), intent(in)                   :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    class(functional), allocatable :: reduced
    class(field), allocatable :: values, measure

    associate (u1 => input_graph); end associate

    if (present(inputs)) then
       if (size(inputs) < 1 .or. size(inputs) > this % num_arguments()) then
          error stop 'operation_reduction: values and an optional measure are bound'
       end if
       if (size(inputs) >= 2) then
          call bound_value(inputs, this % argument(1), values)
          call bound_value(inputs, this % argument(2), measure)
          call reduce_measured(this, values, measure, reduced)
       else
          call bound_value(inputs, this % argument(1), values)
          call this % reduce(values, reduced)
       end if
    else
       call this % initialize(reduced)
    end if

    if (allocated(output)) deallocate(output)
    allocate(output, source=reduced)

  end subroutine reduction_apply

  ! A separate procedure so the measure is passed to a required dummy
  ! argument: gfortran crashes (gfc_get_descriptor_field; verified on
  ! 11.4, 13.4, 15.2, and 16.0 trunk) when an optional class dummy is
  ! passed a polymorphic array element. Delete when the compiler
  ! compiles the direct call without crashing.
  subroutine reduce_measured(this, u, v, image)

    class(reduction), intent(in)   :: this
    class(field), intent(in) :: u
    class(field), intent(in) :: v
    class(functional), allocatable, intent(inout) :: image

    call this % reduce(u, image, measure=v)

  end subroutine reduce_measured

  !===================================================================!
  ! The broadcast's operation interface: the transpose of the
  ! reduction's. The one input field must be a functional; the filled
  ! field is returned on the graph's vertices.
  !===================================================================!

  pure function broadcast_name(this) result(name)

    class(broadcast), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'broadcast'

  end function broadcast_name
  subroutine broadcast_apply(this, input_graph, inputs, output)

    class(broadcast), intent(in)                   :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field) :: out
    class(field), allocatable :: value

    out = stored_field('broadcast', input_graph % vertex_set(), input_graph % num_vertices())

    if (present(inputs)) then
       call bound_value(inputs, this % argument(1), value)
       select type (f => value)
       class is (functional)
          call this % broadcast(f, out)
       class default
          error stop 'broadcast: the operation interface requires a functional'
       end select
    end if

    call emit(out, output)

  end subroutine broadcast_apply

  !===================================================================!
  ! Fill every stored value of the field from the functional's one.
  ! Copy assigns each value J; share assigns each value J divided by
  ! the count of stored values, so a subsequent sum returns J. A real
  ! J fills a real field; a complex J fills a complex field and
  ! transports a complex-step seed; any other kind fills zeros,
  ! following the value-kind rule the fields declare.
  !===================================================================!

  pure subroutine broadcast_functional(this, scalar, values)

    class(broadcast)       , intent(in)    :: this
    class(functional), intent(in)    :: scalar
    class(field)     , intent(inout) :: values

    real(dp)    :: value
    complex(dp) :: complex_value
    integer     :: n, i

    n = values % num_entries() * max(values % num_components(), 1)

    if (scalar % value_kind() == FIELD_COMPLEX) then
       call scalar % complex_value(complex_value)
       if (this % rule == BROADCAST_SHARE .and. n > 0) then
          complex_value = complex_value / real(n, dp)
       end if
       call values % set_complex_vector([(complex_value, i = 1, n)])
    else
       call scalar % real_value(value)
       if (this % rule == BROADCAST_SHARE .and. n > 0) then
          value = value / real(n, dp)
       end if
       call values % set_real_vector([(value, i = 1, n)])
    end if

  end subroutine broadcast_functional

end module operation_reduction
