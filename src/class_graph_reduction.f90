!=====================================================================!
! Concrete graph reductions.
!
! A reduction turns a field on a support into one number. There is a
! single concrete type here carrying a rule, rather than one class per
! rule, for the same reason there is one concrete field per side: a
! caller can then hold reductions in a plain array, and adding a new
! rule costs a case rather than a class.
!
!=====================================================================!
!
!                          WHY FOUR STEPS
!
! Adding numbers up looks like it needs one procedure. It needs four,
! and the reason is that the graph may be in pieces:
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
! Part one averages to 2 and part two to 7. Averaging those gives 4.5,
! which is wrong; the answer is 4, because twenty divided by five is
! four. The sum and the count have to travel together and the division
! has to happen once, at the end. Finish early on each part and a
! parallel run quietly disagrees with a serial one.
!
! Minimum and maximum are safe to finish early; average and norm are
! not. The four steps cost nothing for the safe ones and save the
! unsafe ones, so every rule uses them.
!
!=====================================================================!
!
!                        THE MEASURE
!
! Pass a measure and a bare sum becomes an integral. Weight each cell
! by its volume, or each face by its area, and the answer stops
! depending on how finely the mesh was cut:
!
!      sum        J = sum q_i
!      integral   J = sum q_i V_i          <- measure is the volume
!      average    J = sum q_i V_i / sum V_i
!      norm       J = ( sum |q_i|^p V_i )^(1/p)
!
! The measure carries one value per entry. A field several components
! wide weights every component of an entry by that entry's measure.
!
! The measure seat is also the inner product's second field: a sum
! reduced with measure v answers the sum of q times v, the product
! <q, v>.
!
!=====================================================================!
!
!                        WHICH KINDS
!
! Summing works on real and on complex values, and the complex case is
! not a curiosity: a complex-step objective is a weighted sum, and its
! derivative is the imaginary part. Lose that and the whole reason the
! functional carries complex is lost with it.
!
! Ordering rules - minimum, maximum, norm - are real only, because
! complex numbers do not order.
!
! All and any work on logical fields. They are what let a question
! such as "is this graph acyclic" come back as true or false rather
! than as a one or a zero.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_reduction

  use iso_fortran_env       , only : dp => REAL64
  use graph_grammar         , only : graph_field
  use graph_grammar         , only : GRAPH_FIELD_REAL, GRAPH_FIELD_COMPLEX
  use graph_grammar         , only : GRAPH_FIELD_LOGICAL
  use graph_calculus        , only : graph_reduction, graph_broadcast
  use graph_calculus        , only : graph_functional
  use class_graph_functional, only : scalar_result => functional

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
  ! One reduction, carrying the rule it follows.
  !===================================================================!

  type, extends(graph_reduction) :: reduction

     integer  :: rule  = REDUCE_SUM
     real(dp) :: power = 2.0_dp     ! which norm, when the rule is a norm

   contains

     !----------------------------------------------------------------!
     ! Start empty, take values in, join two parts, finish once.
     !----------------------------------------------------------------!

     procedure :: initialize
     procedure :: accumulate
     procedure :: combine
     procedure :: finalize

     !----------------------------------------------------------------!
     ! All four at once, for a caller holding the whole thing.
     !----------------------------------------------------------------!

     procedure :: reduce

  end type reduction

  !===================================================================!
  ! Constructor. Name the rule; the power matters only for a norm.
  !===================================================================!

  interface reduction
     module procedure create
  end interface reduction

  !===================================================================!
  ! The reduction's pair: one value fills a field. Copy is the
  ! transpose of a sum, share the transpose of an average, and the
  ! round trip reduce(broadcast(J)) = J pins them. The rule
  ! component is public, so the intrinsic structure constructor
  ! broadcast(BROADCAST_SHARE) builds one.
  !===================================================================!

  type, extends(graph_broadcast) :: broadcast

     integer :: rule = BROADCAST_COPY

   contains

     procedure :: broadcast => broadcast_functional

  end type broadcast

contains

  !===================================================================!
  ! Build a reduction that follows one rule.
  !===================================================================!

  pure type(reduction) function create(rule, power) result(this)

    integer , intent(in)           :: rule
    real(dp), intent(in), optional :: power

    this % rule = rule

    if (present(power)) this % power = power

  end function create

  !===================================================================!
  ! The empty answer a reduction starts from. Zero for a sum, the largest
  ! number there is for a minimum, true for an "all".
  !
  ! Summing starts real. If the first field turns out to be complex,
  ! accumulate promotes it - a complex zero and a real zero are the
  ! same number, so nothing is lost by starting either way.
  !===================================================================!

  pure subroutine initialize(this, state)

    class(reduction), intent(in)                            :: this
    class(graph_functional), allocatable, intent(inout)     :: state

    if (allocated(state)) deallocate(state)
    allocate(scalar_result :: state)

    select type (state)
    type is (scalar_result)

       state % tally  = 0.0_dp
       state % weight = 0.0_dp

       select case (this % rule)
       case (REDUCE_MINIMUM)
          call set_real(state, huge(1.0_dp))
       case (REDUCE_MAXIMUM)
          call set_real(state, -huge(1.0_dp))
       case (REDUCE_ALL)
          call set_logical(state, .true.)
       case (REDUCE_ANY)
          call set_logical(state, .false.)
       case default
          call set_real(state, 0.0_dp)
       end select

    end select

  end subroutine initialize

  !===================================================================!
  ! Add one part's values into the running answer.
  !===================================================================!

  pure subroutine accumulate(this, field, state, measure)

    class(reduction)       , intent(in)    :: this
    class(graph_field)     , intent(in)    :: field
    class(graph_functional), intent(inout) :: state
    class(graph_field)     , intent(in), optional :: measure

    real(dp)   , allocatable :: v(:), m(:)
    complex(dp), allocatable :: cv(:)
    logical    , allocatable :: lv(:)

    real(dp)    :: acc
    complex(dp) :: cacc
    logical     :: lacc
    integer     :: i, c, k, ncomp, nentry

    ncomp  = field % num_components()
    nentry = field % num_entries()

    call weights_of(measure, nentry, m)

    select case (this % rule)

    case (REDUCE_ALL, REDUCE_ANY)

       call field % get_logical_vector(lv)
       call get_logical(state, lacc)
       if (this % rule == REDUCE_ALL) then
          lacc = lacc .and. all(lv)
       else
          lacc = lacc .or. any(lv)
       end if
       call set_logical(state, lacc)

    case (REDUCE_COUNT)

       select type (state)
       type is (scalar_result)
          state % tally = state % tally + real(nentry, dp)
       end select

    case default

       ! A complex field takes the complex road; everything else the
       ! real one. Only summing is defined for complex values, since
       ! ordering them has no meaning.
       if (field % value_kind() == GRAPH_FIELD_COMPLEX) then

          call field % get_complex_vector(cv)
          call get_complex(state, cacc)
          do i = 1, nentry
             do c = 1, ncomp
                k = (i - 1) * ncomp + c
                if (k <= size(cv)) cacc = cacc + cv(k) * m(i)
             end do
          end do
          call set_complex(state, cacc)

       else

          call field % get_real_vector(v)

          select case (this % rule)

          case (REDUCE_SUM)
             call get_real(state, acc)
             do i = 1, nentry
                do c = 1, ncomp
                   k = (i - 1) * ncomp + c
                   if (k <= size(v)) acc = acc + v(k) * m(i)
                end do
             end do
             call set_real(state, acc)

          case (REDUCE_AVERAGE, REDUCE_NORM)
             ! Both carry a running total and a running weight, and
             ! both divide or root only at the very end.
             select type (state)
             type is (scalar_result)
                do i = 1, nentry
                   do c = 1, ncomp
                      k = (i - 1) * ncomp + c
                      if (k <= size(v)) then
                         if (this % rule == REDUCE_AVERAGE) then
                            state % tally = state % tally + v(k) * m(i)
                         else
                            state % tally = state % tally + abs(v(k))**this % power * m(i)
                         end if
                         state % weight = state % weight + m(i)
                      end if
                   end do
                end do
             end select

          case (REDUCE_MINIMUM)
             call get_real(state, acc)
             do k = 1, size(v)
                acc = min(acc, v(k))
             end do
             call set_real(state, acc)

          case (REDUCE_MAXIMUM)
             call get_real(state, acc)
             do k = 1, size(v)
                acc = max(acc, v(k))
             end do
             call set_real(state, acc)

          end select

       end if

    end select

  end subroutine accumulate

  !===================================================================!
  ! One weight per entry: the measure if there is one, otherwise one.
  !===================================================================!

  pure subroutine weights_of(measure, nentry, m)

    class(graph_field), intent(in), optional :: measure
    integer           , intent(in)           :: nentry
    real(dp), allocatable, intent(out)       :: m(:)

    real(dp), allocatable :: raw(:)

    allocate(m(max(nentry, 1)))
    m = 1.0_dp

    if (present(measure)) then
       call measure % get_real_vector(raw)
       if (size(raw) >= nentry) m(1:nentry) = raw(1:nentry)
    end if

  end subroutine weights_of

  !===================================================================!
  ! Join two part answers. The result must not depend on the order
  ! the parts arrive in; otherwise a parallel run would depend on
  ! which image finishes first.
  !===================================================================!

  pure subroutine combine(this, left, right, combined)

    class(reduction)       , intent(in)    :: this
    class(graph_functional), intent(in)    :: left
    class(graph_functional), intent(in)    :: right
    class(graph_functional), allocatable, intent(inout) :: combined

    real(dp)    :: a, b
    complex(dp) :: ca, cb
    logical     :: la, lb

    if (allocated(combined)) deallocate(combined)
    allocate(scalar_result :: combined)

    select case (this % rule)

    case (REDUCE_ALL, REDUCE_ANY)

       call get_logical(left, la)
       call get_logical(right, lb)
       if (this % rule == REDUCE_ALL) then
          call set_logical(combined, la .and. lb)
       else
          call set_logical(combined, la .or. lb)
       end if

    case (REDUCE_AVERAGE, REDUCE_NORM, REDUCE_COUNT)

       ! The running totals join; the division and the root still wait.
       select type (combined)
       type is (scalar_result)
          select type (left)
          type is (scalar_result)
             select type (right)
             type is (scalar_result)
                combined % tally  = left % tally  + right % tally
                combined % weight = left % weight + right % weight
             end select
          end select
       end select

    case (REDUCE_MINIMUM, REDUCE_MAXIMUM)

       call get_real(left, a)
       call get_real(right, b)
       if (this % rule == REDUCE_MINIMUM) then
          call set_real(combined, min(a, b))
       else
          call set_real(combined, max(a, b))
       end if

    case default

       ! Summing, on whichever road the parts travelled.
       if (left % value_kind() == GRAPH_FIELD_COMPLEX .or. &
            & right % value_kind() == GRAPH_FIELD_COMPLEX) then
          call get_complex(left, ca)
          call get_complex(right, cb)
          call set_complex(combined, ca + cb)
       else
          call get_real(left, a)
          call get_real(right, b)
          call set_real(combined, a + b)
       end if

    end select

  end subroutine combine

  !===================================================================!
  ! Finish, once, after every part has been joined in. This is where
  ! an average divides and a norm takes its root - and doing either
  ! any earlier is exactly the bug the four steps exist to prevent.
  !===================================================================!

  pure subroutine finalize(this, state, functional)

    class(reduction)       , intent(in) :: this
    class(graph_functional), intent(in) :: state
    class(graph_functional), allocatable, intent(inout) :: functional

    if (allocated(functional)) deallocate(functional)
    allocate(functional, source=state)

    select case (this % rule)

    case (REDUCE_AVERAGE)

       select type (state)
       type is (scalar_result)
          if (state % weight > 0.0_dp) then
             call set_real(functional, state % tally / state % weight)
          else
             call set_real(functional, 0.0_dp)
          end if
       end select

    case (REDUCE_NORM)

       select type (state)
       type is (scalar_result)
          if (state % tally > 0.0_dp) then
             call set_real(functional, state % tally**(1.0_dp / this % power))
          else
             call set_real(functional, 0.0_dp)
          end if
       end select

    case (REDUCE_COUNT)

       select type (state)
       type is (scalar_result)
          call functional % set_integer_vector([nint(state % tally)])
       end select

    end select

  end subroutine finalize

  !===================================================================!
  ! All four steps for a caller holding the whole graph.
  !
  ! Not pure. This is the one place a reduction spread
  ! across images is allowed to talk to the other images, and a
  ! distributed reduction would sum here before finalizing.
  !===================================================================!

  subroutine reduce(this, field, functional, measure)

    class(reduction)    , intent(in) :: this
    class(graph_field)  , intent(in) :: field
    class(graph_functional), allocatable, intent(inout) :: functional
    class(graph_field)  , intent(in), optional :: measure

    class(graph_functional), allocatable :: state

    call this % initialize(state)
    call this % accumulate(field, state, measure)
    call this % finalize(state, functional)

  end subroutine reduce

  !===================================================================!
  ! Fill every stored value of the field from the functional's one.
  ! Copy hands each value J; share hands each value J over the count
  ! of stored values, so a later sum returns J. A real J fills a
  ! real field; a complex J fills a complex field and carries a
  ! complex-step seed; any other kind fills zeros, following the
  ! value-kind rule the fields state.
  !===================================================================!

  pure subroutine broadcast_functional(this, functional, field)

    class(broadcast)       , intent(in)    :: this
    class(graph_functional), intent(in)    :: functional
    class(graph_field)     , intent(inout) :: field

    real(dp)    :: value
    complex(dp) :: complex_value
    integer     :: n, i

    n = field % num_entries() * max(field % num_components(), 1)

    if (functional % value_kind() == GRAPH_FIELD_COMPLEX) then
       call get_complex(functional, complex_value)
       if (this % rule == BROADCAST_SHARE .and. n > 0) then
          complex_value = complex_value / real(n, dp)
       end if
       call field % set_complex_vector([(complex_value, i = 1, n)])
    else
       call get_real(functional, value)
       if (this % rule == BROADCAST_SHARE .and. n > 0) then
          value = value / real(n, dp)
       end if
       call field % set_real_vector([(value, i = 1, n)])
    end if

  end subroutine broadcast_functional

  !===================================================================!
  ! One value in and out of any functional, through the contract's
  ! vector adapters. A wrong-kind read answers the zero of the asked
  ! kind, matching the adapters' zero-length signal.
  !===================================================================!

  pure subroutine get_real(f, x)

    class(graph_functional), intent(in) :: f
    real(dp), intent(out) :: x

    real(dp), allocatable :: t(:)

    call f % get_real_vector(t)
    if (size(t) >= 1) then
       x = t(1)
    else
       x = 0.0_dp
    end if

  end subroutine get_real

  pure subroutine set_real(f, x)

    class(graph_functional), intent(inout) :: f
    real(dp), intent(in) :: x

    call f % set_real_vector([x])

  end subroutine set_real

  pure subroutine get_complex(f, x)

    class(graph_functional), intent(in) :: f
    complex(dp), intent(out) :: x

    complex(dp), allocatable :: t(:)

    call f % get_complex_vector(t)
    if (size(t) >= 1) then
       x = t(1)
    else
       x = (0.0_dp, 0.0_dp)
    end if

  end subroutine get_complex

  pure subroutine set_complex(f, x)

    class(graph_functional), intent(inout) :: f
    complex(dp), intent(in) :: x

    call f % set_complex_vector([x])

  end subroutine set_complex

  pure subroutine get_logical(f, x)

    class(graph_functional), intent(in) :: f
    logical, intent(out) :: x

    logical, allocatable :: t(:)

    call f % get_logical_vector(t)
    if (size(t) >= 1) then
       x = t(1)
    else
       x = .false.
    end if

  end subroutine get_logical

  pure subroutine set_logical(f, x)

    class(graph_functional), intent(inout) :: f
    logical, intent(in) :: x

    call f % set_logical_vector([x])

  end subroutine set_logical

end module class_graph_reduction
