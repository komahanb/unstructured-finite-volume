!=====================================================================!
! The operation prime: the verb within a graph, (graph, fields) ->
! field. Three symbols. name says what it is; domain says which
! member set the result lives on and how many entries it has; apply
! does the work.
!
! apply writes its result into the output argument and never adds to
! what was there. The argument is intent(inout) only so a caller
! already holding a buffer of the right shape can lend it and save
! an allocation; lending changes the cost of the call, not its
! meaning.
!
! A concrete operation is handed the fields it reads when it is
! constructed - a coefficient, a measure, a geometry field arrives as
! an argument the compiler checks - so apply fetches nothing by name.
!
! ARGUMENTS. An operation F(x_1, ..., x_m) reads its inputs by
! position, and the positions form the operation's argument space.
! The space is declared once, by the constructor, with
! declare_arguments(m): m is the number of readable positions, not a
! required input length - an operation may accept fewer inputs than
! it declares, as it always could. An argument is an opaque ordinal
! in one operation's space; two arguments match only when they name
! the same position of the same space, so an argument of another
! operation can never stand for one of this operation's, whatever its
! position. Arguments are obtained from the operation that owns
! them, by argument(k); no caller constructs one.
!
! VARIATIONS. A partial directional derivative
!
!      D_{a_1} ... D_{a_k} F(x) [v_1, ..., v_k]
!
! is parameterized by the factors (a_j, v_j): one argument and one
! direction each. A variation is one such factor, and partial_action
! takes a list of them; a list of arguments beside a separate list of
! directions is not representable.
!
! Every operation reports max_degree, the highest order of exact
! partial action it computes (0 unless overridden), and
! partial_action, one mixed partial - differentiated once per
! variation, contracted against that variation's direction, on its
! own domain. The default stops the program, since an operation with
! max_degree 0 declares no partials.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_action

  use iso_fortran_env, only : dp => REAL64
  use view_directed , only : directed_graph
  use graph_fractal       , only : graph
  use field_calculus, only : field
  use field_stored  , only : stored_field
  use token_identity, only : token, next_token

  implicit none

  private

  public :: operation
  public :: argument
  public :: variation

  !===================================================================!
  ! One argument of one operation: the operation's argument space
  ! and a position in it. Both components are private; matches is
  ! the only comparison, and it is false across spaces.
  !===================================================================!

  type :: argument

     type(token), private :: space
     integer    , private :: ordinal = 0

   contains

     procedure :: matches  => argument_matches
     procedure :: is_named => argument_is_named

  end type argument

  !===================================================================!
  ! One differentiation factor: the argument differentiated and the
  ! direction the derivative is contracted against. The direction is
  ! held by value, a stored field, behind the accessor: a later
  ! borrowed representation changes no caller.
  !===================================================================!

  type :: variation

     type(argument)    , private :: wrt
     type(stored_field), private :: along

   contains

     procedure :: argument_is   => variation_argument_is
     procedure :: argument      => variation_argument
     procedure :: direction     => variation_direction
     procedure :: domain        => variation_domain
     procedure :: field         => variation_field
     procedure :: with_argument => variation_with_argument

  end type variation

  interface variation
     module procedure create_variation
  end interface variation

  type, abstract :: operation

     type(token), private :: arguments_space
     integer    , private :: declared_arguments = 0

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain
     procedure(operation_apply_interface) , deferred :: apply

     procedure :: max_degree     => operation_max_degree
     procedure :: partial_action => operation_partial_action

     procedure :: declare_arguments
     procedure :: num_arguments
     procedure :: argument => operation_argument
     procedure :: owns
     procedure :: require_owned

  end type operation

  abstract interface

     pure function operation_name_interface(this) result(name)
       import :: operation
       class(operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     !---------------------------------------------------------------!
     ! Where the result lives: WHICH set, and HOW MANY entries it
     ! has. The count travels beside the identity because every
     ! caller wants exactly those two things - to check the domain
     ! matches, and to size a field.
     !---------------------------------------------------------------!

     subroutine operation_domain_interface(this, input_graph, domain, &
          & num_entries)
       import :: operation, directed_graph, graph
       class(operation), intent(in)  :: this
       class(directed_graph)          , intent(in)  :: input_graph
       type(graph)       , intent(out) :: domain
       integer               , intent(out) :: num_entries
     end subroutine operation_domain_interface

     subroutine operation_apply_interface(this, input_graph, input_data, output)
       import :: operation, directed_graph, field
       class(operation), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       class(field), intent(in), optional :: input_data(:)
       class(field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

  end interface

contains

  !===================================================================!
  ! Same space, same position. An undeclared space matches nothing,
  ! including itself.
  !===================================================================!

  pure logical function argument_matches(this, other) result(same)

    class(argument), intent(in) :: this
    type(argument) , intent(in) :: other

    same = this % space % matches(other % space) .and. &
         & this % ordinal == other % ordinal

  end function argument_matches

  pure logical function argument_is_named(this) result(named)

    class(argument), intent(in) :: this

    named = this % space % declared() .and. this % ordinal > 0

  end function argument_is_named

  !===================================================================!
  ! Declare the argument space: mint it once, on the first call, and
  ! record how many positions are readable. A later call changes the
  ! count only, so arguments handed out earlier still belong to the
  ! same space. A negative count stops the program.
  !===================================================================!

  subroutine declare_arguments(this, n)

    class(operation), intent(inout) :: this
    integer         , intent(in)    :: n

    if (n < 0) then
       error stop 'operation: the argument count is nonnegative'
    end if

    if (.not. this % arguments_space % declared()) then
       this % arguments_space = next_token()
    end if

    this % declared_arguments = n

  end subroutine declare_arguments

  pure integer function num_arguments(this)

    class(operation), intent(in) :: this

    num_arguments = this % declared_arguments

  end function num_arguments

  !===================================================================!
  ! The k-th argument of this operation. An undeclared space or a
  ! position outside 1..num_arguments() stops the program: there is
  ! no such argument to name.
  !===================================================================!

  pure function operation_argument(this, k) result(a)

    class(operation), intent(in) :: this
    integer         , intent(in) :: k
    type(argument) :: a

    if (.not. this % arguments_space % declared()) then
       error stop 'operation: the argument space is declared before an argument is named'
    end if
    if (k < 1 .or. k > this % declared_arguments) then
       error stop 'operation: the argument is declared'
    end if

    a % space   = this % arguments_space
    a % ordinal = k

  end function operation_argument

  !===================================================================!
  ! Whether an argument belongs to this operation's space.
  !===================================================================!

  pure logical function owns(this, a)

    class(operation), intent(in) :: this
    type(argument)  , intent(in) :: a

    owns = this % arguments_space % declared() .and. &
         & a % space % matches(this % arguments_space) .and. &
         & a % ordinal >= 1 .and. a % ordinal <= this % declared_arguments

  end function owns

  !===================================================================!
  ! Refuse a variation list that names an argument of another
  ! operation. Every partial_action calls this first.
  !===================================================================!

  pure subroutine require_owned(this, variations)

    class(operation), intent(in) :: this
    type(variation) , intent(in) :: variations(:)

    integer :: j

    do j = 1, size(variations)
       if (.not. this % owns(variations(j) % wrt)) then
          error stop 'operation: a variation names an argument of the operation'
       end if
    end do

  end subroutine require_owned

  !===================================================================!
  ! Variations.
  !===================================================================!

  function create_variation(wrt, along) result(this)

    type(argument)    , intent(in) :: wrt
    type(stored_field), intent(in) :: along
    type(variation) :: this

    this % wrt   = wrt
    this % along = along

  end function create_variation

  pure logical function variation_argument_is(this, a)

    class(variation), intent(in) :: this
    type(argument)  , intent(in) :: a

    variation_argument_is = this % wrt % matches(a)

  end function variation_argument_is

  pure function variation_argument(this) result(a)

    class(variation), intent(in) :: this
    type(argument) :: a

    a = this % wrt

  end function variation_argument

  pure subroutine variation_direction(this, values)

    class(variation), intent(in)          :: this
    real(dp), allocatable, intent(out) :: values(:)

    call this % along % real_vector(values)

  end subroutine variation_direction

  function variation_domain(this) result(domain)

    class(variation), intent(in) :: this
    type(graph) :: domain

    domain = this % along % domain()

  end function variation_domain

  function variation_field(this) result(along)

    class(variation), intent(in) :: this
    type(stored_field) :: along

    along = this % along

  end function variation_field

  !===================================================================!
  ! The same direction on another argument: how an operation that
  ! wraps another restates a variation in the wrapped operation's
  ! argument space before delegating.
  !===================================================================!

  function variation_with_argument(this, a) result(moved)

    class(variation), intent(in) :: this
    type(argument)  , intent(in) :: a
    type(variation) :: moved

    moved % wrt   = a
    moved % along = this % along

  end function variation_with_argument

  !===================================================================!
  ! The default: no exact partial action of any order.
  !===================================================================!

  pure function operation_max_degree(this) result(degree)

    class(operation), intent(in) :: this
    integer :: degree

    associate (u1 => this); end associate

    degree = 0

  end function operation_max_degree

  !===================================================================!
  ! The default refuses every request, because max_degree is 0. A
  ! concrete type that declares a positive max_degree overrides both
  ! bindings; the order requested must not exceed its max_degree.
  !===================================================================!

  subroutine operation_partial_action(this, input_graph, input_data, &
       & variations, output)

    class(operation), intent(in)             :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in)                 :: input_data(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    associate (u1 => this, u2 => input_graph, u3 => input_data, &
         & u4 => variations); end associate
    if (allocated(output)) deallocate(output)

    error stop 'operation: the requested order is within max_degree'

  end subroutine operation_partial_action

end module operation_action
