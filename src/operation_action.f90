!=====================================================================!
! The operation prime: the map within a graph, (graph, fields) ->
! field. Three symbols. name states the operation's name; domain
! states which member set the result is defined on and how many
! entries it has; apply computes the result.
!
! apply writes its result into the output argument and never adds to
! what was there. The argument is intent(inout) only so a caller
! already storing a buffer of the right shape can pass it and avoid
! an allocation; passing a buffer changes the cost of the call, not
! its result.
!
! A concrete operation receives the fields it reads at construction
! - a coefficient, a measure, a geometry field is passed as an
! argument the compiler checks - so apply retrieves nothing by name.
!
! ARGUMENTS. An operation F(x_1, ..., x_m) declares an argument
! space. The declaration has m slots; the slots supply operation-owned
! argument identities, not caller-owned names. An argument is an
! opaque ordinal in one operation's space; two arguments match only
! when they name the same position of the same space, so an argument
! of another operation can never denote one of this operation's,
! whatever its position. Arguments are obtained from the operation
! that owns them, by argument(k); no caller constructs one.
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
! partial action it computes (0 unless declared), and
! partial_action, one mixed partial - differentiated once per
! variation, contracted against that variation's direction, on its
! own domain. The default stops the program, since an operation with
! max_degree 0 declares no partials.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_action

  use util_precision  , only : dp
  use view_directed , only : directed_graph
  use graph_fractal       , only : graph
  use field_calculus, only : field, FIELD_NONE
  use field_stored  , only : stored_field, typed_field_domain
  use token_identity, only : token, next_token
  use util_derivative_terms, only : derivative_terms

  implicit none

  private

  public :: operation
  public :: emit, emit_real
  public :: seeded_argument
  public :: argument
  public :: contract
  public :: variation
  public :: binding, moved_binding
  public :: is_bound, bound_on, bound_value, bound_real_vector, bound_integer_vector
  public :: applied, varied
  public :: design_partial, jacobian_of
  public :: dense_of_triples

  !===================================================================!
  ! One argument of one operation: the operation's argument space
  ! and a position in it. Both components are private; matches is
  ! the only comparison, and it is false across spaces.
  !===================================================================!

  !===================================================================!
  ! What a bound field must be: one of the value kinds, and the
  ! component count when the operation fixes it. A count of zero
  ! leaves the count open for an operation that reads the shape from
  ! the field.
  !===================================================================!

  type :: contract

     integer, allocatable, private :: value_kinds(:)
     integer, private :: components = 0

   contains

     procedure :: accepts => contract_accepts

  end type contract

  interface contract
     module procedure create_contract
     module procedure create_contract_set
  end interface contract

  type :: argument

     type(token), private :: space
     integer    , private :: ordinal = 0
     type(contract), private :: required

   contains

     procedure :: matches  => argument_matches
     procedure :: is_named => argument_is_named
     procedure :: contract => argument_contract

  end type argument

  !===================================================================!
  ! One differentiation factor: the argument differentiated and the
  ! direction the derivative is contracted against. The direction is
  ! stored by value, as a stored field, behind the accessor: a later
  ! non-owning representation changes no caller.
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

  !===================================================================!
  ! One argument bound to one field. The binding owns the field and
  ! checks the argument's contract at construction; an apply reads
  ! its inputs by argument identity from an array of these. The
  ! driver binds a rule's argument k to the datum at the k-th vertex
  ! the rule reads; a direct caller binds by position through bind.
  !===================================================================!

  type :: binding

     type(argument)            , private :: to
     class(field), allocatable , private :: value

   contains

     procedure :: argument_is => binding_argument_is

  end type binding

  interface binding
     module procedure create_binding
  end interface binding

  type, abstract :: operation

     type(token), private :: arguments_space
     integer    , private :: declared_arguments = 0
     type(contract), allocatable, private :: argument_contracts(:)

     ! the name reported, and the highest order of exact partial
     ! action, both declared with the arguments; a concretion that
     ! computes either overrides the accessor
     character(len=:), allocatable, private :: label
     integer, private :: exact_degree = 0

     ! THE VERSION. A statement that is unchanged between two solves
     ! has the same version, and a direct solver retains its factors
     ! while the version it last factorised is unchanged. Zero is no
     ! version: the default, and always factorised again.
     integer    , private :: version_number = 0
     logical    , private :: is_transposed = .false.

   contains

     procedure :: name   => operation_name
     procedure :: domain => operation_domain
     procedure(operation_apply_interface) , deferred :: apply

     procedure :: max_degree       => operation_max_degree
     procedure :: defined_at_zero  => operation_defined_at_zero
     procedure :: partial_action   => operation_partial_action
     procedure :: value_by_partial_action
     procedure :: explicit_tangent => operation_explicit_tangent

     procedure :: declare_arguments
     procedure :: versioned
     procedure :: version
     procedure :: transpose_version
     procedure :: num_arguments
     procedure :: argument => operation_argument
     procedure :: contracts
     procedure :: owns
     procedure :: require_owned
     procedure :: require_variations
     procedure, private :: bind_fields, bind_bindings
     generic   :: bind => bind_fields, bind_bindings

  end type operation

  abstract interface

     !---------------------------------------------------------------!
     ! The domain of the result: WHICH set, and HOW MANY entries it
     ! has. The count is returned beside the identity because every
     ! caller requires exactly those two quantities - to check the
     ! domain matches, and to size a field.
     !---------------------------------------------------------------!


     subroutine operation_apply_interface(this, input_graph, inputs, output)
       import :: operation, directed_graph, field, binding
       class(operation), intent(in) :: this
       class(directed_graph), intent(in) :: input_graph
       type(binding), intent(in), optional :: inputs(:)
       class(field), allocatable, intent(inout) :: output
     end subroutine operation_apply_interface

  end interface

contains

  pure function create_contract(value_kind, components) result(this)

    integer, intent(in)           :: value_kind
    integer, intent(in), optional :: components
    type(contract) :: this

    if (value_kind == FIELD_NONE) then
       error stop 'operation: a contract requires a value kind'
    end if

    this % value_kinds = [value_kind]
    call fix_components(this, components)

  end function create_contract

  pure function create_contract_set(value_kinds, components) result(this)

    integer, intent(in)           :: value_kinds(:)
    integer, intent(in), optional :: components
    type(contract) :: this

    if (size(value_kinds) < 1 .or. any(value_kinds == FIELD_NONE)) then
       error stop 'operation: a contract requires value kinds'
    end if

    this % value_kinds = value_kinds
    call fix_components(this, components)

  end function create_contract_set

  pure subroutine fix_components(this, components)

    type(contract), intent(inout)        :: this
    integer       , intent(in), optional :: components

    this % components = 0
    if (present(components)) then
       if (components < 1) then
          error stop 'operation: a contract requires a positive component count'
       end if
       this % components = components
    end if

  end subroutine fix_components

  pure logical function contract_accepts(this, value) result(accepted)

    class(contract), intent(in) :: this
    class(field)   , intent(in) :: value

    if (.not. allocated(this % value_kinds)) then
       accepted = .true.
       return
    end if

    accepted = any(value % value_kind() == this % value_kinds)
    if (this % components > 0) then
       accepted = accepted .and. value % num_components() == this % components
    end if

  end function contract_accepts

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

  pure logical function argument_is_named(this) result(is_defined)

    class(argument), intent(in) :: this

    is_defined = this % space % declared() .and. this % ordinal > 0

  end function argument_is_named

  pure function argument_contract(this) result(required)

    class(argument), intent(in) :: this
    type(contract) :: required

    if (.not. this % is_named()) then
       error stop 'operation: a contract belongs to a named argument'
    end if

    required = this % required

  end function argument_contract

  !===================================================================!
  ! Declare the argument space: allocate its token once, on the first
  ! call, and record how many positions are readable. A later call
  ! changes the count only, so arguments returned earlier still
  ! belong to the same space. A negative count stops the program.
  ! The label and the highest exact degree are recorded when given
  ! and retained otherwise.
  !===================================================================!

  subroutine declare_arguments(this, n, contracts, label, max_degree)

    class(operation), intent(inout) :: this
    integer         , intent(in)    :: n
    type(contract)  , intent(in), optional :: contracts(:)
    character(len=*), intent(in), optional :: label
    integer         , intent(in), optional :: max_degree

    if (n < 0) then
       error stop 'operation: the argument count is nonnegative'
    end if
    if (present(contracts)) then
       if (size(contracts) /= n) then
          error stop 'operation: every argument has one contract'
       end if
    end if

    if (.not. this % arguments_space % declared()) then
       this % arguments_space = next_token()
    end if

    this % declared_arguments = n
    if (allocated(this % argument_contracts)) deallocate(this % argument_contracts)
    allocate(this % argument_contracts(n))
    if (present(contracts))  this % argument_contracts = contracts
    if (present(label))      this % label = label
    if (present(max_degree)) this % exact_degree = max_degree

  end subroutine declare_arguments

  pure integer function num_arguments(this)

    class(operation), intent(in) :: this

    num_arguments = this % declared_arguments

  end function num_arguments

  !===================================================================!
  ! The contracts in argument order, for an operation that declares
  ! the argument space of one it composes.
  !===================================================================!

  pure function contracts(this) result(required)

    class(operation), intent(in) :: this
    type(contract), allocatable :: required(:)

    allocate(required(this % declared_arguments))
    if (allocated(this % argument_contracts)) required = this % argument_contracts

  end function contracts

  !===================================================================!
  ! The declared label, or 'operation' when none was declared.
  !===================================================================!

  pure function operation_name(this) result(name)

    class(operation), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'operation'
    if (allocated(this % label)) name = this % label

  end function operation_name

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
    a % required = this % argument_contracts(k)

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
  ! Reject a variation list that names an argument of another
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
  ! Reject a variation list that names an argument of another
  ! operation, or that is longer than the highest exact degree. A
  ! partial_action computed over derivative terms calls this first.
  !===================================================================!

  pure subroutine require_variations(this, variations)

    class(operation), intent(in) :: this
    type(variation) , intent(in) :: variations(:)

    call this % require_owned(variations)
    if (size(variations) > this % max_degree()) then
       error stop 'operation: the requested order is within max_degree'
    end if

  end subroutine require_variations

  !===================================================================!
  ! Argument k read from the bindings as derivative terms over the
  ! variations: entry j is x_j with direction i set to the i-th
  ! variation's j-th entry for every variation naming argument k, and
  ! consumed counts those variations. A direction of another length
  ! than the argument stops the program. The caller rejects a
  ! variation naming another argument by comparing consumed with the
  ! count of variations.
  !===================================================================!

  subroutine seeded_argument(this, inputs, variations, k, terms, consumed)

    class(operation), intent(in) :: this
    type(binding)   , intent(in) :: inputs(:)
    type(variation) , intent(in) :: variations(:)
    integer         , intent(in) :: k
    type(derivative_terms), allocatable, intent(out) :: terms(:)
    integer         , intent(out) :: consumed

    real(dp), allocatable :: x(:), v(:)
    type(argument) :: a
    integer :: i, j

    a = this % argument(k)
    call bound_real_vector(inputs, a, x)

    allocate(terms(size(x)))
    do j = 1, size(x)
       terms(j) = derivative_terms(x(j), size(variations))
    end do

    consumed = 0
    do i = 1, size(variations)
       if (.not. variations(i) % argument_is(a)) cycle
       call variations(i) % direction(v)
       if (size(v) /= size(x)) then
          error stop 'operation: a direction has one entry per value of the argument it varies'
       end if
       do j = 1, size(x)
          call terms(j) % set_direction(i, v(j))
       end do
       consumed = consumed + 1
    end do

  end subroutine seeded_argument

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
  ! The highest exact degree declared with the arguments: zero, no
  ! exact partial action of any order, unless declared.
  !===================================================================!

  pure function operation_max_degree(this) result(degree)

    class(operation), intent(in) :: this
    integer :: degree

    degree = this % exact_degree

  end function operation_max_degree

  !===================================================================!
  ! WHETHER THE ZERO STATE LIES IN THE OPERATION'S DOMAIN. A(0) is the
  ! constant part a linear solver subtracts to read the linear part of
  ! an affine statement, and a minimizer evaluates it when the
  ! operation is stated. An operation built from sums, products and
  ! non-negative powers has a value there, and the default declares
  ! so. An operation that divides by, or takes a negative power, a
  ! logarithm or a square root of, a quantity read from the state
  ! overrides this: the zero state is outside its domain and is never
  ! evaluated.
  !===================================================================!

  pure logical function operation_defined_at_zero(this) result(defined)

    class(operation), intent(in) :: this

    associate (u1 => this); end associate
    defined = .true.

  end function operation_defined_at_zero

  !===================================================================!
  ! The default rejects every request, because max_degree is 0. A
  ! concrete type that declares a positive max_degree overrides both
  ! bindings; the order requested must not exceed its max_degree.
  !===================================================================!

  subroutine versioned(this, version_number, transposed)

    class(operation), intent(inout) :: this
    integer         , intent(in)    :: version_number
    logical         , intent(in), optional :: transposed

    this % version_number   = version_number
    this % is_transposed = .false.
    if (present(transposed)) this % is_transposed = transposed

  end subroutine versioned

  !===================================================================!
  ! Whether the versioned statement is the transpose of the one the
  ! version names: read from its pattern where the statement is
  ! formed, and stored with the version so that a solver storing the
  ! factors of the one substitutes them transposed for the other.
  !===================================================================!

  pure logical function transpose_version(this)

    class(operation), intent(in) :: this

    transpose_version = this % is_transposed

  end function transpose_version

  pure integer function version(this) result(version_number)

    class(operation), intent(in) :: this

    version_number = this % version_number

  end function version

  !===================================================================!
  ! THE EXPLICIT TANGENT. A statement that can express its own
  ! tangent in one argument as triples - row, column, weight -
  ! reports so here, and a minimizer governing it may then state the
  ! explicit operator instead of forming the tangent by matvecs. The
  ! default tangent_defined is false; the arrays
  ! are then not assigned. Nothing here is a matvec: a statement that
  ! reports its tangent explicitly stores its own structure.
  !===================================================================!

  subroutine operation_explicit_tangent(this, input_graph, inputs, which, &
       & rows, columns, weights, tangent_defined)

    class(operation)     , intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(binding)        , intent(in)  :: inputs(:)
    integer              , intent(in)  :: which
    integer , allocatable, intent(out) :: rows(:), columns(:)
    real(dp), allocatable, intent(out) :: weights(:)
    logical              , intent(out) :: tangent_defined

    associate (u1 => this, u2 => input_graph, u3 => inputs, u4 => which); end associate
    allocate(rows(0), columns(0), weights(0))
    tangent_defined = .false.

  end subroutine operation_explicit_tangent

  subroutine operation_partial_action(this, input_graph, inputs, &
       & variations, output)

    class(operation), intent(in)             :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding), intent(in)                :: inputs(:)
    type(variation), intent(in)              :: variations(:)
    class(field), allocatable, intent(inout) :: output

    associate (u1 => this, u2 => input_graph, u3 => inputs, &
         & u4 => variations); end associate
    if (allocated(output)) deallocate(output)

    error stop 'operation: the requested order is within max_degree'

  end subroutine operation_partial_action

  !===================================================================!
  ! The value as the partial action of order zero, for an operation
  ! computed over derivative terms: apply is this call. Absent inputs
  ! stop the program.
  !===================================================================!

  subroutine value_by_partial_action(this, input_graph, inputs, output)

    class(operation)     , intent(in)        :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding)        , intent(in), optional :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(variation) :: none(0)

    if (.not. present(inputs)) then
       error stop 'operation: the arguments are bound'
    end if

    call this % partial_action(input_graph, inputs, none, output)

  end subroutine value_by_partial_action

  subroutine applied(action, context, inputs, y)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: context
    type(stored_field)   , intent(in) :: inputs(:)
    real(dp), allocatable, intent(out) :: y(:)

    class(field), allocatable :: out

    call action % apply(context, action % bind(inputs), out)
    call out % real_vector(y)

  end subroutine applied

  subroutine varied(action, context, inputs, which, domain, v, y, which2, domain2, v2)

    class(operation)     , intent(in) :: action
    class(directed_graph), intent(in) :: context
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: which
    type(graph)          , intent(in) :: domain
    real(dp)             , intent(in) :: v(:)
    real(dp), allocatable, intent(out) :: y(:)
    integer    , intent(in), optional :: which2
    type(graph), intent(in), optional :: domain2
    real(dp)   , intent(in), optional :: v2(:)

    type(stored_field) :: direction, second
    type(typed_field_domain) :: first_domain, second_domain
    class(field), allocatable :: out

    first_domain = typed_field_domain(domain, size(v))
    direction    = first_domain % direction(v)

    if (present(which2)) then
       second_domain = typed_field_domain(domain2, size(v2))
       second        = second_domain % direction(v2)
       call action % partial_action(context, action % bind(inputs), &
            & [variation(action % argument(which), direction), &
            &  variation(action % argument(which2), second)], out)
    else
       call action % partial_action(context, action % bind(inputs), &
            & [variation(action % argument(which), direction)], out)
    end if

    call out % real_vector(y)

  end subroutine varied

  subroutine design_partial(rows, unknowns, inputs, n, design_domain, d)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: n
    type(graph)          , intent(in) :: design_domain
    real(dp), allocatable, intent(out) :: d(:)

    call varied(rows, unknowns, inputs, 2, design_domain, spread(1.0_dp, 1, n), d)

  end subroutine design_partial

  subroutine jacobian_of(rows, unknowns, inputs, num_unknowns, state_domain, a)

    class(operation)     , intent(in) :: rows
    class(directed_graph), intent(in) :: unknowns
    type(stored_field)   , intent(in) :: inputs(:)
    integer              , intent(in) :: num_unknowns
    type(graph)          , intent(in) :: state_domain
    real(dp), allocatable, intent(out) :: a(:,:)

    real(dp), allocatable :: v(:), column(:), w(:)
    integer , allocatable :: r(:), c(:)
    logical :: tangent_defined
    integer :: j

    call rows % explicit_tangent(unknowns, rows % bind(inputs), 1, r, c, w, tangent_defined)
    if (tangent_defined) then
       call dense_of_triples(num_unknowns, r, c, w, a)
       return
    end if

    allocate(a(num_unknowns, num_unknowns), source=0.0_dp)
    allocate(v(num_unknowns), source=0.0_dp)
    do j = 1, num_unknowns
       v    = 0.0_dp
       v(j) = 1.0_dp
       call varied(rows, unknowns, inputs, 1, state_domain, v, column)
       a(:, j) = column
    end do

  end subroutine jacobian_of

  !===================================================================!
  ! The dense n by n matrix of (row, column, weight) triples: equal
  ! pairs sum, in the order given.
  !===================================================================!

  pure subroutine dense_of_triples(n, rows, columns, weights, a)

    integer , intent(in) :: n, rows(:), columns(:)
    real(dp), intent(in) :: weights(:)
    real(dp), allocatable, intent(out) :: a(:,:)

    integer :: e

    allocate(a(n, n), source=0.0_dp)
    do e = 1, size(rows)
       a(rows(e), columns(e)) = a(rows(e), columns(e)) + weights(e)
    end do

  end subroutine dense_of_triples

  !===================================================================!
  ! Bindings.
  !===================================================================!

  function create_binding(to, value) result(this)

    type(argument), intent(in) :: to
    class(field)  , intent(in) :: value
    type(binding) :: this

    call require_bindable(to, value)
    this % to = to
    allocate(this % value, source=value)

  end function create_binding

  !===================================================================!
  ! The same binding taking ownership of the field rather than
  ! copying it; the driver passes each datum it read to the binding
  ! this way.
  !===================================================================!

  function moved_binding(to, stored) result(this)

    type(argument)           , intent(in)    :: to
    class(field), allocatable, intent(inout) :: stored
    type(binding) :: this

    if (.not. allocated(stored)) then
       error stop 'operation: a binding contains a field'
    end if
    call require_bindable(to, stored)
    this % to = to
    call move_alloc(stored, this % value)

  end function moved_binding

  !===================================================================!
  ! An unnamed argument, or a field outside the argument's contract,
  ! stops the program.
  !===================================================================!

  subroutine require_bindable(to, value)

    type(argument), intent(in) :: to
    class(field)  , intent(in) :: value
    type(contract) :: required

    if (.not. to % is_named()) then
       error stop 'operation: a binding names an argument'
    end if
    required = to % contract()
    if (.not. required % accepts(value)) then
       error stop 'operation: a bound field satisfies its argument contract'
    end if

  end subroutine require_bindable

  pure logical function binding_argument_is(this, a)

    class(binding), intent(in) :: this
    type(argument), intent(in) :: a

    binding_argument_is = this % to % matches(a)

  end function binding_argument_is

  !===================================================================!
  ! Bind fields to this operation's arguments by position: field k
  ! to argument k. Fewer fields than arguments leave the rest unbound;
  ! more than declared is an error.
  !===================================================================!

  function bind_fields(this, fields) result(bound)

    class(operation), intent(in) :: this
    class(field)    , intent(in) :: fields(:)
    type(binding), allocatable :: bound(:)

    integer :: k

    if (size(fields) > this % num_arguments()) then
       error stop 'operation: every bound field names a declared argument'
    end if

    allocate(bound(size(fields)))
    do k = 1, size(fields)
       bound(k) = binding(this % argument(k), fields(k))
    end do

  end function bind_fields

  !===================================================================!
  ! Rebind another operation's bindings to this operation's arguments
  ! by position, for an operation that passes its inputs to one it
  ! composes.
  !===================================================================!

  function bind_bindings(this, others) result(bound)

    class(operation), intent(in) :: this
    type(binding)   , intent(in) :: others(:)
    type(binding), allocatable :: bound(:)

    integer :: k

    if (size(others) > this % num_arguments()) then
       error stop 'operation: every bound field names a declared argument'
    end if

    allocate(bound(size(others)))
    do k = 1, size(others)
       if (.not. allocated(others(k) % value)) then
          error stop 'operation: a binding contains a field'
       end if
       bound(k) = binding(this % argument(k), others(k) % value)
    end do

  end function bind_bindings

  !===================================================================!
  ! The position of an argument's binding, zero when it has none. An
  ! unnamed argument or a repeated binding stops the program.
  !===================================================================!

  pure integer function binding_position(bound, a) result(found)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a

    integer :: k

    if (.not. a % is_named()) then
       error stop 'operation: a bound value is named by an argument'
    end if

    found = 0
    do k = 1, size(bound)
       if (bound(k) % argument_is(a)) then
          if (found /= 0) then
             error stop 'operation: an argument is bound once'
          end if
          found = k
       end if
    end do

  end function binding_position

  !===================================================================!
  ! The one binding of an argument: absent stops the program.
  !===================================================================!

  integer function bound_index(bound, a) result(found)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a

    found = binding_position(bound, a)
    if (found == 0) then
       error stop 'operation: the argument is bound'
    end if
    if (.not. allocated(bound(found) % value)) then
       error stop 'operation: a binding contains a field'
    end if

  end function bound_index

  !===================================================================!
  ! Whether an argument has a binding: the driver leaves unbound the
  ! argument of a vertex that no operation has written.
  !===================================================================!

  pure logical function is_bound(bound, a)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a

    integer :: k

    k = binding_position(bound, a)
    is_bound = .false.
    if (k > 0) is_bound = allocated(bound(k) % value)

  end function is_bound

  !===================================================================!
  ! Whether the bound field of an argument is defined on the domain
  ! with the stated number of values (entries times components): the
  ! check read from the binding itself, without a copy of the field.
  !===================================================================!

  logical function bound_on(bound, a, domain, num_values) result(defined)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a
    type(graph)   , intent(in) :: domain
    integer       , intent(in) :: num_values

    integer :: found

    found   = bound_index(bound, a)
    defined = bound(found) % value % defined_on(domain) .and. &
         & bound(found) % value % num_entries() * bound(found) % value % num_components() == num_values

  end function bound_on

  subroutine bound_value(bound, a, value)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a
    class(field), allocatable, intent(inout) :: value

    integer :: found

    found = bound_index(bound, a)
    if (allocated(value)) deallocate(value)
    allocate(value, source=bound(found) % value)

  end subroutine bound_value

  subroutine bound_real_vector(bound, a, values)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a
    real(dp), allocatable, intent(out) :: values(:)

    integer :: found

    found = bound_index(bound, a)
    call bound(found) % value % real_vector(values)

  end subroutine bound_real_vector

  subroutine bound_integer_vector(bound, a, values)

    type(binding) , intent(in) :: bound(:)
    type(argument), intent(in) :: a
    integer, allocatable, intent(out) :: values(:)

    integer :: found

    found = bound_index(bound, a)
    call bound(found) % value % integer_vector(values)

  end subroutine bound_integer_vector

  !===================================================================!
  ! The domain of an operation's result, unless the operation
  ! overrides it: one entry per vertex of the graph it is applied on.
  ! An operation whose result is defined on its edges, or on a domain
  ! of its own, overrides this.
  !===================================================================!

  subroutine operation_domain(this, input_graph, domain, num_entries)

    class(operation)     , intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % vertex_set()
    num_entries = input_graph % num_vertices()

  end subroutine operation_domain

  !===================================================================!
  ! The one way an operation returns its result: a supplied buffer
  ! is overwritten, never added to.
  !===================================================================!

  subroutine emit(out, output)

    class(field)             , intent(in)    :: out
    class(field), allocatable, intent(inout) :: output

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine emit

  !===================================================================!
  ! Real values emitted as a stored field named name on the domain,
  ! with n entries of num_components components each.
  !===================================================================!

  subroutine emit_real(name, domain, n, values, output, num_components)

    character(len=*), intent(in) :: name
    type(graph)     , intent(in) :: domain
    integer         , intent(in) :: n
    real(dp)        , intent(in) :: values(:)
    class(field), allocatable, intent(inout) :: output
    integer, intent(in), optional :: num_components

    type(stored_field) :: out

    out = stored_field(name, domain, n, num_components)
    call out % set_real_vector(values)
    call emit(out, output)

  end subroutine emit_real

end module operation_action
