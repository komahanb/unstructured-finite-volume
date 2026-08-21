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
! ARGUMENTS. An operation F(x_1, ..., x_m) declares its argument
! space once, by the constructor, with declare_arguments(m). An
! argument is an opaque ordinal in one operation's space; two
! arguments match only when they name the same position of the same
! space, so an argument of another operation can never stand for one
! of this operation's. Arguments are obtained from the operation that
! owns them, by argument(k); no caller constructs one.
!
! BINDINGS. A field and the argument it supplies are one object:
!
!      binding = (argument, field)
!
! the same pairing variation already makes of (argument, direction).
! A list of fields beside an implicit list of arguments is not the
! meaning; the bindings are.
!
! APPLICATIONS. A constituted application is
!
!      application = (operation identity, host carriers, bindings)
!
! It records the operation's argument space and the host's vertex and
! edge carrier identities - not the operation, not the host, and no
! pointer to either. It is validated at construction and complete by
! construction: every declared argument is bound, exactly once, by an
! argument this operation owns. There is no partial application and
! no arity to infer. A caller that means zero binds a zero field.
!
! THE GATE. apply and partial_action are concrete and NOT
! overridable. They check that the operation is the one the
! application was constituted for and that the host carries the
! recorded carriers and counts, and only then reach the operation's
! own act. Every concrete operation implements act and partial_act
! PRIVATELY, so nothing outside its module can reach past the gate;
! check_operation_contract.sh enforces that across both trees.
!
! The host check binds the carrier frame and the counts. It does not
! bind incidence, because the graph's representation is public; that
! is a defect of the graph, recorded separately, not a promise made
! here.
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
  public :: binding
  public :: application
  public :: constitute
  public :: reindex
  public :: bind_moved

  !-------------------------------------------------------------------!
  ! DIAGNOSTIC, temporary. constitutions counts applications built;
  ! payload_copies counts fields copied into a binding. Both are
  ! module state and belong to no design; they come out before any
  ! amendment is proposed.
  !-------------------------------------------------------------------!
  integer, public :: constitutions   = 0
  integer, public :: payload_copies  = 0

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

  !===================================================================!
  ! One supplied input: the argument it supplies and the field
  ! supplying it. The payload is class(field), not stored_field,
  ! because a functional is a lawful payload and does not extend
  ! stored_field. It is an owned snapshot: taken by value at
  ! construction, so a later write by the caller moves nothing here.
  !===================================================================!

  type :: binding

     type(argument)           , private :: to
     class(field), allocatable, private :: supplied

   contains

     procedure :: argument_of => binding_argument
     procedure :: argument_is => binding_argument_is

  end type binding

  interface binding
     module procedure create_binding
  end interface binding

  !===================================================================!
  ! One constituted application. It holds identities, not objects:
  ! the operation's argument space and the host's two carriers, so a
  ! copied host and a copied operation are both admitted and an
  ! unrelated graph is not. The bindings are owned.
  !===================================================================!

  type :: application

     type(token), private :: space
     type(graph), private :: vertices
     type(graph), private :: edges
     integer    , private :: nv = 0
     integer    , private :: ne = 0

     !----------------------------------------------------------------!
     ! An application either OWNS its bindings, or is a REINDEXED VIEW
     ! onto one that does. A view stores the map r from its own
     ! operation's arguments to the parent's:
     !
     !      view % field_for(a)  =  parent % field_for(r(a))
     !
     ! and copies nothing. It borrows for the length of one call and
     ! must not outlive its parent or be stored; reindex is reached
     ! only from an operation composing another, and the view is a
     ! local there.
     !----------------------------------------------------------------!
     type(binding), allocatable, private :: bound(:)

     type(application), pointer  , private :: parent => null()
     type(argument)   , allocatable, private :: outer(:)

   contains

     procedure :: field_for     => application_field_for
     procedure :: num_bindings  => application_num_bindings

  end type application

  type, abstract :: operation

     type(token), private :: arguments_space
     integer    , private :: declared_arguments = 0

   contains

     procedure(operation_name_interface)  , deferred :: name
     procedure(operation_domain_interface), deferred :: domain

     !----------------------------------------------------------------!
     ! What a concrete operation writes, and nothing outside its own
     ! module may call. The gate below is the only road in.
     !----------------------------------------------------------------!
     procedure(operation_act_interface), deferred, private :: act
     procedure, private :: partial_act => operation_partial_act

     !----------------------------------------------------------------!
     ! THE GATE. Concrete, not overridable, and the only caller of act.
     !----------------------------------------------------------------!
     procedure, non_overridable :: apply_empty
     procedure, non_overridable :: apply_fields
     procedure, non_overridable :: apply_application
     generic :: apply => apply_empty, apply_fields, apply_application

     procedure, non_overridable :: partial_action_fields
     procedure, non_overridable :: partial_action_application
     generic :: partial_action => partial_action_fields, &
          &                       partial_action_application

     procedure :: max_degree     => operation_max_degree

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

     !---------------------------------------------------------------!
     ! The operation's own work, on a constituted application. The
     ! application is a target so an operation may read a payload by
     ! pointer without copying it; that pointer must not outlive this
     ! call.
     !---------------------------------------------------------------!

     subroutine operation_act_interface(this, host, app, output)
       import :: operation, directed_graph, application, field
       class(operation)     , intent(in)         :: this
       class(directed_graph), intent(in)         :: host
       type(application)    , intent(in), target :: app
       class(field), allocatable, intent(inout)  :: output
     end subroutine operation_act_interface

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
  ! concrete type that declares a positive max_degree overrides this
  ! binding, privately; the order requested must not exceed its
  ! max_degree.
  !===================================================================!

  subroutine operation_partial_act(this, host, app, variations, output)

    class(operation)     , intent(in)         :: this
    class(directed_graph), intent(in)         :: host
    type(application)    , intent(in), target :: app
    type(variation)      , intent(in)         :: variations(:)
    class(field), allocatable, intent(inout)  :: output

    associate (u1 => this, u2 => host, u3 => app % nv, &
         & u4 => variations); end associate
    if (allocated(output)) deallocate(output)

    error stop 'operation: the requested order is within max_degree'

  end subroutine operation_partial_act

  !===================================================================!
  ! BINDINGS.
  !===================================================================!

  function create_binding(to, supplied) result(this)

    type(argument), intent(in) :: to
    class(field)  , intent(in) :: supplied
    type(binding) :: this

    this % to = to
    allocate(this % supplied, source=supplied)
    payload_copies = payload_copies + 1

  end function create_binding

  !===================================================================!
  ! OWNERSHIP TRANSFER. A field a caller has just created, and holds
  ! no other reference to, is MOVED into its binding rather than
  ! copied. The caller's variable is left deallocated, so there is
  ! exactly one owner afterwards and nothing to keep in step.
  !
  ! This is for freshly allocated payloads only. A field the caller
  ! still owns - a parameter it will use again - is bound by binding()
  ! and copied, as before.
  !===================================================================!

  subroutine bind_moved(to, supplied, this)

    type(argument), intent(in)                :: to
    class(field)  , allocatable, intent(inout) :: supplied
    type(binding) , intent(out)               :: this

    if (.not. allocated(supplied)) then
       error stop 'binding: the field moved into a binding is allocated'
    end if

    this % to = to
    call move_alloc(supplied, this % supplied)

  end subroutine bind_moved

  pure function binding_argument(this) result(a)

    class(binding), intent(in) :: this
    type(argument) :: a

    a = this % to

  end function binding_argument

  pure logical function binding_argument_is(this, a)

    class(binding), intent(in) :: this
    type(argument), intent(in) :: a

    binding_argument_is = this % to % matches(a)

  end function binding_argument_is

  !===================================================================!
  ! APPLICATIONS. The field bound to an argument, by identity and
  ! never by position, as a pointer: reading an input copies nothing.
  ! An application is complete by construction, so an argument this
  ! application does not bind is a defect and stops the program
  ! rather than answering null for a caller to branch on.
  !===================================================================!

  recursive function application_field_for(this, a) result(f)

    class(application), intent(in), target :: this
    type(argument)    , intent(in)         :: a
    class(field), pointer :: f

    integer :: k

    f => null()

    if (.not. a % space % matches(this % space)) then
       error stop 'application: the argument asked for is one this application binds'
    end if

    if (associated(this % parent)) then
       if (a % ordinal < 1 .or. a % ordinal > size(this % outer)) then
          error stop 'application: the argument asked for is one this application binds'
       end if
       f => this % parent % field_for(this % outer(a % ordinal))
       return
    end if

    do k = 1, size(this % bound)
       if (this % bound(k) % to % matches(a)) then
          f => this % bound(k) % supplied
          return
       end if
    end do

    error stop 'application: the argument asked for is one this application binds'

  end function application_field_for

  pure integer function application_num_bindings(this)

    class(application), intent(in), target :: this

    application_num_bindings = 0
    if (associated(this % parent)) then
       application_num_bindings = size(this % outer)
    else if (allocated(this % bound)) then
       application_num_bindings = size(this % bound)
    end if

  end function application_num_bindings

  !===================================================================!
  ! Constitute an application, and refuse anything less than a whole
  ! one. Each check stops the program:
  !
  !    foreign     a binding naming an argument of another operation
  !    duplicate   two bindings on one argument
  !    missing     a declared argument left unbound
  !    arity       a list longer than the argument list
  !
  ! The builder is consumed: the caller's list is moved in, not
  ! copied, and is left deallocated so no second reference survives.
  !===================================================================!

  subroutine constitute(of, host, bound, app)

    class(operation)     , intent(in)    :: of
    class(directed_graph), intent(in)    :: host
    type(binding), allocatable, intent(inout) :: bound(:)
    type(application)    , intent(out)   :: app

    integer :: j, k
    logical :: found

    if (.not. allocated(bound)) then
       error stop 'application: the binding list is given'
    end if
    if (associated(app % parent)) then
       error stop 'application: a reindexed view does not own bindings'
    end if

    do j = 1, size(bound)
       if (.not. of % owns(bound(j) % to)) then
          error stop 'application: a binding names an argument of the operation'
       end if
       do k = j + 1, size(bound)
          if (bound(j) % to % matches(bound(k) % to)) then
             error stop 'application: an argument is bound once'
          end if
       end do
    end do

    do k = 1, of % declared_arguments
       found = .false.
       do j = 1, size(bound)
          if (bound(j) % to % matches(of % argument(k))) found = .true.
       end do
       if (.not. found) then
          error stop 'application: every declared argument is bound'
       end if
    end do

    if (size(bound) /= of % declared_arguments) then
       error stop 'application: the binding list is exactly the argument list'
    end if

    constitutions = constitutions + 1

    app % space    = of % arguments_space
    app % vertices = host % vertex_set()
    app % edges    = host % edge_set()
    app % nv       = host % num_vertices()
    app % ne       = host % num_edges()

    call move_alloc(bound, app % bound)

  end subroutine constitute

  !===================================================================!
  ! A REINDEXED VIEW. An operation composing another states, for each
  ! of the inner operation's arguments, which of its own supplies it.
  ! Nothing is copied: the view answers through the parent. Checks,
  ! each stopping the program: the map is exactly the inner argument
  ! list, and every entry names an argument of the parent's space.
  !
  ! The view borrows. It is valid while the parent lives, must stay
  ! local to the call that made it, and cannot be constituted into an
  ! owning application afterwards.
  !===================================================================!

  subroutine reindex(of, host, parent, outer, view)

    class(operation)     , intent(in)         :: of
    class(directed_graph), intent(in)         :: host
    type(application)    , intent(in), target :: parent
    type(argument)       , intent(in)         :: outer(:)
    type(application)    , intent(out)        :: view

    integer :: k

    if (size(outer) /= of % declared_arguments) then
       error stop 'application: the reindex map is exactly the argument list'
    end if

    do k = 1, size(outer)
       if (.not. outer(k) % space % matches(parent % space)) then
          error stop 'application: a reindex entry names an argument of the parent'
       end if
    end do

    view % space    = of % arguments_space
    view % vertices = host % vertex_set()
    view % edges    = host % edge_set()
    view % nv       = host % num_vertices()
    view % ne       = host % num_edges()

    view % parent => parent
    view % outer  =  outer

  end subroutine reindex

  !===================================================================!
  ! THE GATE. The operation must be the one the application was
  ! constituted for, and the host must carry the recorded carriers
  ! and counts. Only then does the operation's own act run.
  !===================================================================!

  subroutine require_admissible(this, host, app)

    class(operation)     , intent(in) :: this
    class(directed_graph), intent(in) :: host
    type(application)    , intent(in), target :: app

    type(graph) :: v, e

    if (.not. this % arguments_space % matches(app % space)) then
       error stop 'operation: the application was constituted for this operation'
    end if

    v = host % vertex_set()
    e = host % edge_set()

    if (.not. app % vertices % same_as(v)) then
       error stop 'operation: the host carries the application''s vertex set'
    end if
    if (.not. app % edges % same_as(e)) then
       error stop 'operation: the host carries the application''s edge set'
    end if
    if (host % num_vertices() /= app % nv .or. &
         & host % num_edges() /= app % ne) then
       error stop 'operation: the host carries the application''s counts'
    end if

  end subroutine require_admissible

  subroutine apply_application(this, host, app, output)

    class(operation)     , intent(in)         :: this
    class(directed_graph), intent(in)         :: host
    type(application)    , intent(in), target :: app
    class(field), allocatable, intent(inout)  :: output

    call require_admissible(this, host, app)
    call this % act(host, app, output)

  end subroutine apply_application

  subroutine partial_action_application(this, host, app, variations, output)

    class(operation)     , intent(in)         :: this
    class(directed_graph), intent(in)         :: host
    type(application)    , intent(in), target :: app
    type(variation)      , intent(in)         :: variations(:)
    class(field), allocatable, intent(inout)  :: output

    call require_admissible(this, host, app)
    call this % require_owned(variations)
    call this % partial_act(host, app, variations, output)

  end subroutine partial_action_application

  !===================================================================!
  ! A zero-argument operation is applied on the empty binding list.
  !===================================================================!

  subroutine apply_empty(this, host, output)

    class(operation)     , intent(in) :: this
    class(directed_graph), intent(in) :: host
    class(field), allocatable, intent(inout) :: output

    type(binding), allocatable :: bound(:)
    type(application) :: app

    allocate(bound(0))
    call constitute(this, host, bound, app)
    call this % apply_application(host, app, output)

  end subroutine apply_empty

  !===================================================================!
  ! THE POSITIONAL ADAPTER. A caller that still holds an ordered
  ! field array hands it over here; position k supplies argument k,
  ! and the list must be exactly the argument list, as everywhere
  ! else. It exists so callers may migrate one at a time, and it is
  ! deleted when the last of them has.
  !===================================================================!

  subroutine apply_fields(this, host, input_data, output)

    class(operation)     , intent(in) :: this
    class(directed_graph), intent(in) :: host
    class(field)         , intent(in) :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(binding), allocatable :: bound(:)
    type(application) :: app

    call bindings_of(this, input_data, bound)
    call constitute(this, host, bound, app)
    call this % apply_application(host, app, output)

  end subroutine apply_fields

  subroutine partial_action_fields(this, host, input_data, variations, output)

    class(operation)     , intent(in) :: this
    class(directed_graph), intent(in) :: host
    class(field)         , intent(in) :: input_data(:)
    type(variation)      , intent(in) :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(binding), allocatable :: bound(:)
    type(application) :: app

    call bindings_of(this, input_data, bound)
    call constitute(this, host, bound, app)
    call this % partial_action_application(host, app, variations, output)

  end subroutine partial_action_fields

  subroutine bindings_of(this, input_data, bound)

    class(operation), intent(in) :: this
    class(field)    , intent(in) :: input_data(:)
    type(binding), allocatable, intent(out) :: bound(:)

    integer :: k

    if (size(input_data) /= this % declared_arguments) then
       error stop 'operation: the input list is exactly the argument list'
    end if

    allocate(bound(size(input_data)))
    do k = 1, size(input_data)
       bound(k) = binding(this % argument(k), input_data(k))
    end do

  end subroutine bindings_of

end module operation_action
