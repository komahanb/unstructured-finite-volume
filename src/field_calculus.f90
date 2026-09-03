!=====================================================================!
! LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS
!
! The level defines one concept: WHAT VALUES ARE DEFINED ON A DOMAIN. A
! field is a function over one finite domain,
!
!      f : A -> V         or        f : S -> V,   S c--> A
!
! and that domain is named by a SET GRAPH: one identity, never a
! union, never a side flag (AGENTS.md 20, CALCULATOR.md 12). A field
! needs a domain; it does not need a graph container.
!
!             WHAT A FIELD STORES, AND WHAT IT REQUIRES
!
! Identity and a frozen count, and nothing else about the domain:
!
!      domain()        WHICH set - a set graph, by value
!      num_entries()   HOW MANY - the count taken at construction
!
! Those are the only two queries every field caller makes, and
! neither is a query about membership, so neither needs a map.
! A caller that requires WHICH members belong, WHERE a member is
! positioned, or WHAT the domain is called reads the set map, the
! inclusion map or the label map - explicitly, with those maps in
! scope, at the call site.
!
! The count is frozen because it always was: a field has stored a COPY
! of its domain since the first version, so num_entries never tracked
! later mutation of the caller's set. Storing the integer states
! explicitly what copying implied, and stops N fields on one domain
! from storing N copies of its extension.
!
! This module is the RELOCATED field ontology: the one abstract
! field and its value-kind constants. They came from the old
! grammar, which re-exported them for a period and was emptied and
! deleted in PR2; every consumer now imports from here, which is
! where they are defined. One abstract field, one concrete field
! (field_stored), before and after the move.
!
! THE SHAPE INVARIANT, now that the domain is mathematically real:
!
!      stored scalars  =  domain % num_members() * num_components()
!
! with the established interleaving
!
!      position = (domain_local_position - 1) * num_components + component
!
! Values are addressed by the DOMAIN'S local position, never by raw
! member value: a field on the subset declared { d a b } stores d's
! value first, whatever integer d is.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module field_calculus

  use util_precision  , only : dp
  use graph_fractal  , only : graph

  implicit none

  private
  public :: field
  public :: functional
  public :: FIELD_INTEGER, FIELD_REAL, FIELD_COMPLEX
  public :: FIELD_LOGICAL, FIELD_CHARACTER, FIELD_NONE

  !===================================================================!
  ! The five value kinds: one absorbed axis.
  !===================================================================!

  integer, parameter :: FIELD_INTEGER   = 1
  integer, parameter :: FIELD_REAL      = 2
  integer, parameter :: FIELD_COMPLEX   = 3
  ! A field that has never been set stores no kind. The absence is a
  ! member of the enumeration, not a default substituted for one.
  integer, parameter :: FIELD_NONE      = 0
  integer, parameter :: FIELD_LOGICAL   = 4
  integer, parameter :: FIELD_CHARACTER = 5

  !===================================================================!
  ! The abstract field: identity, domain, shape, and the
  ! plain-vector adapters - read once, compute in arrays, write back
  ! once. num_entries is the construction snapshot, not a query.
  !===================================================================!

  type, abstract :: field

     ! THE VALUES, stored once for every field in the tower: one kind
     ! at a time, of whichever kind was last set, in the order the
     ! domain lists its members with the components of one member
     ! adjacent. Only this module reads or writes them.
     class(*), allocatable, private :: values(:)

     ! THE DESCRIPTION, stored once for every field in the tower: the
     ! field's name, the unit of its values, WHICH set they are
     ! defined on and HOW MANY entries that set had at construction,
     ! and how many components each entry contains. A concretion
     ! states all of it through describe, once.
     character(len=:), allocatable, private :: label
     character(len=:), allocatable, private :: unit_name
     type(graph), private :: graph
     integer    , private :: ne = 0
     integer    , private :: nc = 1

   contains

     procedure :: name           => field_name
     procedure :: units          => field_units
     procedure :: domain         => field_domain
     procedure :: defined_on     => field_defined_on
     procedure :: inner_product  => field_inner_product
     procedure :: num_components => field_num_components
     procedure :: num_entries    => field_num_entries
     procedure :: value_kind     => field_value_kind
     procedure :: describe

     ! The ten adapters: read once, compute in arrays, write back once.
     ! A getter of the wrong kind returns a zero-length array; any
     ! setter replaces the values and the kind together, and must
     ! fill the domain exactly.
     procedure :: integer_vector       => field_integer_vector
     procedure :: set_integer_vector   => field_set_integer_vector
     procedure :: real_vector          => field_real_vector
     procedure :: real_values          => field_real_values
     procedure :: set_real_vector      => field_set_real_vector
     procedure :: complex_vector       => field_complex_vector
     procedure :: set_complex_vector   => field_set_complex_vector
     procedure :: logical_vector       => field_logical_vector
     procedure :: set_logical_vector   => field_set_logical_vector
     procedure :: character_vector     => field_character_vector
     procedure :: set_character_vector => field_set_character_vector
     procedure, private :: store

     !---------------------------------------------------------------!
     ! PLACE THIS VALUE AT A LOCATION OF THE SAME TYPE. Intrinsic
     ! assignment to a polymorphic location is barred unless that
     ! location is allocatable, and an element of an array never is.
     ! So a caller with an array of class(field) locations cannot
     ! fill one without naming a type. A concretion names the type
     ! once, here, and no caller needs to name it. A location of any
     ! other type is an error and stops the program with a message.
     !---------------------------------------------------------------!
     procedure(field_assign_in_interface), deferred :: assign_in

  end type field

  abstract interface

     subroutine field_assign_in_interface(this, location)
       import :: field
       class(field), intent(in)    :: this
       class(field), intent(inout) :: location
     end subroutine field_assign_in_interface

  end interface

  !===================================================================!
  ! FUNCTIONAL. The field at domain size one: a single
  ! value with the whole inherited interface. The type exists so
  ! an argument may require the one-entry case at compile time - a
  ! reduction returns a functional, not a field of size one. The
  ! invariant num_entries() == 1 is checked by the test suite,
  ! because the type system cannot state it.
  !===================================================================!

  type, abstract, extends(field) :: functional

   contains

     !----------------------------------------------------------------!
     ! The invariant of the type, stated here rather than stored: a
     ! functional is one entry of one component, however it was
     ! allocated.
     !----------------------------------------------------------------!

     procedure :: num_components => functional_one
     procedure :: num_entries    => functional_one

     !----------------------------------------------------------------!
     ! The scalar adapters: the vector adapters at length one, since
     ! a one-entry field and a scalar are the same value. Written
     ! here once for every functional, through the vector adapters
     ! only - no storage is defined at this level.
     !----------------------------------------------------------------!

     procedure :: integer_value       => functional_integer_value
     procedure :: set_integer_value   => functional_set_integer_value
     procedure :: real_value          => functional_real_value
     procedure :: set_real_value      => functional_set_real_value
     procedure :: complex_value       => functional_complex_value
     procedure :: set_complex_value   => functional_set_complex_value
     procedure :: logical_value       => functional_logical_value
     procedure :: set_logical_value   => functional_set_logical_value
     procedure :: character_value     => functional_character_value
     procedure :: set_character_value => functional_set_character_value

  end type functional

contains

  !===================================================================!
  ! The description, stated once by a concretion's constructor. An
  ! undeclared domain, or a negative entry count, stops the program.
  !===================================================================!

  subroutine describe(this, label, domain, num_entries, num_components, unit_name)

    class(field)    , intent(inout)        :: this
    character(len=*), intent(in)           :: label
    type(graph)     , intent(in)           :: domain
    integer         , intent(in)           :: num_entries
    integer         , intent(in), optional :: num_components
    character(len=*), intent(in), optional :: unit_name

    if (.not. domain % same_as(domain)) then
       error stop 'field: a field requires a declared domain'
    end if
    if (num_entries < 0) then
       error stop 'field: the entry count of a domain is not negative'
    end if

    this % label = label
    this % graph = domain
    this % ne    = num_entries
    this % nc    = 1
    if (present(num_components)) this % nc = num_components
    this % unit_name = '-'
    if (present(unit_name)) this % unit_name = unit_name

  end subroutine describe

  !===================================================================!
  ! The field's name, and the unit of its values: an empty name when
  ! no label was given, a dash when no unit was given.
  !===================================================================!

  pure function field_name(this) result(name)

    class(field), intent(in) :: this
    character(len=:), allocatable :: name

    name = ''
    if (allocated(this % label)) name = this % label

  end function field_name

  pure function field_units(this) result(units)

    class(field), intent(in) :: this
    character(len=:), allocatable :: units

    units = '-'
    if (allocated(this % unit_name)) units = this % unit_name

  end function field_units

  !===================================================================!
  ! WHICH set the values are defined on, by value: a copy of a set
  ! graph stores its token, so the result IS the domain - same_as
  ! decides, and no reference is returned. The counts were frozen at
  ! construction.
  !===================================================================!

  type(graph) function field_domain(this) result(domain)

    class(field), intent(in) :: this

    domain = this % graph

  end function field_domain

  pure integer function field_num_components(this)

    class(field), intent(in) :: this

    field_num_components = this % nc

  end function field_num_components

  pure integer function field_num_entries(this)

    class(field), intent(in) :: this

    field_num_entries = this % ne

  end function field_num_entries

  !===================================================================!
  ! Whether this field is defined on that domain: the same set by
  ! identity, never by extent. Every check that a state, a history
  ! state, a direction, a right-hand side or an action's result is
  ! defined on the domain a calculation expects evaluates this one
  ! predicate; what is rejected, and why, is stated at the call site.
  !===================================================================!

  logical function field_defined_on(this, domain) result(defined)

    class(field), intent(in) :: this
    type(graph) , intent(in) :: domain

    defined = this % graph % same_as(domain)

  end function field_defined_on

  !===================================================================!
  ! THE INNER PRODUCT sum_i this(i) other(i): the reduction
  ! dot_product performs on real arrays, taken here over two fields -
  ! two graphs with values. Fields on different domains, or of a kind
  ! other than real, stop the program rather than pair mismatched or
  ! non-numeric values as a silent zero.
  !===================================================================!

  real(dp) function field_inner_product(this, other) result(prod)

    class(field), intent(in) :: this
    class(field), intent(in) :: other

    real(dp), allocatable :: u(:), v(:)

    if (.not. this % defined_on(other % domain())) then
       error stop 'field: an inner product pairs fields on the same domain'
    end if
    if (this % value_kind() /= FIELD_REAL .or. other % value_kind() /= FIELD_REAL) then
       error stop 'field: an inner product pairs real-valued fields'
    end if

    call this  % real_vector(u)
    call other % real_vector(v)
    prod = dot_product(u, v)

  end function field_inner_product

  !===================================================================!
  ! The kind stored: read from the values themselves. A field that
  ! stores nothing yet reads as FIELD_NONE.
  !===================================================================!

  pure integer function field_value_kind(this) result(kind)

    class(field), intent(in) :: this

    kind = FIELD_NONE
    if (.not. allocated(this % values)) return

    kind = FIELD_REAL

    select type (stored => this % values)
    type is (integer)
       kind = FIELD_INTEGER
    type is (real(dp))
       kind = FIELD_REAL
    type is (complex(dp))
       kind = FIELD_COMPLEX
    type is (logical)
       kind = FIELD_LOGICAL
    type is (character(len=*))
       kind = FIELD_CHARACTER
    end select

  end function field_value_kind

  !===================================================================!
  ! The one setter behind the five: the values must fill the domain
  ! exactly - num_entries times num_components - or the program
  ! stops; they replace whatever kind was stored.
  !===================================================================!

  pure subroutine store(this, values)

    class(field), intent(inout) :: this
    class(*)    , intent(in)    :: values(:)

    if (size(values) /= this % num_entries() * this % num_components()) then
       error stop 'field: a value vector must fill its domain exactly'
    end if

    ! a store of the same kind and length is written in place: the
    ! deallocate-allocate pair performs two allocator operations, and a
    ! polymorphic copy moves the values one at a time through the
    ! type's own copy
    if (allocated(this % values)) then
       if (size(this % values) == size(values)) then
          select type (stored => this % values)
          type is (real(dp))
             select type (values)
             type is (real(dp))
                stored = values
                return
             end select
          type is (integer)
             select type (values)
             type is (integer)
                stored = values
                return
             end select
          type is (complex(dp))
             select type (values)
             type is (complex(dp))
                stored = values
                return
             end select
          type is (logical)
             select type (values)
             type is (logical)
                stored = values
                return
             end select
          end select
       end if
       deallocate(this % values)
    end if
    allocate(this % values, source=values)

  end subroutine store

  !===================================================================!
  ! The adapters, one pair per kind. A getter of another kind returns
  ! a zero-length array: no conversion, no inference, and a pure
  ! procedure has no error path, so the zero length is the indicator.
  !===================================================================!

  pure subroutine field_integer_vector(this, values)

    class(field), intent(in)          :: this
    integer, allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (integer)
          values = stored
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_integer_vector

  pure subroutine field_set_integer_vector(this, values)

    class(field), intent(inout) :: this
    integer     , intent(in)    :: values(:)

    call this % store(values)

  end subroutine field_set_integer_vector

  !===================================================================!
  ! The stored reals themselves, without a copy. A field of another
  ! kind returns null, which is the same indicator the zero-length
  ! getter gives. The result points into the field, so it is valid
  ! only while the field is unchanged.
  !===================================================================!

  function field_real_values(this) result(stored_values)

    class(field), intent(in), target :: this
    real(dp), pointer :: stored_values(:)

    stored_values => null()
    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (real(dp))
          stored_values => stored
       end select
    end if

  end function field_real_values

  pure subroutine field_real_vector(this, values)

    class(field), intent(in)           :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (real(dp))
          values = stored
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_real_vector

  pure subroutine field_set_real_vector(this, values)

    class(field), intent(inout) :: this
    real(dp)    , intent(in)    :: values(:)

    call this % store(values)

  end subroutine field_set_real_vector

  pure subroutine field_complex_vector(this, values)

    class(field), intent(in)              :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (complex(dp))
          values = stored
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_complex_vector

  pure subroutine field_set_complex_vector(this, values)

    class(field), intent(inout) :: this
    complex(dp) , intent(in)    :: values(:)

    call this % store(values)

  end subroutine field_set_complex_vector

  pure subroutine field_logical_vector(this, values)

    class(field), intent(in)          :: this
    logical, allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (logical)
          values = stored
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_logical_vector

  pure subroutine field_set_logical_vector(this, values)

    class(field), intent(inout) :: this
    logical     , intent(in)    :: values(:)

    call this % store(values)

  end subroutine field_set_logical_vector

  pure subroutine field_character_vector(this, values)

    class(field), intent(in)                   :: this
    character(len=:), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (stored => this % values)
       type is (character(len=*))
          values = stored
          return
       end select
    end if
    allocate(character(len=1) :: values(0))

  end subroutine field_character_vector

  pure subroutine field_set_character_vector(this, values)

    class(field)    , intent(inout) :: this
    character(len=*), intent(in)    :: values(:)

    call this % store(values)

  end subroutine field_set_character_vector

  pure integer function functional_one(this)

    class(functional), intent(in) :: this

    associate (u1 => this); end associate
    functional_one = 1

  end function functional_one

  !===================================================================!
  ! The scalar adapters of a functional. A getter reads the vector
  ! adapter of its kind and takes the one entry; when the functional
  ! stores another kind the adapter returns zero-length and the getter
  ! returns the zero of the requested kind. A setter passes the one
  ! value to the vector adapter, which replaces the value and the kind.
  !===================================================================!

  pure subroutine functional_integer_value(this, value)

    class(functional), intent(in) :: this
    integer, intent(out) :: value

    integer, allocatable :: t(:)

    call this % integer_vector(t)
    value = 0
    if (size(t) >= 1) value = t(1)

  end subroutine functional_integer_value

  pure subroutine functional_set_integer_value(this, value)

    class(functional), intent(inout) :: this
    integer, intent(in) :: value

    call this % set_integer_vector([value])

  end subroutine functional_set_integer_value

  pure subroutine functional_real_value(this, value)

    class(functional), intent(in) :: this
    real(dp), intent(out) :: value

    real(dp), allocatable :: t(:)

    call this % real_vector(t)
    value = 0.0_dp
    if (size(t) >= 1) value = t(1)

  end subroutine functional_real_value

  pure subroutine functional_set_real_value(this, value)

    class(functional), intent(inout) :: this
    real(dp), intent(in) :: value

    call this % set_real_vector([value])

  end subroutine functional_set_real_value

  pure subroutine functional_complex_value(this, value)

    class(functional), intent(in) :: this
    complex(dp), intent(out) :: value

    complex(dp), allocatable :: t(:)

    call this % complex_vector(t)
    value = (0.0_dp, 0.0_dp)
    if (size(t) >= 1) value = t(1)

  end subroutine functional_complex_value

  pure subroutine functional_set_complex_value(this, value)

    class(functional), intent(inout) :: this
    complex(dp), intent(in) :: value

    call this % set_complex_vector([value])

  end subroutine functional_set_complex_value

  pure subroutine functional_logical_value(this, value)

    class(functional), intent(in) :: this
    logical, intent(out) :: value

    logical, allocatable :: t(:)

    call this % logical_vector(t)
    value = .false.
    if (size(t) >= 1) value = t(1)

  end subroutine functional_logical_value

  pure subroutine functional_set_logical_value(this, value)

    class(functional), intent(inout) :: this
    logical, intent(in) :: value

    call this % set_logical_vector([value])

  end subroutine functional_set_logical_value

  pure subroutine functional_character_value(this, value)

    class(functional), intent(in) :: this
    character(len=:), allocatable, intent(out) :: value

    character(len=:), allocatable :: t(:)

    call this % character_vector(t)
    value = ''
    if (size(t) >= 1) value = t(1)

  end subroutine functional_character_value

  pure subroutine functional_set_character_value(this, value)

    class(functional), intent(inout) :: this
    character(len=*), intent(in) :: value

    call this % set_character_vector([value])

  end subroutine functional_set_character_value

end module field_calculus
