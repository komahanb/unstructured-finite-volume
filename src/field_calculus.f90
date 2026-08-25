!=====================================================================!
! LEVEL 5 OF THE NEW TOWER . THE FIELD CALCULUS
!
! The level answers one question: WHAT VALUES LIVE ON A DOMAIN. A
! field is a function over one finite domain,
!
!      f : A -> V         or        f : S -> V,   S c--> A
!
! and that domain is named by a SET GRAPH: one identity, never a
! union, never a side flag (AGENTS.md 20, CALCULATOR.md 12). A field
! needs a domain; it does not need a graph container.
!
!             WHAT A FIELD KEEPS, AND WHAT IT ASKS FOR
!
! Identity and a frozen count, and nothing else about the domain:
!
!      domain()        WHICH set - a set graph, by value
!      num_entries()   HOW MANY - the count taken at construction
!
! Those are the only two questions every field caller asks, and
! neither is a question about membership, so neither needs a map.
! Whoever wants to know WHO belongs, WHERE a member stands, or WHAT
! the domain is called asks the set map, the inclusion map or the
! label map - explicitly, holding them, at the call site.
!
! The count is frozen because it always was: a field has held a COPY
! of its domain since the first version, so num_entries never tracked
! later mutation of the caller's set. Storing the integer says out
! loud what copying said in private, and stops N fields on one domain
! from holding N copies of its extension.
!
! This module is the REHOMED field ontology: the one abstract
! field and its value-kind constants. They came from the old
! grammar, which re-exported them for a while and was drained and
! deleted in PR2; every consumer now asks here, which is where they
! are defined. One abstract field, one concrete field
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
! value first, whatever integer d happens to be.
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
  public :: FIELD_LOGICAL, FIELD_CHARACTER

  !===================================================================!
  ! The five value kinds: one absorbed axis, as ever.
  !===================================================================!

  integer, parameter :: FIELD_INTEGER   = 1
  integer, parameter :: FIELD_REAL      = 2
  integer, parameter :: FIELD_COMPLEX   = 3
  integer, parameter :: FIELD_LOGICAL   = 4
  integer, parameter :: FIELD_CHARACTER = 5

  !===================================================================!
  ! The abstract field: identity, domain, shape, and the
  ! plain-vector adapters - fetch once, work in arrays, write back
  ! once. num_entries is the construction snapshot, not a query.
  !===================================================================!

  type, abstract :: field

     ! THE VALUES, held once for every field in the tower: one kind
     ! at a time, of whichever kind was last set, in the order the
     ! domain lists its members with the components of one member
     ! next to each other. Only this module reads or writes them.
     class(*), allocatable, private :: values(:)

   contains

     procedure(field_name_interface), deferred :: name
     procedure(field_name_interface), deferred :: units

     procedure(field_domain_interface), deferred :: domain
     procedure :: defined_on => field_defined_on
     procedure(field_count_interface) , deferred :: num_components
     procedure(field_count_interface) , deferred :: num_entries
     procedure :: value_kind => field_value_kind

     ! The ten adapters: fetch once, work in arrays, write back once.
     ! A getter of the wrong kind answers a zero-length array; any
     ! setter replaces the values and the kind together, and must
     ! fill the domain exactly.
     procedure :: integer_vector       => field_integer_vector
     procedure :: set_integer_vector   => field_set_integer_vector
     procedure :: real_vector          => field_real_vector
     procedure :: set_real_vector      => field_set_real_vector
     procedure :: complex_vector       => field_complex_vector
     procedure :: set_complex_vector   => field_set_complex_vector
     procedure :: logical_vector       => field_logical_vector
     procedure :: set_logical_vector   => field_set_logical_vector
     procedure :: character_vector     => field_character_vector
     procedure :: set_character_vector => field_set_character_vector
     procedure, private :: hold

  end type field

  !===================================================================!
  ! FUNCTIONAL. The field at domain size one: a single
  ! value with the whole inherited interface. The type exists so
  ! an argument may demand the one-entry case at compile time - a
  ! reduction returns a functional, not a field that happens to be
  ! small. The law num_entries() == 1 is held by the test suite,
  ! because the type system cannot state it.
  !===================================================================!

  type, abstract, extends(field) :: functional

   contains

     !----------------------------------------------------------------!
     ! The scalar adapters: the vector adapters at length one, since
     ! a one-entry field and a scalar are the same value. Written
     ! here once for every functional, through the vector adapters
     ! only - no storage is known at this level.
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

  abstract interface

     pure function field_name_interface(this) result(name)
       import :: field
       class(field), intent(in) :: this
       character(len=:), allocatable :: name
     end function field_name_interface

     !--------------------------------------------------------------!
     ! WHICH set the values live on, by value. A copy of a set graph
     ! carries its token, so the answer IS the domain - same_as
     ! decides, and nothing is lent.
     !--------------------------------------------------------------!

     type(graph) function field_domain_interface(this)
       import :: field, graph
       class(field), intent(in) :: this
     end function field_domain_interface

     pure integer function field_count_interface(this)
       import :: field
       class(field), intent(in) :: this
     end function field_count_interface

  end interface

contains

  !===================================================================!
  ! Whether this field is defined on that domain: the same set by
  ! identity, never by extent. Every check that a state, a history
  ! state, a direction, a right-hand side or an action's result lives
  ! where a calculation expects it asks this one question; what is
  ! refused, and why, is said at the call site.
  !===================================================================!

  logical function field_defined_on(this, domain) result(defined)

    class(field), intent(in) :: this
    type(graph) , intent(in) :: domain

    type(graph) :: on

    on      = this % domain()
    defined = on % same_as(domain)

  end function field_defined_on

  !===================================================================!
  ! The kind held: read off the values themselves. A field that holds
  ! nothing yet reads as real, the kind a field is born to.
  !===================================================================!

  pure integer function field_value_kind(this) result(kind)

    class(field), intent(in) :: this

    kind = FIELD_REAL
    if (.not. allocated(this % values)) return

    select type (held => this % values)
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
  ! stops; they replace whatever kind was held.
  !===================================================================!

  pure subroutine hold(this, values)

    class(field), intent(inout) :: this
    class(*)    , intent(in)    :: values(:)

    if (size(values) /= this % num_entries() * this % num_components()) then
       error stop 'field: a value vector must fill its domain exactly'
    end if

    if (allocated(this % values)) deallocate(this % values)
    allocate(this % values, source=values)

  end subroutine hold

  !===================================================================!
  ! The adapters, one pair per kind. A getter of another kind answers
  ! a zero-length array: no conversion, no inference, and a pure
  ! procedure has no error path, so the zero length is the signal.
  !===================================================================!

  pure subroutine field_integer_vector(this, values)

    class(field), intent(in)          :: this
    integer, allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (held => this % values)
       type is (integer)
          values = held
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_integer_vector

  pure subroutine field_set_integer_vector(this, values)

    class(field), intent(inout) :: this
    integer     , intent(in)    :: values(:)

    call this % hold(values)

  end subroutine field_set_integer_vector

  pure subroutine field_real_vector(this, values)

    class(field), intent(in)           :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (held => this % values)
       type is (real(dp))
          values = held
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_real_vector

  pure subroutine field_set_real_vector(this, values)

    class(field), intent(inout) :: this
    real(dp)    , intent(in)    :: values(:)

    call this % hold(values)

  end subroutine field_set_real_vector

  pure subroutine field_complex_vector(this, values)

    class(field), intent(in)              :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (held => this % values)
       type is (complex(dp))
          values = held
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_complex_vector

  pure subroutine field_set_complex_vector(this, values)

    class(field), intent(inout) :: this
    complex(dp) , intent(in)    :: values(:)

    call this % hold(values)

  end subroutine field_set_complex_vector

  pure subroutine field_logical_vector(this, values)

    class(field), intent(in)          :: this
    logical, allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (held => this % values)
       type is (logical)
          values = held
          return
       end select
    end if
    allocate(values(0))

  end subroutine field_logical_vector

  pure subroutine field_set_logical_vector(this, values)

    class(field), intent(inout) :: this
    logical     , intent(in)    :: values(:)

    call this % hold(values)

  end subroutine field_set_logical_vector

  pure subroutine field_character_vector(this, values)

    class(field), intent(in)                   :: this
    character(len=:), allocatable, intent(out) :: values(:)

    if (allocated(this % values)) then
       select type (held => this % values)
       type is (character(len=*))
          values = held
          return
       end select
    end if
    allocate(character(len=1) :: values(0))

  end subroutine field_character_vector

  pure subroutine field_set_character_vector(this, values)

    class(field)    , intent(inout) :: this
    character(len=*), intent(in)    :: values(:)

    call this % hold(values)

  end subroutine field_set_character_vector

  !===================================================================!
  ! The scalar adapters of a functional. A getter reads the vector
  ! adapter of its kind and takes the one entry; when the functional
  ! holds another kind the adapter answers zero-length and the getter
  ! answers the zero of the asked kind. A setter hands the one value
  ! to the vector adapter, which replaces the value and the kind.
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
