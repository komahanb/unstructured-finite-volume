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

  use iso_fortran_env, only : dp => REAL64
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

   contains

     procedure(field_name_interface), deferred :: name
     procedure(field_name_interface), deferred :: units

     procedure(field_domain_interface), deferred :: domain
     procedure :: defined_on => field_defined_on
     procedure(field_count_interface) , deferred :: num_components
     procedure(field_count_interface) , deferred :: num_entries
     procedure(field_count_interface) , deferred :: value_kind

     procedure(field_integer_vector_interface)  , deferred :: integer_vector
     procedure(field_set_integer_vector_interface)  , deferred :: set_integer_vector
     procedure(field_real_vector_interface)     , deferred :: real_vector
     procedure(field_set_real_vector_interface)     , deferred :: set_real_vector
     procedure(field_complex_vector_interface)  , deferred :: complex_vector
     procedure(field_set_complex_vector_interface)  , deferred :: set_complex_vector
     procedure(field_logical_vector_interface)  , deferred :: logical_vector
     procedure(field_set_logical_vector_interface)  , deferred :: set_logical_vector
     procedure(field_character_vector_interface), deferred :: character_vector
     procedure(field_set_character_vector_interface), deferred :: set_character_vector

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

     pure subroutine field_integer_vector_interface(this, values)
       import :: field
       class(field), intent(in) :: this
       integer, allocatable, intent(out) :: values(:)
     end subroutine field_integer_vector_interface

     pure subroutine field_set_integer_vector_interface(this, values)
       import :: field
       class(field), intent(inout) :: this
       integer, intent(in) :: values(:)
     end subroutine field_set_integer_vector_interface

     pure subroutine field_real_vector_interface(this, values)
       import :: field, dp
       class(field), intent(in) :: this
       real(dp), allocatable, intent(out) :: values(:)
     end subroutine field_real_vector_interface

     pure subroutine field_set_real_vector_interface(this, values)
       import :: field, dp
       class(field), intent(inout) :: this
       real(dp), intent(in) :: values(:)
     end subroutine field_set_real_vector_interface

     pure subroutine field_complex_vector_interface(this, values)
       import :: field, dp
       class(field), intent(in) :: this
       complex(dp), allocatable, intent(out) :: values(:)
     end subroutine field_complex_vector_interface

     pure subroutine field_set_complex_vector_interface(this, values)
       import :: field, dp
       class(field), intent(inout) :: this
       complex(dp), intent(in) :: values(:)
     end subroutine field_set_complex_vector_interface

     pure subroutine field_logical_vector_interface(this, values)
       import :: field
       class(field), intent(in) :: this
       logical, allocatable, intent(out) :: values(:)
     end subroutine field_logical_vector_interface

     pure subroutine field_set_logical_vector_interface(this, values)
       import :: field
       class(field), intent(inout) :: this
       logical, intent(in) :: values(:)
     end subroutine field_set_logical_vector_interface

     pure subroutine field_character_vector_interface(this, values)
       import :: field
       class(field), intent(in) :: this
       character(len=:), allocatable, intent(out) :: values(:)
     end subroutine field_character_vector_interface

     pure subroutine field_set_character_vector_interface(this, values)
       import :: field
       class(field), intent(inout) :: this
       character(len=*), intent(in) :: values(:)
     end subroutine field_set_character_vector_interface

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
