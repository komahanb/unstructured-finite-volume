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
! graph_field and its value-kind constants, moved here from the old
! grammar, which now re-exports them for its remaining consumers
! but no longer owns them. One abstract field, one concrete field
! (class_graph_field), before and after the move.
!
! THE SHAPE INVARIANT, now that the domain is mathematically real:
!
!      stored scalars  =  domain % size() * num_components()
!
! with the established interleaving
!
!      position = (domain_local_position - 1) * ncomp + component
!
! Values are addressed by the DOMAIN'S local position, never by raw
! member value: a field on the subset declared { d a b } stores d's
! value first, whatever integer d happens to be.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_field_calculus

  use iso_fortran_env, only : dp => REAL64
  use fractal_graph  , only : set_graph => graph

  implicit none

  private
  public :: graph_field
  public :: set_graph
  public :: GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL, GRAPH_FIELD_COMPLEX
  public :: GRAPH_FIELD_LOGICAL, GRAPH_FIELD_CHARACTER

  !===================================================================!
  ! The five value kinds: one absorbed axis, as ever.
  !===================================================================!

  integer, parameter :: GRAPH_FIELD_INTEGER   = 1
  integer, parameter :: GRAPH_FIELD_REAL      = 2
  integer, parameter :: GRAPH_FIELD_COMPLEX   = 3
  integer, parameter :: GRAPH_FIELD_LOGICAL   = 4
  integer, parameter :: GRAPH_FIELD_CHARACTER = 5

  !===================================================================!
  ! The abstract field: identity, domain, shape, and the
  ! plain-vector adapters - fetch once, work in arrays, write back
  ! once. num_entries is the construction snapshot, not a query.
  !===================================================================!

  type, abstract :: graph_field

   contains

     procedure(field_name_interface), deferred :: name
     procedure(field_name_interface), deferred :: units

     procedure(field_domain_interface), deferred :: domain
     procedure(field_count_interface) , deferred :: num_components
     procedure(field_count_interface) , deferred :: num_entries
     procedure(field_count_interface) , deferred :: value_kind

     procedure(field_get_integer_interface)  , deferred :: get_integer_vector
     procedure(field_set_integer_interface)  , deferred :: set_integer_vector
     procedure(field_get_real_interface)     , deferred :: get_real_vector
     procedure(field_set_real_interface)     , deferred :: set_real_vector
     procedure(field_get_complex_interface)  , deferred :: get_complex_vector
     procedure(field_set_complex_interface)  , deferred :: set_complex_vector
     procedure(field_get_logical_interface)  , deferred :: get_logical_vector
     procedure(field_set_logical_interface)  , deferred :: set_logical_vector
     procedure(field_get_character_interface), deferred :: get_character_vector
     procedure(field_set_character_interface), deferred :: set_character_vector

  end type graph_field

  abstract interface

     pure function field_name_interface(this) result(name)
       import :: graph_field
       class(graph_field), intent(in) :: this
       character(len=:), allocatable :: name
     end function field_name_interface

     !--------------------------------------------------------------!
     ! WHICH set the values live on, by value. A copy of a set graph
     ! carries its token, so the answer IS the domain - same_as
     ! decides, and nothing is lent.
     !--------------------------------------------------------------!

     type(set_graph) function field_domain_interface(this)
       import :: graph_field, set_graph
       class(graph_field), intent(in) :: this
     end function field_domain_interface

     pure integer function field_count_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_count_interface

     pure subroutine field_get_integer_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       integer, allocatable, intent(out) :: values(:)
     end subroutine field_get_integer_interface

     pure subroutine field_set_integer_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       integer, intent(in) :: values(:)
     end subroutine field_set_integer_interface

     pure subroutine field_get_real_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(in) :: this
       real(dp), allocatable, intent(out) :: values(:)
     end subroutine field_get_real_interface

     pure subroutine field_set_real_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(inout) :: this
       real(dp), intent(in) :: values(:)
     end subroutine field_set_real_interface

     pure subroutine field_get_complex_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(in) :: this
       complex(dp), allocatable, intent(out) :: values(:)
     end subroutine field_get_complex_interface

     pure subroutine field_set_complex_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(inout) :: this
       complex(dp), intent(in) :: values(:)
     end subroutine field_set_complex_interface

     pure subroutine field_get_logical_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       logical, allocatable, intent(out) :: values(:)
     end subroutine field_get_logical_interface

     pure subroutine field_set_logical_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       logical, intent(in) :: values(:)
     end subroutine field_set_logical_interface

     pure subroutine field_get_character_interface(this, values)
       import :: graph_field
       class(graph_field), intent(in) :: this
       character(len=:), allocatable, intent(out) :: values(:)
     end subroutine field_get_character_interface

     pure subroutine field_set_character_interface(this, values)
       import :: graph_field
       class(graph_field), intent(inout) :: this
       character(len=*), intent(in) :: values(:)
     end subroutine field_set_character_interface

  end interface

end module graph_field_calculus
