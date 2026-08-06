!=====================================================================!
! Values carried on a graph.
!
! A field holds values over the members of a domain. A functional is
! the case where the domain has one member: a single answer, with the
! whole battery of a field still available to read it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module data_graph_field

  use iso_fortran_env, only : dp => REAL64
  use structure_graph, only : graph

  implicit none

  private

  public :: graph_field
  public :: graph_functional
  public :: GRAPH_FIELD_INTEGER
  public :: GRAPH_FIELD_REAL
  public :: GRAPH_FIELD_COMPLEX
  public :: GRAPH_FIELD_LOGICAL
  public :: GRAPH_FIELD_CHARACTER

  integer, parameter :: GRAPH_FIELD_INTEGER   = 1
  integer, parameter :: GRAPH_FIELD_REAL      = 2
  integer, parameter :: GRAPH_FIELD_COMPLEX   = 3
  integer, parameter :: GRAPH_FIELD_LOGICAL   = 4
  integer, parameter :: GRAPH_FIELD_CHARACTER = 5

  !===================================================================!
  ! GRAPH_FIELD. The carrier of values.
  !
  ! A field is a function: values over a domain, and the domain is a
  ! graph - a member set of some host. One value kind per field,
  ! num_components values per entry, laid out by the formula in the
  ! header.
  !
  !      field  ---get--->  [ v1 v2 v3 v4 ]  ---> a solver, a file
  !                                               writer, an outside
  !      field  <--set----  [ v1 v2 v3 v4 ]  <--- library; the answer
  !                                               coming back
  !
  ! Fetch once, work in plain arrays where the arithmetic is free,
  ! write back once. Scaling and adding are not graph theory and have
  ! no procedures here on purpose. One get/set pair per value kind:
  ! the five roads of one absorbed axis.
  !
  ! A field whose domain has ONE entry is a single number. Level 1
  ! names that case the functional; nothing in this role depends on
  ! the size of the domain, which is why the name can wait.
  !
  ! num_entries repeats the size of the domain. That repetition is a
  ! priced convenience for call sites, recorded here so nobody
  ! mistakes it for a generator.
  !===================================================================!

  type, abstract :: graph_field

   contains

     ! Identity: what it is called and what unit it carries.
     procedure(field_name_interface), deferred :: name
     procedure(field_name_interface), deferred :: units

     ! Where it lives and how it is shaped.
     procedure(field_domain_interface), deferred :: domain
     procedure(field_count_interface) , deferred :: num_components
     procedure(field_count_interface) , deferred :: num_entries
     procedure(field_count_interface) , deferred :: value_kind

     ! The plain-vector adapters, one pair per value kind.
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

  !===================================================================!
  ! GRAPH_FUNCTIONAL. The grammar's field at domain size one: a
  ! single value, with the whole inherited battery intact.
  !
  !      a field    [ v1 v2 v3 v4 v5 v6 ]     values over members
  !      a functional      [ J ]              one value: the answer
  !
  ! Nothing is declared here because nothing new can be asked: name,
  ! units, kind, and the get/set roads all mean exactly what they
  ! meant one level down, read at one entry. The type exists so that
  ! an argument may DEMAND the one-entry case at compile time - a
  ! reduction returns a functional, not a field that happens to be
  ! small. The suite holds the law the compiler cannot:
  !
  !      num_entries() == 1
  !===================================================================!

  type, abstract, extends(graph_field) :: graph_functional
  end type graph_functional

  abstract interface

     pure function field_name_interface(this) result(name)
       import :: graph_field
       class(graph_field), intent(in) :: this
       character(len=:), allocatable :: name
     end function field_name_interface
     subroutine field_domain_interface(this, domain)
       import :: graph_field, graph
       class(graph_field), intent(in) :: this
       class(graph), allocatable, intent(out) :: domain
     end subroutine field_domain_interface
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

end module data_graph_field
