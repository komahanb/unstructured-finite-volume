!=====================================================================!
! Concrete graph fields.
!
! A field is one value for every entry of a support. There is exactly
! one concrete vertex field here and exactly one concrete edge field,
! and each holds every kind of value the contract allows. Because
! there is only one of each, a graph holds its data in a plain array:
!
!      type(vertex_field), allocatable :: vdata(:)
!
! A Fortran array holds a single dynamic type, so a plain array of
! fields requires a single concrete field type.
!
!=====================================================================!
!
!                        THE VALUE-KIND RULE
!
! A field holds one kind of value at a time. Only the matching store
! is ever allocated:
!
!      value_kind() = GRAPH_FIELD_REAL     ->   only rvals is live
!      value_kind() = GRAPH_FIELD_INTEGER  ->   only ivals is live
!
! From that, three rules that hold for all ten accessors:
!
!      ask first     a caller checks value_kind() before reaching for
!                    a vector
!
!      wrong getter  returns a zero-length array. No conversion and
!                    no inference happens, and a pure procedure has
!                    no error path. The zero-length result is the
!                    signal
!
!      any setter    replaces both the values and the kind. Setting
!                    reals onto a field that held integers makes it a
!                    real field
!
! No conversion happens anywhere. A field that holds boundary names
! will not hand them back as numbers.
!
!=====================================================================!
!
!                        WHERE A VALUE SITS
!
! A field stores its values in the order the support lists its cells,
! and keeps the components of one cell next to each other:
!
!      support indices     7        7        3        3
!      component       1        2        1        2
!                   +--------+--------+--------+--------+
!      values       |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
!                   +--------+--------+--------+--------+
!
!      position = (entry_position - 1) * num_components + component
!
! Everything that reads a flat vector out of a field depends on this -
! a linear solver, a file writer, a matrix adapter. It is the reason
! this library needs no degree-of-freedom bookkeeping of its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_field

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_support
  use abstract_graph_types, only : graph_vertex_support, graph_edge_support
  use abstract_graph_types, only : graph_vertex_field, graph_edge_field
  use abstract_graph_types, only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use abstract_graph_types, only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use abstract_graph_types, only : GRAPH_FIELD_CHARACTER
  use class_graph_support , only : vertex_support, edge_support

  implicit none

  private
  public :: vertex_field, edge_field

  !===================================================================!
  ! Values on a set of vertices.
  !===================================================================!

  type, extends(graph_vertex_field) :: vertex_field

     character(len=:), allocatable :: label
     character(len=:), allocatable :: unit_name

     type(vertex_support) :: on

     integer :: ncomp = 1
     integer :: vkind = GRAPH_FIELD_REAL

     integer         , allocatable :: ivals(:)
     real(dp)        , allocatable :: rvals(:)
     complex(dp)     , allocatable :: cvals(:)
     logical         , allocatable :: lvals(:)
     character(len=:), allocatable :: svals(:)

   contains

     !----------------------------------------------------------------!
     ! Identity, and where the values sit.
     !----------------------------------------------------------------!

     procedure :: name           => vertex_field_name
     procedure :: units          => vertex_field_units
     procedure :: support        => vertex_field_support
     procedure :: vertex_support => vertex_field_vertex_support

     !----------------------------------------------------------------!
     ! Shape and kind.
     !----------------------------------------------------------------!

     procedure :: num_components => vertex_field_num_components
     procedure :: num_entries    => vertex_field_num_entries
     procedure :: value_kind     => vertex_field_value_kind

     !----------------------------------------------------------------!
     ! The plain-vector adapters, one pair per kind.
     !----------------------------------------------------------------!

     procedure :: get_integer_vector   => vertex_field_get_integer_vector
     procedure :: set_integer_vector   => vertex_field_set_integer_vector
     procedure :: get_real_vector      => vertex_field_get_real_vector
     procedure :: set_real_vector      => vertex_field_set_real_vector
     procedure :: get_complex_vector   => vertex_field_get_complex_vector
     procedure :: set_complex_vector   => vertex_field_set_complex_vector
     procedure :: get_logical_vector   => vertex_field_get_logical_vector
     procedure :: set_logical_vector   => vertex_field_set_logical_vector
     procedure :: get_character_vector => vertex_field_get_character_vector
     procedure :: set_character_vector => vertex_field_set_character_vector

  end type vertex_field

  !===================================================================!
  ! Values on a set of edges. The same story, told about faces.
  !===================================================================!

  type, extends(graph_edge_field) :: edge_field

     character(len=:), allocatable :: label
     character(len=:), allocatable :: unit_name

     type(edge_support) :: on

     integer :: ncomp = 1
     integer :: vkind = GRAPH_FIELD_REAL

     integer         , allocatable :: ivals(:)
     real(dp)        , allocatable :: rvals(:)
     complex(dp)     , allocatable :: cvals(:)
     logical         , allocatable :: lvals(:)
     character(len=:), allocatable :: svals(:)

   contains

     !----------------------------------------------------------------!
     ! Identity, and where the values sit.
     !----------------------------------------------------------------!

     procedure :: name         => edge_field_name
     procedure :: units        => edge_field_units
     procedure :: support      => edge_field_support
     procedure :: edge_support => edge_field_edge_support

     !----------------------------------------------------------------!
     ! Shape and kind.
     !----------------------------------------------------------------!

     procedure :: num_components => edge_field_num_components
     procedure :: num_entries    => edge_field_num_entries
     procedure :: value_kind     => edge_field_value_kind

     !----------------------------------------------------------------!
     ! The plain-vector adapters, one pair per kind.
     !----------------------------------------------------------------!

     procedure :: get_integer_vector   => edge_field_get_integer_vector
     procedure :: set_integer_vector   => edge_field_set_integer_vector
     procedure :: get_real_vector      => edge_field_get_real_vector
     procedure :: set_real_vector      => edge_field_set_real_vector
     procedure :: get_complex_vector   => edge_field_get_complex_vector
     procedure :: set_complex_vector   => edge_field_set_complex_vector
     procedure :: get_logical_vector   => edge_field_get_logical_vector
     procedure :: set_logical_vector   => edge_field_set_logical_vector
     procedure :: get_character_vector => edge_field_get_character_vector
     procedure :: set_character_vector => edge_field_set_character_vector

  end type edge_field

  !===================================================================!
  ! Constructors. Name the field, say where it lives, say how wide
  ! each entry is. The values arrive afterwards through a setter,
  ! which is also what fixes the kind.
  !===================================================================!

  interface vertex_field
     module procedure create_vertex_field
  end interface vertex_field

  interface edge_field
     module procedure create_edge_field
  end interface edge_field

contains

  !===================================================================!
  ! Build an empty vertex field on a set of vertices.
  !===================================================================!

  pure type(vertex_field) function create_vertex_field(label, on, ncomp, unit_name) result(this)

    character(len=*)    , intent(in)           :: label
    type(vertex_support), intent(in)           :: on
    integer             , intent(in), optional :: ncomp
    character(len=*)    , intent(in), optional :: unit_name

    this % label = label
    this % on    = on

    if (present(ncomp)) then
       this % ncomp = ncomp
    else
       this % ncomp = 1
    end if

    if (present(unit_name)) then
       this % unit_name = unit_name
    else
       this % unit_name = '-'
    end if

  end function create_vertex_field

  !===================================================================!
  ! Build an empty edge field on a set of edges.
  !===================================================================!

  pure type(edge_field) function create_edge_field(label, on, ncomp, unit_name) result(this)

    character(len=*)  , intent(in)           :: label
    type(edge_support), intent(in)           :: on
    integer           , intent(in), optional :: ncomp
    character(len=*)  , intent(in), optional :: unit_name

    this % label = label
    this % on    = on

    if (present(ncomp)) then
       this % ncomp = ncomp
    else
       this % ncomp = 1
    end if

    if (present(unit_name)) then
       this % unit_name = unit_name
    else
       this % unit_name = '-'
    end if

  end function create_edge_field

  !===================================================================!
  ! Vertex field: identity and support.
  !===================================================================!

  pure function vertex_field_name(this) result(name)

    class(vertex_field), intent(in) :: this
    character(len=:), allocatable   :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function vertex_field_name

  !===================================================================!
  ! What the values are measured in; a dash when nobody said.
  !===================================================================!

  pure function vertex_field_units(this) result(units)

    class(vertex_field), intent(in) :: this
    character(len=:), allocatable   :: units

    if (allocated(this % unit_name)) then
       units = this % unit_name
    else
       units = '-'
    end if

  end function vertex_field_units

  !===================================================================!
  ! The set the values live on, answered as a plain support.
  !===================================================================!

  subroutine vertex_field_support(this, support)

    class(vertex_field), intent(in)                :: this
    class(graph_support), allocatable, intent(out) :: support

    allocate(support, source=this % on)

  end subroutine vertex_field_support

  !===================================================================!
  ! The same set, answered with its vertex type on show.
  !===================================================================!

  subroutine vertex_field_vertex_support(this, support)

    class(vertex_field), intent(in)                       :: this
    class(graph_vertex_support), allocatable, intent(out) :: support

    allocate(support, source=this % on)

  end subroutine vertex_field_vertex_support

  !===================================================================!
  ! Vertex field: shape and kind. The entry count comes from the
  ! support, so a field cannot disagree with the set it lives on.
  !===================================================================!

  pure integer function vertex_field_num_components(this)

    class(vertex_field), intent(in) :: this

    vertex_field_num_components = this % ncomp

  end function vertex_field_num_components

  !===================================================================!
  ! One entry per vertex of the support.
  !===================================================================!

  pure integer function vertex_field_num_entries(this)

    class(vertex_field), intent(in) :: this

    vertex_field_num_entries = this % on % size()

  end function vertex_field_num_entries

  !===================================================================!
  ! Which of the five kinds is live.
  !===================================================================!

  pure integer function vertex_field_value_kind(this)

    class(vertex_field), intent(in) :: this

    vertex_field_value_kind = this % vkind

  end function vertex_field_value_kind

  !===================================================================!
  ! Vertex field: the adapters. A getter for the wrong kind answers
  ! with a zero-length array; a setter replaces the values and the
  ! kind together.
  !===================================================================!

  pure subroutine vertex_field_get_integer_vector(this, values)

    class(vertex_field), intent(in)   :: this
    integer, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_INTEGER .and. allocated(this % ivals)) then
       values = this % ivals
    else
       allocate(values(0))
    end if

  end subroutine vertex_field_get_integer_vector

  !===================================================================!
  ! Take integers in; the live kind becomes integer.
  !===================================================================!

  pure subroutine vertex_field_set_integer_vector(this, values)

    class(vertex_field), intent(inout) :: this
    integer             , intent(in)   :: values(:)

    this % ivals = values
    this % vkind = GRAPH_FIELD_INTEGER

  end subroutine vertex_field_set_integer_vector

  !===================================================================!
  ! The values, when the live kind is real; zero-length otherwise.
  !===================================================================!

  pure subroutine vertex_field_get_real_vector(this, values)

    class(vertex_field), intent(in)    :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_REAL .and. allocated(this % rvals)) then
       values = this % rvals
    else
       allocate(values(0))
    end if

  end subroutine vertex_field_get_real_vector

  !===================================================================!
  ! Take reals in; the live kind becomes real.
  !===================================================================!

  pure subroutine vertex_field_set_real_vector(this, values)

    class(vertex_field), intent(inout) :: this
    real(dp)           , intent(in)    :: values(:)

    this % rvals = values
    this % vkind = GRAPH_FIELD_REAL

  end subroutine vertex_field_set_real_vector

  !===================================================================!
  ! The values, when the live kind is complex; zero-length otherwise.
  !===================================================================!

  pure subroutine vertex_field_get_complex_vector(this, values)

    class(vertex_field), intent(in)       :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_COMPLEX .and. allocated(this % cvals)) then
       values = this % cvals
    else
       allocate(values(0))
    end if

  end subroutine vertex_field_get_complex_vector

  !===================================================================!
  ! Take complex values in; the live kind becomes complex. This is
  ! the seat a complex-step perturbation rides into the engine.
  !===================================================================!

  pure subroutine vertex_field_set_complex_vector(this, values)

    class(vertex_field), intent(inout) :: this
    complex(dp)        , intent(in)    :: values(:)

    this % cvals = values
    this % vkind = GRAPH_FIELD_COMPLEX

  end subroutine vertex_field_set_complex_vector

  !===================================================================!
  ! The values, when the live kind is logical; zero-length otherwise.
  !===================================================================!

  pure subroutine vertex_field_get_logical_vector(this, values)

    class(vertex_field), intent(in)   :: this
    logical, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_LOGICAL .and. allocated(this % lvals)) then
       values = this % lvals
    else
       allocate(values(0))
    end if

  end subroutine vertex_field_get_logical_vector

  !===================================================================!
  ! Take logicals in; the live kind becomes logical.
  !===================================================================!

  pure subroutine vertex_field_set_logical_vector(this, values)

    class(vertex_field), intent(inout) :: this
    logical            , intent(in)    :: values(:)

    this % lvals = values
    this % vkind = GRAPH_FIELD_LOGICAL

  end subroutine vertex_field_set_logical_vector

  !===================================================================!
  ! The words, when the live kind is character; zero-length otherwise.
  !===================================================================!

  pure subroutine vertex_field_get_character_vector(this, values)

    class(vertex_field), intent(in)                :: this
    character(len=:), allocatable, intent(out)     :: values(:)

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % svals)) then
       values = this % svals
    else
       allocate(character(len=1) :: values(0))
    end if

  end subroutine vertex_field_get_character_vector

  !===================================================================!
  ! Take words in; the live kind becomes character.
  !===================================================================!

  pure subroutine vertex_field_set_character_vector(this, values)

    class(vertex_field), intent(inout) :: this
    character(len=*)   , intent(in)    :: values(:)

    this % svals = values
    this % vkind = GRAPH_FIELD_CHARACTER

  end subroutine vertex_field_set_character_vector

  !===================================================================!
  ! Edge field: identity and support.
  !===================================================================!

  pure function edge_field_name(this) result(name)

    class(edge_field), intent(in)  :: this
    character(len=:), allocatable  :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function edge_field_name

  !===================================================================!
  ! What the values are measured in; a dash when nobody said.
  !===================================================================!

  pure function edge_field_units(this) result(units)

    class(edge_field), intent(in)  :: this
    character(len=:), allocatable  :: units

    if (allocated(this % unit_name)) then
       units = this % unit_name
    else
       units = '-'
    end if

  end function edge_field_units

  !===================================================================!
  ! The set the values live on, answered as a plain support.
  !===================================================================!

  subroutine edge_field_support(this, support)

    class(edge_field), intent(in)                  :: this
    class(graph_support), allocatable, intent(out) :: support

    allocate(support, source=this % on)

  end subroutine edge_field_support

  !===================================================================!
  ! The same set, answered with its edge type on show.
  !===================================================================!

  subroutine edge_field_edge_support(this, support)

    class(edge_field), intent(in)                       :: this
    class(graph_edge_support), allocatable, intent(out) :: support

    allocate(support, source=this % on)

  end subroutine edge_field_edge_support

  !===================================================================!
  ! Edge field: shape and kind.
  !===================================================================!

  pure integer function edge_field_num_components(this)

    class(edge_field), intent(in) :: this

    edge_field_num_components = this % ncomp

  end function edge_field_num_components

  !===================================================================!
  ! One entry per edge of the support.
  !===================================================================!

  pure integer function edge_field_num_entries(this)

    class(edge_field), intent(in) :: this

    edge_field_num_entries = this % on % size()

  end function edge_field_num_entries

  !===================================================================!
  ! Which of the five kinds is live.
  !===================================================================!

  pure integer function edge_field_value_kind(this)

    class(edge_field), intent(in) :: this

    edge_field_value_kind = this % vkind

  end function edge_field_value_kind

  !===================================================================!
  ! Edge field: the adapters, under the same rule.
  !===================================================================!

  pure subroutine edge_field_get_integer_vector(this, values)

    class(edge_field), intent(in)     :: this
    integer, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_INTEGER .and. allocated(this % ivals)) then
       values = this % ivals
    else
       allocate(values(0))
    end if

  end subroutine edge_field_get_integer_vector

  !===================================================================!
  ! Take integers in; the live kind becomes integer.
  !===================================================================!

  pure subroutine edge_field_set_integer_vector(this, values)

    class(edge_field), intent(inout) :: this
    integer          , intent(in)    :: values(:)

    this % ivals = values
    this % vkind = GRAPH_FIELD_INTEGER

  end subroutine edge_field_set_integer_vector

  !===================================================================!
  ! The values, when the live kind is real; zero-length otherwise.
  !===================================================================!

  pure subroutine edge_field_get_real_vector(this, values)

    class(edge_field), intent(in)      :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_REAL .and. allocated(this % rvals)) then
       values = this % rvals
    else
       allocate(values(0))
    end if

  end subroutine edge_field_get_real_vector

  !===================================================================!
  ! Take reals in; the live kind becomes real.
  !===================================================================!

  pure subroutine edge_field_set_real_vector(this, values)

    class(edge_field), intent(inout) :: this
    real(dp)         , intent(in)    :: values(:)

    this % rvals = values
    this % vkind = GRAPH_FIELD_REAL

  end subroutine edge_field_set_real_vector

  !===================================================================!
  ! The values, when the live kind is complex; zero-length otherwise.
  !===================================================================!

  pure subroutine edge_field_get_complex_vector(this, values)

    class(edge_field), intent(in)         :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_COMPLEX .and. allocated(this % cvals)) then
       values = this % cvals
    else
       allocate(values(0))
    end if

  end subroutine edge_field_get_complex_vector

  !===================================================================!
  ! Take complex values in; the live kind becomes complex. This is
  ! the seat a complex-step perturbation rides into the engine.
  !===================================================================!

  pure subroutine edge_field_set_complex_vector(this, values)

    class(edge_field), intent(inout) :: this
    complex(dp)      , intent(in)    :: values(:)

    this % cvals = values
    this % vkind = GRAPH_FIELD_COMPLEX

  end subroutine edge_field_set_complex_vector

  !===================================================================!
  ! The values, when the live kind is logical; zero-length otherwise.
  !===================================================================!

  pure subroutine edge_field_get_logical_vector(this, values)

    class(edge_field), intent(in)     :: this
    logical, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_LOGICAL .and. allocated(this % lvals)) then
       values = this % lvals
    else
       allocate(values(0))
    end if

  end subroutine edge_field_get_logical_vector

  !===================================================================!
  ! Take logicals in; the live kind becomes logical.
  !===================================================================!

  pure subroutine edge_field_set_logical_vector(this, values)

    class(edge_field), intent(inout) :: this
    logical          , intent(in)    :: values(:)

    this % lvals = values
    this % vkind = GRAPH_FIELD_LOGICAL

  end subroutine edge_field_set_logical_vector

  !===================================================================!
  ! The words, when the live kind is character; zero-length otherwise.
  !===================================================================!

  pure subroutine edge_field_get_character_vector(this, values)

    class(edge_field), intent(in)              :: this
    character(len=:), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % svals)) then
       values = this % svals
    else
       allocate(character(len=1) :: values(0))
    end if

  end subroutine edge_field_get_character_vector

  !===================================================================!
  ! Take words in; the live kind becomes character.
  !===================================================================!

  pure subroutine edge_field_set_character_vector(this, values)

    class(edge_field), intent(inout) :: this
    character(len=*) , intent(in)    :: values(:)

    this % svals = values
    this % vkind = GRAPH_FIELD_CHARACTER

  end subroutine edge_field_set_character_vector

end module class_graph_field
