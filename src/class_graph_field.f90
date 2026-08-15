!=====================================================================!
! The concrete field: values over a domain.
!
! One concrete type serves every field in the tower. Its domain is a
! set - an ambient set or a subset subobject - and the
! domain's identity is the only thing that ever distinguishes a
! cell field from a face field; the field carries no side flag. Because there is
! exactly one concrete field, a plain Fortran array can hold a
! collection of them.
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
! From that, three rules that hold for all ten adapters:
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
! A field stores its values in the order the domain lists its
! members, and keeps the components of one member next to each other:
!
!      member          7        7        3        3
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

  use iso_fortran_env    , only : dp => REAL64
  use graph_set       , only : set
  use graph_field_calculus, only : graph_field
  use graph_field_calculus, only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use graph_field_calculus, only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use graph_field_calculus, only : GRAPH_FIELD_CHARACTER

  implicit none

  private
  public :: field

  !===================================================================!
  ! One field: a name, a unit, a domain, a width, and one live store.
  !===================================================================!

  type, extends(graph_field) :: field

     character(len=:), allocatable :: label
     character(len=:), allocatable :: unit_name

     ! The one domain: an ambient set or a subset subobject,
     ! copied at construction so its identity travels. Private, so
     ! consumers ask through domain() rather than inspecting on.
     class(set), allocatable, private :: on

     integer :: ncomp = 1
     integer :: vkind = GRAPH_FIELD_REAL

     integer         , allocatable :: ivals(:)
     real(dp)        , allocatable :: rvals(:)
     complex(dp)     , allocatable :: cvals(:)
     logical         , allocatable :: lvals(:)
     character(len=:), allocatable :: svals(:)

   contains

     !----------------------------------------------------------------!
     ! Identity, and where the values live.
     !----------------------------------------------------------------!

     procedure :: name   => field_name
     procedure :: units  => field_units
     procedure :: domain => field_domain

     !----------------------------------------------------------------!
     ! Shape and kind.
     !----------------------------------------------------------------!

     procedure :: num_components => field_num_components
     procedure :: num_entries    => field_num_entries
     procedure :: value_kind     => field_value_kind

     !----------------------------------------------------------------!
     ! The plain-vector adapters, one pair per kind.
     !----------------------------------------------------------------!

     procedure :: get_integer_vector   => field_get_integer_vector
     procedure :: set_integer_vector   => field_set_integer_vector
     procedure :: get_real_vector      => field_get_real_vector
     procedure :: set_real_vector      => field_set_real_vector
     procedure :: get_complex_vector   => field_get_complex_vector
     procedure :: set_complex_vector   => field_set_complex_vector
     procedure :: get_logical_vector   => field_get_logical_vector
     procedure :: set_logical_vector   => field_set_logical_vector
     procedure :: get_character_vector => field_get_character_vector
     procedure :: set_character_vector => field_set_character_vector

  end type field

  !===================================================================!
  ! Constructor. Name the field, say where it lives, say how wide
  ! each entry is. The values arrive afterwards through a setter,
  ! which is also what fixes the kind.
  !===================================================================!

  interface field
     module procedure create
  end interface field

contains

  !===================================================================!
  ! Build an empty field on a domain. The domain's identity says
  ! whether this is a cell field or a face field; the field does not
  ! repeat the fact.
  !===================================================================!

  type(field) function create(label, on, ncomp, unit_name) result(this)

    character(len=*), intent(in)           :: label
    class(set), intent(in)          :: on
    integer         , intent(in), optional :: ncomp
    character(len=*), intent(in), optional :: unit_name

    if (.not. on % equals(on)) then
       error stop 'class_graph_field: a field needs a declared domain'
    end if

    this % label = label
    allocate(this % on, source=on)

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

  end function create

  !===================================================================!
  ! What the field is called; an empty name when nobody named it.
  !===================================================================!

  pure function field_name(this) result(name)

    class(field), intent(in)      :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function field_name

  !===================================================================!
  ! What the values are measured in; a dash when nobody said.
  !===================================================================!

  pure function field_units(this) result(units)

    class(field), intent(in)      :: this
    character(len=:), allocatable :: units

    if (allocated(this % unit_name)) then
       units = this % unit_name
    else
       units = '-'
    end if

  end function field_units

  !===================================================================!
  ! The set domain the values live on.
  !===================================================================!

  subroutine field_domain(this, domain)

    class(field), intent(in)                    :: this
    class(set), allocatable, intent(out) :: domain

    allocate(domain, source=this % on)

  end subroutine field_domain

  !===================================================================!
  ! Shape and kind. The entry count comes from the domain, so a
  ! field cannot disagree with the set it lives on.
  !===================================================================!

  pure integer function field_num_components(this)

    class(field), intent(in) :: this

    field_num_components = this % ncomp

  end function field_num_components

  pure integer function field_num_entries(this)

    class(field), intent(in) :: this

    field_num_entries = this % on % size()

  end function field_num_entries

  pure integer function field_value_kind(this)

    class(field), intent(in) :: this

    field_value_kind = this % vkind

  end function field_value_kind

  !===================================================================!
  ! The adapters. A getter for the wrong kind answers with a
  ! zero-length array; a setter replaces the values and the kind
  ! together.
  !===================================================================!

  pure subroutine field_get_integer_vector(this, values)

    class(field), intent(in)          :: this
    integer, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_INTEGER .and. allocated(this % ivals)) then
       values = this % ivals
    else
       allocate(values(0))
    end if

  end subroutine field_get_integer_vector

  pure subroutine field_set_integer_vector(this, values)

    class(field), intent(inout) :: this
    integer     , intent(in)    :: values(:)

    if (size(values) /= this % on % size() * this % ncomp) then
       error stop 'class_graph_field: a value vector must fill its domain exactly'
    end if
    this % ivals = values
    this % vkind = GRAPH_FIELD_INTEGER

  end subroutine field_set_integer_vector

  pure subroutine field_get_real_vector(this, values)

    class(field), intent(in)           :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_REAL .and. allocated(this % rvals)) then
       values = this % rvals
    else
       allocate(values(0))
    end if

  end subroutine field_get_real_vector

  pure subroutine field_set_real_vector(this, values)

    class(field), intent(inout) :: this
    real(dp)    , intent(in)    :: values(:)

    if (size(values) /= this % on % size() * this % ncomp) then
       error stop 'class_graph_field: a value vector must fill its domain exactly'
    end if
    this % rvals = values
    this % vkind = GRAPH_FIELD_REAL

  end subroutine field_set_real_vector

  !===================================================================!
  ! The complex pair is the seat a complex-step perturbation rides
  ! into the engine.
  !===================================================================!

  pure subroutine field_get_complex_vector(this, values)

    class(field), intent(in)              :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_COMPLEX .and. allocated(this % cvals)) then
       values = this % cvals
    else
       allocate(values(0))
    end if

  end subroutine field_get_complex_vector

  pure subroutine field_set_complex_vector(this, values)

    class(field), intent(inout) :: this
    complex(dp) , intent(in)    :: values(:)

    if (size(values) /= this % on % size() * this % ncomp) then
       error stop 'class_graph_field: a value vector must fill its domain exactly'
    end if
    this % cvals = values
    this % vkind = GRAPH_FIELD_COMPLEX

  end subroutine field_set_complex_vector

  pure subroutine field_get_logical_vector(this, values)

    class(field), intent(in)          :: this
    logical, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_LOGICAL .and. allocated(this % lvals)) then
       values = this % lvals
    else
       allocate(values(0))
    end if

  end subroutine field_get_logical_vector

  pure subroutine field_set_logical_vector(this, values)

    class(field), intent(inout) :: this
    logical     , intent(in)    :: values(:)

    if (size(values) /= this % on % size() * this % ncomp) then
       error stop 'class_graph_field: a value vector must fill its domain exactly'
    end if
    this % lvals = values
    this % vkind = GRAPH_FIELD_LOGICAL

  end subroutine field_set_logical_vector

  pure subroutine field_get_character_vector(this, values)

    class(field), intent(in)                   :: this
    character(len=:), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % svals)) then
       values = this % svals
    else
       allocate(character(len=1) :: values(0))
    end if

  end subroutine field_get_character_vector

  pure subroutine field_set_character_vector(this, values)

    class(field), intent(inout)  :: this
    character(len=*), intent(in) :: values(:)

    if (size(values) /= this % on % size() * this % ncomp) then
       error stop 'class_graph_field: a value vector must fill its domain exactly'
    end if
    this % svals = values
    this % vkind = GRAPH_FIELD_CHARACTER

  end subroutine field_set_character_vector

end module class_graph_field
