!=====================================================================!
! The concrete functional: a field whose domain has one entry.
!
! A functional is one value reduced from a field - a total, an
! objective, a norm, an answer to a yes-or-no question.
!
!      +----+----+----+----+                          +-----+
!      | q1 | q2 | q3 | q4 |  ---- reduce ---->       |  J  |
!      +----+----+----+----+                          +-----+
!
! By the tower's reading this is not a new kind of thing: it is the
! field at domain size one, and this type answers the whole field
! contract at that size. The vector adapters move arrays of length
! one; the scalar pairs below them are conveniences for the callers
! who hold the concrete type and want the value without the array.
!
! Complex is here so a complex-step derivative survives a reduction.
! The derivative is the imaginary part, and a real-only functional
! throws it away:
!
!      (2.0, 1e-20) + (3.0, 3e-20) = (5.0, 4e-20)
!                                           \
!                                            the number we were after
!
! Logical is here so a question such as "is this graph acyclic" comes
! back as true or false rather than as a one or a zero.
!
! This type is also what a reduction uses to carry a part-way answer.
! A running average needs a sum and a count, so both live here; the
! value it reports is the quotient. A reduction that needed something
! else again would write its own functional and keep it private.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_functional

  use iso_fortran_env    , only : dp => REAL64
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use graph_field_calculus, only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use graph_field_calculus, only : GRAPH_FIELD_CHARACTER
  use graph_field_calculus     , only : graph_functional
  use graph_fractal      , only : set_graph => graph

  implicit none

  private
  public :: functional

  !===================================================================!
  ! One value, of whichever kind was last set.
  !===================================================================!

  type, extends(graph_functional) :: functional

     ! The one-entry home, declared at construction so domain()
     ! answers one stable identity for the life of the functional.
     type(set_graph), private :: home

     character(len=:), allocatable :: label
     character(len=:), allocatable :: unit_name

     integer :: vkind = GRAPH_FIELD_REAL

     integer                       :: ival = 0
     real(dp)                      :: rval = 0.0_dp
     complex(dp)                   :: cval = (0.0_dp, 0.0_dp)
     logical                       :: lval = .false.
     character(len=:), allocatable :: sval

     !----------------------------------------------------------------!
     ! Work carried while a reduction is still running. A sum alone
     ! needs none of it; an average needs the tally, a norm needs the
     ! power it is gathering. These are not part of the answer and no
     ! caller outside a reduction should read them.
     !----------------------------------------------------------------!

     real(dp) :: tally  = 0.0_dp
     real(dp) :: weight = 0.0_dp

   contains

     !----------------------------------------------------------------!
     ! The field contract, answered at one entry.
     !----------------------------------------------------------------!

     procedure :: name           => functional_name
     procedure :: units          => functional_units
     procedure :: domain         => functional_domain
     procedure :: num_components => functional_num_components
     procedure :: num_entries    => functional_num_entries
     procedure :: value_kind     => functional_value_kind

     procedure :: get_integer_vector   => functional_get_integer_vector
     procedure :: set_integer_vector   => functional_set_integer_vector
     procedure :: get_real_vector      => functional_get_real_vector
     procedure :: set_real_vector      => functional_set_real_vector
     procedure :: get_complex_vector   => functional_get_complex_vector
     procedure :: set_complex_vector   => functional_set_complex_vector
     procedure :: get_logical_vector   => functional_get_logical_vector
     procedure :: set_logical_vector   => functional_set_logical_vector
     procedure :: get_character_vector => functional_get_character_vector
     procedure :: set_character_vector => functional_set_character_vector

     !----------------------------------------------------------------!
     ! Scalar conveniences beyond the contract, for callers holding
     ! the concrete type: the same ten doors without the array.
     !----------------------------------------------------------------!

     procedure :: get_integer_value
     procedure :: set_integer_value
     procedure :: get_real_value
     procedure :: set_real_value
     procedure :: get_complex_value
     procedure :: set_complex_value
     procedure :: get_logical_value
     procedure :: set_logical_value
     procedure :: get_character_value
     procedure :: set_character_value

  end type functional

  !===================================================================!
  ! Constructor. Name it; the value arrives through a setter, which
  ! is also what fixes the kind.
  !===================================================================!

  interface functional
     module procedure create
  end interface functional

contains

  !===================================================================!
  ! Build a functional that holds nothing yet.
  !===================================================================!

  type(functional) function create(label, unit_name) result(this)

    character(len=*), intent(in), optional :: label
    character(len=*), intent(in), optional :: unit_name

    if (present(label)) then
       this % label = label
    else
       this % label = ''
    end if

    if (present(unit_name)) then
       this % unit_name = unit_name
    else
       this % unit_name = '-'
    end if


    call this % home % declare()
  end function create

  !===================================================================!
  ! Identity.
  !===================================================================!

  pure function functional_name(this) result(name)

    class(functional), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function functional_name

  pure function functional_units(this) result(units)

    class(functional), intent(in) :: this
    character(len=:), allocatable :: units

    if (allocated(this % unit_name)) then
       units = this % unit_name
    else
       units = '-'
    end if

  end function functional_units

  !===================================================================!
  ! The terminal domain: one member, no edges. A single value has no
  ! side of its own; the terminal support answers vertex by
  ! convention, and nothing downstream reads it.
  !===================================================================!

  type(set_graph) function functional_domain(this) result(domain)

    class(functional), intent(in) :: this

    domain = this % home

  end function functional_domain

  !===================================================================!
  ! Shape: one entry, one component, one live kind.
  !===================================================================!

  pure integer function functional_num_components(this)

    class(functional), intent(in) :: this

    associate (u1 => this); end associate

    functional_num_components = 1

  end function functional_num_components

  pure integer function functional_num_entries(this)

    class(functional), intent(in) :: this

    associate (u1 => this); end associate

    functional_num_entries = 1

  end function functional_num_entries

  pure integer function functional_value_kind(this)

    class(functional), intent(in) :: this

    functional_value_kind = this % vkind

  end function functional_value_kind

  !===================================================================!
  ! The vector adapters, at length one. The same three rules as any
  ! field: ask first, a wrong getter answers zero-length, any setter
  ! replaces the value and the kind together. A setter handed an
  ! empty array changes nothing.
  !===================================================================!

  pure subroutine functional_get_integer_vector(this, values)

    class(functional), intent(in)     :: this
    integer, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_INTEGER) then
       values = [this % ival]
    else
       allocate(values(0))
    end if

  end subroutine functional_get_integer_vector

  pure subroutine functional_set_integer_vector(this, values)

    class(functional), intent(inout) :: this
    integer          , intent(in)    :: values(:)

    if (size(values) >= 1) then
       this % ival  = values(1)
       this % vkind = GRAPH_FIELD_INTEGER
    end if

  end subroutine functional_set_integer_vector

  pure subroutine functional_get_real_vector(this, values)

    class(functional), intent(in)      :: this
    real(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_REAL) then
       values = [this % rval]
    else
       allocate(values(0))
    end if

  end subroutine functional_get_real_vector

  pure subroutine functional_set_real_vector(this, values)

    class(functional), intent(inout) :: this
    real(dp)         , intent(in)    :: values(:)

    if (size(values) >= 1) then
       this % rval  = values(1)
       this % vkind = GRAPH_FIELD_REAL
    end if

  end subroutine functional_set_real_vector

  pure subroutine functional_get_complex_vector(this, values)

    class(functional), intent(in)         :: this
    complex(dp), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_COMPLEX) then
       values = [this % cval]
    else
       allocate(values(0))
    end if

  end subroutine functional_get_complex_vector

  pure subroutine functional_set_complex_vector(this, values)

    class(functional), intent(inout) :: this
    complex(dp)      , intent(in)    :: values(:)

    if (size(values) >= 1) then
       this % cval  = values(1)
       this % vkind = GRAPH_FIELD_COMPLEX
    end if

  end subroutine functional_set_complex_vector

  pure subroutine functional_get_logical_vector(this, values)

    class(functional), intent(in)     :: this
    logical, allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_LOGICAL) then
       values = [this % lval]
    else
       allocate(values(0))
    end if

  end subroutine functional_get_logical_vector

  pure subroutine functional_set_logical_vector(this, values)

    class(functional), intent(inout) :: this
    logical          , intent(in)    :: values(:)

    if (size(values) >= 1) then
       this % lval  = values(1)
       this % vkind = GRAPH_FIELD_LOGICAL
    end if

  end subroutine functional_set_logical_vector

  pure subroutine functional_get_character_vector(this, values)

    class(functional), intent(in)              :: this
    character(len=:), allocatable, intent(out) :: values(:)

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % sval)) then
       allocate(character(len=len(this % sval)) :: values(1))
       values(1) = this % sval
    else
       allocate(character(len=1) :: values(0))
    end if

  end subroutine functional_get_character_vector

  pure subroutine functional_set_character_vector(this, values)

    class(functional), intent(inout) :: this
    character(len=*) , intent(in)    :: values(:)

    if (size(values) >= 1) then
       this % sval  = values(1)
       this % vkind = GRAPH_FIELD_CHARACTER
    end if

  end subroutine functional_set_character_vector

  !===================================================================!
  ! The scalar doors. A getter for the wrong kind answers the zero
  ! of the asked kind; a setter replaces the value and the kind.
  !===================================================================!

  pure subroutine get_integer_value(this, value)

    class(functional), intent(in) :: this
    integer, intent(out) :: value

    if (this % vkind == GRAPH_FIELD_INTEGER) then
       value = this % ival
    else
       value = 0
    end if

  end subroutine get_integer_value

  pure subroutine set_integer_value(this, value)

    class(functional), intent(inout) :: this
    integer, intent(in) :: value

    this % ival  = value
    this % vkind = GRAPH_FIELD_INTEGER

  end subroutine set_integer_value

  pure subroutine get_real_value(this, value)

    class(functional), intent(in) :: this
    real(dp), intent(out) :: value

    if (this % vkind == GRAPH_FIELD_REAL) then
       value = this % rval
    else
       value = 0.0_dp
    end if

  end subroutine get_real_value

  pure subroutine set_real_value(this, value)

    class(functional), intent(inout) :: this
    real(dp), intent(in) :: value

    this % rval  = value
    this % vkind = GRAPH_FIELD_REAL

  end subroutine set_real_value

  pure subroutine get_complex_value(this, value)

    class(functional), intent(in) :: this
    complex(dp), intent(out) :: value

    if (this % vkind == GRAPH_FIELD_COMPLEX) then
       value = this % cval
    else
       value = (0.0_dp, 0.0_dp)
    end if

  end subroutine get_complex_value

  pure subroutine set_complex_value(this, value)

    class(functional), intent(inout) :: this
    complex(dp), intent(in) :: value

    this % cval  = value
    this % vkind = GRAPH_FIELD_COMPLEX

  end subroutine set_complex_value

  pure subroutine get_logical_value(this, value)

    class(functional), intent(in) :: this
    logical, intent(out) :: value

    if (this % vkind == GRAPH_FIELD_LOGICAL) then
       value = this % lval
    else
       value = .false.
    end if

  end subroutine get_logical_value

  pure subroutine set_logical_value(this, value)

    class(functional), intent(inout) :: this
    logical, intent(in) :: value

    this % lval  = value
    this % vkind = GRAPH_FIELD_LOGICAL

  end subroutine set_logical_value

  pure subroutine get_character_value(this, value)

    class(functional), intent(in) :: this
    character(len=:), allocatable, intent(out) :: value

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % sval)) then
       value = this % sval
    else
       value = ''
    end if

  end subroutine get_character_value

  pure subroutine set_character_value(this, value)

    class(functional), intent(inout) :: this
    character(len=*), intent(in) :: value

    this % sval  = value
    this % vkind = GRAPH_FIELD_CHARACTER

  end subroutine set_character_value

end module class_graph_functional
