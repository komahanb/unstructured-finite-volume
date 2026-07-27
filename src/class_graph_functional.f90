!=====================================================================!
! A concrete graph functional.
!
! A functional is one value reduced from a field - a total, an
! objective, a norm, an answer to a yes-or-no question.
!
!      +----+----+----+----+                          +-----+
!      | q1 | q2 | q3 | q4 |  ---- reduce ---->        |  J  |
!      +----+----+----+----+                          +-----+
!
! It carries every kind of value a field can, under the same rule
! stated in class_graph_field: only the matching store is live, a
! getter for the wrong kind answers with the zero of that kind, and
! any setter replaces the value and the kind together.
!
! Complex is here so a complex-step derivative survives the fold. The
! derivative is the imaginary part, and a real-only functional throws
! it away:
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

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_functional
  use abstract_graph_types, only : GRAPH_FIELD_INTEGER, GRAPH_FIELD_REAL
  use abstract_graph_types, only : GRAPH_FIELD_COMPLEX, GRAPH_FIELD_LOGICAL
  use abstract_graph_types, only : GRAPH_FIELD_CHARACTER

  implicit none

  private
  public :: functional

  !===================================================================!
  ! One value, of whichever kind was last set.
  !===================================================================!

  type, extends(graph_functional) :: functional

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
     ! Identity.
     !----------------------------------------------------------------!

     procedure :: name  => fn_name
     procedure :: units => fn_units

     !----------------------------------------------------------------!
     ! Which kind of value is in here.
     !----------------------------------------------------------------!

     procedure :: value_kind => fn_value_kind

     !----------------------------------------------------------------!
     ! One value in and out, per kind.
     !----------------------------------------------------------------!

     procedure :: get_integer_value   => fn_get_integer
     procedure :: set_integer_value   => fn_set_integer
     procedure :: get_real_value      => fn_get_real
     procedure :: set_real_value      => fn_set_real
     procedure :: get_complex_value   => fn_get_complex
     procedure :: set_complex_value   => fn_set_complex
     procedure :: get_logical_value   => fn_get_logical
     procedure :: set_logical_value   => fn_set_logical
     procedure :: get_character_value => fn_get_character
     procedure :: set_character_value => fn_set_character

  end type functional

  !===================================================================!
  ! Constructor. Name it; the value arrives through a setter, which is
  ! also what fixes the kind.
  !===================================================================!

  interface functional
     module procedure create_functional
  end interface functional

contains

  !===================================================================!
  ! Build a functional that holds nothing yet.
  !===================================================================!

  pure type(functional) function create_functional(label, unit_name) result(this)

    character(len=*), intent(in)           :: label
    character(len=*), intent(in), optional :: unit_name

    this % label = label

    if (present(unit_name)) then
       this % unit_name = unit_name
    else
       this % unit_name = '-'
    end if

  end function create_functional

  !===================================================================!
  ! Identity.
  !===================================================================!

  pure function fn_name(this) result(name)

    class(functional), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function fn_name

  pure function fn_units(this) result(units)

    class(functional), intent(in) :: this
    character(len=:), allocatable :: units

    if (allocated(this % unit_name)) then
       units = this % unit_name
    else
       units = '-'
    end if

  end function fn_units

  pure integer function fn_value_kind(this)

    class(functional), intent(in) :: this

    fn_value_kind = this % vkind

  end function fn_value_kind

  !===================================================================!
  ! The accessors. A getter for the wrong kind answers with the zero
  ! of that kind rather than converting or guessing; it is pure and
  ! has no way to complain. A caller checks value_kind() first.
  !===================================================================!

  pure subroutine fn_get_integer(this, value)

    class(functional), intent(in) :: this
    integer          , intent(out) :: value

    if (this % vkind == GRAPH_FIELD_INTEGER) then
       value = this % ival
    else
       value = 0
    end if

  end subroutine fn_get_integer

  pure subroutine fn_set_integer(this, value)

    class(functional), intent(inout) :: this
    integer          , intent(in)    :: value

    this % ival  = value
    this % vkind = GRAPH_FIELD_INTEGER

  end subroutine fn_set_integer

  pure subroutine fn_get_real(this, value)

    class(functional), intent(in)  :: this
    real(dp)         , intent(out) :: value

    if (this % vkind == GRAPH_FIELD_REAL) then
       value = this % rval
    else
       value = 0.0_dp
    end if

  end subroutine fn_get_real

  pure subroutine fn_set_real(this, value)

    class(functional), intent(inout) :: this
    real(dp)         , intent(in)    :: value

    this % rval  = value
    this % vkind = GRAPH_FIELD_REAL

  end subroutine fn_set_real

  pure subroutine fn_get_complex(this, value)

    class(functional), intent(in)  :: this
    complex(dp)      , intent(out) :: value

    if (this % vkind == GRAPH_FIELD_COMPLEX) then
       value = this % cval
    else
       value = (0.0_dp, 0.0_dp)
    end if

  end subroutine fn_get_complex

  pure subroutine fn_set_complex(this, value)

    class(functional), intent(inout) :: this
    complex(dp)      , intent(in)    :: value

    this % cval  = value
    this % vkind = GRAPH_FIELD_COMPLEX

  end subroutine fn_set_complex

  pure subroutine fn_get_logical(this, value)

    class(functional), intent(in)  :: this
    logical          , intent(out) :: value

    if (this % vkind == GRAPH_FIELD_LOGICAL) then
       value = this % lval
    else
       value = .false.
    end if

  end subroutine fn_get_logical

  pure subroutine fn_set_logical(this, value)

    class(functional), intent(inout) :: this
    logical          , intent(in)    :: value

    this % lval  = value
    this % vkind = GRAPH_FIELD_LOGICAL

  end subroutine fn_set_logical

  pure subroutine fn_get_character(this, value)

    class(functional), intent(in)              :: this
    character(len=:), allocatable, intent(out) :: value

    if (this % vkind == GRAPH_FIELD_CHARACTER .and. allocated(this % sval)) then
       value = this % sval
    else
       value = ''
    end if

  end subroutine fn_get_character

  pure subroutine fn_set_character(this, value)

    class(functional), intent(inout) :: this
    character(len=*) , intent(in)    :: value

    this % sval  = value
    this % vkind = GRAPH_FIELD_CHARACTER

  end subroutine fn_set_character

end module class_graph_functional
