!=====================================================================!
! Module that states the concrete functional: a field whose domain has
! one entry.
!
! A functional is one value reduced from a field: a total, an
! objective, a norm, the value of a predicate.
!
!      +----+----+----+----+                          +-----+
!      | q1 | q2 | q3 | q4 |  ---- reduce ---->       |  J  |
!      +----+----+----+----+                          +-----+
!
! By the tower's definition it is not a new kind of object: it is the
! field at domain size one, and this type implements the whole field
! contract at that size. The vector adapters move arrays of length
! one; the scalar pairs below them are conveniences for the callers
! that reference the concrete type and want the value without the
! array.
!
! Complex is here so that a complex-step derivative is preserved
! through a reduction. The derivative is the imaginary part, and a
! real-only functional discards it:
!
!      (2.0, 1e-20) + (3.0, 3e-20) = (5.0, 4e-20)
!                                           \
!                                            the derivative
!
! Logical is here so that a predicate, such as "is this graph
! acyclic", is returned as true or false rather than as one or zero.
!
! This type is also what a reduction uses to store an intermediate
! result. A running average needs a sum and a count, so both are
! stored here, and the value it reports is the quotient. A reduction
! that needs other intermediate data defines its own functional and
! keeps it private.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module field_functional

  use util_precision, only : dp
  use view_directed , only : directed_graph
  use field_calculus, only : functional, field
  use graph_fractal , only : graph

  implicit none

  private

  public :: stored_functional

  !-------------------------------------------------------------------!
  ! Type that models one value, of whichever kind was last set, in
  ! the store every field inherits; the adapters are field_calculus's,
  ! at length one
  !-------------------------------------------------------------------!

  type, extends(functional) :: stored_functional

     ! The intermediate data of a running reduction. A sum needs none
     ! of it; an average needs the tally, a norm the accumulated power.
     ! They are not part of the result, and no caller outside a
     ! reduction reads them.
     real(dp) :: tally  = 0.0_dp
     real(dp) :: weight = 0.0_dp

   contains

     procedure :: assign_in => stored_functional_assign_in

  end type stored_functional

  ! The constructor sets the label; the value is set through a setter,
  ! which also fixes its kind
  interface stored_functional
     module procedure create
  end interface stored_functional

contains

  !===================================================================!
  ! Creates a functional that stores no value yet, with its label and
  ! unit name ('-' when none is given)
  !===================================================================!

  type(stored_functional) function create(label, unit_name) result(this)

    character(len=*), intent(in), optional :: label
    character(len=*), intent(in), optional :: unit_name

    type(graph)                   :: domain
    character(len=:), allocatable :: field_label, field_units

    field_label = ''
    if (present(label)) field_label = label
    field_units = '-'
    if (present(unit_name)) field_units = unit_name

    ! A new domain of one entry
    call domain % declare()
    call this % describe(field_label, domain, 1, 1, field_units)

  end function create

  !===================================================================!
  ! Places this value at a location that is a stored functional. A
  ! location of any other type is an error: the caller requested a copy
  ! the location cannot store.
  !===================================================================!

  subroutine stored_functional_assign_in(this, location)

    class(stored_functional), intent(in)    :: this
    class(field)            , intent(inout) :: location

    select type (location)
    type is (stored_functional)
       location = this
    class default
       error stop 'field_functional: location''s dynamic type is not type(stored_functional); &
            &a stored functional can only be placed at a stored functional'
    end select

  end subroutine stored_functional_assign_in

end module field_functional
