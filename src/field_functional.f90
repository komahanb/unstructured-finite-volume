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

module field_functional

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use field_calculus     , only : functional
  use graph_fractal      , only : graph

  implicit none

  private
  public :: stored_functional

  !===================================================================!
  ! One value, of whichever kind was last set, held in the store every
  ! field inherits; the adapters are field_calculus's, at length one.
  !===================================================================!

  type, extends(functional) :: stored_functional

     ! The one-entry home, declared at construction so domain()
     ! answers one stable identity for the life of the functional.
     type(graph), private :: home

     character(len=:), allocatable :: label
     character(len=:), allocatable :: unit_name



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


  end type stored_functional

  !===================================================================!
  ! Constructor. Name it; the value arrives through a setter, which
  ! is also what fixes the kind.
  !===================================================================!

  interface stored_functional
     module procedure create
  end interface stored_functional

contains

  !===================================================================!
  ! Build a functional that holds nothing yet.
  !===================================================================!

  type(stored_functional) function create(label, unit_name) result(this)

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

    class(stored_functional), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function functional_name

  pure function functional_units(this) result(units)

    class(stored_functional), intent(in) :: this
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

  type(graph) function functional_domain(this) result(domain)

    class(stored_functional), intent(in) :: this

    domain = this % home

  end function functional_domain

  !===================================================================!
  ! Shape: one entry, one component, one live kind.
  !===================================================================!

  pure integer function functional_num_components(this)

    class(stored_functional), intent(in) :: this

    associate (u1 => this); end associate

    functional_num_components = 1

  end function functional_num_components

  pure integer function functional_num_entries(this)

    class(stored_functional), intent(in) :: this

    associate (u1 => this); end associate

    functional_num_entries = 1

  end function functional_num_entries

end module field_functional
