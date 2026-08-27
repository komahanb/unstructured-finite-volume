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

     !----------------------------------------------------------------!
     ! Work carried while a reduction is still running. A sum alone
     ! needs none of it; an average needs the tally, a norm needs the
     ! power it is gathering. These are not part of the answer and no
     ! caller outside a reduction should read them.
     !----------------------------------------------------------------!

     real(dp) :: tally  = 0.0_dp
     real(dp) :: weight = 0.0_dp

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

    type(graph) :: home
    character(len=:), allocatable :: called, measured

    called = ''
    if (present(label)) called = label
    measured = '-'
    if (present(unit_name)) measured = unit_name

    call home % declare()
    call this % describe(called, home, 1, 1, measured)

  end function create

end module field_functional
