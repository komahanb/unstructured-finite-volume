!=====================================================================!
! The concrete field: values over a domain.
!
! One concrete type serves every field in the tower. Its domain is a
! set GRAPH, and the domain's identity is the only thing that ever
! distinguishes a cell field from a face field; the field carries no
! side flag. Because there is exactly one concrete field, a plain
! Fortran array can hold a collection of them.
!
!                  WHAT THE FIELD KEEPS OF ITS DOMAIN
!
!      type(graph) :: on          which set        O(1)
!      integer     :: num_entries    how many         O(1)
!
! and nothing else. It once kept a COPY of the whole domain object,
! which for a listed domain meant a copy of the member roll: 40 fields
! on a 200 000-member domain carried 28.7 MB of duplicated extension,
! measured, against 30.5 MB predicted for exactly that duplication.
! The extension now lives once, in whatever set map the caller holds.
!
! Nothing was lost with it. A field only ever answered two questions
! about its domain - WHICH and HOW MANY - and both survive by value.
! The copy already froze the count, so freezing it explicitly changes
! no behaviour; it only stops the copy from being O(N_extent).
!
!=====================================================================!
!
!                        THE VALUE-KIND RULE
!
! A field holds one kind of value at a time, in the one store every
! field inherits from field_calculus, where the ten adapters are
! written once. From that, three rules that hold for all of them:
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

module field_stored

  use util_precision  , only : dp
  use field_calculus, only : field
  use graph_fractal , only : graph

  implicit none

  private
  public :: stored_field

  !===================================================================!
  ! One field: the description and the store the abstract field
  ! holds, stated by the constructor and nothing more.
  !===================================================================!

  type, extends(field) :: stored_field
  end type stored_field

  !===================================================================!
  ! Constructor. Name the field, say where it lives, say how wide
  ! each entry is. The values arrive afterwards through a setter,
  ! which is also what fixes the kind.
  !===================================================================!

  interface stored_field
     module procedure create
  end interface stored_field

contains

  !===================================================================!
  ! Build an empty field on a domain. The domain's identity says
  ! whether this is a cell field or a face field; the field does not
  ! repeat the fact.
  !===================================================================!

  type(stored_field) function create(label, on, num_entries, num_components, unit_name) &
       & result(this)

    character(len=*), intent(in)           :: label
    type(graph) , intent(in)           :: on
    integer         , intent(in)           :: num_entries
    integer         , intent(in), optional :: num_components
    character(len=*), intent(in), optional :: unit_name

    call this % describe(label, on, num_entries, num_components, unit_name)

  end function create

end module field_stored
