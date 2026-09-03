!=====================================================================!
! The concrete field: values over a domain.
!
! One concrete type serves every field in the tower. Its domain is a
! set GRAPH, and the domain's identity is the only datum that ever
! distinguishes a cell field from a face field; the field stores no
! side flag. Because there is exactly one concrete field, a plain
! Fortran array can store a collection of them.
!
!                  WHAT THE FIELD STORES OF ITS DOMAIN
!
!      type(graph) :: graph       which set        O(1)
!      integer     :: num_entries    how many         O(1)
!
! and nothing else. An earlier version stored a COPY of the whole
! domain object, which for a listed domain meant a copy of the member
! list: 40 fields on a 200 000-member domain stored 28.7 MB of
! duplicated extension, measured, against 30.5 MB predicted for
! exactly that duplication. The extension is now stored once, in
! whatever set map the caller owns.
!
! No information was lost. A field only ever returned two values
! about its domain - WHICH and HOW MANY - and both are retained by
! value. The copy already froze the count, so freezing it explicitly
! changes no behaviour; it only stops the copy from costing
! O(N_extent).
!
!=====================================================================!
!
!                        THE VALUE-KIND RULE
!
! A field stores one kind of value at a time, in the one store every
! field inherits from field_calculus, where the ten adapters are
! written once. From that, three rules that hold for all of them:
!
!      check first   a caller checks value_kind() before reading a
!                    vector
!
!      wrong getter  returns a zero-length array. No conversion and
!                    no inference happens, and a pure procedure has
!                    no error path. The zero-length result is the
!                    indicator
!
!      any setter    replaces both the values and the kind. Setting
!                    reals onto a field that stored integers makes it
!                    a real field
!
! No conversion happens anywhere. A field that stores boundary names
! does not return them as numbers.
!
!=====================================================================!
!
!                       WHERE A VALUE IS STORED
!
! A field stores its values in the order the domain lists its
! members, and stores the components of one member contiguously:
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
! this library needs no degree-of-freedom index map of its own.
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
  ! contains, stated by the constructor and nothing more.
  !===================================================================!

  type, extends(field) :: stored_field
   contains
     procedure :: assign_in => stored_field_assign_in
  end type stored_field

  !===================================================================!
  ! Constructor. Name the field, state its domain, state how many
  ! components each entry has. The values are set afterwards through
  ! a setter, which is also what fixes the kind.
  !===================================================================!

  interface stored_field
     module procedure create
  end interface stored_field

contains

  !===================================================================!
  ! Build an empty field on a domain. The domain's identity states
  ! whether this is a cell field or a face field; the field does not
  ! store the fact a second time.
  !===================================================================!

  type(stored_field) function create(label, domain, num_entries, num_components, unit_name) &
       & result(this)

    character(len=*), intent(in)           :: label
    type(graph) , intent(in)           :: domain
    integer         , intent(in)           :: num_entries
    integer         , intent(in), optional :: num_components
    character(len=*), intent(in), optional :: unit_name

    call this % describe(label, domain, num_entries, num_components, unit_name)

  end function create

  !===================================================================!
  ! Place this value at a location that is a stored field. A location
  ! of any other type is an error: the caller requested a copy the
  ! location cannot store.
  !===================================================================!

  subroutine stored_field_assign_in(this, location)

    class(stored_field), intent(in)    :: this
    class(field)       , intent(inout) :: location

    select type (location)
    type is (stored_field)
       location = this
    class default
       error stop 'field_stored: a stored field is placed at a stored field'
    end select

  end subroutine stored_field_assign_in

end module field_stored
