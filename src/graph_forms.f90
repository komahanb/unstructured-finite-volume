!=====================================================================!
! LEVEL 1 . THE FORMS
!
! A form is a family of functions of position - a basis shape - and
! the tower now says what that IS: a support whose members are basis
! functions. Concretion continued in the graph direction, abstraction
! reopened in the evaluation direction - the zoom, done deliberately:
!
!      graph -> graph_support -> support -> form -> polynomial | wave
!
! Everything a roster once did, membership does: the standing basis
! members ARE the support's members, indices into the concretion's
! own table of functions. Pruning a form is building a smaller
! member set - a graph act, owned by the form sector one level up.
! No active(:) array survives; a set does not need a second list to
! say who belongs to it.
!
! What the form adds beyond membership is only its evaluation
! symbols, read over the FULL table, membership saying who stands:
!
!      size_of                  the table's width
!      values(x, at)            each table entry, evaluated at x,
!                               reckoned about the point `at`
!      slopes(x, at, n)         each entry's derivative along n
!
! and one act of its own: restrict, which sets that membership. It
! is here rather than at the caller because a citizen's structure is
! its own business - whoever decides a member should go says so, and
! the form does it. When the form sector becomes a transform the
! restriction will hand back a NEW form and this verb becomes the
! constructor it calls.
!
! Evaluating a form at a point is calculus; choosing its
! coefficients is minimization and lives one level up. Polynomials
! are one concretion, waves another; a fit holds a form the way an
! operator holds coefficients - as data about shape.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_forms

  use iso_fortran_env    , only : dp => REAL64
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use graph_set, only : subset, set

  implicit none

  private
  public :: form

  type, abstract, extends(subset) :: form

   contains

     procedure(form_size_interface)  , deferred :: size_of
     procedure(form_values_interface), deferred :: values
     procedure(form_slopes_interface), deferred :: slopes

     procedure :: restrict

  end type form

  abstract interface

     pure integer function form_size_interface(this)
       import :: form
       class(form), intent(in) :: this
     end function form_size_interface

     pure subroutine form_values_interface(this, x, at, phi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(3), at(3)
       real(dp), intent(out) :: phi(:)
     end subroutine form_values_interface

     pure subroutine form_slopes_interface(this, x, at, direction, dphi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(3), at(3), direction(3)
       real(dp), intent(out) :: dphi(:)
     end subroutine form_slopes_interface

  end interface

contains

  !===================================================================!
  ! Stand only these table entries. The kept indices name members of
  ! the concretion's own table, and membership is the whole roster -
  ! a set says who belongs by holding them.
  !===================================================================!

  subroutine restrict(this, kept)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: kept(:)

    class(set), allocatable :: table

    table = this % ambient()
    this % subset = subset('basis', table, kept)

  end subroutine restrict

end module graph_forms
