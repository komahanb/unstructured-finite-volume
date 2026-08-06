!=====================================================================!
! LEVEL 1 . THE FORMS
!
! A graph_form is a family of functions of position - a basis shape - and
! the tower now says what that IS: a support whose members are basis
! functions. Concretion continued in the graph direction, abstraction
! reopened in the evaluation direction - the zoom, done deliberately:
!
!      graph -> graph_support -> support -> graph_form -> polynomial | wave
!
! Everything a roster once did, membership does: the standing basis
! members ARE the support's members, indices into the concretion's
! own table of functions. Pruning a graph_form is building a smaller
! member set - a graph act, owned by the graph_form sector one level up.
! No active(:) array survives; a set does not need a second list to
! say who belongs to it.
!
! What the graph_form adds beyond membership is only its evaluation
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
! the graph_form does it. When the graph_form sector becomes a transform the
! restriction will hand back a NEW graph_form and this verb becomes the
! constructor it calls.
!
! Evaluating a graph_form at a point is calculus; choosing its
! coefficients is minimization and lives one level up. Polynomials
! are one concretion, waves another; a fit holds a graph_form the way an
! operator holds coefficients - as data about shape.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module structure_graph_form

  use iso_fortran_env    , only : dp => REAL64
  use structure_graph, only : GRAPH_SIDE_VERTEX
  use structure_support, only : support

  implicit none

  private
  public :: graph_form

  type, abstract, extends(support) :: graph_form

   contains

     procedure(form_size_interface)  , deferred :: size_of
     procedure(form_values_interface), deferred :: values
     procedure(form_slopes_interface), deferred :: slopes

     procedure :: restrict

  end type graph_form

  abstract interface

     pure integer function form_size_interface(this)
       import :: graph_form
       class(graph_form), intent(in) :: this
     end function form_size_interface

     pure subroutine form_values_interface(this, x, at, phi)
       import :: graph_form, dp
       class(graph_form), intent(in) :: this
       real(dp), intent(in)  :: x(3), at(3)
       real(dp), intent(out) :: phi(:)
     end subroutine form_values_interface

     pure subroutine form_slopes_interface(this, x, at, direction, dphi)
       import :: graph_form, dp
       class(graph_form), intent(in) :: this
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

  pure subroutine restrict(this, kept)

    class(graph_form), intent(inout) :: this
    integer    , intent(in)    :: kept(:)

    this % support = support(GRAPH_SIDE_VERTEX, kept)

  end subroutine restrict

end module structure_graph_form
