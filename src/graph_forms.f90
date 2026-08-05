!=====================================================================!
! LEVEL 1 . THE FORMS
!
! A form is a family of functions of position - a basis shape - and
! nothing more: no goal, no fitting, no solve. Evaluating a form at
! a point is calculus; choosing its coefficients is minimization and
! lives one level up; adapting which members stand active is the
! form sector's governance, above that. This module holds only the
! first: what the functions ARE.
!
!      values(x, at)            each member, evaluated at x,
!                               reckoned about the point `at`
!      slopes(x, at, n)         each member's derivative along n
!      active                   the roster the form sector writes
!                               and every fit honours
!
! Polynomials are one concretion, waves another; a fit holds a form
! the way an operator holds coefficients - as data about shape,
! owned one level below the act that uses it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_forms

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: form

  type, abstract :: form

     logical, allocatable :: active(:)

   contains

     procedure(form_size_interface)  , deferred :: size_of
     procedure(form_values_interface), deferred :: values
     procedure(form_slopes_interface), deferred :: slopes

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

end module graph_forms
