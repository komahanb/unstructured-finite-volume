!=====================================================================!
! Module that states the identity shared by every level of the tower.
!
! It is not a level but a service every identified object uses.
! Carriers have identity, relations have identity, and the relational
! graph will be the third to use it. What they share is not
! mathematics but the rules of identity, and these are defined here
! once:
!
!      creation   new, unrepeatable, contents not selectable
!      copying    whole-object only, and a copy IS the identity
!      matching   the one comparison; serial zero matches nothing
!
! The token is opaque. Its parts are private, so no caller can compose
! one with chosen contents: a token comes to exist only by creation or
! by copying, which is the difference between declaring a domain and
! being one. A token is now an (image, serial) pair, so two coarray
! images can never create the same token:
!
!   image 1:  (1, 1)  (1, 2)  (1, 3)  ...
!   image 2:  (2, 1)  (2, 2)  ...           never equal to image 1's
!   (0, 0)    the undeclared token, which matches nothing
!
! Nothing outside this module can read the parts, so the
! representation can be extended for distributed use.
!
! serial_number() is a read-only diagnostic for messages and tests,
! local to one image by design; matches is the one comparison.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module token_identity

  implicit none

  private

  public :: token, next_token

  ! The serial counter of this image. Serial zero is reserved for the
  ! undeclared: a default-initialized token is no identity at all.
  integer, save :: last_serial = 0

  !-------------------------------------------------------------------!
  ! Type that models an identity: an (image, serial) pair, opaque to
  ! every caller
  !-------------------------------------------------------------------!

  type :: token

     integer, private :: image  = 0
     integer, private :: serial = 0

   contains

     ! Queries
     procedure :: matches
     procedure :: declared
     procedure :: serial_number

  end type token

contains

  !===================================================================!
  ! Returns the next token of this image: new, unrepeatable, contents
  ! not selectable.
  !
  ! The increment and the read of the counter are one atomic capture
  ! under OpenMP, so two threads never receive one serial. Without
  ! -fopenmp the sentinel lines are comments and the statements are the
  ! serial ones.
  !===================================================================!

  type(token) function next_token()

    integer :: serial

    !$omp atomic capture
    last_serial = last_serial + 1
    serial      = last_serial
    !$omp end atomic

    next_token % serial  = serial
    next_token % image   = this_image()

  end function next_token

  !===================================================================!
  ! Returns whether two tokens are the same identity: both declared,
  ! of the same image and the same serial
  !===================================================================!

  pure logical function matches(this, other)

    class(token), intent(in) :: this
    type(token) , intent(in) :: other

    matches = (this % serial /= 0)              .and. &
         &    (this % serial == other % serial) .and. &
         &    (this % image  == other % image)

  end function matches

  !===================================================================!
  ! Returns whether the token was created, its serial not zero
  !===================================================================!

  pure logical function declared(this)

    class(token), intent(in) :: this

    declared = this % serial /= 0

  end function declared

  !===================================================================!
  ! Returns the token's serial, a diagnostic local to its image
  !===================================================================!

  pure integer function serial_number(this)

    class(token), intent(in) :: this

    serial_number = this % serial

  end function serial_number

end module token_identity
