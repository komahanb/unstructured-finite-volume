!=====================================================================!
! THE IDENTITY . INFRASTRUCTURE BENEATH THE TOWER
!
! Not a level: a service every identified citizen draws on. Carriers
! have identity; relations have identity; the relational graph will
! be the third to sign. What they share is not mathematics but LAW,
! and the law lives here once:
!
!      minting    fresh, unrepeatable, contents unchoosable
!      copying    whole-object only - and a copy IS the identity
!      matching   the one comparison; serial zero matches nothing
!
! The token is OPAQUE. Its parts are private, so no caller can
! compose one with chosen contents - the only ways a token comes to
! exist are minting and copying, which is precisely the difference
! between declaring a domain and being one. Today a token is an
! (image, serial) pair, so two coarray images can never mint the
! same stamp; nothing outside this module can read the parts, so
! the representation stays free to grow distributed-safe further.
!
! serial_number() is a read-only diagnostic for messages and tests,
! local to one image BY DESIGN - matches is the law, and the only
! one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module token_identity

  implicit none

  private
  public :: token, mint_token

  !===================================================================!
  ! The stamp roll of this image. Serial zero is reserved for the
  ! undeclared: a default-initialized token is no identity at all.
  !===================================================================!

  integer, save :: last_serial = 0

  type :: token

     integer, private :: image  = 0
     integer, private :: serial = 0

   contains

     procedure :: matches
     procedure :: declared
     procedure :: serial_number

  end type token

contains

  !===================================================================!
  ! mint_token hands out the next stamp of this image: fresh,
  ! unrepeatable, contents unchoosable.
  !===================================================================!

  type(token) function mint_token()

    last_serial          = last_serial + 1
    mint_token % serial  = last_serial
    mint_token % image   = this_image()

  end function mint_token

  !===================================================================!
  ! The token's three answers.
  !===================================================================!

  pure logical function matches(this, other)

    class(token), intent(in) :: this
    type(token) , intent(in) :: other

    matches = (this % serial /= 0)              .and. &
         &    (this % serial == other % serial) .and. &
         &    (this % image  == other % image)

  end function matches

  pure logical function declared(this)

    class(token), intent(in) :: this

    declared = this % serial /= 0

  end function declared

  pure integer function serial_number(this)

    class(token), intent(in) :: this

    serial_number = this % serial

  end function serial_number

end module token_identity
