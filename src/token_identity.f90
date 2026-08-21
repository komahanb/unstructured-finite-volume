!=====================================================================!
! THE IDENTITY . INFRASTRUCTURE BENEATH THE TOWER
!
! Not a level: a service every identified object draws on. Carriers
! have identity; relations have identity; the relational graph will
! be the third to sign. What they share is not mathematics but the
! identity rules, and those live here once:
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
! local to one image BY DESIGN - matches is the one comparison,
! and the only one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module token_identity

  implicit none

  private
  public :: token, next_token
  public :: index_of

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
  ! next_token hands out the next stamp of this image: fresh,
  ! unrepeatable, contents unchoosable.
  !===================================================================!

  type(token) function next_token()

    last_serial          = last_serial + 1
    next_token % serial  = last_serial
    next_token % image   = this_image()

  end function next_token

  !===================================================================!
  ! The index of key in keys, or zero if no element matches. The one
  ! written form of the identity scan: every map that stores tokens
  ! finds them here, so identity comparison over a collection has
  ! exactly one body.
  !===================================================================!

  pure integer function index_of(keys, key) result(at)

    type(token), intent(in) :: keys(:)
    type(token), intent(in) :: key

    integer :: k

    at = 0
    do k = 1, size(keys)
       if (keys(k) % matches(key)) then
          at = k
          return
       end if
    end do

  end function index_of

  !===================================================================!
  ! The token's three queries.
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
