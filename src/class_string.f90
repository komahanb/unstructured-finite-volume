!=====================================================================!
! Contains a derived type 'string' and implemented procedures
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_string

  use iso_fortran_env, only : dp => REAL64, error_unit

  implicit none
  
  private
  public :: string

  !-------------------------------------------------------------------!
  ! Derived type for string
  !-------------------------------------------------------------------!
  
  type :: string

     character(:), allocatable :: str ! character array
     type(integer) :: count ! length

   contains

     ! Override
     procedure :: print
     procedure :: equals
     procedure :: tokenize
     procedure :: asinteger
     procedure :: asreal

     ! Destructor
     final :: destroy

  end type string

  !-------------------------------------------------------------------!
  ! Interface to construct a string
  !-------------------------------------------------------------------!

  interface string
     module procedure create
  end interface string

contains
  
  !===================================================================!
  ! Tokenize the string object and return an array of tokens 
  !===================================================================!
  
  pure subroutine tokenize(this, delimiter, tokens, num_tokens)

    ! Arguments
    class(string)    , intent(in)               :: this
    character(len=*) , intent(in)               :: delimiter
    type(string)     , intent(out), allocatable :: tokens(:)
    integer          , intent(out)              :: num_tokens
    
    ! Locals
    integer , allocatable :: tidx(:,:)
    integer :: sidx, eidx   
    integer :: token_idx, token_ctr
    integer :: i

    num_tokens = 0
    if (len(delimiter) .eq. 0) return ! doesnt match
    if (index(this % str, delimiter) .eq. 0) return ! doesnt match

    ! Lower and upper index of tokens
    allocate(tidx(2, this % count))
    
    ! Initialize
    sidx      = 1
    eidx      = len(this % str)
    token_ctr = 0 
    parse: do while (len(this % str(sidx:eidx)) .gt. 0)

       ! Get the -th index of delimiter
       token_idx = index(this % str(sidx:eidx), delimiter)

       if (token_idx .ne. 0) then

          token_ctr = token_ctr + 1
          tidx(:,token_ctr) = [sidx, token_idx - 1 + sidx]

          ! We found the match record the index
          sidx = sidx + token_idx

       else

          ! Check if its the last substring
          if (token_ctr .gt. 1) then

             ! Yes, this is a token
             token_ctr = token_ctr + 1

             token_idx = 1

             tidx(:,token_ctr) = [sidx, eidx]

          end if

          exit parse

       end if

    end do parse

    ! Set the return arguments
    num_tokens = token_ctr
    allocate(tokens(num_tokens))
    do i = 1, num_tokens
       tokens(i) = string(this % str(tidx(1,i):tidx(2,i)))
    end do

  end subroutine tokenize

  !===================================================================!
  ! Construct a string object from the supplied literal, find its
  ! length, initialize its hashcode as zero.
  !===================================================================!

  pure type(string) function create(str) result (this)

    type(character(*)), intent(in) :: str
    
    allocate(this % str, source=str) ! source copies, mold does not
    
    this % count = len(str)    

  end function create
  
  !===================================================================!
  ! Destructor for string object
  !===================================================================!
  
  pure subroutine destroy(this)

    type(string), intent(inout) :: this

    if(allocated(this % str)) deallocate(this % str)

  end subroutine destroy

  !===================================================================!
  ! Overridden string equality logic. Based on comparison of entries
  !===================================================================!
  
  pure elemental type(logical) function equals(this, element)
    
    class(string), intent(in) :: this
    class(string), intent(in) :: element 

    ! string objects are equal if their values are equal
    equals = (element % str .eq. this % str)

  end function equals

  !===================================================================!
  ! Returns the string representation of the object
  !===================================================================!
  
  impure elemental subroutine print(this)
    
    class(string), intent(in) :: this
    
    print *, "string : ", this % str
    
  end subroutine print

  !===================================================================!  
  ! Return the integer evaluation of string
  !===================================================================!
  
  pure elemental type(integer)  function asinteger(this)

    class(string), intent(in) :: this

    read (this % str,*) asinteger

  end function asinteger

  !===================================================================!  
  ! Get the real number from string
  !===================================================================!
  
  pure elemental type(real(dp)) function asreal(this)

    class(string), intent(in) :: this

    read (this % str,*) asreal
    
  end function asreal
  
end module class_string
