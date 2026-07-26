!=====================================================================!
! This module contains a derived type 'string' and its implemented
! procedures.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_string

  use iso_fortran_env, only : dp => REAL64

  implicit none
  
  private
  public :: string

  !-------------------------------------------------------------------!
  ! A derived type for a string.
  !-------------------------------------------------------------------!
  
  type :: string

     character(:), allocatable :: str ! The character array.
     type(integer) :: count ! The length.

   contains

     ! Overridden procedures.
     procedure :: print
     procedure :: equals
     procedure :: tokenize
     procedure :: asinteger
     procedure :: asreal

     ! The destructor.
     final :: destroy

  end type string

  !-------------------------------------------------------------------!
  ! The constructor interface for a string.
  !-------------------------------------------------------------------!

  interface string
     module procedure create
  end interface string

contains
  
  !===================================================================!
  ! Tokenize the string object and return an array of tokens.
  !===================================================================!
  
  pure subroutine tokenize(this, delimiter, num_tokens, tokens)

    ! Arguments.
    class(string)    , intent(in)  :: this
    character(len=*) , intent(in)  :: delimiter
    integer          , intent(out) :: num_tokens    
    type(string)     , intent(out), allocatable, optional :: tokens(:)
    
    ! Locals.
    integer , allocatable :: tidx(:,:)
    integer :: sidx, eidx   
    integer :: token_idx, token_ctr
    integer :: i

    num_tokens = 0
    if (len(delimiter) .eq. 0) return ! The delimiter does not match.
    if (index(this % str, delimiter) .eq. 0) return ! The delimiter does not match.

    ! The table holds the lower and upper index of each token.
    allocate(tidx(2, this % count))
    
    ! Initialize the indices and the token counter.
    sidx      = 1
    eidx      = len(this % str)
    token_ctr = 0 
    parse: do while (len(this % str(sidx:eidx)) .ge. 0)

       ! Find the next index of the delimiter.
       token_idx = index(this % str(sidx:eidx), delimiter)

       if (token_idx .ne. 0) then

          token_ctr = token_ctr + 1
          tidx(:,token_ctr) = [sidx, token_idx - 1 + sidx]

          ! A match was found; record it and advance the start index.
          sidx = sidx + token_idx

       else

          ! Check whether this is the last substring.
          if (token_ctr .ge. 1) then

             ! The trailing substring is itself a token.
             token_ctr = token_ctr + 1
                       
             token_idx = 1

             tidx(:,token_ctr) = [sidx, eidx]

          end if

          exit parse

       end if

    end do parse

    ! Set the return arguments.
    num_tokens = token_ctr
    if (present(tokens)) then
       if (allocated(tokens)) deallocate(tokens)
       allocate(tokens(num_tokens))
       do concurrent (i=1:num_tokens)
          tokens(i) = string(this % str(tidx(1,i):tidx(2,i)))
       end do
    end if
    
  end subroutine tokenize

  !===================================================================!
  ! Construct a string object from the supplied literal, find its
  ! length, and initialize its hashcode as zero.
  !===================================================================!

  pure elemental type(string) function create(str) result (this)

    type(character(*)), intent(in) :: str
    
    allocate(this % str, source=str) ! Source copies the value; mold does not.
    
    this % count = len(str)    

  end function create
  
  !===================================================================!
  ! The destructor for a string object.
  !===================================================================!
  
  pure subroutine destroy(this)

    type(string), intent(inout) :: this

    if(allocated(this % str)) deallocate(this % str)

  end subroutine destroy

  !===================================================================!
  ! Overridden string equality logic, based on a comparison of entries.
  !===================================================================!
  
  pure elemental type(logical) function equals(this, element)
    
    class(string), intent(in) :: this
    class(string), intent(in) :: element 

    ! Two string objects are equal when their values are equal.
    equals = (element % str .eq. this % str)

  end function equals

  !===================================================================!
  ! Print the string representation of the object.
  !===================================================================!
  
  impure elemental subroutine print(this, fmt)
    
    class(string)   , intent(in) :: this
    character(len=*), intent(in), optional :: fmt

    if (present(fmt)) then
       
       if (allocated(this % str)) then
          write(*,fmt) this % str
       else
          print *, "string : ", "NULL"
       end if
       
       return
       
    end if
    
    if (allocated(this % str)) then
       print *, "string : ", this % str
    else
       print *, "string : ", "NULL"
    end if
    
  end subroutine print

  !===================================================================!  
  ! Return the integer evaluation of the string.
  !===================================================================!
  
  pure elemental type(integer)  function asinteger(this)

    class(string), intent(in) :: this

    read (this % str,*) asinteger

  end function asinteger

  !===================================================================!  
  ! Return the real number parsed from the string.
  !===================================================================!
  
  pure elemental type(real(dp)) function asreal(this)

    class(string), intent(in) :: this

    read (this % str,*) asreal
    
  end function asreal
  
end module class_string
