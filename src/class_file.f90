module class_file

  use class_string, only : string

  implicit none
  
  private
  public :: file

  !-------------------------------------------------------------------!
  ! Derived type for file
  !-------------------------------------------------------------------!
  
  type :: file

     character(:), allocatable :: filename ! filename
     type(integer) :: file_unit            ! file unit number

   contains

     ! Override
     procedure :: open
     procedure :: close
     procedure :: get_unit
     procedure :: read_line
     procedure :: read_lines
     procedure :: get_num_lines

  end type file

  !-------------------------------------------------------------------!
  ! Interface to construct a file
  !-------------------------------------------------------------------!

  interface file
     module procedure create
  end interface file

contains

  !===================================================================!
  ! Construct a file object with filename
  !===================================================================!
  
  type(file) function create(filename) result (this)

    type(character(*)), intent(in) :: filename
    logical :: ok
    integer :: i

    ! Set the file name
    allocate(this % filename, source=filename)

    ! Use an available handle for opening
    i = 99
    check_unit: do 
       i = i + 1 
       inquire(unit=i, opened=ok)
       if(ok .eqv. .false.) then
          this % file_unit = i
          exit check_unit
       end if
    end do check_unit

  end function create
  
  !===================================================================!
  ! Open the file
  !===================================================================!
  
  subroutine open(this)

    class(file), intent(in) :: this

    open(unit = this % file_unit, file = this % filename)

  end subroutine open

  !===================================================================!
  ! Close the file
  !===================================================================!

  subroutine close(this)

    class(file), intent(in) :: this

    close(unit = this % file_unit)

  end subroutine close

  !===================================================================!
  ! return the file unit number
  !===================================================================!
  
  pure type(integer) function get_unit(this)

    class(file), intent(in) :: this

    get_unit = this % file_unit

  end function get_unit

    !===================================================================!
    ! Utility function for get number of lines in mesh file
    !===================================================================!
    
    type(integer) function get_num_lines(this) result(nlines)

      class(file) , intent(in) :: this
      integer :: stat

      nlines = 0 
      call this % open()
      do
         read(this % get_unit(),*, iostat=stat)
         if (stat .ne. 0) exit
         nlines = nlines + 1
      end do
      call this % close()

    end function get_num_lines

    !=================================================================!
    ! Read all lines and return a string array
    !=================================================================!
    
    subroutine read_lines(this, lines)

      ! Arguments
      class(file)  , intent(in)            :: this
      type(string) , allocatable, intent(out) :: lines(:)

      ! Locals
      integer :: num_lines
      integer :: iline

      ! Get number of lines in file to allocate space (uses a separate
      ! handle) fix!
      num_lines = this % get_num_lines()
      allocate(lines(num_lines))

      ! Loop through each line and read and store into lines
      call this % open()
      do iline = 1, num_lines
         call this % read_line(lines(iline))
      end do
      call this % close()

    end subroutine read_lines
    
    !=================================================================!
    ! Read one line and return a string object
    !=================================================================!
!!$    
!!$    subroutine read_line(this, line)
!!$
!!$      ! Arguments
!!$      class(file)  , intent(in)  :: this
!!$      type(string) , intent(out) :: line
!!$      character(len=1024) :: buffer(:)
!!$      integer :: iostat, isize
!!$
!!$      print *, 'reading line from unit', this % get_unit()
!!$
!!$      read( this % get_unit(), &
!!$           & iostat  = iostat, &
!!$           & fmt     = '(a)' , &
!!$           & advance = 'no'  , &
!!$           & size    = isize) buffer
!!$      
!!$      print *, 'size of line', isize, 'iostat', iostat, buffer
!!$
!!$      line = string(buffer)
!!$
!!$    end subroutine read_line
!!$    
    subroutine read_line(this, linestr)

      ! Arguments
      class(file)  , intent(in)    :: this
      type(string) , intent(out)   :: linestr

      ! Locals
      character(len=:) , allocatable :: line
      integer          , parameter   :: buflen=1024
      character(len=buflen)          :: buffer
      integer                        :: last
      integer                        :: isize, ier

      isize = 0
      line = ''
      ier = 0

      ! read characters from line and append to result
      INFINITE: do
         ! read next buffer (an improvement might be to use stream I/O
         ! for files other than stdin so system line limit is not
         ! limiting)
         read(this % get_unit(),iostat=ier,fmt='(a)',advance='no',size=isize) buffer
         ! append what was read to result
         if(isize.gt.0)line=line//buffer(:isize)
         ! if hit EOR reading is complete unless backslash ends the line
         if(is_iostat_eor(ier))then
            last=len(line)
            if(last.ne.0)then
               ! if line ends in backslash it is assumed a continued line
               if(line(last:last).eq."\")then
                  ! remove backslash
                  line=line(:last-1)
                  ! continue on and read next line and append to result
                  cycle INFINITE
               endif
            endif
            ! hitting end of record is not an error for this routine
            ier=0
            ! end of reading line
            exit INFINITE
            ! end of file or error
         elseif(ier.ne.0)then
            exit INFINITE
         endif
      enddo INFINITE

      line = trim(line)
      linestr = string(line)

    end subroutine read_line

end module class_file
