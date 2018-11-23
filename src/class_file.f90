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
     integer       :: buffer_size

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
  
  type(file) function create(filename, line_width) result (this)

    type(character(*)), intent(in)           :: filename
    type(integer)     , intent(in), optional :: line_width
    logical :: ok
    integer :: i

    ! Set the file name
    allocate(this % filename, source=filename)

    ! Line width
    if (present(line_width)) then
       this % buffer_size = line_width
    else
       this % buffer_size = 100
    end if

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
    logical :: file_exists

    inquire(file=this % filename, exist=file_exists)
    if (file_exists  .eqv. .false.) then
       print *, 'file does not exist ', this % filename
       error stop
    end if
    
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
  ! Read one line and return a string object
  !=================================================================!

  subroutine read_line(this, line)

    ! Arguments
    class(file)  , intent(in)  :: this
    type(string) , intent(out) :: line

    ! Locals
    character(len=this % buffer_size) :: buffer

    read(this % get_unit(), fmt = '(a)') buffer
    line = string(trim(buffer))

  end subroutine read_line

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
   
end module class_file
