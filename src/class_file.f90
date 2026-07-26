!=====================================================================!
! A file, seen the way we see everything here: a chain. Lines are the
! vertices, "next line" is the only edge, and the reader can do just
! one thing - take the edge forward:
!
!     (line 1)──▶(line 2)──▶(line 3)──▶ ...
!
! Open and close bracket the walk, read_line takes one step,
! read_lines walks the whole chain into memory, and get_num_lines
! measures its length.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_file

  use class_string, only : string

  implicit none
  
  private
  public :: file

  !-------------------------------------------------------------------!
  ! A derived type for a file.
  !-------------------------------------------------------------------!
  
  type :: file

     character(:), allocatable :: filename ! The file name.
     type(integer) :: file_unit            ! The file unit number.
     integer       :: buffer_size          ! The line buffer width.

   contains

     ! Overridden procedures.
     procedure :: open
     procedure :: close
     procedure :: get_unit
     procedure :: read_line
     procedure :: read_lines
     procedure :: get_num_lines

     ! The destructor.
     final :: destroy

  end type file

  !-------------------------------------------------------------------!
  ! The constructor interface for a file.
  !-------------------------------------------------------------------!

  interface file
     module procedure create
  end interface file

contains

  !===================================================================!
  ! Construct a file object with the given file name.
  !===================================================================!
  
  impure type(file) function create(filename, line_width) result (this)

    type(character(*)), intent(in)           :: filename
    type(integer)     , intent(in), optional :: line_width
    logical :: ok
    integer :: i

    ! Set the file name.
    allocate(this % filename, source=filename)

    ! Set the line width.
    if (present(line_width)) then
       this % buffer_size = line_width
    else
       this % buffer_size = 100
    end if

    ! Use an available handle for opening.
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
  ! The destructor for a file object.
  !===================================================================!
  
  pure subroutine destroy(this)

    type(file), intent(inout) :: this

    if(allocated(this % filename)) deallocate(this % filename)

  end subroutine destroy

  !===================================================================!
  ! Open the file.
  !===================================================================!

  impure subroutine open(this)

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
  ! Close the file.
  !===================================================================!

  impure subroutine close(this)

    class(file), intent(in) :: this

    close(unit = this % file_unit)

  end subroutine close

  !===================================================================!
  ! Return the file unit number.
  !===================================================================!

  pure type(integer) function get_unit(this)

    class(file), intent(in) :: this

    get_unit = this % file_unit

  end function get_unit

  !===================================================================!
  ! A utility function that counts the number of lines in the mesh file.
  !===================================================================!

  impure type(integer) function get_num_lines(this) result(nlines)

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
  ! Read one line and return a string object.
  !=================================================================!

  impure subroutine read_line(this, line)

    ! Arguments.
    class(file)  , intent(in)  :: this
    type(string) , intent(out) :: line

    ! Locals.
    character(len=this % buffer_size) :: buffer

    read(this % get_unit(), fmt = '(a)') buffer
    line = string(trim(buffer))

  end subroutine read_line

  !=================================================================!
  ! Read all lines and return a string array.
  !=================================================================!

  impure subroutine read_lines(this, lines)

    ! Arguments.
    class(file)  , intent(in)            :: this
    type(string) , allocatable, intent(out) :: lines(:)

    ! Locals.
    integer :: num_lines
    integer :: iline

    !-----------------------------------------------------------------!
    ! Count the lines in the file so that space can be allocated.
    ! This measurement walks the chain on a separate handle, which
    ! needs a fix.
    !-----------------------------------------------------------------!

    num_lines = this % get_num_lines()
    allocate(lines(num_lines))

    ! Walk the chain and store each line into the lines array.
    call this % open()
    do iline = 1, num_lines
       call this % read_line(lines(iline))
    end do
    call this % close()

  end subroutine read_lines
   
end module class_file
