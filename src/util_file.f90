!=====================================================================!
! Module that represents a file as a chain graph. Lines are the
! vertices, "next line" is the only edge, and the reader has exactly
! one operation, to traverse the edge forward:
!
!     (line 1) --> (line 2) --> (line 3) --> ...
!
! Open and close delimit the traversal, read_line takes one step,
! read_lines reads the whole chain into memory, and num_lines measures
! its length.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_file

  use util_string, only : string

  implicit none

  private

  public :: file

  !-------------------------------------------------------------------!
  ! Type that models a file: its name, its unit and the width of a
  ! line's buffer
  !-------------------------------------------------------------------!

  type :: file

     character(:), allocatable :: filename    ! the file name
     type(integer)             :: file_unit   ! the file unit number
     integer                   :: buffer_size ! the line buffer width

   contains

     ! Provided procedures
     procedure :: open
     procedure :: close
     procedure :: unit
     procedure :: read_line
     procedure :: read_lines
     procedure :: num_lines

     ! The destructor
     final :: destroy

  end type file

  !-------------------------------------------------------------------!
  ! The constructor interface for a file
  !-------------------------------------------------------------------!

  interface file
     module procedure create
  end interface file

contains

  !===================================================================!
  ! Creates a file object of the given name, with lines of the given
  ! width (100 when none is given), on the first unopened unit from 100
  !===================================================================!

  impure type(file) function create(filename, line_width) result (this)

    type(character(*)), intent(in)           :: filename
    type(integer)     , intent(in), optional :: line_width

    logical :: valid
    integer :: i

    ! The file name
    allocate(this % filename, source=filename)

    ! The line width
    if (present(line_width)) then
       this % buffer_size = line_width
    else
       this % buffer_size = 100
    end if

    ! The first unopened unit number
    i = 99
    check_unit: do
       i = i + 1
       inquire(unit=i, opened=valid)
       if(valid .eqv. .false.) then
          this % file_unit = i
          exit check_unit
       end if
    end do check_unit

  end function create

  !===================================================================!
  ! Destroys a file object, releasing its name
  !===================================================================!

  pure subroutine destroy(this)

    type(file), intent(inout) :: this

    if(allocated(this % filename)) deallocate(this % filename)

  end subroutine destroy

  !===================================================================!
  ! Opens the file; a file that does not exist stops the program
  !===================================================================!

  impure subroutine open(this)

    class(file), intent(in) :: this

    logical :: file_exists

    inquire(file=this % filename, exist=file_exists)
    if (file_exists  .eqv. .false.) then
       error stop 'file: does not exist - ' // this % filename
    end if

    open(unit = this % file_unit, file = this % filename)

  end subroutine open

  !===================================================================!
  ! Closes the file
  !===================================================================!

  impure subroutine close(this)

    class(file), intent(in) :: this

    close(unit = this % file_unit)

  end subroutine close

  !===================================================================!
  ! Returns the file unit number
  !===================================================================!

  pure type(integer) function unit(this)

    class(file), intent(in) :: this

    unit = this % file_unit

  end function unit

  !===================================================================!
  ! Returns the number of lines in the file, counted by one traversal
  !===================================================================!

  impure type(integer) function num_lines(this) result(nlines)

    class(file), intent(in) :: this

    integer :: stat

    nlines = 0
    call this % open()
    do
       read(this % unit(),*, iostat=stat)
       if (stat .ne. 0) exit
       nlines = nlines + 1
    end do
    call this % close()

  end function num_lines

  !===================================================================!
  ! Reads one line and returns it as a string, trailing blanks trimmed
  !===================================================================!

  impure subroutine read_line(this, line)

    ! Arguments
    class(file) , intent(in)  :: this
    type(string), intent(out) :: line

    ! Locals
    character(len=this % buffer_size) :: buffer

    read(this % unit(), fmt = '(a)') buffer
    line = string(trim(buffer))

  end subroutine read_line

  !===================================================================!
  ! Reads all lines and returns them as an array of strings.
  !
  ! The lines are counted first, so that the array can be allocated;
  ! that count opens and traverses the file a second time, which needs
  ! a fix:
  !
  !   num_lines:   open --> line 1 --> ... --> line n --> close
  !   read:        open --> line 1 --> ... --> line n --> close
  !===================================================================!

  impure subroutine read_lines(this, lines)

    ! Arguments
    class(file) , intent(in)               :: this
    type(string), allocatable, intent(out) :: lines(:)

    ! Locals
    integer :: num_lines
    integer :: iline

    num_lines = this % num_lines()
    allocate(lines(num_lines))

    ! The chain traversed, each line stored
    call this % open()
    do iline = 1, num_lines
       call this % read_line(lines(iline))
    end do
    call this % close()

  end subroutine read_lines

end module util_file
