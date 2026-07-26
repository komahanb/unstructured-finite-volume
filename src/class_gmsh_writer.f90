!=====================================================================!
! Write a gmsh post-processing file: the input mesh copied verbatim with
! the cell-centred solution appended as an $ElementData block, keyed by
! the original gmsh element tags. Open the result with
!
!     gmsh solution.msh
!
! and gmsh renders the field natively - no $Entities synthesis, and the
! tags are guaranteed consistent because they are the very tags the
! loader read from $Elements (grid % cell_numbers). The paraview .vtu
! path is unaffected; this is the gmsh-native sibling.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_gmsh_writer

  use iso_fortran_env , only : dp => real64
  use class_file      , only : file
  use class_string    , only : string

  implicit none

  private
  public :: gmsh_writer

  type :: gmsh_writer

     ! The writer copies the geometry and topology from this 4.1
     ! input mesh file.
     character(len=:), allocatable :: meshfile

   contains

     procedure :: write => write_solution
     procedure :: write_time_series

  end type gmsh_writer

  interface gmsh_writer
     module procedure create
  end interface gmsh_writer

contains

  !===================================================================!
  ! Constructor: remember which mesh file the fields belong to.
  !===================================================================!

  pure type(gmsh_writer) function create(meshfile) result(this)

    character(len=*), intent(in) :: meshfile

    this % meshfile = meshfile

  end function create

  !===================================================================!
  ! Copy the input mesh and append the cell field as $ElementData.
  !===================================================================!

  impure subroutine write_solution(this, filename, cell_numbers, values, label)

    class(gmsh_writer), intent(in) :: this
    character(len=*)  , intent(in) :: filename
    integer           , intent(in) :: cell_numbers(:)
    real(dp)          , intent(in) :: values(:)
    character(len=*)  , intent(in) :: label

    type(file)                :: src
    type(string), allocatable :: lines(:)
    integer                   :: unit, i

    ! Copy the input mesh verbatim; long lines demand a wide buffer.
    src = file(this % meshfile, 4096)
    call src % read_lines(lines)

    open(newunit = unit, file = filename, action = 'write', status = 'replace')

    do i = 1, size(lines)
       write(unit, '(a)') trim(lines(i) % str)
    end do

    ! Append the cell-centred field as a gmsh view.
    write(unit, '(a)') "$ElementData"
    write(unit, '(a)') "1"                 ! One string tag.
    write(unit, '(a)') '"'//label//'"'     ! The view name.
    write(unit, '(a)') "1"                 ! One real tag.
    write(unit, '(a)') "0.0"               ! The time.
    write(unit, '(a)') "3"                 ! Three integer tags.
    write(unit, '(a)') "0"                 ! The time step.
    write(unit, '(a)') "1"                 ! One component (a scalar).
    write(unit, '(i0)') size(values)       ! The number of element values.

    do i = 1, size(values)
       write(unit, '(i0,1x,es22.15)') cell_numbers(i), values(i)
    end do

    write(unit, '(a)') "$EndElementData"

    close(unit)

  end subroutine write_solution

  !===================================================================!
  ! Copy the input mesh and append several named cell fields over a
  ! time series. The values array is shaped (ncell, nfield, nstep);
  ! each field becomes a gmsh view (named by names) and each step an
  ! $ElementData block carrying its time and step index, so gmsh
  ! groups same-named blocks into one animated view. Steady export is
  ! just nstep = 1.
  !===================================================================!

  impure subroutine write_time_series(this, filename, cell_numbers, names, times, values)

    class(gmsh_writer), intent(in) :: this
    character(len=*)  , intent(in) :: filename
    integer           , intent(in) :: cell_numbers(:)
    character(len=*)  , intent(in) :: names(:)        ! (nfield)
    real(dp)          , intent(in) :: times(:)        ! (nstep)
    real(dp)          , intent(in) :: values(:,:,:)   ! (ncell, nfield, nstep)

    type(file)                :: src
    type(string), allocatable :: lines(:)
    integer                   :: unit, i, ifield, istep, ncell

    ncell = size(values, 1)

    ! Copy the input mesh verbatim; long lines demand a wide buffer.
    src = file(this % meshfile, 4096)
    call src % read_lines(lines)

    open(newunit = unit, file = filename, action = 'write', status = 'replace')

    do i = 1, size(lines)
       write(unit, '(a)') trim(lines(i) % str)
    end do

    !-----------------------------------------------------------------!
    ! One $ElementData block is written per (field, step). Blocks
    ! sharing a view name across steps form a single animated view in
    ! gmsh.
    !-----------------------------------------------------------------!

    fields: do ifield = 1, size(names)
       steps: do istep = 1, size(times)

          write(unit, '(a)') "$ElementData"
          write(unit, '(a)') "1"                          ! One string tag.
          write(unit, '(a)') '"'//trim(names(ifield))//'"'! The view name.
          write(unit, '(a)') "1"                          ! One real tag.
          write(unit, '(es22.15)') times(istep)           ! The time.
          write(unit, '(a)') "3"                           ! Three integer tags.
          write(unit, '(i0)') istep - 1                    ! The time step index.
          write(unit, '(a)') "1"                           ! One component (a scalar).
          write(unit, '(i0)') ncell                        ! The number of element values.

          do i = 1, ncell
             write(unit, '(i0,1x,es22.15)') cell_numbers(i), values(i, ifield, istep)
          end do

          write(unit, '(a)') "$EndElementData"

       end do steps
    end do fields

    close(unit)

  end subroutine write_time_series

end module class_gmsh_writer
