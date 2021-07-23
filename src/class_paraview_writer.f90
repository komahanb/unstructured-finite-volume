!=====================================================================!
! Module that implements output capabilities for Paraview
! visualization.
!
! UnstructuredGrid ( .vtu ) is the supported format
!
! Author: Komahan Boopathy (komibuddy@gmail.com)
!=====================================================================!

module class_paraview_writer

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  use class_mesh     , only : mesh_t => mesh
  use class_string   , only : string
  implicit none

  !===================================================================!
  ! Paraview writer datatype
  !===================================================================!

  type :: paraview_writer

     ! attributes
     class(mesh_t), allocatable :: mesh

   contains

     ! type bound procedures
     procedure :: write
    
  end type paraview_writer
  
  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!

  interface paraview_writer
     module procedure construct
  end interface paraview_writer
  
contains

  !===================================================================!
  ! Constructor for paraview wrtier
  !===================================================================!
  
  pure type(paraview_writer) function construct(mesh) result (this)

    class(mesh_t), intent(in) :: mesh

    allocate(this % mesh, source = mesh)

  end function construct

  !===================================================================!
  ! Constructor for paraview wrtier
  !===================================================================!
  
  subroutine write(this, filename, phic, solution_labels)

    ! arguments
    class(paraview_writer) , intent(in)           :: this
    character(len=*)       , intent(in)           :: filename
    real(dp)               , intent(in), optional :: phic(:,:) ! (icell, ivar)
    type(string)           , optional             :: solution_labels(:)

    ! locals
    integer                       :: ierr
    integer, parameter            :: fhandle = 90

    open(unit=fhandle, file=trim(filename), iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') trim(filename)
       return
    end if

    write(fhandle, *) '<VTKFile type="UnstructuredGrid">'
    write(fhandle, *) '<UnstructuredGrid>'

    write(fhandle, *) '<Piece NumberOfPoints="', this % mesh % num_vertices, &
         & '" NumberOfCells="', this % mesh % num_cells, '">'

    write(fhandle, *) '</Piece>'
    write(fhandle, *) '</UnstructuredGrid>'
    write(fhandle, *) '</VTKFile>'

    close(unit=fhandle)
    
  end subroutine write
  
end module class_paraview_writer


