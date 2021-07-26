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
  use iso_fortran_env, only : dp => real64, int32
  use class_mesh     , only : mesh_t => mesh
  use class_string   , only : string
  implicit none

  !===================================================================!
  ! Paraview writer datatype
  !===================================================================!

  type :: paraview_writer

     ! attributes
     class(mesh_t), allocatable :: mesh

     ! support binary and ascii
     ! deduce cell types from num cell vertices
     ! write cell and point data

   contains

     ! type bound procedures
     procedure :: write
    
  end type paraview_writer


  type :: linear_cell_type

     integer(kind=int32) :: VTK_EMPTY_CELL        = 0
     integer(kind=int32) :: VTK_VERTEX            = 1
     integer(kind=int32) :: VTK_POLY_VERTEX       = 2
     integer(kind=int32) :: VTK_LINE              = 3
     integer(kind=int32) :: VTK_POLY_LINE         = 4
     integer(kind=int32) :: VTK_TRIANGLE          = 5
     integer(kind=int32) :: VTK_TRIANGLE_STRIP    = 6
     integer(kind=int32) :: VTK_POLYGON           = 7
     integer(kind=int32) :: VTK_PIXEL             = 8
     integer(kind=int32) :: VTK_QUAD              = 9
     integer(kind=int32) :: VTK_TETRA             = 10
     integer(kind=int32) :: VTK_VOXEL             = 11
     integer(kind=int32) :: VTK_HEXAHEDRON        = 12
     integer(kind=int32) :: VTK_WEDGE             = 13
     integer(kind=int32) :: VTK_PYRAMID           = 14
     integer(kind=int32) :: VTK_PENTAGONAL_PRISM  = 15
     integer(kind=int32) :: VTK_HEXAGONAL_PRISM   = 16

  end type linear_cell_type
  
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

    
    open(unit=fhandle, file=trim(filename), iostat= ierr, action = 'write', form = 'formatted')
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') trim(filename)
       return
    end if

    write(fhandle, '(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(fhandle, '(a)') '<UnstructuredGrid>'
    write(fhandle, '(a,i0,a,i0,a)') '<Piece NumberOfPoints="', this % mesh % num_vertices, &
         & '" NumberOfCells="', this % mesh % num_cells, '">'
    
    write_points: block

      integer :: ivertex, jdim

      write(fhandle, '(a)') '<Points>'
      write(fhandle, '(a)') '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'

      do ivertex = 1, this % mesh % num_vertices
         
         write(fhandle, *) (this % mesh % vertices(jdim, ivertex), jdim = 1, 3)
         
      end do

      write(fhandle, '(a)') '</DataArray>'
      write(fhandle, '(a)') '</Points>'

    end block write_points

    write_cells: block

      integer :: icell, jvertex
      integer, allocatable :: cell_types(:)

      integer :: offset
      
      write(fhandle, '(a)') '<Cells>'

      !---------------------------------------------------------------!
      ! write cell vertex connectivities
      !---------------------------------------------------------------!
      
      write(fhandle, '(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
      
      do icell = 1, this % mesh % num_cells

         ! correct for paraview's 0 based numbering
         write(fhandle, *) (this % mesh % cell_vertices(jvertex, icell) - 1, &
              & jvertex = 1, this % mesh % num_cell_vertices(icell))
         
      end do
      
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! write cell to vertex connectivity offsets
      !---------------------------------------------------------------!
      
 
      write(fhandle, '(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'

      offset = 0
      
      do icell = 1, this % mesh % num_cells

         write(fhandle, '(i0)') offset + this % mesh % num_cell_vertices(icell)
         
         offset = offset + this % mesh % num_cell_vertices(icell)

      end do
      
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! write cell types
      !---------------------------------------------------------------!
      
      write(fhandle, *) '<DataArray type="UInt8" Name="types" format="ascii">'
      
      do icell = 1, this % mesh % num_cells
         write(fhandle, '(i0)') 7 ! polyhedrals
      end do
      write(fhandle, '(a)') '</DataArray>'
      
      write(fhandle, '(a)') '</Cells>'

      write(fhandle, '(a)') '<PointData></PointData>'
      write(fhandle, '(a)') '<CellData></CellData>'
  
      if (allocated(cell_types))   deallocate(cell_types)
      
    end block write_cells

    ! close the opened tags
    write(fhandle, '(a)') '</Piece>'
    write(fhandle, '(a)') '</UnstructuredGrid>'
    write(fhandle, '(a)') '</VTKFile>'

    close(unit=fhandle)
    
  end subroutine write
  
end module class_paraview_writer

