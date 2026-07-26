!=====================================================================!
! This module implements output capabilities for Paraview
! visualization.
!
! UnstructuredGrid (.vtu) is the supported format.
!
! Author: Komahan Boopathy (komibuddy@gmail.com)
!=====================================================================!

module class_paraview_writer

  ! Import the dependencies.
  use iso_fortran_env, only : dp => real64, int32
  use class_mesh     , only : mesh_t => mesh
  use class_string   , only : string

  implicit none

  !===================================================================!
  ! This datatype enumerates the Paraview cell types.
  !===================================================================!

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
     integer(kind=int32) :: VTK_POLYHEDRON        = 42

   contains

     procedure :: get_element_type

  end type linear_cell_type

  !===================================================================!
  ! This datatype writes a mesh and its cell fields to Paraview.
  !===================================================================!

  type :: paraview_writer

     ! These are the attributes.
     class(mesh_t), allocatable :: mesh

     type(linear_cell_type)     :: cell_type

     ! Support binary and ascii forms.
     ! Write cell and point data.

   contains

     ! These are the type-bound procedures.
     procedure :: write

  end type paraview_writer

  !===================================================================!
  ! This interface admits multiple constructors.
  !===================================================================!

  interface paraview_writer
     module procedure construct
  end interface paraview_writer

contains

  !===================================================================!
  ! This function maps gmsh element numbers to paraview cell types.
  !===================================================================!

  pure elemental type(integer) function get_element_type(this, gmsh_type) &
       & result (paraview_type)

    class(linear_cell_type), intent(in) :: this
    integer                , intent(in) :: gmsh_type

    select case (gmsh_type)
    case (1) ! A 2-node line.
       paraview_type = this % VTK_LINE
    case (2) ! A 3-node triangle.
       paraview_type = this % VTK_TRIANGLE
    case (3) ! A 4-node quadrangle.
       paraview_type = this % VTK_QUAD
    case (4) ! A 4-node tetrahedron.
       paraview_type = this % VTK_TETRA
    case (5) ! An 8-node hexahedron.
       paraview_type = this % VTK_HEXAHEDRON
    case (6) ! A 6-node prism.
       paraview_type = this % VTK_WEDGE
    case (7) ! A 5-node prism (a pyramid).
       paraview_type = this % VTK_PYRAMID
    case (-1) ! An agglomerated polygon - our own convention; gmsh has no such element.
       paraview_type = this % VTK_POLYGON
    case default
       paraview_type = this % VTK_POLYHEDRON
    end select

  end function get_element_type

  !===================================================================!
  ! This is the constructor for the paraview writer.
  !===================================================================!

  pure type(paraview_writer) function construct(mesh) result (this)

    class(mesh_t), intent(in) :: mesh

    allocate(this % mesh, source = mesh)

  end function construct

  !===================================================================!
  ! Write the mesh and the optional cell fields to a .vtu file in the
  ! UnstructuredGrid format.
  !===================================================================!

  impure subroutine write(this, filename, phic, solution_labels)

    ! These are the arguments.
    class(paraview_writer) , intent(in)           :: this
    character(len=*)       , intent(in)           :: filename
    real(dp)               , intent(in), optional :: phic(:,:) ! (icell, ivar)
    type(string)           , optional             :: solution_labels(:)
    type(linear_cell_type) :: paraview_type

    ! These are the locals.
    integer                       :: ierr
    integer, parameter            :: fhandle = 90
    integer                       :: iresult

    ! Open the output file for formatted writing.
    open(unit=fhandle, file=trim(filename), iostat= ierr, action = 'write', form = 'formatted')
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') trim(filename)
       return
    end if

    !-----------------------------------------------------------------!
    ! Write the basic header information.
    !-----------------------------------------------------------------!

    write(fhandle, '(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(fhandle, '(a)') '<UnstructuredGrid>'
    write(fhandle, '(a,i0,a,i0,a)') '<Piece NumberOfPoints="', this % mesh % num_points, &
         & '" NumberOfCells="', this % mesh % num_cells, '">'

    !-----------------------------------------------------------------!
    ! Write the vertices.
    !-----------------------------------------------------------------!

    write_points: block

      integer :: ivertex, jdim

      write(fhandle, '(a)') '<Points>'
      write(fhandle, '(a)') '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do ivertex = 1, this % mesh % num_points
         write(fhandle, *) (this % mesh % coordinates(jdim, ivertex), jdim = 1, 3)
      end do
      write(fhandle, '(a)') '</DataArray>'
      write(fhandle, '(a)') '</Points>'

    end block write_points

    write_cells: block

      integer :: icell, jvertex

      integer :: offset

      write(fhandle, '(a)') '<Cells>'

      !---------------------------------------------------------------!
      ! Write the cell-to-vertex connectivities.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
      do icell = 1, this % mesh % num_cells
         ! Correct for paraview's 0-based numbering.
         write(fhandle, *) (this % mesh % cell_vertices(jvertex, icell) - 1, &
              & jvertex = 1, this % mesh % num_cell_vertices(icell))
      end do
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! Write the cell-to-vertex connectivity offsets.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
      offset = 0
      do icell = 1, this % mesh % num_cells
         write(fhandle, '(i0)') offset + this % mesh % num_cell_vertices(icell)
         offset = offset + this % mesh % num_cell_vertices(icell)
      end do
      write(fhandle, '(a)') '</DataArray>'

      !---------------------------------------------------------------!
      ! Write the cell types.
      !---------------------------------------------------------------!

      write(fhandle, *) '<DataArray type="UInt8" Name="types" format="ascii">'
      do icell = 1, this % mesh % num_cells
         write(fhandle, '(i0)')  paraview_type % get_element_type (this % mesh % cell_types(icell))
      end do
      write(fhandle, '(a)') '</DataArray>'

      write(fhandle, '(a)') '</Cells>'

      write(fhandle, '(a)') '<PointData></PointData>'

      !---------------------------------------------------------------!
      ! Write the cell data.
      !---------------------------------------------------------------!

      write(fhandle, '(a)') '<CellData>'

      ! Export the cell volumes.
      write(fhandle, *) '<DataArray type="Float32" Name="volume" format="ascii">'
      do icell = 1, this % mesh % num_cells
         write(fhandle, *)  this % mesh % cell_volumes(icell)
      end do
      write(fhandle, '(a)') '</DataArray>'

      if (present(solution_labels) .and. present(phic)) then
         do iresult = 1, size(solution_labels)
            write(fhandle, '(a,a,a)') '<DataArray type="Float32" Name="', solution_labels(iresult) % str, '" format="ascii">'
            do icell = 1, this % mesh % num_cells
               write(fhandle, *) phic(icell,iresult)
            end do
            write(fhandle, '(a)') '</DataArray>'
         end do
      end if

      write(fhandle, '(a)') '</CellData>'

    end block write_cells

    ! Close the opened tags.
    write(fhandle, '(a)') '</Piece>'
    write(fhandle, '(a)') '</UnstructuredGrid>'
    write(fhandle, '(a)') '</VTKFile>'

    close(unit=fhandle)

  end subroutine write

end module class_paraview_writer
