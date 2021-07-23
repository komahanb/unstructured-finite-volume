module class_tecplot_writer

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  use class_mesh, only : mesh_t => mesh
  use class_assembler, only : assembler_t => assembler

  implicit none
  
  include 'tecio.f90'

  private

contains
  

  !===================================================================!
  ! Write solution to file
  !===================================================================!
  
  subroutine write_solution(assembler, mesh, filename, phic)

    class(assembler_t), intent(in)  :: assembler
    class(mesh_t)     , intent(in)  :: mesh
    character(len=*), intent(in)  :: filename
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
    real(dp)        , intent(in)  :: phic(:)
    integer                       ::  i, ierr
    real(dp)        , allocatable :: phiv(:)

    ! Compute vertex values by interpolating cell center values
    call assembler % evaluate_vertex_flux(phiv, phic)
       
!!$    open_file: block
!!$
!!$      integer(4) :: tecflag, file_format, file_type, debug, is_double
!!$
!!$      file_type   =  0 ! all, grid, mesh
!!$      file_format = 0
!!$      debug       = 0
!!$      is_double   = 1
!!$
!!$      tecflag = tecini142("simple dataset"//char(0), &
!!$           & "X Y T"//char(0), &
!!$           & "solution.plt"//char(0), &
!!$           & "."//char(0), &
!!$           & file_format, &
!!$           & file_type, &
!!$           & debug, &
!!$           & is_double)
!!$      if (tecflag .ne. 0) error stop "TECINI142 failed"
!!$
!!$    end block open_file
!!$
!!$    ! a. Cell-based finite element zones: FELINESEG, FETRIANGLE, FEQUADRILATERAL, FETETRAHEDRON, and FEBRICK.
!!$    ! b. Face-based finite element zones: FEPOLYGON and FEPOLYHEDRON.
!!$    add_zone: block
!!$      
!!$      character(len=128)           :: ZoneTitle      
!!$      integer(kind=4)              :: ZoneType
!!$      integer(kind=4)              :: NumPts, NumElements, NumFaces
!!$      integer(kind=4)              :: ICellMax, JCellMax, KCellMax
!!$      integer(kind=4)              :: StrandID
!!$      integer(kind=4)              :: ParentZone
!!$      integer(kind=4)              :: IsBlock
!!$      integer(kind=4)              :: NumFaceConnections = 0
!!$      integer(kind=4)              :: FaceNeighborMode = 0
!!$      integer(kind=4)              :: TotalNumFaceNodes = 0
!!$      integer(kind=4)              :: NumConnectedBoundaryFaces = 0
!!$      integer(kind=4)              :: TotalNumBoundaryConnections = 0
!!$      integer(kind=4)              :: ShareConnectivityFromZone = 0      
!!$      integer(kind=4), allocatable :: PassiveVarList(:)
!!$      integer(kind=4), allocatable :: ValueLocation(:) 
!!$      integer(kind=4), allocatable :: ShareVarFromZone(:)      
!!$      integer(kind=4)              :: status
!!$      real(kind=dp)                :: SolutionTime
!!$      
!!$      allocate(PassiveVarList(this % num_variables))
!!$      PassiveVarList = 0 ! all variables are active
!!$      
!!$      allocate(ValueLocation(this % num_variables))
!!$      ValueLocation = 1 ! 0 = cell, 1 = node (point format supports only nodal values)
!!$
!!$      allocate(ShareVarFromZone(this % num_variables))
!!$      ShareVarFromZone = 0 ! not shared     
!!$      
!!$      ZoneTitle     = "Temperature Grid" // char(0)
!!$      ZoneType      = 3 ! will use node duplication trick to accomodate FETRIANGLE
!!$      NumPts        = mesh % num_vertices
!!$      NumElements   = mesh % num_cells
!!$      NumFaces      = 0 ! not used except for polyhedral types
!!$      ICellMax      = 0
!!$      JCellMax      = 0
!!$      KCellMax      = 0
!!$      SolutionTime  = 0.0_dp
!!$      StrandID      = 0
!!$      ParentZone    = 0
!!$      IsBlock       = 1
!!$
!!$      status = TECZNE142(ZoneTitle, ZoneType, &
!!$           & NumPts, NumElements, NumFaces, ICellMax, JCellMax, KCellMax, &
!!$           & SolutionTime, StrandID, &
!!$           & ParentZone, IsBlock, &
!!$           & NumFaceConnections, FaceNeighborMode, TotalNumFaceNodes, &
!!$           & NumConnectedBoundaryFaces, TotalNumBoundaryConnections, PassiveVarList, &
!!$           & ValueLocation, &
!!$           & ShareVarFromZone, ShareConnectivityFromZone)
!!$      
!!$      if (status .ne. 0) error stop "TECZNE142 failed"
!!$
!!$      ! if (allocated(ShareVarFromZone)) deallocate(ShareVarFromZone)
!!$
!!$      ! if (allocated(ValueLocation)) deallocate(ValueLocation)
!!$
!!$      ! if (allocated(ShareVarFromZone)) deallocate(ShareVarFromZone)
!!$       
!!$      
!!$    end block add_zone
!!$
!!$
!!$    add_data: block
!!$
!!$      integer(kind=4) :: IsDouble = 0
!!$      integer(kind=4) :: status
!!$
!!$      ! Write vertices
!!$      do i = 1, 3
!!$         status = TECDAT142(mesh % num_vertices, real(mesh % vertices(i,:)), IsDouble)
!!$         if (status .ne. 0) error stop "TECDAT142 failed for grid"
!!$      end do
!!$
!!$      status = TECDAT142(mesh % num_vertices, real(phiv), IsDouble)
!!$      if (status .ne. 0) error stop "TECDAT142 failed for solution"
!!$
!!$      ! Write the connectivities of all cells
!!$      status = TECNOD142(mesh % cell_vertices)
!!$      if (status .ne. 0) error stop "TECNOD142 failed for connectivities"
!!$
!!$        
!!$    end block add_data
!!$    
!!$    close_file: block
!!$
!!$      integer(kind=4) :: status
!!$
!!$      status = TECEND142()
!!$      
!!$      if (status .ne. 0) error stop "TECEND142 failed"
!!$      
!!$    end block close_file

    
    ! Open resource
    path = trim(filename)

    open(unit=90, file=path, iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    ! Write header
    write(90, *) 'TITLE = "FVM-Laplace"'
    write(90, *) 'VARIABLES = "x" "y"  "T"'

    

    !-----------------------------------------------------------------!
    ! Write Triangles/Quads (works only for homogeneous elements)
    !-----------------------------------------------------------------!

!!    if ( maxval(this % grid % num_cell_vertices) .eq. 4 ) then
!!       write(90, *) &
!!            & 'ZONE T="Temperature", N=', this % grid % num_vertices, &
!!            & ', E=', this % grid % num_cells, &
!!            & ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
!!    else
!!       write(90, *) &
!!            & 'ZONE T="Temperature", N=', this % grid % num_vertices, &
!!            & ', E=', this % grid % num_cells, &
!!            & ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
!!    end if

    write(90, *) &
         & 'ZONE T="Temperature", N=', mesh % num_vertices, &
         & ', E=', mesh % num_cells, &
         & ', DATAPACKING=BLOCK, ZONETYPE=FEPOLYGON'
    
    ! Write vertices
    do i = 1, mesh % num_vertices
       write(90,*) mesh % vertices(1:2,i), phiv(i)
    end do
    
    ! Write cell connectivities
    do i = 1, mesh % num_cells
       write(90,*) mesh % cell_vertices(1 : mesh % num_cell_vertices(i), i)
    end do

    !-----------------------------------------------------------------!
    ! Write Other types of elements too
    !-----------------------------------------------------------------!

    ! Close resource
    close(90)
    
    if (allocated(path)) deallocate(path)
    if (allocated(new_name)) deallocate(new_name)
    if (allocated(phiv)) deallocate(phiv)

  end subroutine write_solution


end module class_tecplot_writer
