module class_gmsh_loader

  ! import dependencies
  use iso_fortran_env       , only : dp => real64
  use interface_mesh_loader , only : mesh_loader
  use class_file            , only : file
  use class_string          , only : string
  use class_set             , only : set

  implicit none

  private
  public :: gmsh_loader

  !-------------------------------------------------------------------!
  ! Interface to construct a mesh_loader for GMSH
  !-------------------------------------------------------------------!
  
  interface gmsh_loader
     module procedure create
  end interface gmsh_loader
  
  type, extends(mesh_loader) :: gmsh_loader

     type(file) :: file ! mesh file

   contains

     ! Implement deferred procedure from interface
     procedure :: get_mesh_data

  end type gmsh_loader

contains

  type(gmsh_loader) function create(filename) result (this)

    type(character(*)), intent(in) :: filename

    this % file = file(filename)

  end function create


  !====================================================================!
  ! Supply all information needed to create a mesh object
  !====================================================================!

  subroutine get_mesh_data(this, &
       & num_vertices, vertex_numbers, vertex_tags , vertices ,  & 
       & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
       & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices   )
    
    ! Arguments
    class(gmsh_loader)  , intent(in)   :: this

    ! Vertices
    integer , intent(out)              :: num_vertices
    integer , intent(out), allocatable :: vertex_numbers(:)
    integer , intent(out), allocatable :: vertex_tags(:)
    real(dp), intent(out), allocatable :: vertices(:,:)

    ! Edges
    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_numbers(:)
    integer, intent(out), allocatable :: edge_tags(:)
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)

    ! Faces    
    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)

    ! Cells
    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)

    ! Local
    type(string), allocatable, dimension(:) :: lines    
    
    ! Mesh tag
    integer :: idx_start_mesh
    integer :: idx_end_mesh

    ! Physical_Names tag    
    integer :: idx_start_physical_names
    integer :: idx_end_physical_names

    ! Nodes tag
    integer :: idx_start_nodes
    integer :: idx_end_nodes

    ! Elements tag
    integer :: idx_start_elements
    integer :: idx_end_elements

    integer :: face_idx, cell_idx, vertex_idx
    
    ! Collection to store face information
    type(set) :: set_face_vertices
    type(set) :: set_face_tags
    type(set) :: set_face_numbers

    integer, allocatable :: tagged_face_tags(:,:)
    integer, allocatable :: tagged_face_numbers(:,:)

    ! Load the mesh into memory
    write(*,'(a,a)') "Loading mesh file :", this % file % filename
    call this % file % read_lines(lines)

    write(*,'(a)') "Identifying tags..."
    call find_tags(lines, &
         & idx_start_mesh           , idx_end_mesh  , &
         & idx_start_physical_names , idx_end_physical_names , &
         & idx_start_nodes          , idx_end_nodes, &
         & idx_start_elements       , idx_end_elements)
    
    write(*,'(a,i8,i8)') "mesh           : " , idx_start_mesh           , idx_end_mesh
    write(*,'(a,i8,i8)') "physical names : " , idx_start_physical_names , idx_end_physical_names
    write(*,'(a,i8,i8)') "nodes          : " , idx_start_nodes          , idx_end_nodes
    write(*,'(a,i8,i8)') "elements       : " , idx_start_elements       , idx_end_elements   

    write(*,'(a)') "Reading vertices..."
    
    ! Process nodes
    process_nodes: block
      
      integer                   :: num_tokens
      type(string), allocatable :: tokens(:)
      integer                   :: ivertex

      associate(vlines => lines(idx_start_nodes+2:idx_end_nodes-1))
      
        ! Set the number of vertices
        num_vertices = size(vlines)

        allocate(vertices(3, num_vertices))
        vertices = 0
        
        allocate(vertex_numbers(num_vertices))
        vertex_numbers = 0        

        allocate(vertex_tags(num_vertices))
        vertex_tags = 0

        ! Parse lines and store vertices
        do concurrent(ivertex=1:num_vertices)
           
           ! Get the numbers of tokens and tokens
           call vlines(ivertex) % tokenize(" ", num_tokens, tokens)
           
           ! First token is the vertex number
           vertex_numbers(ivertex) = ivertex !tokens(1) % asinteger()
           
           ! Second, third and fourth tokens are the coordinates
           vertices(:,ivertex) = tokens(2:4) % asreal()
           
        end do
   
        if (allocated(tokens)) deallocate(tokens)

      end associate

    end block process_nodes
    
    write(*,'(a)') "Reading elements..."

    elements : block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens    
      integer                   :: num_lines, iline

      associate(elines => lines(idx_start_elements+2:idx_end_elements-1))

        ! Extract start and end indices of different mesh tags used by
        ! GMSH
        num_lines = size(elines)

        ! Zero counters for elements
        num_edges = 0
        num_faces = 0
        num_cells = 0

        ! Count the number of cells present in elements
        do iline = 1, num_lines

           call elines(iline) % tokenize(" ", num_tokens, tokens)

           if (tokens(2) % asinteger() .eq. 2) then
              ! Triangular element
              num_cells = num_cells + 1
           else if (tokens(2) % asinteger() .eq. 3) then
              ! Quadrilateral element
              num_cells = num_cells + 1
           end if

        end do

        ! Allocate space based on counted num of cells
        allocate(cell_numbers(num_cells))
        cell_numbers = 0

        allocate(num_cell_vertices(num_cells))
        num_cell_vertices = 0

        allocate(cell_vertices(4,num_cells))
        cell_vertices = 0

        allocate(cell_tags(num_cells))
        cell_tags = 0

        ! Create space for processing face information
        set_face_vertices = set(2, 4*num_cells)
        set_face_tags = set(1, 4*num_cells)
        set_face_numbers = set(1, 4*num_cells)
        
        face_idx   = 0
        cell_idx   = 0
        vertex_idx = 0
        
        do iline = 1, num_lines

           call elines(iline) % tokenize(" ", num_tokens, tokens)

           ! Line element
           if (tokens(2) % asinteger() .eq. 1) then

              ! Carry out processing of physically tagged faces
              face_idx = face_idx + 1
              
              call set_face_vertices     % add_entry([tokens(6:7) % asinteger()])
              call set_face_tags         % add_entry([tokens(5) % asinteger()])
              call set_face_numbers      % add_entry([face_idx])

              ! Triangular element
           else if (tokens(2) % asinteger() .eq. 2) then

              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(5) % asinteger()
              cell_vertices(1:3,cell_idx)  = tokens(6:8) % asinteger()
              num_cell_vertices(cell_idx)  = 3

              ! Quadrilateral element
           else if (tokens(2) % asinteger() .eq. 3) then

              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(5) % asinteger()
              cell_vertices(1:4,cell_idx)  = tokens(6:9) % asinteger()
              num_cell_vertices(cell_idx)  = 4

              ! Node
           else if (tokens(2) % asinteger() .eq. 15) then

              ! Complete vertex processing with tagging
              vertex_idx = vertex_idx + 1
              vertex_tags(vertex_idx) = tokens(5) % asinteger()
              
           else

              print *, 'unknown element number', tokens(2) % asinteger()
              error stop

           end if

           if (allocated(tokens)) deallocate(tokens)

        end do

      end associate

    end block elements

    ! Extract the remaining faces based on cell vertices
    process_faces: block

      type(integer) :: idx(2)
      type(integer) :: icell, iverpair, iface
      
      ! Make ordered pair of vertices as faces in 2D
      do icell = 1, num_cells
         do iverpair = 1, num_cell_vertices(icell)
            if (iverpair .eq. num_cell_vertices(icell)) then
               idx = [iverpair, 1]
            else
               idx = [iverpair, iverpair+1]
            end if
            ! Add ordered pair of integers into set
            call set_face_vertices % add_entry(cell_vertices(idx, icell))
         end do
      end do
      
      ! Get the relevant entires into memory
      call set_face_vertices % get_entries(face_vertices)
      num_faces = size(face_vertices,dim=2)

      ! Set face numbers
      call set_face_numbers % get_entries(tagged_face_numbers)

      allocate(face_numbers(num_faces))
      do concurrent(iface=1:num_faces)
         face_numbers(iface) = iface
      end do
      print *, face_numbers
      !face_numbers(1:face_idx) = tagged_face_numbers(:,1)

      print *, 'num_faces', num_faces
      
      ! Set face tags
      call set_face_tags % get_entries(tagged_face_tags)
      allocate(face_tags(num_faces))
      do concurrent(iface=1:num_faces)
         face_tags(iface) = iface
      end do
      !face_tags(1:face_idx) = tagged_face_tags(:,1)

      ! Linear face has two vertices
      allocate(num_face_vertices(num_faces))
      num_face_vertices = 2

    end block process_faces

    deallocate(lines)
    
    process_edges : block

      ! Edges are not in 2D
      allocate(edge_numbers(num_edges))
      edge_numbers = 0

      allocate(edge_vertices(4,num_edges))
      edge_vertices = 0

      allocate(num_edge_vertices(num_edges))
      num_edge_vertices = 0

      allocate(edge_tags(num_edges))
      edge_tags = 0

    end block process_edges


    if (num_faces .gt. 4*num_cells) error stop
  end subroutine get_mesh_data

  pure subroutine find_tags(lines, &
       & idx_start_mesh           , idx_end_mesh  , &
       & idx_start_physical_names , idx_end_physical_names , &
       & idx_start_nodes          , idx_end_nodes, &
       & idx_start_elements       , idx_end_elements)
    
    ! Arguments
    type(string)       , intent(in) :: lines(:)

    ! Mesh tag
    integer, intent(out) :: idx_start_mesh
    integer, intent(out) :: idx_end_mesh

    ! Physical_Names tag    
    integer, intent(out) :: idx_start_physical_names
    integer, intent(out) :: idx_end_physical_names

    ! Nodes tag
    integer, intent(out) :: idx_start_nodes
    integer, intent(out) :: idx_end_nodes

    ! Elements tag
    integer, intent(out) :: idx_start_elements
    integer, intent(out) :: idx_end_elements

    character(len=*), parameter :: BEGIN_MESH           = "$MeshFormat"
    character(len=*), parameter :: END_MESH             = "$EndMeshFormat"  
    character(len=*), parameter :: BEGIN_PHYSICAL_NAMES = "$PhysicalNames"
    character(len=*), parameter :: END_PHYSICAL_NAMES   = "$EndPhysicalNames"    
    character(len=*), parameter :: BEGIN_NODES          = "$Nodes"
    character(len=*), parameter :: END_NODES            = "$EndNodes"  
    character(len=*), parameter :: BEGIN_ELEMENTS       = "$Elements"
    character(len=*), parameter :: END_ELEMENTS         = "$EndElements"  

    integer :: num_lines, iline
    
    ! Extract start and end indices of different mesh tags used by
    ! GMSH
    num_lines = size(lines)
    do concurrent(iline = 1:num_lines)

       ! Find mesh start and end
       if (index(lines(iline) % str, BEGIN_MESH) .eq. 1) then
          idx_start_mesh = iline
       end if
       if (index(lines(iline) % str, END_MESH) .eq. 1) then
          idx_end_mesh = iline
       end if

       ! Find physical_names start and end
       if (index(lines(iline) % str, BEGIN_PHYSICAL_NAMES) .eq. 1) then
          idx_start_physical_names = iline
       end if
       if (index(lines(iline) % str, END_PHYSICAL_NAMES) .eq. 1) then
          idx_end_physical_names = iline
       end if

       ! Find nodes start and end
       if (index(lines(iline) % str, BEGIN_NODES) .eq. 1) then
          idx_start_nodes = iline
       end if
       if (index(lines(iline) % str, END_NODES) .eq. 1) then
          idx_end_nodes = iline
       end if

       ! Find elements start and end
       if (index(lines(iline) % str, BEGIN_ELEMENTS) .eq. 1) then
          idx_start_elements = iline
       end if
       if (index(lines(iline) % str, END_ELEMENTS) .eq. 1) then
          idx_end_elements = iline
       end if

    end do

  end subroutine find_tags
  
end module class_gmsh_loader
