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

     ! Helper routines
     procedure :: find_tags
     procedure :: process_vertices
     procedure :: process_elements

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

    ! Load the mesh into memory
    write(*,'(a,a)') "Loading mesh file : ", this % file % filename
    call this % file % read_lines(lines)
    ! call lines % print()  

    write(*,'(a)') "Identifying tags..."
    call this % find_tags(lines, &
         & idx_start_mesh           , idx_end_mesh  , &
         & idx_start_physical_names , idx_end_physical_names , &
         & idx_start_nodes          , idx_end_nodes, &
         & idx_start_elements       , idx_end_elements)
    write(*,'(a,i8,i8)') "mesh           : " , idx_start_mesh           , idx_end_mesh
    write(*,'(a,i8,i8)') "physical names : " , idx_start_physical_names , idx_end_physical_names
    write(*,'(a,i8,i8)') "nodes          : " , idx_start_nodes          , idx_end_nodes
    write(*,'(a,i8,i8)') "elements       : " , idx_start_elements       , idx_end_elements   

    write(*,'(a)') "Reading vertices... "        
    call this % process_vertices(lines(idx_start_nodes+2:idx_end_nodes-1), &
         & num_vertices, vertices, vertex_numbers, vertex_tags)

    ! How to find vertex tags?   
    write(*,'(a)') "Reading elements... "
    call this % process_elements(lines(idx_start_elements+2:idx_end_elements-1), &
         & num_edges, edge_numbers, edge_tags, edge_vertices, num_edge_vertices, &              
         & num_faces, face_numbers, face_tags, face_vertices, num_face_vertices, &
         & num_cells, cell_numbers, cell_tags, cell_vertices, num_cell_vertices  )

    deallocate(lines)
    
    face_finder: block         

      type(set)     :: faces
      type(integer) :: idx(2)
      type(integer) :: icell, iverpair, iface

      ! Create space for as many faces possible (rehash maybe?)
      faces = set(maxval(num_cell_vertices)*num_cells)

      ! Make ordered pair of vertices as faces in 2D
      do icell = 1, num_cells
         do iverpair = 1, num_cell_vertices(icell)
            if (iverpair .eq. num_cell_vertices(icell)) then           
               idx = [iverpair, 1]
            else           
               idx = [iverpair, iverpair+1]
            end if
            ! Add ordered pair of integers into set
            call faces % add_entry(cell_vertices(idx, icell))
         end do
      end do

      ! Get the relevant entires into memory
      call faces % get_entries(face_vertices)
      num_faces = size(face_vertices,dim=2)

      ! Set face numbers
      allocate(face_numbers(num_faces))
      do concurrent(iface=1:num_faces)
         face_numbers(iface) = iface
      end do

      ! Linear face has two vertices
      allocate(num_face_vertices(num_faces))
      num_face_vertices = 2

    end block face_finder

    tag_find: block 

      ! Tag vertices, faces and cells
      allocate(vertex_tags(num_vertices)); vertex_tags = 0
      !allocate(cell_tags(num_cells)); cell_tags = 0
      allocate(face_tags(num_faces)); face_tags = 0
      allocate(edge_tags(num_edges)); edge_tags = 0

      ! Here the principle is that a cell with only one neighbour is
      ! tagged as boundary

      ! Then we loop through all such tagged cells as 'say 1', extract
      ! the three faces

      ! Maybe have a default option for user to tag things, if not we
      ! can do some processing internally after all the inverse
      ! information is obtained and tag these separately than interior
      ! nodes. But for advanced appliations, the user might want to
      ! supply tags externally.

      ! With this information it is not possible to Use geometry and
      ! x,y coords to find vertices, then faces, then cells, edges?

    end block tag_find
    
    num_edges = 0    
    
  end subroutine get_mesh_data

  pure subroutine find_tags(this, lines, &
       & idx_start_mesh           , idx_end_mesh  , &
       & idx_start_physical_names , idx_end_physical_names , &
       & idx_start_nodes          , idx_end_nodes, &
       & idx_start_elements       , idx_end_elements)
    
    ! Arguments
    class(gmsh_loader) , intent(in) :: this
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
  
  subroutine process_elements(this, &
       & lines, &
       & num_edges, edge_numbers, edge_tags, edge_vertices, num_edge_vertices, &              
       & num_faces, face_numbers, face_tags, face_vertices, num_face_vertices, &
       & num_cells, cell_numbers, cell_tags, cell_vertices, num_cell_vertices  &       
       & )

    class(gmsh_loader) , intent(in)   :: this
    type(string)       , intent(in)   :: lines(:)
    
    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_numbers(:)
    integer, intent(out), allocatable :: edge_tags(:)        
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)
    
    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)        
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    
    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)        
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)

    type(string), allocatable :: tokens(:)
    integer                   :: num_tokens    
    integer                   :: num_lines, iline
    integer                   :: face_idx, edge_idx, cell_idx, node_idx
    
    ! Extract start and end indices of different mesh tags used by
    ! GMSH
    num_lines = size(lines)

    !num_edges = 0
    !num_faces = 0
    num_cells = 0
    
    do iline = 1, num_lines
       
       call lines(iline) % tokenize(" ", num_tokens, tokens)
   
       if (tokens(2) % asinteger() .eq. 1) then
          ! Line element
     !     num_faces = num_faces + 1
       else if (tokens(2) % asinteger() .eq. 2) then
          ! Triangular element
          num_cells = num_cells + 1
       else if (tokens(2) % asinteger() .eq. 3) then
          ! Quadrilateral element
          num_cells = num_cells + 1
       else if (tokens(2) % asinteger() .eq. 15) then
          print *, 'skip node'
       else
          error stop
       end if
       
    end do

!!$    allocate(edge_numbers(num_edges))
!!$    allocate(face_numbers(num_faces))
    allocate(cell_numbers(num_cells))
!!$    edge_numbers = 0
!!$    face_numbers = 0
    cell_numbers = 0
!!$    
!!$    allocate(edge_tags(num_edges))
!!$    allocate(face_tags(num_faces))
    allocate(cell_tags(num_cells))
!!$    edge_tags = 0
!!$    face_tags = 0
    cell_tags = 0
    
!!$    allocate(edge_vertices(4,num_edges)) ! upto quads
!!$    allocate(face_vertices(2,num_faces)) ! linear faces
    allocate(cell_vertices(4,num_cells)) ! upto quads
!!$    edge_vertices = 0
!!$    face_vertices = 0
    cell_vertices = 0

!!$    allocate(num_edge_vertices(num_edges))
!!$    allocate(num_face_vertices(num_faces))
    allocate(num_cell_vertices(num_cells))

!!$    face_idx = 0
    cell_idx = 0
    !vertex_idx = 0
    do iline = 1, num_lines
       
       call lines(iline) % tokenize(" ", num_tokens, tokens)
       
       ! Line element
       if (tokens(2) % asinteger() .eq. 1) then
          
      !    face_idx = face_idx + 1
!!$
!!$          face_numbers(face_idx)      = face_idx
!!$          face_tags(face_idx)         = tokens(5) % asinteger()
!!$          face_vertices(:,face_idx)   = tokens(6:7) % asinteger()
!!$          num_face_vertices(face_idx) = 2
!!$          
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

!!$          vertex_idx = vertex_idx + 1
!!$          
!!$          cell_numbers(vertex_idx )      = iline
!!$          cell_tags(vertex_idx)          = tokens(5) % asinteger()
!!$          cell_vertices(1:4,vertex_idx)  = tokens(6:9) % asinteger()
!!$          num_cell_vertices(vertex_idx)  = 4

       else

          print *, 'unknown element number', tokens(2) % asinteger()
          error stop

       end if

       if (allocated(tokens)) deallocate(tokens)

    end do

  end subroutine process_elements
  
  pure subroutine process_vertices(this, lines, &
       & num_vertices, vertices, vertex_numbers, vertex_tags)
    
    ! Arguments
    class(gmsh_loader) , intent(in)               :: this
    type(string)       , intent(in)               :: lines(:)
    integer            , intent(out)              :: num_vertices
    real(dp)           , intent(out), allocatable :: vertices(:,:)
    integer            , intent(out), allocatable :: vertex_numbers(:)
    integer            , intent(out), allocatable :: vertex_tags(:)

    ! Process nodes
    process_nodes: block
      
      integer                   :: iline
      integer                   :: num_tokens
      integer                   :: ivertex
      type(string), allocatable :: tokens(:)

      ! Set the number of vertices
      num_vertices = size(lines)
      allocate(vertices(3, num_vertices))
      allocate(vertex_numbers(num_vertices))
!!$      allocate(vertex_tags(num_vertices))
      vertices       = 0
      vertex_numbers = 0
  !!$    vertex_tags    = 0
      
      ! Parse lines and store vertices
      do concurrent(ivertex=1:num_vertices)

         ! Get the numbers of tokens and tokens
         call lines(ivertex) % tokenize(" ", num_tokens, tokens)

         ! First token is the vertex number
         vertex_numbers(ivertex) = tokens(1) % asinteger()

         ! Second, third and fourth token are the coordinates
         vertices(:,ivertex) = tokens(2:) % asreal()

      end do

      if (allocated(tokens)) deallocate(tokens)

      ! Determine tags?
      
    end block process_nodes

  end subroutine process_vertices

end module class_gmsh_loader
