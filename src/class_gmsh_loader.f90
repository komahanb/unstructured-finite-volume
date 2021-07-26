module class_gmsh_loader

  ! import dependencies
  use iso_fortran_env       , only : dp => real64
  use interface_mesh_loader , only : mesh_loader
  use class_file            , only : file
  use class_string          , only : string
  use class_set             , only : set
  use class_list            , only : list

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
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
       & num_tags    , tag_numbers , tag_info )
    
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
    
    ! Tagging boundaries and domain with integers/strings
    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    type(string), allocatable, intent(out) :: tag_info(:)

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

    integer :: face_idx, cell_idx, edge_idx
    
    ! Collection to store face information
    type(set)  :: set_face_vertices
    type(set)  :: set_face_numbers
    type(list) :: list_face_tags

    ! Load the mesh into memory
    write(*,'(a,a)') "Loading mesh file :", this % file % filename
    call this % file % read_lines(lines)

    write(*,'(a)') "Identifying tags..."
    call find_tags(lines, &
         & idx_start_mesh           , idx_end_mesh  , &
         & idx_start_physical_names , idx_end_physical_names , &
         & idx_start_nodes          , idx_end_nodes, &
         & idx_start_elements       , idx_end_elements)

    write(*,*) "mesh           : " , idx_start_mesh           , idx_end_mesh
    write(*,*) "physical names : " , idx_start_physical_names , idx_end_physical_names
    write(*,*) "nodes          : " , idx_start_nodes          , idx_end_nodes
    write(*,*) "elements       : " , idx_start_elements       , idx_end_elements   
       
    process_mesh_version: block
      
      integer                   :: num_tokens
      type(string), allocatable :: tokens(:)
      
      write(*,'(a)') "Reading mesh information..."

      associate(mlines => lines(idx_start_mesh+1:idx_start_mesh+1))
        call mlines(1) % tokenize(" ", num_tokens, tokens)
        if (floor(tokens(1) % asreal()) > 2) then
           print *, "unsupported version of gmsh version ", tokens(1) % str
           stop
        end if
      end associate

      if (allocated(tokens)) deallocate(tokens)

      write(*,'(a)') "Reading mesh information completed..."
          
    end block process_mesh_version


    process_tags: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens    
      integer                   :: iline, length

      write(*,'(a)') "Reading physical tags..."

      associate(tag_lines => lines(idx_start_physical_names+1:idx_end_physical_names-1))

        ! Set the intent(out) variable for number of tags present 
        num_tags = tag_lines(1) % asinteger()
        
        ! Allocate space for other two return variables
        allocate(tag_info(num_tags))        
        allocate(tag_numbers(num_tags))
        tag_numbers = 0
        
        do concurrent(iline = 1: num_tags)

           ! Tokenize based on delimited space
           call tag_lines(iline+1) % tokenize(" ", num_tokens, tokens)
           
           ! Second tag is the tag number
           tag_numbers(iline) = tokens(2) % asinteger()

           ! Remove quotes on third tag
           length = len(tokens(3) % str)
           tag_info(iline) = string(tokens(3) % str(2:length-1))

        end do

      end associate

      if (allocated(tokens)) deallocate(tokens)

      write(*,'(4x,a,i0)') "num physical tags : ", num_tags

      write(*,'(4x,a)') "physical tags are : "

      write(*,*) tag_numbers

      call tag_info % print('(8x,a)')

      write(*,'(a)') "Reading physical tags completed..."

    end block process_tags
      

    ! Process nodes
    process_nodes: block
      
      integer                   :: num_tokens
      type(string), allocatable :: tokens(:)
      integer                   :: ivertex

      write(*,'(a)') "Reading vertices..."

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
           vertex_numbers(ivertex) = ivertex

           ! check for consistency of local node numebrs with gmsh
           if (ivertex .ne. tokens(1) % asinteger()) then
              error stop "inconsistent node numbering"
           end if           
           
           ! Second, third and fourth tokens are the coordinates
           vertices(:,ivertex) = tokens(2:4) % asreal()
           
        end do
   
        if (allocated(tokens)) deallocate(tokens)

      end associate

      if (allocated(tokens)) deallocate(tokens)

      write(*,'(4x,a,i0)') "num vertices   : ", num_vertices
      
      write(*,'(a)') "Reading vertices completed..."
      
    end block process_nodes

    elements : block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens    
      integer                   :: num_lines, iline
      type(logical)             :: added
      
      write(*,'(a)') "Reading elements..."

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
          
           if (tokens(2) % asinteger() .eq. 1) then

              ! 2-node line. 
              num_edges = num_edges + 1
              
           else if (tokens(2) % asinteger() .eq. 2) then

              ! 3-node triangle. 
              num_faces = num_faces + 1

           else if (tokens(2) % asinteger() .eq. 3) then

              ! 4-node quadrangle. 
              num_faces = num_faces + 1

           else if (tokens(2) % asinteger() .eq. 4) then
              
              ! 4-node tetrahedron.
              num_cells = num_cells + 1

           else if (tokens(2) % asinteger() .eq. 5) then
              
              ! 8-node hexahedron. 
              num_cells = num_cells + 1

           else if (tokens(2) % asinteger() .eq. 6) then

              ! 6-node prism
              num_cells = num_cells + 1

           else if (tokens(2) % asinteger() .eq. 7) then
                 
              ! 5-node prism
              num_cells = num_cells + 1

           else if (tokens(2) % asinteger() .eq. 15) then
                 
              ! 1-node point (skip)
                            
           else

              call tokens(2) % print()
              
              error stop "unsupported GMSH mesh element type"
           
           end if

        end do

        write(*,'(4x,a,i0)') "num edges      : ", num_edges
        write(*,'(4x,a,i0)') "num faces      : ", num_faces
        write(*,'(4x,a,i0)') "num cells      : ", num_cells       

        ! Allocate space for cells
        allocate(cell_numbers(num_cells))
        cell_numbers = 0

        allocate(num_cell_vertices(num_cells))
        num_cell_vertices = 0

        allocate(cell_vertices(8,num_cells)) ! max is 8noded tetrahedron
        cell_vertices = 0

        allocate(cell_tags(num_cells))
        cell_tags = 0
        
        ! Allocate space for faces
        allocate(face_numbers(num_faces))
        face_numbers = 0

        allocate(num_face_vertices(num_faces))
        num_face_vertices = 0

        allocate(face_vertices(4,num_faces)) ! max is 4 noded quadrilateral
        face_vertices = 0

        allocate(face_tags(num_faces))
        face_tags = 0

        ! Allocate space for edges
        allocate(edge_numbers(num_edges))
        edge_numbers = 0

        allocate(edge_vertices(2,num_edges))
        edge_vertices = 0

        allocate(num_edge_vertices(num_edges))
        num_edge_vertices = 0

        allocate(edge_tags(num_edges))
        edge_tags = 0

!!$        ! Create space for processing face information
!!$        set_face_vertices = set(2, 4*num_cells)
!!$        set_face_numbers  = set(1, 4*num_cells)
!!$        list_face_tags    = list(1, 4*num_cells)

        edge_idx   = 0
        face_idx   = 0
        cell_idx   = 0
        
        do iline = 1, num_lines

           ! elm-number elm-type number-of-tags < tag > … node-number-list
           call elines(iline) % tokenize(" ", num_tokens, tokens)

           if (tokens(2) % asinteger() .eq. 1) then

              !!$              ! Carry out processing of physically tagged faces
!!$              added = set_face_vertices % insert([tokens(6:7) % asinteger()])             
!!$              if (added .eqv. .true.) then               
!!$                 face_idx = face_idx + 1
!!$                 call set_face_numbers % add_entry([face_idx])
!!$                 call list_face_tags % add_entry([tokens(4) % asinteger()])                 
!!$              end if
!!$
!!$              error stop

              ! 2-node line.
              edge_idx = edge_idx + 1

              edge_numbers(edge_idx)       = edge_idx
              edge_tags(edge_idx)          = tokens(4) % asinteger()
              edge_vertices(1:2,edge_idx)  = tokens(6:6+2-1) % asinteger()
              num_edge_vertices(edge_idx)  = 2
              
           else if (tokens(2) % asinteger() .eq. 2) then

              ! 3-node triangle.
              face_idx = face_idx + 1

              face_numbers(face_idx)       = face_idx
              face_tags(face_idx)          = tokens(4) % asinteger()
              face_vertices(1:3,face_idx)  = tokens(6:6+3-1) % asinteger()
              num_face_vertices(face_idx)  = 3

           else if (tokens(2) % asinteger() .eq. 3) then

              ! 4-node quadrangle.
              face_idx = face_idx + 1

              face_numbers(face_idx)       = face_idx
              face_tags(face_idx)          = tokens(4) % asinteger()
              face_vertices(1:4,face_idx)  = tokens(6:6+4-1) % asinteger()
              num_face_vertices(face_idx)  = 4

           else if (tokens(2) % asinteger() .eq. 4) then
              
              ! 4-node tetrahedron.
              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(4) % asinteger()
              cell_vertices(1:4,cell_idx)  = tokens(6:6+4-1) % asinteger()
              num_cell_vertices(cell_idx)  = 4

           else if (tokens(2) % asinteger() .eq. 5) then
              
              ! 8-node hexahedron.
              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(4) % asinteger()
              cell_vertices(1:8,cell_idx)  = tokens(6:6+8-1) % asinteger()
              num_cell_vertices(cell_idx)  = 8
              

           else if (tokens(2) % asinteger() .eq. 6) then

              ! 6-node prism
              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(4) % asinteger()
              cell_vertices(1:6,cell_idx)  = tokens(6:6+6-1) % asinteger()
              num_cell_vertices(cell_idx)  = 6

           else if (tokens(2) % asinteger() .eq. 7) then
                 
              ! 5-node prism
              cell_idx = cell_idx + 1

              cell_numbers(cell_idx )      = cell_idx
              cell_tags(cell_idx)          = tokens(4) % asinteger()
              cell_vertices(1:5,cell_idx)  = tokens(6:6+5-1) % asinteger()
              num_cell_vertices(cell_idx)  = 5

           else if (tokens(2) % asinteger() .eq. 15) then

              ! 1-node point (skip)
              
           else

              call tokens(2) % print()
              
              error stop "unsupported GMSH mesh element type"
           
           end if

           if (allocated(tokens)) deallocate(tokens)

        end do

      end associate
      
      deallocate(lines)

      write(*,'(a)') "Reading elements completed..."

    end block elements

!!$    error stop
!!$
!!$    write(*,'(a)') "Finding faces..."
!!$    
!!$    ! Extract the remaining faces based on cell vertices
!!$    process_faces: block
!!$
!!$      type(integer) :: idx(2)
!!$      type(integer) :: icell, iverpair
!!$      type(logical) :: added
!!$      integer, allocatable :: lface_numbers(:,:), lface_tags(:,:)
!!$
!!$      ! Make ordered pair of vertices as faces in 2D
!!$      do icell = 1, num_cells
!!$         do iverpair = 1, num_cell_vertices(icell)
!!$
!!$            ! Figure out a vertex pair that makes a face
!!$            if (iverpair .eq. num_cell_vertices(icell)) then
!!$               idx = [iverpair, 1]
!!$            else
!!$               idx = [iverpair, iverpair+1]
!!$            end if
!!$
!!$            ! Add ordered pair of integers into set
!!$            added = set_face_vertices % insert(cell_vertices(idx, icell))
!!$            if (added .eqv. .true.) then               
!!$               face_idx = face_idx + 1
!!$               call set_face_numbers % add_entry([face_idx])
!!$               call list_face_tags % add_entry([0])
!!$            end if
!!$
!!$         end do
!!$      end do
!!$
!!$      num_faces = face_idx
!!$
!!$      allocate(num_face_vertices(num_faces))
!!$      num_face_vertices = 2
!!$
!!$      call set_face_vertices % get_entries(face_vertices)
!!$
!!$      ! Set face numbers
!!$      call set_face_numbers % get_entries(lface_numbers)
!!$      allocate(face_numbers(num_faces))
!!$      face_numbers = lface_numbers(1,:)
!!$      deallocate(lface_numbers)
!!$
!!$      ! Set face tags
!!$      call list_face_tags % get_entries(lface_tags)
!!$      allocate(face_tags(num_faces))
!!$      face_tags = lface_tags(1,:)
!!$      deallocate(lface_tags)
!!$
!!$    end block process_faces
!!$    
!!$    write(*,'(a)') "Processing edges..."
!!$
!!$    process_edges : block
!!$
!!$      ! Edges are not in 2D
!!$      allocate(edge_numbers(num_edges))
!!$      edge_numbers = 0
!!$
!!$      allocate(edge_vertices(2,num_edges))
!!$      edge_vertices = 0
!!$
!!$      allocate(num_edge_vertices(num_edges))
!!$      num_edge_vertices = 0
!!$
!!$      allocate(edge_tags(num_edges))
!!$      edge_tags = 0
!!$
!!$    end block process_edges
!!$
!!$    if (num_faces .gt. 4*num_cells) error stop
    
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
