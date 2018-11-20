module class_test_mesh_loader

  ! import dependencies
  use iso_fortran_env       , only : dp => real64, error_unit
  use interface_mesh_loader , only : mesh_loader

  ! Helper classes for file handing and string manipulation
  use class_file            , only : file
  use class_string          , only : string
  use class_set             , only : set

  implicit none

  private
  public :: test_mesh_loader

  !-------------------------------------------------------------------!
  ! Interface to construct a mesh loader for test mesh
  !-------------------------------------------------------------------!
  
  interface test_mesh_loader
     module procedure create
  end interface test_mesh_loader
  
  type, extends(mesh_loader) :: test_mesh_loader

     type(file) :: coord_file
     type(file) :: elem_file

   contains

     ! Implement deferred procedure from interface
     procedure :: get_mesh_data

  end type test_mesh_loader

contains

  type(test_mesh_loader) function create(coord_filename, elem_filename) &
       & result (this)

    type(character(*)), intent(in) :: coord_filename
    type(character(*)), intent(in) :: elem_filename
    
    this % coord_file = file(coord_filename)
    this % elem_file  = file(elem_filename)

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
    class(test_mesh_loader)  , intent(in)   :: this

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

    load_elements: block

      ! Local variables
      type(string), allocatable, dimension(:) :: lines    
      type(string), allocatable               :: tokens(:)
      integer                                 :: num_tokens
      integer                                 :: icell = 0
      integer                                 :: num_lines = 0

      ! Allocate memory based on prediction and zero as necessary
      ! No edges in 2D. Edges are faces.
      num_edges = 0
      allocate(edge_numbers(num_edges))
      allocate(edge_tags(num_edges))
      allocate(edge_vertices(0,num_edges))
      allocate(num_edge_vertices(num_edges))
      edge_numbers      = 0
      edge_tags         = 0
      edge_vertices     = 0
      num_edge_vertices = 0

      ! Read the elements into memory and predict the number of entities
      write(*,'(a,a)') "Loading elements file : ", this % elem_file % filename
      call this % elem_file % read_lines(lines)
      num_lines = size(lines)

      num_cells = num_lines      
      allocate(cell_numbers(num_cells))
      !allocate(cell_tags(num_cells))
      allocate(cell_vertices(3,num_cells))
      allocate(num_cell_vertices(num_cells))
      cell_numbers      = 0
      !cell_tags         = 0
      cell_vertices     = 0
      num_cell_vertices = 0

      do concurrent (icell=1:num_lines)

         call lines(icell) % tokenize(",", num_tokens, tokens)
         
         if (num_tokens .gt. 0) then

            cell_numbers(icell)       = icell
            cell_vertices(1:3,icell)  = tokens % asinteger()
            num_cell_vertices(icell)  = num_tokens ! hopefully 3

         end if

      end do
          
      if (allocated(tokens)) deallocate(tokens)
      if (allocated(lines)) deallocate(lines)

    end block load_elements

    ! Process nodes
    process_nodes: block

      ! Local variables
      type(string), allocatable, dimension(:) :: lines
      type(string), allocatable               :: tokens(:)
      integer                                 :: num_tokens
      integer                                 :: ivertex

      ! Read the coordinate into memory
      write(*,'(a,a)') "Loading coordinate file : ", this % coord_file % filename
      call this % coord_file % read_lines(lines)

      ! Set the number of vertices
      num_vertices = size(lines)
      allocate(vertices(3, num_vertices))
      allocate(vertex_numbers(num_vertices))
      allocate(vertex_tags(num_vertices))
      vertices       = 0
      vertex_numbers = 0
      vertex_tags    = 0

      ! Parse lines and store vertices
      do concurrent(ivertex=1:num_vertices)

         ! Get the numbers of tokens and tokens
         call lines(ivertex) % tokenize(",", num_tokens, tokens)

         if (num_tokens .gt. 0) then
            ! We create our own vertex number
            vertex_numbers(ivertex) = ivertex

            ! First and second tokens are the x and y coordinates
            vertices(1:2,ivertex) = tokens(:) % asreal()
         end if

      end do

      ! Tag vertices if they are a part of boundary tag (later) after
      ! reading the elements

      if (allocated(tokens)) deallocate(tokens)
      if (allocated(lines)) deallocate(lines)

    end block process_nodes

    face_finder: block         

      type(set)     :: faces
      type(integer) :: idx(2)
      type(integer) :: icell, iverpair, iface

      ! Create space for as many faces possible (rehash maybe?)
      faces = set(3*num_cells)

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
    
    ! Tag vertices, faces and cells
    !allocate(face_tags(num_faces))
    !face_tags         = 0
    
    ! All 16 inputs OK?

!!$
!!$    write(*,'(a)') "Identifying tags..."
!!$    call this % find_tags(lines, &
!!$         & idx_start_mesh           , idx_end_mesh  , &
!!$         & idx_start_physical_names , idx_end_physical_names , &
!!$         & idx_start_nodes          , idx_end_nodes, &
!!$         & idx_start_elements       , idx_end_elements)
!!$    write(*,'(a,i8,i8)') "mesh           : " , idx_start_mesh           , idx_end_mesh
!!$    write(*,'(a,i8,i8)') "physical names : " , idx_start_physical_names , idx_end_physical_names
!!$    write(*,'(a,i8,i8)') "nodes          : " , idx_start_nodes          , idx_end_nodes
!!$    write(*,'(a,i8,i8)') "elements       : " , idx_start_elements       , idx_end_elements   
!!$
!!$    write(*,'(a)') "Reading vertices... "        
!!$    call this % process_vertices(lines(idx_start_nodes+2:idx_end_nodes-1), &
!!$         & num_vertices, vertices, vertex_numbers, vertex_tags)
!!$    write(*,'(a,i8)') "number of vertices", num_vertices
!!$
!!$    ! How to find vertex tags?   
!!$    write(*,'(a)') "Reading elements... "
!!$    call this % process_elements(lines(idx_start_elements+2:idx_end_elements-1), &
!!$         & num_edges, edge_numbers, edge_tags, edge_vertices, num_edge_vertices, &              
!!$         & num_faces, face_numbers, face_tags, face_vertices, num_face_vertices, &
!!$         & num_cells, cell_numbers, cell_tags, cell_vertices, num_cell_vertices  )
!!$    write(*,'(a,i8)') "number of cells :", num_cells
!!$    write(*,'(a,i8)') "number of faces :", num_faces
!!$    write(*,'(a,i8)') "number of edges :", num_edges

  end subroutine get_mesh_data

end module class_test_mesh_loader
