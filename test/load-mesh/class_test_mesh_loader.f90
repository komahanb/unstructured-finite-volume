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
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
       & num_tags    , tag_numbers   , tag_info &
       & )

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

    ! Tagging boundaries and domain with integers/strings
    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    type(string), allocatable, intent(out) :: tag_info(:)

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
      allocate(edge_vertices(0,num_edges))
      allocate(num_edge_vertices(num_edges))
      edge_numbers      = 0
      edge_vertices     = 0
      num_edge_vertices = 0

      ! Read the elements into memory and predict the number of entities
      write(*,'(a,a)') "Loading elements file : ", this % elem_file % filename
      call this % elem_file % read_lines(lines)
      num_lines = size(lines)

      ! Allocate space
      num_cells = num_lines      

      allocate(num_cell_vertices(num_cells))
      num_cell_vertices = 0

      allocate(cell_numbers(num_cells))
      cell_numbers = 0

      ! Do an intial read to know the kind of cells present
      do concurrent(icell=1:num_lines)
         call lines(icell) % tokenize(",", num_tokens)
         num_cell_vertices(icell)  = num_tokens
      end do

      allocate(cell_vertices(maxval(num_cell_vertices),num_cells))
      cell_vertices = 0
      
      ! Now read the data
      do concurrent (icell=1:num_lines)

         call lines(icell) % tokenize(",", num_tokens, tokens)
         
         if (num_tokens .ne. 0) then
            cell_numbers(icell) = icell
            cell_vertices(1:num_cell_vertices(icell),icell) = tokens % asinteger() 
         else
            error stop
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
      vertices       = 0
      vertex_numbers = 0

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
      faces = set(2,maxval(num_cell_vertices)*num_cells)

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
      allocate(cell_tags(num_cells)); cell_tags = 0
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

  end subroutine get_mesh_data

end module class_test_mesh_loader
