!=====================================================================!
! Module that reads gmsh mesh files. Currently supports only version 2.
!
! The format information is present here:
! http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
!=====================================================================!

module class_gmsh_loader

  ! import dependencies
  use iso_fortran_env       , only : dp => real64
  use interface_mesh_loader , only : mesh_loader
  use class_file            , only : file
  use class_string          , only : string
  use module_mesh_utils     , only : find, elem_type_face_count

  implicit none

  private
  public :: gmsh_loader

  !-------------------------------------------------------------------!
  ! Interface to construct a mesh_loader for GMSH
  !-------------------------------------------------------------------!

  interface gmsh_loader
     module procedure create
  end interface gmsh_loader

  !-------------------------------------------------------------------!
  ! gmsh_reader datatype
  !-------------------------------------------------------------------!
  
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
       & cell_types  , face_types    , edge_types  , &
       & num_tags    , tag_numbers   , tag_physical_dimensions, tag_info )

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
    integer, intent(out), allocatable :: edge_types(:)

    ! Faces
    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    integer, intent(out), allocatable :: face_types(:)

    ! Cells
    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)
    integer, intent(out), allocatable :: cell_types(:)

    ! Tagging boundaries and domain with integers/strings
    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    integer     , allocatable, intent(out) :: tag_physical_dimensions(:)
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
        allocate(tag_physical_dimensions(num_tags))
        tag_physical_dimensions = 0

        do concurrent(iline = 1: num_tags)

           ! Tokenize based on delimited space
           call tag_lines(iline+1) % tokenize(" ", num_tokens, tokens)

           ! First token is the physical dim
           tag_physical_dimensions(iline) = tokens(1) % asinteger()

           ! Second token is the tag number
           tag_numbers(iline) = tokens(2) % asinteger()

           ! Remove quotes on third token
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
           vertex_numbers(ivertex) = tokens(1) % asinteger() ! ivertex

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
      integer                   :: vloc
      integer                   :: ivertex, num_vertices

      write(*,'(a)') "Reading elements..."

      !---------------------------------------------------------------!
      ! $Elements
      ! number-of-elements
      ! elm-number elm-type number-of-tags < tag > … node-number-list
      ! …
      ! $EndElements
      !---------------------------------------------------------------!

      associate(elines => lines(idx_start_elements+2:idx_end_elements-1))

        ! Extract start and end indices of different mesh tags used by
        ! GMSH
        num_lines = size(elines)

        ! Zero counters for elements
        num_edges = 0
        num_faces = 0
        num_cells = 0

        ! gives the number of integer tags that follow for the n-th
        ! element. By default, the first tag is the tag of the
        ! physical entity to which the element belongs; the second is
        ! the tag of the elementary model entity to which the element
        ! belongs; the third is the number of mesh partitions to which
        ! the element belongs, followed by the partition ids (negative
        ! partition ids indicate ghost cells). A zero tag is
        ! equivalent to no tag. Gmsh and most codes using the MSH 2
        ! format require at least the first two tags (physical and
        ! elementary tags).

        ! Count the number of cells present in elements https://gmsh.info/doc/texinfo/gmsh.html
        do iline = 1, num_lines

           call elines(iline) % tokenize(" ", num_tokens, tokens)

           select case (tokens(2) % asinteger())
           case (1)
              ! 2-node line.
              num_edges = num_edges + 1
           case (2:3)
              ! 3-node triangle, 4-node quadrangle
              num_faces = num_faces + 1
           case (4:7)
              ! 4-node tetrahedron, 8-node hexahedron, 6-node prism, 5-node prism (pyramid)
              num_cells = num_cells + 1
           case (15)
              ! 1-node point (skip)
           case default
              call tokens(2) % print()
              error stop "unsupported GMSH mesh element type"
           end select

        end do

        write(*,'(4x,a,i0)')  "num edges                 : ", num_edges
        write(*,'(4x,a,2i0)') "num faces [boundary]      : ", num_faces
        write(*,'(4x,a,i0)')  "num cells                 : ", num_cells

        ! Allocate space for cells
        allocate(cell_numbers(num_cells))
        cell_numbers = 0

        allocate(num_cell_vertices(num_cells))
        num_cell_vertices = 0

        allocate(cell_vertices(8,num_cells)) ! max is 8 noded tetrahedron
        cell_vertices = 0

        allocate(cell_tags(num_cells))
        cell_tags = 0

        allocate(cell_types(num_cells))
        cell_types = 0

        ! Allocate space for faces
        allocate(face_numbers(num_faces))
        face_numbers = 0

        allocate(num_face_vertices(num_faces))
        num_face_vertices = 0

        allocate(face_vertices(4,num_faces)) ! max is 4 noded quadrilateral
        face_vertices = 0

        allocate(face_tags(num_faces))
        face_tags = 0

        allocate(face_types(num_faces))
        face_types = 0

        ! Allocate space for edges
        allocate(edge_numbers(num_edges))
        edge_numbers = 0

        allocate(edge_vertices(2,num_edges))
        edge_vertices = 0

        allocate(num_edge_vertices(num_edges))
        num_edge_vertices = 0

        allocate(edge_tags(num_edges))
        edge_tags = 0

        allocate(edge_types(num_edges))
        edge_types = 0

        edge_idx   = 0
        face_idx   = 0
        cell_idx   = 0

        do iline = 1, num_lines

           ! elm-number elm-type number-of-tags < tag > … node-number-list
           call elines(iline) % tokenize(" ", num_tokens, tokens)

           if (tokens(2) % asinteger() .eq. 1) then

              ! 2-node line.
              edge_idx                               = edge_idx + 1
              num_vertices                           = 2
              edge_numbers(edge_idx)                 = tokens(1) % asinteger()
              edge_types(edge_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              edge_tags(edge_idx)                    = tokens(4) % asinteger()
              edge_vertices(1:num_vertices,edge_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_edge_vertices(edge_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 edge_vertices(ivertex, edge_idx) = find(vertex_numbers, edge_vertices(ivertex,edge_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 2) then

              ! 3-node triangle.
              face_idx                               = face_idx + 1
              num_vertices                           = 3
              face_numbers(face_idx)                 = tokens(1) % asinteger()
              face_types(face_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              face_tags(face_idx)                    = tokens(4) % asinteger()
              face_vertices(1:num_vertices,face_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_face_vertices(face_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 face_vertices(ivertex, face_idx) = find(vertex_numbers, face_vertices(ivertex,face_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 3) then

              ! 4-node quadrangle.
              face_idx                               = face_idx + 1
              num_vertices                           = 4
              face_numbers(face_idx)                 = tokens(1) % asinteger()
              face_types(face_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              face_tags(face_idx)                    = tokens(4) % asinteger()
              face_vertices(1:num_vertices,face_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_face_vertices(face_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 face_vertices(ivertex, face_idx) = find(vertex_numbers, face_vertices(ivertex,face_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 4) then

              ! 4-node tetrahedron.
              cell_idx                               = cell_idx + 1
              num_vertices                           = 4
              cell_numbers(cell_idx)                 = tokens(1) % asinteger()
              cell_types(cell_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              cell_tags(cell_idx)                    = tokens(4) % asinteger()
              cell_vertices(1:num_vertices,cell_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_cell_vertices(cell_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 cell_vertices(ivertex, cell_idx) = find(vertex_numbers, cell_vertices(ivertex,cell_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 5) then

              ! 8-node hexahedron.
              cell_idx                                = cell_idx + 1
              num_vertices                            = 8
              cell_numbers(cell_idx)                  = tokens(1) % asinteger()
              cell_types(cell_idx)                    = tokens(2) % asinteger()
              vloc                                    = 4 + tokens(3) % asinteger()
              cell_tags(cell_idx)                     = tokens(4) % asinteger()
              cell_vertices(1:num_vertices, cell_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger() ! global numbers in parallel
              num_cell_vertices(cell_idx)             = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 cell_vertices(ivertex, cell_idx) = find(vertex_numbers, cell_vertices(ivertex,cell_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 6) then

              ! 6-node prism
              cell_idx                               = cell_idx + 1
              num_vertices                           = 6
              cell_numbers(cell_idx)                 = tokens(1) % asinteger()
              cell_types(cell_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              cell_tags(cell_idx)                    = tokens(4) % asinteger()
              cell_vertices(1:num_vertices,cell_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_cell_vertices(cell_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 cell_vertices(ivertex, cell_idx) = find(vertex_numbers, cell_vertices(ivertex,cell_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 7) then

              ! 5-node prism
              cell_idx                               = cell_idx + 1
              num_vertices                           = 5
              cell_numbers(cell_idx)                 = tokens(1) % asinteger()
              cell_types(cell_idx)                   = tokens(2) % asinteger()
              vloc                                   = 4 + tokens(3) % asinteger()
              cell_tags(cell_idx)                    = tokens(4) % asinteger()
              cell_vertices(1:num_vertices,cell_idx) = tokens(vloc:vloc+num_vertices-1) % asinteger()
              num_cell_vertices(cell_idx)            = num_vertices

              do concurrent (ivertex = 1 : num_vertices)
                 cell_vertices(ivertex, cell_idx) = find(vertex_numbers, cell_vertices(ivertex,cell_idx))
              end do

           else if (tokens(2) % asinteger() .eq. 15) then

              ! 1-node point (skip)

           else

              call tokens(2) % print()

              error stop "unsupported GMSH mesh element type"

           end if

           if (allocated(tokens)) deallocate(tokens)

        end do

      end associate

      if (count(face_tags .ne. 0) .eq.0 ) then
         write(*,*) "untagged faces exist - mesh not useful for simulation"
      end if

      deallocate(lines)

      write(*,'(a)') "Reading elements completed..."

    end block elements

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
