!=====================================================================!
! Module that reads gmsh mesh files in the MSH 4.1 format (gmsh 4.x,
! e.g. 4.15). The format is documented here:
! http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
!
! 4.1 differs from the legacy 2.2 in three ways the parser handles:
!   - $Entities: an element no longer carries a physical tag per line;
!     it belongs to a geometric entity that carries the tag. so the
!     per-element tag is a lookup  element -> (entityDim, entityTag) ->
!     physical tag.
!   - $Nodes is block-structured: per entity block, all node tags then
!     all coordinates (node tags may be sparse - mapped through find).
!   - $Elements is block-structured: per entity block, element lines
!     with no per-element tag list.
!
! 2.2 files are no longer read - regenerate with python meshgen/generate.py.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_gmsh_loader

  ! import dependencies
  use iso_fortran_env       , only : dp => real64
  use interface_mesh_loader , only : mesh_loader
  use class_file            , only : file
  use class_string          , only : string
  use module_mesh_utils     , only : find, elem_type_dimension, &
       & elem_type_vertex_count
  use module_verbosity      , only : verbosity

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

  impure type(gmsh_loader) function create(filename) result (this)

    type(character(*)), intent(in) :: filename

    ! msh 4.1 lines are long (entity bounding boxes, coordinate triples),
    ! well past the default 100-char buffer - read them whole.
    this % file = file(filename, 4096)

  end function create

  !====================================================================!
  ! Supply all information needed to create a mesh object
  !====================================================================!

  impure subroutine get_mesh_data(this, &
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

    ! Section markers
    integer :: idx_start_mesh     , idx_end_mesh
    integer :: idx_start_physical_names, idx_end_physical_names
    integer :: idx_start_entities , idx_end_entities
    integer :: idx_start_nodes    , idx_end_nodes
    integer :: idx_start_elements , idx_end_elements

    ! Entity -> physical-tag maps, one (tag, phys) pair list per dimension
    integer, allocatable :: ent_tag0(:), ent_phys0(:)   ! points
    integer, allocatable :: ent_tag1(:), ent_phys1(:)   ! curves
    integer, allocatable :: ent_tag2(:), ent_phys2(:)   ! surfaces
    integer, allocatable :: ent_tag3(:), ent_phys3(:)   ! volumes

    integer :: mesh_dim

    ! Load the mesh into memory
    if (verbosity .ge. 1) write(*,'(a,a)') "Loading mesh file :", this % file % filename
    call this % file % read_lines(lines)

    if (verbosity .ge. 1) write(*,'(a)') "Identifying tags..."
    call find_tags(lines, &
         & idx_start_mesh           , idx_end_mesh  , &
         & idx_start_physical_names , idx_end_physical_names , &
         & idx_start_entities       , idx_end_entities , &
         & idx_start_nodes          , idx_end_nodes, &
         & idx_start_elements       , idx_end_elements)

    if (verbosity .ge. 1) then
       write(*,*) "mesh           : " , idx_start_mesh           , idx_end_mesh
       write(*,*) "physical names : " , idx_start_physical_names , idx_end_physical_names
       write(*,*) "entities       : " , idx_start_entities       , idx_end_entities
       write(*,*) "nodes          : " , idx_start_nodes          , idx_end_nodes
       write(*,*) "elements       : " , idx_start_elements       , idx_end_elements
    end if

    process_mesh_version: block

      integer                   :: num_tokens
      type(string), allocatable :: tokens(:)

      if (verbosity .ge. 1) write(*,'(a)') "Reading mesh information..."

      associate(mlines => lines(idx_start_mesh+1:idx_start_mesh+1))
        call mlines(1) % tokenize(" ", num_tokens, tokens)
        if (floor(tokens(1) % asreal()) .ne. 4) then
           print *, "mesh format ", tokens(1) % str, &
                & " is not msh 4.x - regenerate with python meshgen/generate.py"
           error stop
        end if
      end associate

      if (allocated(tokens)) deallocate(tokens)

      if (verbosity .ge. 1) write(*,'(a)') "Reading mesh information completed..."

    end block process_mesh_version

    process_tags: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: iline

      if (verbosity .ge. 1) write(*,'(a)') "Reading physical tags..."

      associate(tag_lines => lines(idx_start_physical_names+1:idx_end_physical_names-1))

        ! Set the intent(out) variable for number of tags present
        num_tags = tag_lines(1) % asinteger()

        ! Allocate space for other two return variables
        allocate(tag_info(num_tags))
        allocate(tag_numbers(num_tags))
        tag_numbers = 0
        allocate(tag_physical_dimensions(num_tags))
        tag_physical_dimensions = 0

        do iline = 1, num_tags

           ! Tokenize based on delimited space
           call tag_lines(iline+1) % tokenize(" ", num_tokens, tokens)

           ! First token is the physical dim, second is the tag number
           tag_physical_dimensions(iline) = tokens(1) % asinteger()
           tag_numbers(iline)             = tokens(2) % asinteger()

           ! The name is the quoted string on the line. Take everything
           ! between the first and last double quote so names with spaces
           ! (e.g. "cylindrical wall") survive.
           associate(s => tag_lines(iline+1) % str)
             tag_info(iline) = string(s(index(s,'"')+1 : index(s,'"',back=.true.)-1))
           end associate

        end do

      end associate

      if (allocated(tokens)) deallocate(tokens)

      if (verbosity .ge. 1) then
         write(*,'(4x,a,i0)') "num physical tags : ", num_tags
         write(*,'(4x,a)') "physical tags are : "
         write(*,*) tag_numbers
         call tag_info % print('(8x,a)')
         write(*,'(a)') "Reading physical tags completed..."
      end if

    end block process_tags

    ! Build the entity -> physical-tag maps from $Entities. an element
    ! gets its physical tag from the entity it belongs to (4.1 change).
    process_entities: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: np, nc, ns, nv, ie, il, nphys

      if (verbosity .ge. 1) write(*,'(a)') "Reading entities..."

      il = idx_start_entities + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      np = tokens(1) % asinteger()
      nc = tokens(2) % asinteger()
      ns = tokens(3) % asinteger()
      nv = tokens(4) % asinteger()

      allocate(ent_tag0(np), ent_phys0(np)); ent_phys0 = 0
      allocate(ent_tag1(nc), ent_phys1(nc)); ent_phys1 = 0
      allocate(ent_tag2(ns), ent_phys2(ns)); ent_phys2 = 0
      allocate(ent_tag3(nv), ent_phys3(nv)); ent_phys3 = 0

      il = il + 1

      ! points: tag x y z numPhys phys... (numPhys at token 5)
      do ie = 1, np
         call lines(il) % tokenize(" ", num_tokens, tokens)
         ent_tag0(ie) = tokens(1) % asinteger()
         nphys = tokens(5) % asinteger()
         if (nphys .gt. 0) ent_phys0(ie) = tokens(6) % asinteger()
         il = il + 1
      end do

      ! curves / surfaces / volumes: tag bbox(6) numPhys phys... (numPhys at token 8)
      do ie = 1, nc
         call lines(il) % tokenize(" ", num_tokens, tokens)
         ent_tag1(ie) = tokens(1) % asinteger()
         nphys = tokens(8) % asinteger()
         if (nphys .gt. 0) ent_phys1(ie) = tokens(9) % asinteger()
         il = il + 1
      end do

      do ie = 1, ns
         call lines(il) % tokenize(" ", num_tokens, tokens)
         ent_tag2(ie) = tokens(1) % asinteger()
         nphys = tokens(8) % asinteger()
         if (nphys .gt. 0) ent_phys2(ie) = tokens(9) % asinteger()
         il = il + 1
      end do

      do ie = 1, nv
         call lines(il) % tokenize(" ", num_tokens, tokens)
         ent_tag3(ie) = tokens(1) % asinteger()
         nphys = tokens(8) % asinteger()
         if (nphys .gt. 0) ent_phys3(ie) = tokens(9) % asinteger()
         il = il + 1
      end do

      if (allocated(tokens)) deallocate(tokens)

      if (verbosity .ge. 1) write(*,'(a)') "Reading entities completed..."

    end block process_entities

    ! $Nodes is block-structured: per block, all node tags then all
    ! coordinates. node tags can be sparse and are mapped through find.
    process_nodes: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: numblocks, il, ib, k, i, ivert

      if (verbosity .ge. 1) write(*,'(a)') "Reading nodes..."

      il = idx_start_nodes + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      numblocks    = tokens(1) % asinteger()
      num_vertices = tokens(2) % asinteger()

      allocate(vertex_numbers(num_vertices)); vertex_numbers = 0
      allocate(vertex_tags(num_vertices))   ; vertex_tags    = 0
      allocate(vertices(3, num_vertices))   ; vertices       = 0.0_dp

      il    = il + 1
      ivert = 0

      do ib = 1, numblocks

         ! entityDim entityTag parametric numNodesInBlock
         call lines(il) % tokenize(" ", num_tokens, tokens)
         k  = tokens(4) % asinteger()
         il = il + 1

         ! the node tags (a bare integer per line - read directly, since
         ! tokenize needs a delimiter and there is none on a single value)
         do i = 1, k
            read(lines(il) % str, *) vertex_numbers(ivert+i)
            il = il + 1
         end do

         ! the matching coordinates
         do i = 1, k
            read(lines(il) % str, *) vertices(1, ivert+i), vertices(2, ivert+i), vertices(3, ivert+i)
            il = il + 1
         end do

         ivert = ivert + k

      end do

      if (allocated(tokens)) deallocate(tokens)

      if (verbosity .ge. 1) then
         write(*,'(4x,a,i0)') "num vertices   : ", num_vertices
         write(*,'(a)') "Reading nodes completed..."
      end if

    end block process_nodes

    ! $Elements is block-structured. the top dimension are cells, one
    ! lower are faces, two lower are edges (so 2d and 3d both work). the
    ! physical tag of each element comes from its block's entity.
    process_elements: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: numblocks, il0, il, ib, k, i, j
      integer                   :: bdim, btag, etype, nv, edim, phys
      integer                   :: cell_idx, face_idx, edge_idx

      if (verbosity .ge. 1) write(*,'(a)') "Reading elements..."

      il = idx_start_elements + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      numblocks = tokens(1) % asinteger()
      il0       = idx_start_elements + 2

      ! pass A: the mesh (top) dimension from the element types present
      mesh_dim = 0
      il = il0
      do ib = 1, numblocks
         call lines(il) % tokenize(" ", num_tokens, tokens)
         etype    = tokens(3) % asinteger()
         k        = tokens(4) % asinteger()
         mesh_dim = max(mesh_dim, elem_type_dimension(etype))
         il = il + 1 + k
      end do

      ! pass B: count cells / faces / edges by relative dimension
      num_cells = 0
      num_faces = 0
      num_edges = 0
      il = il0
      do ib = 1, numblocks
         call lines(il) % tokenize(" ", num_tokens, tokens)
         etype = tokens(3) % asinteger()
         k     = tokens(4) % asinteger()
         edim  = elem_type_dimension(etype)
         if (edim .eq. mesh_dim) then
            num_cells = num_cells + k
         else if (edim .eq. mesh_dim - 1) then
            num_faces = num_faces + k
         else if (edim .eq. mesh_dim - 2) then
            num_edges = num_edges + k
         end if
         il = il + 1 + k
      end do

      if (verbosity .ge. 1) then
         write(*,'(4x,a,i0)')  "num edges                 : ", num_edges
         write(*,'(4x,a,2i0)') "num faces [boundary]      : ", num_faces
         write(*,'(4x,a,i0)')  "num cells                 : ", num_cells
      end if

      ! Allocate space for cells (max 8 noded hexahedron)
      allocate(cell_numbers(num_cells))      ; cell_numbers      = 0
      allocate(num_cell_vertices(num_cells)) ; num_cell_vertices = 0
      allocate(cell_vertices(8,num_cells))   ; cell_vertices     = 0
      allocate(cell_tags(num_cells))         ; cell_tags         = 0
      allocate(cell_types(num_cells))        ; cell_types        = 0

      ! Allocate space for faces (max 4 noded quadrilateral)
      allocate(face_numbers(num_faces))      ; face_numbers      = 0
      allocate(num_face_vertices(num_faces)) ; num_face_vertices = 0
      allocate(face_vertices(4,num_faces))   ; face_vertices     = 0
      allocate(face_tags(num_faces))         ; face_tags         = 0
      allocate(face_types(num_faces))        ; face_types        = 0

      ! Allocate space for edges
      allocate(edge_numbers(num_edges))      ; edge_numbers      = 0
      allocate(edge_vertices(2,num_edges))   ; edge_vertices     = 0
      allocate(num_edge_vertices(num_edges)) ; num_edge_vertices = 0
      allocate(edge_tags(num_edges))         ; edge_tags         = 0
      allocate(edge_types(num_edges))        ; edge_types        = 0

      ! pass C: fill, taking the physical tag from the block's entity
      cell_idx = 0
      face_idx = 0
      edge_idx = 0
      il = il0
      do ib = 1, numblocks

         ! entityDim entityTag elementType numElementsInBlock
         call lines(il) % tokenize(" ", num_tokens, tokens)
         bdim  = tokens(1) % asinteger()
         btag  = tokens(2) % asinteger()
         etype = tokens(3) % asinteger()
         k     = tokens(4) % asinteger()

         edim  = elem_type_dimension(etype)
         nv    = elem_type_vertex_count(etype)
         phys  = entity_phys(bdim, btag)

         il = il + 1

         do i = 1, k

            ! elementTag node1 node2 ...
            call lines(il) % tokenize(" ", num_tokens, tokens)

            if (edim .eq. mesh_dim) then

               cell_idx = cell_idx + 1
               cell_numbers(cell_idx)      = tokens(1) % asinteger()
               cell_types(cell_idx)        = etype
               cell_tags(cell_idx)         = phys
               num_cell_vertices(cell_idx) = nv
               do j = 1, nv
                  cell_vertices(j,cell_idx) = find(vertex_numbers, tokens(1+j) % asinteger())
               end do

            else if (edim .eq. mesh_dim - 1) then

               face_idx = face_idx + 1
               face_numbers(face_idx)      = tokens(1) % asinteger()
               face_types(face_idx)        = etype
               face_tags(face_idx)         = phys
               num_face_vertices(face_idx) = nv
               do j = 1, nv
                  face_vertices(j,face_idx) = find(vertex_numbers, tokens(1+j) % asinteger())
               end do

            else if (edim .eq. mesh_dim - 2) then

               edge_idx = edge_idx + 1
               edge_numbers(edge_idx)      = tokens(1) % asinteger()
               edge_types(edge_idx)        = etype
               edge_tags(edge_idx)         = phys
               num_edge_vertices(edge_idx) = nv
               do j = 1, nv
                  edge_vertices(j,edge_idx) = find(vertex_numbers, tokens(1+j) % asinteger())
               end do

            end if

            il = il + 1

         end do

      end do

      if (allocated(tokens)) deallocate(tokens)

      if (count(face_tags .ne. 0) .eq. 0) then
         write(*,*) "untagged faces exist - mesh not useful for simulation"
      end if

      if (verbosity .ge. 1) write(*,'(a)') "Reading elements completed..."

    end block process_elements

    deallocate(lines)

  contains

    !==================================================================!
    ! Physical tag of the entity (edim, etag), 0 if it has no group
    !==================================================================!

    pure integer function entity_phys(edim, etag)

      integer, intent(in) :: edim, etag

      select case (edim)
      case (0)
         entity_phys = ent_phys0(find(ent_tag0, etag))
      case (1)
         entity_phys = ent_phys1(find(ent_tag1, etag))
      case (2)
         entity_phys = ent_phys2(find(ent_tag2, etag))
      case (3)
         entity_phys = ent_phys3(find(ent_tag3, etag))
      case default
         entity_phys = 0
      end select

    end function entity_phys

  end subroutine get_mesh_data

  !====================================================================!
  ! Scan the file for the start/end line of each section we need
  !====================================================================!

  pure subroutine find_tags(lines, &
       & idx_start_mesh           , idx_end_mesh  , &
       & idx_start_physical_names , idx_end_physical_names , &
       & idx_start_entities       , idx_end_entities , &
       & idx_start_nodes          , idx_end_nodes, &
       & idx_start_elements       , idx_end_elements)

    ! Arguments
    type(string)       , intent(in) :: lines(:)

    integer, intent(out) :: idx_start_mesh           , idx_end_mesh
    integer, intent(out) :: idx_start_physical_names , idx_end_physical_names
    integer, intent(out) :: idx_start_entities       , idx_end_entities
    integer, intent(out) :: idx_start_nodes          , idx_end_nodes
    integer, intent(out) :: idx_start_elements       , idx_end_elements

    character(len=*), parameter :: BEGIN_MESH           = "$MeshFormat"
    character(len=*), parameter :: END_MESH             = "$EndMeshFormat"
    character(len=*), parameter :: BEGIN_PHYSICAL_NAMES = "$PhysicalNames"
    character(len=*), parameter :: END_PHYSICAL_NAMES   = "$EndPhysicalNames"
    character(len=*), parameter :: BEGIN_ENTITIES       = "$Entities"
    character(len=*), parameter :: END_ENTITIES         = "$EndEntities"
    character(len=*), parameter :: BEGIN_NODES          = "$Nodes"
    character(len=*), parameter :: END_NODES            = "$EndNodes"
    character(len=*), parameter :: BEGIN_ELEMENTS       = "$Elements"
    character(len=*), parameter :: END_ELEMENTS         = "$EndElements"

    integer :: num_lines, iline

    ! Sections absent unless found
    idx_start_entities = 0
    idx_end_entities   = 0

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

       ! Find entities start and end
       if (index(lines(iline) % str, BEGIN_ENTITIES) .eq. 1) then
          idx_start_entities = iline
       end if
       if (index(lines(iline) % str, END_ENTITIES) .eq. 1) then
          idx_end_entities = iline
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
