!=====================================================================!
! This module reads gmsh mesh files in the MSH 4.1 format (gmsh 4.x,
! e.g. 4.15). The format is documented here:
! http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
!
! 4.1 differs from the legacy 2.2 in three ways the parser handles:
!   - $Entities: an element line no longer lists a physical tag; the
!     element belongs to a geometric entity that stores the tag. So the
!     per-element tag is a lookup:
!
!        element -> (entityDim, entityTag) -> physical tag
!
!   - $Nodes is block-structured: per entity block, all node tags then
!     all coordinates (node tags may be sparse - mapped through findloc).
!   - $Elements is block-structured: per entity block, element lines
!     with no per-element tag list.
!
! Legacy 2.2 files are no longer read; regenerate the mesh with
! python meshgen/generate.py.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_gmsh_loader

  ! Import the dependencies.
  use util_precision  , only : dp
  use view_mesh_loader , only : mesh_loader
  use util_file            , only : file
  use util_string          , only : string
  use view_mesh_geometry   , only : element_dimension, &
       & element_num_vertices, widest_element
  use util_verbosity      , only : verbosity

  implicit none

  private
  public :: gmsh_loader

  !-------------------------------------------------------------------!
  ! The sections the parser reads, in the order they are numbered.
  !-------------------------------------------------------------------!

  integer, parameter :: num_sections = 5
  integer, parameter :: MESH = 1, PHYSICAL_NAMES = 2, ENTITIES = 3, NODES = 4, ELEMENTS = 5
  character(len=13), parameter :: sections(num_sections) = &
       & [character(len=13) :: 'MeshFormat', 'PhysicalNames', 'Entities', 'Nodes', 'Elements']

  !-------------------------------------------------------------------!
  ! The interface to construct a mesh_loader for GMSH.
  !-------------------------------------------------------------------!

  interface gmsh_loader
     module procedure create
  end interface gmsh_loader

  !-------------------------------------------------------------------!
  ! The gmsh_loader datatype.
  !-------------------------------------------------------------------!

  type, extends(mesh_loader) :: gmsh_loader

     type(file) :: file ! The mesh file.

   contains

     ! Implement the deferred procedure from the interface.
     procedure :: mesh_data

  end type gmsh_loader

contains

  !===================================================================!
  ! Constructor: point at the .msh file and size the line buffer.
  !===================================================================!

  impure type(gmsh_loader) function create(filename) result (this)

    type(character(*)), intent(in) :: filename

    !----------------------------------------------------------------!
    ! MSH 4.1 lines are long (entity bounding boxes, coordinate
    ! triples), well past the default 100-character buffer, so the
    ! file reads whole lines.
    !----------------------------------------------------------------!

    this % file = file(filename, 4096)

  end function create

  !====================================================================!
  ! Supply all the information needed to create a mesh object.
  !====================================================================!

  impure subroutine mesh_data(this, &
       & num_vertices, vertex_numbers, vertex_tags , vertices ,  &
       & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
       & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
       & cell_types  , face_types    , edge_types  , &
       & num_tags    , tag_numbers   , tag_physical_dimensions, tag_info )

    ! The arguments.
    class(gmsh_loader)  , intent(in)   :: this

    ! The vertices.
    integer , intent(out)              :: num_vertices
    integer , intent(out), allocatable :: vertex_numbers(:)
    integer , intent(out), allocatable :: vertex_tags(:)
    real(dp), intent(out), allocatable :: vertices(:,:)

    ! The edges.
    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_numbers(:)
    integer, intent(out), allocatable :: edge_tags(:)
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)
    integer, intent(out), allocatable :: edge_types(:)

    ! The faces.
    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    integer, intent(out), allocatable :: face_types(:)

    ! The cells.
    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)
    integer, intent(out), allocatable :: cell_types(:)

    ! These tag the boundaries and the domain with integers and strings.
    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    integer     , allocatable, intent(out) :: tag_physical_dimensions(:)
    type(string), allocatable, intent(out) :: tag_info(:)

    ! Local variables.
    type(string), allocatable, dimension(:) :: lines

    ! The section markers: the line each section starts and ends on.
    integer :: section_start(num_sections), section_end(num_sections), s

    ! The entity table: each entity's dimension, tag and physical tag.
    integer, allocatable :: ent_dim(:), ent_tag(:), ent_phys(:)

    integer :: mesh_dim

    ! Load the mesh into memory.
    if (verbosity .ge. 1) write(*,'(a,a)') "Loading mesh file :", this % file % filename
    call this % file % read_lines(lines)

    if (verbosity .ge. 1) write(*,'(a)') "Identifying tags..."
    call find_tags(lines, section_start, section_end)

    if (verbosity .ge. 1) then
       do s = 1, num_sections
          write(*,*) sections(s), " : ", section_start(s), section_end(s)
       end do
    end if

    process_mesh_version: block

      integer                   :: num_tokens
      type(string), allocatable :: tokens(:)

      if (verbosity .ge. 1) write(*,'(a)') "Reading mesh information..."

      ! The first line of $MeshFormat contains the version number.
      associate(mlines => lines(section_start(MESH)+1:section_start(MESH)+1))
        call mlines(1) % tokenize(" ", num_tokens, tokens)
        if (floor(tokens(1) % as_real()) .ne. 4) then
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

      associate(tag_lines => lines(section_start(PHYSICAL_NAMES)+1:section_end(PHYSICAL_NAMES)-1))

        ! Set the intent(out) variable for the number of tags present.
        num_tags = tag_lines(1) % as_integer()

        ! Allocate space for the other two return variables.
        allocate(tag_info(num_tags))
        allocate(tag_numbers(num_tags), tag_physical_dimensions(num_tags), source=0)

        do iline = 1, num_tags

           ! Tokenize on the space delimiter.
           call tag_lines(iline+1) % tokenize(" ", num_tokens, tokens)

           ! The first token is the physical dimension; the second is the tag number.
           tag_physical_dimensions(iline) = tokens(1) % as_integer()
           tag_numbers(iline)             = tokens(2) % as_integer()

           !----------------------------------------------------------!
           ! The name is the quoted string on the line. Take everything
           ! between the first and last double quote so names with
           ! spaces (e.g. "outer boundary") are preserved.
           !----------------------------------------------------------!

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

    !----------------------------------------------------------------!
    ! Build the entity table from $Entities. An element takes its
    ! physical tag from the entity it belongs to (a 4.1 change). The
    ! header counts the points, curves, surfaces and volumes, and the
    ! lines follow in that order; a point line reads tag x y z numPhys
    ! phys..., the others tag bbox(6) numPhys phys..., so numPhys is
    ! token 5 of a point and token 8 otherwise.
    !----------------------------------------------------------------!

    process_entities: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: counts(0:3), d, i, ie, il, at

      if (verbosity .ge. 1) write(*,'(a)') "Reading entities..."

      il = section_start(ENTITIES) + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      counts = tokens(1:4) % as_integer()

      allocate(ent_dim(sum(counts)), ent_tag(sum(counts)), ent_phys(sum(counts)), source=0)

      ie = 0
      do d = 0, 3
         at = merge(5, 8, d == 0)
         do i = 1, counts(d)
            ie = ie + 1
            il = il + 1
            call lines(il) % tokenize(" ", num_tokens, tokens)
            ent_dim(ie) = d
            ent_tag(ie) = tokens(1) % as_integer()
            if (tokens(at) % as_integer() .gt. 0) ent_phys(ie) = tokens(at + 1) % as_integer()
         end do
      end do

      if (allocated(tokens)) deallocate(tokens)

      if (verbosity .ge. 1) write(*,'(a)') "Reading entities completed..."

    end block process_entities

    !----------------------------------------------------------------!
    ! $Nodes is block-structured: per block, all node tags then all
    ! coordinates. Node tags can be sparse and are mapped through
    ! findloc.
    !----------------------------------------------------------------!

    process_nodes: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: numblocks, il, ib, k, i, ivert

      if (verbosity .ge. 1) write(*,'(a)') "Reading nodes..."

      il = section_start(NODES) + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      numblocks    = tokens(1) % as_integer()
      num_vertices = tokens(2) % as_integer()

      allocate(vertex_numbers(num_vertices), vertex_tags(num_vertices), source=0)
      allocate(vertices(3, num_vertices), source=0.0_dp)

      il    = il + 1
      ivert = 0

      do ib = 1, numblocks

         ! The block header reads: entityDim entityTag parametric numNodesInBlock.
         call lines(il) % tokenize(" ", num_tokens, tokens)
         k  = tokens(4) % as_integer()
         il = il + 1

         !------------------------------------------------------------!
         ! The node tags are listed one integer per line, so they are
         ! read directly; tokenize requires a delimiter and a single
         ! value contains none.
         !------------------------------------------------------------!

         do i = 1, k
            read(lines(il) % str, *) vertex_numbers(ivert+i)
            il = il + 1
         end do

         ! Read the matching coordinates.
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

    !----------------------------------------------------------------!
    ! $Elements is block-structured. The top dimension contains the
    ! cells, one lower the faces, two lower the edges (so 2d and 3d
    ! are both handled). The physical tag of each element comes from
    ! its block's entity.
    !----------------------------------------------------------------!

    process_elements: block

      type(string), allocatable :: tokens(:)
      integer                   :: num_tokens
      integer                   :: numblocks, il0, il, ib, k, i
      integer                   :: bdim, btag, etype, nv, edim, phys
      integer                   :: filled(0:2)

      if (verbosity .ge. 1) write(*,'(a)') "Reading elements..."

      il = section_start(ELEMENTS) + 1
      call lines(il) % tokenize(" ", num_tokens, tokens)
      numblocks = tokens(1) % as_integer()
      il0       = section_start(ELEMENTS) + 2

      ! Pass A finds the mesh (top) dimension from the element types present.
      mesh_dim = 0
      il = il0
      do ib = 1, numblocks
         call lines(il) % tokenize(" ", num_tokens, tokens)
         etype    = tokens(3) % as_integer()
         k        = tokens(4) % as_integer()
         mesh_dim = max(mesh_dim, element_dimension(etype))
         il = il + 1 + k
      end do

      ! Pass B counts the cells, faces and edges by relative dimension.
      num_cells = 0
      num_faces = 0
      num_edges = 0
      il = il0
      do ib = 1, numblocks
         call lines(il) % tokenize(" ", num_tokens, tokens)
         etype = tokens(3) % as_integer()
         k     = tokens(4) % as_integer()
         edim  = element_dimension(etype)
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

      ! Allocate space for the cells, the faces and the edges, each as
      ! wide as the widest element of its dimension.
      allocate(cell_numbers(num_cells), num_cell_vertices(num_cells), cell_tags(num_cells), &
           & cell_types(num_cells), cell_vertices(widest_element(3), num_cells), source=0)
      allocate(face_numbers(num_faces), num_face_vertices(num_faces), face_tags(num_faces), &
           & face_types(num_faces), face_vertices(widest_element(2), num_faces), source=0)
      allocate(edge_numbers(num_edges), num_edge_vertices(num_edges), edge_tags(num_edges), &
           & edge_types(num_edges), edge_vertices(widest_element(1), num_edges), source=0)

      ! Pass C fills the arrays, taking the physical tag from the block's entity.
      filled = 0
      il = il0
      do ib = 1, numblocks

         ! The block header reads: entityDim entityTag elementType numElementsInBlock.
         call lines(il) % tokenize(" ", num_tokens, tokens)
         bdim  = tokens(1) % as_integer()
         btag  = tokens(2) % as_integer()
         etype = tokens(3) % as_integer()
         k     = tokens(4) % as_integer()

         edim  = element_dimension(etype)
         nv    = element_num_vertices(etype)
         phys  = entity_phys(bdim, btag)

         il = il + 1

         do i = 1, k

            ! An element line reads: elementTag node1 node2 ...
            call lines(il) % tokenize(" ", num_tokens, tokens)

            select case (mesh_dim - edim)
            case (0)
               filled(0) = filled(0) + 1
               call record(tokens, etype, phys, nv, filled(0), &
                    & cell_numbers, cell_types, cell_tags, num_cell_vertices, cell_vertices)
            case (1)
               filled(1) = filled(1) + 1
               call record(tokens, etype, phys, nv, filled(1), &
                    & face_numbers, face_types, face_tags, num_face_vertices, face_vertices)
            case (2)
               filled(2) = filled(2) + 1
               call record(tokens, etype, phys, nv, filled(2), &
                    & edge_numbers, edge_types, edge_tags, num_edge_vertices, edge_vertices)
            end select

            il = il + 1

         end do

      end do

      if (allocated(tokens)) deallocate(tokens)

      ! Warn when no face has a physical tag.
      if (count(face_tags .ne. 0) .eq. 0) then
         write(*,*) "no face has a physical tag - the mesh cannot be used for simulation"
      end if

      if (verbosity .ge. 1) write(*,'(a)') "Reading elements completed..."

    end block process_elements

    deallocate(lines)

  contains

    !==================================================================!
    ! Return the physical tag of the entity (edim, etag), or 0 if it
    ! belongs to no physical group.
    !==================================================================!

    pure integer function entity_phys(edim, etag)

      integer, intent(in) :: edim, etag

      integer :: at

      entity_phys = 0
      at = findloc(ent_dim == edim .and. ent_tag == etag, .true., dim=1)
      if (at >= 1) entity_phys = ent_phys(at)

    end function entity_phys

    !==================================================================!
    ! Record one element line - its number, then its nodes - at
    ! position `at` of one of the three element lists, the nodes
    ! mapped from file tags to positions in vertex_numbers.
    !==================================================================!

    pure subroutine record(tokens, etype, phys, nv, at, numbers, types, tags, counts, vertices)

      type(string), intent(in)    :: tokens(:)
      integer     , intent(in)    :: etype, phys, nv, at
      integer     , intent(inout) :: numbers(:), types(:), tags(:), counts(:), vertices(:,:)

      integer :: j

      numbers(at) = tokens(1) % as_integer()
      types(at)   = etype
      tags(at)    = phys
      counts(at)  = nv
      do j = 1, nv
         vertices(j, at) = findloc(vertex_numbers, tokens(1 + j) % as_integer(), dim=1)
      end do

    end subroutine record

  end subroutine mesh_data

  !====================================================================!
  ! Scan the file for the start and end line of each section the
  ! parser reads: $Name opens a section and $EndName closes it. A
  ! section absent from the file has both lines zero.
  !====================================================================!

  pure subroutine find_tags(lines, section_start, section_end)

    type(string), intent(in)  :: lines(:)
    integer     , intent(out) :: section_start(num_sections), section_end(num_sections)

    integer :: iline, s

    section_start = 0
    section_end   = 0

    do iline = 1, size(lines)
       do s = 1, num_sections
          if (index(lines(iline) % str, '$' // trim(sections(s))) .eq. 1) then
             section_start(s) = iline
          end if
          if (index(lines(iline) % str, '$End' // trim(sections(s))) .eq. 1) then
             section_end(s) = iline
          end if
       end do
    end do

  end subroutine find_tags

end module view_gmsh_loader
