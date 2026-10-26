!=====================================================================!
! Module that reads gmsh mesh files in the MSH 4.1 format (gmsh 4.x,
! e.g. 4.15). The format is documented here:
! http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
!
! The sections the parser reads, each opened by $Name and closed by
! $EndName:
!
!      $MeshFormat      version  file-type  data-size
!      $PhysicalNames   count, then  dimension  tag  "name"  per line
!      $Entities        numPoints numCurves numSurfaces numVolumes,
!                       then one line per entity:
!                         point     tag x y z numPhys phys ...
!                         other     tag bbox(6) numPhys phys ...
!      $Nodes           numBlocks numNodes minTag maxTag, then per block
!                         entityDim entityTag parametric numNodesInBlock
!                         one node tag per line
!                         one x y z per line
!      $Elements        numBlocks numElements minTag maxTag, then per block
!                         entityDim entityTag elementType numElementsInBlock
!                         elementTag node1 node2 ...  per line
!      $Periodic        numLinks, then per link
!                         entityDim entityTag entityTagMaster
!                         numAffine affine(16)
!                         numPairs, then  nodeTag nodeTagMaster  per line
!
! 4.1 differs from the legacy 2.2 in three ways the parser handles:
!
!   - $Entities: an element line no longer lists a physical tag; the
!     element belongs to a geometric entity that stores the tag. So the
!     per-element tag is a lookup:
!
!        element -> (entityDim, entityTag) -> physical tag
!
!   - $Nodes is block-structured: per entity block, all node tags then
!     all coordinates. Node tags may be sparse, and are mapped through a
!     table from number to position.
!   - $Elements is block-structured: per entity block, element lines
!     with no per-element tag list.
!
! Legacy 2.2 files are no longer read; regenerate the mesh with
! python meshgen/generate.py.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_gmsh_loader

  use util_precision    , only : dp
  use view_mesh_loader  , only : mesh_loader
  use util_file         , only : file
  use util_string       , only : string
  use view_mesh_geometry, only : element_dimension, &
       & element_num_vertices, widest_element
  use util_verbosity    , only : verbosity

  implicit none

  private

  public :: gmsh_loader

  ! The sections the parser reads, in the order they are numbered
  integer, parameter :: num_sections = 6
  integer, parameter :: MESH = 1, PHYSICAL_NAMES = 2, ENTITIES = 3, NODES = 4, ELEMENTS = 5, PERIODIC = 6
  character(len=13), parameter :: sections(num_sections) = &
       & [character(len=13) :: 'MeshFormat', 'PhysicalNames', 'Entities', 'Nodes', 'Elements', 'Periodic']

  ! The interface to construct a mesh_loader for gmsh
  interface gmsh_loader
     module procedure create
  end interface gmsh_loader

  !-------------------------------------------------------------------!
  ! Type that models the gmsh loader: the mesh file it reads
  !-------------------------------------------------------------------!

  type, extends(mesh_loader) :: gmsh_loader

     type(file) :: file ! The mesh file

   contains

     ! Implement the deferred procedures from the interface
     procedure :: mesh_data
     procedure :: periodic_links

  end type gmsh_loader

contains

  !===================================================================!
  ! Creates a loader of a .msh file. MSH 4.1 lines are long (entity
  ! bounding boxes, coordinate triples), well past the default
  ! 100-character buffer, so the file reads whole lines of up to 4096.
  !===================================================================!

  impure type(gmsh_loader) function create(filename) result (this)

    type(character(*)), intent(in) :: filename

    this % file = file(filename, 4096)

  end function create

  !===================================================================!
  ! Supplies all the information needed to create a mesh object. The
  ! elements are classified by their dimension relative to the
  ! highest present:
  !
  !      element dimension    is a        in 3d          in 2d
  !      the highest          cell        tetrahedron    triangle
  !      one lower            face        triangle       line
  !      two lower            edge        line           point
  !===================================================================!

  impure subroutine mesh_data(this, &
       & num_vertices, vertex_numbers, vertex_tags , vertices ,  &
       & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
       & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
       & cell_types  , face_types    , edge_types  , &
       & num_tags    , tag_numbers   , tag_physical_dimensions, tag_info )

    ! The arguments
    class(gmsh_loader), intent(in) :: this

    ! The vertices
    integer , intent(out)              :: num_vertices
    integer , intent(out), allocatable :: vertex_numbers(:)
    integer , intent(out), allocatable :: vertex_tags(:)
    real(dp), intent(out), allocatable :: vertices(:,:)

    ! The edges
    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_numbers(:)
    integer, intent(out), allocatable :: edge_tags(:)
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)
    integer, intent(out), allocatable :: edge_types(:)

    ! The faces
    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    integer, intent(out), allocatable :: face_types(:)

    ! The cells
    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)
    integer, intent(out), allocatable :: cell_types(:)

    ! These tag the boundaries and the domain with integers and strings
    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    integer     , allocatable, intent(out) :: tag_physical_dimensions(:)
    type(string), allocatable, intent(out) :: tag_info(:)

    ! The lines of the file
    type(string), allocatable, dimension(:) :: lines

    ! The section markers: the line each section starts and ends on
    integer :: section_start(num_sections), section_end(num_sections), s

    ! The entity table: each entity's dimension, tag and physical tag
    integer, allocatable :: ent_dim(:), ent_tag(:), ent_phys(:)

    integer :: mesh_dim

    ! The position in vertex_numbers of each node number
    integer, allocatable :: position(:)

    ! Load the mesh into memory
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

      ! The first line of $MeshFormat contains the version number
      associate(mlines => lines(section_start(MESH)+1:section_start(MESH)+1))
        call mlines(1) % tokenize(" ", num_tokens, tokens)
        if (floor(tokens(1) % as_real()) .ne. 4) then
           error stop "mesh format " // tokens(1) % str // &
                & " is not msh 4.x - regenerate with python meshgen/generate.py"
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

        ! Set the intent(out) variable for the number of tags present
        num_tags = tag_lines(1) % as_integer()

        ! Allocate space for the other two return variables
        allocate(tag_info(num_tags))
        allocate(tag_numbers(num_tags), tag_physical_dimensions(num_tags), source=0)

        do iline = 1, num_tags

           ! Tokenize on the space delimiter
           call tag_lines(iline+1) % tokenize(" ", num_tokens, tokens)

           ! The first token is the physical dimension; the second is
           ! the tag number
           tag_physical_dimensions(iline) = tokens(1) % as_integer()
           tag_numbers(iline)             = tokens(2) % as_integer()

           ! The name is the quoted string on the line. Take everything
           ! between the first and last double quote so names with
           ! spaces (e.g. "outer boundary") are preserved.
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

    !-----------------------------------------------------------------!
    ! Build the entity table from $Entities. An element takes its
    ! physical tag from the entity it belongs to (a 4.1 change). The
    ! header counts the points, curves, surfaces and volumes, and the
    ! lines follow in that order; a point line reads tag x y z numPhys
    ! phys..., the others tag bbox(6) numPhys phys..., so numPhys is
    ! token 5 of a point and token 8 otherwise:
    !
    !      point    tag  x  y  z  numPhys  phys
    !      token    1    2  3  4  5        6
    !      other    tag  xmin ymin zmin xmax ymax zmax  numPhys  phys
    !      token    1    2    3    4    5    6    7     8        9
    !
    ! An entity of no physical group has physical tag zero; of several,
    ! the first.
    !-----------------------------------------------------------------!

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

    !-----------------------------------------------------------------!
    ! $Nodes is block-structured: per block, all node tags then all
    ! coordinates. Node tags can be sparse and are mapped through a
    ! table from number to position, formed once.
    !-----------------------------------------------------------------!

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

         ! The block header reads: entityDim entityTag parametric
         ! numNodesInBlock
         call lines(il) % tokenize(" ", num_tokens, tokens)
         k  = tokens(4) % as_integer()
         il = il + 1

         ! The node tags are listed one integer per line, so they are
         ! read directly; tokenize requires a delimiter and a single
         ! value contains none
         do i = 1, k
            read(lines(il) % str, *) vertex_numbers(ivert+i)
            il = il + 1
         end do

         ! Read the matching coordinates
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

    position = position_of_number(vertex_numbers)

    !-----------------------------------------------------------------!
    ! $Elements is block-structured. The top dimension contains the
    ! cells, one lower the faces, two lower the edges (so 2d and 3d are
    ! both handled). The physical tag of each element comes from its
    ! block's entity. Three passes over the blocks:
    !
    !      A    the mesh's top dimension, from the element types present
    !      B    the cells, faces and edges counted, and their arrays
    !           allocated as wide as the widest element of each
    !      C    each element recorded, its nodes mapped to positions
    !-----------------------------------------------------------------!

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

      ! Pass A finds the mesh (top) dimension from the element types
      ! present
      mesh_dim = 0
      il = il0
      do ib = 1, numblocks
         call lines(il) % tokenize(" ", num_tokens, tokens)
         etype    = tokens(3) % as_integer()
         k        = tokens(4) % as_integer()
         mesh_dim = max(mesh_dim, element_dimension(etype))
         il = il + 1 + k
      end do

      ! Pass B counts the cells, faces and edges by relative dimension
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
      ! wide as the widest element of its dimension
      allocate(cell_numbers(num_cells), num_cell_vertices(num_cells), cell_tags(num_cells), &
           & cell_types(num_cells), cell_vertices(widest_element(3), num_cells), source=0)
      allocate(face_numbers(num_faces), num_face_vertices(num_faces), face_tags(num_faces), &
           & face_types(num_faces), face_vertices(widest_element(2), num_faces), source=0)
      allocate(edge_numbers(num_edges), num_edge_vertices(num_edges), edge_tags(num_edges), &
           & edge_types(num_edges), edge_vertices(widest_element(1), num_edges), source=0)

      ! Pass C fills the arrays, taking the physical tag from the
      ! block's entity
      filled = 0
      il = il0
      do ib = 1, numblocks

         ! The block header reads: entityDim entityTag elementType
         ! numElementsInBlock
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

      ! Warn when no face has a physical tag
      if (count(face_tags .ne. 0) .eq. 0) then
         write(*,*) "no face has a physical tag - the mesh cannot be used for simulation"
      end if

      if (verbosity .ge. 1) write(*,'(a)') "Reading elements completed..."

    end block process_elements

    deallocate(lines)

  contains

    !-----------------------------------------------------------------!
    ! Returns the physical tag of the entity (edim, etag), or 0 if it
    ! belongs to no physical group
    !-----------------------------------------------------------------!

    pure integer function entity_phys(edim, etag)

      integer, intent(in) :: edim, etag

      integer :: at

      entity_phys = 0
      at = findloc(ent_dim == edim .and. ent_tag == etag, .true., dim=1)
      if (at >= 1) entity_phys = ent_phys(at)

    end function entity_phys

    !-----------------------------------------------------------------!
    ! Records one element line, its number and then its nodes, at
    ! position at of one of the three element lists, the nodes mapped
    ! from file tags to positions in vertex_numbers
    !-----------------------------------------------------------------!

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
         vertices(j, at) = located(position, tokens(1 + j) % as_integer())
      end do

    end subroutine record

  end subroutine mesh_data

  !===================================================================!
  ! Returns the periodic links of the file's $Periodic section, when
  ! present: each link identifies one entity with a master entity under
  ! a translation, master + translation = image, and lists the vertex
  ! pairs (image, master) it identifies. The vertices are positions in
  ! vertex_numbers, the array mesh_data returns. A link whose affine
  ! transformation is not a translation is invalid input: the mesh then
  ! identifies entities under a rotation, which the periodic face shift
  ! does not represent.
  !
  !      first(l), last(l)   the pairs of link l, positions in image/master
  !      translation(:, l)   the three components of the translation
  !
  ! The affine matrix is stated by rows, so the translation is its last
  ! column, entries 4, 8 and 12; the rest must be the identity:
  !
  !      | 1  0  0  t_x |
  !      | 0  1  0  t_y |
  !      | 0  0  1  t_z |
  !      | 0  0  0  1   |
  !===================================================================!

  impure subroutine periodic_links(this, vertex_numbers, first, last, image, master, translation)

    class(gmsh_loader)   , intent(in)  :: this
    integer              , intent(in)  :: vertex_numbers(:)
    integer , allocatable, intent(out) :: first(:), last(:), image(:), master(:)
    real(dp), allocatable, intent(out) :: translation(:,:)

    type(string), allocatable :: lines(:), tokens(:)
    integer                   :: section_start(num_sections), section_end(num_sections)
    integer                   :: num_links, num_tokens, il, l, num_affine, num_pairs, i, at, tag_image, tag_master
    real(dp)    , allocatable :: affine(:)
    integer     , allocatable :: position(:)
    real(dp)    , parameter   :: identity(16) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         &                                       0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
         &                                       0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         &                                       0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    character(len=250)        :: message

    call this % file % read_lines(lines)
    call find_tags(lines, section_start, section_end)

    if (section_start(PERIODIC) == 0) then
       allocate(first(0), last(0), image(0), master(0), translation(3, 0))
       return
    end if

    position = position_of_number(vertex_numbers)
    il = section_start(PERIODIC) + 1
    read(lines(il) % str, *) num_links
    allocate(first(num_links), last(num_links), translation(3, num_links))

    ! One pass counts the pairs, the second records them
    num_pairs = 0
    il = il + 1
    do l = 1, num_links
       il = il + 1
       read(lines(il) % str, *) num_affine
       il = il + 1
       read(lines(il) % str, *) i
       num_pairs = num_pairs + i
       il = il + 1 + i
    end do
    allocate(image(num_pairs), master(num_pairs))

    at = 0
    il = section_start(PERIODIC) + 2
    do l = 1, num_links

       ! The link's line reads: entityDim entityTag entityTagMaster
       call lines(il) % tokenize(" ", num_tokens, tokens)
       tag_image  = tokens(2) % as_integer()
       tag_master = tokens(3) % as_integer()
       il = il + 1

       ! The affine line reads: numAffine, then the 4 x 4 matrix by rows
       call lines(il) % tokenize(" ", num_tokens, tokens)
       num_affine = tokens(1) % as_integer()
       if (num_affine /= 16) then
          write(message,'(a,i0,a,i0,a,i0)') 'view_gmsh_loader: a periodic link is a translation, &
               &stated as a 4 x 4 affine matrix; entity ', tag_image, ' onto ', tag_master, &
               & ' has affine values = ', num_affine
          error stop trim(message)
       end if
       affine = tokens(2:17) % as_real()
       translation(:, l) = affine([4, 8, 12])
       affine([4, 8, 12]) = 0.0_dp
       if (any(abs(affine - identity) > 1.0e-12_dp)) then
          write(message,'(a,i0,a,i0)') 'view_gmsh_loader: a periodic link is a translation; entity ', &
               & tag_image, ' onto ', tag_master, ' is identified under a rotation or a scaling'
          error stop trim(message)
       end if
       il = il + 1

       ! The pairs: nodeTag nodeTagMaster, one per line
       read(lines(il) % str, *) num_pairs
       il = il + 1
       first(l) = at + 1
       do i = 1, num_pairs
          call lines(il) % tokenize(" ", num_tokens, tokens)
          at = at + 1
          image(at)  = located(position, tokens(1) % as_integer())
          master(at) = located(position, tokens(2) % as_integer())
          if (image(at) == 0 .or. master(at) == 0) then
             error stop 'view_gmsh_loader: a periodic pair names a node absent from $Nodes'
          end if
          il = il + 1
       end do
       last(l) = at

    end do

  end subroutine periodic_links

  !===================================================================!
  ! Returns the start and end line of each section the parser reads:
  ! $Name opens a section and $EndName closes it. A section absent from
  ! the file has both lines zero.
  !===================================================================!

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

  !===================================================================!
  ! Returns the position of each node number in the list of numbers,
  ! the first where a number repeats: a table indexed by number, so
  ! that reading an element's nodes is linear in the elements.
  !
  !      numbers           [ 7  3  12 ]
  !      position(1:12)    [ 0 0 2 0 0 0 1 0 0 0 0 3 ]
  !===================================================================!

  pure function position_of_number(numbers) result(position)
    integer, intent(in)  :: numbers(:)
    integer, allocatable :: position(:)
    integer              :: i
    allocate(position(max(0, maxval(numbers))), source=0)
    do i = size(numbers), 1, -1
       if (numbers(i) >= 1) position(numbers(i)) = i
    end do
  end function position_of_number

  !===================================================================!
  ! Returns a number's position in the table, zero for a number the
  ! list does not contain
  !===================================================================!

  pure integer function located(position, number)
    integer, intent(in) :: position(:), number
    located = 0
    if (number >= 1 .and. number <= size(position)) located = position(number)
  end function located

end module view_gmsh_loader
