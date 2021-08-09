!=====================================================================!
! Unstructured mesh handler.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_mesh

  use iso_fortran_env       , only : dp => REAL64, error_unit
  use interface_mesh_loader , only : mesh_loader
  use class_string          , only : string
  use module_mesh_utils     , only : find, distance, elem_type_face_count, &
       & cross_product, form_cell_faces, transpose_connectivities, &
       & order_face_vertices, sparse_transpose_matmul

  implicit none

  private
  public :: mesh

  ! Constructor
  interface mesh
     module procedure create_mesh
  end interface mesh

  type :: face_group

     ! Faces
     type(string) :: group_name
     integer :: group_number
     integer :: group_num_faces
     integer, allocatable :: group_face_numbers(:)
     integer, allocatable :: group_face_vertices(:,:)
     integer, allocatable :: group_num_face_vertices(:)
     integer, allocatable :: group_face_types(:)

   contains

     procedure :: print
     final     :: destroy_face_group

  end type face_group

  ! constructor function for face_group
  interface face_group
     module procedure create_face_group
  end interface face_group

  !-------------------------------------------------------------------!
  ! Mesh datatype. A collection of vertices, cells and faces.
  !-------------------------------------------------------------------!

  type :: mesh ! rename as topology?

     integer :: max_print = 20
     logical :: initialized = .false.

     ! Based on PhysicalNames?
     integer                   :: num_tags
     integer     , allocatable :: tag_numbers(:)
     integer     , allocatable :: tag_physical_dimensions(:)
     type(string), allocatable :: tag_info(:)

     !================================================================!
     ! Basic Topology information
     !================================================================!

     ! Fundamental vertex info
     ! type(vertex_group), allocatable :: vertex_groups(:)

     integer :: num_vertices
     real(dp) , allocatable :: vertices(:,:)     ! [[x,y,z],1:num_vertices]
     integer  , allocatable :: vertex_numbers(:) ! global
     integer  , allocatable :: vertex_tags(:)    ! not set

     ! Fundamental face info
     ! type(edge_group), allocatable :: face_groups(:)

     integer :: num_edges
     integer  , allocatable :: edge_numbers(:)
     integer  , allocatable :: edge_tags(:)
     integer  , allocatable :: edge_types(:)
     integer  , allocatable :: edge_vertices(:,:)        ! [[v1,v2],1:nedges]
     integer  , allocatable :: num_edge_vertices(:)

     ! Fundamental face info
     type(face_group), allocatable :: face_groups(:)

     integer :: num_faces
     integer  , allocatable :: face_numbers(:)
     integer  , allocatable :: face_tags(:)            ! one face can be a part of two tags?
     integer  , allocatable :: face_types(:)
     integer  , allocatable :: face_vertices(:,:)        ! [[v1,v2],1:nfaces]
     integer  , allocatable :: num_face_vertices(:)

     ! Fundamental cell info
     integer :: num_cells
     integer  , allocatable :: cell_numbers(:)
     integer  , allocatable :: cell_tags(:)
     integer  , allocatable :: cell_types(:)
     integer  , allocatable :: cell_vertices(:,:)   ! [[v1,v2,v3], 1:ncells]
     integer  , allocatable :: num_cell_vertices(:) ! [1:ncells]

     !================================================================!
     ! Derived Topology information
     !================================================================!

     ! Inverse cell information
     integer  , allocatable :: vertex_cells(:,:)    ! [[c1,c2,c3],[1:nvertices]]
     integer  , allocatable :: num_vertex_cells(:)  ! [1:nvertices]

     ! Inverse face information
     integer  , allocatable :: vertex_faces(:,:)    ! [[f1,f2,f3],[1:nfaces]]
     integer  , allocatable :: num_vertex_faces(:)  ! [1:nfaces]

     ! Inverse edge information
     integer  , allocatable :: vertex_edges(:,:)    ! [[e1,e2,e3],[1:nedges]]
     integer  , allocatable :: num_vertex_edges(:)  ! [1:nedges]

     ! Intermidiate connectivities and their inverse
     integer  , allocatable :: num_cell_faces(:)         ! [1:ncells]
     integer  , allocatable :: cell_faces(:,:)           ! [[f1,f2,f3..],1:ncells]
     integer  , allocatable :: cell_faces_type(:,:)      ! [[tri,qua,tri..],1:ncells]

     integer  , allocatable :: num_face_cells(:)         ! [1:nfaces]
     integer  , allocatable :: face_cells(:,:)           ! [[c1,c2...],1:nfaces]
     integer  , allocatable :: face_cells_type(:,:)      ! [[hex,tet,..],1:nfaces]

     ! Intermidiate connectivities and their inverse
     integer  , allocatable :: num_face_edges(:)         ! [1:nfaces]
     integer  , allocatable :: face_edges(:,:)           ! [[e1,e2,e3..],1:nfaces]
     integer  , allocatable :: num_edge_faces(:)         ! [1:nedges]
     integer  , allocatable :: edge_faces(:,:)           ! [[f1,f2...],1:nedges]

     !================================================================!
     ! Derived Geometry information
     !================================================================!

     ! Derived cell info
     real(dp) , allocatable :: cell_centers(:,:)         ! [[x,y,z] , 1:ncells]
     real(dp) , allocatable :: cell_volumes(:)           ! [1:ncells]

     ! Derived vertex info
     real(dp) , allocatable :: face_centers(:,:)         ! [[x,y,z],1:nfaces]
     real(dp) , allocatable :: face_areas(:)             ! [1:nfaces]
     real(dp) , allocatable :: face_deltas(:)             ! [1:nfaces]
     real(dp) , allocatable :: lvec(:,:)                 ! [[lx,ly,lz],1:nfaces]

     ! maybe unify t1, t2 and n as cell_face_triad ?  also assuming
     ! normal flux is always predominant can stall convergence.
     real(dp) , allocatable :: cell_face_tangents(:,:,:,:) ! [[tx,ty,tz], t1 or t2 , [f1,f2,f3..] 1:ncells]
     real(dp) , allocatable :: cell_face_normals(:,:,:)  ! [[nx,ny,nz], [f1,f2,f3..] 1:ncells]

     real(dp) , allocatable :: vertex_cell_weights(:,:)  ! [[wc1,wc2...],1:vertices]
     real(dp) , allocatable :: face_cell_weights(:,:)    ! [[wc1,wc2],1:nfaces]

   contains

     ! Type bound procedures
     procedure :: to_string
     procedure :: initialize

     ! Evaluation routines
     procedure :: evaluate_cell_centers
     procedure :: evaluate_cell_volumes
     procedure :: evaluate_face_centers_areas
     procedure :: evaluate_face_tangents_normals
     procedure :: evaluate_centroidal_vector
     procedure :: evaluate_face_deltas
     procedure :: evaluate_face_weight
     procedure :: evaluate_vertex_weight

     procedure :: evaluate_face_centers_areas_2d
     procedure :: evaluate_face_tangents_normals_2d

     ! Mesh info
     procedure :: get_num_elems
     procedure :: create_face_groups
     procedure :: invert_connectivities

     ! Destructor
     final   :: destroy

  end type mesh

contains

  impure type(face_group) function create_face_group( &
       & group_name, group_number, &
       & num_faces, face_numbers, face_vertices, &
       & num_face_vertices, face_types ) &
       & result(this)

    type(string) , intent(in) :: group_name
    integer      , intent(in) :: group_number
    integer      , intent(in) :: num_faces
    integer      , intent(in) :: face_numbers(:)
    integer      , intent(in) :: face_vertices(:,:)
    integer      , intent(in) :: num_face_vertices(:)
    integer      , intent(in) :: face_types(:)

    this % group_name      = group_name
    this % group_number    = group_number
    this % group_num_faces = num_faces

    allocate(this % group_face_numbers     , source = face_numbers)
    allocate(this % group_face_vertices    , source = face_vertices)
    allocate(this % group_num_face_vertices, source = num_face_vertices)
    allocate(this % group_face_types       , source = face_types)

    ! sanity checks
    if (this % group_num_faces .ne. size(face_numbers)) &
         & error stop "face group creation failed: inconsistent face numbers"
    if (this % group_num_faces .ne. size(face_vertices, dim=2)) &
         & error stop "face group creation failed: inconsistent face vertices"
    if (this % group_num_faces .ne. size(face_types)) &
         & error stop "face group creation failed : inconsistent face types"
    if (this % group_num_faces .ne. size(num_face_vertices)) &
         & error stop "face group creation failed : inconsistent # face vertices"

  end function create_face_group

  pure subroutine destroy_face_group(this)

    type(face_group), intent(inout) :: this

    if (allocated(this % group_name % str))        deallocate(this % group_name % str)
    if (allocated(this % group_face_numbers))      deallocate(this % group_face_numbers)
    if (allocated(this % group_face_vertices))     deallocate(this % group_face_vertices)
    if (allocated(this % group_num_face_vertices)) deallocate(this % group_num_face_vertices)
    if (allocated(this % group_face_types))        deallocate(this % group_face_types)

  end subroutine destroy_face_group

  impure elemental subroutine print(this)

    class(face_group), intent(in) :: this

    print *, "face group number : ", this % group_number
    print *, "face group number of faces : ", this % group_num_faces
    print *, "face group name : ", this % group_name % str

  end subroutine print

  !================================================================!
  ! Constructor for mesh object using mesh loader
  !================================================================!

  type(mesh) function create_mesh(loader) result(me)

    class(mesh_loader), intent(in)  :: loader
    integer           , allocatable :: bface_numbers(:)
    integer           , allocatable :: bface_tags(:)
    integer           , allocatable :: bface_types(:)
    integer           , allocatable :: bface_vertices(:,:)
    integer           , allocatable :: bnum_face_vertices(:)
    integer                         :: bnum_faces

    !-----------------------------------------------------------------!
    ! Get the fundamental information needed
    !-----------------------------------------------------------------!

    call loader % get_mesh_data( &
         & me % num_vertices, me % vertex_numbers, me % vertex_tags , me % vertices ,  &
         & me % num_edges   , me % edge_numbers  , me % edge_tags   , me % edge_vertices , me % num_edge_vertices , &
         & bnum_faces       , bface_numbers      , bface_tags       , bface_vertices     , bnum_face_vertices , &
         & me % num_cells   , me % cell_numbers  , me % cell_tags   , me % cell_vertices , me % num_cell_vertices , &
         & me % cell_types  , bface_types       , me % edge_types  , &
         & me % num_tags    , me % tag_numbers   , me % tag_physical_dimensions, me % tag_info )

    ! determine the number of faces algebraically
    me % num_faces = (sum(elem_type_face_count(me % cell_types)) - bnum_faces)/2 + bnum_faces

    ! allocate the mesh's face data
    allocate(me % face_numbers(me % num_faces))
    allocate(me % face_tags(me % num_faces))
    allocate(me % face_types(me % num_faces))
    allocate(me % face_vertices(4, me % num_faces)) ! max is quadrilateral face
    allocate(me % num_face_vertices(me % num_faces))

    me % face_numbers = 0
    me % face_tags = 0
    me % face_types = 0
    me % face_vertices = 0
    me % num_face_vertices = 0

    ! add the remaining boundary faces to the array
    add_boundary_faces : block

      integer :: iface

      do iface = 1, bnum_faces

         me % face_numbers(iface)      = iface
         me % face_tags(iface)         = bface_tags(iface)
         me % face_types(iface)        = bface_types(iface)
         me % num_face_vertices(iface) = bnum_face_vertices(iface)
         me % face_vertices(1:bnum_face_vertices(iface),iface) = bface_vertices(1:bnum_face_vertices(iface), iface)

      end do

    end block add_boundary_faces

    form_internal_faces: block

      integer   :: num_faces
      integer   :: icell, jcell, ivertex
      integer   :: shared_node_count
      integer   :: shared_vertices(4)

      ! This logic will naturally avoid forming the boundary faces,
      ! because we are dotting only with a different interior cell.
      ! Therefore we can simply append these new internal faces to the
      ! existing (boundary) faces provided by the mesh

      ! may need to adjust a bit for partitioned mesh (may be ignore ghost cells?)
      ! face numbers need to be global in assignment?

      ! counting can be challenging for higher order mesh (uff)

      ! can apply the same logic for CP x PC?

      ! becareful about :, counts and 0

      num_faces = bnum_faces ! me % num_boundary_faces

      do icell = 1, me % num_cells

         shared_vertices = 0

         ! pick only upper triangle to avoid double counting
         do jcell = icell + 1, me % num_cells

            ! sparsely carry out dot product between two rows to yield
            ! the count of shared vertices
            shared_node_count = 0
            do ivertex = 1, me % num_cell_vertices(icell)
               if (any (me % cell_vertices(1 : me % num_cell_vertices(jcell),jcell) .eq. &
                    & me % cell_vertices(ivertex, icell) &
                    & ) .eqv. .true.) then
                  shared_node_count = shared_node_count + 1
                  shared_vertices(shared_node_count) = me % cell_vertices(ivertex, icell)
               end if
            end do

            ! post-process: how is icell related to jcell?
            select case (shared_node_count)
            case(0)
               ! no vertices are shared between the two cells
            case (1:2)
               ! a vertex or an edge is shared between the two cells
            case (3)
               ! a 3-noded triangle is shared between the two cells
               num_faces = num_faces + 1
               me % face_numbers(num_faces)      = num_faces
               me % face_tags(num_faces)         = me % cell_tags(icell)
               me % face_types(num_faces)        = 2
               me % num_face_vertices(num_faces) = 3
               call order_face_vertices(me % cell_types(jcell), me % cell_vertices(:,jcell), shared_vertices(1:3))
               me % face_vertices(1:3,num_faces) = shared_vertices(1:3)
            case (4)
               ! a 4-node quadrangle is shared between the two cells
               num_faces = num_faces + 1
               me % face_numbers(num_faces)      = num_faces
               me % face_tags(num_faces)         = me % cell_tags(icell)
               me % face_types(num_faces)        = 3
               me % num_face_vertices(num_faces) = 4
               call order_face_vertices(me % cell_types(jcell), me % cell_vertices(:,jcell), shared_vertices(1:4))
               me % face_vertices(1:4,num_faces) = shared_vertices(1:4)
            case default
               print *, shared_node_count
               error stop "confused in shared node counting forming internal faces"
            end select

         end do

      end do

      if (num_faces .ne. me % num_faces) then
         print *, "numfaces don't match", num_faces, me % num_faces
      end if

      if (maxval(me % face_vertices) .ne. me % num_vertices) then
         print *, "check face vertices"
         error stop
      end if


    end block form_internal_faces

    ! Sanity check (make sure numbering is continuous), although it
    ! may not start from one. Applicable only for unpartitioned mesh
!!$    if (num_images() .eq. 1) then
!!$       if (me % num_vertices .gt. 0 .and. &
!!$            & maxval(me % vertex_numbers) -  minval(me % vertex_numbers) + 1 .ne. me % num_vertices) &
!!$            & error stop
!!$       if (me % num_edges    .gt. 0 .and. &
!!$            & maxval(me % edge_numbers  ) -  minval(me % edge_numbers  ) + 1 .ne. me % num_edges   ) &
!!$            & error stop
!!$       if (me % num_faces    .gt. 0 .and. &
!!$            & maxval(me % face_numbers  ) -  minval(me % face_numbers  ) + 1 .ne. me % num_faces   ) &
!!$            & error stop
!!$       if (me % num_cells    .gt. 0 .and. &
!!$            & maxval(me % cell_numbers  ) -  minval(me % cell_numbers  ) + 1 .ne. me % num_cells   ) &
!!$            & error stop
!!$    end if

    ! Perform initialization tasks and store the resulting flag
    me % initialized = me % initialize()
    if (me % initialized .eqv. .false.) then
       write(error_unit,*) "Mesh.Construct: failed"
       error stop
    end if

    ! call me % create_face_groups()

  end function create_mesh

  impure elemental subroutine create_face_groups(this)

    class(mesh), intent(inout) :: this

    integer, allocatable :: face_tag_numbers(:)
    logical, allocatable :: group_mask(:)

    integer :: ifacegroup, face_group_number
    integer :: num_node_tags, num_edge_tags, num_face_tags, num_cell_tags

    num_node_tags = count(this % tag_physical_dimensions .eq. 0)
    num_edge_tags = count(this % tag_physical_dimensions .eq. 1)
    num_face_tags = count(this % tag_physical_dimensions .eq. 2)
    num_cell_tags = count(this % tag_physical_dimensions .eq. 3)

    write (*, '(1x, a, 1x, i0)') "number of node tags :", num_node_tags
    write (*, '(1x, a, 1x, i0)') "number of edge tags :", num_edge_tags
    write (*, '(1x, a, 1x, i0)') "number of face tags :", num_face_tags
    write (*, '(1x, a, 1x, i0)') "number of cell tags :", num_cell_tags

    allocate(face_tag_numbers(num_face_tags))
    face_tag_numbers = pack(this % tag_numbers, mask = this % tag_physical_dimensions .eq. 2)

    allocate(this % face_groups(num_face_tags))

    ! create face groups
    loop_face_groups: do ifacegroup = 1, num_face_tags

       face_group_number = face_tag_numbers(ifacegroup)

       ! allocate group mask
       allocate(group_mask, source = this % face_tags .eq. face_group_number)

       associate(&
            & gname              => this % tag_info(face_group_number), &
            & gnumber            => face_group_number, &
            & gnum_faces         => count(group_mask), &
            & gface_numbers      => pack(this % face_numbers, mask = group_mask), &
            & gface_vertices     => pack(this % face_vertices, mask = spread(group_mask, 1, size(this % face_vertices, 1))), &
            & gnum_face_vertices => pack(this % num_face_vertices, mask = group_mask), &
            & gface_types        => pack(this % face_types, mask = group_mask) &
            & )

         ! todo: may not need to use the first dimension of
         ! 'face_vertices', instead extract from 'face_ypes' in the
         ! group, and the corresponding # vertices

         this % face_groups(ifacegroup) = face_group(gname, gnumber, gnum_faces, &
              & gface_numbers, &
              & reshape(gface_vertices, [size(this % face_vertices, 1), gnum_faces]), &
              & gnum_face_vertices, gface_types)

       end associate

     ! deallocate group mask
     if (allocated(group_mask)) deallocate(group_mask)

    end do loop_face_groups

    if (allocated(face_tag_numbers)) deallocate(face_tag_numbers)

  end subroutine create_face_groups

  !===================================================================!
  ! Destructor for file object
  !===================================================================!

  pure subroutine destroy(this)

    type(mesh), intent(inout) :: this

    if (allocated(this % tag_numbers)) deallocate(this % tag_numbers)
    if (allocated(this % tag_physical_dimensions)) deallocate(this % tag_physical_dimensions)
    if (allocated(this % tag_info)) deallocate(this % tag_info)

    if (allocated(this % vertices)) deallocate(this % vertices)
    if (allocated(this % vertex_numbers)) deallocate(this % vertex_numbers)

    if (allocated(this % vertex_tags)) deallocate(this % vertex_tags)

    if (allocated(this % edge_numbers)) deallocate(this % edge_numbers)
    if (allocated(this % edge_tags)) deallocate(this % edge_tags)
    if (allocated(this % edge_vertices)) deallocate(this % edge_vertices)
    if (allocated(this % num_edge_vertices)) deallocate( this % num_edge_vertices)

    if (allocated(this % face_numbers)) deallocate(this % face_numbers)
    if (allocated(this % face_tags)) deallocate(this % face_tags)
    if (allocated(this % face_vertices)) deallocate(this % face_vertices)
    if (allocated(this % num_face_vertices)) deallocate(this % num_face_vertices)

    if (allocated(this % cell_numbers)) deallocate(this % cell_numbers)
    if (allocated(this % cell_tags)) deallocate(this % cell_tags)
    if (allocated(this % cell_vertices)) deallocate(this % cell_vertices)
    if (allocated(this % num_cell_vertices)) deallocate(this % num_cell_vertices)

    if (allocated(this % vertex_cells)) deallocate(this % vertex_cells)
    if (allocated(this % num_vertex_cells)) deallocate(this % num_vertex_cells)

    if (allocated(this % vertex_faces)) deallocate(this % vertex_faces)
    if (allocated(this % num_vertex_faces)) deallocate(this % num_vertex_faces)

    if (allocated(this % vertex_edges)) deallocate(this % vertex_edges)
    if (allocated(this % num_vertex_edges)) deallocate(this % num_vertex_edges)

    if (allocated(this % num_cell_faces)) deallocate(this % num_cell_faces)
    if (allocated(this % cell_faces)) deallocate(this % cell_faces)
    if (allocated(this % num_face_cells)) deallocate(this % num_face_cells)
    if (allocated(this % face_cells)) deallocate(this % face_cells)

    if (allocated(this % num_face_edges)) deallocate(this % num_face_edges)
    if (allocated(this % face_edges)) deallocate(this % face_edges)
    if (allocated(this % num_edge_faces)) deallocate(this % num_edge_faces)
    if (allocated(this % edge_faces)) deallocate(this % edge_faces)

    if (allocated(this % cell_centers)) deallocate(this % cell_centers)
    if (allocated(this % cell_volumes)) deallocate(this % cell_volumes)

    if (allocated(this % face_centers)) deallocate(this % face_centers)
    if (allocated(this % face_areas)) deallocate(this % face_areas)
    if (allocated(this % face_deltas)) deallocate(this % face_deltas)
    if (allocated(this % lvec)) deallocate(this % lvec)

    if (allocated(this % cell_face_tangents)) deallocate(this % cell_face_tangents)
    if (allocated(this % cell_face_normals)) deallocate(this % cell_face_normals)
    if (allocated(this % vertex_cell_weights)) deallocate(this % vertex_cell_weights)
    if (allocated(this % face_cell_weights)) deallocate(this % face_cell_weights)

  end subroutine destroy

  subroutine invert_connectivities(this)

    class(mesh), intent(inout) :: this

    !-----------------------------------------------------------------!
    ! Find VertexCell conn. by inverting CellVertex conn.
    !-----------------------------------------------------------------!

    vertex_cell: block

      integer :: ivertex

      write(*,'(a)') "Inverting CellVertex Map..."

      call transpose_connectivities( &
           & this % cell_vertices, &
           & this % num_cell_vertices, &
           & this % vertex_cells, &
           & this % num_vertex_cells)

      write(*,'(a)') "Inverting CellVertex Map complete ..."

      if (allocated(this % vertex_cells) .and. size(this % vertex_cells, dim = 2) .gt. 0) then

         write(*,'(a,i8,a,i8)') &
              & "Vertex to cell info for", min(this % max_print,this % num_vertices), &
              & " vertices out of ", this % num_vertices

         do ivertex = 1, min(this % max_print,this % num_vertices)
            write(*,*) &
                 & 'vertex [', this % vertex_numbers(ivertex), ']', &
                 & 'num cells [', this % num_vertex_cells(ivertex), ']',&
                 & 'cells [', this % vertex_cells(1:this % num_vertex_cells(ivertex),ivertex), ']'
         end do

         ! Sanity check
         if (minval(this % num_vertex_cells) .lt. 1) then
            write(error_unit, *) 'Error: There are vertices not mapped to a cell'
            error stop
         end if

      else

         write(*,'(a)') "Vertex to cell info not computed"

      end if

    end block vertex_cell

    !-----------------------------------------------------------------!
    ! Find VertexFace conn. by inverting FaceVertex conn.
    !-----------------------------------------------------------------!
!!$
!!$    vertex_face: block
!!$
!!$      integer :: ivertex
!!$
!!$      write(*,'(a)') "Inverting FaceVertex Map..."
!!$
!!$      call transpose_connectivities( &
!!$           & this % face_vertices, &
!!$           & this % num_face_vertices, &
!!$           & this % vertex_faces, &
!!$           & this % num_vertex_faces)
!!$
!!$      write(*,'(a)') "Inverting FaceVertex Map complete..."
!!$
!!$      if (allocated(this % vertex_faces) .and. size(this % vertex_faces, dim = 2) .gt. 0) then
!!$
!!$         write(*,'(a,i8,a,i8)') &
!!$              & "Vertex to face info for", min(this % max_print,this % num_vertices), &
!!$              & " vertices out of ", this % num_vertices
!!$
!!$         do ivertex = 1, min(this % max_print,this % num_vertices)
!!$            write(*,*) &
!!$                 & 'vertex [', this % vertex_numbers(ivertex), ']',&
!!$                 & 'num_vertex_faces [', this % num_vertex_faces(ivertex), ']',&
!!$                 & 'faces [', this % vertex_faces( &
!!$                 & 1:this % num_vertex_faces(ivertex),ivertex), ']'
!!$         end do
!!$
!!$         ! Sanity check
!!$         if (minval(this % num_vertex_faces) .lt. 1) then
!!$            write(error_unit, *) 'Error: There are vertices not mapped to a face'
!!$            error stop
!!$         end if
!!$
!!$      else
!!$
!!$         write(*,'(a)') "Vertex to face info not computed"
!!$
!!$      end if
!!$
!!$    end block vertex_face

    !-----------------------------------------------------------------!
    ! Find VertexEdge conn. by inverting EdgeVertex conn.
    !-----------------------------------------------------------------!

!!$    vertex_edge: block
!!$
!!$      integer :: ivertex
!!$
!!$      write(*,'(a)') "Inverting EdgeVertex Map..."
!!$
!!$      call transpose_connectivities( &
!!$           & this % edge_vertices, &
!!$           & this % num_edge_vertices, &
!!$           & this % vertex_edges, &
!!$           & this % num_vertex_edges)
!!$
!!$      write(*,'(a)') "Inverting EdgeVertex Map complete..."
!!$
!!$      if (allocated(this % vertex_edges) .and. size(this % vertex_edges, dim = 2) .gt. 0) then
!!$
!!$         write(*,'(a,i8,a,i8)') &
!!$              & "Vertex to edge info for", min(this % max_print,this % num_vertices), &
!!$              & " vertices out of ", this % num_vertices
!!$
!!$         do ivertex = 1, min(this % max_print,this % num_vertices)
!!$            write(*,*) &
!!$                 & 'vertex [', this % vertex_numbers(ivertex), ']',&
!!$                 & 'num_vertex_edges [', this % num_vertex_edges(ivertex) , ']',&
!!$                 & 'edges [', this % vertex_edges(&
!!$                 & 1:this % num_vertex_edges(ivertex),ivertex), ']'
!!$         end do
!!$
!!$         ! Sanity check
!!$         if (minval(this % num_vertex_edges) .lt. 1) then
!!$            write(error_unit, *) 'Error: There are vertices not mapped to a edge'
!!$            error stop
!!$         end if
!!$
!!$      else
!!$
!!$         write(*,'(a)') "Vertex to edge info not computed"
!!$
!!$      end if
!!$
!!$    end block vertex_edge

  end subroutine invert_connectivities

  type(logical) function initialize(this)

    class(mesh), intent(inout) :: this

    call this % invert_connectivities()

    !-----------------------------------------------------------------!
    ! Find Cell Face conn. by combining two maps CF = CV x VF
    !-----------------------------------------------------------------!

    cell_face_face_cell: block

      integer :: icell, iface

      write(*,*) "Forming CellFace and FaceCell connectivities..."

      call sparse_transpose_matmul(&
           & this % num_cell_vertices, this % cell_types, this % cell_tags, this % cell_vertices, &
           & this % num_face_vertices, this % face_types, this % face_tags, this % face_vertices, &
           & this % num_cell_faces, this % cell_faces, this % cell_faces_type, &
           & this % num_face_cells, this % face_cells, this % face_cells_type)

      if (allocated(this % cell_faces)) then

         write(*,'(a,i8,a,i8)') &
              & "Cell to face info for", min(this % max_print,this % num_cells), &
              & " cells out of ", this % num_cells

         do icell = 1, min(this % max_print,this % num_cells)
            write(*,*) &
                 & 'cell ['  , icell, ']',&
                 & 'nfaces [', this % num_cell_faces(icell), ']',&
                 & 'faces [' , this % cell_faces(&
                 & 1:this % num_cell_faces(icell),icell), ']'
         end do

      else

         write(*,'(a)') "Cell to face info not computed"

      end if

      if (allocated(this % face_cells)) then

         do iface = 1, min(this % max_print,this % num_faces)
            print *, &
                 & 'face [', this % face_numbers(iface), ']', &
                 & 'cells [', this % face_cells(1 : this % num_face_cells(iface),iface), ']'
         end do

         if (minval(this % num_face_cells) .lt. 1) then
            write(error_unit, *) 'Error: There are faces not mapped to a cell'
         end if

      else

         write(*,'(a)') "Face to cell info not computed"

      end if

    end block cell_face_face_cell

    !-----------------------------------------------------------------!
    ! Evaluate all geometric quantities needed for FVM assembly
    !-----------------------------------------------------------------!

    geometric_quantities : block

      write(*,*) 'Calculating mesh geometry information'

      call this % evaluate_cell_centers()
      call this % evaluate_face_centers_areas()
      call this % evaluate_face_tangents_normals()
      call this % evaluate_cell_volumes()
      call this % evaluate_centroidal_vector()
      call this % evaluate_face_deltas()
      call this % evaluate_face_weight()
      call this % evaluate_vertex_weight()

    end block geometric_quantities

    ! Signal that all tasks are complete
    initialize = .true.

  end function initialize

  subroutine evaluate_vertex_weight(this)

    class(mesh), intent(inout) :: this
    integer,  allocatable :: cells(:)
    real(dp) :: total, dcell
    integer  :: icell, ivertex

    write(*,*) 'Evaluating face weights for interpolation from cells to vertex'

    allocate(cells(maxval(this % num_vertex_cells)))

    allocate( &
         & this % vertex_cell_weights( &
         & 1:maxval(this % num_vertex_cells), &
         & this % num_faces) &
         & )
    this % vertex_cell_weights = 0

    do concurrent (ivertex = 1 : this % num_vertices)

       ! actual cells numbers
       cells(1:this % num_vertex_cells(ivertex)) = &
            & this % vertex_cells(1:this % num_vertex_cells(ivertex), &
            & ivertex)

       total = real(0,dp)

       do icell = 1, this % num_vertex_cells(ivertex)

          dcell = distance(this % cell_centers(:,cells(icell)), this % vertices(:,ivertex))

          this % vertex_cell_weights(icell,ivertex) = real(1,dp)/dcell

          total = total + this % vertex_cell_weights(icell,ivertex)

       end do

       this % vertex_cell_weights(1:this % num_vertex_cells(ivertex),ivertex) = &
            & this % vertex_cell_weights(1:this % num_vertex_cells(ivertex),ivertex)/total

    end do

    do ivertex = 1, min(this % max_print,this % num_vertices)
       write(*,*) &
            & "vertex [", this % vertex_numbers(ivertex), ']', &
            & "weights [", this % vertex_cell_weights(&
            & 1:this % num_vertex_cells(ivertex),ivertex), ']'
    end do

    deallocate(cells)

  end subroutine evaluate_vertex_weight

  subroutine evaluate_face_weight(this)

    class(mesh), intent(inout) :: this
    integer  :: iface
    integer  :: cellindex1, cellindex2
    real(dp) :: xcellcenter1(3), xcellcenter2(3), xfacecenter(3)
    real(dp) :: d1, d2
    real(dp) :: dinv1, dinv2
    real(dp) :: weight

    write(*, *) 'Evaluating face weights for interpolation from cells to face'
    allocate(this % face_cell_weights(2, this % num_faces))
    do concurrent (iface = 1 : this % num_faces)

       ! first cell is found for all faces
       cellindex1   = this % face_cells(1, iface)
       xcellcenter1 = this % cell_centers(:, cellindex1)
       xfacecenter  = this % face_centers(:,iface)
       d1           = distance(xcellcenter1, xfacecenter)
       dinv1        = 1.0_dp/d1

       ! Extract the second cell if this is not a boundary face/hole
       if (this % num_face_cells(iface) .ne. 1) then
          cellindex2   = this % face_cells(2, iface)
          xcellcenter2 = this % cell_centers(:, cellindex2)
          d2           = distance(xcellcenter2, xfacecenter)
          dinv2        = 1.0_dp/d2
       else
          dinv2        = 0.0_dp
       end if

       weight       = dinv1/(dinv1+dinv2)

       this % face_cell_weights(1:2,iface) = [weight,1.0_dp - weight]

    end do

    do iface = 1, min(this % max_print,this % num_faces)
       write(*,*) &
            & "face [", iface, "] ",&
            & "weight [", this % face_cell_weights(1:2,iface), "] "
    end do

  end subroutine evaluate_face_weight

!!$  subroutine evaluate_face_deltas(this)
!!$
!!$    class(mesh), intent(inout) :: this
!!$    integer  :: gface, gcell, lface
!!$    real(dp) :: fn(3)
!!$
!!$    write(*,*) "Evaluating face deltas"
!!$    allocate(this % face_deltas(this % num_faces))
!!$
!!$    do concurrent (gface=1:this % num_faces)
!!$
!!$       ! First cell belonging to the face
!!$       gcell = this % face_cells(1, gface)
!!$
!!$       ! Face number in local numbering
!!$       ! avoid this??
!!$       lface = find(this % cell_faces(:,gcell), gface)
!!$
!!$       ! Index into normal array
!!$       fn =  this % cell_face_normals(:, lface, gcell)
!!$
!!$       ! Take absolute value of dot product
!!$       this % face_deltas(gface) = abs(dot_product(this % lvec(1:3,gface), fn))
!!$
!!$       print *, "face [", gface, "] ",&
!!$            & "delta [", this % face_deltas(gface), "] ",&
!!$            & "skewness t.l [", &
!!$            & dot_product(this % lvec(1:3,gface)                    , this % cell_face_tangents(:, lface, gcell)), "] ", &
!!$            & "orthogonality t.n [", &
!!$            & dot_product(this % cell_face_tangents(:, lface, gcell), this % cell_face_normals(:, lface, gcell)), "] ", &
!!$            & this % cell_face_normals(:, lface, gcell)
!!$
!!$    end do
!!$
!!$    ! Check for negative volumes
!!$    if (abs(minval(this % face_deltas)) < 1.0d-10) then
!!$       print *, 'collinear faces/bad cell?'
!!$       error stop
!!$    end if
!!$
!!$  end subroutine evaluate_face_deltas

  subroutine evaluate_centroidal_vector(this)

    class(mesh), intent(inout) :: this
    integer :: iface

    write (*,*) "Evaluating centroidal vector..."

    allocate(this % lvec(3,this % num_faces))
    this % lvec = real(0,dp)

    do concurrent (iface = 1 : this % num_faces)

       associate(face_cells => this % face_cells(1:this % num_face_cells(iface),iface))

         ! centroidal vector joining the centroid of two cells
         check_boundary: if (this % num_face_cells(iface) .eq. 1) then
            ! boundary faces
            this % lvec(:,iface) = this % face_centers(:,iface) - this % cell_centers(:,face_cells(1))
         else
            ! internal faces
            this % lvec(:,iface) = this % cell_centers(:,face_cells(2)) - this % cell_centers(:,face_cells(1))
         end if check_boundary

       end associate

    end do

  end subroutine evaluate_centroidal_vector

  subroutine evaluate_cell_volumes(this)

    class(mesh), intent(inout) :: this

    ! Use divergence theorem to find volumes
    integer :: lcell, lface, gface

    write (*,*) "Evaluating cell volumes..."

    allocate(this % cell_volumes (this % num_cells))
    this % cell_volumes = real(0,dp)

    ! V = \sum_f nx_f \times  xmid_f \times A_f
    do concurrent (lcell = 1 : this % num_cells)
       do lface = 1, this % num_cell_faces(lcell)
          ! Global face index
          gface = this % cell_faces(lface, lcell)
          associate( &
               & face_center => this % face_centers(:,gface)/real(3,dp), &
               & face_normal => this % cell_face_normals(:,lface,lcell),&
               & face_area => this % face_areas(gface))
            this % cell_volumes(lcell) = this % cell_volumes(lcell) + &
                 & face_area*dot_product(face_normal, face_center)
          end associate
       end do
    end do

    ! Check for negative volumes
    associate(minvol => minval(this % cell_volumes))
      if (minvol .lt. 0) then
         print *, 'negative volume encountered', minvol
         error stop
      end if
    end associate

  end subroutine evaluate_cell_volumes

  subroutine evaluate_face_deltas(this)

    class(mesh), intent(inout) :: this

    ! Use divergence theorem to find volumes
    integer :: lcell, lface, gface

    write (*,*) "Evaluating face delta..."

    allocate(this % face_deltas (this % num_faces))
    this % face_deltas = real(0,dp)

    do concurrent (lcell = 1 : this % num_cells)

       do lface = 1, this % num_cell_faces(lcell)

          gface = this % cell_faces(lface, lcell)

          this % face_deltas(lface) = &
               & abs(dot_product(this % lvec(:,gface), this % cell_face_normals(:, lface, lcell)))

       end do

    end do

  end subroutine evaluate_face_deltas

  subroutine evaluate_face_centers_areas(this)

    class(mesh), intent(inout) :: this

    ! Currently the length as its a 1D face
    type(integer) :: iface
    real(dp)      :: n(3)

    write(*, *) 'Evaluating face centers and areas'

    allocate(this % face_areas(this % num_faces))
    this % face_areas = real(0,dp)

    allocate(this % face_centers(3,this % num_faces))
    this % face_centers = real(0,dp)

    face_loop: do concurrent(iface = 1 : this % num_faces)

       associate(facenodes => this % face_vertices(1:this % num_face_vertices(iface),iface))

         ! Compute the coordinates of face centers
         associate(num_vertices => real(this % num_face_vertices(iface), kind=dp))
           this % face_centers(:,iface) = sum(this % vertices(:,facenodes),dim=2)/num_vertices
         end associate

         ! triangle
         associate(&
              & t12 => this % vertices(:,facenodes(2)) - this % vertices(:,facenodes(1)), &
              & t13 => this % vertices(:,facenodes(3)) - this % vertices(:,facenodes(1)) &
              & )

           call cross_product(t12,t13,n)

           this % face_areas(iface) = norm2(n)/real(2,dp)

           ! if quadrilateral add the second triangle
           if (this % num_face_vertices(iface) .gt. 3) then

              associate (t14 => this % vertices(:,facenodes(4)) - this % vertices(:,facenodes(1)))

                call cross_product(t13,t14,n)

                ! add the area of second triangle
                this % face_areas(iface) = this % face_areas(iface) + norm2(n)/real(2,dp)

              end associate

           end if

         end associate

       end associate

    end do face_loop

    ! Check for negative areas
    associate(minarea => minval(this % face_areas))
      if (minarea .lt. 0) then
         print *, 'negative area encountered', minarea
         error stop
      end if
    end associate

  end subroutine evaluate_face_centers_areas

  subroutine evaluate_face_centers_areas_2d(this)

    class(mesh), intent(inout) :: this

    ! Currently the length as its a 1D face
    type(integer) :: iface

    write(*, *) 'Evaluating face centers and areas'

    allocate(this % face_areas(this % num_faces))
    allocate(this % face_centers(3,this % num_faces))

    do concurrent(iface = 1 : this % num_faces)

       ! Area calculation is complicated in 3D
       associate(facenodes => this % face_vertices(:,iface))

         ! Compute the coordinates of face centers
         this % face_centers(1:3, iface) = &
              & sum(this % vertices(1:3, facenodes),dim=2)/&
              & real(2,kind=dp) ! this face has 2 edges

         ! Compute face areas
         associate(v1 => this % vertices(:,facenodes(1)), &
              & v2 => this % vertices(:,facenodes(2))  )
           this % face_areas(iface) = distance(v1, v2)
         end associate

       end associate

    end do

    ! Check for zero areas
    if (abs(minval(this % face_areas)) .lt. 10.0d0*tiny(1.0d0)) then
       print *, 'same points/bad face?'
       error stop
    end if

  end subroutine evaluate_face_centers_areas_2d

  subroutine evaluate_cell_centers(this)

    class(mesh), intent(inout) :: this

    ! Find cell centers O = (A + B + C + ...) /count(vertices)
    type(integer) :: icell

    write(*,*) 'Evaluating cell centers'

    allocate(this % cell_centers(3, this % num_cells))

    do concurrent(icell=1:this % num_cells)
       associate(&
            & num_vertices => real(this % num_cell_vertices(icell), kind=dp), &
            & vids => this % cell_vertices(1:this % num_cell_vertices(icell), icell) &
            & )
         this % cell_centers(:, icell) = sum(this % vertices(:,vids),dim=2)/num_vertices
       end associate
    end do

  end subroutine evaluate_cell_centers

  subroutine evaluate_face_tangents_normals(this)

    class(mesh), intent(inout) :: this
    integer  :: icell, iface, gface
    real(dp) :: tmp(3)

    write(*,*) 'Evaluating face tangents normals'

    allocate(this % cell_face_normals(3, maxval(this % num_cell_faces), this % num_cells))
    this % cell_face_normals = real(0,dp)

    ! 3d, 2 tangents per face
    allocate(this % cell_face_tangents(3, 2, maxval(this % num_cell_faces), this % num_cells))
    this % cell_face_tangents = real(0,dp)

    ! loop cells
    cell_loop: do concurrent (icell = 1 : this % num_cells)

       ! loop faces of each cell
       face_loop: do iface = 1, this % num_cell_faces(icell)

          ! gface = this % face_numbers
          gface = this % cell_faces(iface, icell)

          associate(ifv => this % face_vertices(1:this % num_face_vertices(gface),gface), &
               & normal => this % cell_face_normals(:,iface,icell) )

            associate(&
                 & t12 => this % vertices(:,ifv(2)) - this % vertices(:,ifv(1)), &
                 & t13 => this % vertices(:,ifv(3)) - this % vertices(:,ifv(1))  &
                 & )

              call cross_product(t12, t13, normal)

              ! if quadrilateral add the second triangle (may not need this at all)
              if (this % num_face_vertices(gface) .gt. 3) then

                 associate (t14 => this % vertices(:,ifv(4)) - this % vertices(:,ifv(1)))

                   call cross_product(t13,t14,tmp)

                   normal = normal + tmp

                 end associate

              end if

              normal = normal/norm2(normal)

              ! determine whether the normal is inward or outward by
              ! projecting it to the vector connecting cell center and
              ! face center
              if (dot_product(normal, this % face_centers(:,gface) - this % cell_centers(:,icell)) &
                   & .lt. real(0,dp)) then

                 ! t1, t2, n forms a local orthonormal coordinate
                 ! system (notice we flip the order of cross pdt to
                 ! account for sign change). That is, the normal will
                 ! be negative only when we do t2 cross t1.

                 normal = - normal
                 this % cell_face_tangents(:, 1, iface, icell) = t13/norm2(t13)
                 this % cell_face_tangents(:, 2, iface, icell) = t12/norm2(t12)

              else

                 ! t1, t2, n forms a local orthonormal coordinate system
                 this % cell_face_tangents(:, 1, iface, icell) = t12/norm2(t12)
                 this % cell_face_tangents(:, 2, iface, icell) = t13/norm2(t13)

              end if

            end associate

          end associate

       end do face_loop

  end do cell_loop

  end subroutine evaluate_face_tangents_normals

  subroutine evaluate_face_tangents_normals_2d(this)

    class(mesh), intent(inout) :: this

    integer  :: icell, iface, gface
    real(dp) :: t(3), n(3), tcn(3) ! all spatial dim
    integer  :: ifv(2)

    write(*,*) 'Evaluating face tangents normals'

    allocate(this % cell_face_normals (3, maxval(this % num_cell_faces), this % num_cells))
    allocate(this % cell_face_tangents(3, 1, maxval(this % num_cell_faces), this % num_cells))

    ! loop cells
    do concurrent (icell = 1 : this % num_cells)

       ! get cell verties
       associate( icv =>  this % cell_vertices(:, icell) )

         ! loop faces of each cell
         do iface = 1, this % num_cell_faces(icell)

            if (iface .eq. this % num_cell_faces(icell)) then
               ifv(1) = icv(iface)
               ifv(2) = icv(1)
            else
               ifv(1) = icv(iface)
               ifv(2) = icv(iface+1)
            end if

            ! find the face vertex in cell order
            gface = this % cell_faces(iface,icell)

            t = this % vertices(:,ifv(2)) - this % vertices(:,ifv(1))
            t = t/norm2(t)

            ! By anticlockwise convention
            n(1) =  t(2)
            n(2) = -t(1)
            n(3) = 0

            ! Sanity check if the normal if facing out of the face
            call cross_product(n,t,tcn)
            if (abs(tcn(3) - 1.0d0) > 1.0d-10) then ! tangent cross normal should give +k vector
               print *, 'face', gface, 'of cell', icell, 'has inward/non-unit normal', tcn(3)
               error stop
            end if

            this % cell_face_normals (:, iface, icell) = n
            this % cell_face_tangents(:, 1, iface, icell) = t

         end do

       end associate

    end do

  end subroutine evaluate_face_tangents_normals_2d

  !===================================================================!
  ! Counts and returns the number of elements with supplied number of
  ! nodes
  !===================================================================!

  pure elemental type(integer) function get_num_elems(this, num_elem_nodes)

    class(mesh)  , intent(in) :: this
    type(integer), intent(in) :: num_elem_nodes

    get_num_elems = count(MASK=(this % num_cell_vertices .eq. num_elem_nodes))

  end function get_num_elems

  !===================================================================!
  ! Constructor for mesh creation
  !===================================================================!

  subroutine to_string(this)

    class(mesh), intent(in) :: this

    integer :: icell, ivertex, iface, iedge, itag

    write(*,*) 'Number of vertices :', this % num_vertices
    write(*,*) 'Number of cells    :', this % num_cells
    write(*,*) 'Number of faces    :', this % num_faces

    do itag = 1, this % num_tags
       write(*,'(1x,a,i0,a,a,a,a)') &
            & "tag number [", this % tag_numbers(itag) , "] ", &
            & "info [", this % tag_info(itag) % str, "] "
    end do

    if (this % num_vertices .gt. 0) then
       write(*,'(a,i8,a,i8)') "Vertex info for ", min(this % max_print,this % num_vertices), &
            & ' vertices out of ', this % num_vertices
       write(*,*) "number tag x y z"
       do ivertex = 1, min(this % max_print,this % num_vertices)
          write(*,'(i8,i2,3ES15.3)') &
               & this % vertex_numbers(ivertex), &
               & this % vertex_tags(ivertex), &
               & this % vertices(:, ivertex)
       end do
    end if

    if (this % num_cells .gt. 0) then
       write(*,'(a,i8,a,i8)') "Cell info for ", min(this % max_print,this % num_cells), &
            & ' cells out of ', this % num_cells
       write(*,*) "cno ctag ncv iverts"
       do icell = 1, min(this % max_print,this % num_cells)
          write(*,'(i8,i2,i2,10i8)') &
               & this % cell_numbers(icell), &
               & this % cell_tags(icell), &
               & this % num_cell_vertices(icell), &
               & this % cell_vertices(1:this % num_cell_vertices(icell), icell)
       end do
    end if

    if (this % num_faces .gt. 0) then
       write(*,'(a,i8,a,i8)') "Face info for ", min(this % max_print,this % num_faces), &
            & ' faces out of ', this % num_faces
       write(*,*) "fno ftag nfv iverts"
       do iface = 1, min(this % max_print,this % num_faces)
          write(*,'(i8,i2,i2,10i8)') &
               & this % face_numbers(iface), &
               & this % face_tags(iface), &
               & this % num_face_vertices(iface), &
               & this % face_vertices(1:this % num_face_vertices(iface), iface)
       end do
    end if

    if (this % num_edges .gt. 0) then
       write(*,'(a,i8,a,i8)') "Edge info for ", min(this % max_print,this % num_edges), &
            & ' edges out of ', this % num_edges
       write(*,*) "eno etag nev iverts"
       do iedge = 1, min(this % max_print,this % num_edges)
          write(*,'(i8,i2,i2,10i8)') &
               & this % edge_numbers(iedge), &
               & this % edge_tags(iedge), &
               & this % num_edge_vertices(iedge), &
               & this % edge_vertices(1:this % num_edge_vertices(iedge), iedge)
       end do
    end if

    if (this % initialized .eqv. .true.) then

       write(*,*) "Cell Geo. Data [index] [center] [volume] "
       do icell = 1, min(this % max_print,this % num_cells)
          write(*,*) &
               & "local number [", this % cell_numbers(icell)   ,"] ", &
               & "center [", this % cell_centers(:,icell) ,"] ", &
               & "volume [", this % cell_volumes(icell)   ,"] "
       end do
       write(*,*) "total volume ", sum(this % cell_volumes)

       write(*,*) "Face Data [index] [center] [area] "
       do iface = 1, min(this % max_print,this % num_faces)
          write(*,*) &
               & "num [",iface,"] ", &
               & "center [",this % face_centers(:, iface),"] ", &
               & "area [",this % face_areas(iface),"] ", &
               & "lvec [",this % lvec(:,iface),"] ", &
               & "delta [",this % face_deltas(iface),"] "
       end do

    end if

  end subroutine to_string

end module class_mesh
