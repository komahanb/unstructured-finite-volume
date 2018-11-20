!=====================================================================!
! Unstructured mesh handler.
! 
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_mesh

  use iso_fortran_env       , only : dp => REAL64, error_unit
  use interface_mesh_loader , only : mesh_loader    
  use module_mesh_utils
  
  implicit none

  private
  public :: mesh
  
  ! Constructor
  interface mesh
     module procedure create_mesh
     module procedure create_mesh_from_file
  end interface mesh

  !-------------------------------------------------------------------! 
  ! Mesh datatype. A collection of vertices, cells and faces.
  !-------------------------------------------------------------------!
  
  type :: mesh ! rename as topology?
     
     !================================================================!
     ! Basic Topology information
     !================================================================!

     ! Fundamental vertex info
     integer :: num_vertices
     real(dp) , allocatable :: vertices(:,:)             ! [[x,y,z],1:nvertices]
     integer  , allocatable :: vertex_numbers(:)
     integer  , allocatable :: vertex_tags(:)

     ! Fundamental face info
     integer :: num_edges
     integer  , allocatable :: edge_numbers(:)
     integer  , allocatable :: edge_tags(:)
     integer  , allocatable :: edge_vertices(:,:)        ! [[v1,v2],1:nedges]
     integer  , allocatable :: num_edge_vertices(:)

     ! Fundamental face info
     integer :: num_faces
     integer  , allocatable :: face_numbers(:)          
     integer  , allocatable :: face_tags(:)          
     integer  , allocatable :: face_vertices(:,:)        ! [[v1,v2],1:nfaces]
     integer  , allocatable :: num_face_vertices(:)

     ! Fundamental cell info
     integer :: num_cells
     integer  , allocatable :: cell_numbers(:)
     integer  , allocatable :: cell_tags(:)
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
     integer  , allocatable :: num_face_cells(:)         ! [1:nfaces]
     integer  , allocatable :: face_cells(:,:)           ! [[c1,c2...],1:nfaces]
     
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
     real(dp) , allocatable :: cell_gamma(:)             ! [1:ncells]

     ! Derived vertex info
     real(dp) , allocatable :: face_centers(:,:)         ! [[x,y,z],1:nfaces]
     real(dp) , allocatable :: face_areas(:)             ! [1:nfaces]
     real(dp) , allocatable :: face_gamma(:)             ! [1:nfaces]
     real(dp) , allocatable :: face_delta(:)             ! [1:nfaces]
     real(dp) , allocatable :: lvec(:,:)                 ! [[lx,ly,lz],1:nfaces]     
  
     real(dp) , allocatable :: cell_face_tangents(:,:,:) ! [[tx,ty,tz], [f1,f2,f3..] 1:ncells]
     real(dp) , allocatable :: cell_face_normals(:,:,:)  ! [[nx,ny,nz], [f1,f2,f3..] 1:ncells]          
    
     real(dp) , allocatable :: vertex_cell_weights(:,:)  ! [[wc1,wc2...],1:vertices]
     real(dp) , allocatable :: face_cell_weights(:,:)    ! [[wc1,wc2],1:nfaces]

     
     ! Generalize to tag numbers. Right now we can't differentiate
     ! lower face from upper face.
     integer  :: num_boundary_faces
     integer  , allocatable :: boundary_face_face(:)      ! [1:num_boundary_faces]
     integer  , allocatable :: is_face_boundary_face(:)   ! [1:num_faces]
     integer  , allocatable :: is_node_boundary_node(:)   ! [1:num_vertices]
     
   contains

     ! Type bound procedures
     procedure :: get_num_lines
     procedure :: to_string
     procedure :: initialize

     ! Evaluation routines
     procedure :: evaluate_tangents_normals
     procedure :: evaluate_cell_centers
     procedure :: evaluate_face_centers_areas
     procedure :: evaluate_cell_volumes
     procedure :: evaluate_centroidal_vector
     procedure :: evaluate_face_delta
     procedure :: evaluate_face_weight
     procedure :: evaluate_vertex_weight

     ! Destructor
     ! final   :: destruct
     
  end type mesh

contains
  
  subroutine initialize(this)

    class(mesh), intent(inout) :: this
    integer :: i
    integer :: lvertex, gvertex, icell, iface

    ! write(*,'(a)') 'Finding inverse mapping of vertex, edge, face, and cell numbers...'
    ! revere mapping from cell vertices to vertex_cells

    write(*,'(a)') 'Finding inverse topologies...'
    write(*,'(a)') '1. vertex to cell connectivites'

    ! revere mapping from cell vertices to vertex_cells
    do i = 1, this % num_cells
       print *, 'cell', i, this % cell_numbers(i), 'num_cell_vertices', this % num_cell_vertices(i) ,'vertices', &
            & this % cell_vertices(1:this % num_cell_vertices(i),i)
    end do

    call reverse_map( &
         & this % cell_vertices, &
         & this % num_cell_vertices, &
         & this % vertex_cells, &
         & this % num_vertex_cells)
    do i = 1, size(this % vertex_cells, dim=2)
       print *, 'vertex', i, 'num_vertex_cells', this % num_vertex_cells(i) ,'cells', &
            & this % vertex_cells(1:this % num_vertex_cells(i),i)
    end do

    write(*,'(a)') '2. vertex to face connectivites'
    do i = 1, this % num_faces
       print *, 'face', i, this % face_numbers(i), 'num_face_vertices', this % num_face_vertices(i) ,'vertices', &
            & this % face_vertices(1:this % num_face_vertices(i),i)
    end do
    call reverse_map( &
         & this % face_vertices, &
         & this % num_face_vertices, &
         & this % vertex_faces, &
         & this % num_vertex_faces)
    do i = 1, size(this % vertex_faces, dim=2)
       print *, 'vertex', i, 'num_vertex_faces', this % num_vertex_faces(i) ,'cells', &
            & this % vertex_faces(1:this % num_vertex_faces(i),i)
    end do
!!$
!!$    write(*,'(a)') '3. vertex to edge connectivites'
!!$    call reverse_map( &
!!$         & this % edge_vertices, &
!!$         & this % num_edge_vertices, &
!!$         & this % vertex_edges, &
!!$         & this % num_vertex_edges)
!!$    do i = 1, size(this % vertex_edges, dim=2)
!!$       print *, 'vertex', i, 'num_vertex_edges', this % num_vertex_edges(i) ,'cells', &
!!$            & this % vertex_edges(1:this % num_vertex_edges(i),i)
!!$    end do

    !=================================================================!
    ! Intermediate Topologies
    !=================================================================!

    write(*,'(a)') 'Finding intermediate topologies...'
    write(*,'(a)') '1.a. face to cell connectivites'
    call get_face_cells( &
       & this % cell_vertices, this % num_cell_vertices, &
       & this % face_vertices, this % num_face_vertices, &
       & this % face_cells   , this % num_face_cells)
    do iface = 1, size(this % face_cells, dim=2)
       print *, 'face', iface,  'num face cells', this%num_face_cells(iface), 'cells',&
            & this % face_cells(1:this%num_face_cells(iface),iface)
    end do

!!$
!!$    write(*,'(a)') '1.b. face to cell connectivites'
!!$    ! Invert cell_faces
!!$    call reverse_map(cell_faces, num_cell_faces, face_cells, num_face_cells)






!!$
!!$    ! Finding the vertex tags based on faces
!!$    do lvertex = 1, this % num_vertices
!!$       gvertex = this % vertex_numbers(lvertex)
!!$       gface   = this % vertex_faces(1,gvertex)
!!$       this % vertex_tags(lvertex) = this % face_tags(lface(gf))
!!$    end do

    stop

    write(*,'(a)') 'Calculating geometry information...'











    print *, 'finding mesh geometry'
    
    allocate(this % cell_gamma(this % num_cells))      
    this % cell_gamma = 1.0d0

    call this % evaluate_cell_centers()
!!$    print *, 'cell center'
!!$    do i = 1, this % num_cells
!!$       print *, i, this % cell_centers(:,i)
!!$    end do

    call this % evaluate_face_centers_areas()      
!!$    print *, 'face'
!!$    do i = 1, this % num_faces
!!$       print *, i, this % face_areas(i), &
!!$            & this % face_centers(:,i)
!!$    end do

    ! Use divergence theorem to find are
    call this % evaluate_tangents_normals()
    call this % evaluate_cell_volumes()
    call this % evaluate_centroidal_vector()
    call this % evaluate_face_delta()
    call this % evaluate_face_weight()   
    call this % evaluate_vertex_weight()    

  end subroutine initialize
  
  subroutine evaluate_vertex_weight(this)

    class(mesh), intent(inout) :: this
    integer,  allocatable :: cells(:)
    real(dp) :: total, dcell
    integer  :: icell, ivertex

    print *, 'face weights for interpolation from cells to vertex'

    allocate(cells(maxval(this % num_vertex_cells)))

    allocate( &
         & this % vertex_cell_weights( &
         & 1:maxval(this % num_vertex_cells), &
         & this % num_faces) &
         & )
    this % vertex_cell_weights = 0

    do ivertex = 1, this % num_vertices

       cells(1:this % num_vertex_cells(ivertex)) = &
            & this % vertex_cells(1:this % num_vertex_cells(ivertex), ivertex)

       total  = 0.0d0

       do icell = 1, this % num_vertex_cells(ivertex)

          dcell = distance(this % cell_centers(:,icell), this % vertices(:,ivertex))

          this % vertex_cell_weights(icell,ivertex) = 1.0_dp/dcell

          total = total + this % vertex_cell_weights(icell,ivertex)

       end do

       this % vertex_cell_weights(:,ivertex) = this % vertex_cell_weights(:,ivertex)/total

       print *, "vertex", ivertex, this % vertex_cell_weights(1:this % num_vertex_cells(ivertex),ivertex)

    end do

  end subroutine evaluate_vertex_weight

  subroutine evaluate_face_weight(this)

    class(mesh), intent(inout) :: this
    integer  :: iface
    integer  :: cellindex1, cellindex2
    real(dp) :: xcellcenter1(3), xcellcenter2(3), xfacecenter(3)
    real(dp) :: d1, d2
    real(dp) :: dinv1, dinv2
    real(dp) :: weight

    print *, 'face weights for interpolation from cells to face'
    allocate(this % face_cell_weights(2, this % num_faces))      

    !do concurrent(iface = 1: this % num_faces)

    do iface = 1, this % num_faces

       cellindex1   = this % face_cells(1, iface)
       xcellcenter1 = this % cell_centers(:, cellindex1)
       xfacecenter  = this % face_centers(:,iface)       
       d1           = distance(xcellcenter1, xfacecenter)
       dinv1        = 1.0_dp/d1

       ! Extract the second cell if this is not a boundary face
       if (this % is_face_boundary_face(iface) .ne. 1) then
          cellindex2   = this % face_cells(2, iface)
          xcellcenter2 = this % cell_centers(:, cellindex2)
          d2           = distance(xcellcenter2, xfacecenter)
          dinv2        = 1.0_dp/d2
       else
          dinv2        = 0.0_dp                    
       end if

       weight       = dinv1/(dinv1+dinv2)

       this % face_cell_weights(1:2,iface) = [weight,1.0_dp - weight]

       print *, "face weight", iface, this % face_cell_weights(1:2,iface)

    end do

  end subroutine evaluate_face_weight

  subroutine evaluate_face_delta(this)

    class(mesh), intent(inout) :: this
    integer  :: gface, gcell, lface
    real(dp) :: fn(3)

    allocate(this % face_delta(this % num_faces))

    do gface = 1, this % num_faces

       ! First cell belonging to the face
       gcell = this % face_cells(1, gface)

       ! Face number in local numbering
       lface = find(this % cell_faces(:,gcell), gface)

       ! Index into normal array
       fn =  this % cell_face_normals(:, lface, gcell)

       ! Take absolute value of dot product
       this % face_delta(gface) = abs(dot_product(this % lvec(1:3,gface), fn))

       print *, "face", gface, "delta", this % face_delta(gface), &
            & "skewness", dot_product(this % lvec(1:3,gface), &
            & this % cell_face_tangents(:, lface, gcell)), &
            & dot_product(this % cell_face_tangents(:, lface, gcell), &
            & this % cell_face_normals(:, lface, gcell))

    end do

  end subroutine evaluate_face_delta

  subroutine evaluate_centroidal_vector(this)

    class(mesh), intent(inout) :: this
    integer :: iface, cells(2)

    allocate(this % lvec(3,this % num_faces))

    do iface = 1, this % num_faces

       cells = this % face_cells(:,iface)

       if (this % is_face_boundary_face(iface) .eq. 1) then
          ! Boundary faces .or. iface is in bfaces
          this % lvec(:,iface) = this % face_centers(:,iface) - this % cell_centers(:,cells(1))         
       else
          ! Interior face; subtract neighbouring cell centers (not sure which orientation)
          this % lvec(:,iface) = this % cell_centers(:,cells(2)) - this % cell_centers(:,cells(1))          
       end if

    end do

  end subroutine evaluate_centroidal_vector

  subroutine evaluate_cell_volumes(this)

    class(mesh), intent(inout) :: this

    ! Use divergence theorem to find volumes
      integer :: lcell, lface, gface

      allocate(this % cell_volumes (this % num_cells))
      this % cell_volumes = 0_dp      

      ! V = \sum_f nx_f \times  xmid_f \times A_f
      do lcell = 1, this % num_cells
         this % cell_volumes(lcell) = 0.0d0
         do lface = 1, this % num_cell_faces(lcell)
            ! Global face index
            gface = this % cell_faces(lface, lcell)
            associate( &
                 & xmid => this % face_centers(1,gface), &
                 & nx   => this % cell_face_normals(1,lface,lcell),&
                 & area => this % face_areas(gface))
            this % cell_volumes(lcell) = this % cell_volumes(lcell) + &
                 & nx*xmid*area              
          end associate
       end do
       print *, 'cell', lcell, 'volume', this % cell_volumes(lcell)       
    end do

end subroutine evaluate_cell_volumes
  
  subroutine evaluate_face_centers_areas(this)

    class(mesh), intent(inout) :: this


    ! Currently the length as its a 1D face
    type(integer) :: iface

    allocate(this % face_areas(this % num_faces))
    allocate(this % face_centers(3,this % num_faces))         

    do concurrent(iface = 1 : this % num_faces)

       ! Area calculation is complicated in 3D
       associate(facenodes => this % face_vertices(:,iface))

         ! Compute the coordinates of face centers
       this % face_centers(1:3, iface) = &
            & sum(this % vertices(1:3, facenodes),dim=2)/&
            & real(2,kind=dp) ! this face has 2 edges

       associate(v1 => this % vertices(:,facenodes(1)), &
            & v2 => this % vertices(:,facenodes(2))  )

         ! Compute face areas
       this % face_areas(iface) = distance(v1, v2)

        end associate
       
      end associate

   end do
   
 end subroutine evaluate_face_centers_areas
 
  subroutine evaluate_cell_centers(this)

    class(mesh), intent(inout) :: this

     ! Find cell centers O = (A + B + C) /3
    type(integer) :: icell

    !print *, 'num_vertices for each cell', this % num_cell_vertices

    allocate(this % cell_centers(3, this % num_cells))
    
    do concurrent(icell = 1 : this % num_cells)
       this % cell_centers(:, icell) = sum(&
            & this % vertices(&
            & :, this % cell_vertices(&
            & 1:this % num_cell_vertices(icell),icell)&
            & ), dim=2)&
            & /real(this % num_cell_vertices(icell), kind=dp)
    end do
    
  end subroutine evaluate_cell_centers
  
  subroutine evaluate_tangents_normals(this)

    class(mesh), intent(inout) :: this

    integer  :: icell, iface, cell1, gface
    real(dp) :: t(3), a(3), b(3), n(3), tcn(3) ! all spatial dim
    integer  :: ifv(2)

    allocate(this % cell_face_normals (3, maxval(this % num_cell_faces), this % num_cells))
    allocate(this % cell_face_tangents(3, maxval(this % num_cell_faces), this % num_cells))

    ! loop cells
    do icell = 1, this % num_cells

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
          if (tcn(3) < 0) then
             print *, 'face', gface, 'of cell', icell, 'has inward normal'
             error stop
          end if

          this % cell_face_normals (:, iface, icell) = n
          this % cell_face_tangents(:, iface, icell) = t
          
          print *, icell, gface, n(1:2), tcn

       end do

       end associate

    end do

  end subroutine evaluate_tangents_normals

  !===================================================================!
  ! Constructor for mesh creation
  !===================================================================!
  
  subroutine to_string(this)

    class(mesh), intent(in) :: this

    integer :: icell, ivertex, iface

    write(*,*) 'Number of vertices :', this % num_vertices
    write(*,*) 'Number of cells    :', this % num_cells
    write(*,*) 'Number of faces    :', this % num_faces

    write(*,*) "Vertex Info:"
    write(*,*) "number tag x y z"
    do ivertex = 1, this % num_vertices
       write(*,'(i6,i2,3E15.6)') &
            & this % vertex_numbers(ivertex), &
            & this % vertex_tags(ivertex), &
            & this % vertices(:, ivertex)
    end do
    
    write(*,*) "Cell Info:"
    write(*,*) "cno ctag ncv iverts"
    do icell = 1, this % num_cells
       write(*,'(i6,i2,i2,10i6)') &
            & this % cell_numbers(icell), &
            & this % cell_tags(icell), &
            & this % num_cell_vertices(icell), &
            & this % cell_vertices(1:this % num_cell_vertices(icell), icell)
    end do
    
    write(*,*) "Face Info:"
    write(*,*) "fno ftag nfv iverts"
    do iface = 1, this % num_faces
       write(*,'(i6,i2,i2,10i6)') &
            & this % face_numbers(iface), &
            & this % face_tags(iface), &
            & this % num_face_vertices(iface), &
            & this % face_vertices(1:this % num_face_vertices(iface), iface)
    end do

!!$    write(*,*) "Face Data [index] [center] [volume]"
!!$    do icell = 1, this % num_cells
!!$       write(*,*) "[",icell,"]", "[",this % cell_centers(:, icell),"]", &
!!$            & "[",this % cell_volumes(icell),"]"
!!$    end do
!!$
!!$    write(*,*) "Face Data [index] [center] [area]"
!!$    do iface = 1, this % num_faces
!!$       write(*,*) "[",iface,"]", "[",this % face_centers(:, iface),"]", &
!!$            & "[",this % face_areas(iface),"]"
!!$    end do
!!$    
!!$    write(*,*) "Face to Face Connectivity"
!!$    do icell = 1, this % num_cells
!!$       write(*,*) icell, this % cell_faces(1:this % num_cell_faces(icell), icell)
!!$    end do
!!$
!!$    write(*,*) "Face to Cell Connectivity"
!!$    do iface = 1, this % num_faces
!!$       write(*,*) iface, this % face_cells(1:this % num_face_cells(iface), iface)
!!$    end do

  end subroutine to_string
    
  type(mesh) function create_mesh_from_file(loader) result(me)

    ! Arguments
    class(mesh_loader), intent(in) :: loader
    
    ! Get the fundamental information needed 
    call loader % get_mesh_data( &
         & me % num_vertices, me % vertex_numbers, me % vertex_tags , me % vertices ,  & 
         & me % num_edges   , me % edge_numbers  , me % edge_tags   , me % edge_vertices , me % num_edge_vertices , &
         & me % num_faces   , me % face_numbers  , me % face_tags   , me % face_vertices , me % num_face_vertices , &
         & me % num_cells   , me % cell_numbers  , me % cell_tags   , me % cell_vertices , me % num_cell_vertices   &
         & )

    ! Check allocations and print error messages and stop

    ! Sanity check (make sure numbering is continuous), although it may not start from one
    if (me % num_vertices .gt. 0 .and. maxval(me % vertex_numbers) -  minval(me % vertex_numbers) + 1 .ne. me % num_vertices) &
         & error stop
    if (me % num_edges    .gt. 0 .and. maxval(me % edge_numbers  ) -  minval(me % edge_numbers  ) + 1 .ne. me % num_edges   ) &
         & error stop
    if (me % num_faces    .gt. 0 .and. maxval(me % face_numbers  ) -  minval(me % face_numbers  ) + 1 .ne. me % num_faces   ) &
         & error stop
    if (me % num_cells    .gt. 0 .and. maxval(me % cell_numbers  ) -  minval(me % cell_numbers  ) + 1 .ne. me % num_cells   ) &
         & error stop

    call me % to_string()

    print *, 'passed initialization check'

    stop

    ! Perform initialization tasks
    call me % initialize()

  end function create_mesh_from_file

  !===================================================================!
  ! Constructor for mesh creation
  !===================================================================!

  type(mesh) function create_mesh(xptsin, facesin, connin, &
       & num_cell_verticesin, &
       & cell_facesin, num_cell_facesin, &
       & face_cellsin, num_face_cellsin, &
       & boundary_face_facein) result(this)
    
    real(dp) , intent(in) :: xptsin(:,:)
    integer  , intent(in) :: facesin(:,:)
    integer  , intent(in) :: connin(:,:)
    integer  , intent(in) :: num_cell_verticesin(:)
    integer  , intent(in) :: cell_facesin(:,:)
    integer  , intent(in) :: face_cellsin(:,:)
    integer  , intent(in) :: num_cell_facesin(:)
    integer  , intent(in) :: num_face_cellsin(:)
    integer  , intent(in) :: boundary_face_facein(:)
    integer :: iface, i, ivertex
    
    allocate(this % vertices          , source = xptsin)
    allocate(this % face_vertices     , source = facesin)
    allocate(this % cell_vertices     , source = connin)
    allocate(this % num_cell_vertices , source = num_cell_verticesin)
    allocate(this % cell_faces        , source = cell_facesin)
    allocate(this % face_cells        , source = face_cellsin)
    allocate(this % num_cell_faces    , source = num_cell_facesin)
    allocate(this % num_face_cells    , source = num_face_cellsin)
    allocate(this % boundary_face_face    , source = boundary_face_facein)

    this % num_vertices       = size(this % vertices       , dim=2)
    this % num_cells          = size(this % cell_vertices  , dim=2)
    this % num_faces          = size(this % face_vertices  , dim=2)
    this % num_boundary_faces = size(this % boundary_face_face , dim=1)

    ! Find if a face is boundary face (Tag faces with index) face_tags
    ! [t1,t2,1:nfaces]
    allocate(this % is_face_boundary_face(this % num_faces))
    do iface = 1, this % num_faces
       if (any(this % boundary_face_face .eq. iface) .eqv. .true.) then
          this % is_face_boundary_face(iface) = 1
       else
          this % is_face_boundary_face(iface) = 0
       end if
    end do

    ! Find if a node is boundary node (node tag)
    allocate(this % is_node_boundary_node(this % num_vertices))
    do ivertex = 1, this % num_vertices
       if (any(this % face_vertices(:,this % boundary_face_face) &
            & .eq. ivertex) .eqv. .true.) then
          this % is_node_boundary_node(ivertex) = 1
       else
          this % is_node_boundary_node(ivertex) = 0
       end if
    end do

    call this % initialize()
    
!!$    ! revere mapping from cell vertices to vertex_cells
!!$    call reverse_map(this % cell_vertices, this % num_cell_vertices, &
!!$         & this % vertex_cells, this % num_vertex_cells)
!!$    do i = 1, size(this % vertex_cells, dim=2)
!!$       print *, 'vertex', i, 'cells', this % vertex_cells(1:this%num_vertex_cells(i),i)
!!$    end do
!!$    
!!$    allocate(this % cell_gamma(this % num_cells))      
!!$    this % cell_gamma = 1.0d0
!!$
!!$    call this % evaluate_cell_centers()
!!$    print *, 'cell center'
!!$    do i = 1, this % num_cells
!!$       print *, i, this % cell_centers(:,i)
!!$    end do
!!$
!!$    call this % evaluate_face_centers_areas()      
!!$    print *, 'face'
!!$    do i = 1, this % num_faces
!!$       print *, i, this % face_areas(i), &
!!$            & this % face_centers(:,i)
!!$    end do
!!$
!!$    ! Use divergence theorem to find are
!!$    call this % evaluate_tangents_normals()
!!$    call this % evaluate_cell_volumes()
!!$    call this % evaluate_centroidal_vector()
!!$    call this % evaluate_face_delta()
!!$    call this % evaluate_face_weight()   
!!$    call this % evaluate_vertex_weight()    

  end function create_mesh
  
  !===================================================================!
  ! Utility function for get number of lines in mesh file
  !===================================================================!

  integer function get_num_lines(this, filename) result(nlines)

    class(mesh)     , intent(in)  :: this 
    character(len=*), intent(in)  :: filename
    integer :: stat

    nlines = 0 
    open (111, file = filename)
    do
       read(111,*,iostat=stat)
       if (stat .ne. 0) exit
       nlines = nlines + 1
    end do
    close (111)

  end function get_num_lines

end module class_mesh

!===================================================================!
! Main program to run basic test of functionalities of mesh module.
!===================================================================!

module mesh_loader

  ! Dependencies
  use iso_fortran_env , only : dp => REAL64
  use module_mesh_utils
  use class_mesh     , only : mesh

  implicit none
  
contains
  
  subroutine create(grid)
    
    type(mesh), intent(out) :: grid

    real(dp) , allocatable :: xpts(:,:)        
    integer  , allocatable :: cell_vertices(:,:)
    integer  , allocatable :: num_cell_vertices(:)

    integer  , allocatable :: vertex_cells(:,:)
    integer  , allocatable :: num_vertex_cells(:)
    integer  , allocatable :: vertex_cell_ptr(:)

    integer  , allocatable :: face_vertices(:,:)
    integer  , allocatable :: num_face_vertices(:)  
    integer  , allocatable :: vertex_faces(:,:)    
    integer  , allocatable :: num_vertex_faces(:)

    integer  , allocatable :: cell_faces(:,:)
    integer  , allocatable :: face_cells(:,:)
    integer  , allocatable :: num_face_cells(:)
    integer  , allocatable :: num_cell_faces(:)
    integer  , allocatable :: boundary_faces(:)

    integer :: i
    integer :: npoints, ncells, nfaces, nbfaces
    integer :: icell, iface
    integer :: v1, v2
    integer :: face_ptr, ctr

!!$  load_mesh : block
    ! Load mesh and extract xpts and cell_vertices  
!!$    npoints = msh % get_num_lines("coordinates_10.input")
!!$    ncells  = msh % get_num_lines("elements_10.input")
!!$
!!$    ! Load veritices
!!$    allocate(xpts(3,npoints))
!!$    xpts = 0.0_dp
!!$    open(unit=10, file="coordinates_10.input")
!!$    do i = 1, npoints
!!$       read(10,*) xpts(1,i), xpts(2,i)
!!$    end do
!!$    close(10)
!!$    
!!$    ! Load cells assuming triangles
!!$    allocate(cell_vertices(3,ncells))
!!$    open(unit=10, file="elements_10.input")
!!$    do i = 1, ncells
!!$       read(10,*) cell_vertices(1,i), cell_vertices(2,i), cell_vertices(3,i)
!!$    end do
!!$    close(10)
!!$            
!!$    npoints = 9
!!$    allocate(xpts(3,npoints))  
!!$    xpts(:,1) = [0.0d0, 0.0d0, 0.0d0]
!!$    xpts(:,2) = [0.5d0, 0.0d0, 0.0d0]
!!$    xpts(:,3) = [1.0d0, 0.0d0, 0.0d0]
!!$    xpts(:,4) = [0.0d0, 0.5d0, 0.0d0]
!!$    xpts(:,5) = [0.5d0, 0.5d0, 0.0d0]
!!$    xpts(:,6) = [1.0d0, 0.5d0, 0.0d0]
!!$    xpts(:,7) = [0.0d0, 1.0d0, 0.0d0]
!!$    xpts(:,8) = [0.5d0, 1.0d0, 0.0d0]
!!$    xpts(:,9) = [1.0d0, 1.0d0, 0.0d0]
!!$
!!$    ncells = 8
!!$    allocate(cell_vertices(3,ncells))
!!$    cell_vertices(:,1) = [1,2,4]
!!$    cell_vertices(:,2) = [2,5,4]
!!$    cell_vertices(:,3) = [2,6,5]
!!$    cell_vertices(:,4) = [2,3,6]    
!!$    cell_vertices(:,5) = [4,8,7]
!!$    cell_vertices(:,6) = [4,5,8]
!!$    cell_vertices(:,7) = [5,9,8]
!!$    cell_vertices(:,8) = [5,6,9]
!!$
!!$    call reverse_map(cell_vertices, vertex_cells, num_vertex_cells)
!!$    do i = 1, npoints
!!$       print *, 'node', i, 'cells', vertex_cells(:,i)       
!!$    end do
!!$   
!!$    deallocate(xpts)
!!$    deallocate(cell_vertices)
!!$
!!$  end block load_mesh

!  manual_unstructured_mesh: block        
    
!!$    call get_triangular_test_mesh(&
!!$         & npoints, xpts, &
!!$         & ncells, cell_vertices, num_cell_vertices, &
!!$         & nfaces, face_vertices, num_face_vertices &
!!$         & )

!!$    call get_rectangular_test_mesh(&
!!$         & npoints, xpts, &
!!$         & ncells, cell_vertices, num_cell_vertices, &
!!$         & nfaces, face_vertices, num_face_vertices &
!!$         & )
    
    call get_heterogenous_test_mesh(&
         & npoints, xpts, &
         & ncells, cell_vertices, num_cell_vertices, &
         & nfaces, face_vertices, num_face_vertices &
         & )
    
    !-----------------------------------------------------------------!
    ! Post process information to extract the rest of mapping
    !-----------------------------------------------------------------!
!!$    
!!$    print *, ''
!!$    ! Invert to vertex_faces
!!$    call reverse_map(cell_vertices, num_cell_vertices, vertex_cells, num_vertex_cells)
!!$    do i = 1, size(vertex_cells, dim=2)
!!$       print *, 'vertex', i, 'cells', vertex_cells(1:num_vertex_cells(i),i)
!!$    end do
    
    print *, ''
    ! Invert to vertex_faces
    call reverse_map(face_vertices, num_face_vertices, vertex_faces, num_vertex_faces)
    do i = 1, size(vertex_faces, dim=2)
       print *, 'vertex', i, 'faces', vertex_faces(1:num_vertex_faces(i),i)
    end do

    ! Combine maps to get cell_faces
    call get_cell_faces(cell_vertices, vertex_faces, num_vertex_faces, &
         & cell_faces, num_cell_faces)    
    do icell = 1, ncells
       print *, 'cell', icell, 'faces', cell_faces(1:num_cell_faces(icell),icell)
    end do
    
    print *, ''

    ! Invert cell_faces
    call reverse_map(cell_faces, num_cell_faces, face_cells, num_face_cells)
    do iface = 1, size(face_cells, dim=2)
       print *, 'face', iface, 'cells', face_cells(1:num_face_cells(iface),iface)
    end do

    ! Form boundary faces from faces with 1 boundary
    call get_boundary_faces(num_face_cells, boundary_faces)
    print *, boundary_faces
    
    grid = mesh(xpts, face_vertices, &
         & cell_vertices, num_cell_vertices, &
         & cell_faces, num_cell_faces, face_cells, num_face_cells,&
         &  boundary_faces)
    
!  end block manual_unstructured_mesh
  end subroutine create
  
  !===================================================================!
  ! Determine if the face is a boundary face based on how many
  ! neighbouring cells it has.
  !===================================================================!
  
  pure subroutine get_boundary_faces(num_face_cells, boundary_faces)

    integer, intent(in) :: num_face_cells(:)
    integer, intent(out), allocatable :: boundary_faces(:)
    integer :: iface, nfaces, nbfaces, ctr

    nfaces = size(num_face_cells, dim=1)

    ! Boundary faces are the faces corresponding to just cell
    nbfaces = 0
    do iface = 1, nfaces
       if (num_face_cells(iface) .eq. 1) then
          nbfaces = nbfaces + 1
       end if
    end do

    allocate(boundary_faces(nbfaces))
    ctr = 0 
    do iface = 1, nfaces
       if (num_face_cells(iface) .eq. 1) then
          ctr = ctr + 1
          boundary_faces(ctr) = iface
       end if
    end do

  end subroutine get_boundary_faces
    subroutine get_triangular_test_mesh(npoints, xpts, ncells, &
       & cell_vertices, num_cell_vertices, &
       & nfaces, face_vertices, num_face_vertices)

    integer, intent(out) :: npoints, ncells, nfaces

    real(dp), allocatable, intent(out) :: xpts(:,:)
    integer , allocatable, intent(out) :: cell_vertices(:,:)
    integer , allocatable, intent(out) :: face_vertices(:,:)
    integer , allocatable, intent(out) :: num_cell_vertices(:)
    integer , allocatable, intent(out) :: num_face_vertices(:)
    
    npoints = 9
    allocate(xpts(3,npoints))
    xpts(:,1) = [0.0d0, 0.0d0, 0.0d0]
    xpts(:,2) = [0.5d0, 0.0d0, 0.0d0]
    xpts(:,3) = [1.0d0, 0.0d0, 0.0d0]
    xpts(:,4) = [0.0d0, 0.5d0, 0.0d0]
    xpts(:,5) = [0.5d0, 0.5d0, 0.0d0]
    xpts(:,6) = [1.0d0, 0.5d0, 0.0d0]
    xpts(:,7) = [0.0d0, 1.0d0, 0.0d0]
    xpts(:,8) = [0.5d0, 1.0d0, 0.0d0]
    xpts(:,9) = [1.0d0, 1.0d0, 0.0d0]

    ! All are triangles
    ncells = 8
    allocate(num_cell_vertices(ncells))
    num_cell_vertices(1:ncells) = 3
    allocate(cell_vertices(maxval(num_cell_vertices),ncells))    
    cell_vertices(:,1) = [1,2,4]
    cell_vertices(:,2) = [2,5,4]
    cell_vertices(:,3) = [2,6,5]
    cell_vertices(:,4) = [2,3,6]
    cell_vertices(:,5) = [4,8,7]
    cell_vertices(:,6) = [4,5,8]
    cell_vertices(:,7) = [5,9,8]
    cell_vertices(:,8) = [5,6,9]

    ! Face vertices (load from file if possible)
    nfaces = 16
    allocate(num_face_vertices(nfaces))
    num_face_vertices = 2
    allocate(face_vertices(maxval(num_face_vertices),nfaces))    
    face_vertices(:,1)  = [1,2]
    face_vertices(:,2)  = [2,3]
    face_vertices(:,3)  = [4,5]
    face_vertices(:,4)  = [5,6]
    face_vertices(:,5)  = [7,8]
    face_vertices(:,6)  = [8,9]  
    face_vertices(:,7)  = [1,4]
    face_vertices(:,8)  = [2,5]
    face_vertices(:,9)  = [3,6]
    face_vertices(:,10) = [4,7]
    face_vertices(:,11) = [5,8]
    face_vertices(:,12) = [6,9]
    face_vertices(:,13) = [2,4]
    face_vertices(:,14) = [2,6]
    face_vertices(:,15) = [4,8]
    face_vertices(:,16) = [5,9]

  end subroutine get_triangular_test_mesh

  subroutine get_rectangular_test_mesh(npoints, xpts, ncells, &
              & cell_vertices, num_cell_vertices, &
              & nfaces, face_vertices, num_face_vertices)

    integer, intent(out) :: npoints, ncells, nfaces

    real(dp), allocatable, intent(out) :: xpts(:,:)
    integer , allocatable, intent(out) :: cell_vertices(:,:)
    integer , allocatable, intent(out) :: face_vertices(:,:)
    integer , allocatable, intent(out) :: num_cell_vertices(:)
    integer , allocatable, intent(out) :: num_face_vertices(:)
    
    npoints = 9
    allocate(xpts(3,npoints))
    xpts(:,1) = [0.0d0, 0.0d0, 0.0d0]
    xpts(:,2) = [0.5d0, 0.0d0, 0.0d0]
    xpts(:,3) = [1.0d0, 0.0d0, 0.0d0]
    xpts(:,4) = [0.0d0, 0.5d0, 0.0d0]
    xpts(:,5) = [0.5d0, 0.5d0, 0.0d0]
    xpts(:,6) = [1.0d0, 0.5d0, 0.0d0]
    xpts(:,7) = [0.0d0, 1.0d0, 0.0d0]
    xpts(:,8) = [0.5d0, 1.0d0, 0.0d0]
    xpts(:,9) = [1.0d0, 1.0d0, 0.0d0]
   
    ! All are rectangles
    ncells = 4        
    allocate(num_cell_vertices(ncells))
    num_cell_vertices(1:ncells) = 4    
    allocate(cell_vertices(maxval(num_cell_vertices),ncells))
    cell_vertices(:,1) = [1,2,5,4]
    cell_vertices(:,2) = [2,3,6,5]
    cell_vertices(:,3) = [4,5,8,7]
    cell_vertices(:,4) = [5,6,9,8]
    
    ! Face vertices (load from file if possible)    
    nfaces = 12
    allocate(num_face_vertices(nfaces))
    num_face_vertices = 2
    allocate(face_vertices(maxval(num_face_vertices),nfaces))    
    face_vertices(:,1)  = [1,2]
    face_vertices(:,2)  = [2,3]
    face_vertices(:,3)  = [4,5]
    face_vertices(:,4)  = [5,6]
    face_vertices(:,5)  = [7,8]
    face_vertices(:,6)  = [8,9]  
    face_vertices(:,7)  = [1,4]
    face_vertices(:,8)  = [2,5]
    face_vertices(:,9)  = [3,6]
    face_vertices(:,10) = [4,7]
    face_vertices(:,11) = [5,8]
    face_vertices(:,12) = [6,9]

  end subroutine get_rectangular_test_mesh

  subroutine get_heterogenous_test_mesh(npoints, xpts, ncells, &
       & cell_vertices, num_cell_vertices, &
       & nfaces, face_vertices, num_face_vertices)

    integer, intent(out) :: npoints, ncells, nfaces

    real(dp), allocatable, intent(out) :: xpts(:,:)
    integer , allocatable, intent(out) :: cell_vertices(:,:)
    integer , allocatable, intent(out) :: face_vertices(:,:)
    integer , allocatable, intent(out) :: num_cell_vertices(:)
    integer , allocatable, intent(out) :: num_face_vertices(:)
    
    npoints = 9
    allocate(xpts(3,npoints))
    xpts(:,1) = [0.0d0, 0.0d0, 0.0d0]
    xpts(:,2) = [0.5d0, 0.0d0, 0.0d0]
    xpts(:,3) = [1.0d0, 0.0d0, 0.0d0]
    xpts(:,4) = [0.0d0, 0.5d0, 0.0d0]
    xpts(:,5) = [0.5d0, 0.5d0, 0.0d0]
    xpts(:,6) = [1.0d0, 0.5d0, 0.0d0]
    xpts(:,7) = [0.0d0, 1.0d0, 0.0d0]
    xpts(:,8) = [0.5d0, 1.0d0, 0.0d0]
    xpts(:,9) = [1.0d0, 1.0d0, 0.0d0]

    ! Triangles and rectangles
    ncells = 6
    allocate(num_cell_vertices(ncells))
    num_cell_vertices(1:2) = 4
    num_cell_vertices(3:ncells) = 3
    allocate(cell_vertices(maxval(num_cell_vertices),ncells))
    cell_vertices(:,1) = [1,2,5,4]
    cell_vertices(:,2) = [2,3,6,5]
    cell_vertices(:,3) = [4,8,7]
    cell_vertices(:,4) = [4,5,8]
    cell_vertices(:,5) = [5,9,8]
    cell_vertices(:,6) = [5,6,9]

    ! Face vertices (load from file if possible)
    nfaces = 14
    allocate(num_face_vertices(nfaces))
    num_face_vertices = 2
    allocate(face_vertices(maxval(num_face_vertices),nfaces))    
    face_vertices(:,1)  = [1,2]
    face_vertices(:,2)  = [2,3]
    face_vertices(:,3)  = [4,5]
    face_vertices(:,4)  = [5,6]
    face_vertices(:,5)  = [7,8]
    face_vertices(:,6)  = [8,9]  
    face_vertices(:,7)  = [1,4]
    face_vertices(:,8)  = [2,5]
    face_vertices(:,9)  = [3,6]
    face_vertices(:,10) = [4,7]
    face_vertices(:,11) = [5,8]
    face_vertices(:,12) = [6,9]
    face_vertices(:,13) = [4,8]
    face_vertices(:,14) = [5,9]

  end subroutine get_heterogenous_test_mesh

end module mesh_loader
