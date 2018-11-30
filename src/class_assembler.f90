module class_assembler

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  use class_mesh, only : mesh
  
  implicit none

  private
  public :: assembler

  !===================================================================!
  ! Class responsible for matrix, right hand side assembly and boundary
  ! conditions
  !===================================================================!

  type :: assembler

     ! set symmetry to .true. for structured grid
     logical :: symmetry = .true.

     ! Mesh object
     class(mesh), pointer :: grid

     ! Number of state varibles 
     integer :: num_state_vars

     ! Flux vector
     real(dp), allocatable :: phi(:)

   contains

     ! Evaluation routines
     procedure :: evaluate_vertex_flux
     procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_source
     procedure :: get_skew_source
     !procedure :: get_jacobian
     procedure :: get_jacobian_vector_product
     procedure :: get_transpose_jacobian_vector_product
     procedure :: write_solution
     
     ! Destructor
     final :: destroy

  end type assembler

  interface assembler
     module procedure construct
  end interface assembler

contains
    
  !===================================================================!
  ! Constructor for physics
  !===================================================================!
  
  type(assembler) function construct(grid) result (this)

    type(mesh), intent(in) :: grid

    print *, "constructing assembler"

    ! Set mesh
    allocate(this % grid, source  = grid)
    call this % grid % to_string()

    ! Non symmetric jacobian
    this % symmetry = .true.

    ! Determine the number of state variables to solve based on the
    ! mesh. In FVM it is the number of cells present.
    this % num_state_vars = this % grid % num_cells

    ! Allocate the flux vector
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

  end function construct

  !===================================================================!
  ! Destructor for file object
  !===================================================================!
  
  pure subroutine destroy(this)
    
    type(assembler), intent(inout) :: this
    
    if(associated(this % grid)) then
       deallocate(this % grid)
       nullify(this % grid)
    end if

    if (allocated(this % phi)) deallocate(this % phi)
    
  end subroutine destroy

!!$  subroutine get_jacobian(this, A, x)
!!$
!!$    class(assembler) , intent(in)    :: this
!!$    real(dp)         , intent(in)    :: x(:)
!!$    real(dp)         , intent(out)   :: A(:,:)
!!$
!!$    ! okay for nonlinear case?
!!$    matdim = size(x)
!!$    x = 
!!$    do icol = 1, matdim
!!$       call this % get_jacobian_vector_product(A(icol,:),x)
!!$    end do
!!$
!!$  end subroutine get_jacobian

  subroutine get_jacobian_vector_product(this, Ax, x)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(out)   :: Ax(:)

!!$    test: block 
!!$      integer               :: n
!!$      real(dp), allocatable :: A(:,:)
!!$      n = size(x,dim=1)
!!$      allocate(A(n,n))
!!$      A(:,1) = [-6.0d0,1.0d0,1.0d0,0.0d0]
!!$      A(:,2) = [1.0d0,-6.0d0,0.0d0,1.0d0]
!!$      A(:,3) = [1.0d0,0.0d0,-6.0d0,1.0d0]
!!$      A(:,4) = [0.0d0,1.0d0,1.0d0,-6.0d0]
!!$
!!$      Ax = matmul(A,x)
!!$      deallocate(A)
!!$    end block test
!!$
!!$    print *, "Ax=",Ax
!!$    
!!$    return
    
    laplace_normal: block
      
      integer :: icell, iface
      integer :: ncell, fcells(2)

      ! Loop cells
      !loop_cells: do icell = 1, this % grid % num_cells
      loop_cells: do concurrent (icell = 1 : this % grid % num_cells)

         ! Get the faces corresponding to this cell
         associate( &
              & faces => this % grid % cell_faces &
              & (1:this % grid % num_cell_faces(icell),icell), &
              & highest_tag => maxval(this % grid % tag_numbers) &
              & )

           !print *, icell, faces, highest_tag
           
           ! Loop faces
           Ax(icell) = 0.0d0

           loop_faces: do iface = 1, this % grid % num_cell_faces(icell)
              
              associate (&
                   & ftag   => this % grid % face_tags(faces(iface))  , &
                   & fdelta => this % grid % face_deltas(faces(iface)), &
                   & farea  => this % grid % face_areas(faces(iface)),  &
                   & nfcells => this % grid % num_face_cells(faces(iface)) &
                   & )

                ! Interpolate to get face gammas
                !print *, iface, faces(iface), ftag, fdelta, farea !, !fgamma
                
              ! Add contribution from internal faces
              if (ftag .eq. highest_tag) then

                 ! Neighbour cell index
                 fcells(1:nfcells) = this % grid % face_cells(1:nfcells,faces(iface))

                 ! Neighbour is the one that has a different cell
                 ! index than current icell
                 if (fcells(1) .eq. icell) then 
                    ncell = fcells(2)
                 else 
                    ncell = fcells(1)
                 end if

                 !print *, "cell=",icell, "face=", faces(iface), "ncell=", ncell

                 ! Interior faces (call tagged physics) (FVM Equation)
                 Ax(icell) = Ax(icell) + farea*(x(ncell) - x(icell))/fdelta

                 !print *, icell, "internal", faces(iface), ftag, fdelta, farea !, !fgamma

!!$              else if (ftag .eq. 1) then
!!$
!!$                 Ax(icell) = Ax(icell) + farea*1.0d0
                 
              else

                 ! Boundary faces (call boundary physics)
                 Ax(icell) = Ax(icell) + farea*(0.0d0 - x(icell))/fdelta
                 !print *, icell, "boundary", faces(iface), ftag, fdelta, farea !, !fgamma

              end if

            end associate

           end do loop_faces

         end associate

      end do loop_cells
      
    end block laplace_normal

  end subroutine get_jacobian_vector_product

  !===================================================================!
  ! Compute vertex values by interpolating cell center values
  !===================================================================!

  subroutine evaluate_vertex_flux(this, phiv, phic)

    class(assembler) , intent(in)               :: this
    real(dp)         , intent(in)               :: phic(:)
    real(dp)         , intent(out), allocatable :: phiv(:)

    integer :: ivertex

    ! Interpolate the supplied cell centered solution to form the
    ! nodal solution
    allocate(phiv(this % grid % num_vertices)); phiv = 0
    do concurrent (ivertex = 1: this % grid % num_vertices)
       associate(&
            & w => this % grid % vertex_cell_weights(&
            & 1:this % grid % num_vertex_cells(ivertex), ivertex&
            & ), &
            & icells => this % grid % vertex_cells(&
            & 1:this % grid % num_vertex_cells(ivertex), ivertex)&
            & )
         phiv(ivertex) = dot_product(phic(icells), w)
       end associate
    end do

  end subroutine evaluate_vertex_flux

  !===================================================================!
  ! Compute face center values by interpolating cell center values
  !===================================================================!

  subroutine evaluate_face_flux(this, phif, phic)

    class(assembler) , intent(in)               :: this
    real(dp)         , intent(in)               :: phic(:)
    real(dp)         , intent(out), allocatable :: phif(:)

    integer :: iface

    ! Interpolate the supplied cell centered solution to form the
    ! face center solution
    allocate(phif(this % grid % num_faces)); phif = 0
    do concurrent (iface = 1: this % grid % num_faces)
       associate(&
            & w => this % grid % face_cell_weights(&
            & 1:this % grid % num_face_cells(iface), iface&
            & ), &
            & icells => this % grid % face_cells(&
            & 1:this % grid % num_face_cells(iface), iface)&
            & )
         phif(iface) = dot_product(phic(icells), w)
       end associate
    end do

  end subroutine evaluate_face_flux

  subroutine get_transpose_jacobian_vector_product(this, Ax, x)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(out)   :: Ax(:)
    real(dp), allocatable :: A(:,:)
    integer               :: n

!!$    n = size(x)
!!$    allocate(A(n,n))
!!$    
!!$    A(:,1) = [-6.0d0,1.0d0,1.0d0,0.0d0]
!!$    A(:,2) = [1.0d0,-6.0d0,0.0d0,1.0d0]
!!$    A(:,3) = [1.0d0,0.0d0,-6.0d0,1.0d0]
!!$    A(:,4) = [0.0d0,1.0d0,1.0d0,-6.0d0]
!!$
!!$    Ax = matmul(transpose(A),x)
!!$
!!$    deallocate(A)

  end subroutine get_transpose_jacobian_vector_product
  
  !===================================================================!
  ! Evaluate internal skew source based on the current cell states and
  ! return
  !===================================================================!
  
  !subroutine get_tangential_flux(this, ss, phic)  
  subroutine get_skew_source(this, ss, phic)

    class(assembler), intent(in)  :: this
    type(real(dp))  , intent(in)  :: phic(:)
    type(real(dp))  , intent(out) :: ss(:)

    ! Evaluation procedure
    evaluate: block
    
    ! Local variables
    real(8) :: scale
    integer :: icell, iface, gface
    type(real(dp)) , allocatable :: phiv(:)

    ! Evaluate nodal values of phi
    call this % evaluate_vertex_flux(phiv, phic)

    ! Make space for skew source terms
    ss = 0;

    loop_cells: do icell = 1, this % grid % num_cells

       ! Get the faces corresponding to this cell
       associate( &
            & faces => this % grid % cell_faces &
            & (1:this % grid % num_cell_faces(icell),icell), &
            & highest_tag => maxval(this % grid % tag_numbers) &
            & )

           loop_faces : do iface = 1, this % grid % num_cell_faces(icell) 

              ! Global face number
              gface = faces(iface)
                 
              ! Ignore boundary faces from skew source evaluation
              if (this % grid % face_tags(gface) .eq. highest_tag) then

                 ! Compute tangent.dot.lvector/delta
                 scale = dot_product(&
                      & this % grid % lvec(1:3,gface), &
                      & this % grid % cell_face_tangents(1:3,iface,icell) &
                      & )/this % grid % face_deltas(gface)

                 ! Get the vertices associated with this face
                 associate(&
                      & fvertices => this % grid % face_vertices( &
                      & 1:this % grid % num_face_vertices(gface), gface)&
                      & )

                 ss(icell) = ss(icell) + scale*(phiv(fvertices(2))-phiv(fvertices(1)))

               end associate
              
              end if

           end do loop_faces
           
         end associate
         
      end do loop_cells

      ! Deallocate the vertex flux values
      deallocate(phiv)

    end block evaluate

  end subroutine get_skew_source
    
  subroutine get_source(this, b)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: b(:)

    real(dp) , parameter :: phib = 0.0d0
!!$    
!!$    block
!!$      b(1) = 4.0d-1
!!$      b(2) = 4.0d-1
!!$      b(3) = 4.0d-1
!!$      b(4) = 4.0d-1
!!$    end block
!!$
!!$    print *, 'source', b
!!$    return
    
    add_boundary_terms: block

      integer :: icell, iface
      integer :: ncell, fcells(2)
      
      ! Loop cells
      !loop_cells: do icell = 1, this % grid % num_cells
      loop_cells: do concurrent (icell = 1:this % grid % num_cells)

         ! Get the faces corresponding to this cell
         associate( &
              & faces => this % grid % cell_faces &
              & (1:this % grid % num_cell_faces(icell),icell), &
              & highest_tag => maxval(this % grid % tag_numbers) &
              & )
           
           !print *, icell, faces, highest_tag

           ! Loop faces
           b(icell) = 0.0d0
         
         loop_faces: do iface = 1, this % grid % num_cell_faces(icell)

            associate (&
                 & ftag   => this % grid % face_tags(faces(iface))  , &
                 & fdelta => this % grid % face_deltas(faces(iface)), &
                 & farea  => this % grid % face_areas(faces(iface))  &
                 & )
              
              ! Interpolate to get face gammas
              !print *, iface, faces(iface), ftag, fdelta, farea !, !fgamma
              
              ! Add contribution from internal faces
              if (ftag .ne. highest_tag) then ! homogenous dirichlet T = 1.0d0
               
               ! Boundary faces (call boundary physics)
                b(icell) = b(icell) + farea*(-phib)/fdelta
                !print *, icell, "boundary", faces(iface), ftag, fdelta, farea !, !fgamma
               
              end if
            
            end associate

         end do loop_faces
         
       end associate

      end do loop_cells
    
    end block add_boundary_terms

    
    cell_source: block

      integer :: icell

      ! Loop cells
      loop_cells: do concurrent (icell = 1 : this % grid % num_cells)
         associate( &
              & x => this % grid % cell_centers(:,icell), &
              & cell_volume => this % grid % cell_volumes(icell))
           b(icell) = b(icell) + evaluate_source(x)*cell_volume
         end associate
      end do loop_cells

    end block cell_source

  end subroutine get_source
  
  pure type(real(dp)) function evaluate_source(x)

    real(dp), intent(in) :: x(3)

    evaluate_source = -1.0d0 !sin(x(1)) + cos(x(2))

  end function evaluate_source
  
  !===================================================================!
  ! Write solution to file
  !===================================================================!
  
  subroutine write_solution(this, filename, phic)

    class(assembler), intent(in)  :: this
    character(len=*), intent(in)  :: filename
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
    real(dp)        , intent(in)  :: phic(:)
    integer                       ::  i, ierr

    real(dp), allocatable :: phiv(:)

    ! Open resource
    path = trim(filename)

    open(unit=90, file=path, iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    ! Compute vertex values by interpolating cell center values
    call this % evaluate_vertex_flux(phiv, phic)
    
    ! Write header
    write(90, *) 'TITLE = "FVM-Laplace"'
    write(90, *) 'VARIABLES = "x" "y"  "T"'

    !-----------------------------------------------------------------!
    ! Write Triangles
    !-----------------------------------------------------------------!

    write(90, *) 'ZONE T="Temperature", N=', this % grid % num_vertices, &
         & ', E=', this % grid % num_cells, &
         & ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
    
    ! Write vertices
    do i = 1, this % grid % num_vertices
       write(90,*) this % grid % vertices(1:2,i), phiv(i)
    end do
    
    ! Write cell connectivities
    do i = 1, this % grid % num_cells
       write(90,*) this % grid % cell_vertices(1:this % grid % num_cell_vertices(i),i)
    end do

    !-----------------------------------------------------------------!
    ! Write Other types of elements too
    !-----------------------------------------------------------------!

    ! Close resource
    close(90)
    
    if (allocated(path)) deallocate(path)
    if (allocated(new_name)) deallocate(new_name)
    if (allocated(phiv)) deallocate(phiv)

  end subroutine write_solution

end module class_assembler
