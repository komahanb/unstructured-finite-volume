module class_assembler

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  use class_mesh, only : mesh
  
  implicit none

  private
  public :: assembler
  public :: solve_conjugate_gradient

  !===================================================================!
  ! Class responsible for matrix, right hand side assembly and boundary
  ! conditions
  !===================================================================!

  type :: assembler

     ! set symmetry to .true. for structured grid
     logical :: symmetry = .true.

     ! Mesh object
     type(mesh), pointer :: grid

     ! Number of state varibles 
     integer :: num_state_vars

     ! Flux vector
     real(dp), allocatable :: phi(:)

   contains

     ! Evaluation routines
     !procedure :: evaluate_vertex_flux
     !procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_source
     !procedure :: add_skew_source
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
    
    laplace: block
      
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

              else

                 ! Boundary faces (call boundary physics)
                 Ax(icell) = Ax(icell) + farea*(0.0d0 - x(icell))/fdelta
                 !print *, icell, "boundary", faces(iface), ftag, fdelta, farea !, !fgamma

              end if

            end associate

           end do loop_faces

         end associate

      end do loop_cells
      
    end block laplace
    
  end subroutine get_jacobian_vector_product
  
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
  
  subroutine get_source(this, b)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: b(:)

    real(dp) , parameter :: phib = -1.0d0
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

    evaluate_source = 0.0d0 !sin(x(1)) + cos(x(2))

  end function evaluate_source

  subroutine solve_conjugate_gradient(oassembler, max_it, max_tol, x)

    ! Arguments
    type(assembler)  , intent(in)    :: oassembler
    integer          , intent(in)    :: max_it
    real(dp)         , intent(in)    :: max_tol
    real(dp)         , intent(inout) :: x(:)

    ! Create local data
    real(dp), allocatable :: p(:), r(:), w(:), Ax(:), tmp(:)
    real(dp), allocatable :: b(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    integer               :: iter
    real(dp)              :: rho(2)

    ! Start the iteration counter
    iter = 1

    ! Memory allocations
    allocate(b,p,r,w,Ax,tmp,mold=x)

    ! Norm of the right hand side
    call oassembler % get_source(tmp)
    if (oassembler % symmetry .eqv. .false.) then
       call oassembler % get_transpose_jacobian_vector_product(b, tmp)
    else
       b = tmp
    end if
    bnorm = norm2(b)

    ! Norm of the initial residual
    call oassembler % get_jacobian_vector_product(tmp, x)
    if (oassembler % symmetry .eqv. .false.) then
       call oassembler % get_transpose_jacobian_vector_product(Ax, tmp)
    else
       Ax = tmp
    end if
    r         = b - Ax ! could directly form this residual using get_residual_call
    rnorm     = norm2(r)
    tol       = rnorm/bnorm
    rho(2)    = rnorm*rnorm
    
    open(13, file='cg.log', action='write', position='append')

    tol = huge(1.0_dp)

    ! Apply Iterative scheme until tolerance is achieved
    do while ((tol .gt. max_tol) .and. (iter .lt. max_it))

       ! step (a) compute the descent direction
       if ( iter .eq. 1) then
          ! steepest descent direction p
          p = r
       else
          ! take a conjugate direction
          beta = rho(2)/rho(1)
          p = r + beta*p
       end if

       ! step (b) compute the solution update
       call oassembler % get_jacobian_vector_product(tmp, p)
       if (oassembler % symmetry .eqv. .false.) then
          call oassembler % get_transpose_jacobian_vector_product(w, tmp)
       else
          w = tmp
       end if

       ! step (c) compute the step size for update
       alpha = rho(2)/dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual
       r = r - alpha*w

       ! step(f) update values before next iteration
       rnorm = norm2(r)
       tol = rnorm/bnorm

       write(13,*) iter, tol
       write(*,*) iter, tol, rnorm, rho ! causes valgrind errors

       iter = iter + 1

       rho(1) = rho(2)
       rho(2) = rnorm*rnorm

    end do

    close(13)

    deallocate(r, p, w, b, Ax, tmp)

  end subroutine solve_conjugate_gradient
  
  !===================================================================!
  ! Write solution to file
  !===================================================================!
  
  subroutine write_solution(this, filename)!, phi)

    class(assembler), intent(in)  :: this
    character(len=*), intent(in)  :: filename
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
!    real(dp)        , intent(in)  :: phi(:)
    integer                       ::  i, ierr

    ! Open resource
    path = trim(filename)

    open(unit=90, file=path, iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    ! Write header
    write(90, *) 'TITLE = "FVM-Laplace"'
    write(90, *) 'VARIABLES = "x" "y"'
    write(90, *) 'ZONE T="Temperature", N=', this % grid % num_vertices, &
         & ', E=', this % grid % num_cells, &
         & ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
    
    ! Write vertices
    do i = 1, this % grid % num_vertices
       write(90,*) this % grid % vertices(1:2,i)!, phi(i)
    end do
    
    ! Write cell connectivities
    do i = 1, this % grid % num_cells
       write(90,*) this % grid % cell_vertices(:,i)
    end do

    ! Close resource
    close(90)
    
    if (allocated(path)) deallocate(path)
    if (allocated(new_name)) deallocate(new_name)

  end subroutine write_solution

end module class_assembler
