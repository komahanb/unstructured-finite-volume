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

     type(mesh), pointer :: grid
     logical :: symmetry = .false.

   contains

     ! Evaluation routines
     !procedure :: evaluate_vertex_flux
     !procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_source
     !procedure :: add_skew_source
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

    ! Non symmetric jacobian
    this % symmetry = .false.

    call this % grid % to_string()

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

  end subroutine destroy

  subroutine get_jacobian_vector_product(this, x, Ax)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(out)   :: Ax(:)

    integer               :: n
    real(dp), allocatable :: A(:,:)

    n = size(x,dim=1)
    allocate(A(n,n))

    A(1,:) = [2.0d0,1.0d0,1.0d0]
    A(2,:) = [1.0d0,2.0d0,1.0d0]
    A(3,:) = [0.0d0,1.0d0,2.0d0]

    Ax = matmul(A,x)

    deallocate(A)

  end subroutine get_jacobian_vector_product
  
  subroutine get_transpose_jacobian_vector_product(this, x, Ax)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(out)   :: Ax(:)
    real(dp), allocatable :: A(:,:)
    integer               :: n

    n = size(x)
    allocate(A(n,n))

    A(1,:) = [2.0d0,1.0d0,1.0d0]
    A(2,:) = [1.0d0,2.0d0,1.0d0]
    A(3,:) = [0.0d0,1.0d0,2.0d0]

    Ax = matmul(transpose(A),x)

    deallocate(A)

  end subroutine get_transpose_jacobian_vector_product
  
  subroutine get_source(this, b)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: b(:)

    b(1) = 1.0d0
    b(2) = 1.0d0
    b(3) = 1.0d0

  end subroutine get_source

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
       call oassembler % get_transpose_jacobian_vector_product(tmp, b)
    end if
    bnorm = norm2(b)

    ! Norm of the initial residual
    call oassembler % get_jacobian_vector_product(x, tmp)
    if (oassembler % symmetry .eqv. .false.) then
       call oassembler % get_transpose_jacobian_vector_product(tmp, Ax)
    end if
    r         = b - Ax
    rnorm     = norm2(r)
    tol       = rnorm/bnorm
    rho(2)    = rnorm*rnorm
    
    open(10, file='cg.log', action='write', position='append')

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
       call oassembler % get_jacobian_vector_product(p, tmp)
       if (oassembler % symmetry .eqv. .false.) then
          call oassembler % get_transpose_jacobian_vector_product(tmp, w)
       end if
       !w = matmul(A,p)

       ! step (c) compute the step size for update
       alpha = rho(2)/dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual
       r = r - alpha*w

       ! step(f) update values before next iteration
       rnorm = norm2(r)
       tol = rnorm/bnorm

       write(10,*) iter, tol
       print *, iter, tol, rnorm, rho

       iter = iter + 1

       rho(1) = rho(2)
       rho(2) = rnorm*rnorm

    end do

    close(10)

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
