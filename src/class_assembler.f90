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
     ! type(mesh), pointer :: grid
     class(mesh)   , allocatable :: grid
     !class(physics), allocatable :: system(:) ! poisson on \Omega, dirichlet on dOmega1 , dirchlet dOmega3 , dirichlet, Neumann dOmega4

     ! Number of state varibles
     integer :: num_state_vars

     ! Number of variables each cell
     integer :: num_variables

     ! Flux vector
     real(dp), allocatable :: phi(:)

     ! Matrix filters
     integer :: DIAGONAL       = 0
     integer :: LOWER_TRIANGLE = -1
     integer :: UPPER_TRIANGLE = 1

   contains

     procedure :: create_vector

     ! Evaluation routines
     procedure :: evaluate_vertex_flux
     procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_source
     procedure :: get_skew_source
     procedure :: get_jacobian
     procedure :: get_transpose_jacobian
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

    ! Number of unknowns in the problem (currently only "T")
    this % num_variables = 1

    ! Determine the number of state variables to solve based on the
    ! mesh. In FVM it is the number of cells present.
    this % num_state_vars = this % grid % num_cells

    ! Allocate the flux vector
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

    this % DIAGONAL       = 0
    this % LOWER_TRIANGLE = -1
    this % UPPER_TRIANGLE = 1

  end function construct

  !===================================================================!
  ! Destructor for file object
  !===================================================================!

  pure subroutine destroy(this)

    type(assembler), intent(inout) :: this

    if (allocated(this % grid)) deallocate(this % grid)
!!$    if(associated(this % grid)) then
!!$       deallocate(this % grid)
!!$       nullify(this % grid)
!!$    end if

    if (allocated(this % phi)) deallocate(this % phi)

  end subroutine destroy

  !===================================================================!
  ! Assemble and return the full jacobian matrix
  !===================================================================!

  subroutine get_jacobian(this, A, filter)

    ! Arguments
    class(assembler)      , intent(in)    :: this
    real(dp), allocatable , intent(out)   :: A(:,:)
    integer , optional    , intent(in)    :: filter

    ! Locals
    real(dp), allocatable :: ex(:)
    integer :: icol

    allocate(A(this % num_state_vars, this % num_state_vars))
    allocate(ex(this % num_state_vars)); ex = 0

    if (present(filter)) then

       ! okay for nonlinear case? the state vectors are required for linearization
       ! Assemble only a part of the matrix (lower(-1), upper(1), or diagonal (0))
       do icol = 1, this % num_state_vars
          ex(icol) = 1.0d0
          call this % get_jacobian_vector_product(A(:,icol), ex, filter)
          ex(icol) = 0.0d0
       end do

    else

       ! Assemble Full Matrix A = L + D + U
       do icol = 1, this % num_state_vars
          ex(icol) = 1.0d0
          call this % get_jacobian_vector_product(A(:,icol),ex)
          ex(icol) = 0.0d0
       end do

    end if

    deallocate(ex)

  end subroutine get_jacobian

  !===================================================================!
  ! Assemble and return the full transpose jacobian matrix
  !===================================================================!

  subroutine get_transpose_jacobian(this, A, filter)

    ! Arguments
    class(assembler)      , intent(in)    :: this
    real(dp), allocatable , intent(out)   :: A(:,:)
    integer , optional    , intent(in)    :: filter

    ! Locals
    real(dp), allocatable :: ex(:)
    integer :: irow

    allocate(A(this % num_state_vars, this % num_state_vars))
    allocate(ex(this % num_state_vars)); ex = 0

    if (present(filter)) then

       ! okay for nonlinear case? the state vectors are required for linearization
       ! Assemble only a part of the matrix (lower(-1), upper(1), or diagonal (0))
       do irow = 1, this % num_state_vars
          ex(irow) = 1.0d0
          call this % get_jacobian_vector_product(A(irow,:), ex, filter)
          ex(irow) = 0.0d0
       end do

    else

       ! Assemble Full Matrix A = L + D + U
       do irow = 1, this % num_state_vars
          ex(irow) = 1.0d0
          call this % get_jacobian_vector_product(A(irow,:),ex)
          ex(irow) = 0.0d0
       end do

    end if

    deallocate(ex)

  end subroutine get_transpose_jacobian

  subroutine get_jacobian_vector_product(this, Aq, q, filter)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: q(:)
    real(dp)         , intent(out)   :: Aq(:)
    integer, optional, intent(in)    :: filter

    ! Finite difference coeff for flux approximation between two cells
    real(dp), parameter  :: alpha(2) = [-1.0_dp, 1.0_dp]

    laplace_normal: block

      integer :: icell, iface
      integer :: ncell, fcells(2)

      ! Loop cells
      loop_cells: do icell = 1 , this % grid % num_cells

         ! Get the faces corresponding to this cell
         associate(faces => this % grid % cell_faces (1:this % grid % num_cell_faces(icell),icell))

           ! Loop faces
           Aq(icell) = 0.0d0

           loop_faces: do iface = 1, this % grid % num_cell_faces(icell)

              associate (&
                   & fdelta => this % grid % face_deltas(faces(iface)), &
                   & farea  => this % grid % face_areas(faces(iface)),  &
                   & num_face_cells => this % grid % num_face_cells(faces(iface)) &
                   & )

                ! Interpolate to get face gammas
                !print *, iface, faces(iface), ftag, fdelta, farea !, !fgamma

                ! Add contribution from internal faces
                domain: if (num_face_cells .eq. 2) then

                   ! Neighbour cell index !? filter here probably?
                   fcells(1:num_face_cells) = this % grid % face_cells(1:num_face_cells,faces(iface))

                   ! Neighbour is the one that has a different cell
                   ! index than current icell
                   if (fcells(1) .eq. icell) then
                      ncell = fcells(2)
                   else
                      ncell = fcells(1)
                   end if

                   ! Interior faces (call tagged physics) (FVM Equation)
                   ! Aq(icell) = Aq(icell) + farea*(x(ncell)-x(icell))/fdelta
                   ! for interior faces: ncell .ne. icell
                   present_filter: if (present(filter)) then

                      !-------------------------------------------------!
                      ! Assemble part of matrix
                      !-------------------------------------------------!

                      apply_filter: if (filter .eq. this % UPPER_TRIANGLE) then

                         ! Activate Upper Trianglular Filter
                         if (ncell .gt. icell) then
                            Aq(icell) = Aq(icell) + farea*(q(ncell))/fdelta
                         end if

                      else if (filter .eq. this % LOWER_TRIANGLE) then

                         ! Activate Lower Trianglular Filter
                         if (ncell .lt. icell) then
                            Aq(icell) = Aq(icell) + farea*(q(ncell))/fdelta
                         end if

                      else if (filter .eq. this % DIAGONAL) then

                         ! Activate DIAGONAL Trianglular Filter (note icell)
                         Aq(icell) = Aq(icell) + farea*(-q(icell))/fdelta

                      end if apply_filter

                   else

                      !-------------------------------------------------!
                      ! Assemble Full of matrix
                      !-------------------------------------------------!

                      Aq(icell) = Aq(icell) + farea*(q(ncell)-q(icell))/fdelta

                   end if present_filter

                else

                   !----------------------------------------------------!
                   ! Diagonal self contributions
                   !----------------------------------------------------!

                   ! Adds more things to diagonal (makes the matrix more
                   ! diagonally dominant) Boundary faces (call boundary
                   ! physics) If the supplied filter is diagonal then add.
                   if (present(filter)) then
                      ! If diagonals are flagged
                      if (filter .eq. this % DIAGONAL) then
                         Aq(icell) = Aq(icell) + farea*(0.0d0 - q(icell))/fdelta
                      end if
                   else
                      ! Add diagonal if filter is not supplied
                      Aq(icell) = Aq(icell) + farea*(0.0d0 - q(icell))/fdelta
                   end if

                   print *, "check on how the boundary faces contribute to jacobian"

                end if domain

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

  subroutine get_transpose_jacobian_vector_product(this, Aq, x)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: x(:)
    real(dp)         , intent(out)   :: Aq(:)
    real(dp), allocatable :: A(:,:)
    integer               :: n

    error stop 'not implemented'

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
    ss = real(0,dp)

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
              domain: if (this % grid % num_cell_faces(gface) .eq. 2) then

                 ! Compute tangent.dot.lvector/delta
                 scale = dot_product(&
                      & this % grid % lvec(1:3,gface), &
                      & this % grid % cell_face_tangents(1:3,1,iface,icell) &
                      & )/this % grid % face_deltas(gface)

                 ! Get the vertices associated with this face
                 associate(&
                      & fvertices => this % grid % face_vertices( &
                      & 1:this % grid % num_face_vertices(gface), gface)&
                      & )

                 ss(icell) = ss(icell) + scale*(phiv(fvertices(2))-phiv(fvertices(1)))

               end associate

            end if domain

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

    ! Homogenous dirichlet boundary conditions
    real(dp) , parameter :: phi_left   = 1.0d0
    real(dp) , parameter :: phi_right  = 0.0d0
    real(dp) , parameter :: phi_top    = 0.0d0
    real(dp) , parameter :: phi_bottom = 0.0d0
    real(dp) , parameter :: phib = 0.0d0

    add_boundary_terms: block

      integer :: icell, iface

      ! Loop cells
      !loop_cells: do icell = 1, this % grid % num_cells
      loop_cells: do concurrent (icell = 1:this % grid % num_cells)

         ! Get the faces corresponding to this cell
         associate( &
              & faces => this % grid % cell_faces &
              & (1:this % grid % num_cell_faces(icell),icell), &
              & highest_tag => maxval(this % grid % tag_numbers) &
              & )

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
              boundary_faces: if (this % grid % num_cell_faces(faces(iface)) .eq. 1) then

                 ! Boundary faces (call boundary physics) (minus as we moved it to rhs)
                 b(icell) = b(icell) + farea*(-phib)/fdelta
                 !print *, icell, "boundary", faces(iface), ftag, fdelta, farea !, !fgamma

              end if boundary_faces

            end associate

         end do loop_faces

       end associate

      end do loop_cells

    end block add_boundary_terms

    cell_source: block

      integer :: icell

      loop_cells: do concurrent (icell = 1 : this % grid % num_cells)

         associate(&
              & x => this % grid % cell_centers(:,icell), &
              & cell_volume => this % grid % cell_volumes(icell)&
              & )

           b(icell) = b(icell) + evaluate_source(x)*cell_volume

         end associate

      end do loop_cells

    end block cell_source

  end subroutine get_source

  pure type(real(dp)) function evaluate_source(x)

    real(dp), intent(in) :: x(3)

    evaluate_source = sin(x(1)) !1.0d0 !-1.0d0 !sin(x(1)) + cos(x(2))

  end function evaluate_source

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine write_solution(this, filename, phic)

    use class_paraview_writer, only  : paraview_writer

    class(assembler), intent(in)  :: this
    character(len=*), intent(in)  :: filename
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
    real(dp)        , intent(in)  :: phic(:)
    integer                       ::  i, ierr
    real(dp)        , allocatable :: phiv(:)

    class(paraview_writer), allocatable :: pwriter

    allocate(pwriter, source = paraview_writer(this % grid))

    call pwriter % write('paraview.vtu'//char(0))

    if (allocated(pwriter)) deallocate(pwriter)

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
    ! Write Triangles/Quads (works only for homogeneous elements)
    !-----------------------------------------------------------------!

    if ( maxval(this % grid % num_cell_vertices) .eq. 4 ) then
       write(90, *) &
            & 'ZONE T="Temperature", N=', this % grid % num_vertices, &
            & ', E=', this % grid % num_cells, &
            & ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
    else
       write(90, *) &
            & 'ZONE T="Temperature", N=', this % grid % num_vertices, &
            & ', E=', this % grid % num_cells, &
            & ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
    end if

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

  !===================================================================!
  ! Create a state vector and sets values if a scalar is supplied
  !===================================================================!

  subroutine create_vector(this, x, scalar)

    class(assembler), intent(in)               :: this
    real(dp)        , intent(out), allocatable :: x(:)
    real(dp)        , intent(in) , optional    :: scalar

    if (allocated(x)) error stop "vector already allocated"
    allocate(x(this % num_state_vars))
    if (present(scalar))  x = scalar

  end subroutine create_vector

end module class_assembler
