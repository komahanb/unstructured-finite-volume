module class_assembler

  ! import dependencies
  use iso_fortran_env         , only : dp => REAL64
  use class_mesh              , only : mesh
  use class_string            , only : string
  use class_boundary_condition, only : boundary_condition, dirichlet, neumann, robin
  use interface_equation      , only : equation
  use class_diffusion         , only : diffusion
  use class_graph             , only : graph

  implicit none

  private
  public :: assembler

  !===================================================================!
  ! Class responsible for matrix, right hand side assembly and boundary
  ! conditions
  !===================================================================!

  type :: assembler

     ! Mesh object
     ! type(mesh), pointer :: grid
     class(mesh)   , allocatable :: grid
     !class(physics), allocatable :: system(:) ! poisson on \Omega, dirichlet on dOmega1 , dirchlet dOmega3 , dirichlet, Neumann dOmega4

     ! Number of state varibles
     integer :: num_state_vars

     ! Number of variables each cell
     integer :: num_variables

     ! Dof / connectivity graph (cells = vertices, interior faces = edges)
     type(graph) :: g

     ! Flux vector
     real(dp), allocatable :: phi(:)

     ! Boundary conditions - one per (face tag, variable), addressed by name
     type(boundary_condition), allocatable :: bcs(:,:)

     ! The pde being discretized (diffusion tensor, source, #variables)
     class(equation), allocatable :: model

     ! Transient (backward euler) support. The operator becomes
     ! (M/dt - A) and the rhs (M/dt) phi_old - b, with M = cell volume.
     logical               :: transient = .false.
     real(dp)              :: dt = 0.0_dp
     real(dp), allocatable :: phi_old(:)

     ! Matrix filters
     integer :: DIAGONAL       = 0
     integer :: LOWER_TRIANGLE = -1
     integer :: UPPER_TRIANGLE = 1

   contains

     procedure :: create_vector

     ! Boundary conditions - addressed by physical group name
     procedure :: set_dirichlet
     procedure :: set_neumann
     procedure :: set_robin
     procedure :: set_dirichlet_tag

     ! The pde
     procedure :: set_equation

     ! Transient marching (backward euler)
     procedure :: set_transient

     ! Evaluation routines
     procedure :: evaluate_vertex_flux
     procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_source
     procedure :: get_skew_source
     procedure :: get_jacobian
     procedure :: get_transpose_jacobian
     procedure :: get_jacobian_vector_product
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
    ! One variable per cell by default ("T"); set_equation bumps this up
    this % num_variables = 1

    ! Build the dof/connectivity graph and size the system from it. For a
    ! single field num_state_vars = num_cells, as before.
    this % g = graph(grid, this % num_variables)
    this % num_state_vars = this % g % num_dofs()

    ! Allocate the flux vector
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

    ! Default every (boundary tag, variable) to homogeneous dirichlet
    ! (phi=0). The driver overrides by name. Interior faces never look
    ! this table up.
    allocate(this % bcs(maxval(this % grid % tag_numbers), this % num_variables))

    ! Default physics: isotropic unit diffusion, no source - recovers the
    ! old laplace operator. The driver overrides with set_equation.
    allocate(this % model, source = diffusion(1.0_dp))

    this % DIAGONAL       = 0
    this % LOWER_TRIANGLE = -1
    this % UPPER_TRIANGLE = 1

  end function construct

  !===================================================================!
  ! Set a dirichlet condition on a named physical group. The optional
  ! var targets a single variable; omitted, it applies to all variables.
  !===================================================================!

  subroutine set_dirichlet(this, name, value, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: value
    integer, optional, intent(in)    :: var

    call set_bc(this, this % grid % find_tag_by_name(name), name, dirichlet(value), var)

  end subroutine set_dirichlet

  !===================================================================!
  ! Set a neumann condition on a named physical group
  !===================================================================!

  subroutine set_neumann(this, name, flux, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: flux
    integer, optional, intent(in)    :: var

    call set_bc(this, this % grid % find_tag_by_name(name), name, neumann(flux), var)

  end subroutine set_neumann

  !===================================================================!
  ! Set a robin condition  a*phi + b*dphi/dn = c  on a named group
  !===================================================================!

  subroutine set_robin(this, name, a, b, c, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: a, b, c
    integer, optional, intent(in)    :: var

    call set_bc(this, this % grid % find_tag_by_name(name), name, robin(a,b,c), var)

  end subroutine set_robin

  !===================================================================!
  ! Convenience for replicating tag-keyed setups without names
  !===================================================================!

  subroutine set_dirichlet_tag(this, tag, value, var)

    class(assembler) , intent(inout) :: this
    integer          , intent(in)    :: tag
    real(dp)         , intent(in)    :: value
    integer, optional, intent(in)    :: var

    call set_bc(this, tag, "", dirichlet(value), var)

  end subroutine set_dirichlet_tag

  !===================================================================!
  ! Set the pde being solved (diffusion tensor, source, #variables).
  ! Call this BEFORE setting boundary conditions - changing the number
  ! of variables rebuilds the (tag,variable) bc table.
  !===================================================================!

  subroutine set_equation(this, model)

    class(assembler), intent(inout) :: this
    class(equation) , intent(in)    :: model

    ! Guard the ordering footgun: set_equation rebuilds the (tag,variable)
    ! bc table, so any boundary conditions set earlier would be lost.
    if (allocated(this % bcs)) then
       if (any(this % bcs % tag .ne. -1)) then
          write(*,*) "warning: set_equation called after boundary conditions; ", &
               & "they have been reset - call set_equation first"
       end if
    end if

    if (allocated(this % model)) deallocate(this % model)
    allocate(this % model, source = model)

    ! Resize the system for this many variables per cell
    this % num_variables     = model % num_variables
    this % g % num_variables = this % num_variables
    this % num_state_vars    = this % g % num_dofs()

    if (allocated(this % phi)) deallocate(this % phi)
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

    if (allocated(this % bcs)) deallocate(this % bcs)
    allocate(this % bcs(maxval(this % grid % tag_numbers), this % num_variables))

    if (this % transient) then
       if (allocated(this % phi_old)) deallocate(this % phi_old)
       allocate(this % phi_old(this % num_state_vars))
       this % phi_old = 0.0_dp
    end if

  end subroutine set_equation

  !===================================================================!
  ! Enable backward-euler transient marching with step dt. The
  ! time_integrator advances phi_old between solves.
  !===================================================================!

  subroutine set_transient(this, dt)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(in)    :: dt

    this % transient = .true.
    this % dt = dt

    if (allocated(this % phi_old)) deallocate(this % phi_old)
    allocate(this % phi_old(this % num_state_vars))
    this % phi_old = 0.0_dp

  end subroutine set_transient

  !===================================================================!
  ! Stamp a resolved bc into the table (module-private helper). With no
  ! variable index the bc applies to every variable on the boundary.
  !===================================================================!

  subroutine set_bc(this, tag, name, bc, var)

    class(assembler)        , intent(inout) :: this
    integer                 , intent(in)    :: tag
    character(len=*)        , intent(in)    :: name
    type(boundary_condition), intent(in)    :: bc
    integer, optional       , intent(in)    :: var

    integer :: v, vlo, vhi

    if (tag .lt. 1 .or. tag .gt. size(this % bcs, 1)) then
       write(*,*) "set_bc: unknown boundary '", name, "' -> tag ", tag
       error stop
    end if

    vlo = 1
    vhi = this % num_variables

    if (present(var)) then
       vlo = var
       vhi = var
    end if

    do v = vlo, vhi
       this % bcs(tag, v)        = bc
       this % bcs(tag, v) % tag  = tag
       this % bcs(tag, v) % name = string(name)
    end do

  end subroutine set_bc

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

      integer  :: icell, iface, ivar, p, n
      integer  :: ncell, fcells(2)
      real(dp) :: nf(3), Kf(3,3), keff

      associate(nv => this % g % num_variables)

      ! Loop cells
      loop_cells: do icell = 1 , this % grid % num_cells

         ! Get the faces corresponding to this cell
         associate(faces => this % grid % cell_faces (1:this % grid % num_cell_faces(icell),icell))

           ! Zero this cell's dofs
           do ivar = 1, nv
              Aq(this % g % dof(icell,ivar)) = 0.0d0
           end do

           loop_faces: do iface = 1, this % grid % num_cell_faces(icell)

              associate (&
                   & fdelta => this % grid % face_deltas(faces(iface)), &
                   & farea  => this % grid % face_areas(faces(iface)),  &
                   & ftag   => this % grid % face_tags(faces(iface)),   &
                   & num_face_cells => this % grid % num_face_cells(faces(iface)) &
                   & )

                ! Effective normal conductivity from the diffusion tensor
                nf   = this % grid % cell_face_normals(1:3, iface, icell)
                Kf   = this % model % diffusion_tensor(this % grid % face_centers(1:3, faces(iface)))
                keff = dot_product(nf, matmul(Kf, nf))

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

                   ! Interior face: couple each variable to the same
                   ! variable in the neighbour cell (decoupled fields).
                   ! The upper/lower split is by cell index, which for a
                   ! fixed variable matches the dof ordering.
                   variables_interior: do ivar = 1, nv

                      p = this % g % dof(icell, ivar)
                      n = this % g % dof(ncell, ivar)

                      present_filter: if (present(filter)) then

                         apply_filter: if (filter .eq. this % UPPER_TRIANGLE) then
                            if (ncell .gt. icell) Aq(p) = Aq(p) + keff*farea*(q(n))/fdelta
                         else if (filter .eq. this % LOWER_TRIANGLE) then
                            if (ncell .lt. icell) Aq(p) = Aq(p) + keff*farea*(q(n))/fdelta
                         else if (filter .eq. this % DIAGONAL) then
                            Aq(p) = Aq(p) + keff*farea*(-q(p))/fdelta
                         end if apply_filter

                      else

                         Aq(p) = Aq(p) + keff*farea*(q(n)-q(p))/fdelta

                      end if present_filter

                   end do variables_interior

                else

                   !----------------------------------------------------!
                   ! Boundary face - the bc contributes to the diagonal
                   !----------------------------------------------------!

                   ! Robin-type bc: flux = lhs_coeff*phi_p + const. Only
                   ! the phi_p (diagonal) part lives in the operator; the
                   ! constant goes to the rhs in get_source. A neumann bc
                   ! has lhs_coeff = 0, so it adds nothing here.
                   variables_boundary: do ivar = 1, nv

                      p = this % g % dof(icell, ivar)

                      if (present(filter)) then
                         if (filter .eq. this % DIAGONAL) then
                            Aq(p) = Aq(p) &
                                 & + this % bcs(ftag,ivar) % lhs_coeff(farea, fdelta, keff)*q(p)
                         end if
                      else
                         Aq(p) = Aq(p) &
                              & + this % bcs(ftag,ivar) % lhs_coeff(farea, fdelta, keff)*q(p)
                      end if

                   end do variables_boundary

                end if domain

              end associate

           end do loop_faces

           ! Transient: the operator is (M/dt - A). Negate the spatial
           ! part assembled above and add the mass term M/dt to the
           ! diagonal (full and DIAGONAL filter only).
           if (this % transient) then
              do ivar = 1, nv
                 p = this % g % dof(icell, ivar)
                 Aq(p) = -Aq(p)
                 if (present(filter)) then
                    if (filter .eq. this % DIAGONAL) &
                         & Aq(p) = Aq(p) + this % grid % cell_volumes(icell)/this % dt*q(p)
                 else
                    Aq(p) = Aq(p) + this % grid % cell_volumes(icell)/this % dt*q(p)
                 end if
              end do
           end if

         end associate

      end do loop_cells

      end associate

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
    integer  :: icell, iface, gface, ivar, jcell
    integer  :: iv, gv, nfv
    real(dp) :: nf(3), dunit(3), kvec(3), gradt(3)
    real(dp) :: t1(3), t2(3), fc(3), rv(3), Kf(3,3)
    real(dp) :: cosang, dnorm
    real(dp) :: av(8), bv(8), pv(8)             ! in-face coords (a,b) and phi per vertex
    real(dp) :: abar, bbar, pbar, da, db, dphi
    real(dp) :: Saa, Sab, Sbb, Sap, Sbp, det, ca, cb
    type(real(dp)) , allocatable :: phiv(:), phic_v(:)

    ! Make space for skew source terms
    ss = real(0,dp)

    allocate(phic_v(this % grid % num_cells))

    ! One field at a time - gather its cell values, interpolate to
    ! vertices, then add its non-orthogonal correction to its own dofs.
    variables: do ivar = 1, this % g % num_variables

    do jcell = 1, this % grid % num_cells
       phic_v(jcell) = phic(this % g % dof(jcell, ivar))
    end do
    call this % evaluate_vertex_flux(phiv, phic_v)

    loop_cells: do icell = 1, this % grid % num_cells

       ! Get the faces corresponding to this cell
       associate( &
            & faces => this % grid % cell_faces &
            & (1:this % grid % num_cell_faces(icell),icell) &
            & )

           loop_faces : do iface = 1, this % grid % num_cell_faces(icell)

              ! Global face number
              gface = faces(iface)

              ! Skew (non-orthogonal) correction applies to interior faces only
              domain: if (this % grid % num_face_cells(gface) .eq. 2) then

                 ! Unit face normal as seen from this cell
                 nf = this % grid % cell_face_normals(1:3,iface,icell)

                 ! Unit centroid-to-centroid vector
                 dnorm  = norm2(this % grid % lvec(1:3,gface))
                 dunit  = this % grid % lvec(1:3,gface)/dnorm
                 cosang = dot_product(dunit, nf)

                 ! Minimum-correction vector: the component of the face
                 ! normal orthogonal to the centroid vector. It is bounded
                 ! (|kvec| = |sin(angle)| <= 1), so the explicit correction
                 ! stays finite on highly skewed/non-orthogonal cells where
                 ! the old (lvec.t)/delta form blew up (delta -> 0).
                 kvec = nf - cosang*dunit

                 ! Least-squares tangential gradient over ALL face vertices.
                 ! Work in the in-face coordinates (a,b) spanned by the
                 ! orthonormal face-tangent pair (t1,t2), fitting
                 !   phi ~ phi0 + ca*a + cb*b
                 ! so grad_t = ca*t1 + cb*t2 lies in the face plane. For a
                 ! triangle (3 pts, 3 unknowns) this is exact; for a quad it
                 ! is a least-squares plane. More accurate than a single-edge
                 ! difference on stretched/warped faces.
                 t1  = this % grid % cell_face_tangents(1:3,1,iface,icell)
                 t2  = this % grid % cell_face_tangents(1:3,2,iface,icell)
                 fc  = this % grid % face_centers(1:3,gface)
                 nfv = this % grid % num_face_vertices(gface)

                 abar = 0.0_dp
                 bbar = 0.0_dp
                 pbar = 0.0_dp

                 do iv = 1, nfv
                    gv = this % grid % face_vertices(iv, gface)
                    rv = this % grid % vertices(1:3,gv) - fc
                    av(iv) = dot_product(rv, t1)
                    bv(iv) = dot_product(rv, t2)
                    pv(iv) = phiv(gv)
                    abar = abar + av(iv)
                    bbar = bbar + bv(iv)
                    pbar = pbar + pv(iv)
                 end do
                 abar = abar/real(nfv,dp)
                 bbar = bbar/real(nfv,dp)
                 pbar = pbar/real(nfv,dp)

                 Saa = 0.0_dp
                 Sab = 0.0_dp
                 Sbb = 0.0_dp
                 Sap = 0.0_dp
                 Sbp = 0.0_dp

                 do iv = 1, nfv
                    da   = av(iv) - abar
                    db   = bv(iv) - bbar
                    dphi = pv(iv) - pbar
                    Saa  = Saa + da*da
                    Sab  = Sab + da*db
                    Sbb  = Sbb + db*db
                    Sap  = Sap + da*dphi
                    Sbp  = Sbp + db*dphi
                 end do

                 ! Solve the 2x2 normal equations for the in-plane slopes
                 det = Saa*Sbb - Sab*Sab
                 if (abs(det) .gt. tiny(1.0_dp)) then
                    ca    = (Sbb*Sap - Sab*Sbp)/det
                    cb    = (Saa*Sbp - Sab*Sap)/det
                    gradt = ca*t1 + cb*t2
                 else
                    gradt = 0.0_dp     ! degenerate face: skip correction
                 end if

                 ! Deferred non-orthogonal correction flux, projected
                 ! through the diffusion tensor: Area*((K grad_t) . kvec)
                 Kf = this % model % diffusion_tensor(fc)
                 associate(p => this % g % dof(icell, ivar))
                 ss(p) = ss(p) &
                      & + this % grid % face_areas(gface)*dot_product(matmul(Kf, gradt), kvec)
                 end associate

            end if domain

           end do loop_faces

         end associate

      end do loop_cells

      end do variables

      ! Deallocate the vertex flux values
      deallocate(phiv)
      deallocate(phic_v)

    end block evaluate

  end subroutine get_skew_source

  subroutine get_source(this, b)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: b(:)

    add_boundary_terms: block

      integer  :: icell, iface, ivar
      real(dp) :: nf(3), Kf(3,3), keff

      associate(nv => this % g % num_variables)

      ! Loop cells
      ! plain do: the per-face temporaries (nf, Kf, keff) are shared, so
      ! this is not a safe do concurrent without f2018 locality specs
      loop_cells: do icell = 1, this % grid % num_cells

         ! Get the faces corresponding to this cell
         associate( &
              & faces => this % grid % cell_faces &
              & (1:this % grid % num_cell_faces(icell),icell) &
              & )

           ! Zero this cell's dofs
           do ivar = 1, nv
              b(this % g % dof(icell,ivar)) = 0.0d0
           end do

         loop_faces: do iface = 1, this % grid % num_cell_faces(icell)

            associate (&
                 & ftag   => this % grid % face_tags(faces(iface))  , &
                 & fdelta => this % grid % face_deltas(faces(iface)), &
                 & farea  => this % grid % face_areas(faces(iface))  &
                 & )

              ! Boundary face: add each variable's bc constant to the rhs
              ! (the phi_p part went to the diagonal in the jvp).
              boundary_faces: if (this % grid % num_face_cells(faces(iface)) .eq. 1) then

                 nf   = this % grid % cell_face_normals(1:3, iface, icell)
                 Kf   = this % model % diffusion_tensor(this % grid % face_centers(1:3, faces(iface)))
                 keff = dot_product(nf, matmul(Kf, nf))

                 do ivar = 1, nv
                    associate(p => this % g % dof(icell,ivar))
                    b(p) = b(p) + this % bcs(ftag,ivar) % rhs_coeff(farea, fdelta, keff)
                    end associate
                 end do

              end if boundary_faces

            end associate

         end do loop_faces

       end associate

      end do loop_cells

      end associate

    end block add_boundary_terms

    cell_source: block

      integer :: icell, ivar

      associate(nv => this % g % num_variables)

      loop_cells: do icell = 1, this % grid % num_cells

         associate(&
              & x => this % grid % cell_centers(:,icell), &
              & cell_volume => this % grid % cell_volumes(icell)&
              & )

           do ivar = 1, nv
              associate(p => this % g % dof(icell,ivar))
              b(p) = b(p) + this % model % source(x)*cell_volume
              end associate
           end do

         end associate

      end do loop_cells

      end associate

    end block cell_source

    ! Transient: rhs becomes (M/dt) phi_old - b_steady
    transient_rhs: block
      integer :: icell, ivar
      if (this % transient) then
         associate(nv => this % g % num_variables)
         do icell = 1, this % grid % num_cells
            do ivar = 1, nv
               associate(p => this % g % dof(icell, ivar))
               b(p) = this % grid % cell_volumes(icell)/this % dt*this % phi_old(p) - b(p)
               end associate
            end do
         end do
         end associate
      end if
    end block transient_rhs

  end subroutine get_source

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine write_solution(this, filename, phic)

    use class_paraview_writer, only : paraview_writer

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: filename
    real(dp)        , intent(in) :: phic(:)

    class(paraview_writer), allocatable :: pwriter
    real(dp)    , allocatable :: cellfields(:,:)   ! (cell, variable)
    type(string), allocatable :: labels(:)
    integer :: icell, ivar, nv
    character(len=16) :: vname

    ! Scatter the flat dof vector into per-variable cell fields and write
    ! them all (3d, any element type) through the paraview writer.
    nv = this % g % num_variables
    allocate(cellfields(this % grid % num_cells, nv))
    allocate(labels(nv))
    do ivar = 1, nv
       do icell = 1, this % grid % num_cells
          cellfields(icell, ivar) = phic(this % g % dof(icell, ivar))
       end do
       write(vname,'(a,i0)') "phi", ivar
       labels(ivar) = string(trim(vname))
    end do

    allocate(pwriter, source = paraview_writer(this % grid))
    call pwriter % write(trim(filename), cellfields, labels)

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
