#include "scalar.fpp"

module class_assembler

  ! import dependencies
  use iso_fortran_env         , only : dp => REAL64
  use interface_assembler     , only : base_assembler => assembler, &
       & DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE
  use class_mesh              , only : mesh
  use class_string            , only : string
  use class_boundary_condition, only : boundary_condition, dirichlet, neumann, robin
  use interface_flux          , only : flux
  use interface_physics       , only : source, point_state
  use class_diffusion_flux    , only : diffusion_flux, constant_source
  use class_fvm_field         , only : fvm_field
  use class_csr               , only : csr_matrix
  use module_verbosity        , only : verbosity

  implicit none

  private
  public :: assembler
  public :: CONVECTION_CENTRAL, CONVECTION_UPWIND
  public :: DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE   ! re-exported filter codes

  ! Convection scheme for the advective face flux
  integer, parameter :: CONVECTION_CENTRAL = 1   ! 2nd-order central differencing
  integer, parameter :: CONVECTION_UPWIND  = 2   ! 1st-order upwind (stable, high Peclet)

  !===================================================================!
  ! Class responsible for matrix, right hand side assembly and boundary
  ! conditions
  !===================================================================!

  type, extends(base_assembler) :: assembler

     ! Mesh object
     ! type(mesh), pointer :: grid
     class(mesh)   , allocatable :: grid
     !class(physics), allocatable :: system(:) ! poisson on \Omega, dirichlet on dOmega1 , dirchlet dOmega3 , dirichlet, Neumann dOmega4

     ! num_state_vars, differential_order and the state S(nvars,order+1)
     ! are inherited from the abstract base_assembler.

     ! Number of variables each cell
     integer :: num_variables

     ! Flux vector
     real(dp), allocatable :: phi(:)

     ! Boundary conditions - one per (face tag, variable), addressed by name
     type(boundary_condition), allocatable :: bcs(:,:)

     ! The law-agnostic seam: the flux operator, the source operator, and
     ! the cell-centred reconstruction field. A concrete law (diffusion,
     ! advection, ...) is set by providing its flux and source.
     class(flux)        , allocatable :: fx     ! flux operator F(q, grad q)
     class(source)      , allocatable :: src    ! source operator S(q, grad q)
     type(fvm_field)    , allocatable :: fld    ! cell-centred reconstruction

     ! Convective face-flux scheme (central by default; upwind for stability)
     integer :: convection_scheme = CONVECTION_CENTRAL

   contains

     procedure :: create_vector

     ! Boundary conditions - addressed by physical group name
     procedure :: set_dirichlet
     procedure :: set_neumann
     procedure :: set_robin
     procedure :: set_dirichlet_tag
     procedure, private :: set_bc

     ! The pde
     procedure :: set_equation
     procedure :: set_convection_scheme

     ! convection-scheme-dependent advective face closures (central/upwind)
     procedure, private :: adv_weights
     procedure, private :: adv_boundary

     ! Evaluation routines
     procedure :: evaluate_vertex_flux
     procedure :: evaluate_face_flux

     ! Assembly Routines
     procedure :: get_residual
     procedure :: get_source
     procedure :: get_skew_source
     procedure :: get_jacobian
     procedure :: get_transpose_jacobian
     procedure :: get_jacobian_vector_product
     ! law-agnostic flux/source assembly (the operator + rhs)
     procedure :: get_jvp_via_flux
     procedure :: get_source_via_flux
     procedure :: get_operator_csr       ! assemble the operator as sparse CSR (for AMG)
     procedure :: write_solution
     procedure :: write_solution_fields
     procedure :: write_gmsh_series

     ! Semi-discrete contract driven by the time integrator: the residual
     ! R = M*udot + A*u - b, its jacobian-vector product and the initial
     ! condition (deferred in base_assembler).
     procedure :: add_residual
     procedure :: add_jacobian_vector_product
     procedure :: add_initial_condition

     ! Design-variable sensitivity support: the design variables live on
     ! the equation (kappa); the residual design partial dR/dx is finite
     ! differenced for now (a stand-in for analytic / complex-step).
     procedure :: get_num_design_vars
     procedure :: set_design_vars
     procedure :: get_design_vars
     procedure :: add_design_residual_transpose_product

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

  impure type(assembler) function construct(grid) result (this)

    type(mesh), intent(in) :: grid

    if (verbosity .ge. 1) print *, "constructing assembler"

    ! Set mesh
    allocate(this % grid, source  = grid)
    if (verbosity .ge. 2) call this % grid % to_string()
    ! One variable per cell by default ("T"); set_equation bumps this up
    this % num_variables = 1

    ! The grid IS the dof/connectivity graph now; size the system from
    ! it. For a single field num_state_vars = num_cells, as before.
    this % grid % num_variables = this % num_variables
    this % num_state_vars = this % grid % num_dofs()

    ! Semi-discrete in time: R(u,udot) = M*udot + A*u - b is first order,
    ! and the jacobian-vector product below is exact (not approximated).
    ! Allocate the integrator state S(num_state_vars, order+1).
    call this % set_differential_order(1)
    this % approximate_jacobian = .false.
    call this % create_state(this % S, 0.0_dp)

    ! Allocate the flux vector
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

    ! Default every (boundary tag, variable) to homogeneous dirichlet
    ! (phi=0). The driver overrides by name. Interior faces never look
    ! this table up.
    allocate(this % bcs(maxval(this % grid % tag_numbers), this % num_variables))

    ! Default physics: isotropic unit diffusion, no source - recovers the
    ! old laplace operator. The driver overrides with set_equation.
    allocate(this % fx , source = diffusion_flux(1.0_dp, this % num_variables))
    allocate(this % src, source = constant_source(0.0_dp, this % num_variables))
    allocate(this % fld, source = fvm_field(this % num_variables, this % num_state_vars))

  end function construct

  !===================================================================!
  ! Set a dirichlet condition on a named physical group. The optional
  ! var targets a single variable; omitted, it applies to all variables.
  !===================================================================!

  impure subroutine set_dirichlet(this, name, value, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: value
    integer, optional, intent(in)    :: var

    call this % set_bc(this % grid % find_tag_by_name(name), name, dirichlet(value), var)

  end subroutine set_dirichlet

  !===================================================================!
  ! Set a neumann condition on a named physical group
  !===================================================================!

  impure subroutine set_neumann(this, name, flux, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: flux
    integer, optional, intent(in)    :: var

    call this % set_bc(this % grid % find_tag_by_name(name), name, neumann(flux), var)

  end subroutine set_neumann

  !===================================================================!
  ! Set a robin condition  a*phi + b*dphi/dn = c  on a named group
  !===================================================================!

  impure subroutine set_robin(this, name, a, b, c, var)

    class(assembler) , intent(inout) :: this
    character(len=*) , intent(in)    :: name
    real(dp)         , intent(in)    :: a, b, c
    integer, optional, intent(in)    :: var

    call this % set_bc(this % grid % find_tag_by_name(name), name, robin(a,b,c), var)

  end subroutine set_robin

  !===================================================================!
  ! Convenience for replicating tag-keyed setups without names
  !===================================================================!

  impure subroutine set_dirichlet_tag(this, tag, value, var)

    class(assembler) , intent(inout) :: this
    integer          , intent(in)    :: tag
    real(dp)         , intent(in)    :: value
    integer, optional, intent(in)    :: var

    call this % set_bc(tag, "", dirichlet(value), var)

  end subroutine set_dirichlet_tag

  !===================================================================!
  ! Set the pde being solved (diffusion tensor, source, #variables).
  ! Call this BEFORE setting boundary conditions - changing the number
  ! of variables rebuilds the (tag,variable) bc table.
  !===================================================================!

  impure subroutine set_equation(this, flux_op, source_op)

    class(assembler), intent(inout) :: this
    class(flux)     , intent(in)    :: flux_op    ! F(q, grad q)
    class(source)   , intent(in)    :: source_op  ! S(q, grad q)

    ! Guard the ordering footgun: set_equation rebuilds the (tag,variable)
    ! bc table, so any boundary conditions set earlier would be lost.
    if (allocated(this % bcs)) then
       if (any(this % bcs % tag .ne. -1)) then
          write(*,*) "warning: set_equation called after boundary conditions; ", &
               & "they have been reset - call set_equation first"
       end if
    end if

    ! Resize the system for this many variables per cell
    this % num_variables     = flux_op % num_components
    this % grid % num_variables = this % num_variables
    this % num_state_vars    = this % grid % num_dofs()

    if (allocated(this % phi)) deallocate(this % phi)
    allocate(this % phi(this % num_state_vars))
    this % phi = 0

    if (allocated(this % bcs)) deallocate(this % bcs)
    allocate(this % bcs(maxval(this % grid % tag_numbers), this % num_variables))

    ! store the flux/source operators and the reconstruction field
    if (allocated(this % fx))  deallocate(this % fx)
    if (allocated(this % src)) deallocate(this % src)
    if (allocated(this % fld)) deallocate(this % fld)
    allocate(this % fx , source = flux_op)
    allocate(this % src, source = source_op)
    allocate(this % fld, source = fvm_field(this % num_variables, this % num_state_vars))

  end subroutine set_equation

  !===================================================================!
  ! Stamp a resolved bc into the table (module-private helper). With no
  ! variable index the bc applies to every variable on the boundary.
  !===================================================================!

  impure subroutine set_bc(this, tag, name, bc, var)

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

  pure subroutine get_jacobian(this, A, filter)

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

  pure subroutine get_transpose_jacobian(this, A, filter)

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

  pure subroutine get_jacobian_vector_product(this, Aq, q, filter)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: q(:)
    real(dp)         , intent(out)   :: Aq(:)
    integer, optional, intent(in)    :: filter

    ! The diffusion operator is assembled law-agnostically through the
    ! flux seam: F = -K grad q, reconstructed at faces and dotted with
    ! the face normals. (Was the hard-coded laplace_normal block.)
    call this % get_jvp_via_flux(Aq, q, filter)

  end subroutine get_jacobian_vector_product

  !===================================================================!
  ! Jacobian-vector product via the flux seam (parallel to the legacy
  ! get_jacobian_vector_product, used to pin the new path against the
  ! old before the switch). The normal diffusivity keff = n^T K n comes
  ! from the flux's dF/d(grad q); the diagonal (own) and neighbour face-
  ! jacobian contributions are formed as separate sub-expressions so the
  ! L/U/D split stays bit-exact.
  !===================================================================!

  pure subroutine get_jvp_via_flux(this, Aq, q, filter)

    class(assembler) , intent(in)    :: this
    real(dp)         , intent(in)    :: q(:)
    real(dp)         , intent(out)   :: Aq(:)
    integer, optional, intent(in)    :: filter

    integer           :: icell, iface, ivar, p, n, ncell, fcells(2), gface
    real(dp)          :: nf(3), keff, vn, wp, wn, alhs, arhs, diag, neigh
    type(point_state) :: st

    associate(nv => this % grid % num_variables)

    allocate(st % q(nv), st % gradq(3, nv))
    st % nv = nv
    st % q  = 0.0d0

    loop_cells: do icell = 1, this % grid % num_cells

       associate(faces => this % grid % cell_faces(1:this % grid % num_cell_faces(icell), icell))

         do ivar = 1, nv
            Aq(this % grid % dof(icell, ivar)) = 0.0d0
         end do

         loop_faces: do iface = 1, this % grid % num_cell_faces(icell)

            gface = faces(iface)

            associate( &
                 & fdelta => this % grid % face_deltas(gface), &
                 & farea  => this % grid % face_areas(gface),  &
                 & ftag   => this % grid % face_tags(gface),   &
                 & nfc    => this % grid % num_face_cells(gface))

              nf        = this % grid % cell_face_normals(1:3, iface, icell)
              st % x    = this % grid % face_centers(1:3, gface)
              st % gradq = 0.0d0

              domain: if (nfc .eq. 2) then

                 fcells(1:nfc) = this % grid % face_cells(1:nfc, gface)
                 if (fcells(1) .eq. icell) then
                    ncell = fcells(2)
                 else
                    ncell = fcells(1)
                 end if

                 do ivar = 1, nv
                    p = this % grid % dof(icell, ivar)
                    n = this % grid % dof(ncell, ivar)
                    keff  = this % fx % normal_diffusivity(st, nf, ivar)
                    vn    = this % fx % normal_speed(st, nf, ivar)
                    call this % adv_weights(vn, wp, wn)
                    ! d(F.n)/dq_p and /dq_n: diffusion (keff/fdelta) + advection
                    diag  = -farea*( keff/fdelta + wp)*q(p)
                    neigh = -farea*(-keff/fdelta + wn)*q(n)
                    if (present(filter)) then
                       if (filter .eq. UPPER_TRIANGLE) then
                          if (ncell .gt. icell) Aq(p) = Aq(p) + neigh
                       else if (filter .eq. LOWER_TRIANGLE) then
                          if (ncell .lt. icell) Aq(p) = Aq(p) + neigh
                       else if (filter .eq. DIAGONAL) then
                          Aq(p) = Aq(p) + diag
                       end if
                    else
                       Aq(p) = Aq(p) + diag + neigh
                    end if
                 end do

              else

                 do ivar = 1, nv
                    p    = this % grid % dof(icell, ivar)
                    keff = this % fx % normal_diffusivity(st, nf, ivar)
                    vn   = this % fx % normal_speed(st, nf, ivar)
                    call this % adv_boundary(this % bcs(ftag,ivar), farea, fdelta, vn, &
                         &            alhs, arhs)
                    associate(bc => this % bcs(ftag,ivar))
                      if (present(filter)) then
                         if (filter .eq. DIAGONAL) &
                              & Aq(p) = Aq(p) + ( bc % lhs_coeff(farea, fdelta, keff) + alhs )*q(p)
                      else
                         Aq(p) = Aq(p) + ( bc % lhs_coeff(farea, fdelta, keff) + alhs )*q(p)
                      end if
                    end associate
                 end do

              end if domain

            end associate

         end do loop_faces

       end associate

    end do loop_cells

    end associate

  end subroutine get_jvp_via_flux

  !===================================================================!
  ! Assemble the steady operator as a sparse CSR matrix - the same matrix
  ! the matrix-free get_jvp_via_flux applies, but with explicit entries so
  ! algebraic multigrid can build a hierarchy on it. Mirrors that routine's
  ! face loop: interior face contributes A(p,p) += -farea*keff/fdelta and
  ! A(p,n) += +farea*keff/fdelta; a boundary face contributes the bc
  ! lhs_coeff to the diagonal. Single-field-or-decoupled (the ivar loop is
  ! the multi-field hook); steady only.
  !===================================================================!

  impure subroutine get_operator_csr(this, A)

    class(assembler), intent(in)  :: this
    type(csr_matrix), intent(out) :: A

    integer , allocatable :: nnz_row(:), row_ptr(:), col_idx(:), cursor(:)
    integer               :: ndof, icell, iface, ivar, p, n, ncell, fcells(2), gface, ftag, nfc
    real(dp)              :: nf(3), keff, vn, wp, wn, alhs, arhs, farea, fdelta
    type(point_state)     :: st

    associate(nv => this % grid % num_variables)

    ndof = this % grid % num_dofs()

    !---------------------------------------------------------------!
    ! symbolic pass: row p has its diagonal + one column per interior
    ! face neighbour (same variable)
    !---------------------------------------------------------------!
    allocate(nnz_row(ndof)); nnz_row = 1
    do icell = 1, this % grid % num_cells
       do iface = 1, this % grid % num_cell_faces(icell)
          gface = this % grid % cell_faces(iface, icell)
          if (this % grid % num_face_cells(gface) .eq. 2) then
             do ivar = 1, nv
                p = this % grid % dof(icell, ivar)
                nnz_row(p) = nnz_row(p) + 1
             end do
          end if
       end do
    end do

    allocate(row_ptr(ndof + 1))
    row_ptr(1) = 1
    do p = 1, ndof
       row_ptr(p+1) = row_ptr(p) + nnz_row(p)
    end do
    allocate(col_idx(row_ptr(ndof+1) - 1))

    ! seed each row's diagonal, then append neighbour columns
    allocate(cursor(ndof))
    do p = 1, ndof
       col_idx(row_ptr(p)) = p
       cursor(p) = row_ptr(p) + 1
    end do
    do icell = 1, this % grid % num_cells
       do iface = 1, this % grid % num_cell_faces(icell)
          gface = this % grid % cell_faces(iface, icell)
          if (this % grid % num_face_cells(gface) .eq. 2) then
             fcells(1:2) = this % grid % face_cells(1:2, gface)
             if (fcells(1) .eq. icell) then
                ncell = fcells(2)
             else
                ncell = fcells(1)
             end if
             do ivar = 1, nv
                p = this % grid % dof(icell, ivar)
                col_idx(cursor(p)) = this % grid % dof(ncell, ivar)
                cursor(p) = cursor(p) + 1
             end do
          end if
       end do
    end do

    A = csr_matrix(ndof, ndof, row_ptr, col_idx)   ! pattern, values zeroed

    !---------------------------------------------------------------!
    ! numeric pass: identical entries to get_jvp_via_flux
    !---------------------------------------------------------------!
    allocate(st % q(nv), st % gradq(3, nv))
    st % nv = nv
    st % q  = 0.0d0

    do icell = 1, this % grid % num_cells
       do iface = 1, this % grid % num_cell_faces(icell)

          gface  = this % grid % cell_faces(iface, icell)
          nf     = this % grid % cell_face_normals(1:3, iface, icell)
          st % x = this % grid % face_centers(1:3, gface)
          st % gradq = 0.0d0
          farea  = this % grid % face_areas(gface)
          fdelta = this % grid % face_deltas(gface)
          ftag   = this % grid % face_tags(gface)
          nfc    = this % grid % num_face_cells(gface)

          if (nfc .eq. 2) then
             fcells(1:2) = this % grid % face_cells(1:2, gface)
             if (fcells(1) .eq. icell) then
                ncell = fcells(2)
             else
                ncell = fcells(1)
             end if
             do ivar = 1, nv
                p    = this % grid % dof(icell, ivar)
                n    = this % grid % dof(ncell, ivar)
                keff = this % fx % normal_diffusivity(st, nf, ivar)
                vn   = this % fx % normal_speed(st, nf, ivar)
                call this % adv_weights(vn, wp, wn)
                call A % add_entry(p, p, -farea*(keff/fdelta + wp))
                call A % add_entry(p, n,  farea*(keff/fdelta - wn))
             end do
          else
             do ivar = 1, nv
                p    = this % grid % dof(icell, ivar)
                keff = this % fx % normal_diffusivity(st, nf, ivar)
                vn   = this % fx % normal_speed(st, nf, ivar)
                call this % adv_boundary(this % bcs(ftag,ivar), farea, fdelta, vn, &
                     &            alhs, arhs)
                call A % add_entry(p, p, this % bcs(ftag,ivar) % lhs_coeff(farea, fdelta, keff) + alhs)
             end do
          end if

       end do
    end do

    end associate

  end subroutine get_operator_csr

  !===================================================================!
  ! Source via the flux seam (parallel to get_source): boundary-condition
  ! constants (keff from the flux) + volumetric source from the source
  ! operator + the transient rhs.
  !===================================================================!

  pure subroutine get_source_via_flux(this, b, boundary_only)

    class(assembler), intent(in)           :: this
    real(dp)        , intent(out)          :: b(:)
    logical         , intent(in), optional :: boundary_only

    logical           :: bnd_only
    integer           :: icell, iface, ivar, p, gface
    real(dp)          :: nf(3), keff, vn, alhs, arhs
    type(point_state) :: st
    type(scalar)      :: Sval(this % grid % num_variables)

    bnd_only = .false.
    if (present(boundary_only)) bnd_only = boundary_only

    associate(nv => this % grid % num_variables)

    allocate(st % q(nv), st % gradq(3, nv))
    st % nv = nv; st % q = 0.0d0; st % gradq = 0.0d0

    ! boundary-condition constants
    do icell = 1, this % grid % num_cells
       associate(faces => this % grid % cell_faces(1:this % grid % num_cell_faces(icell), icell))
         do ivar = 1, nv
            b(this % grid % dof(icell, ivar)) = 0.0d0
         end do
         do iface = 1, this % grid % num_cell_faces(icell)
            gface = faces(iface)
            associate( &
                 & fdelta => this % grid % face_deltas(gface), &
                 & farea  => this % grid % face_areas(gface),  &
                 & ftag   => this % grid % face_tags(gface))
              if (this % grid % num_face_cells(gface) .eq. 1) then
                 nf     = this % grid % cell_face_normals(1:3, iface, icell)
                 st % x = this % grid % face_centers(1:3, gface)
                 do ivar = 1, nv
                    p    = this % grid % dof(icell, ivar)
                    keff = this % fx % normal_diffusivity(st, nf, ivar)
                    vn   = this % fx % normal_speed(st, nf, ivar)
                    call this % adv_boundary(this % bcs(ftag,ivar), farea, fdelta, vn, &
                         &            alhs, arhs)
                    b(p) = b(p) + this % bcs(ftag,ivar) % rhs_coeff(farea, fdelta, keff) + arhs
                 end do
              end if
            end associate
         end do
       end associate
    end do

    ! volumetric source
    if (.not. bnd_only) then
       do icell = 1, this % grid % num_cells
          st % x = this % grid % cell_centers(:, icell)
          Sval   = this % src % value(st)
          do ivar = 1, nv
             p    = this % grid % dof(icell, ivar)
             b(p) = b(p) + real(Sval(ivar), dp)*this % grid % cell_volumes(icell)
          end do
       end do
    end if

    end associate

  end subroutine get_source_via_flux

  !===================================================================!
  ! Select the convective face-flux scheme (CONVECTION_CENTRAL / _UPWIND)
  !===================================================================!

  pure subroutine set_convection_scheme(this, scheme)
    class(assembler), intent(inout) :: this
    integer         , intent(in)    :: scheme
    this % convection_scheme = scheme
  end subroutine set_convection_scheme

  !===================================================================!
  ! Interior advective face weights: the face flux is
  !   Phi_adv = farea * (wp*q_p + wn*q_n).
  ! central -> wp = wn = vn/2 (2nd order, skew-symmetric); upwind -> the
  ! upstream cell carries it (wp = max(vn,0), wn = min(vn,0); 1st order,
  ! unconditionally stable - adds |vn|/2 numerical diffusion).
  !===================================================================!

  pure subroutine adv_weights(this, vn, wp, wn)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(in)  :: vn
    real(dp)        , intent(out) :: wp, wn

    if (this % convection_scheme .eq. CONVECTION_UPWIND) then
       wp = max(vn, 0.0_dp)
       wn = min(vn, 0.0_dp)
    else
       wp = 0.5_dp*vn
       wn = 0.5_dp*vn
    end if

  end subroutine adv_weights

  !===================================================================!
  ! Boundary advective contribution: diagonal coefficient (on q_p) and the
  ! rhs constant. central and upwind-inflow use the bc face value (the same
  ! robin elimination as the diffusive closure); upwind-outflow (vn >= 0)
  ! upwinds the interior value q_p, so the flux farea*vn*q_p is all diagonal.
  !===================================================================!

  pure subroutine adv_boundary(this, bc, farea, fdelta, vn, alhs, arhs)

    class(assembler)        , intent(in)  :: this
    type(boundary_condition), intent(in)  :: bc
    real(dp)                , intent(in)  :: farea, fdelta, vn
    real(dp)                , intent(out) :: alhs, arhs

    if (this % convection_scheme .eq. CONVECTION_UPWIND .and. vn .ge. 0.0_dp) then
       alhs = -farea*vn
       arhs = 0.0_dp
    else
       alhs = bc % adv_lhs_coeff(farea, fdelta, vn)
       arhs = bc % adv_rhs_coeff(farea, fdelta, vn)
    end if

  end subroutine adv_boundary

  !===================================================================!
  ! Compute vertex values by interpolating cell center values
  !===================================================================!

  pure subroutine evaluate_vertex_flux(this, phiv, phic)

    class(assembler) , intent(in)               :: this
    real(dp)         , intent(in)               :: phic(:)
    real(dp)         , intent(out), allocatable :: phiv(:)

    integer :: ivertex

    ! Interpolate the supplied cell centered solution to form the
    ! nodal solution
    allocate(phiv(this % grid % num_points)); phiv = 0
    do concurrent (ivertex = 1: this % grid % num_points)
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

  pure subroutine evaluate_face_flux(this, phif, phic)

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
  ! The residual at state x, in this discretization's vocabulary: the
  ! constant source plus the solution-dependent skew correction, minus
  ! the operator action - composed here, on the layer that owns the
  ! words, and seen by solvers only through the deferred query.
  !===================================================================!

  impure subroutine get_residual(this, r, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: r(:)
    real(dp)        , intent(in)  :: x(:)

    real(dp), allocatable :: b(:), s(:), ax(:)

    allocate(b, s, ax, mold = x)

    call this % get_source(b)
    call this % get_skew_source(s, x)
    call this % get_jacobian_vector_product(ax, x)

    r = (b + s) - ax

  end subroutine get_residual

  !===================================================================!
  ! Evaluate internal skew source based on the current cell states and
  ! return
  !===================================================================!

  !subroutine get_tangential_flux(this, ss, phic)
  pure subroutine get_skew_source(this, ss, phic)

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
    variables: do ivar = 1, this % grid % num_variables

    do jcell = 1, this % grid % num_cells
       phic_v(jcell) = phic(this % grid % dof(jcell, ivar))
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
                    rv = this % grid % coordinates(1:3,gv) - fc
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
                 ! through the diffusion tensor: Area*((K grad_t) . kvec).
                 ! (Diffusion-specific correction; reads the tensor from the
                 ! diffusion flux.)
                 Kf = 0.0_dp
                 select type (fxp => this % fx)
                 type is (diffusion_flux)
                    Kf = fxp % kmat
                 end select
                 associate(p => this % grid % dof(icell, ivar))
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

  pure subroutine get_source(this, b, boundary_only)

    class(assembler), intent(in)           :: this
    real(dp)        , intent(out)          :: b(:)
    logical         , intent(in), optional :: boundary_only

    ! The source (boundary-condition constants + volumetric source + the
    ! transient rhs) is assembled law-agnostically through the flux seam.
    ! boundary_only drops the volumetric source (used by the adjoint
    ! design partial). (Was the add_boundary_terms / cell_source blocks.)
    call this % get_source_via_flux(b, boundary_only)

  end subroutine get_source

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  impure subroutine write_solution(this, filename, phic)

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
    nv = this % grid % num_variables
    allocate(cellfields(this % grid % num_cells, nv))
    allocate(labels(nv))
    do ivar = 1, nv
       do icell = 1, this % grid % num_cells
          cellfields(icell, ivar) = phic(this % grid % dof(icell, ivar))
       end do
       write(vname,'(a,i0)') "phi", ivar
       labels(ivar) = string(trim(vname))
    end do

    allocate(pwriter, source = paraview_writer(this % grid))
    call pwriter % write(trim(filename), cellfields, labels)

  end subroutine write_solution

  !===================================================================!
  ! Write several named flat-dof fields (e.g. the state and the adjoint
  ! state) as cell data in one paraview file. fields is (ndof, nfield)
  ! and labels names each. For multi-variable systems each field expands
  ! to one column per variable (suffixed _vN).
  !===================================================================!

  impure subroutine write_solution_fields(this, filename, fields, labels)

    use class_paraview_writer, only : paraview_writer

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: filename
    real(dp)        , intent(in) :: fields(:,:)   ! (ndof, nfield)
    character(len=*), intent(in) :: labels(:)     ! (nfield)

    class(paraview_writer), allocatable :: pwriter
    real(dp)    , allocatable :: cellfields(:,:)  ! (cell, column)
    type(string), allocatable :: colnames(:)
    integer :: icell, ivar, nv, ifield, nfield, col
    character(len=64) :: cname

    nv     = this % grid % num_variables
    nfield = size(fields, 2)

    allocate(cellfields(this % grid % num_cells, nfield*nv))
    allocate(colnames(nfield*nv))

    col = 0
    do ifield = 1, nfield
       do ivar = 1, nv
          col = col + 1
          do icell = 1, this % grid % num_cells
             cellfields(icell, col) = fields(this % grid % dof(icell, ivar), ifield)
          end do
          if (nv .eq. 1) then
             colnames(col) = string(trim(labels(ifield)))
          else
             write(cname, '(a,a,i0)') trim(labels(ifield)), "_v", ivar
             colnames(col) = string(trim(cname))
          end if
       end do
    end do

    allocate(pwriter, source = paraview_writer(this % grid))
    call pwriter % write(trim(filename), cellfields, colnames)

  end subroutine write_solution_fields

  !===================================================================!
  ! Export named flat-dof fields over a time series to a gmsh post file
  ! (overrides base_assembler % write_gmsh_series). Scatters each field
  ! to its cells (first variable), keys by the gmsh element tags and
  ! writes through the gmsh writer. fields is (ndof, nfield, nstep).
  !===================================================================!

  impure subroutine write_gmsh_series(this, meshfile, filename, fields, names, times)

    use class_gmsh_writer, only : gmsh_writer

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: meshfile, filename
    real(dp)        , intent(in) :: fields(:,:,:)   ! (ndof, nfield, nstep)
    character(len=*), intent(in) :: names(:)        ! (nfield)
    real(dp)        , intent(in) :: times(:)        ! (nstep)

    type(gmsh_writer)     :: gw
    real(dp), allocatable :: cellvals(:,:,:)        ! (ncell, nfield, nstep)
    integer , allocatable :: cell_numbers(:)
    integer               :: icell, ifield, istep, ncell, nfield, nstep

    ncell  = this % grid % num_cells
    nfield = size(fields, 2)
    nstep  = size(fields, 3)

    allocate(cellvals(ncell, nfield, nstep), cell_numbers(ncell))

    ! scatter the flat dof fields to cell values (first variable) and
    ! collect the gmsh element tags the loader read
    do istep = 1, nstep
       do ifield = 1, nfield
          do icell = 1, ncell
             cellvals(icell, ifield, istep) = fields(this % grid % dof(icell, 1), ifield, istep)
          end do
       end do
    end do
    do icell = 1, ncell
       cell_numbers(icell) = this % grid % cell_numbers(icell)
    end do

    gw = gmsh_writer(meshfile)
    call gw % write_time_series(trim(filename), cell_numbers, names, times, cellvals)

    deallocate(cellvals, cell_numbers)

  end subroutine write_gmsh_series

  !===================================================================!
  ! Create a state vector and sets values if a value is supplied
  ! (overrides base_assembler % create_vector)
  !===================================================================!

  impure subroutine create_vector(this, x, val)

    class(assembler), intent(in)               :: this
    real(dp)        , intent(out), allocatable :: x(:)
    real(dp)        , intent(in) , optional    :: val

    if (allocated(x)) error stop "vector already allocated"
    allocate(x(this % num_state_vars))
    if (present(val))  x = val

  end subroutine create_vector

  !===================================================================!
  ! Semi-discrete residual  R = M*udot + A*u - b, assembled from the
  ! inherited state S (S(:,1) = u, S(:,2) = udot). The spatial operator
  ! A*u and the source b are the steady FVM pieces
  ! (get_jacobian_vector_product and get_source with the legacy transient
  ! flag off); M = cell volume is the mass term. The time integrator owns
  ! the marching, so the assembler's own transient flag stays off here.
  !
  ! Sign: R = M*udot - A*u + b, so the steady root (udot = 0) recovers the
  ! laplace solve A*u = b, and dR/du = -A makes the newton operator
  ! beta*M - alpha*A symmetric positive definite (CG-friendly).
  !===================================================================!

  pure subroutine add_residual(this, residual, filter)

    class(assembler), intent(in)           :: this
    type(scalar)    , intent(inout)        :: residual(:)
    integer         , intent(in), optional :: filter

    real(dp), allocatable :: Au(:), b(:)
    integer               :: icell, ivar, p

    allocate(Au(this % num_state_vars))
    allocate(b (this % num_state_vars))

    ! Steady spatial operator A*u and source b
    call this % get_jacobian_vector_product(Au, this % S(:,1), filter)
    call this % get_source(b)

    do icell = 1, this % grid % num_cells
       do ivar = 1, this % grid % num_variables
          p = this % grid % dof(icell, ivar)
          residual(p) = residual(p) &
               & + this % grid % cell_volumes(icell)*this % S(p,2) &
               & - Au(p) + b(p)
       end do
    end do

    deallocate(Au, b)

  end subroutine add_residual

  !===================================================================!
  ! Jacobian-vector product for the integrator
  !   pdt += [alpha dR/du + beta dR/dudot] vec = [beta*M - alpha*A] vec
  ! with scalars = [alpha, beta] (the linearization coefficients). The
  ! spatial operator action A*vec is the steady jacobian-vector product;
  ! M*vec is the cell-volume mass term.
  !===================================================================!

  pure subroutine add_jacobian_vector_product(this, pdt, vec, scalars, filter)

    class(assembler), intent(in)           :: this
    type(scalar)    , intent(inout)        :: pdt(:)
    type(scalar)    , intent(in)           :: vec(:)
    type(scalar)    , intent(in)           :: scalars(:)
    integer         , intent(in), optional :: filter

    real(dp), allocatable :: Av(:)
    integer               :: icell, ivar, p

    allocate(Av(this % num_state_vars))

    ! Steady spatial operator action A*vec
    call this % get_jacobian_vector_product(Av, vec, filter)

    do icell = 1, this % grid % num_cells
       do ivar = 1, this % grid % num_variables
          p = this % grid % dof(icell, ivar)
          pdt(p) = pdt(p) &
               & + scalars(2)*this % grid % cell_volumes(icell)*vec(p) &
               & - scalars(1)*Av(p)
       end do
    end do

    deallocate(Av)

  end subroutine add_jacobian_vector_product

  !===================================================================!
  ! Initial condition: u(0) from the stored field phi (a driver may prime
  ! it; defaults to zero); udot(0) = 0, since the first BDF step recovers
  ! the time derivative from the stencil.
  !===================================================================!

  pure subroutine add_initial_condition(this, U)

    class(assembler), intent(in)    :: this
    type(scalar)    , intent(inout) :: U(:,:)

    U(:,1) = this % phi
    U(:,2) = 0.0_dp

  end subroutine add_initial_condition

  !===================================================================!
  ! Design variables live on the flux operator (the isotropic
  ! conductivity kappa); the assembler forwards to it. The operator now
  ! reads kappa from fx, so perturbing a design variable must update fx
  ! (not the legacy model) for the forward solve to respond.
  !===================================================================!

  pure integer function get_num_design_vars(this)

    class(assembler), intent(in) :: this

    get_num_design_vars = this % fx % num_design_vars()

  end function get_num_design_vars

  pure subroutine set_design_vars(this, x)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(in)    :: x(:)

    call this % fx % set_design_vars(x)

  end subroutine set_design_vars

  pure subroutine get_design_vars(this, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: x(:)

    call this % fx % get_design_vars(x)

  end subroutine get_design_vars

  !===================================================================!
  ! Adjoint design contribution  dfdx(k) += psi^T dR/dx_k, assembled
  ! through the flux seam (state S held fixed):
  !
  !   dR/dx_k = ( + integral_faces (dF/dx_k . n)        )   interior flux
  !             ( + integral_volume (dS/dx_k)           )   source
  !             ( + boundary-condition design derivative )
  !
  ! The interior contribution is law-agnostic: the flux operator supplies
  ! its own design partial dF/dx_k (for diffusion dF/dkappa = -grad q),
  ! and the assembler integrates it - no assumption of linearity in the
  ! design variable. The boundary uses the (diffusion-shaped) Robin
  ! closure, whose coefficients scale with keff ~ kappa, so its design
  ! derivative is the boundary residual divided by kappa (generalising
  ! with the boundary flux is deferred). For constant diffusion this
  ! reproduces the exact (b_bc - A*u)/kappa.
  !===================================================================!

  pure subroutine add_design_residual_transpose_product(this, dfdx, psi)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(inout) :: dfdx(:)
    type(scalar)    , intent(in)    :: psi(:)

    real(dp), allocatable :: kappa(:), dRdk(:)
    type(point_state)     :: st
    integer               :: ndv, n, k, icell, iface, ivar, p, gface

    ndv = this % get_num_design_vars()
    if (ndv .eq. 0) return

    n = this % num_state_vars
    allocate(kappa(ndv), dRdk(n))
    call this % get_design_vars(kappa)

    associate(nv => this % grid % num_variables)

    design_vars: do k = 1, ndv

       dRdk = 0.0_dp

       do icell = 1, this % grid % num_cells
          associate(faces => this % grid % cell_faces(1:this % grid % num_cell_faces(icell), icell))
          do iface = 1, this % grid % num_cell_faces(icell)

             gface = faces(iface)
             associate( &
                  & farea  => this % grid % face_areas(gface),  &
                  & fdelta => this % grid % face_deltas(gface), &
                  & ftag   => this % grid % face_tags(gface))

               ! reconstruct (q, grad q) at the face from the (fixed) state
               call this % fld % face_state(this % grid, this % grid, &
                    & this % S(:,1), icell, iface, gface, st)

               interior: if (this % grid % num_face_cells(gface) .eq. 2) then

                  ! integral of (dF/dx_k . n) over the face - the flux's own
                  ! design partial, integrated (law-agnostic)
                  block
                    type(scalar) :: dFk(3, nv)
                    real(dp)     :: nf(3)
                    nf  = this % grid % cell_face_normals(1:3, iface, icell)
                    dFk = this % fx % dflux_ddesign(st, k)
                    do ivar = 1, nv
                       p = this % grid % dof(icell, ivar)
                       dRdk(p) = dRdk(p) + farea*real(dot_product(dFk(:,ivar), nf), dp)
                    end do
                  end block

               else

                  ! boundary: Robin closure scales with keff ~ kappa, so the
                  ! design derivative is the boundary residual over kappa
                  block
                    real(dp) :: nf(3), keff, lhs, rhs
                    nf = this % grid % cell_face_normals(1:3, iface, icell)
                    do ivar = 1, nv
                       p    = this % grid % dof(icell, ivar)
                       keff = this % fx % normal_diffusivity(st, nf, ivar)
                       lhs  = this % bcs(ftag,ivar) % lhs_coeff(farea, fdelta, keff)
                       rhs  = this % bcs(ftag,ivar) % rhs_coeff(farea, fdelta, keff)
                       dRdk(p) = dRdk(p) + (rhs - lhs*real(this % S(p,1), dp))/kappa(k)
                    end do
                  end block

               end if interior

             end associate

          end do
          end associate
       end do

       dfdx(k) = dfdx(k) + real(dot_product(psi, dRdk), dp)

    end do design_vars

    end associate

    deallocate(kappa, dRdk)

  end subroutine add_design_residual_transpose_product



end module class_assembler
