#include "scalar.fpp"

!=====================================================================!
! Module that contains class that assembles residual and performs
! matrix vector products on all elements in mesh
!
! Author: Komahan Boopathy
!=====================================================================!

module interface_assembler

  use iso_fortran_env  , only : dp => REAL64
  use class_csr        , only : csr_matrix
  use module_solve_mode, only : FORWARD, REVERSE, &
       &                        WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE
  implicit none

  private
  public :: assembler
  public :: WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE   ! re-exported tags

  !===================================================================!
  ! Assembler for the physical system
  !===================================================================!
  
  type, abstract :: assembler

     type(character(len=:)), allocatable :: description

     type(integer) :: num_state_vars
     type(logical) :: approximate_jacobian
     type(integer) :: differential_order

     real(dp), allocatable :: S(:,:)

   contains  

     ! Deferred procedures
     procedure(add_residual_interface)               , deferred :: add_residual
     procedure(add_jacobian_vector_product_interface), deferred :: add_jacobian_vector_product
     procedure(add_initial_condition_interface)      , deferred :: add_initial_condition

     ! The essentials a solver may ask of the system: the imbalance at a
     ! given answer, the ONE product (jacobian on a residual-descended
     ! vector; mode = traversal direction, part = subgraph), and the
     ! inner product that gives the space its geometry (owned by the
     ! system because distribution is the system's business - a
     ! partitioned system reduces across images here).
     procedure :: get_residual
     procedure :: get_jacobian_residual_product
     procedure :: inner_product

     ! Machine-precision audits of the product (verify-before witnesses;
     ! the finite-difference twin is only the coarse third witness)
     procedure :: verify_transpose_consistency
     procedure :: verify_parts_honesty

     ! steady transpose action behind the one product's REVERSE direction
     procedure, private :: transpose_product

     ! The pieces the essentials are composed from. The base defaults are
     ! trivial (zero); a spatial assembler overrides them. Solvers call
     ! the essentials above, never these directly.
     procedure :: get_source
     procedure :: get_skew_source
     procedure :: get_jacobian_vector_product
     procedure :: get_operator_csr   ! assembled entries - ONLY for building
                                     ! algebraic preconditioners, never to iterate
     
     ! Assembler knows the size of state array
     procedure :: create_vector
     procedure :: create_state

     ! Provided procedures
     procedure :: get_num_state_vars, set_num_state_vars
     procedure :: get_description   , set_description
  
     ! Defined procedures
     procedure :: get_differential_order
     procedure :: set_differential_order

     ! Adjoint support (provided defaults; physics overrides as needed)
     procedure :: add_jacobian_vector_product_transpose
     procedure :: get_num_design_vars
     procedure :: set_design_vars
     procedure :: get_design_vars
     procedure :: add_design_residual_transpose_product

     ! Post-processing (provided default: no-op; an assembler that owns a
     ! mesh overrides these to export the named fields)
     procedure :: write_solution_fields
     procedure :: write_gmsh_series

  end type assembler

  !===================================================================!
  ! Interfaces for deferred procedures
  !===================================================================!

  abstract interface

     !================================================================!
     ! Interface for residual assembly R(U,xi)
     !================================================================!

     impure subroutine add_residual_interface(this, residual, filter)

       import :: assembler

       class(assembler), intent(in)           :: this
       type(scalar)    , intent(inout)        :: residual(:)
       type(integer)   , intent(in), optional :: filter

     end subroutine add_residual_interface

     !================================================================!
     ! Routine to return the product of jacobian matrix with a compati-
     ! ble vector pdt <---- [scalar(i)*dR(U,X)/dU(i)]*vec
     !================================================================!
     
     impure subroutine add_jacobian_vector_product_interface(this, pdt, vec, scalars, filter)

       import :: assembler

       class(assembler) , intent(in)    :: this
       type(scalar)     , intent(inout) :: pdt(:)
       type(scalar)     , intent(in)    :: vec(:)
       type(scalar)     , intent(in)    :: scalars(:)
       type(integer)    , intent(in), optional :: filter

     end subroutine add_jacobian_vector_product_interface

     !================================================================!
     ! Supplying the initial condition to march in time
     !================================================================!

     impure subroutine add_initial_condition_interface(this, U)

       import :: assembler

       class(assembler), intent(in)    :: this
       type(scalar)    , intent(inout) :: U(:,:)

     end subroutine add_initial_condition_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the number of state variables in the physical system
  !===================================================================!
  
  pure type(integer) function get_num_state_vars(this)

    class(assembler), intent(in) :: this

    get_num_state_vars = this % num_state_vars

  end function get_num_state_vars

  !===================================================================!
  ! Sets the number of state variables in the physical system
  !===================================================================!
  
  pure subroutine set_num_state_vars(this, num_state_vars)

    class(assembler), intent(inout) :: this
    type(integer)  , intent(in)   :: num_state_vars

    this % num_state_vars  = num_state_vars

  end subroutine set_num_state_vars
  
  !===================================================================!
  ! Returns the description set for the physical system
  !===================================================================!
  
  pure type(character) function get_description(this)

    class(assembler), intent(in) :: this

    get_description = this % description

  end function get_description

  !===================================================================!
  ! Sets the description for physical system
  !===================================================================!

  pure subroutine set_description(this, description)

    class(assembler), intent(inout) :: this
    type(character(len=*)), intent(in) :: description
    
    allocate(this % description, source = trim(description))

  end subroutine set_description

  !===================================================================!
  ! Returns the highest order of time derivative in the assembler
  !===================================================================!
  
  pure type(integer) function get_differential_order(this)

    class(assembler), intent(in) :: this

    get_differential_order = this % differential_order

  end function get_differential_order

  !===================================================================!
  ! Sets the highest order of time derivative in the assembler
  !===================================================================!

  pure subroutine set_differential_order(this, order)

    class(assembler), intent(inout) :: this
    type(integer) , intent(in)    :: order

    this % differential_order = order

  end subroutine set_differential_order

  !===================================================================!
  ! Create a state vector and sets values if a scalar is supplied
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
  ! Create a state vector and sets values if a scalar is supplied
  !===================================================================!

  impure subroutine create_state(this, S, val)
    
    class(assembler), intent(in)               :: this
    real(dp)        , intent(out), allocatable :: S(:,:)
    real(dp)        , intent(in) , optional    :: val
    
    if (allocated(S)) error stop "vector already allocated"
    allocate( &
         & S( &
         & this % num_state_vars, &
         & this % get_differential_order() + 1 &
         & ))
    if (present(val))  S = val

  end subroutine create_state

  !===================================================================!
  ! Transpose jacobian-vector product  pdt += [scalar(i) dR/dU(i)]^T vec.
  ! Default: assume a symmetric jacobian (A^T = A), so the transpose
  ! action equals the forward one. Physics with a non-symmetric operator
  ! (e.g. advection) overrides this with a true transpose.
  !===================================================================!

  impure subroutine add_jacobian_vector_product_transpose(this, pdt, vec, scalars, filter)

    class(assembler), intent(in)           :: this
    type(scalar)    , intent(inout)        :: pdt(:)
    type(scalar)    , intent(in)           :: vec(:)
    type(scalar)    , intent(in)           :: scalars(:)
    type(integer)   , intent(in), optional :: filter

    call this % add_jacobian_vector_product(pdt, vec, scalars, filter)

  end subroutine add_jacobian_vector_product_transpose

  !===================================================================!
  ! Number of design variables x the residual depends on. Default: none
  ! (no design dependence), so the adjoint total derivative is just the
  ! function's explicit df/dx. Physics with design variables overrides.
  !===================================================================!

  pure type(integer) function get_num_design_vars(this)

    class(assembler), intent(in) :: this

    get_num_design_vars = 0

  end function get_num_design_vars

  !===================================================================!
  ! Set / get the design variables. Default: no design dependence, so
  ! these are no-ops. Physics carrying design variables overrides them.
  !===================================================================!

  pure subroutine set_design_vars(this, x)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(in)    :: x(:)

  end subroutine set_design_vars

  pure subroutine get_design_vars(this, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: x(:)

    x = 0.0_dp

  end subroutine get_design_vars

  !===================================================================!
  ! Accumulate the adjoint design contribution  dfdx += psi^T dR/dx,
  ! the product of the adjoint variables with the residual's design
  ! jacobian. Default: no design dependence, so this adds nothing.
  ! Physics with design variables overrides it (analytic, or finite
  ! differenced as a temporary stand-in).
  !===================================================================!

  impure subroutine add_design_residual_transpose_product(this, dfdx, psi)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(inout) :: dfdx(:)
    type(scalar)    , intent(in)    :: psi(:)

  end subroutine add_design_residual_transpose_product


  !===================================================================!
  ! Export named flat-dof fields (state, adjoint state, ...) for post-
  ! processing. fields is (num_state_vars, nfield), labels names each.
  ! Default: no-op (an abstract / mesh-less system has nothing to write);
  ! a mesh-backed assembler overrides this to write a real file.
  !===================================================================!

  impure subroutine write_solution_fields(this, filename, fields, labels)

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: filename
    real(dp)        , intent(in) :: fields(:,:)
    character(len=*), intent(in) :: labels(:)

  end subroutine write_solution_fields

  !===================================================================!
  ! Export named flat-dof fields over a time series as a gmsh post file.
  ! fields is (num_state_vars, nfield, nstep); names labels each field
  ! (a gmsh view), times gives the time of each step. meshfile is the
  ! source mesh copied verbatim and keyed by its element tags. Default:
  ! no-op; a mesh-backed assembler overrides it. (Steady = nstep 1.)
  !===================================================================!

  impure subroutine write_gmsh_series(this, meshfile, filename, fields, names, times)

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: meshfile, filename
    real(dp)        , intent(in) :: fields(:,:,:)
    character(len=*), intent(in) :: names(:)
    real(dp)        , intent(in) :: times(:)

  end subroutine write_gmsh_series

  !===================================================================!
  ! The imbalance at answer x:  r = R(x). For the linear system this is
  ! the constant source plus the solution-dependent correction, minus
  ! the operator action - composed here so a solver never names the
  ! pieces. Forward only: the adjoint's residual needs the functional
  ! and stays on the linearized path until that unification lands.
  !===================================================================!

  pure subroutine get_residual(this, r, x)

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
  ! The ONE product: the jacobian's action on a residual-descended
  ! vector. mode picks the traversal direction (FORWARD = J v,
  ! REVERSE = J^T v); part picks the subgraph (WHOLE, DIAGONAL,
  ! LOWER_TRIANGLE, UPPER_TRIANGLE). Defaults: FORWARD, WHOLE.
  !
  ! The REVERSE direction currently inherits the symmetric default
  ! (J^T = J) of add_jacobian_vector_product_transpose; a system with a
  ! non-symmetric operator must override that transpose before REVERSE
  ! products may be trusted - verify_transpose_consistency is the audit.
  !===================================================================!

  pure subroutine get_jacobian_residual_product(this, w, v, mode, part)

    class(assembler), intent(in)           :: this
    real(dp)        , intent(out)          :: w(:)
    real(dp)        , intent(in)           :: v(:)
    integer         , intent(in), optional :: mode
    integer         , intent(in), optional :: part

    integer :: dir, sub

    dir = FORWARD
    if (present(mode)) dir = mode
    sub = WHOLE
    if (present(part)) sub = part

    if (dir .eq. REVERSE) then
       ! transpose action at the steady linearization; symmetric systems
       ! inherit J^T = J, non-symmetric ones override the transpose
       call this % transpose_product(w, v, sub)
    else
       if (sub .eq. WHOLE) then
          call this % get_jacobian_vector_product(w, v)
       else
          call this % get_jacobian_vector_product(w, v, filter = sub)
       end if
    end if

  end subroutine get_jacobian_residual_product

  !===================================================================!
  ! Steady transpose action w = J^T v (helper for the one product).
  ! Default: the symmetric assumption J^T = J, matching the default of
  ! add_jacobian_vector_product_transpose.
  !===================================================================!

  pure subroutine transpose_product(this, w, v, sub)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: w(:)
    real(dp)        , intent(in)  :: v(:)
    integer         , intent(in)  :: sub

    if (sub .eq. WHOLE) then
       call this % get_jacobian_vector_product(w, v)
    else
       call this % get_jacobian_vector_product(w, v, filter = sub)
    end if

  end subroutine transpose_product

  !===================================================================!
  ! The inner product of two vectors of this system's space. The system
  ! owns it because distribution is the system's business: this serial
  ! default is the plain dot; a partitioned system sums its own rows and
  ! reduces across images. The volume weighting attaches here later.
  !===================================================================!

  pure real(dp) function inner_product(this, a, b)

    class(assembler), intent(in) :: this
    real(dp)        , intent(in) :: a(:)
    real(dp)        , intent(in) :: b(:)

    inner_product = dot_product(a, b)

  end function inner_product

  !===================================================================!
  ! Audit: <w, J v> = <J^T w, v> for deterministic pseudo-random v, w,
  ! per part. Two products and two inner products - an analytic witness
  ! at machine precision, no truncation error. Returns the largest
  ! relative defect over the parts.
  !===================================================================!

  impure real(dp) function verify_transpose_consistency(this) result(defect)

    class(assembler), intent(in) :: this

    real(dp), allocatable :: v(:), w(:), jv(:), jtw(:)
    real(dp) :: lhs, rhs, scale
    integer  :: parts(4), i, n

    n = this % num_state_vars
    allocate(v(n), w(n), jv(n), jtw(n))
    call fill_deterministic(v, 17)
    call fill_deterministic(w, 31)

    parts  = [WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE]
    defect = 0.0_dp

    do i = 1, size(parts)
       call this % get_jacobian_residual_product(jv,  v, mode = FORWARD, part = parts(i))
       call this % get_jacobian_residual_product(jtw, w, mode = REVERSE, part = parts(i))
       ! transpose of the lower triangle is the upper triangle: pair them
       if (parts(i) .eq. LOWER_TRIANGLE) then
          call this % get_jacobian_residual_product(jtw, w, mode = REVERSE, part = UPPER_TRIANGLE)
       else if (parts(i) .eq. UPPER_TRIANGLE) then
          call this % get_jacobian_residual_product(jtw, w, mode = REVERSE, part = LOWER_TRIANGLE)
       end if
       lhs   = this % inner_product(w, jv)
       rhs   = this % inner_product(jtw, v)
       scale = max(abs(lhs), abs(rhs), 1.0_dp)
       defect = max(defect, abs(lhs - rhs)/scale)
    end do

  end function verify_transpose_consistency

  !===================================================================!
  ! Audit: (diagonal + lower + upper) v = whole v for a deterministic
  ! pseudo-random v. Three part-products against one whole-product;
  ! catches every part-implementation bug at machine precision. Returns
  ! the relative defect.
  !===================================================================!

  impure real(dp) function verify_parts_honesty(this) result(defect)

    class(assembler), intent(in) :: this

    real(dp), allocatable :: v(:), wd(:), wl(:), wu(:), wf(:)
    integer :: n

    n = this % num_state_vars
    allocate(v(n), wd(n), wl(n), wu(n), wf(n))
    call fill_deterministic(v, 7)

    call this % get_jacobian_residual_product(wd, v, part = DIAGONAL)
    call this % get_jacobian_residual_product(wl, v, part = LOWER_TRIANGLE)
    call this % get_jacobian_residual_product(wu, v, part = UPPER_TRIANGLE)
    call this % get_jacobian_residual_product(wf, v, part = WHOLE)

    defect = norm2((wd + wl + wu) - wf)/max(norm2(wf), 1.0_dp)

  end function verify_parts_honesty

  !===================================================================!
  ! Deterministic pseudo-random fill (linear congruential), so the
  ! audits are reproducible run to run.
  !===================================================================!

  pure subroutine fill_deterministic(v, seed)

    real(dp), intent(out) :: v(:)
    integer , intent(in)  :: seed

    integer :: i, s

    s = seed
    do i = 1, size(v)
       s    = mod(s*1103515245 + 12345, 2147483647)
       v(i) = real(mod(s, 10000), dp)/10000.0_dp - 0.5_dp
    end do

  end subroutine fill_deterministic

  !===================================================================!
  ! Default right-hand side: none (a spatial assembler with boundary
  ! conditions / a volumetric source overrides this).
  !===================================================================!

  pure subroutine get_source(this, b, boundary_only)

    class(assembler), intent(in)           :: this
    real(dp)        , intent(out)          :: b(:)
    logical         , intent(in), optional :: boundary_only

    b = 0.0_dp

  end subroutine get_source

  !===================================================================!
  ! Default non-orthogonal skew correction: none (orthogonal default).
  !===================================================================!

  pure subroutine get_skew_source(this, ss, phic)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(in)  :: phic(:)
    real(dp)        , intent(out) :: ss(:)

    ss = 0.0_dp

  end subroutine get_skew_source

  !===================================================================!
  ! Default operator action Aq = A q via the deferred jacobian-vector
  ! product at the steady linearization (dR/du only).
  !===================================================================!

  pure subroutine get_jacobian_vector_product(this, Aq, q, filter)

    class(assembler) , intent(in)           :: this
    real(dp)         , intent(out)          :: Aq(:)
    real(dp)         , intent(in)           :: q(:)
    integer          , intent(in), optional :: filter

    Aq = 0.0_dp

  end subroutine get_jacobian_vector_product

  !===================================================================!
  ! Default: this assembler does not assemble a sparse operator. A
  ! spatial assembler (class_assembler) overrides this with the real csr.
  !===================================================================!

  impure subroutine get_operator_csr(this, A)

    class(assembler), intent(in)  :: this
    type(csr_matrix), intent(out) :: A

    error stop "get_operator_csr: not implemented for this assembler"

  end subroutine get_operator_csr

end module interface_assembler
