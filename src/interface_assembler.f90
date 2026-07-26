#include "scalar.fpp"

!=====================================================================!
! This module holds the abstract assembler, the system every solver
! questions: it assembles the residual and performs matrix-vector
! products on all elements in the mesh.
!
! Author: Komahan Boopathy
!=====================================================================!

module interface_assembler

  use iso_fortran_env  , only : dp => REAL64
  use class_csr        , only : csr_matrix
  use module_solve_mode, only : FORWARD, REVERSE, &
       &                        WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE, &
       &                        is_valid_mode, is_valid_part
  implicit none

  private
  public :: assembler
  public :: WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE   ! These tags are re-exported.

  !===================================================================!
  ! This type is the assembler for the physical system.
  !===================================================================!
  
  type, abstract :: assembler

     type(character(len=:)), allocatable :: description

     type(integer) :: num_state_vars
     type(logical) :: approximate_jacobian
     type(integer) :: differential_order

     ! These flags hold the transpose claim and its verification.
     ! Every REVERSE seat -
     ! transpose_product behind the unified product AND the linearized
     ! adjoint seat add_jacobian_vector_product_transpose - refuses
     ! unless the configured instance either declares its operator
     ! symmetric (an instance property, never a class property - one
     ! class serves symmetric diffusion and non-symmetric advection
     ! through its flux objects) or overrides the seat with a genuine
     ! transpose. The entry gate (converge) verifies the claim, once per
     ! system, and caches the verdict here.
     logical :: operator_is_symmetric = .false.
     logical :: transpose_verified    = .false.

     real(dp), allocatable :: S(:,:)

     ! The frozen linearization is set by linearize and cleared by
     ! clear_linearization. When the coefficients are set, the system
     ! answers the two solver questions AS the linear system - no
     ! solver ever holds linearization state:
     !
     !                        nonlinear self          frozen self
     !                        ---------------         ------------------
     !    get_residual   -->  R(x)                    b - J x
     !    product        -->  dR/du v                 J v
     !
     !    J = sum_n lin_coeff(n) dR/dU(n)   (J^T when lin_transpose -
     !        an adjoint solve is then a plain FORWARD march of the
     !        transposed frozen system; no REVERSE tag travels)
     !    b = lin_rhs when the equation brings its own (the adjoint's
     !        -df/du), else -R at the current state.
     type(scalar), allocatable :: lin_coeff(:)
     real(dp)    , allocatable :: lin_rhs(:)
     logical                   :: lin_transpose = .false.

   contains

     ! The concrete assembler must provide these deferred procedures.
     procedure(add_residual_interface)               , deferred :: add_residual
     procedure(add_jacobian_vector_product_interface), deferred :: add_jacobian_vector_product
     procedure(add_initial_condition_interface)      , deferred :: add_initial_condition

     ! These are the queries a solver makes of the system: the residual at a
     ! given state, the jacobian-vector product (mode selects forward
     ! or transpose, part selects the whole operator or a sub-part),
     ! and the inner product (provided by the system because the data
     ! distribution is the system's concern - a partitioned system
     ! reduces across images here). get_residual is PROVIDED: it lets
     ! the frozen linearization answer first, then asks the deferred
     ! state_residual - so no implementor can forget the freeze, and
     ! discretization vocabulary stays on the discretization layer.
     procedure :: get_residual
     procedure(state_residual_interface), deferred :: state_residual
     procedure :: get_jacobian_residual_product
     procedure :: inner_product

     ! These bindings freeze and thaw the linearization the two
     ! queries answer for, and provide the frozen residual that
     ! get_residual asks first.
     procedure :: linearize
     procedure :: clear_linearization
     procedure :: linearized_residual

     ! These are analytic consistency checks of the product, exact to
     ! machine precision; finite differences remain only as a coarse
     ! independent check.
     procedure :: verify_transpose_consistency
     procedure :: verify_parts_consistency

     ! The steady transpose action sits behind the REVERSE direction;
     ! the binding is public so a system with a non-symmetric operator
     ! can override it with a genuine transpose, the stated remedy of
     ! its refusal.
     procedure :: transpose_product

     ! This is the forward routing of the deferred add_ mechanism. It
     ! is public because system-layer concretes (the partitioned
     ! assembler's replicated fallback) and matrix-free checks consume
     ! it; its consolidation into the one product is a tracked
     ! deferral
     ! (ROADMAP: tracked deferrals). The base default is trivial (zero);
     ! a spatial assembler overrides it.
     procedure :: get_jacobian_vector_product
     procedure :: get_operator_csr   ! Assembled entries serve ONLY to build
                                     ! algebraic preconditioners, never to iterate.
     
     ! The assembler knows the size of the state array.
     procedure :: create_vector
     procedure :: create_state

     ! These accessors read and write the state count and the label.
     procedure :: get_num_state_vars, set_num_state_vars
     procedure :: get_description   , set_description

     ! These accessors read and write the differential order.
     procedure :: get_differential_order
     procedure :: set_differential_order

     ! The lumped mass of the semi-discrete residual is provided with
     ! a refusing default; a system that owns a mass overrides it.
     procedure :: get_lumped_mass

     ! Adjoint support carries provided defaults; the physics
     ! overrides them as needed.
     procedure :: add_jacobian_vector_product_transpose
     procedure :: get_num_design_vars
     procedure :: set_design_vars
     procedure :: get_design_vars
     procedure :: add_design_residual_transpose_product

     ! Post-processing defaults to a no-op; an assembler that owns a
     ! mesh overrides these to export the named fields.
     procedure :: write_solution_fields
     procedure :: write_gmsh_series

  end type assembler

  !===================================================================!
  ! This block declares the interfaces for the deferred procedures.
  !===================================================================!

  abstract interface

     !================================================================!
     ! This is the interface for residual assembly R(U,xi).
     !================================================================!

     impure subroutine add_residual_interface(this, residual, filter)

       import :: assembler

       class(assembler), intent(in)           :: this
       type(scalar)    , intent(inout)        :: residual(:)
       type(integer)   , intent(in), optional :: filter

     end subroutine add_residual_interface

     !================================================================!
     ! This routine returns the product of the jacobian matrix with a
     ! compatible vector:
     !
     !    pdt <---- [scalar(i)*dR(U,X)/dU(i)]*vec
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
     ! This interface supplies the initial condition to march in time.
     !================================================================!

     impure subroutine add_initial_condition_interface(this, U)

       import :: assembler

       class(assembler), intent(in)    :: this
       type(scalar)    , intent(inout) :: U(:,:)

     end subroutine add_initial_condition_interface

     !================================================================!
     ! The residual at state x is r = R(x). It is deferred: each
     ! system states its own residual in its own vocabulary, and the
     ! solver layer sees only this query. It is forward only: the
     ! adjoint right-hand side also needs the functional, and remains
     ! on the linearized path.
     !================================================================!

     ! This is the discretization's own residual at state x - the seat
     ! each system implements. Solvers reach it only through the
     ! provided get_residual, which lets a frozen linearization answer
     ! first.
     impure subroutine state_residual_interface(this, r, x)

       import :: assembler, dp

       class(assembler), intent(in)  :: this
       real(dp)        , intent(out) :: r(:)
       real(dp)        , intent(in)  :: x(:)

     end subroutine state_residual_interface

  end interface

contains
  
  !===================================================================!
  ! Return the number of state variables in the physical system.
  !===================================================================!
  
  pure type(integer) function get_num_state_vars(this)

    class(assembler), intent(in) :: this

    get_num_state_vars = this % num_state_vars

  end function get_num_state_vars

  !===================================================================!
  ! Set the number of state variables in the physical system.
  !===================================================================!
  
  pure subroutine set_num_state_vars(this, num_state_vars)

    class(assembler), intent(inout) :: this
    type(integer)  , intent(in)   :: num_state_vars

    this % num_state_vars  = num_state_vars

  end subroutine set_num_state_vars
  
  !===================================================================!
  ! Return the description set for the physical system.
  !===================================================================!
  
  pure type(character) function get_description(this)

    class(assembler), intent(in) :: this

    get_description = this % description

  end function get_description

  !===================================================================!
  ! Set the description for the physical system.
  !===================================================================!

  pure subroutine set_description(this, description)

    class(assembler), intent(inout) :: this
    type(character(len=*)), intent(in) :: description
    
    allocate(this % description, source = trim(description))

  end subroutine set_description

  !===================================================================!
  ! Return the highest order of the time derivative in the assembler.
  !===================================================================!
  
  pure type(integer) function get_differential_order(this)

    class(assembler), intent(in) :: this

    get_differential_order = this % differential_order

  end function get_differential_order

  !===================================================================!
  ! Set the highest order of the time derivative in the assembler.
  !===================================================================!

  pure subroutine set_differential_order(this, order)

    class(assembler), intent(inout) :: this
    type(integer) , intent(in)    :: order

    this % differential_order = order

  end subroutine set_differential_order

  !===================================================================!
  ! Create a state-sized vector and set its values when a scalar is
  ! supplied.
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
  ! Create the state array over the differential order and set its
  ! values when a scalar is supplied.
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
  ! The transpose jacobian-vector product
  !
  !    pdt += [scalar(i) dR/dU(i)]^T vec
  !
  ! is the seat the linearized adjoint route (newton, bdf) drives.
  ! The law: names must not lie, on every REVERSE seat. It is refused
  ! unless the configured instance declares its operator symmetric
  ! (then the forward product serves the transpose as an explicit
  ! claim) or a subclass overrides this with a genuine transpose. The
  ! genuine non-symmetric transpose inside the spatial assembler is a reserved
  ! decision, tracked in the register.
  !===================================================================!

  impure subroutine add_jacobian_vector_product_transpose(this, pdt, vec, scalars, filter)

    class(assembler), intent(in)           :: this
    type(scalar)    , intent(inout)        :: pdt(:)
    type(scalar)    , intent(in)           :: vec(:)
    type(scalar)    , intent(in)           :: scalars(:)
    type(integer)   , intent(in), optional :: filter

    if (.not. this % operator_is_symmetric) then
       error stop "assembler % add_jacobian_vector_product_transpose: no transpose " // &
            & "available - override with a genuine transpose, or set " // &
            & "operator_is_symmetric on the configured instance"
    end if

    ! The symmetric identity J^T = J is an explicit per-instance claim.
    call this % add_jacobian_vector_product(pdt, vec, scalars, filter)

  end subroutine add_jacobian_vector_product_transpose

  !===================================================================!
  ! The lumped mass of the semi-discrete residual carries one entry
  ! per dof - the M in R = M*udot - A*u + b - handed out whole so an
  ! explicit step can divide by it. The default refuses loudly; a
  ! system that owns a mass (the fvm assembler's cell volumes)
  ! overrides this.
  !===================================================================!

  pure subroutine get_lumped_mass(this, m)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: m(:)

    m = 0.0_dp
    error stop "assembler: this system provides no lumped mass"

  end subroutine get_lumped_mass

  !===================================================================!
  ! Return the number of design variables x the residual depends on.
  ! The default is none (no design dependence), so the adjoint total
  ! derivative is just the function's explicit df/dx. Physics with
  ! design variables overrides this.
  !===================================================================!

  pure type(integer) function get_num_design_vars(this)

    class(assembler), intent(in) :: this

    get_num_design_vars = 0

  end function get_num_design_vars

  !===================================================================!
  ! Set the design variables. The default assumes no design
  ! dependence, so this is a no-op. Physics carrying design variables
  ! overrides it.
  !===================================================================!

  pure subroutine set_design_vars(this, x)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(in)    :: x(:)

  end subroutine set_design_vars

  !===================================================================!
  ! Get the design variables. The default assumes no design
  ! dependence, so this returns zeros. Physics carrying design
  ! variables overrides it.
  !===================================================================!

  pure subroutine get_design_vars(this, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: x(:)

    x = 0.0_dp

  end subroutine get_design_vars

  !===================================================================!
  ! Accumulate the adjoint design contribution  dfdx += psi^T dR/dx,
  ! the product of the adjoint variables with the residual's design
  ! jacobian. The default assumes no design dependence, so this adds
  ! nothing. Physics with design variables overrides it, analytically
  ! or finite differenced as a temporary stand-in.
  !===================================================================!

  impure subroutine add_design_residual_transpose_product(this, dfdx, psi)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(inout) :: dfdx(:)
    type(scalar)    , intent(in)    :: psi(:)

  end subroutine add_design_residual_transpose_product


  !===================================================================!
  ! Export named flat-dof fields (state, adjoint state, ...) for post-
  ! processing. fields is (num_state_vars, nfield), and labels names
  ! each field. The default is a no-op, because an abstract, mesh-less
  ! system has nothing to write; a mesh-backed assembler overrides
  ! this to write a real file.
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
  ! (a gmsh view), and times gives the time of each step. meshfile is
  ! the source mesh copied verbatim and keyed by its element tags. The
  ! default is a no-op; a mesh-backed assembler overrides it. A steady
  ! export uses an nstep of one.
  !===================================================================!

  impure subroutine write_gmsh_series(this, meshfile, filename, fields, names, times)

    class(assembler), intent(in) :: this
    character(len=*), intent(in) :: meshfile, filename
    real(dp)        , intent(in) :: fields(:,:,:)
    character(len=*), intent(in) :: names(:)
    real(dp)        , intent(in) :: times(:)

  end subroutine write_gmsh_series

  !===================================================================!
  ! This is the unified jacobian-vector product. mode selects the
  ! direction (FORWARD = J v, REVERSE = J^T v); part selects the
  ! operator part (WHOLE, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE).
  ! The defaults are FORWARD and WHOLE.
  !
  ! The law of the (mode, part) pair: part names the part of the operator
  ! selected by mode - a part of J under FORWARD, a part of J^T under
  ! REVERSE. Because the lower triangle of J^T is the transposed upper
  ! triangle of J, the transpose of a forward triangle product is the
  ! REVERSE product of the OPPOSITE triangle:
  !     (forward LOWER_TRIANGLE)^T == reverse UPPER_TRIANGLE
  ! - exactly the pairing verify_transpose_consistency checks.
  !
  ! The REVERSE direction refuses unless the instance declares its
  ! operator symmetric or overrides transpose_product with a genuine
  ! transpose; the entry gate verifies either claim before the first
  ! REVERSE march.
  !===================================================================!

  impure subroutine get_jacobian_residual_product(this, w, v, mode, part)

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

    !-----------------------------------------------------------------!
    ! A wrong tag dies at the door with its name and is never silently
    ! reinterpreted; the mode and part ranges are disjoint.
    !-----------------------------------------------------------------!

    if (.not. is_valid_mode(dir)) then
       write(*,'(1x,a,i0)') "get_jacobian_residual_product: invalid mode tag ", dir
       error stop "get_jacobian_residual_product: mode must be FORWARD or REVERSE"
    end if
    if (.not. is_valid_part(sub)) then
       write(*,'(1x,a,i0)') "get_jacobian_residual_product: invalid part tag ", sub
       error stop "get_jacobian_residual_product: part must be WHOLE, DIAGONAL, " // &
            & "LOWER_TRIANGLE or UPPER_TRIANGLE"
    end if

    !-----------------------------------------------------------------!
    ! The frozen linearization answers first: w = J v with the
    ! declared coefficients. The freeze holds the direction, and a
    ! REVERSE request composes with it - the transpose of a
    ! transposed freeze is the forward action.
    !-----------------------------------------------------------------!

    if (allocated(this % lin_coeff)) then

       frozen: block

         type(scalar), allocatable :: Jv(:)
         logical :: trans

         allocate(Jv(size(w)))
         Jv    = 0.0d0
         trans = (this % lin_transpose .neqv. (dir .eq. REVERSE))

         if (sub .eq. WHOLE) then
            if (trans) then
               call this % add_jacobian_vector_product_transpose(Jv, v, this % lin_coeff)
            else
               call this % add_jacobian_vector_product(Jv, v, this % lin_coeff)
            end if
         else
            if (trans) then
               call this % add_jacobian_vector_product_transpose(Jv, v, this % lin_coeff, filter = sub)
            else
               call this % add_jacobian_vector_product(Jv, v, this % lin_coeff, filter = sub)
            end if
         end if

         w = real(Jv, dp)

       end block frozen

       return
    end if

    if (dir .eq. REVERSE) then
       !--------------------------------------------------------------!
       ! The transpose action at the steady linearization is refused
       ! unless the instance declares symmetry or overrides the
       ! transpose seat.
       !--------------------------------------------------------------!

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
  ! The residual query is provided once for every system:
  !
  !    get_residual ──▶ frozen self answers?  b - J x   (linearized)
  !                 └─▶ else the deferred     R(x)      (state_residual)
  !
  ! The freeze is structurally unforgettable - no implementor can
  ! bypass it, because solvers only ever call this seat.
  !===================================================================!

  impure subroutine get_residual(this, r, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: r(:)
    real(dp)        , intent(in)  :: x(:)

    if (this % linearized_residual(r, x)) return

    call this % state_residual(r, x)

  end subroutine get_residual

  !===================================================================!
  ! Freeze the linearization: declare the coefficients (and the
  ! external right-hand side, when the equation brings its own).
  ! transpose = .true. freezes J^T, so an adjoint solve is a plain
  ! forward march of the transposed frozen system. Freezing declares
  ! the weights, not a snapshot - the products still read the live
  ! state.
  !===================================================================!

  pure subroutine linearize(this, coeff, rhs, transpose)

    class(assembler), intent(inout)        :: this
    type(scalar)    , intent(in)           :: coeff(:)
    real(dp)        , intent(in), optional :: rhs(:)
    logical         , intent(in), optional :: transpose

    if (allocated(this % lin_coeff)) deallocate(this % lin_coeff)
    if (allocated(this % lin_rhs))   deallocate(this % lin_rhs)

    this % lin_coeff = coeff
    if (present(rhs)) this % lin_rhs = rhs

    this % lin_transpose = .false.
    if (present(transpose)) this % lin_transpose = transpose

  end subroutine linearize

  !===================================================================!
  ! Thaw: return to the nonlinear self.
  !===================================================================!

  pure subroutine clear_linearization(this)

    class(assembler), intent(inout) :: this

    if (allocated(this % lin_coeff)) deallocate(this % lin_coeff)
    if (allocated(this % lin_rhs))   deallocate(this % lin_rhs)
    this % lin_transpose = .false.

  end subroutine clear_linearization

  !===================================================================!
  ! The frozen system's residual r = b - J x is answered without
  ! touching the state. It returns .true. when it answered; every
  ! get_residual implementation asks this first, so the freeze works
  ! for any discretization.
  !===================================================================!

  impure logical function linearized_residual(this, r, x) result(answered)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: r(:)
    real(dp)        , intent(in)  :: x(:)

    type(scalar), allocatable :: res(:), Jx(:)

    answered = allocated(this % lin_coeff)
    if (.not. answered) return

    if (allocated(this % lin_rhs)) then
       r = this % lin_rhs
    else
       allocate(res(size(r)))
       res = 0.0d0
       call this % add_residual(res)
       r = -real(res, dp)
    end if

    allocate(Jx(size(r)))
    Jx = 0.0d0
    if (this % lin_transpose) then
       call this % add_jacobian_vector_product_transpose(Jx, x, this % lin_coeff)
    else
       call this % add_jacobian_vector_product(Jx, x, this % lin_coeff)
    end if
    r = r - real(Jx, dp)

  end function linearized_residual

  !===================================================================!
  ! The steady transpose action w = J^T v sits behind the REVERSE
  ! direction. The law: names must not lie. The base refuses unless
  ! the configured instance declares its operator symmetric - then
  ! J^T = J is an explicit claim and the forward product serves it.
  ! A system with a
  ! non-symmetric operator must override this with a genuine transpose.
  ! Either claim is verified by the entry gate (converge) through
  ! verify_transpose_consistency before the first REVERSE march.
  !===================================================================!

  impure subroutine transpose_product(this, w, v, sub)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: w(:)
    real(dp)        , intent(in)  :: v(:)
    integer         , intent(in)  :: sub

    if (.not. is_valid_part(sub)) then
       write(*,'(1x,a,i0)') "transpose_product: invalid part tag ", sub
       error stop "transpose_product: part must be WHOLE, DIAGONAL, " // &
            & "LOWER_TRIANGLE or UPPER_TRIANGLE"
    end if

    if (.not. this % operator_is_symmetric) then
       error stop "assembler % transpose_product: no transpose available - " // &
            & "override transpose_product with a genuine transpose, or set " // &
            & "operator_is_symmetric on the configured instance; the REVERSE " // &
            & "entry gate verifies whichever claim is made"
    end if

    ! The symmetric identity J^T = J is an explicit per-instance claim.
    if (sub .eq. WHOLE) then
       call this % get_jacobian_vector_product(w, v)
    else
       call this % get_jacobian_vector_product(w, v, filter = sub)
    end if

  end subroutine transpose_product

  !===================================================================!
  ! The system owns the inner product of two vectors of its space,
  ! because distribution is the system's business: this serial default
  ! is the plain dot; a partitioned system sums its own rows and
  ! reduces across images. The volume weighting attaches here later.
  !===================================================================!

  impure real(dp) function inner_product(this, a, b)

    class(assembler), intent(in) :: this
    real(dp)        , intent(in) :: a(:)
    real(dp)        , intent(in) :: b(:)

    inner_product = dot_product(a, b)

  end function inner_product

  !===================================================================!
  ! This consistency check tests <w, J v> = <J^T w, v> for
  ! deterministic pseudo-random v and w, per part. Two products and
  ! two inner products suffice - exact to machine precision, with no
  ! truncation error. It returns the largest relative defect over the
  ! parts.
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
       ! The transpose of the lower triangle is the upper triangle, so pair them.
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
  ! This consistency check tests (diagonal + lower + upper) v = whole
  ! v for a deterministic pseudo-random v. Three part-products stand
  ! against one whole-product, catching part-implementation errors at
  ! machine precision. It returns the relative defect.
  !===================================================================!

  impure real(dp) function verify_parts_consistency(this) result(defect)

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

  end function verify_parts_consistency

  !===================================================================!
  ! Fill deterministically with a linear congruential pseudo-random
  ! sequence, so the checks are reproducible run to run.
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
  ! The default operator action Aq = A q goes through the deferred
  ! jacobian-vector product at the steady linearization (dR/du only).
  !===================================================================!

  pure subroutine get_jacobian_vector_product(this, Aq, q, filter)

    class(assembler) , intent(in)           :: this
    real(dp)         , intent(out)          :: Aq(:)
    real(dp)         , intent(in)           :: q(:)
    integer          , intent(in), optional :: filter

    Aq = 0.0_dp

  end subroutine get_jacobian_vector_product

  !===================================================================!
  ! By default, this assembler does not assemble a sparse operator. A
  ! spatial assembler (class_assembler) overrides this with the real
  ! csr.
  !===================================================================!

  impure subroutine get_operator_csr(this, A)

    class(assembler), intent(in)  :: this
    type(csr_matrix), intent(out) :: A

    error stop "get_operator_csr: not implemented for this assembler"

  end subroutine get_operator_csr

end module interface_assembler
