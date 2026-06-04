#include "scalar.fpp"

!=====================================================================!
! Module that contains class that assembles residual and performs
! matrix vector products on all elements in mesh
!
! Author: Komahan Boopathy
!=====================================================================!

module interface_assembler

  use iso_fortran_env  , only : dp => REAL64
  implicit none
  
  private
  public :: assembler
 
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

  subroutine create_vector(this, x, val)
    
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

  subroutine create_state(this, S, val)
    
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

  subroutine set_design_vars(this, x)

    class(assembler), intent(inout) :: this
    real(dp)        , intent(in)    :: x(:)

  end subroutine set_design_vars

  subroutine get_design_vars(this, x)

    class(assembler), intent(in)  :: this
    real(dp)        , intent(out) :: x(:)

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

end module interface_assembler
