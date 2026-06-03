!=====================================================================!
! Abstract interface for a field u(t, x, y, z).
!
! A field is the continuous solution over the spatial domain x; time t
! evolves it, the random vector y makes it a random field, and the
! design vector z parametrises a family of fields. Discretising the
! spatial domain turns the field into its state - the finite vector of
! coefficients U in a basis {phi_i} of the resolved space V_h:
!
!     u  =  R U  +  e         (reconstruct the state, plus a remainder)
!        =  Pi u + (I - Pi) u
!
!   discretize  : u    -> U           Delta : V   -> R^N   (restriction)
!   reconstruct : U    -> u_h = R U   R     : R^N -> V_h   (prolongation)
!   project     : u    -> Pi u        Pi    : V   -> V_h   (orthogonal)
!   remainder   : u, U -> e = (I-Pi)u                      (in V_h^perp)
!
! Concrete fields (cell-centred fvm, nodal fem, spectral, ...) extend
! this and supply the operators.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_field

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: field

  !-------------------------------------------------------------------!
  ! Abstract field type
  !-------------------------------------------------------------------!

  type, abstract :: field

     ! Codomain: m = 1 scalar field, m > 1 vector field
     integer :: num_components

     ! Sizes of the (discretised) domains
     integer :: num_state_vars    ! N : spatial dofs of the state U
     integer :: num_random_dim    ! dim(y) : stochastic domain
     integer :: num_design_dim    ! dim(z) : design domain

   contains

     ! u(t, x, y, z) - the field as a function of its arguments
     procedure(evaluate_interface)     , deferred :: evaluate

     ! Field <-> state
     procedure(discretize_interface)   , deferred :: discretize   ! Delta : u -> U
     procedure(reconstruct_interface)  , deferred :: reconstruct  ! R     : U -> u
     procedure(project_interface)      , deferred :: project      ! Pi    : u -> Pi u
     procedure(remainder_interface)    , deferred :: remainder    ! e = (I - Pi) u
     procedure(basis_interface)        , deferred :: basis        ! phi_i(x)

     ! Differentiation with respect to each argument
     procedure(gradient_interface)     , deferred :: grad         ! d u / d x_idim (spatial)
     procedure(time_deriv_interface)   , deferred :: ddt          ! d u / d t      (temporal)
     procedure(sensitivity_interface)  , deferred :: drandom      ! d u / d y_k    (stochastic)
     procedure(sensitivity_interface)  , deferred :: ddesign      ! d u / d z_k    (design)

     ! L2 structure (gives the norm and Galerkin orthogonality e _|_ V_h)
     procedure(inner_product_interface), deferred :: inner_product
     procedure(norm_interface)         , deferred :: norm

  end type field

  abstract interface

     !================================================================!
     ! Evaluate u(t, x, y, z); returns the m components of the field
     !================================================================!

     pure function evaluate_interface(this, t, x, y, z) result(u)
       import :: field, dp
       class(field), intent(in) :: this
       real(dp)    , intent(in) :: t
       real(dp)    , intent(in) :: x(:)    ! spatial
       real(dp)    , intent(in) :: y(:)    ! stochastic
       real(dp)    , intent(in) :: z(:)    ! design
       real(dp)                 :: u(this % num_components)
     end function evaluate_interface

     !================================================================!
     ! Delta : restrict the field to its state coefficients U in R^N
     !================================================================!

     subroutine discretize_interface(this, U)
       import :: field, dp
       class(field)         , intent(in)  :: this
       real(dp), allocatable, intent(out) :: U(:)
     end subroutine discretize_interface

     !================================================================!
     ! R : reconstruct the field from a state, u_h = sum_i U_i phi_i
     !================================================================!

     subroutine reconstruct_interface(this, U)
       import :: field, dp
       class(field), intent(inout) :: this
       real(dp)    , intent(in)    :: U(:)
     end subroutine reconstruct_interface

     !================================================================!
     ! Pi : orthogonal projection of the field onto V_h (resolved state)
     !================================================================!

     subroutine project_interface(this, U)
       import :: field, dp
       class(field)         , intent(in)  :: this
       real(dp), allocatable, intent(out) :: U(:)
     end subroutine project_interface

     !================================================================!
     ! e = (I - Pi) u : the unresolved remainder for the state U
     !================================================================!

     subroutine remainder_interface(this, U, e)
       import :: field, dp
       class(field)             , intent(in)  :: this
       real(dp)                 , intent(in)  :: U(:)
       class(field), allocatable, intent(out) :: e
     end subroutine remainder_interface

     !================================================================!
     ! phi_i : the i-th basis function of V_h at a spatial point
     !================================================================!

     pure function basis_interface(this, i, x) result(phi)
       import :: field, dp
       class(field), intent(in) :: this
       integer     , intent(in) :: i
       real(dp)    , intent(in) :: x(:)
       real(dp)                 :: phi
     end function basis_interface

     !================================================================!
     ! d u / d x_idim : spatial gradient component (a field)
     !================================================================!

     subroutine gradient_interface(this, idim, du)
       import :: field
       class(field)             , intent(in)  :: this
       integer                  , intent(in)  :: idim
       class(field), allocatable, intent(out) :: du
     end subroutine gradient_interface

     !================================================================!
     ! d u / d t : temporal rate (a field)
     !================================================================!

     subroutine time_deriv_interface(this, du)
       import :: field
       class(field)             , intent(in)  :: this
       class(field), allocatable, intent(out) :: du
     end subroutine time_deriv_interface

     !================================================================!
     ! Sensitivity to the k-th random (y_k) or design (z_k) variable
     !================================================================!

     subroutine sensitivity_interface(this, k, du)
       import :: field
       class(field)             , intent(in)  :: this
       integer                  , intent(in)  :: k
       class(field), allocatable, intent(out) :: du
     end subroutine sensitivity_interface

     !================================================================!
     ! L2 inner product <this, other>; norm follows as sqrt(<u,u>)
     !================================================================!

     pure function inner_product_interface(this, other) result(ip)
       import :: field, dp
       class(field), intent(in) :: this
       class(field), intent(in) :: other
       real(dp)                 :: ip
     end function inner_product_interface

     pure function norm_interface(this) result(nrm)
       import :: field, dp
       class(field), intent(in) :: this
       real(dp)                 :: nrm
     end function norm_interface

  end interface

end module interface_field
