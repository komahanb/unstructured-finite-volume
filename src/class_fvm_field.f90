#include "scalar.fpp"

!=====================================================================!
! The cell-centred finite-volume field is the discretization that turns
! a flat state vector U (one cell-average per variable per cell) into
! the pointwise (q, grad q) the physics needs. It implements
! interface_field; the only method on the FV solve path is face_state,
! which reconstructs the face value and the FACE-NORMAL gradient
! between two cells:
!
!     grad q(:,j) = ((q_n - q_p)/fdelta) * n
!
! so that a flux F = -K grad q dotted with n gives exactly the two-point
! normal diffusive flux keff*(q_n - q_p)/fdelta. The tangential / non-
! orthogonal part is a separate deferred correction in the assembler
! (get_skew_source), so this reconstruction is deliberately the normal
! component only - a full least-squares gradient here would change the
! operator and break the operator-split and order-of-accuracy oracles.
!
! The remaining interface_field operators (basis, project, ddt, random/
! design sensitivities, ...) are FEM/UQ extensions, not used by the FV
! solve, and are stubbed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_fvm_field

  use iso_fortran_env  , only : dp => REAL64
  use interface_field  , only : field
  use interface_physics, only : point_state
  use class_mesh       , only : mesh
  use interface_graph  , only : graph

  implicit none

  private
  public :: fvm_field

  type, extends(field) :: fvm_field
   contains
     ! This is the FV reconstruction the assembler uses.
     procedure :: face_state

     ! The interface_field deferred set is stubbed except where trivial.
     procedure :: evaluate
     procedure :: discretize
     procedure :: reconstruct
     procedure :: project
     procedure :: remainder
     procedure :: basis
     procedure :: grad
     procedure :: ddt
     procedure :: drandom
     procedure :: ddesign
     procedure :: inner_product
     procedure :: norm
  end type fvm_field

  interface fvm_field
     module procedure create
  end interface fvm_field

contains

  !===================================================================!
  ! Construct an FV field of m components over n state dofs.
  !===================================================================!

  pure type(fvm_field) function create(num_components, num_state_vars) result(this)
    integer, intent(in) :: num_components
    integer, intent(in) :: num_state_vars
    this % num_components = num_components
    this % num_state_vars = num_state_vars
    this % num_random_dim = 0
    this % num_design_dim = 0
  end function create

  !===================================================================!
  ! Reconstruct (q, grad q) at a face from the cell-average state U.
  ! The grid and the graph are passed in (the assembler owns them) to
  ! avoid any aliasing when the assembler is cloned into integrators
  ! and adjoints. Here iface is the LOCAL face index of icell (for the
  ! normal); gface is the GLOBAL face id (for area/delta/centre/
  ! neighbour).
  !===================================================================!

  pure subroutine face_state(this, grid, g, U, icell, iface, gface, st)

    class(fvm_field) , intent(in)  :: this
    class(mesh)      , intent(in)  :: grid
    class(graph)     , intent(in)  :: g
    real(dp)         , intent(in)  :: U(:)
    integer          , intent(in)  :: icell, iface, gface
    type(point_state), intent(out) :: st

    integer  :: nv, j, ncell, p, n
    real(dp) :: nf(3), fdelta, dphidn

    nv = this % num_components

    allocate(st % q(nv), st % gradq(3, nv))
    st % nv    = nv
    st % q     = 0.0_dp
    st % gradq = 0.0_dp

    nf     = grid % cell_face_normals(1:3, iface, icell)
    fdelta = grid % face_deltas(gface)
    st % x = grid % face_centers(1:3, gface)

    interior: if (grid % num_face_cells(gface) .eq. 2) then

       ncell = grid % across(icell, gface)

       do j = 1, nv
          p = g % dof(icell, j)
          n = g % dof(ncell, j)
          dphidn         = (U(n) - U(p))/fdelta
          st % gradq(:,j) = dphidn*nf            ! The normal-component gradient.
          st % q(j)       = 0.5_dp*(U(p) + U(n)) ! The face value (midpoint).
       end do

    else

       !--------------------------------------------------------------!
       ! At a boundary face the bc closure (lhs/rhs_coeff) supplies
       ! the flux; return the owner value with no reconstructed
       ! normal gradient.
       !--------------------------------------------------------------!

       do j = 1, nv
          st % q(j)       = U(g % dof(icell, j))
          st % gradq(:,j) = 0.0_dp
       end do

    end if interior

  end subroutine face_state

  !===================================================================!
  ! The remaining interface_field deferred procedures are FEM/UQ
  ! extensions, not on the FV solve path. They are stubbed: the pure
  ! ones return zero and the rest error stop.
  !===================================================================!

  !===================================================================!
  ! The evaluate stub returns zero.
  !===================================================================!

  pure function evaluate(this, t, x, y, z) result(u)
    class(fvm_field), intent(in) :: this
    real(dp)        , intent(in) :: t, x(:), y(:), z(:)
    real(dp)                     :: u(this % num_components)
    u = 0.0_dp
  end function evaluate

  !===================================================================!
  ! The discretize stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine discretize(this, U)
    class(fvm_field)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: U(:)
    error stop "fvm_field % discretize: not used by the fv solve"
  end subroutine discretize

  !===================================================================!
  ! The reconstruct stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine reconstruct(this, U)
    class(fvm_field), intent(inout) :: this
    real(dp)        , intent(in)    :: U(:)
    error stop "fvm_field % reconstruct: not used by the fv solve"
  end subroutine reconstruct

  !===================================================================!
  ! The project stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine project(this, U)
    class(fvm_field)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: U(:)
    error stop "fvm_field % project: not used by the fv solve"
  end subroutine project

  !===================================================================!
  ! The remainder stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine remainder(this, U, e)
    class(fvm_field)         , intent(in)  :: this
    real(dp)                 , intent(in)  :: U(:)
    class(field), allocatable, intent(out) :: e
    error stop "fvm_field % remainder: not used by the fv solve"
  end subroutine remainder

  !===================================================================!
  ! The basis stub returns zero.
  !===================================================================!

  pure function basis(this, i, x) result(phi)
    class(fvm_field), intent(in) :: this
    integer         , intent(in) :: i
    real(dp)        , intent(in) :: x(:)
    real(dp)                     :: phi
    phi = 0.0_dp
  end function basis

  !===================================================================!
  ! The grad stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine grad(this, idim, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: idim
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % grad: not used by the fv solve"
  end subroutine grad

  !===================================================================!
  ! The ddt stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine ddt(this, du)
    class(fvm_field)         , intent(in)  :: this
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % ddt: not used by the fv solve"
  end subroutine ddt

  !===================================================================!
  ! The drandom stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine drandom(this, k, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: k
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % drandom: not used by the fv solve"
  end subroutine drandom

  !===================================================================!
  ! The ddesign stub stops; it is not used by the FV solve.
  !===================================================================!

  impure subroutine ddesign(this, k, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: k
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % ddesign: not used by the fv solve"
  end subroutine ddesign

  !===================================================================!
  ! The inner_product stub returns zero.
  !===================================================================!

  pure function inner_product(this, other) result(ip)
    class(fvm_field), intent(in) :: this
    class(field)    , intent(in) :: other
    real(dp)                     :: ip
    ip = 0.0_dp
  end function inner_product

  !===================================================================!
  ! The norm stub returns zero.
  !===================================================================!

  pure function norm(this) result(nrm)
    class(fvm_field), intent(in) :: this
    real(dp)                     :: nrm
    nrm = 0.0_dp
  end function norm

end module class_fvm_field
