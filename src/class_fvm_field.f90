#include "scalar.fpp"

!=====================================================================!
! Cell-centred finite-volume field: the discretization that turns a
! flat state vector U (one cell-average per variable per cell) into the
! pointwise (q, grad q) the physics needs. Implements interface_field;
! the only method on the FV solve path is face_state, which reconstructs
! the face value and the FACE-NORMAL gradient between two cells:
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
     ! the FV reconstruction the assembler uses
     procedure :: face_state

     ! interface_field deferred set (stubbed except where trivial)
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
  ! Construct an FV field of m components over n state dofs
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
  ! grid/graph are passed in (the assembler owns them) to avoid any
  ! aliasing when the assembler is cloned into integrators/adjoints.
  ! iface is the LOCAL face index of icell (for the normal); gface is
  ! the GLOBAL face id (for area/delta/centre/neighbour).
  !===================================================================!

  pure subroutine face_state(this, grid, g, U, icell, iface, gface, st)

    class(fvm_field) , intent(in)  :: this
    class(mesh)      , intent(in)  :: grid
    class(graph)     , intent(in)  :: g
    real(dp)         , intent(in)  :: U(:)
    integer          , intent(in)  :: icell, iface, gface
    type(point_state), intent(out) :: st

    integer  :: nv, j, ncell, fcells(2), p, n
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

       fcells(1:2) = grid % face_cells(1:2, gface)
       if (fcells(1) .eq. icell) then
          ncell = fcells(2)
       else
          ncell = fcells(1)
       end if

       do j = 1, nv
          p = g % dof(icell, j)
          n = g % dof(ncell, j)
          dphidn         = (U(n) - U(p))/fdelta
          st % gradq(:,j) = dphidn*nf            ! normal-component gradient
          st % q(j)       = 0.5_dp*(U(p) + U(n)) ! face value (midpoint)
       end do

    else

       ! boundary face: the bc closure (lhs/rhs_coeff) supplies the flux;
       ! return the owner value with no reconstructed normal gradient.
       do j = 1, nv
          st % q(j)       = U(g % dof(icell, j))
          st % gradq(:,j) = 0.0_dp
       end do

    end if interior

  end subroutine face_state

  !===================================================================!
  ! interface_field deferred procedures - FEM/UQ extensions, not on the
  ! FV solve path. Stubbed (pure ones return zero; the rest error stop).
  !===================================================================!

  pure function evaluate(this, t, x, y, z) result(u)
    class(fvm_field), intent(in) :: this
    real(dp)        , intent(in) :: t, x(:), y(:), z(:)
    real(dp)                     :: u(this % num_components)
    u = 0.0_dp
  end function evaluate

  impure subroutine discretize(this, U)
    class(fvm_field)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: U(:)
    error stop "fvm_field % discretize: not used by the fv solve"
  end subroutine discretize

  impure subroutine reconstruct(this, U)
    class(fvm_field), intent(inout) :: this
    real(dp)        , intent(in)    :: U(:)
    error stop "fvm_field % reconstruct: not used by the fv solve"
  end subroutine reconstruct

  impure subroutine project(this, U)
    class(fvm_field)     , intent(in)  :: this
    real(dp), allocatable, intent(out) :: U(:)
    error stop "fvm_field % project: not used by the fv solve"
  end subroutine project

  impure subroutine remainder(this, U, e)
    class(fvm_field)         , intent(in)  :: this
    real(dp)                 , intent(in)  :: U(:)
    class(field), allocatable, intent(out) :: e
    error stop "fvm_field % remainder: not used by the fv solve"
  end subroutine remainder

  pure function basis(this, i, x) result(phi)
    class(fvm_field), intent(in) :: this
    integer         , intent(in) :: i
    real(dp)        , intent(in) :: x(:)
    real(dp)                     :: phi
    phi = 0.0_dp
  end function basis

  impure subroutine grad(this, idim, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: idim
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % grad: not used by the fv solve"
  end subroutine grad

  impure subroutine ddt(this, du)
    class(fvm_field)         , intent(in)  :: this
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % ddt: not used by the fv solve"
  end subroutine ddt

  impure subroutine drandom(this, k, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: k
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % drandom: not used by the fv solve"
  end subroutine drandom

  impure subroutine ddesign(this, k, du)
    class(fvm_field)         , intent(in)  :: this
    integer                  , intent(in)  :: k
    class(field), allocatable, intent(out) :: du
    error stop "fvm_field % ddesign: not used by the fv solve"
  end subroutine ddesign

  pure function inner_product(this, other) result(ip)
    class(fvm_field), intent(in) :: this
    class(field)    , intent(in) :: other
    real(dp)                     :: ip
    ip = 0.0_dp
  end function inner_product

  pure function norm(this) result(nrm)
    class(fvm_field), intent(in) :: this
    real(dp)                     :: nrm
    nrm = 0.0_dp
  end function norm

end module class_fvm_field
