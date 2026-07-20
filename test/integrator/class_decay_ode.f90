#include "scalar.fpp"

!=====================================================================!
! Toy first-order system  R = qdot + lambda*q = 0,  q(0) = 1, whose
! exact solution is q(t) = exp(-lambda t). A minimal class(assembler)
! used to verify the time integrator: it holds its state in the
! inherited S(nvars, order+1) and supplies the residual, the
! matrix-free jacobian-vector product and the initial condition.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_decay_ode

  use interface_assembler, only : assembler

  implicit none

  private
  public :: decay_ode

  type, extends(assembler) :: decay_ode

     real(8) :: lambda = 1.0d0

   contains

     procedure :: state_residual
     procedure :: add_residual
     procedure :: add_jacobian_vector_product
     procedure :: add_initial_condition

  end type decay_ode

  interface decay_ode
     module procedure create
  end interface decay_ode

contains

  !===================================================================!
  ! Construct the decay system with the given rate lambda
  !===================================================================!

  type(decay_ode) function create(lambda) result(this)

    real(8), intent(in) :: lambda

    this % lambda = lambda

    call this % set_num_state_vars(1)
    call this % set_differential_order(1)

    allocate(this % S(1, 2))
    this % S = 0.0d0

  end function create

  !===================================================================!
  ! The steady residual at state x (qdot = 0):  r = -lambda * x
  ! (the deferred seat behind the provided get_residual)
  !===================================================================!

  impure subroutine state_residual(this, r, x)

    class(decay_ode), intent(in)  :: this
    real(8)         , intent(out) :: r(:)
    real(8)         , intent(in)  :: x(:)

    r = -this % lambda * x

  end subroutine state_residual

  !===================================================================!
  ! Residual  R = qdot + lambda*q   (S(1,1) = q, S(1,2) = qdot)
  !===================================================================!

  subroutine add_residual(this, residual, filter)

    class(decay_ode), intent(in)           :: this
    type(scalar)    , intent(inout)        :: residual(:)
    integer         , intent(in), optional :: filter

    residual(1) = residual(1) + this % S(1,2) + this % lambda*this % S(1,1)

  end subroutine add_residual

  !===================================================================!
  ! pdt += [scalars(1) dR/dq + scalars(2) dR/dqdot] vec
  !      =  [scalars(1)*lambda + scalars(2)] vec
  !===================================================================!

  subroutine add_jacobian_vector_product(this, pdt, vec, scalars, filter)

    class(decay_ode), intent(in)           :: this
    type(scalar)    , intent(inout)        :: pdt(:)
    type(scalar)    , intent(in)           :: vec(:)
    type(scalar)    , intent(in)           :: scalars(:)
    integer         , intent(in), optional :: filter

    pdt(1) = pdt(1) + (scalars(1)*this % lambda + scalars(2))*vec(1)

  end subroutine add_jacobian_vector_product

  !===================================================================!
  ! Initial condition  q(0) = 1,  qdot(0) = -lambda*q(0)
  !===================================================================!

  subroutine add_initial_condition(this, U)

    class(decay_ode), intent(in)    :: this
    type(scalar)    , intent(inout) :: U(:,:)

    U(1,1) = 1.0d0
    U(1,2) = -this % lambda

  end subroutine add_initial_condition

end module class_decay_ode
