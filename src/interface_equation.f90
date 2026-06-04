!=====================================================================!
! The physics behind the discretization. The assembler is told nothing
! about the pde - it asks the equation, per cell, for
!
!   - the diffusion tensor K(x)   [3x3, anisotropic]
!   - the volumetric source s(x)
!   - how many variables live on a cell
!
! Concrete equations (diffusion, advection-diffusion, ...) extend this.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_equation

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: equation

  type, abstract :: equation

     integer :: num_variables = 1

   contains

     procedure(diffusion_tensor_interface), deferred :: diffusion_tensor
     procedure(source_interface)          , deferred :: source

     ! Design variables for sensitivity analysis (provided defaults; a
     ! concrete equation overrides them to expose its parameters)
     procedure :: get_num_design_vars
     procedure :: set_design_vars
     procedure :: get_design_vars

  end type equation

  interface

     ! Diffusion tensor K at a spatial point
     pure function diffusion_tensor_interface(this, x) result(K)
       import :: equation, dp
       class(equation), intent(in) :: this
       real(dp)       , intent(in) :: x(3)
       real(dp)                    :: K(3,3)
     end function diffusion_tensor_interface

     ! Volumetric source at a spatial point
     pure function source_interface(this, x) result(s)
       import :: equation, dp
       class(equation), intent(in) :: this
       real(dp)       , intent(in) :: x(3)
       real(dp)                    :: s
     end function source_interface

  end interface

contains

  !===================================================================!
  ! Number of design variables this equation exposes. Default: none.
  !===================================================================!

  pure integer function get_num_design_vars(this)

    class(equation), intent(in) :: this

    get_num_design_vars = 0

  end function get_num_design_vars

  !===================================================================!
  ! Set / get the design variables. Default: no design dependence.
  !===================================================================!

  subroutine set_design_vars(this, x)

    class(equation), intent(inout) :: this
    real(dp)       , intent(in)    :: x(:)

  end subroutine set_design_vars

  subroutine get_design_vars(this, x)

    class(equation), intent(in)  :: this
    real(dp)       , intent(out) :: x(:)

  end subroutine get_design_vars

end module interface_equation
