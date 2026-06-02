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

end module interface_equation
