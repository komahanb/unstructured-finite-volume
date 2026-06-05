#include "scalar.fpp"

!=====================================================================!
! Continuous, pointwise, discretization-agnostic physics operators.
!
! A conservation law is written as
!
!     dq/dt + div F(q, grad q) = S(q, grad q)
!
! and a function of interest as  J = integral f(q, grad q) dV. The three
! operators - flux F, source S, objective f - are pointwise functions of
! the state q and its gradient grad q (first derivatives only; the
! divergence theorem / integration by parts removes the second). They
! provide their value and their state/design partials. They know nothing
! about meshes or cells: a finite-volume assembler integrates F over
! faces and S over the volume; a finite-element assembler would weight
! the same operators differently. This is the seam that keeps the
! framework law-agnostic (and discretization-agnostic).
!
! flux, source and objective extend a common base `physics`. Fortran has
! single inheritance, so a concrete law supplies a flux object and a
! source object separately (composed at setup), not one type that is
! both.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_physics

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: physics, source, objective, point_state

  !-------------------------------------------------------------------!
  ! Pointwise evaluation context: the state and its gradient at a point
  ! (a face centroid for the flux, a cell centroid for the source). The
  ! state is type(scalar) so complex-step partials drop in later; the
  ! spatial point x is geometry (real).
  !-------------------------------------------------------------------!

  type :: point_state
     integer                   :: nv = 1       ! number of variables
     type(scalar), allocatable :: q(:)         ! q        (nv)
     type(scalar), allocatable :: gradq(:,:)   ! grad q   (3, nv)
     real(dp)                  :: x(3) = 0.0_dp ! spatial point
  end type point_state

  !-------------------------------------------------------------------!
  ! Common base: number of components it acts on + design variables.
  !-------------------------------------------------------------------!

  type, abstract :: physics
     integer :: num_components = 1
   contains
     procedure :: num_design_vars => physics_num_design_vars  ! default 0
     procedure :: set_design_vars => physics_set_design_vars  ! default no-op
     procedure :: get_design_vars => physics_get_design_vars  ! default no-op
  end type physics

  !-------------------------------------------------------------------!
  ! Volumetric source S(q, grad q): value -> S(nv) + partials.
  !-------------------------------------------------------------------!

  type, extends(physics), abstract :: source
   contains
     procedure(source_value_interface), deferred :: value
     procedure :: dsource_dq      => source_dq_zero           ! (nv,nv)
     procedure :: dsource_dgradq  => source_dgradq_zero       ! (nv,3,nv)
     procedure :: dsource_ddesign => source_ddesign_zero      ! (nv) for design var k
  end type source

  !-------------------------------------------------------------------!
  ! Function-of-interest integrand f(q, grad q): value -> scalar + parts.
  !-------------------------------------------------------------------!

  type, extends(physics), abstract :: objective
   contains
     procedure(objective_value_interface), deferred :: value
     procedure :: dobj_dq      => objective_dq_zero           ! (nv)
     procedure :: dobj_dgradq  => objective_dgradq_zero       ! (3,nv)
     procedure :: dobj_ddesign => objective_ddesign_zero      ! scalar for design var k
  end type objective

  !-------------------------------------------------------------------!
  ! Deferred interfaces
  !-------------------------------------------------------------------!

  abstract interface

     pure function source_value_interface(this, st) result(S)
       import :: source, point_state
       class(source)    , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: S(this % num_components)
     end function source_value_interface

     pure function objective_value_interface(this, st) result(f)
       import :: objective, point_state
       class(objective) , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: f
     end function objective_value_interface

  end interface

contains

  !===================================================================!
  ! Design-variable defaults (no design dependence)
  !===================================================================!

  pure integer function physics_num_design_vars(this) result(n)
    class(physics), intent(in) :: this
    n = 0
  end function physics_num_design_vars

  subroutine physics_set_design_vars(this, x)
    class(physics), intent(inout) :: this
    real(dp)      , intent(in)    :: x(:)
  end subroutine physics_set_design_vars

  subroutine physics_get_design_vars(this, x)
    class(physics), intent(in)  :: this
    real(dp)      , intent(out) :: x(:)
  end subroutine physics_get_design_vars

  !===================================================================!
  ! Default partials: zero (overridden by laws that depend on them)
  !===================================================================!

  pure function source_dq_zero(this, st) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: dS(this % num_components, this % num_components)
    dS = 0.0_dp
  end function source_dq_zero

  pure function source_dgradq_zero(this, st) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: dS(this % num_components, 3, this % num_components)
    dS = 0.0_dp
  end function source_dgradq_zero

  pure function source_ddesign_zero(this, st, k) result(dS)
    class(source)    , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: dS(this % num_components)
    dS = 0.0_dp
  end function source_ddesign_zero

  pure function objective_dq_zero(this, st) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: df(this % num_components)
    df = 0.0_dp
  end function objective_dq_zero

  pure function objective_dgradq_zero(this, st) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    type(scalar)                  :: df(3, this % num_components)
    df = 0.0_dp
  end function objective_dgradq_zero

  pure function objective_ddesign_zero(this, st, k) result(df)
    class(objective) , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: df
    df = 0.0_dp
  end function objective_ddesign_zero

end module interface_physics
