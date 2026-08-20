!=====================================================================!
! The conduction law: the material's answer to a gradient.
!
! LEVEL 3 OF THE STRATIFICATION. A conduction law holds one tensor,
!
!      K = | kxx kxy kxz |
!          | kyx kyy kyz |         isotropic k is K = k * I
!          | kzx kzy kzz |
!
! and supplies two coefficient arrays, one entry per face, from the
! mesh's own normals and areas:
!
!      normal_conductivity   keff_e = n^T K n    what the old flux
!                                                called the normal
!                                                diffusivity
!      edge_coefficients     keff_e * area_e     the dictionary's
!                                                interior diffusion
!                                                coefficient, zero on
!                                                the headless faces -
!                                                a boundary face gets
!                                                its coefficient from
!                                                its condition, not
!                                                from the material
!
! It owns no operator and no balance; it computes numbers.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conduction

  use iso_fortran_env  , only : dp => REAL64
  use field_stored, only : stored_field
  use view_mesh , only : mesh

  implicit none

  private
  public :: conduction

  type :: conduction

     real(dp) :: k(3, 3) = 0.0_dp

   contains

     procedure :: normal_conductivity
     procedure :: edge_coefficients

  end type conduction

  interface conduction
     module procedure isotropic, tensor
  end interface conduction

contains

  !===================================================================!
  ! The two ways to state the law: one number for an isotropic
  ! material, the full tensor otherwise.
  !===================================================================!

  pure type(conduction) function isotropic(k) result(this)

    real(dp), intent(in) :: k

    integer :: i

    this % k = 0.0_dp
    do i = 1, 3
       this % k(i, i) = k
    end do

  end function isotropic

  pure type(conduction) function tensor(k) result(this)

    real(dp), intent(in) :: k(3, 3)

    this % k = k

  end function tensor

  !===================================================================!
  ! keff_e = n^T K n, one entry per face, from the mesh's normals.
  !===================================================================!

  subroutine normal_conductivity(this, m, values)

    class(conduction), intent(in)      :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: values(:)

    type(stored_field) :: fn
    real(dp), allocatable :: normals(:)
    real(dp) :: n(3)
    integer :: ne, e

    fn = m % face_normal()
    call fn % get_real_vector(normals)

    ne = m % num_edges()
    allocate(values(ne))

    do e = 1, ne
       n = normals(3 * e - 2 : 3 * e)
       values(e) = dot_product(n, matmul(this % k, n))
    end do

  end subroutine normal_conductivity

  !===================================================================!
  ! The dictionary's interior coefficient: keff_e * area_e, zero on
  ! the headless faces.
  !===================================================================!

  subroutine edge_coefficients(this, m, values)

    class(conduction), intent(in)      :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: values(:)

    type(stored_field) :: fa
    real(dp), allocatable :: areas(:)
    integer :: e

    call this % normal_conductivity(m, values)

    fa = m % face_area()
    call fa % get_real_vector(areas)

    do e = 1, size(values)
       if (m % edge_has_head(e)) then
          values(e) = values(e) * areas(e)
       else
          values(e) = 0.0_dp
       end if
    end do

  end subroutine edge_coefficients

end module operation_conduction
