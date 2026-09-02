!=====================================================================!
! The conduction law: the flux a material produces from a gradient.
!
! LEVEL 3 OF THE STRATIFICATION. A conduction law stores one tensor,
!
!      K = | kxx kxy kxz |
!          | kyx kyy kyz |         isotropic k is K = k * I
!          | kzx kzy kzz |
!
! and supplies two coefficient arrays, one entry per face, from the
! mesh's own normals and areas:
!
!      normal_conductivity   keff_e = n^T K n    the normal
!                                                diffusivity of the
!                                                earlier flux code
!      edge_coefficients     keff_e * area_e     the dictionary's
!                                                interior diffusion
!                                                coefficient, zero on
!                                                the headless faces -
!                                                a boundary face takes
!                                                its coefficient from
!                                                its condition, not
!                                                from the material
!
! The law stores no operator and no balance; it computes coefficients.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conduction

  use util_precision  , only : dp
  use field_stored, only : stored_field
  use view_mesh , only : mesh

  implicit none

  private
  public :: conduction

  type :: conduction

     ! one number for an isotropic material, or the tensor, with the
     ! dimension of the mesh's space
     real(dp) :: scalar = 0.0_dp
     real(dp), allocatable :: k(:,:)

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

    this % scalar = k

  end function isotropic

  pure type(conduction) function tensor(k) result(this)

    real(dp), intent(in) :: k(:,:)

    if (size(k, 1) /= size(k, 2)) error stop 'conduction: the conductivity tensor is square'
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
    real(dp), allocatable :: normals(:), n(:)
    integer :: ne, e, d

    fn = m % face_normal()
    call fn % real_vector(normals)

    ne = m % num_edges()
    d  = m % dimension
    allocate(values(ne), n(d))

    if (allocated(this % k)) then
       if (size(this % k, 1) /= d) then
          error stop 'conduction: the conductivity tensor has the dimension of the space'
       end if
    end if

    do e = 1, ne
       n = normals(d * (e - 1) + 1 : d * e)
       if (allocated(this % k)) then
          values(e) = dot_product(n, matmul(this % k, n))
       else
          values(e) = this % scalar * dot_product(n, n)
       end if
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
    call fa % real_vector(areas)

    do e = 1, size(values)
       if (m % edge_has_head(e)) then
          values(e) = values(e) * areas(e)
       else
          values(e) = 0.0_dp
       end if
    end do

  end subroutine edge_coefficients

end module operation_conduction
