!=====================================================================!
! The advection law: the flow's answer to a state.
!
! LEVEL 3 OF THE STRATIFICATION. An advection law holds one velocity
! and supplies two coefficient arrays, one entry per face, from the
! mesh's own normals and areas:
!
!      normal_speed        vn_e = v . n           the speed through
!                                                 each face
!      edge_coefficients   vn_e * area_e          the dictionary's
!                                                 advection
!                                                 coefficient, zero
!                                                 on the headless
!                                                 faces - a boundary
!                                                 face gets its
!                                                 advective closure
!                                                 from its condition
!
! The scheme is not the law's business. Handed to the calculus with
! one_sided true, the signed coefficient upwinds by itself - the
! sign of vn picks the end the flow leaves; handed with one_sided
! false, the term is the central average. The old assembler's
! weights are exactly these two cases:
!
!      upwind    wp = max(vn, 0)   wn = min(vn, 0)
!      central   wp = wn = vn / 2
!
! It owns no operator and no balance; it computes numbers.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_advection

  use iso_fortran_env  , only : dp => REAL64
  use class_graph_field, only : field
  use class_graph_mesh , only : mesh

  implicit none

  private
  public :: advection

  type :: advection

     real(dp) :: velocity(3) = 0.0_dp

   contains

     procedure :: normal_speed
     procedure :: edge_coefficients

  end type advection

  interface advection
     module procedure create
  end interface advection

contains

  pure type(advection) function create(velocity) result(this)

    real(dp), intent(in) :: velocity(3)

    this % velocity = velocity

  end function create

  !===================================================================!
  ! vn_e = v . n, one entry per face, from the mesh's normals.
  !===================================================================!

  subroutine normal_speed(this, m, values)

    class(advection), intent(in)       :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: values(:)

    type(field) :: fn
    real(dp), allocatable :: normals(:)
    integer :: ne, e

    fn = m % face_normal()
    call fn % get_real_vector(normals)

    ne = m % num_edges()
    allocate(values(ne))

    do e = 1, ne
       values(e) = dot_product(this % velocity, normals(3 * e - 2 : 3 * e))
    end do

  end subroutine normal_speed

  !===================================================================!
  ! The dictionary's advection coefficient: vn_e * area_e, zero on
  ! the headless faces.
  !===================================================================!

  subroutine edge_coefficients(this, m, values)

    class(advection), intent(in)       :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: values(:)

    type(field) :: fa
    real(dp), allocatable :: areas(:)
    integer :: e

    call this % normal_speed(m, values)

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

end module class_advection
