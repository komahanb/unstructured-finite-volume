!=====================================================================!
! Module that states the face coefficient of a material or a
! velocity: a form in the face normal.
!
! This is level 3 of the stratification. One law stores one form in
! the unit normal n of a face, of degree zero, one or two. The rank of
! the stored coefficient is the number of times n enters, n.n = 1 for
! the unit normal contracting it away at degree zero:
!
!             n
!             ^
!     cell    |    cell        conduction(k)    k          isotropic, K = k * I
!   ----------+----------      advection(v)     v.n        the speed through the face
!            face              conduction(K)    n^T K n    the tensor material
!
! It supplies two coefficient arrays, one entry per face, from the
! mesh's own normals and areas:
!
!      normal_value        the form at the face normal
!      edge_coefficients   the form times the face area, zero on the
!                          headless faces when requested: a boundary
!                          face takes its closure from its condition
!
! The scheme is not part of the law. An advection coefficient passed
! to the calculus with one_sided true upwinds, the sign of v.n
! selecting the upstream end; with one_sided false the term is the
! central average:
!
!      upwind    wp = max(vn, 0)   wn = min(vn, 0)
!      central   wp = wn = vn / 2
!
! The law stores no operator and no balance; it computes coefficients.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_conduction

  use util_precision, only : dp
  use view_mesh     , only : mesh, values_of

  implicit none

  private

  public :: conduction, advection

  !-------------------------------------------------------------------!
  ! Type that models a face coefficient: a number, a vector or a
  ! tensor, contracted with the face normal
  !-------------------------------------------------------------------!

  type :: conduction

     ! Exactly one of the three is the form: the number for an
     ! isotropic material, the vector for a velocity, and otherwise the
     ! tensor, of the dimension of the mesh's space
     real(dp)              :: scalar = 0.0_dp
     real(dp), allocatable :: vector(:)
     real(dp), allocatable :: tensor(:,:)

   contains

     procedure :: normal_value
     procedure :: edge_coefficients

  end type conduction

  interface conduction
     module procedure isotropic, tensor
  end interface conduction

contains

  !===================================================================!
  ! Creates the law of an isotropic material, one number
  !===================================================================!

  pure type(conduction) function isotropic(k) result(this)

    real(dp), intent(in) :: k

    this % scalar = k

  end function isotropic

  !===================================================================!
  ! Creates the law of a velocity
  !===================================================================!

  pure type(conduction) function advection(velocity) result(this)

    real(dp), intent(in) :: velocity(:)

    this % vector = velocity

  end function advection

  !===================================================================!
  ! Creates the law of a material by its full tensor, which must be
  ! square
  !===================================================================!

  pure type(conduction) function tensor(k) result(this)

    real(dp), intent(in) :: k(:,:)

    if (size(k, 1) /= size(k, 2)) error stop 'conduction: k must be square, but size(k,1) /= size(k,2)'
    this % tensor = k

  end function tensor

  !===================================================================!
  ! Returns the form at the normal, one entry per face. Invalid input:
  ! a vector or a tensor whose extent is not the mesh's dimension.
  !===================================================================!

  subroutine normal_value(this, m, values)

    class(conduction)    , intent(in)  :: this
    type(mesh)           , intent(in)  :: m
    real(dp), allocatable, intent(out) :: values(:)

    real(dp), allocatable :: normals(:), n(:)
    integer               :: ne, e, d
    character(len=250)    :: message

    call values_of(m % face_normal(), normals)

    ne = m % num_edges()
    d  = m % dimension
    allocate(values(ne), n(d))

    if (allocated(this % vector)) then
       if (size(this % vector) /= d) then
          write(message,'(a,i0,a,i0)') 'advection: velocity must have one component per space &
               &dimension; size(velocity) = ', size(this % vector), ', dimension = ', d
          error stop trim(message)
       end if
    end if
    if (allocated(this % tensor)) then
       if (size(this % tensor, 1) /= d) then
          write(message,'(a,i0,a,i0)') 'conduction: the conductivity tensor must have the &
               &dimension of the space; size(tensor,1) = ', size(this % tensor, 1), &
               & ', dimension = ', d
          error stop trim(message)
       end if
    end if

    ! The normals stored face after face, d components each
    do e = 1, ne
       n = normals(d * (e - 1) + 1 : d * e)
       if (allocated(this % tensor)) then
          values(e) = dot_product(n, matmul(this % tensor, n))
       else if (allocated(this % vector)) then
          values(e) = dot_product(this % vector, n)
       else
          values(e) = this % scalar
       end if
    end do

  end subroutine normal_value

  !===================================================================!
  ! Returns the coefficient the calculus reads: the form times the face
  ! area, and zero on the headless faces when headless_zero is true
  !===================================================================!

  subroutine edge_coefficients(this, m, headless_zero, values)

    class(conduction)    , intent(in)  :: this
    type(mesh)           , intent(in)  :: m
    logical              , intent(in)  :: headless_zero
    real(dp), allocatable, intent(out) :: values(:)

    real(dp), allocatable :: areas(:)
    integer               :: e

    call this % normal_value(m, values)
    call values_of(m % face_area(), areas)

    do e = 1, size(values)
       if (m % edge_has_head(e) .or. .not. headless_zero) then
          values(e) = values(e) * areas(e)
       else
          values(e) = 0.0_dp
       end if
    end do

  end subroutine edge_coefficients

end module operation_conduction
