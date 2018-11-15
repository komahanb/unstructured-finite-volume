module module_mesh_utils

  use iso_fortran_env , only : dp => REAL64, error_unit
  
  implicit none

  interface distance
     module procedure distanceX
     module procedure distanceAB
  end interface distance

contains
  
  !===================================================================!
  ! Cross product for area computations
  ! ===================================================================!
  
  pure subroutine cross_product(a, b, pdt)

    real(dp), intent(in)  :: a(3), b(3)
    real(dp), intent(out) :: pdt(3)

    ! use skew form and generalize to n-dimensions?
    pdt(1) = a(2) * b(3) - a(3) * b(2)
    pdt(2) = a(3) * b(1) - a(1) * b(3)
    pdt(3) = a(1) * b(2) - a(2) * b(1)

  end subroutine cross_product

  !===================================================================!
  ! Geometric distance between two points
  !===================================================================!

  pure real(dp) function distanceX(X)

    real(dp), intent(in)  :: X(:,:) ! [[x,y,z], [1:2]]

    distanceX = sqrt(sum((X(:,1)-X(:,2))**2))

  end function distanceX

  pure real(dp) function distanceAB(x, y)

    real(dp), intent(in)  :: X(:), y(:) ! [[x,y,z], [1:2]]

    distanceAB = sqrt(sum((x-y)**2))

  end function distanceAB

  pure subroutine reverse_map(map, num_map_vals, inverse, num_inverse_vals)

    ! Arguments
    integer, intent(in) :: map(:,:)
    integer, intent(in) :: num_map_vals(:)

    integer, allocatable, intent(out) :: inverse(:,:)
    integer, allocatable, intent(out) :: num_inverse_vals(:)

    ! Locals
    integer              :: key, value
    integer              :: i, j              ! loop indices
    integer              :: nkeysin, nkeysout ! input map size
    integer              :: nvalsin, nvalsout ! output map size    
    integer, allocatable :: ptr(:)

    ! Forward mapping size 
    nkeysin = size(map, dim = 2)
    nvalsin = size(map, dim = 1)

    ! Find the maximum size of the inverse map values
    nkeysout = maxval(map)
    allocate(num_inverse_vals(nkeysout))
    num_inverse_vals = 0
    do i = 1, nkeysin
       do j = 1, num_map_vals(i)
          value = map(j,i)
          num_inverse_vals(value) = num_inverse_vals(value) + 1
       end do
    end do
    nvalsout = maxval(num_inverse_vals)

    ! Allocate the inverse map based on determined sizes
    allocate(inverse(nvalsout, nkeysout))
    inverse = 0

    ! Point into the next available slot for each inverse_key
    allocate(ptr(nkeysout))
    ptr = 0
    do i = 1, nkeysin
       do j = 1, num_map_vals(i)
          value = map(j,i)
          ptr(value) = ptr(value) + 1
          inverse(ptr(value), value) = i
       end do
    end do
    deallocate(ptr)

  end subroutine reverse_map

  pure subroutine intersection(a, b, c)

    ! Arguments
    integer, intent(in)  :: a(:)
    integer, intent(in)  :: b(:)
    integer, intent(out) :: c(:)

    ! Locals
    integer :: i, j
    integer :: sizea, sizeb
    integer :: ctr

    sizea = size(a)
    sizeb = size(b)

    ctr = 0
    do i = 1, sizea
       do j = 1, sizeb
          ! Copy entry to new list if equal
          if (a(i) .eq. b(j)) then
             ctr = ctr + 1
             c(ctr) = a(i)
             return
          end if
       end do
    end do

  end subroutine intersection

  pure type(integer) function find(array, target_value)

    integer, intent(in) :: array(:)
    integer, intent(in) :: target_value
    integer :: i, nentries

    nentries = size(array, dim=1)

    do i = 1, nentries
       if (array(i) .eq. target_value) then
          find = i
          return
       endif
    end do

    find = -1

  end function find

end module module_mesh_utils
