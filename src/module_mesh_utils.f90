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


  subroutine get_face_cells( &
       & cell_vertices, num_cell_vertices, &
       & face_vertices, num_face_vertices, &
       & face_cells   , num_face_cells)

    integer, intent(in)  :: cell_vertices(:,:)
    integer, intent(in)  :: num_cell_vertices(:)
    integer, intent(in)  :: face_vertices(:,:)
    integer, intent(in)  :: num_face_vertices(:)
    integer, allocatable, intent(out) :: face_cells(:,:)
    integer, allocatable, intent(out) :: num_face_cells(:)

    integer :: num_faces, iface
    integer :: num_cells, icell

    num_faces = size(face_vertices, dim = 2)
    num_cells = size(cell_vertices, dim = 2)

    allocate(face_cells(2, num_faces))
    face_cells = 0

    allocate(num_face_cells(num_faces))
    num_face_cells = 0

    loop_faces : do iface = num_faces, num_faces

       ! Get the face vertices
       associate( fvs => face_vertices(1:num_face_vertices(iface),iface) ) 
     
         print *, "face vertices:", fvs

         loop_cells: do icell = 1, num_cells

            ! Get cell vertices and face vertices and check if it is a
            ! subset
            print *, "cell vertices:", icell, cell_vertices(1:num_cell_vertices(icell),icell)

            if (is_subset(fvs, &
                 & cell_vertices(:,icell)) &
                 & .eqv. .true.) then

               print *, "found for face", iface, icell

               ! Found a cell
               num_face_cells(iface) = num_face_cells(iface) + 1

               face_cells(num_face_cells(iface), iface) = icell

            end if

            if (num_face_cells(iface) .eq. 2) exit loop_cells

         end do loop_cells

     end associate

  end do loop_faces

  end subroutine get_face_cells

  !===================================================================!
  ! Forms the cell faces from a pair of vertices belonging to cell.
  !===================================================================!
!!$  
!!$  subroutine get_face_cells( &
!!$       & cell_vertices, num_cell_vertices, &
!!$       & vertex_faces , num_vertex_faces, &
!!$       & face_cells   , num_face_cells)
!!$
!!$    integer, intent(in)  :: cell_vertices(:,:)
!!$    integer, intent(in)  :: num_cell_vertices(:)
!!$    integer, intent(in)  :: vertex_cells(:,:)
!!$    integer, intent(in)  :: num_vertex_cells(:)
!!$    integer, allocatable, intent(out) :: face_cells(:,:)
!!$    integer, allocatable, intent(out) :: num_face_cells(:)
!!$
!!$    integer :: num_faces, iface
!!$    integer :: num_cells, icell
!!$
!!$    num_faces = size(face_vertices, dim = 2)
!!$    num_cells = size(cell_vertices, dim = 2)
!!$
!!$    allocate(face_cells(2, num_faces))
!!$    face_cells = 0
!!$
!!$    allocate(num_face_cells(num_faces))
!!$    num_face_cells = 0
!!$
!!$    loop_faces : do iface = 1, num_faces
!!$
!!$       ! Get the face vertices
!!$       associate( fvs => face_vertices(1:num_face_vertices(iface),iface)) 
!!$     
!!$         print *, "face vertices:", fvs
!!$
!!$       loop_cells: do icell = 1, num_cells
!!$
!!$          ! Get cell vertices and face vertices and check if it is a
!!$          ! subset
!!$          print *, "cell vertices:", cell_vertices(1:num_cell_vertices(icell),icell)
!!$          if (is_subset(fvs, cell_vertices(1:num_cell_vertices(icell),icell)) .eqv. .true.) then
!!$
!!$             ! Found a cell
!!$             num_face_cells(iface) = num_face_cells(iface) + 1
!!$
!!$             face_cells(num_face_cells(iface), iface) = icell
!!$
!!$          end if
!!$
!!$          if (num_face_cells(iface) .eq. 2) exit loop_cells
!!$
!!$       end do loop_cells
!!$
!!$     end associate
!!$
!!$  end do loop_faces
!!$
!!$  end subroutine get_face_cells
  
  !===================================================================!
  ! Checks if the first argument is a subset of the second argument
  !===================================================================!
  
  pure type(logical) function is_subset(small, big)

    integer, intent(in) :: small(:)
    integer, intent(in) :: big(:)

    integer, allocatable :: sub(:)
    integer, allocatable :: set(:)

    integer :: lensub, i, j

    is_subset = .false.

    if (size(small) .gt. size(big)) return

    ! Create local copy of arrays
    allocate(sub, source = small)
    allocate(set, source = big)

    ! Sort two arrays
    call isort(sub)
    call isort(set)    
    lensub = size(sub)

    ! Check if all entries are equal upto the length of the smallest
    ! array
    is_subset = .true.
    do i = 1, lensub
       if (any(set .eq. sub(i)) .eqv. .false.) then
          is_subset = .false.
          exit
       end if
    end do

    deallocate(sub,set)    

  end function is_subset

  pure subroutine isort(array)

    integer, intent(inout) :: array(:)
    integer :: temp , j , k
    integer :: n

    n = ubound(array,1)

    do j = 1 , n
       do k = j + 1 , n
          if(array(j) > array(k)) then
             temp     = array(k)
             array(k) = array(j)
             array(j) = temp
          end if
       end do
    end do

  end subroutine isort

  !===================================================================!
  ! Forms the cell faces from a pair of vertices belonging to cell.
  !===================================================================!
  
  subroutine get_cell_faces( cell_vertices, &
       & vertex_faces, num_vertex_faces, &
       & cell_faces, num_cell_faces ) 
    
    integer, intent(in)  :: cell_vertices(:,:)
    integer, intent(in)  :: vertex_faces(:,:)
    integer, intent(in)  :: num_vertex_faces(:)

    integer, allocatable, intent(out) :: cell_faces(:,:)
    integer, allocatable, intent(out) :: num_cell_faces(:)

    integer :: icell, iface
    integer :: v1, v2
    integer :: face_ptr, ctr
    integer :: ncells, nvertices, nfaces

    nvertices = size(cell_vertices,dim=1)
    nfaces = nvertices
    ncells = size(cell_vertices,dim=2)

    ! find how many faces are there based on nodes
    allocate(num_cell_faces(ncells))
    do icell = 1, ncells      
       ctr = 0
       do iface = 1, nvertices
          if (cell_vertices(iface, icell) .ne. 0) then
             ctr = ctr + 1
          end if
       end do
       num_cell_faces(icell) = ctr
    end do

    ! Cell to face cell_vertices
    allocate(cell_faces(maxval(num_cell_faces),ncells))
    cell_faces = 0
    do icell = 1, ncells
       face_ptr = 0
       do iface = 1, num_cell_faces(icell)
          ! Get the first two vertices
          if (iface .eq. num_cell_faces(icell)) then
             v1 = cell_vertices(iface,icell)
             v2 = cell_vertices(1,icell)
          else 
             v1 = cell_vertices(iface,icell)
             v2 = cell_vertices(iface+1,icell)
          end if

          face_ptr = face_ptr + 1
          call intersection( &
               & vertex_faces(1:num_vertex_faces(v2), v2), &
               & vertex_faces(1:num_vertex_faces(v1), v1), &
               & cell_faces(face_ptr:face_ptr,icell))
       end do
    end do

  end subroutine get_cell_faces

end module module_mesh_utils
