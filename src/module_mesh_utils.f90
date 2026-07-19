module module_mesh_utils

  use iso_fortran_env , only : dp => REAL64

  implicit none

  interface distance
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

  pure real(dp) function distanceAB(x, y)

    real(dp), intent(in)  :: X(:), y(:) ! [[x,y,z], [1:2]]

    distanceAB = sqrt(sum((x-y)**2))

  end function distanceAB

  !===================================================================!
  ! Index of a target value if present in the array
  !===================================================================!

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

  !===================================================================!
  ! Checks if the first smaller array is a subset of the second larger
  ! array
  !===================================================================!

  pure type(logical) function is_subset(small, big)

    integer, intent(in) :: small(:)
    integer, intent(in) :: big(:)

    integer, allocatable :: sub(:)
    integer, allocatable :: set(:)

    integer :: lensub, i

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

  !===================================================================!
  ! Sort an integer array ! move elsewhere?
  !===================================================================!

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

  ! generalize the name to return the number of lower dimensional entities
  pure elemental integer function elem_type_face_count(cell_type) result (num_faces)

    integer, intent(in) :: cell_type

    select case (cell_type)
    case (1)
       ! 2-node line.
       num_faces = 2 ! vertex
    case (2)
       ! 3-node triangle
       num_faces = 3  ! 3 edges
    case (3)
       ! 4-node quadrangle
       num_faces = 4 ! 4 edges
    case (4)
       ! 4-node tetrahedron
       num_faces = 4
    case (5)
       ! 8-node hexahedron
       num_faces = 6
    case (6)
       !  6-node prism
       num_faces = 5
    case(7)
       ! 5-node prism (pyramid)
       num_faces = 5
    case default
       num_faces = 0
    end select

  end function elem_type_face_count

  ! Spatial dimension of a gmsh element type (point=0, line=1,
  ! tri/quad=2, tet/hex/prism/pyramid=3). Used to classify elements into
  ! cells / faces / edges by dimension rather than by type.
  pure elemental integer function elem_type_dimension(elem_type) result (dim)

    integer, intent(in) :: elem_type

    select case (elem_type)
    case (15)
       ! 1-node point
       dim = 0
    case (1)
       ! 2-node line
       dim = 1
    case (2:3)
       ! triangle, quadrangle
       dim = 2
    case (4:7)
       ! tet, hex, prism, pyramid
       dim = 3
    case default
       dim = -1
    end select

  end function elem_type_dimension

  ! generalize the name to return the number of lower dimensional entities
  pure elemental integer function elem_type_vertex_count(elem_type) result (num_vertices)

    integer, intent(in) :: elem_type

    select case (elem_type)
    case (1)
       ! 2-node line.
       num_vertices = 2
    case (2)
       ! 3-node triangle
       num_vertices = 3
    case (3)
       ! 4-node quadrangle
       num_vertices = 4
    case (4)
       ! 4-node tetrahedron
       num_vertices = 4
    case (5)
       ! 8-node hexahedron
       num_vertices = 8
    case (6)
       !  6-node prism
       num_vertices = 6
    case(7)
       ! 5-node prism (pyramid)
       num_vertices = 5
    case default
       num_vertices = 0
    end select

  end function elem_type_vertex_count

  impure subroutine order_face_vertices(cell_type, cell_vertices, face_vertices_unordered)

    integer, intent(in)    :: cell_type
    integer, intent(in)    :: cell_vertices(:)
    integer, intent(inout) :: face_vertices_unordered(:)

    integer, allocatable :: face_vertices(:,:)
    integer :: num_face_vertices, num_cell_faces
    integer :: iface, ivertex
    integer :: match_count

    select case (cell_type)
    case (2)
       ! 3-node triangle (2d cell) - faces are its 3 edges
       num_face_vertices = 2
       num_cell_faces = 3

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3)]
       face_vertices(:,3) = [cell_vertices(3), cell_vertices(1)]

    case (3)
       ! 4-node quadrangle (2d cell) - faces are its 4 edges
       num_face_vertices = 2
       num_cell_faces = 4

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3)]
       face_vertices(:,3) = [cell_vertices(3), cell_vertices(4)]
       face_vertices(:,4) = [cell_vertices(4), cell_vertices(1)]

    case (4)
       ! 4-node tetrahedron
       num_face_vertices = 3
       num_cell_faces = 4

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       ! normals must point outward
       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2), cell_vertices(3)]
       face_vertices(:,2) = [cell_vertices(3), cell_vertices(2), cell_vertices(4)]
       face_vertices(:,3) = [cell_vertices(1), cell_vertices(3), cell_vertices(4)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(4), cell_vertices(2)]

    case (5)
       ! 8-node hexahedron
       num_face_vertices = 4
       num_cell_faces = 6

       if (num_face_vertices .ne. size(face_vertices_unordered)) &
            & error stop "inconstent vertices"

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(5), cell_vertices(8), cell_vertices(4)]
       face_vertices(:,2) = [cell_vertices(2), cell_vertices(3), cell_vertices(7), cell_vertices(6)]
       face_vertices(:,3) = [cell_vertices(5), cell_vertices(6), cell_vertices(7), cell_vertices(8)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(4), cell_vertices(3), cell_vertices(2)]
       face_vertices(:,5) = [cell_vertices(7), cell_vertices(3), cell_vertices(4), cell_vertices(8)]
       face_vertices(:,6) = [cell_vertices(1), cell_vertices(2), cell_vertices(6), cell_vertices(5)]

    case (6)
       !  6-node prism (5 faces)
       num_cell_faces = 5
       num_face_vertices = 4 ! take the max

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(4), cell_vertices(5), cell_vertices(6)]
       face_vertices(:,2) = [cell_vertices(1), cell_vertices(3), cell_vertices(6), cell_vertices(4)]
       face_vertices(:,3) = [cell_vertices(2), cell_vertices(5), cell_vertices(6), cell_vertices(3)]
       face_vertices(:,4) = [cell_vertices(1), cell_vertices(2), cell_vertices(3)]
       face_vertices(:,5) = [cell_vertices(4), cell_vertices(6), cell_vertices(5)]

    case(7)
       ! 5-node prism (pyramid) 5 faces
       num_cell_faces = 5
       num_face_vertices = 4 ! take the max

       allocate(face_vertices(num_face_vertices, num_cell_faces))
       face_vertices = 0

       face_vertices(:,1) = [cell_vertices(1), cell_vertices(2), cell_vertices(3), cell_vertices(4)]
       face_vertices(:,2) = [cell_vertices(1), cell_vertices(2), cell_vertices(5)]
       face_vertices(:,3) = [cell_vertices(2), cell_vertices(3), cell_vertices(5)]
       face_vertices(:,4) = [cell_vertices(3), cell_vertices(4), cell_vertices(5)]
       face_vertices(:,5) = [cell_vertices(1), cell_vertices(5), cell_vertices(4)]

    case default
       print *, "unknown type"
       error stop
    end select

    ! compare the input with face 1,2,3,4 and see which has the
    ! matching count equal to input vertex count

    match_face: do iface = 1, num_cell_faces

       match_count = 0
       match_vertex: do ivertex = 1, num_face_vertices
          if (any (face_vertices_unordered(:) .eq. face_vertices(ivertex,iface) ) .eqv. .true.) then
             ! exit match_vertex ! skip and go to next face
             ! else
             match_count = match_count + 1
          end if
       end do match_vertex

       ! found the face with correct ordering of vertices
       if (match_count .eq. num_face_vertices) then
          face_vertices_unordered = face_vertices(:, iface)
          exit match_face
       end if

    end do match_face

    if (allocated(face_vertices)) deallocate(face_vertices)

  end subroutine order_face_vertices

end module module_mesh_utils
