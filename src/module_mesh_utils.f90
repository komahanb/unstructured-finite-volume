module module_mesh_utils

  use iso_fortran_env , only : dp => REAL64, error_unit
  use class_set             , only : set
  use class_list            , only : list

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

  pure real(dp) function distanceX(X)

    real(dp), intent(in)  :: X(:,:) ! [[x,y,z], [1:2]]

    distanceX = sqrt(sum((X(:,1)-X(:,2))**2))

  end function distanceX

  pure real(dp) function distanceAB(x, y)

    real(dp), intent(in)  :: X(:), y(:) ! [[x,y,z], [1:2]]

    distanceAB = sqrt(sum((x-y)**2))

  end function distanceAB

  !===================================================================!
  ! Find the intersection of two arrays (move elsewhere)?
  !===================================================================!

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

  subroutine transpose_connectivities(cell_vertices, num_cell_vertices, &
       & vertex_cells, num_vertex_cells)

    ! Arguments
    integer, intent(in) :: cell_vertices(:,:)
    integer, intent(in) :: num_cell_vertices(:)

    integer, allocatable, intent(out) :: vertex_cells(:,:)
    integer, allocatable, intent(out) :: num_vertex_cells(:)

    ! Locals
    integer, allocatable :: A(:,:)
    integer, allocatable :: cell_indices(:)
    integer :: num_cells, num_vertices, icell, ivertex

    num_cells    = size(num_cell_vertices, dim=1)
    num_vertices = maxval(cell_vertices) ! assume continuity starting from 1 ... num_vertices
    if (num_images().gt.1) error stop "only serial"

    ! form the matrix from forward connectivities
    allocate(A(num_vertices, num_cells))
    A = 0
    do icell = 1, num_cells
       do ivertex = 1, num_cell_vertices(icell)
          A(cell_vertices(ivertex, icell), icell) = 1
       end do
    end do

    ! count the number of nonzeros referring to cells
    allocate(num_vertex_cells(num_vertices))
    num_vertex_cells = 0
    do ivertex = 1, num_vertices
       num_vertex_cells(ivertex) = count(A(ivertex,:) .ne. 0)
    end do

    allocate(cell_indices(num_cells))
    forall(icell=1:num_cells) cell_indices(icell) = icell

    ! Fill the vertex cells based on A
    allocate(vertex_cells(maxval(num_vertex_cells), num_vertices))
    vertex_cells = 0
    do ivertex = 1, num_vertices
       vertex_cells(1:num_vertex_cells(ivertex),ivertex) = pack(cell_indices, A(ivertex,:) .ne. 0)
    end do

    if (allocated(A)) deallocate(A)
    if (allocated(cell_indices)) deallocate(cell_indices)

    if (sum(num_cell_vertices) .ne. sum(num_vertex_cells)) then
       error stop "inconsistent connectivity inversion"
    end if

  end subroutine transpose_connectivities

  subroutine sparse_transpose_matmul( &
       & num_cell_vertices, cell_types, cell_numbers, cell_tags, cell_vertices, &
       & num_face_vertices, face_types, face_numbers, face_tags, face_vertices, &
       & num_cell_faces   , cell_faces, cell_faces_type, &
       & num_face_cells   , face_cells, face_cells_type )

    ! arguments
    integer, intent(in)               :: num_cell_vertices(:)
    integer, intent(in)               :: cell_types(:)
    integer, intent(in)               :: cell_numbers(:)
    integer, intent(in)               :: cell_tags(:)
    integer, intent(in)               :: cell_vertices(:,:)
    integer, intent(in)               :: num_face_vertices(:)
    integer, intent(in)               :: face_types(:)
    integer, intent(in)               :: face_numbers(:)
    integer, intent(in)               :: face_tags(:)
    integer, intent(in)               :: face_vertices(:,:)

    integer, allocatable, intent(out) :: num_cell_faces(:)
    integer, allocatable, intent(out) :: cell_faces(:,:)
    integer, allocatable, intent(out) :: cell_faces_type(:,:)
    integer, allocatable, intent(out) :: num_face_cells(:)
    integer, allocatable, intent(out) :: face_cells(:,:)
    integer, allocatable, intent(out) :: face_cells_type(:,:)

    ! locals
    integer :: icell, iface, ivertex
    integer :: ncells, nfaces
    integer :: max_faces_in_cell, max_cells_in_face
    integer :: shared_vertex_count
    integer :: shared_vertices(maxval(elem_type_vertex_count(face_types)))
    integer :: inner_pdt_count

    max_faces_in_cell = maxval(elem_type_face_count(cell_types))
    max_cells_in_face = 2

    ncells = ubound(num_cell_vertices,1)
    nfaces = ubound(num_face_vertices,1)

    ! allocate return variables
    allocate(num_cell_faces(ncells))
    num_cell_faces = 0

    allocate(cell_faces(max_faces_in_cell, ncells))
    cell_faces = 0

    allocate(cell_faces_type(max_faces_in_cell, ncells))
    cell_faces_type = 0

    allocate(num_face_cells(nfaces))
    num_face_cells = 0

    allocate(face_cells(max_cells_in_face, nfaces))
    face_cells = 0

    allocate(face_cells_type(max_cells_in_face, nfaces))
    face_cells_type = 0

    inner_pdt_count = 0
    do concurrent (icell=1:ncells, iface=1:nfaces) ! keep track of how many faces of this cell is found

!       face_search: do iface = 1, nfaces

          ! filter out the boundary faces
          if (num_cell_faces(icell) .eq. elem_type_face_count(cell_types(icell))) then
             !            exit face_search
             cycle
          end if

          if ((face_tags(iface) .ne. cell_tags(icell)) &
               & .and. &
               & (num_face_cells(iface) .eq. 1)) then
             cycle
          end if

          if ((face_tags(iface) .eq. cell_tags(icell)) &
               & .and. &
               & (num_face_cells(iface) .eq. 2)) then
             cycle
          end if

          inner_pdt_count = inner_pdt_count + 1

          shared_vertex_count = 0

          shared_vertices = 0

          ! we can probably match faces if that's simpler
          ! each cell has only a given number of faces
          dot_pdt: do ivertex = 1, num_cell_vertices(icell)

             ! dot both functions with inner product in the space of
             ! vertices
             if (any (face_vertices(1:num_face_vertices(iface),iface) .eq. &
                  & cell_vertices(ivertex, icell) &
                  & ) .eqv. .true.) then

                shared_vertex_count = shared_vertex_count + 1

                shared_vertices(shared_vertex_count) = cell_vertices(ivertex, icell)

             end if

          end do dot_pdt

          select case(shared_vertex_count)
          case(3:4)
             ! triangle face or quadrangle
             num_cell_faces(icell)                         = num_cell_faces(icell) + 1
             cell_faces(num_cell_faces(icell), icell)      = face_numbers(iface) ! or iface
             cell_faces_type(num_cell_faces(icell), icell) = face_types(iface)

             num_face_cells(iface)                         = num_face_cells(iface) + 1
             face_cells(num_face_cells(iface), iface)      = cell_numbers(icell)
             face_cells_type(num_face_cells(iface), iface) = cell_types(icell)
          case(5:)
             print *, "more than 4 vertices do not make a face"
             error stop
          case default
             ! no vertex, one vertex, edge may be
          end select

  !     end do face_search

    end do

  end subroutine sparse_transpose_matmul

     !===================================================================!
  ! Invert a map. The keys of 'map' are values of 'inverse'. The
  ! values of 'map' are the keys of 'inverse'
  ! ===================================================================!

  subroutine reverse_map(map, num_map_vals, inverse, num_inverse_vals)

    ! Arguments
    integer, intent(in) :: map(:,:)
    integer, intent(in) :: num_map_vals(:)

    integer, allocatable, intent(out) :: inverse(:,:)
    integer, allocatable, intent(out) :: num_inverse_vals(:)

    ! Locals
    integer              :: value
    integer              :: i, j              ! loop indices
    integer              :: nkeysin, nkeysout ! input map size
    integer              :: nvalsin, nvalsout ! output map size
    integer, allocatable :: ptr(:)

    if (size(num_map_vals).eq.0) return

    ! Forward mapping size
    nkeysin = size(map, dim = 2)
    nvalsin = size(map, dim = 1)

    ! Nothing to do (probably empty map)!
    if (nkeysin.eq.0) return

    ! Find the maximum size of the inverse map values
    nkeysout = maxval(map) !dangerous
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

  !===================================================================!
  ! Forms the cell faces from a pair of vertices belonging to cell.
  !===================================================================!

  subroutine form_cell_faces( cell_vertices, &
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

  end subroutine form_cell_faces

  !===================================================================!
  ! Determine if the face is a boundary face based on how many
  ! neighbouring cells it has.
  !===================================================================!

  pure subroutine get_boundary_faces(num_face_cells, boundary_faces)

    integer, intent(in)               :: num_face_cells(:)
    integer, intent(out), allocatable :: boundary_faces(:)
    integer                           :: iface, nfaces, nbfaces, ctr

    nfaces = size(num_face_cells, dim=1)

    ! Boundary faces are the faces corresponding to just one cell
    nbfaces = 0
    do iface = 1, nfaces
       if (num_face_cells(iface) .eq. 1) then
          nbfaces = nbfaces + 1
       end if
    end do

    allocate(boundary_faces(nbfaces))
    ctr = 0
    do iface = 1, nfaces
       if (num_face_cells(iface) .eq. 1) then
          ctr = ctr + 1
          boundary_faces(ctr) = iface
       end if
    end do

  end subroutine get_boundary_faces

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

  pure type(logical) function same_entity(one, two)

    integer, intent(in) :: one(:), two(:)
    integer :: tmp(size(two))
    integer :: i, n

    n = size(two)
    tmp = two
    same_entity = .false.
    loop: do i = 0, n - 1
       same_entity = all(one .eq. cshift(tmp, i))
       if (same_entity .eqv. .true.) exit loop
    end do loop

  end function same_entity

  impure subroutine order_face_vertices(cell_type, cell_vertices, face_vertices_unordered)

    integer, intent(in)    :: cell_type
    integer, intent(in)    :: cell_vertices(:)
    integer, intent(inout) :: face_vertices_unordered(:)

    integer, allocatable :: face_vertices(:,:)
    integer :: num_face_vertices, num_cell_faces
    integer :: iface, ivertex
    integer :: match_count

    select case (cell_type)
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
