!=====================================================================!
! The mesh suite: rung 2's acceptance.
!
! A hand mesh of three cells in a row, two interior faces and one
! wall face, every measurement chosen small and unequal so nothing
! cancels by accident:
!
!        wall                          e1              e2
!      o------ (1)              (1) --------> (2) --------> (3)
!        e3                     vol  2         3            4
!
!      areas   [1.5, 2.5, 1.0]      deltas  [0.5, 0.25, 1.0]
!      k       [2.0, 4.0,  - ]      q       [1, 3, 6]
!
! Three things are proven here. The mesh IS a graph: the inherited
! structure answers. The measurements come back through compiled
! names, exactly as given. And the dictionary of
! geometry-to-operator-mapping.md holds: the order-2 operator built
! from the mesh's own fields reproduces the old assembler's two-point
! rows entry for entry,
!
!      y_p = - sum over faces of  farea * keff/fdelta * (q_p - q_n)
!
! first with unit measures (the raw row), then seated on the cell
! volumes (the row over vol_p).
!=====================================================================!

program test_graph_mesh

  use iso_fortran_env, only : dp => REAL64
  use structure_graph, only : graph
  use data_graph_field, only : graph_field
  use structure_graph, only : GRAPH_SIDE_VERTEX, GRAPH_SIDE_EDGE
  use structure_support, only : support
  use data_field  , only : field
  use structure_mesh   , only : mesh
  use class_mesh_builder , only : mesh_from_gmsh
  use operation_differential, only : differential
  use operation_differential, only : vertex_differential_operator

  implicit none

  integer :: nfail

  nfail = 0

  call check_mesh_is_a_graph(nfail)
  call check_measurements_return(nfail)
  call check_dictionary_rows(nfail)
  call check_gmsh_path(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all mesh checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' mesh checks FAILED'
     error stop 1
  end if

contains

  !===================================================================!
  ! One line per check, one counter for the verdict.
  !===================================================================!

  subroutine report(ok, message, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! The hand mesh, one place, so every check reads the same numbers.
  !===================================================================!

  type(mesh) function hand_mesh() result(m)

    m = mesh(3, tails=[1, 2, 1], heads=[2, 3, 0], &
         & volumes      = [2.0_dp, 3.0_dp, 4.0_dp], &
         & cell_centers = [0.5_dp, 0.0_dp, 0.0_dp, &
         &                 1.5_dp, 0.0_dp, 0.0_dp, &
         &                 2.5_dp, 0.0_dp, 0.0_dp], &
         & areas        = [1.5_dp, 2.5_dp, 1.0_dp], &
         & deltas       = [0.5_dp, 0.25_dp, 1.0_dp], &
         & normals      = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 1.0_dp, 0.0_dp, 0.0_dp, &
         &                -1.0_dp, 0.0_dp, 0.0_dp], &
         & face_centers = [1.0_dp, 0.0_dp, 0.0_dp, &
         &                 2.0_dp, 0.0_dp, 0.0_dp, &
         &                 0.0_dp, 0.0_dp, 0.0_dp], &
         & weights      = [0.5_dp, 0.5_dp, 1.0_dp], &
         & etags        = [character(len=4) :: '', '', 'wall'])

  end function hand_mesh

  !===================================================================!
  ! The mesh answers as the graph it is.
  !===================================================================!

  subroutine check_mesh_is_a_graph(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    class(graph), allocatable :: wall
    integer, allocatable :: nbrs(:)

    m = hand_mesh()

    call report(m % num_vertices() == 3, 'three cells stand as vertices', nfail)
    call report(m % num_edges() == 3, 'three faces stand as edges', nfail)
    call report(.not. m % edge_has_head(3), 'the wall face has no head', nfail)

    call m % adjacent_vertices(2, nbrs)
    call report(size(nbrs) == 2, 'the middle cell has two neighbours', nfail)

    call m % tagged_edges('wall', wall)
    call report(wall % num_vertices() == 1, &
         & 'the wall tag names one face', nfail)
    call report(wall % global_vertex_index(1) == 3, &
         & 'and it is the third face', nfail)

  end subroutine check_mesh_is_a_graph

  !===================================================================!
  ! The measurements come back through compiled names, as given.
  !===================================================================!

  subroutine check_measurements_return(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(field) :: f
    real(dp), allocatable :: v(:)

    m = hand_mesh()

    f = m % cell_volume()
    call f % get_real_vector(v)
    call report(size(v) == 3 .and. all(abs(v - [2.0_dp, 3.0_dp, 4.0_dp]) < 1.0d-14), &
         & 'cell_volume returns one volume per cell, as given', nfail)

    f = m % face_area()
    call f % get_real_vector(v)
    call report(all(abs(v - [1.5_dp, 2.5_dp, 1.0_dp]) < 1.0d-14), &
         & 'face_area returns one area per face, as given', nfail)

    f = m % face_delta()
    call f % get_real_vector(v)
    call report(all(abs(v - [0.5_dp, 0.25_dp, 1.0_dp]) < 1.0d-14), &
         & 'face_delta returns the spacings, as given', nfail)

    f = m % face_normal()
    call f % get_real_vector(v)
    call report(size(v) == 9 .and. abs(v(7) + 1.0_dp) < 1.0d-14, &
         & 'face_normal is three wide and keeps its signs', nfail)

    f = m % cell_center()
    call report(f % num_components() == 3, &
         & 'cell_center is three wide', nfail)

    f = m % face_weights()
    call f % get_real_vector(v)
    call report(all(abs(v - [0.5_dp, 0.5_dp, 1.0_dp]) < 1.0d-14), &
         & 'face_weights returns the interpolation weights', nfail)

  end subroutine check_measurements_return

  !===================================================================!
  ! THE DICTIONARY CHECK. Build the order-2 operator from the mesh's
  ! own measurements - c_e = keff*farea, h_e = fdelta - and compare
  ! its rows against the old assembler's two-point formula, computed
  ! by hand:
  !
  !      y_1 = -1.5*2/0.5  *(1-3)            =   12
  !      y_2 = -12 - 2.5*4/0.25*(3-6)        =  108
  !      y_3 = -2.5*4/0.25*(6-3)             = -120
  !
  ! The wall face carries no conductivity here, so it contributes
  ! nothing; its condition is rung 3's business. Then the same rows
  ! seated on the cell volumes: y over [2, 3, 4].
  !===================================================================!

  subroutine check_dictionary_rows(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(differential) :: second
    type(field) :: state
    type(support) :: cells
    class(graph_field), allocatable :: y
    real(dp), allocatable :: farea(:), fdelta(:), vol(:), got(:)
    real(dp) :: k(3), c(3)
    type(field) :: fa, fd, fv
    integer :: v

    m = hand_mesh()

    ! The mesh's own measurements feed the constructor, per the
    ! dictionary. The wall face has no conductivity.
    fa = m % face_area()
    call fa % get_real_vector(farea)
    fd = m % face_delta()
    call fd % get_real_vector(fdelta)
    fv = m % cell_volume()
    call fv % get_real_vector(vol)

    k = [2.0_dp, 4.0_dp, 0.0_dp]
    c = k * farea

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, 3)])
    state = field('q', cells)
    call state % set_real_vector([1.0_dp, 3.0_dp, 6.0_dp])

    ! The raw rows: unit measures.
    second = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=fdelta)
    call second % apply(m, [state], y)
    call y % get_real_vector(got)

    call report(size(got) == 3, 'the operator answers one row per cell', nfail)
    call report(all(abs(got - [12.0_dp, 108.0_dp, -120.0_dp]) < 1.0d-11), &
         & 'the rows match the old two-point stencil entry for entry', nfail)

    ! The same rows seated on the cell volumes.
    second = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=fdelta, measures=vol)
    call second % apply(m, [state], y)
    call y % get_real_vector(got)

    call report(all(abs(got - [6.0_dp, 36.0_dp, -30.0_dp]) < 1.0d-11), &
         & 'seated on the volumes, each row divides by its cell', nfail)

  end subroutine check_dictionary_rows

  !===================================================================!
  ! THE GMSH PATH. A real file walks the whole road: parsed and
  ! measured by the bridge, seated on the tower, then proven by the
  ! dictionary's conservation law - with conductivity on the interior
  ! faces only, the interior rows exchange and cancel, so the rows
  ! sum to zero over the whole mesh to machine precision.
  !===================================================================!

  subroutine check_gmsh_path(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(differential) :: second
    type(field) :: state
    type(support) :: cells
    class(graph_field), allocatable :: y
    type(field) :: fa, fd, fv
    real(dp), allocatable :: farea(:), fdelta(:), vol(:), got(:), c(:), q(:)
    integer :: nv, ne, e, v, nwall

    m = mesh_from_gmsh('../square-10.msh')

    nv = m % num_vertices()
    ne = m % num_edges()

    call report(nv > 0 .and. ne > nv, &
         & 'the file arrives: more faces than cells, as a mesh has', nfail)

    fv = m % cell_volume()
    call fv % get_real_vector(vol)
    call report(all(vol > 0.0_dp), 'every cell volume is positive', nfail)

    fa = m % face_area()
    call fa % get_real_vector(farea)
    call report(all(farea > 0.0_dp), 'every face area is positive', nfail)

    fd = m % face_delta()
    call fd % get_real_vector(fdelta)
    call report(all(fdelta > 0.0_dp), 'every face delta is positive', nfail)

    nwall = 0
    do e = 1, ne
       if (.not. m % edge_has_head(e)) nwall = nwall + 1
    end do
    call report(nwall > 0, 'the boundary faces arrived without heads', nfail)

    ! Conservation: conductivity on interior faces only, an uneven
    ! state, and the rows must sum to zero - what every face gives
    ! one cell it takes from the other.
    allocate(c(ne))
    do e = 1, ne
       c(e) = farea(e)
       if (.not. m % edge_has_head(e)) c(e) = 0.0_dp
    end do

    allocate(q(nv))
    do v = 1, nv
       q(v) = real(mod(7 * v, 11), dp)
    end do

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    state = field('q', cells)
    call state % set_real_vector(q)

    second = vertex_differential_operator(order=2, coefficients=c, &
         & spacings=fdelta)
    call second % apply(m, [state], y)
    call y % get_real_vector(got)

    call report(abs(sum(got)) < 1.0d-9 * sum(abs(got)) + 1.0d-12, &
         & 'the interior rows conserve: they sum to zero on the real mesh', nfail)

  end subroutine check_gmsh_path

end program test_graph_mesh
