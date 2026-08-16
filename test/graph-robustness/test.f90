!=====================================================================!
! The robustness suite: numerics that do not ask the mesh's pardon.
!
! A deliberately bad mesh: four cells whose centers wander off the
! face normals, every face skewed, nothing orthogonal anywhere,
!
!        (3) ~~~~ (4)          c1 = (0.00, 0.00)
!         |        |           c2 = (1.00, 0.45)
!        (1) ~~~~ (2)          c3 = (0.15, 1.10)
!                              c4 = (1.30, 1.60)
!
! and one exact linear field, q = 2x + 3y. Three verdicts:
!
!   1. the two-point flux MISSES, and by how much is measured - it
!      reads the derivative along the center line, not the normal
!   2. the exactness weights land the flux to machine precision on
!      the same mesh, because they use the geometry instead of
!      assuming it
!   3. the solver on the exact stencil recovers the exact field -
!      the numerics owe nothing to the mesh's quality
!=====================================================================!

program test_graph_robustness

  use class_graph, only : stored_graph
  use iso_fortran_env, only : dp => REAL64
  use graph_ordinary_view, only : graph
  use graph_field_calculus, only : graph_field
  use graph_calculus , only : GRAPH_SIDE_VERTEX
  use fractal_graph  , only : set_graph => graph
  use class_graph_field  , only : field
  use class_graph_mesh   , only : mesh
  use class_graph_differential_operator, only : edge_differential_operator
  use class_graph_differential_operator, only : differential_operator
  use class_graph_stencil, only : stencil_operator
  use class_fitted_balance, only : fitted_balance_stencil
  use graph_fitting        , only : fit
  use class_polynomial_form, only : polynomial_form
  use class_harmonic_form  , only : harmonic_form
  use class_form_pruner    , only : pruner
  use class_graph_gmres  , only : gmres

  implicit none

  integer :: nfail

  nfail = 0

  call check_derived_stencils(nfail)
  call check_the_form_is_free(nfail)
  call check_two_point_misses(nfail)
  call check_exactness_lands(nfail)
  call check_solver_owes_nothing(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all robustness checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' robustness checks FAILED'
     error stop 1
  end if

contains

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
  ! The bad mesh, in one place. Interior deltas are the true
  ! center-to-center distances; every one of them leans away from
  ! its normal.
  !===================================================================!

  type(mesh) function skewed_mesh() result(m)

    real(dp), parameter :: c1(2) = [0.00_dp, 0.00_dp]
    real(dp), parameter :: c2(2) = [1.00_dp, 0.45_dp]
    real(dp), parameter :: c3(2) = [0.15_dp, 1.10_dp]
    real(dp), parameter :: c4(2) = [1.30_dp, 1.60_dp]

    integer :: kf

    m = mesh(4, &
         & tails=[1, 3, 1, 2,  1, 3, 1, 2, 2, 4, 3, 4], &
         & heads=[2, 4, 3, 4,  0, 0, 0, 0, 0, 0, 0, 0], &
         & volumes      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
         & cell_centers = [c1(1), c1(2), 0.0_dp, &
         &                 c2(1), c2(2), 0.0_dp, &
         &                 c3(1), c3(2), 0.0_dp, &
         &                 c4(1), c4(2), 0.0_dp], &
         & areas        = [(1.0_dp, kf = 1, 12)], &
         & deltas       = [norm2(c2 - c1), norm2(c4 - c3), &
         &                 norm2(c3 - c1), norm2(c4 - c2), &
         &                 (0.5_dp, kf = 1, 8)], &
         & normals      = [ 1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                  0.0_dp, 1.0_dp, 0.0_dp, &
         &                  0.0_dp, 1.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                 -1.0_dp, 0.0_dp, 0.0_dp, &
         &                  0.0_dp,-1.0_dp, 0.0_dp, &
         &                  0.0_dp,-1.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                  1.0_dp, 0.0_dp, 0.0_dp, &
         &                  0.0_dp, 1.0_dp, 0.0_dp, &
         &                  0.0_dp, 1.0_dp, 0.0_dp], &
         & face_centers = [0.5_dp, 0.2_dp, 0.0_dp, &
         &                 0.7_dp, 1.4_dp, 0.0_dp, &
         &                 0.1_dp, 0.6_dp, 0.0_dp, &
         &                 1.2_dp, 1.0_dp, 0.0_dp, &
         &                -0.5_dp, 0.1_dp, 0.0_dp, &
         &                -0.4_dp, 1.2_dp, 0.0_dp, &
         &                 0.1_dp,-0.5_dp, 0.0_dp, &
         &                 1.1_dp,-0.2_dp, 0.0_dp, &
         &                 1.6_dp, 0.5_dp, 0.0_dp, &
         &                 1.9_dp, 1.7_dp, 0.0_dp, &
         &                 0.2_dp, 1.7_dp, 0.0_dp, &
         &                 1.5_dp, 2.2_dp, 0.0_dp], &
         & weights      = [(0.5_dp, kf = 1, 12)])

  end function skewed_mesh

  !===================================================================!
  ! The exact field and its values wherever a point stands.
  !===================================================================!

  pure real(dp) function exact_at(x, y)

    real(dp), intent(in) :: x, y

    exact_at = 2.0_dp * x + 3.0_dp * y

  end function exact_at

  !===================================================================!
  ! VERDICT ZERO. The classical kernels are theorems of the
  ! exactness operation: two collinear points at degree one force
  ! the two-point stencil -+1/delta, derived rather than assumed.
  !===================================================================!

  subroutine check_derived_stencils(nfail)

    integer, intent(inout) :: nfail

    type(fit) :: fitting
    type(stored_graph) :: pair
    type(field)   :: positions
    class(graph_field), allocatable :: answer
    real(dp), allocatable :: w(:)

    pair = stored_graph(2, tails=[integer ::], heads=[integer ::])
    positions = field('positions', pair % vertex_set(), pair % num_vertices(), ncomp=3)
    call positions % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp, &
         &                            0.5_dp, 0.0_dp, 0.0_dp])

    fitting = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call fitting % apply(pair, [positions], answer)
    call answer % get_real_vector(w)

    call report(all(abs(w - [-2.0_dp, 2.0_dp]) < 1.0d-12), &
         & 'two collinear points force the two-point stencil: a theorem', nfail)

    fitting = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp], scale=3.0_dp)
    call fitting % apply(pair, [positions], answer)
    call answer % get_real_vector(w)

    call report(all(abs(w - [-6.0_dp, 6.0_dp]) < 1.0d-12), &
         & 'and the conductivity rides the scale, as the dictionary says', nfail)

  end subroutine check_derived_stencils

  !===================================================================!
  ! VERDICT ZERO-B. The basis is truly the free seat, and the form
  ! is governed. A harmonic fit differentiates its own waves
  ! exactly, where polynomials only approximate them; and the
  ! pruner strikes the members two collinear points cannot see,
  ! after which the polynomial fit still lands its theorem.
  !===================================================================!

  subroutine check_the_form_is_free(nfail)

    integer, intent(inout) :: nfail

    type(fit)    :: wave
    type(fit)    :: poly
    type(pruner) :: gardener
    type(stored_graph) :: trio, pair
    type(field)   :: positions
    class(graph_field), allocatable :: answer
    real(dp), allocatable :: w(:)
    real(dp) :: pts(9), sampled(3), got, expected
    integer :: j
    integer, allocatable :: kept(:)

    ! Three points on a line, one wave through them.
    pts = [0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.4_dp, 0.0_dp, 0.0_dp, &
         & 0.8_dp, 0.0_dp, 0.0_dp]

    trio = stored_graph(3, tails=[integer ::], heads=[integer ::])
    positions = field('positions', trio % vertex_set(), trio % num_vertices(), ncomp=3)
    call positions % set_real_vector(pts)

    wave = fit(harmonic_form([2.5_dp, 0.0_dp, 0.0_dp]), &
         & at=[0.4_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call wave % apply(trio, [positions], answer)
    call answer % get_real_vector(w)

    do j = 1, 3
       sampled(j) = sin(2.5_dp * pts(3 * j - 2))
    end do

    got      = sum(w * sampled)
    expected = 2.5_dp * cos(2.5_dp * 0.4_dp)

    call report(abs(got - expected) < 1.0d-10, &
         & 'the harmonic fit differentiates its own wave exactly', nfail)

    ! Two collinear points: the pruner strikes what they cannot see,
    ! and the fit still lands the two-point theorem.
    pair = stored_graph(2, tails=[integer ::], heads=[integer ::])
    positions = field('positions', pair % vertex_set(), pair % num_vertices(), ncomp=3)
    call positions % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp, &
         &                            0.5_dp, 0.0_dp, 0.0_dp])

    poly = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call gardener % adapt(poly % shape, [0.0_dp, 0.0_dp, 0.0_dp, &
         &                              0.5_dp, 0.0_dp, 0.0_dp])

    call poly % shape % members(kept)
    call report(size(kept) == 2 .and. all(kept == [1, 2]), &
         & 'the pruner strikes the members the points cannot see', nfail)

    call poly % apply(pair, [positions], answer)
    call answer % get_real_vector(w)
    call report(all(abs(w - [-2.0_dp, 2.0_dp]) < 1.0d-12), &
         & 'and the governed fit still lands the theorem', nfail)

  end subroutine check_the_form_is_free

  !===================================================================!
  ! VERDICT ONE. The two-point flux on the skewed mesh, measured
  ! against the exact normal derivative it claims to be.
  !===================================================================!

  subroutine check_two_point_misses(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(differential_operator) :: slope
    type(field) :: state, fd
    type(set_graph) :: cells
    class(graph_field), allocatable :: z
    real(dp), allocatable :: got(:), deltas(:), q(:)
    real(dp) :: exact_flux(4), worst
    integer :: v, e

    m = skewed_mesh()

    cells = m % vertex_set()
    state = field('q', cells, m % num_vertices())
    q = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]
    call state % set_real_vector(q)

    fd = m % face_delta()
    call fd % get_real_vector(deltas)

    slope = edge_differential_operator(order=1, spacings=deltas)
    call slope % apply(m, [state], z)
    call z % get_real_vector(got)

    ! What the flux claims to be: keff*area*dq/dn, here dq/dn along
    ! each interior normal.
    exact_flux = [2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp]

    worst = 0.0_dp
    do e = 1, 4
       worst = max(worst, abs(got(e) - exact_flux(e)) / abs(exact_flux(e)))
    end do

    write(*, '(a, f6.1, a)') '        measured: the two-point flux is off by up to ', &
         & 100.0_dp * worst, ' percent on this mesh'

    call report(worst > 0.25_dp, &
         & 'the two-point flux misses on the skewed mesh, and it is measured', nfail)

  end subroutine check_two_point_misses

  !===================================================================!
  ! VERDICT TWO. The exactness weights on the same mesh, same field:
  ! every row must balance to machine precision, because the flux is
  ! exact for linears by construction and the scatter conserves.
  !===================================================================!

  subroutine check_exactness_lands(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(stencil_operator) :: op
    type(field) :: state
    type(set_graph) :: cells
    class(graph_field), allocatable :: y
    real(dp), allocatable :: got(:), q(:), vb(:), centers(:)
    type(field) :: fc
    integer :: v, e

    m = skewed_mesh()

    ! The wall values: the exact field, read at each wall's center.
    fc = m % face_center()
    call fc % get_real_vector(centers)

    allocate(vb(12))
    vb = 0.0_dp
    do e = 5, 12
       vb(e) = exact_at(centers(3 * e - 2), centers(3 * e - 1))
    end do

    block
      type(field) :: fa
      real(dp), allocatable :: farea(:)
      fa = m % face_area()
      call fa % get_real_vector(farea)
      op = fitted_balance_stencil(m, polynomial_form(), farea, &
           & boundary_values=vb)
    end block

    cells = m % vertex_set()
    state = field('q', cells, m % num_vertices())
    q = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]
    call state % set_real_vector(q)

    call op % apply(m, [state], y)
    call y % get_real_vector(got)

    call report(maxval(abs(got)) < 1.0d-10, &
         & 'the exactness weights balance the same field to machine zero', nfail)

  end subroutine check_exactness_lands

  !===================================================================!
  ! VERDICT THREE. Solve on the bad mesh with the exact stencil: the
  ! recovered field is the exact field. The solver owes nothing to
  ! the mesh.
  !===================================================================!

  subroutine check_solver_owes_nothing(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(stencil_operator) :: op
    type(gmres) :: gm
    type(field) :: fc
    real(dp), allocatable :: vb(:), centers(:), g(:), rhs(:), x(:)
    real(dp) :: achieved
    real(dp) :: q_exact(4)
    integer :: v, e

    m = skewed_mesh()

    fc = m % face_center()
    call fc % get_real_vector(centers)

    allocate(vb(12))
    vb = 0.0_dp
    do e = 5, 12
       vb(e) = exact_at(centers(3 * e - 2), centers(3 * e - 1))
    end do

    block
      type(field) :: fa
      real(dp), allocatable :: farea(:)
      fa = m % face_area()
      call fa % get_real_vector(farea)
      op = fitted_balance_stencil(m, polynomial_form(), farea, &
           & boundary_values=vb)
    end block

    call gm % attach(op, m, m % vertex_set(), m % num_vertices())
    gm % tolerance = 1.0d-12

    call gm % constant(g)
    rhs = -g

    allocate(x(4))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)

    q_exact = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]

    call report(achieved < 1.0d-10, &
         & 'the solve closes on the bad mesh', nfail)
    call report(all(abs(x - q_exact) < 1.0d-8), &
         & 'and recovers the exact field: the numerics owe the mesh nothing', nfail)

  end subroutine check_solver_owes_nothing

end program test_graph_robustness
