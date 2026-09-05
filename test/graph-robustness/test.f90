!=====================================================================!
! The robustness suite: numerics independent of mesh quality.
!
! A deliberately bad mesh: four cells whose centres wander off the
! face normals, every face skewed, nothing orthogonal anywhere,
!
!        (3) ~~~~ (4)          c1 = (0.00, 0.00)
!         |        |           c2 = (1.00, 0.45)
!        (1) ~~~~ (2)          c3 = (0.15, 1.10)
!                              c4 = (1.30, 1.60)
!
! and one exact linear field, q = 2x + 3y. Three checks:
!
!   1. the two-point flux is in error, and by how much is measured -
!      it reads the derivative along the centre line, not the normal
!   2. the exactness weights land the flux to machine precision on
!      the same mesh, because they use the geometry instead of
!      assuming it
!   3. the solver on the exact stencil recovers the exact field -
!      the numerics do not depend on the mesh's quality
!=====================================================================!

program test_graph_robustness

  use view_directed_stored, only : stored_directed_graph
  use iso_fortran_env, only : dp => REAL64
  use field_calculus, only : field
  use graph_fractal  , only : graph
  use field_stored  , only : stored_field
  use view_mesh   , only : mesh
  use operation_differential, only : gradient
  use operation_differential, only : differential_operator
  use operation_stencil, only : stencil
  use operation_fitted_balance, only : fitted_balance_stencil
  use operation_fitting        , only : fit
  use field_forms       , only : polynomial_form
  use field_forms   , only : harmonic_form
  use operation_fitting , only : pruner
  use operation_gmres  , only : gmres

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

  subroutine report(passed, message, nfail)

    logical         , intent(in)    :: passed
    character(len=*), intent(in)    :: message
    integer         , intent(inout) :: nfail

    if (passed) then
       write(*, '(a)') ' PASS : ' // message
    else
       write(*, '(a)') ' FAIL : ' // message
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! The skewed mesh, in one place. Interior deltas are the true
  ! centre-to-centre distances; none is parallel to its normal.
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
         & cell_centres = [c1(1), c1(2), 0.0_dp, &
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
         & face_centres = [0.5_dp, 0.2_dp, 0.0_dp, &
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
  ! The exact field, evaluated at a point.
  !===================================================================!

  pure real(dp) function exact_at(x, y)

    real(dp), intent(in) :: x, y

    exact_at = 2.0_dp * x + 3.0_dp * y

  end function exact_at

  !===================================================================!
  ! CHECK ZERO. The classical kernels are theorems of the
  ! exactness operation: two collinear points at degree one force
  ! the two-point stencil -+1/delta, derived rather than assumed.
  !===================================================================!

  subroutine check_derived_stencils(nfail)

    integer, intent(inout) :: nfail

    type(fit) :: fitting
    type(stored_directed_graph) :: pair
    type(stored_field)   :: positions
    class(field), allocatable :: fitted
    real(dp), allocatable :: w(:)

    pair = stored_directed_graph(2, tails=[integer ::], heads=[integer ::])
    positions = stored_field('positions', pair % vertex_set(), pair % num_vertices(), num_components=3)
    call positions % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp, &
         &                            0.5_dp, 0.0_dp, 0.0_dp])

    fitting = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call fitting % apply(pair, fitting % bind([positions]), fitted)
    call fitted % real_vector(w)

    call report(all(abs(w - [-2.0_dp, 2.0_dp]) < 1.0d-12), &
         & 'two collinear points force the two-point stencil: a theorem', nfail)

    fitting = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp], scale=3.0_dp)
    call fitting % apply(pair, fitting % bind([positions]), fitted)
    call fitted % real_vector(w)

    call report(all(abs(w - [-6.0_dp, 6.0_dp]) < 1.0d-12), &
         & 'and the conductivity multiplies the scale, as the dictionary states', nfail)

  end subroutine check_derived_stencils

  !===================================================================!
  ! CHECK ZERO-B. The basis is a free choice, and the form is
  ! restricted. A harmonic fit differentiates its own waves
  ! exactly, where polynomials only approximate them; and the
  ! pruner removes the members two collinear points cannot
  ! determine, after which the polynomial fit still reproduces
  ! its theorem.
  !===================================================================!

  subroutine check_the_form_is_free(nfail)

    integer, intent(inout) :: nfail

    type(fit)    :: wave
    type(fit)    :: poly
    type(pruner) :: restriction
    type(stored_directed_graph) :: trio, pair
    type(stored_field)   :: positions
    class(field), allocatable :: fitted
    real(dp), allocatable :: w(:)
    real(dp) :: pts(9), sampled(3), computed, expected
    integer :: j
    integer, allocatable :: retained(:)

    ! Three points on a line, one wave through them.
    pts = [0.0_dp, 0.0_dp, 0.0_dp, &
         & 0.4_dp, 0.0_dp, 0.0_dp, &
         & 0.8_dp, 0.0_dp, 0.0_dp]

    trio = stored_directed_graph(3, tails=[integer ::], heads=[integer ::])
    positions = stored_field('positions', trio % vertex_set(), trio % num_vertices(), num_components=3)
    call positions % set_real_vector(pts)

    wave = fit(harmonic_form([2.5_dp, 0.0_dp, 0.0_dp]), &
         & at=[0.4_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call wave % apply(trio, wave % bind([positions]), fitted)
    call fitted % real_vector(w)

    do j = 1, 3
       sampled(j) = sin(2.5_dp * pts(3 * j - 2))
    end do

    computed      = sum(w * sampled)
    expected = 2.5_dp * cos(2.5_dp * 0.4_dp)

    call report(abs(computed - expected) < 1.0d-10, &
         & 'the harmonic fit differentiates its own wave exactly', nfail)

    ! Two collinear points: the pruner removes what they cannot
    ! determine, and the fit still reproduces the two-point theorem.
    pair = stored_directed_graph(2, tails=[integer ::], heads=[integer ::])
    positions = stored_field('positions', pair % vertex_set(), pair % num_vertices(), num_components=3)
    call positions % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp, &
         &                            0.5_dp, 0.0_dp, 0.0_dp])

    poly = fit(polynomial_form(), at=[0.25_dp, 0.0_dp, 0.0_dp], &
         & direction=[1.0_dp, 0.0_dp, 0.0_dp])
    call restriction % adapt(poly % shape, [0.0_dp, 0.0_dp, 0.0_dp, &
         &                              0.5_dp, 0.0_dp, 0.0_dp])

    call poly % shape % members(retained)
    call report(size(retained) == 2 .and. all(retained == [1, 2]), &
         & 'the pruner removes the members the points cannot determine', nfail)

    call poly % apply(pair, poly % bind([positions]), fitted)
    call fitted % real_vector(w)
    call report(all(abs(w - [-2.0_dp, 2.0_dp]) < 1.0d-12), &
         & 'and the restricted fit still reproduces the theorem', nfail)

  end subroutine check_the_form_is_free

  !===================================================================!
  ! CHECK ONE. The two-point flux on the skewed mesh, measured
  ! against the exact normal derivative it claims to be.
  !===================================================================!

  subroutine check_two_point_misses(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(differential_operator) :: slope
    type(stored_field) :: state, fd
    type(graph) :: cells
    class(field), allocatable :: z
    real(dp), allocatable :: computed(:), deltas(:), q(:)
    real(dp) :: exact_flux(4), largest_error
    integer :: e

    m = skewed_mesh()

    cells = m % vertex_set()
    state = stored_field('q', cells, m % num_vertices())
    q = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]
    call state % set_real_vector(q)

    fd = m % face_delta()
    call fd % real_vector(deltas)

    slope = gradient(spacings=deltas)
    call slope % apply(m, slope % bind([state]), z)
    call z % real_vector(computed)

    ! What the flux claims to be: keff*area*dq/dn, here dq/dn along
    ! each interior normal.
    exact_flux = [2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp]

    largest_error = 0.0_dp
    do e = 1, 4
       largest_error = max(largest_error, abs(computed(e) - exact_flux(e)) / abs(exact_flux(e)))
    end do

    write(*, '(a, f6.1, a)') '        measured: the two-point flux is in error by up to ', &
         & 100.0_dp * largest_error, ' percent on this mesh'

    call report(largest_error > 0.25_dp, &
         & 'the two-point flux is in error on the skewed mesh, and the error is measured', nfail)

  end subroutine check_two_point_misses

  !===================================================================!
  ! CHECK TWO. The exactness weights on the same mesh, same field:
  ! every row must balance to machine precision, because the flux is
  ! exact for linears by construction and the scatter conserves.
  !===================================================================!

  subroutine check_exactness_lands(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(stencil) :: op
    type(stored_field) :: state
    type(graph) :: cells
    class(field), allocatable :: y
    real(dp), allocatable :: computed(:), q(:), vb(:), centres(:)
    type(stored_field) :: fc
    integer :: e

    m = skewed_mesh()

    ! The wall values: the exact field, read at each wall's centre.
    fc = m % face_centre()
    call fc % real_vector(centres)

    allocate(vb(12))
    vb = 0.0_dp
    do e = 5, 12
       vb(e) = exact_at(centres(3 * e - 2), centres(3 * e - 1))
    end do

    block
      type(stored_field) :: fa
      real(dp), allocatable :: farea(:)
      fa = m % face_area()
      call fa % real_vector(farea)
      op = fitted_balance_stencil(m, polynomial_form(), farea, &
           & boundary_values=vb)
    end block

    cells = m % vertex_set()
    state = stored_field('q', cells, m % num_vertices())
    q = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]
    call state % set_real_vector(q)

    call op % apply(m, op % bind([state]), y)
    call y % real_vector(computed)

    call report(maxval(abs(computed)) < 1.0d-10, &
         & 'the exactness weights balance the same field to machine zero', nfail)

  end subroutine check_exactness_lands

  !===================================================================!
  ! CHECK THREE. Solve on the skewed mesh with the exact stencil:
  ! the recovered field is the exact field. The solver does not
  ! depend on the mesh.
  !===================================================================!

  subroutine check_solver_owes_nothing(nfail)

    integer, intent(inout) :: nfail

    type(mesh) :: m
    type(stencil) :: op
    type(gmres) :: gm
    type(stored_field) :: fc
    real(dp), allocatable :: vb(:), centres(:), g(:), rhs(:), x(:)
    real(dp) :: achieved
    real(dp) :: q_exact(4)
    integer :: e

    m = skewed_mesh()

    fc = m % face_centre()
    call fc % real_vector(centres)

    allocate(vb(12))
    vb = 0.0_dp
    do e = 5, 12
       vb(e) = exact_at(centres(3 * e - 2), centres(3 * e - 1))
    end do

    block
      type(stored_field) :: fa
      real(dp), allocatable :: farea(:)
      fa = m % face_area()
      call fa % real_vector(farea)
      op = fitted_balance_stencil(m, polynomial_form(), farea, &
           & boundary_values=vb)
    end block

    call gm % state(op, m, m % vertex_set(), m % num_vertices())
    gm % tolerance = 1.0d-12

    g = gm % affine
    rhs = -g

    allocate(x(4))
    x = 0.0_dp
    call gm % solve(rhs, x, achieved)

    q_exact = [exact_at(0.0_dp, 0.0_dp), exact_at(1.0_dp, 0.45_dp), &
         & exact_at(0.15_dp, 1.1_dp), exact_at(1.3_dp, 1.6_dp)]

    call report(achieved < 1.0d-10, &
         & 'the solve converges on the skewed mesh', nfail)
    call report(all(abs(x - q_exact) < 1.0d-8), &
         & 'and recovers the exact field: the numerics do not depend on the mesh', nfail)

  end subroutine check_solver_owes_nothing

end program test_graph_robustness
