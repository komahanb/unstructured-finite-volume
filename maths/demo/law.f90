! The adjoint law on a space graph, through the one compose and the
! one transpose the marcher also uses.
program law
  use iso_fortran_env, only : dp => REAL64
  use field_calculus, only : field
  use view_directed_stored, only : stored_directed_graph
  use field_stored, only : stored_field
  use operation_stencil, only : stencil, compose
  use operation_differential, only : gradient, divergence, laplacian, &
       & vertex_derivative, stencil_of
  implicit none

  type(stored_directed_graph) :: g
  type(stencil) :: grad, div, lap, lap_t, chain, chain_t, flag
  type(stored_field) :: q, z, p
  class(field), allocatable :: out
  real(dp), allocatable :: y(:), w(:), a(:), b(:)
  real(dp) :: lhs, rhs
  integer :: nfail, nv, ne, i

  nfail = 0

  ! a chain of five cells with a boundary face at each end:
  ! edges 1..4 interior, edges 5 and 6 leave the chain (no head)
  nv = 5
  g = stored_directed_graph(nv, tails=[1,2,3,4,1,5], heads=[2,3,4,5,6,7])
  ne = g % num_edges()

  q = stored_field('q', g % vertex_set(), nv); call q % set_real_vector([1.0_dp, 4.0_dp, 9.0_dp, 16.0_dp, 25.0_dp])
  z = stored_field('z', g % edge_set(), ne);   call z % set_real_vector([(0.5_dp * i, i = 1, ne)])
  p = stored_field('p', g % vertex_set(), nv); call p % set_real_vector([3.0_dp, -1.0_dp, 2.0_dp, 0.5_dp, 7.0_dp])

  ! 1. the edge landing compiles to a rectangular stencil
  grad = stencil_of(gradient(spacing=0.5_dp, boundary_value=2.0_dp), g)
  call check(grad % num_rows == ne .and. grad % num_columns == nv, 'gradient is edges x vertices')

  ! 2. rectangular adjointness: <G q, z>_E = <q, G^T z>_V on the linear part
  call grad % apply(g, [q], out); call out % real_vector(y)
  lhs = sum((y - grad % constants) * [(0.5_dp * i, i = 1, ne)])
  flag = grad % transpose()
  call flag % apply(g, [z], out); call out % real_vector(w)
  rhs = sum([1.0_dp, 4.0_dp, 9.0_dp, 16.0_dp, 25.0_dp] * w)
  call check(abs(lhs - rhs) < 1.0d-12, '<G q, z> = <q, G^T z> across edges x vertices')
  call check(flag % num_rows == nv .and. flag % num_columns == ne, 'the transpose is vertices x edges')

  ! 3. (D o G)^T = G^T o D^T: both sides through compose and transpose
  div   = stencil_of(divergence(measure=0.25_dp), g, on_edge_field=.true.)
  chain = compose(div, grad)
  chain_t = chain % transpose()
  flag = compose(grad % transpose(), div % transpose())
  call check(div % num_rows == nv .and. div % num_columns == ne, 'the incidence step is vertices x edges')
  call check(same_map(chain_t, flag), '(D o G)^T = G^T o D^T')

  ! 4. the composed chain is the laplacian the operator compiles
  lap = stencil_of(laplacian(spacing=0.5_dp, measure=0.25_dp, boundary_value=2.0_dp), g)
  call check(same_map(chain, lap), 'D o G is the compiled laplacian, constant included')

  ! 5. the adjoint flag and the transpose are one map
  lap_t = stencil_of(vertex_derivative(2, spacing=0.5_dp, measure=0.25_dp, boundary_value=2.0_dp, adjoint=.true.), g)
  call check(same_map(lap % transpose(), lap_t), 'adjoint=.true. is the transpose of the compiled chain')

  ! 6. L^T^T = L on the linear part
  flag = lap % transpose()
  flag = flag % transpose()
  call check(same_linear(flag, lap), 'L^T^T = L on the linear part')

  if (nfail == 0) then
     print '(a)', ' PASS : the adjoint law holds on the space graph through the one compose and transpose'
  else
     print '(a,i0,a)', ' FAIL : ', nfail, ' checks'
     error stop
  end if

contains

  subroutine check(ok, what)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: what
    if (ok) then
       print '(a)', ' PASS : ' // what
    else
       print '(a)', ' FAIL : ' // what
       nfail = nfail + 1
    end if
  end subroutine check

  ! equality of maps by action on a basis: apply both to every basis
  ! vector of the column side and to zero
  logical function same_linear(x, y) result(same)
    type(stencil), intent(in) :: x, y
    integer :: j
    same = x % num_rows == y % num_rows .and. x % num_columns == y % num_columns
    if (.not. same) return
    do j = 1, x % num_columns
       call column(x, j, a); call column(y, j, b)
       if (maxval(abs(a - b)) > 1.0d-12) same = .false.
    end do
  end function same_linear

  logical function same_map(x, y) result(same)
    type(stencil), intent(in) :: x, y
    same = same_linear(x, y)
    if (same) same = maxval(abs(x % constants - y % constants)) < 1.0d-12
  end function same_map

  subroutine column(x, j, col)
    type(stencil), intent(in) :: x
    integer, intent(in) :: j
    real(dp), allocatable, intent(out) :: col(:)
    real(dp), allocatable :: e(:)
    integer :: k
    allocate(col(x % num_rows), e(x % num_columns))
    e = 0.0_dp; e(j) = 1.0_dp
    col = 0.0_dp
    do k = 1, size(x % rows)
       col(x % rows(k)) = col(x % rows(k)) + x % weights(k) * e(x % columns(k))
    end do
  end subroutine column

end program law
