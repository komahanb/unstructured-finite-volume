!=====================================================================!
! The marching suite: rung 6's acceptance, with a fractal for a
! birth certificate.
!
! Three verdicts. TIME IS A GRAPH: the marcher's instants stand as
! a chain - one vertex per instant, one edge per step, walked in
! order. THE STEP IS EXACT ABOUT ITSELF: euler on decay dq/dt = -q
! lands on (1 - h)^n to machine precision, euler's own discrete
! truth. And THE MAP IS THE MARCH: z -> z^2 + c is forward euler
! with step one on S = z - z^2 - c - an identity, not an
! approximation - so the escape-time fractal falls out of the
! marcher with known points landing where sixty years of arithmetic
! say they land:
!
!      c = 0        never escapes: the origin is fixed
!      c = -1       never escapes: a two-cycle, 0 -> -1 -> 0
!      c = 1        escapes at step three: 0, 1, 2, 5
!      c = -2       never escapes: walks to the fixed point 2
!      c = 2i       escapes at step two: -4 + 2i leaves the circle
!=====================================================================!

program test_graph_marching

  use iso_fortran_env, only : dp => REAL64
  use graph_grammar  , only : graph, graph_field, graph_operation
  use graph_calculus , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use class_graph        , only : stored_graph
  use class_graph_differential_operator, only : vertex_differential_operator
  use class_graph_differential_operator, only : differential_operator
  use class_graph_marcher, only : marcher
  use mandelbrot_law_fixture, only : mandelbrot_law

  implicit none

  integer :: nfail

  nfail = 0

  call check_time_is_a_graph(nfail)
  call check_euler_is_exact_about_itself(nfail)
  call check_the_map_is_the_march(nfail)

  write(*, '(a)') ' ============================================='
  if (nfail == 0) then
     write(*, '(a)') ' all marching checks passed'
  else
     write(*, '(a, i0, a)') ' ', nfail, ' marching checks FAILED'
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
  ! VERDICT ONE. The instants stand as a chain.
  !===================================================================!

  subroutine check_time_is_a_graph(nfail)

    integer, intent(inout) :: nfail

    type(marcher) :: clock
    type(stored_graph) :: chain
    logical :: ordered
    integer :: e

    call clock % instants(10, chain)

    call report(chain % num_vertices() == 11 .and. chain % num_edges() == 10, &
         & 'eleven instants, ten steps: the chain stands', nfail)

    ordered = .true.
    do e = 1, 10
       if (chain % edge_tail(e) /= e .or. chain % edge_head(e) /= e + 1) &
            & ordered = .false.
    end do
    call report(ordered, 'and every step leads to the next instant', nfail)

  end subroutine check_time_is_a_graph

  !===================================================================!
  ! VERDICT TWO. Euler's discrete truth, hit exactly: on S = q the
  ! walk gives q_n = (1 - h)^n q_0.
  !===================================================================!

  subroutine check_euler_is_exact_about_itself(nfail)

    integer, intent(inout) :: nfail

    type(marcher) :: clock
    type(support) :: cells
    type(stored_graph) :: lone
    type(differential_operator) :: decay
    real(dp) :: q(1), expected
    integer :: v

    lone = stored_graph(1, tails=[integer ::], heads=[integer ::])
    associate (u1 => cells, u2 => v)
    end associate

    decay = vertex_differential_operator(order=0, coefficient=1.0_dp)

    clock % step = 0.125_dp
    q = [3.0_dp]
    call clock % march(decay, lone, q, 20)

    expected = 3.0_dp * (1.0_dp - 0.125_dp)**20

    call report(abs(q(1) - expected) < 1.0d-14, &
         & 'euler lands on (1-h)^n exactly: its own discrete truth', nfail)

  end subroutine check_euler_is_exact_about_itself

  !===================================================================!
  ! VERDICT THREE. Five points of the complex plane, marched by
  ! z -> z^2 + c with the escape watched. The law is a level-3
  ! fixture; the marcher neither knows nor cares that the picture
  ! is a fractal.
  !===================================================================!

  subroutine check_the_map_is_the_march(nfail)

    integer, intent(inout) :: nfail

    type(marcher) :: clock
    type(mandelbrot_law) :: law
    type(stored_graph) :: points
    type(support) :: cells
    type(field) :: escape_field
    real(dp), allocatable :: q(:)
    integer, allocatable :: escape(:)
    integer :: v, n
    integer, parameter :: nv = 5, nmax = 30

    ! The five points, as a graph of lone cells.
    points = stored_graph(nv, tails=[integer ::], heads=[integer ::])

    law % creal = [0.0_dp, -1.0_dp, 1.0_dp, -2.0_dp, 0.0_dp]
    law % cimag = [0.0_dp,  0.0_dp, 0.0_dp,  0.0_dp, 2.0_dp]

    allocate(q(2 * nv), escape(nv))
    q      = 0.0_dp
    escape = 0

    clock % step = 1.0_dp

    do n = 1, nmax
       call clock % march(law, points, q, 1)
       do v = 1, nv
          if (escape(v) == 0 .and. &
               & q(2 * v - 1)**2 + q(2 * v)**2 > 4.0_dp) then
             ! The verdict is in; the orbit is frozen so the
             ! remaining march cannot overflow on its way out.
             escape(v) = n
             law % creal(v) = 0.0_dp
             law % cimag(v) = 0.0_dp
             q(2 * v - 1 : 2 * v) = 0.0_dp
          end if
       end do
    end do

    ! The counts leave as an integer field: the tower's native kind
    ! for a picture of whole numbers.
    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    escape_field = field('escape time', cells)
    call escape_field % set_integer_vector(escape)

    call report(escape(1) == 0, 'c = 0 never escapes: the origin holds', nfail)
    call report(escape(2) == 0, 'c = -1 never escapes: the two-cycle turns', nfail)
    call report(escape(3) == 3, 'c = 1 escapes at step three: 0, 1, 2, 5', nfail)
    call report(escape(4) == 0, 'c = -2 never escapes: pinned at two', nfail)
    call report(escape(5) == 2, &
         & 'c = 2i escapes at step two: -4 + 2i leaves the circle', nfail)

  end subroutine check_the_map_is_the_march

end program test_graph_marching
