!=====================================================================!
! THE PRINTED FORMS OF A LAGRANGIAN. The objective and every equation
! of every term of a continuous residual are the summands of one
! expression_view, the view of each paired with the component of its
! multiplier, and that view is written in each of its forms: the
! tree, the diagram, the brackets, the digraph, and the formula, each
! summand on a line of its own,
!
!      L = integral over the time factor of (u*u + dt(u)**2)/2.0
!        + lambda_omega*(dt(dt(u)) - nu*(1.0 - u*u)*dt(u) + u)
!        + lambda_ic(1)*(u - 2.0)
!
! An unknown is written by the name its manifold declares, with the
! component's number for an unknown of several components; a design
! coordinate by its name; a multiplier by its name, with the number
! of the equation it pairs for a term of several equations. An
! equation summed without a multiplier is written alone; an
! inequality as its pairing, the sign of its multiplier not written.
!
! tree(M_h) prints the discrete minimizer M_h paired with a discrete
! Lagrangian: the march over the blocks, and the minimizer of a block
! with its inner minimizers (operation_minimization).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_residual

  use util_precision        , only : dp
  use token_identity        , only : token
  use operation_field       , only : continuous_field
  use operation_field       , only : WHOLE, TIME_FACE, TIME_FACTOR, SPACE_FACTOR, DESIGN_FACTOR
  use operation_lagrangian  , only : continuous_lagrangian
  use operation_minimization, only : discrete_minimizer
  use view_expression       , only : expression_view, operator(+), operator(*)

  implicit none

  private
  public :: formula, tree, diagram, brackets, digraph

  interface tree
     module procedure lagrangian_tree, minimizer_tree
  end interface tree

  ! the length of a written name: a manifold's name of 32 characters
  ! and the number of its component
  integer, parameter :: MAX_LABEL = 44

contains

  function formula(L) result(text)

    type(continuous_lagrangian), intent(in) :: L
    character(len=:), allocatable :: text

    type(expression_view), allocatable :: s(:)
    character(len=:), allocatable :: line
    integer :: k

    s    = summands(L)
    text = 'L = '
    do k = 1, size(s)
       line = s(k) % formula()
       if (k > 1) text = text // new_line('a') // '  + '
       text = text // line
    end do

  end function formula

  function lagrangian_tree(L) result(text)

    type(continuous_lagrangian), intent(in) :: L
    character(len=:), allocatable :: text

    type(expression_view) :: view

    view = lagrangian_view(L)
    text = view % tree()

  end function lagrangian_tree

  function minimizer_tree(M_h) result(text)

    type(discrete_minimizer), intent(in) :: M_h
    character(len=:), allocatable :: text

    text = M_h % tree()

  end function minimizer_tree

  function diagram(L) result(text)

    type(continuous_lagrangian), intent(in) :: L
    character(len=:), allocatable :: text

    type(expression_view) :: view

    view = lagrangian_view(L)
    text = view % diagram()

  end function diagram

  function brackets(L) result(text)

    type(continuous_lagrangian), intent(in) :: L
    character(len=:), allocatable :: text

    type(expression_view) :: view

    view = lagrangian_view(L)
    text = view % brackets()

  end function brackets

  function digraph(L) result(text)

    type(continuous_lagrangian), intent(in) :: L
    character(len=:), allocatable :: text

    type(expression_view) :: view

    view = lagrangian_view(L)
    text = view % digraph()

  end function digraph

  ! the sum of the summands, in their order
  function lagrangian_view(L) result(view)

    type(continuous_lagrangian), intent(in) :: L
    type(expression_view) :: view

    type(expression_view), allocatable :: s(:)
    integer :: k

    s    = summands(L)
    view = s(1)
    do k = 2, size(s)
       view = view + s(k)
    end do

  end function lagrangian_view

  !===================================================================!
  ! THE SUMMANDS: the objective, then each equation of each term in
  ! the order of the terms, times the component of the term's
  ! multiplier where the term is paired.
  !===================================================================!

  function summands(L) result(s)

    type(continuous_lagrangian), intent(in) :: L
    type(expression_view), allocatable :: s(:)

    character(len=MAX_LABEL), allocatable :: unknowns(:), designs(:)
    type(expression_view) :: summand
    integer :: k, i

    allocate(s(0))
    if (allocated(L % objective)) then
       call names_on(L, L % objective % on, unknowns, designs)
       s = [s, field_view(L % objective, unknowns, designs)]
    end if
    do k = 1, size(L % term)
       associate (t => L % term(k))
         do i = 1, size(t % equation)
            call names_on(L, t % equation(i) % on, unknowns, designs)
            summand = field_view(t % equation(i), unknowns, designs)
            if (t % paired) summand = expression_view(component_name(t % multiplier % name, i, &
                 & size(t % equation))) * summand
            s = [s, summand]
         end do
       end associate
    end do

  end function summands

  !===================================================================!
  ! THE VIEW OF A FIELD: its components, separated by commas; the
  ! integral over the region added to them; and, for an integral, the
  ! factor it is taken over, or for a smooth maximum the factor and
  ! the weight, written before them.
  !===================================================================!

  function field_view(f, unknowns, designs) result(view)

    type(continuous_field), intent(in) :: f
    character(len=*)      , intent(in) :: unknowns(:), designs(:)
    type(expression_view) :: view

    character(len=32) :: weight
    integer :: c

    view = expression_view(f % component(1), unknowns, designs)
    do c = 2, size(f % component)
       view = expression_view(',', view, expression_view(f % component(c), unknowns, designs))
    end do
    if (allocated(f % integrand)) then
       view = view + expression_view('integral over the ' // part_name(SPACE_FACTOR) // ' of', &
            & expression_view(f % integrand(1), unknowns, designs))
    end if
    if (.not. f % is_integral()) return
    if (f % weight > 0.0_dp) then
       write(weight, '(g0.6)') f % weight
       view = expression_view('smooth maximum of weight ' // trim(adjustl(weight)) // ' over the ' // &
            & part_name(f % integrated_part) // ' of', view)
    else
       view = expression_view('integral over the ' // part_name(f % integrated_part) // ' of', view)
    end if

  end function field_view

  ! the part of a manifold a field is integrated over
  function part_name(part) result(name)

    integer, intent(in) :: part
    character(len=:), allocatable :: name

    select case (part)
    case (WHOLE);         name = 'manifold'
    case (TIME_FACE);     name = 'face of time'
    case (TIME_FACTOR);   name = 'time factor'
    case (SPACE_FACTOR);  name = 'space factor'
    case (DESIGN_FACTOR); name = 'design factor'
    case default;         name = 'face of space'
    end select

  end function part_name

  !===================================================================!
  ! THE NAMES of the unknowns and of the design coordinates of the
  ! manifold a field is on, read from the first term on that
  ! manifold. Without such a term no names are returned, and the
  ! leaves are written q<f> and nu.
  !===================================================================!

  subroutine names_on(L, on, unknowns, designs)

    type(continuous_lagrangian), intent(in) :: L
    type(token)              , intent(in) :: on
    character(len=MAX_LABEL), allocatable, intent(out) :: unknowns(:), designs(:)

    integer :: k, f

    allocate(unknowns(0), designs(0))
    do k = 1, size(L % term)
       associate (m => L % term(k) % manifold)
         if (.not. m % same_support(on)) cycle
         deallocate(unknowns)
         allocate(unknowns(size(m % unknown_name)))
         do f = 1, size(unknowns)
            unknowns(f) = component_name(m % unknown_name(f), f - findloc(m % unknown_name, m % unknown_name(f), 1) &
                 & + 1, count(m % unknown_name == m % unknown_name(f)))
         end do
         if (allocated(m % design_parameter)) designs = [character(len=MAX_LABEL) :: m % design_parameter % name]
         return
       end associate
    end do

  end subroutine names_on

  ! the name of component k of a field of n components: the name
  ! alone for one component, name(k) for several
  function component_name(name, k, n) result(text)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: k, n
    character(len=:), allocatable :: text

    character(len=12) :: buffer

    text = trim(name)
    if (n == 1) return
    write(buffer, '(i0)') k
    text = text // '(' // trim(buffer) // ')'

  end function component_name

end module view_residual
