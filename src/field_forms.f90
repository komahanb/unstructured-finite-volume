!=====================================================================!
! LEVEL 1 . THE FORMS
!
! A form is a family of functions of position - a basis. A form is
! composed of two independent objects, and HAS them rather than
! being either:
!
!      evaluation      num_members, the table's width recorded where
!                      the basis is declared; values, slopes - the
!                      concretion's own table of functions, read whole
!      active basis    WHICH table entries are active, as a declared set:
!                      an identity, and a representation listing them
!
! The type once EXTENDED subset_set, which stated that a form IS a
! set of basis functions. That inheritance provided one method -
! members() - and required the whole carrier contract: a form
! implemented has(), local_index() and ambient() that nothing called,
! and could not be a set of anything else without becoming a
! different type. Composition provides the same method and requires
! nothing.
!
! Everything a member array once did, the representation does: the
! active basis members ARE the listed representation's members,
! indices into the concretion's own table. Restricting a form is
! relisting them. No second active(:) array remains, for the same
! reason as before - a set does not need two lists to record its
! members.
!
! What the form adds beyond membership is only its evaluation
! procedures, read over the FULL table, membership recording which
! entries are active:
!
!      num_members                  the table's width
!      values(x, at)            each table entry, evaluated at x,
!                               reckoned about the point `at`
!      slopes(x, at, n)         each entry's derivative along n
!
! and one operation of its own: restrict, which sets that membership.
! restrict is here rather than at the caller because a form owns its
! structure - the caller decides which member is removed, and the
! form removes it. When the form sector becomes a transform the
! restriction will return a NEW form and this procedure becomes the
! constructor it calls.
!
! Evaluating a form at a point is calculus; choosing its
! coefficients is minimization and is defined one level up.
! Polynomials are one concretion, waves another; a fit stores a form
! as data describing the basis, as an operator stores coefficients.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module field_forms

  use util_precision  , only : dp
  use view_directed     , only : SIDE_VERTEX
  use graph_fractal      , only : graph
  use map_set_representation, only : listed_set_representation

  implicit none

  private
  public :: form

  type, abstract :: form

     !----------------------------------------------------------------!
     ! WHICH basis, and WHICH entries are active. The identity is declared once
     ! by the concretion; the representation is what restrict replaces.
     !----------------------------------------------------------------!

     type(graph)                , private :: basis
     type(listed_set_representation), private :: active

     ! the table's width, stated once where the basis is declared;
     ! restriction narrows the active set, never the table
     integer, private :: width = 0

   contains

     procedure :: num_members
     procedure(form_values_interface), deferred :: values
     procedure(form_slopes_interface), deferred :: slopes
     procedure(form_derivatives_interface), deferred :: derivatives
     procedure(form_count_interface) , deferred :: dimension

     procedure :: declare_basis
     procedure :: basis_set
     procedure :: members
     procedure :: restrict

  end type form

  abstract interface

     pure subroutine form_values_interface(this, x, at, phi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(:), at(:)
       real(dp), intent(out) :: phi(:)
     end subroutine form_values_interface

     pure subroutine form_slopes_interface(this, x, at, direction, dphi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(:), at(:), direction(:)
       real(dp), intent(out) :: dphi(:)
     end subroutine form_slopes_interface
     ! the derivative of each entry of the multi-index orders, one
     ! order per coordinate, at x reckoned about at
     pure subroutine form_derivatives_interface(this, x, at, orders, dphi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(:), at(:)
       integer , intent(in)  :: orders(:)
       real(dp), intent(out) :: dphi(:)
     end subroutine form_derivatives_interface

     ! how many coordinates the form reads
     pure integer function form_count_interface(this)
       import :: form
       class(form), intent(in) :: this
     end function form_count_interface

  end interface

!=====================================================================!
! The polynomial form: every monomial in the coordinates, as many as
! the space has,
! reckoned about the point of interest, up to a degree - the Taylor
! basis at that degree, whose span is every polynomial field of it.
! Degree one is the constant and the three coordinates, and is the
! default when no degree is given. The monomials are ordered by
! total degree, and within a degree with the first coordinate's
! power falling, so degree one is 1, x, y, z in that order.
!
!=====================================================================!

  public :: polynomial_form

  type, extends(form) :: polynomial_form

     integer, allocatable, private :: power(:,:)

   contains

     procedure :: values    => polynomial_values
     procedure :: slopes    => polynomial_slopes
     procedure :: derivatives => polynomial_derivatives
     procedure :: dimension => polynomial_dimension
     procedure :: pure_members

  end type polynomial_form

  interface polynomial_form
     module procedure create_polynomial
  end interface polynomial_form

!=====================================================================!
! The harmonic form: one wave and the constant,
!
!      { 1,  sin(k . (x - at)),  cos(k . (x - at)) }
!
! whose span contains every wave of that wavenumber, whatever its
! phase. A fit over this form differentiates such waves exactly,
! where a polynomial of any finite degree only approximates them.
!
!=====================================================================!

  public :: harmonic_form

  type, extends(form) :: harmonic_form
     real(dp), allocatable :: wavenumber(:)

   contains

     procedure :: values    => harmonic_values
     procedure :: slopes    => harmonic_slopes
     procedure :: derivatives => harmonic_derivatives
     procedure :: dimension => harmonic_dimension

  end type harmonic_form

  interface harmonic_form
     module procedure create_harmonic
  end interface harmonic_form

contains

  !===================================================================!
  ! A concretion declares its basis once, with every entry of its
  ! table active. The identity is declared here so no concretion has
  ! to declare it; the active set starts full because an unrestricted
  ! form has every entry active.
  !===================================================================!

  subroutine declare_basis(this, width)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: width

    integer :: m

    call this % basis % declare()
    this % width  = width
    this % active = listed_set_representation([(m, m = 1, width)])

  end subroutine declare_basis

  !===================================================================!
  ! The table's width, as declared. Restriction does not change it:
  ! a restricted form still evaluates every entry of its table and
  ! has only some active.
  !===================================================================!

  pure integer function num_members(this)

    class(form), intent(in) :: this

    num_members = this % width

  end function num_members

  !===================================================================!
  ! WHICH basis this form's active members belong to. The identity
  ! survives restriction: restricting a form narrows the active set,
  ! and does not make the form a different basis.
  !===================================================================!

  type(graph) function basis_set(this) result(b)

    class(form), intent(in) :: this

    b = this % basis

  end function basis_set

  !===================================================================!
  ! The active members, in declaration order.
  !===================================================================!

  pure subroutine members(this, active_members)

    class(form)         , intent(in)  :: this
    integer, allocatable, intent(out) :: active_members(:)

    call this % active % members(active_members)

  end subroutine members

  !===================================================================!
  ! Activate only these table entries. The active indices name entries
  ! of the concretion's own table, and the active set is the whole
  ! record of membership.
  !
  ! restrict MUTATES, as before. Making restriction functional - a
  ! new form, a new basis identity - is a separate transformation and
  ! is not included in a type change.
  !===================================================================!

  subroutine restrict(this, active_members)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: active_members(:)

    this % active = listed_set_representation(active_members)

  end subroutine restrict


  ! Constructed with every table entry active: the members are the four.
  type(polynomial_form) function create_polynomial(degree, dimension) result(this)

    integer, intent(in), optional :: degree, dimension

    integer :: p, d, deg, m, width, k
    integer, allocatable :: alpha(:)
    character(len=250) :: message

    p = 1
    if (present(degree)) p = degree
    if (p < 0) then
       write(message,'(a,i0)') 'field_forms: a polynomial degree must be zero or above; degree = ', p
       error stop trim(message)
    end if
    d = 3
    if (present(dimension)) d = dimension
    if (d < 1) then
       write(message,'(a,i0)') 'field_forms: a polynomial must read at least one coordinate; &
            &dimension = ', d
       error stop trim(message)
    end if

    ! every multi-index of d powers summing to at most p: C(d + p, p)
    ! of them, by total degree rising and, within a degree, the first
    ! power falling, then the second, and so on
    width = 1
    do k = 1, p
       width = width * (d + k) / k
    end do
    allocate(this % power(d, width), alpha(d))

    m = 0
    do deg = 0, p
       call multi_indices(deg, 1, alpha, this % power, m)
    end do
    if (m /= width) then
       write(message,'(a,i0,a,i0)') 'field_forms: the multi-indices must fill the table; m = ', &
            & m, ', width = ', width
       error stop trim(message)
    end if

    call this % declare_basis(width)

  end function create_polynomial

  !-------------------------------------------------------------------!
  ! Every way of writing total on the powers from position first on,
  ! the earlier power falling first, each written into the next column
  ! of the table.
  !-------------------------------------------------------------------!

  pure recursive subroutine multi_indices(total, first, alpha, table, m)

    integer, intent(in)    :: total, first
    integer, intent(inout) :: alpha(:), table(:,:), m

    integer :: k

    if (first == size(alpha)) then
       alpha(first) = total
       m = m + 1
       table(:, m) = alpha
       return
    end if

    do k = total, 0, -1
       alpha(first) = k
       call multi_indices(total - k, first + 1, alpha, table, m)
    end do

  end subroutine multi_indices

  pure integer function polynomial_dimension(this)

    class(polynomial_form), intent(in) :: this

    polynomial_dimension = size(this % power, 1)

  end function polynomial_dimension

  pure subroutine polynomial_values(this, x, at, phi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:)
    real(dp), intent(out) :: phi(:)

    real(dp) :: r(size(this % power, 1))
    integer :: m, c

    r = x(1:size(r)) - at(1:size(r))
    do m = 1, size(phi)
       phi(m) = 1.0_dp
       do c = 1, size(r)
          phi(m) = phi(m) * monomial(r(c), this % power(c, m))
       end do
    end do

  end subroutine polynomial_values

  !-------------------------------------------------------------------!
  ! The derivative of each monomial along the direction: the sum over
  ! the coordinates of the direction's component times the power as
  ! a factor.
  !-------------------------------------------------------------------!

  pure subroutine polynomial_slopes(this, x, at, direction, dphi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:), direction(:)
    real(dp), intent(out) :: dphi(:)

    real(dp) :: r(size(this % power, 1)), term
    integer :: m, c, o

    r = x(1:size(r)) - at(1:size(r))
    do m = 1, size(dphi)
       dphi(m) = 0.0_dp
       do c = 1, size(r)
          if (this % power(c, m) == 0) cycle
          term = direction(c) * real(this % power(c, m), dp) &
               & * monomial(r(c), this % power(c, m) - 1)
          do o = 1, size(r)
             if (o == c) cycle
             term = term * monomial(r(o), this % power(o, m))
          end do
          dphi(m) = dphi(m) + term
       end do
    end do

  end subroutine polynomial_slopes

  !-------------------------------------------------------------------!
  ! The derivative of each monomial of the multi-index orders: for
  ! each coordinate the falling factorial of the power times the
  ! monomial of the power less the order, zero where the order
  ! exceeds the power.
  !-------------------------------------------------------------------!

  pure subroutine polynomial_derivatives(this, x, at, orders, dphi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:)
    integer , intent(in)  :: orders(:)
    real(dp), intent(out) :: dphi(:)

    real(dp) :: r(size(this % power, 1))
    integer :: m, c, k

    r = x(1:size(r)) - at(1:size(r))
    do m = 1, size(dphi)
       dphi(m) = 1.0_dp
       do c = 1, size(r)
          if (orders(c) > this % power(c, m)) then
             dphi(m) = 0.0_dp
             exit
          end if
          do k = 0, orders(c) - 1
             dphi(m) = dphi(m) * real(this % power(c, m) - k, dp)
          end do
          dphi(m) = dphi(m) * monomial(r(c), this % power(c, m) - orders(c))
       end do
    end do

  end subroutine polynomial_derivatives

  !-------------------------------------------------------------------!
  ! The members that are powers of one coordinate, the constant among
  ! them: the form restricted to these is the compact one, whose
  ! derivatives on a cell and its face neighbours are the central
  ! differences.
  !-------------------------------------------------------------------!

  pure subroutine pure_members(this, list)

    class(polynomial_form), intent(in)  :: this
    integer, allocatable  , intent(out) :: list(:)

    integer :: m

    list = [integer ::]
    do m = 1, size(this % power, 2)
       if (count(this % power(:, m) > 0) <= 1) list = [list, m]
    end do

  end subroutine pure_members

  pure real(dp) function monomial(r, power) result(v)

    real(dp), intent(in) :: r
    integer , intent(in) :: power

    if (power == 0) then
       v = 1.0_dp
    else
       v = r ** power
    end if

  end function monomial



  ! Constructed with every table entry active: the members are the three.
  type(harmonic_form) function create_harmonic(wavenumber) result(this)

    real(dp), intent(in) :: wavenumber(:)

    integer :: m

    this % wavenumber = wavenumber
    call this % declare_basis(3)

  end function create_harmonic

  pure subroutine harmonic_values(this, x, at, phi)

    class(harmonic_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:)
    real(dp), intent(out) :: phi(:)

    real(dp) :: phase

    phase = dot_product(this % wavenumber, x - at)

    phi(1) = 1.0_dp
    phi(2) = sin(phase)
    phi(3) = cos(phase)

  end subroutine harmonic_values

  !===================================================================!
  ! d/dn of a wave: the chain rule gives the factor k . n.
  !===================================================================!

  pure subroutine harmonic_slopes(this, x, at, direction, dphi)

    class(harmonic_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:), direction(:)
    real(dp), intent(out) :: dphi(:)

    real(dp) :: phase, kn

    phase = dot_product(this % wavenumber, x - at)
    kn    = dot_product(this % wavenumber, direction)

    dphi(1) = 0.0_dp
    dphi(2) =  kn * cos(phase)
    dphi(3) = -kn * sin(phase)

  end subroutine harmonic_slopes


  !===================================================================!
  ! The derivative of the multi-index orders of a wave: each order
  ! along a coordinate multiplies by that wavenumber and advances the
  ! phase by a quarter turn, so the total order n turns sin into
  ! sin, cos, -sin, -cos as n runs mod 4.
  !===================================================================!

  pure subroutine harmonic_derivatives(this, x, at, orders, dphi)

    class(harmonic_form), intent(in) :: this
    real(dp), intent(in)  :: x(:), at(:)
    integer , intent(in)  :: orders(:)
    real(dp), intent(out) :: dphi(:)

    real(dp) :: phase, factor, s, c
    integer :: n, i

    phase  = dot_product(this % wavenumber, x - at)
    factor = 1.0_dp
    do i = 1, size(orders)
       factor = factor * this % wavenumber(i) ** orders(i)
    end do
    n = sum(orders)
    select case (mod(n, 4))
    case (0); s =  sin(phase); c =  cos(phase)
    case (1); s =  cos(phase); c = -sin(phase)
    case (2); s = -sin(phase); c = -cos(phase)
    case default; s = -cos(phase); c =  sin(phase)
    end select
    dphi(1) = merge(1.0_dp, 0.0_dp, n == 0)
    dphi(2) = factor * s
    dphi(3) = factor * c

  end subroutine harmonic_derivatives

  pure integer function harmonic_dimension(this)

    class(harmonic_form), intent(in) :: this

    harmonic_dimension = size(this % wavenumber)

  end function harmonic_dimension

end module field_forms
