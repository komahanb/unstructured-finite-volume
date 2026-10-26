!=====================================================================!
! Module that states the forms, level 1 of the fields.
!
! A form is a family of functions of position: a basis. It is composed
! of two independent objects, and has them rather than being either:
!
!      evaluation      the concretion's own table of functions, read
!                      whole: its width num_members, recorded where the
!                      basis is declared, and its values and derivatives
!      active basis    which entries of the table are active, as a
!                      declared set: an identity, and a representation
!                      listing them
!
!      the table      phi_1   phi_2   phi_3   phi_4      num_members = 4
!                       1       x       y       z
!      active set       *               *       *        restrict([1, 3, 4])
!
! The type once extended subset_set, which stated that a form is a set
! of basis functions. That inheritance provided one method, members(),
! and required the whole carrier contract: a form implemented has(),
! local_index() and ambient(), which nothing called, and could not be a
! set of anything else without becoming a different type. Composition
! provides the same method and requires nothing.
!
! Everything a member array once did, the representation does: the
! active basis members are the listed representation's members,
! indices into the concretion's own table. Restricting a form is
! relisting them. No second active(:) array remains, because a set
! does not need two lists to record its members.
!
! What the form adds beyond membership is its evaluation procedures,
! read over the full table, membership recording which entries are
! active:
!
!      num_members              the table's width
!      values(x, at)            each table entry, evaluated at x,
!                               reckoned about the point at
!      derivatives(x, at, o)    each entry's derivative of the
!                               multi-index o
!
! and one operation of its own: restrict, which sets that membership.
! restrict is here rather than at the caller because a form owns its
! structure: the caller decides which member is removed, and the form
! removes it. When the form sector becomes a transform, restriction
! will return a new form and this procedure becomes the constructor it
! calls.
!
! Evaluating a form at a point is calculus; choosing its coefficients
! is minimization and is defined one level up. Polynomials are one
! concretion, waves another; a fit stores a form as data describing
! the basis, as an operator stores coefficients.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module field_forms

  use util_precision        , only : dp
  use view_directed         , only : SIDE_VERTEX
  use graph_fractal         , only : graph
  use map_set_representation, only : listed_set_representation

  implicit none

  private

  public :: form

  !-------------------------------------------------------------------!
  ! Type that models a form: the identity of its basis, the active
  ! entries of its table, and the table's width. The concretion states
  ! the table itself through values, derivatives and dimension.
  !-------------------------------------------------------------------!

  type, abstract :: form

     !----------------------------------------------------------------!
     ! Which basis, and which entries are active. The identity is
     ! declared once by the concretion; the representation is what
     ! restrict replaces.
     !----------------------------------------------------------------!

     type(graph)                    , private :: basis
     type(listed_set_representation), private :: active

     ! The table's width, stated once where the basis is declared;
     ! restriction narrows the active set, never the table
     integer, private :: width = 0

   contains

     ! The table, stated by the concretion
     procedure                                      :: num_members
     procedure(form_values_interface)     , deferred :: values
     procedure(form_derivatives_interface), deferred :: derivatives
     procedure(form_count_interface)      , deferred :: dimension

     ! The basis and its active members
     procedure :: declare_basis
     procedure :: basis_set
     procedure :: members
     procedure :: restrict

  end type form

  abstract interface

     !----------------------------------------------------------------!
     ! Interface for the value of each table entry at x, reckoned
     ! about at
     !----------------------------------------------------------------!

     pure subroutine form_values_interface(this, x, at, phi)
       import :: form, dp
       class(form), intent(in)  :: this
       real(dp)   , intent(in)  :: x(:), at(:)
       real(dp)   , intent(out) :: phi(:)
     end subroutine form_values_interface

     !----------------------------------------------------------------!
     ! Interface for the derivative of each entry of the multi-index
     ! orders, one order per coordinate, at x reckoned about at
     !----------------------------------------------------------------!

     pure subroutine form_derivatives_interface(this, x, at, orders, dphi)
       import :: form, dp
       class(form), intent(in)  :: this
       real(dp)   , intent(in)  :: x(:), at(:)
       integer    , intent(in)  :: orders(:)
       real(dp)   , intent(out) :: dphi(:)
     end subroutine form_derivatives_interface

     !----------------------------------------------------------------!
     ! Interface for how many coordinates the form reads
     !----------------------------------------------------------------!

     pure integer function form_count_interface(this)
       import :: form
       class(form), intent(in) :: this
     end function form_count_interface

  end interface

  !-------------------------------------------------------------------!
  ! Type that models the polynomial form: every monomial in the
  ! coordinates, as many as the space has, reckoned about the point of
  ! interest, up to a degree. It is the Taylor basis at that degree,
  ! whose span is every polynomial field of it. Degree one is the
  ! constant and the three coordinates, and is the default when no
  ! degree is given.
  !
  ! The monomials are ordered by total degree, and within a degree with
  ! the first coordinate's power falling, so degree one is 1, x, y, z in
  ! that order. The table power(d, width) stores one column per
  ! monomial:
  !
  !      degree 2 in two coordinates, C(4, 2) = 6 monomials
  !
  !      column        1     2     3     4      5      6
  !      monomial      1     x     y     x^2    xy     y^2
  !      power(1,:)    0     1     0     2      1      0
  !      power(2,:)    0     0     1     0      1      2
  !-------------------------------------------------------------------!

  public :: polynomial_form

  type, extends(form) :: polynomial_form

     integer, allocatable, private :: power(:,:)

   contains

     procedure :: values      => polynomial_values
     procedure :: derivatives => polynomial_derivatives
     procedure :: dimension   => polynomial_dimension
     procedure :: pure_members

  end type polynomial_form

  interface polynomial_form
     module procedure create_polynomial
  end interface polynomial_form

  !-------------------------------------------------------------------!
  ! Type that models the harmonic form: one wave and the constant,
  !
  !      { 1,  sin(k . (x - at)),  cos(k . (x - at)) }
  !
  ! whose span contains every wave of that wavenumber, whatever its
  ! phase. A fit over this form differentiates such waves exactly,
  ! where a polynomial of any finite degree only approximates them.
  !-------------------------------------------------------------------!

  public :: harmonic_form

  type, extends(form) :: harmonic_form

     real(dp), allocatable :: wavenumber(:)

   contains

     procedure :: values      => harmonic_values
     procedure :: derivatives => harmonic_derivatives
     procedure :: dimension   => harmonic_dimension

  end type harmonic_form

  interface harmonic_form
     module procedure create_harmonic
  end interface harmonic_form

contains

  !===================================================================!
  ! Declares the basis of a concretion once, with every entry of its
  ! table active. The identity is declared here so that no concretion
  ! has to declare it; the active set starts full because an
  ! unrestricted form has every entry active.
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
  ! Returns the table's width, as declared. Restriction does not change
  ! it: a restricted form still evaluates every entry of its table and
  ! has only some active.
  !===================================================================!

  pure integer function num_members(this)

    class(form), intent(in) :: this

    num_members = this % width

  end function num_members

  !===================================================================!
  ! Returns which basis this form's active members belong to. The
  ! identity survives restriction: restricting a form narrows the
  ! active set, and does not make the form a different basis.
  !===================================================================!

  type(graph) function basis_set(this) result(b)

    class(form), intent(in) :: this

    b = this % basis

  end function basis_set

  !===================================================================!
  ! Returns the active members, in declaration order
  !===================================================================!

  pure subroutine members(this, active_members)

    class(form)         , intent(in)  :: this
    integer, allocatable, intent(out) :: active_members(:)

    call this % active % members(active_members)

  end subroutine members

  !===================================================================!
  ! Activates only these table entries. The active indices name entries
  ! of the concretion's own table, and the active set is the whole
  ! record of membership.
  !
  ! restrict mutates, as before. Making restriction functional, a new
  ! form with a new basis identity, is a separate transformation and is
  ! not included in a type change.
  !===================================================================!

  subroutine restrict(this, active_members)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: active_members(:)

    this % active = listed_set_representation(active_members)

  end subroutine restrict

  !===================================================================!
  ! Creates the polynomial form of a degree, one by default, in a number
  ! of coordinates, three by default, with every table entry active:
  ! C(d + p, p) members, the four 1, x, y, z at the defaults.
  !===================================================================!

  type(polynomial_form) function create_polynomial(degree, dimension) result(this)

    integer, intent(in), optional :: degree, dimension

    integer              :: p, d, deg, m, width, k
    integer, allocatable :: alpha(:)
    character(len=250)   :: message

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

    ! Every multi-index of d powers summing to at most p: C(d + p, p)
    ! of them, by total degree rising and, within a degree, the first
    ! power falling, then the second, and so on. The width is formed as
    ! the product of (d + k) / k for k = 1 .. p, each partial product
    ! the binomial C(d + k, k) and so exact in integers.
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

  !===================================================================!
  ! Writes every way of distributing total over the powers from
  ! position first on, the earlier power falling first, each into the
  ! next column of the table. For total = 2 over two powers:
  !
  !      alpha(1) = 2  ->  (2, 0)
  !      alpha(1) = 1  ->  (1, 1)
  !      alpha(1) = 0  ->  (0, 2)
  !===================================================================!

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

  !===================================================================!
  ! Returns how many coordinates the polynomial form reads
  !===================================================================!

  pure integer function polynomial_dimension(this)

    class(polynomial_form), intent(in) :: this

    polynomial_dimension = size(this % power, 1)

  end function polynomial_dimension

  !===================================================================!
  ! Returns each monomial at x, reckoned about at: the product over the
  ! coordinates of r_c ** power(c, m), with r = x - at
  !===================================================================!

  pure subroutine polynomial_values(this, x, at, phi)

    class(polynomial_form), intent(in)  :: this
    real(dp)              , intent(in)  :: x(:), at(:)
    real(dp)              , intent(out) :: phi(:)

    real(dp) :: r(size(this % power, 1))
    integer  :: m, c

    r = x(1:size(r)) - at(1:size(r))
    do m = 1, size(phi)
       phi(m) = 1.0_dp
       do c = 1, size(r)
          phi(m) = phi(m) * monomial(r(c), this % power(c, m))
       end do
    end do

  end subroutine polynomial_values

  !===================================================================!
  ! Returns the derivative of each monomial of the multi-index orders:
  ! for each coordinate, the falling factorial of the power times the
  ! monomial of the power less the order, zero where the order exceeds
  ! the power,
  !
  !      d^o / dr^o  r^a  =  a (a - 1) ... (a - o + 1)  r^(a - o)
  !===================================================================!

  pure subroutine polynomial_derivatives(this, x, at, orders, dphi)

    class(polynomial_form), intent(in)  :: this
    real(dp)              , intent(in)  :: x(:), at(:)
    integer               , intent(in)  :: orders(:)
    real(dp)              , intent(out) :: dphi(:)

    real(dp) :: r(size(this % power, 1))
    integer  :: m, c, k

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

  !===================================================================!
  ! Returns the members that are powers of one coordinate, the constant
  ! among them. The form restricted to these is the compact one, whose
  ! derivatives on a cell and its face neighbours are the central
  ! differences. At degree 2 in two coordinates:
  !
  !      monomial      1     x     y     x^2    xy     y^2
  !      pure          *     *     *     *             *
  !===================================================================!

  pure subroutine pure_members(this, list)

    class(polynomial_form), intent(in)  :: this
    integer, allocatable  , intent(out) :: list(:)

    integer :: m

    list = [integer ::]
    do m = 1, size(this % power, 2)
       if (count(this % power(:, m) > 0) <= 1) list = [list, m]
    end do

  end subroutine pure_members

  !===================================================================!
  ! Returns r ** power, one for the power zero whatever r is
  !===================================================================!

  pure real(dp) function monomial(r, power) result(v)

    real(dp), intent(in) :: r
    integer , intent(in) :: power

    if (power == 0) then
       v = 1.0_dp
    else
       v = r ** power
    end if

  end function monomial

  !===================================================================!
  ! Creates the harmonic form of a wavenumber, with every table entry
  ! active: the three members 1, sin and cos
  !===================================================================!

  type(harmonic_form) function create_harmonic(wavenumber) result(this)

    real(dp), intent(in) :: wavenumber(:)

    integer :: m

    this % wavenumber = wavenumber
    call this % declare_basis(3)

  end function create_harmonic

  !===================================================================!
  ! Returns the constant and the wave's sine and cosine at x, the phase
  ! k . (x - at)
  !===================================================================!

  pure subroutine harmonic_values(this, x, at, phi)

    class(harmonic_form), intent(in)  :: this
    real(dp)            , intent(in)  :: x(:), at(:)
    real(dp)            , intent(out) :: phi(:)

    real(dp) :: phase

    phase = dot_product(this % wavenumber, x - at)

    phi(1) = 1.0_dp
    phi(2) = sin(phase)
    phi(3) = cos(phase)

  end subroutine harmonic_values

  !===================================================================!
  ! Returns the derivative of the multi-index orders of the wave. Each
  ! order along a coordinate multiplies by that wavenumber and advances
  ! the phase by a quarter turn, so the total order n turns the pair as
  ! n runs mod 4:
  !
  !      n mod 4        0        1        2        3
  !      sin   ->      sin      cos     -sin     -cos
  !      cos   ->      cos     -sin     -cos      sin
  !
  ! The constant's derivative is zero for every order above zero.
  !===================================================================!

  pure subroutine harmonic_derivatives(this, x, at, orders, dphi)

    class(harmonic_form), intent(in)  :: this
    real(dp)            , intent(in)  :: x(:), at(:)
    integer             , intent(in)  :: orders(:)
    real(dp)            , intent(out) :: dphi(:)

    real(dp) :: phase, factor, s, c
    integer  :: n, i

    phase  = dot_product(this % wavenumber, x - at)
    factor = 1.0_dp
    do i = 1, size(orders)
       factor = factor * this % wavenumber(i) ** orders(i)
    end do
    n = sum(orders)
    select case (mod(n, 4))
    case (0)    ; s =  sin(phase); c =  cos(phase)
    case (1)    ; s =  cos(phase); c = -sin(phase)
    case (2)    ; s = -sin(phase); c = -cos(phase)
    case default; s = -cos(phase); c =  sin(phase)
    end select
    dphi(1) = merge(1.0_dp, 0.0_dp, n == 0)
    dphi(2) = factor * s
    dphi(3) = factor * c

  end subroutine harmonic_derivatives

  !===================================================================!
  ! Returns how many coordinates the harmonic form reads: the size of
  ! its wavenumber
  !===================================================================!

  pure integer function harmonic_dimension(this)

    class(harmonic_form), intent(in) :: this

    harmonic_dimension = size(this % wavenumber)

  end function harmonic_dimension

end module field_forms
