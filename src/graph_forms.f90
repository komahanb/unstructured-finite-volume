!=====================================================================!
! LEVEL 1 . THE FORMS
!
! A form is a family of functions of position - a basis shape. It is
! two independent things held together, and it HAS them rather than
! being either:
!
!      evaluation      size_of, values, slopes - the concretion's
!                      own table of functions, read whole
!      active basis    WHICH table entries stand, as a declared set:
!                      an identity, and a representation listing them
!
! It once EXTENDED subset_set, which said a form IS a set of basis
! functions. That inheritance bought one method - members() - and
! charged for the whole carrier contract: a form answered has(),
! local_index() and ambient() that nothing asked, and could not be a
! set of anything else without becoming a different type. Composition
! buys the same method and charges for nothing.
!
! Everything a roster once did, the representation does: the standing
! basis members ARE the listed representation's members, indices into
! the concretion's own table. Pruning a form is relisting them. No
! second active(:) array survives, for the same reason as before - a
! set does not need two lists to say who belongs to it.
!
! What the form adds beyond membership is only its evaluation
! symbols, read over the FULL table, membership saying who stands:
!
!      size_of                  the table's width
!      values(x, at)            each table entry, evaluated at x,
!                               reckoned about the point `at`
!      slopes(x, at, n)         each entry's derivative along n
!
! and one act of its own: restrict, which sets that membership. It
! is here rather than at the caller because a citizen's structure is
! its own business - whoever decides a member should go says so, and
! the form does it. When the form sector becomes a transform the
! restriction will hand back a NEW form and this verb becomes the
! constructor it calls.
!
! Evaluating a form at a point is calculus; choosing its
! coefficients is minimization and lives one level up. Polynomials
! are one concretion, waves another; a fit holds a form the way an
! operator holds coefficients - as data about shape.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_forms

  use iso_fortran_env    , only : dp => REAL64
  use graph_directed_view     , only : GRAPH_SIDE_VERTEX
  use fractal_graph      , only : set_graph => graph
  use graph_set_representation, only : listed_set_representation

  implicit none

  private
  public :: form

  type, abstract :: form

     !----------------------------------------------------------------!
     ! WHICH basis, and WHO stands in it. The identity is declared once
     ! by the concretion; the representation is what restrict replaces.
     !----------------------------------------------------------------!

     type(set_graph)                , private :: basis
     type(listed_set_representation), private :: active

   contains

     procedure(form_size_interface)  , deferred :: size_of
     procedure(form_values_interface), deferred :: values
     procedure(form_slopes_interface), deferred :: slopes

     procedure :: declare_basis
     procedure :: basis_set
     procedure :: members
     procedure :: restrict

  end type form

  abstract interface

     pure integer function form_size_interface(this)
       import :: form
       class(form), intent(in) :: this
     end function form_size_interface

     pure subroutine form_values_interface(this, x, at, phi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(3), at(3)
       real(dp), intent(out) :: phi(:)
     end subroutine form_values_interface

     pure subroutine form_slopes_interface(this, x, at, direction, dphi)
       import :: form, dp
       class(form), intent(in) :: this
       real(dp), intent(in)  :: x(3), at(3), direction(3)
       real(dp), intent(out) :: dphi(:)
     end subroutine form_slopes_interface

  end interface

!=====================================================================!
! The polynomial form: the constant and the three coordinates,
! reckoned about the point of interest - the Taylor shape at degree
! one, whose span is every linear field.
!
!=====================================================================!

  public :: polynomial_form

  type, extends(form) :: polynomial_form

   contains

     procedure :: size_of => polynomial_size
     procedure :: values  => polynomial_values
     procedure :: slopes  => polynomial_slopes

  end type polynomial_form

  interface polynomial_form
     module procedure create_polynomial
  end interface polynomial_form

!=====================================================================!
! The harmonic form: one wave and the constant,
!
!      { 1,  sin(k . (x - at)),  cos(k . (x - at)) }
!
! whose span holds every wave of that wavenumber, whatever its
! phase. A fit over this form differentiates such waves exactly,
! where a polynomial of any finite degree only approximates them.
!
!=====================================================================!

  public :: harmonic_form

  type, extends(form) :: harmonic_form

     real(dp) :: wavenumber(3) = [1.0_dp, 0.0_dp, 0.0_dp]

   contains

     procedure :: size_of => harmonic_size
     procedure :: values  => harmonic_values
     procedure :: slopes  => harmonic_slopes

  end type harmonic_form

  interface harmonic_form
     module procedure create_harmonic
  end interface harmonic_form

contains

  !===================================================================!
  ! A concretion declares its basis once, standing every entry of its
  ! table. The identity is minted here so no concretion has to
  ! remember to; the roster starts full because an unrestricted form
  ! stands whole.
  !===================================================================!

  subroutine declare_basis(this, width)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: width

    integer :: m

    call this % basis % declare()
    this % active = listed_set_representation([(m, m = 1, width)])

  end subroutine declare_basis

  !===================================================================!
  ! WHICH basis this form's standing members belong to. The identity
  ! survives restriction: restricting a form narrows who stands, and
  ! does not make it a different basis.
  !===================================================================!

  type(set_graph) function basis_set(this) result(b)

    class(form), intent(in) :: this

    b = this % basis

  end function basis_set

  !===================================================================!
  ! Who stands, in declaration order.
  !===================================================================!

  pure subroutine members(this, standing)

    class(form)         , intent(in)  :: this
    integer, allocatable, intent(out) :: standing(:)

    call this % active % members(standing)

  end subroutine members

  !===================================================================!
  ! Stand only these table entries. The kept indices name entries of
  ! the concretion's own table, and the roster is the whole statement
  ! of who belongs.
  !
  ! It MUTATES, as it always has. Making restriction functional - a
  ! new form, a new basis identity - is a separate transformation and
  ! is not smuggled in behind a type change.
  !===================================================================!

  subroutine restrict(this, kept)

    class(form), intent(inout) :: this
    integer    , intent(in)    :: kept(:)

    this % active = listed_set_representation(kept)

  end subroutine restrict


  ! Born with every table entry standing: the members are the four.
  type(polynomial_form) function create_polynomial() result(this)

    call this % declare_basis(4)

  end function create_polynomial

  pure integer function polynomial_size(this)

    class(polynomial_form), intent(in) :: this

    associate (u1 => this); end associate

    polynomial_size = 4

  end function polynomial_size

  pure subroutine polynomial_values(this, x, at, phi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3)
    real(dp), intent(out) :: phi(:)

    associate (u1 => this); end associate

    phi(1)   = 1.0_dp
    phi(2:4) = x - at

  end subroutine polynomial_values

  pure subroutine polynomial_slopes(this, x, at, direction, dphi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3), direction(3)
    real(dp), intent(out) :: dphi(:)

    associate (u1 => this, u2 => x, u3 => at); end associate

    dphi(1)   = 0.0_dp
    dphi(2:4) = direction

  end subroutine polynomial_slopes



  ! Born with every table entry standing: the members are the three.
  type(harmonic_form) function create_harmonic(wavenumber) result(this)

    real(dp), intent(in) :: wavenumber(3)

    integer :: m

    this % wavenumber = wavenumber
    call this % declare_basis(3)

  end function create_harmonic

  pure integer function harmonic_size(this)

    class(harmonic_form), intent(in) :: this

    associate (u1 => this); end associate

    harmonic_size = 3

  end function harmonic_size

  pure subroutine harmonic_values(this, x, at, phi)

    class(harmonic_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3)
    real(dp), intent(out) :: phi(:)

    real(dp) :: phase

    phase = dot_product(this % wavenumber, x - at)

    phi(1) = 1.0_dp
    phi(2) = sin(phase)
    phi(3) = cos(phase)

  end subroutine harmonic_values

  !===================================================================!
  ! d/dn of a wave: the chain rule brings down k . n.
  !===================================================================!

  pure subroutine harmonic_slopes(this, x, at, direction, dphi)

    class(harmonic_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3), direction(3)
    real(dp), intent(out) :: dphi(:)

    real(dp) :: phase, kn

    phase = dot_product(this % wavenumber, x - at)
    kn    = dot_product(this % wavenumber, direction)

    dphi(1) = 0.0_dp
    dphi(2) =  kn * cos(phase)
    dphi(3) = -kn * sin(phase)

  end subroutine harmonic_slopes


end module graph_forms
