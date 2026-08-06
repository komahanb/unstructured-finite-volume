!=====================================================================!
! The harmonic form: one wave and the constant,
!
!      { 1,  sin(k . (x - at)),  cos(k . (x - at)) }
!
! whose span holds every wave of that wavenumber, whatever its
! phase. A fit over this form differentiates such waves exactly,
! where a polynomial of any finite degree only approximates them.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_harmonic_form

  use iso_fortran_env    , only : dp => REAL64
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use graph_forms        , only : form

  implicit none

  private
  public :: harmonic_form

  type, extends(form) :: harmonic_form

     real(dp) :: wavenumber(3) = [1.0_dp, 0.0_dp, 0.0_dp]

   contains

     procedure :: size_of => harmonic_size
     procedure :: values  => harmonic_values
     procedure :: slopes  => harmonic_slopes

  end type harmonic_form

  interface harmonic_form
     module procedure create
  end interface harmonic_form

contains

  ! Born with every table entry standing: the members are the three.
  pure type(harmonic_form) function create(wavenumber) result(this)

    real(dp), intent(in) :: wavenumber(3)

    integer :: m

    this % wavenumber = wavenumber
    this % support = support(GRAPH_SIDE_VERTEX, [(m, m = 1, 3)])

  end function create

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

end module class_harmonic_form
