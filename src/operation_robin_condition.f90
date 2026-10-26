!=====================================================================!
! Module that states the robin condition: one tagged boundary
! statement.
!
! This is level 3 of the stratification. The condition specifies the
! boundary behaviour of the material as coefficients: this type
! computes numbers and passes them to the calculus, and it owns no
! operator, no balance and no solve.
!
! Every boundary condition is one statement,
!
!      a*phi + b*dphi/dn = c        on the faces of one tag
!
! Dirichlet is a = 1, b = 0; Neumann is a = 0, b = 1; anything mixed
! is itself. The face value is eliminated one-sidedly,
!
!      phi_b = (c + (b/delta)*phi_p) / (a + b/delta)
!
!   cell centre p                    boundary face b
!        o------------ delta ------------|
!      phi_p                           phi_b
!
! and the face flux splits into a phi_p coefficient and a constant.
! With denom = a + b/delta, every number per face is one formula,
!
!      sign * scale * area * (coefficient / delta^q) / (delta^p * denom)
!
! and the five coefficients are its table (matching
! class_boundary_condition of the previous implementation, checked to
! machine precision in the test suite):
!
!      lhs     = -kappa*area*a/(delta*denom)      multiplies phi_p
!      rhs     = -kappa*area*c/(delta*denom)      the constant
!      adv lhs = -vn*area*(b/delta)/denom         multiplies phi_p
!      adv rhs =  vn*area*c/denom                 the constant
!      operator = kappa*area*a/denom              the edge coefficient
!
! The string enters once. The condition stores its tag, and resolves
! it through tagged_edges when coefficients are requested. No later
! stage stores the string; every later stage stores arrays.
!
! The operator path, for a > 0: the diffusive part uses the calculus
! directly. An edge coefficient kappa*area*a/denom with spacing delta
! and the stored value c/a in the operator's boundary argument
! reproduces the eliminated flux with the sign of the previous row, so
! a boundary increases the diagonal exactly as the previous
! assembler's lhs does; the test suite checks the row. For a = 0 the
! constant enters as a source in the balance instead, because a zero
! coefficient contributes nothing.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_robin_condition

  use util_precision, only : dp
  use view_directed , only : directed_graph
  use graph_fractal , only : graph
  use map_set_store , only : set_store
  use view_mesh     , only : mesh, values_of

  implicit none

  private

  ! The five coefficients of the table
  integer, parameter, public :: COEFFICIENT_LHS           = 1
  integer, parameter, public :: COEFFICIENT_RHS           = 2
  integer, parameter, public :: COEFFICIENT_ADVECTION_LHS = 3
  integer, parameter, public :: COEFFICIENT_ADVECTION_RHS = 4
  integer, parameter, public :: COEFFICIENT_OPERATOR      = 5

  public :: robin_condition
  public :: robin, dirichlet, neumann

  !-------------------------------------------------------------------!
  ! Type that models one condition: a tag and three numbers
  !-------------------------------------------------------------------!

  type :: robin_condition

     character(len=:), allocatable :: tag

     real(dp) :: a = 1.0_dp
     real(dp) :: b = 0.0_dp
     real(dp) :: c = 0.0_dp

   contains

     procedure :: faces
     procedure :: coefficient_values
     procedure :: boundary_values
     procedure :: boundary_relation

  end type robin_condition

contains

  !===================================================================!
  ! Creates the general statement, a*phi + b*dphi/dn = c on the faces
  ! of a tag. Dirichlet and Neumann below are its two standard
  ! specializations.
  !===================================================================!

  pure type(robin_condition) function robin(tag, a, b, c) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: a, b, c

    this % tag = tag
    this % a   = a
    this % b   = b
    this % c   = c

  end function robin

  !===================================================================!
  ! Creates the Dirichlet condition, phi = value
  !===================================================================!

  pure type(robin_condition) function dirichlet(tag, value) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: value

    this = robin(tag, 1.0_dp, 0.0_dp, value)

  end function dirichlet

  !===================================================================!
  ! Creates the Neumann condition, dphi/dn = flux
  !===================================================================!

  pure type(robin_condition) function neumann(tag, flux) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: flux

    this = robin(tag, 0.0_dp, 1.0_dp, flux)

  end function neumann

  !===================================================================!
  ! Resolves the tag, once: which edges this condition covers.
  !
  ! The result is a declared subset, so it is a new declared domain,
  ! and the caller's set store records its members, its label and its
  ! ambient set. The set store is an argument because the result
  ! outlives this call; a set the caller cannot interpret would be
  ! unusable.
  !===================================================================!

  subroutine faces(this, m, sets, members)

    class(robin_condition), intent(in)    :: this
    type(mesh)            , intent(in)    :: m
    type(set_store)       , intent(inout) :: sets
    type(graph)           , intent(out)   :: members

    call m % tagged_edges(this % tag, sets, members)

  end subroutine faces

  !===================================================================!
  ! Returns one of the five coefficients, one entry per tagged face, in
  ! the member set's order: the formula's table, selected by which.
  ! Invalid input: a selector outside the five.
  !===================================================================!

  subroutine coefficient_values(this, m, scale, which, values)

    class(robin_condition), intent(in)  :: this
    type(mesh)            , intent(in)  :: m
    real(dp)              , intent(in)  :: scale
    integer               , intent(in)  :: which
    real(dp), allocatable , intent(out) :: values(:)

    type(graph)           :: members
    type(set_store)       :: sets
    real(dp), allocatable :: area(:), delta(:)
    real(dp)              :: sign_factor, coefficient
    integer               :: q, p
    character(len=250)    :: message

    ! The row of the table
    select case (which)
    case (COEFFICIENT_LHS)
       sign_factor = -1.0_dp; coefficient = this % a; q = 0; p = 1
    case (COEFFICIENT_RHS)
       sign_factor = -1.0_dp; coefficient = this % c; q = 0; p = 1
    case (COEFFICIENT_ADVECTION_LHS)
       sign_factor = -1.0_dp; coefficient = this % b; q = 1; p = 0
    case (COEFFICIENT_ADVECTION_RHS)
       sign_factor =  1.0_dp; coefficient = this % c; q = 0; p = 0
    case (COEFFICIENT_OPERATOR)
       sign_factor =  1.0_dp; coefficient = this % a; q = 0; p = 0
    case default
       write(message,'(a,i0)') 'operation_robin_condition: which does not name one of the five &
            &coefficient projections; which = ', which
       error stop trim(message)
    end select

    call this % faces(m, sets, members)
    call measures_at(m, sets, members, delta, area)
    call face_formula(this, delta, coefficient, q, p, values, prefactor=sign_factor * scale * area)

  end subroutine coefficient_values

  !===================================================================!
  ! Returns the stored value of the operator path, c/a, one entry per
  ! tagged face. A pure Neumann condition has no substitute value; its
  ! constant enters as a source in the balance instead.
  !===================================================================!

  subroutine boundary_values(this, m, values)

    class(robin_condition), intent(in)  :: this
    type(mesh)            , intent(in)  :: m
    real(dp), allocatable , intent(out) :: values(:)

    type(graph)     :: members
    type(set_store) :: sets

    call this % faces(m, sets, members)
    allocate(values(sets % num_members_of(members)))

    if (abs(this % a) > 0.0_dp) then
       values = this % c / this % a
    else
       values = 0.0_dp
    end if

  end subroutine boundary_values

  !===================================================================!
  ! Returns the whole boundary condition in two numbers per face.
  ! Eliminating the face value from
  !
  !      a*phi_b + b*(phi_b - phi_p)/delta = c
  !
  ! gives, with denom = a + b/delta,
  !
  !      phi_b = (1 - w)*phi_p + v        w = a/denom, v = c/denom
  !
  ! and that one line is the entire condition:
  !
  !   w = 1                    0 < w < 1          w = 0
  !   Dirichlet, phi_b = c     mixed              Neumann, phi_b = phi_p + c*delta
  !
  ! Neumann is the boundary that specifies a gradient rather than a
  ! value. A caller that takes only v can express Dirichlet and nothing
  ! else, which is why both numbers are returned. The members are the
  ! caller's, resolved once through faces.
  !===================================================================!

  subroutine boundary_relation(this, m, sets, members, weights, values)

    class(robin_condition), intent(in)  :: this
    type(mesh)            , intent(in)  :: m
    type(set_store)       , intent(in)  :: sets
    type(graph)           , intent(in)  :: members
    real(dp), allocatable , intent(out) :: weights(:)
    real(dp), allocatable , intent(out) :: values(:)

    real(dp), allocatable :: delta(:)

    call measures_at(m, sets, members, delta)
    call face_formula(this, delta, this % a, 0, 0, weights)
    call face_formula(this, delta, this % c, 0, 0, values)

  end subroutine boundary_relation

  !===================================================================!
  ! Returns the one formula, per face,
  !
  !      prefactor * (coefficient / delta^q) / (delta^p * denom)
  !
  ! with q and p zero or one, and the prefactor one when absent
  !===================================================================!

  pure subroutine face_formula(this, delta, coefficient, q, p, values, prefactor)

    class(robin_condition), intent(in)           :: this
    real(dp)              , intent(in)           :: delta(:)
    real(dp)              , intent(in)           :: coefficient
    integer               , intent(in)           :: q, p
    real(dp), allocatable , intent(out)          :: values(:)
    real(dp)              , intent(in), optional :: prefactor(:)

    real(dp) :: numerator
    integer  :: f

    allocate(values(size(delta)))

    do f = 1, size(delta)
       numerator = coefficient / delta(f) ** q
       if (present(prefactor)) numerator = prefactor(f) * numerator
       values(f) = numerator / (delta(f) ** p * denom(this, delta(f)))
    end do

  end subroutine face_formula

  !===================================================================!
  ! Returns the mesh's delta, and the area when requested, at the
  ! member faces, in member order
  !===================================================================!

  subroutine measures_at(m, sets, members, delta, area)

    type(mesh)           , intent(in)            :: m
    type(set_store)      , intent(in)            :: sets
    type(graph)          , intent(in)            :: members
    real(dp), allocatable, intent(out)           :: delta(:)
    real(dp), allocatable, intent(out), optional :: area(:)

    real(dp), allocatable :: all_areas(:), all_deltas(:)
    integer               :: f, e, n

    n = sets % num_members_of(members)

    call values_of(m % face_delta(), all_deltas)
    allocate(delta(n))
    if (present(area)) then
       call values_of(m % face_area(), all_areas)
       allocate(area(n))
    end if

    do f = 1, n
       e = sets % member_of(members, f)
       delta(f) = all_deltas(e)
       if (present(area)) area(f) = all_areas(e)
    end do

  end subroutine measures_at

  !===================================================================!
  ! Returns the shared denominator of every formula, a + b/delta
  !===================================================================!

  pure real(dp) function denom(this, delta)

    class(robin_condition), intent(in) :: this
    real(dp)              , intent(in) :: delta

    denom = this % a + this % b / delta

  end function denom

end module operation_robin_condition
