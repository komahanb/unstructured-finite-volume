!=====================================================================!
! The robin condition: one tagged boundary statement.
!
! LEVEL 3 OF THE STRATIFICATION. The condition says what the
! material does at a boundary, and it says it as COEFFICIENTS - this
! class computes numbers and hands them to the calculus; it owns no
! operator, no balance, no solve.
!
! Every boundary condition is one statement,
!
!      a*phi + b*dphi/dn = c        on the faces of one tag
!
! dirichlet is a = 1, b = 0; neumann is a = 0, b = 1; anything mixed
! is itself. The face value is eliminated one-sidedly,
!
!      phi_b = (c + (b/delta)*phi_p) / (a + b/delta)
!
! and the face flux splits into a phi_p coefficient and a constant.
! With denom = a + b/delta, the four numbers per face (matching
! class_boundary_condition of the old world, checked to machine
! precision in the suite):
!
!      lhs     = -kappa*area*a/(delta*denom)      multiplies phi_p
!      rhs     = -kappa*area*c/(delta*denom)      the constant
!      adv lhs = -vn*area*(b/delta)/denom         multiplies phi_p
!      adv rhs =  vn*area*c/denom                 the constant
!
! THE STRING ENTERS ONCE. The condition holds its tag, and resolves
! it through tagged_edges at the moment coefficients are asked for.
! Nothing downstream holds the string; everything downstream holds
! arrays.
!
! THE OPERATOR PATH, for a > 0: the diffusive part uses the
! calculus directly. An edge coefficient kappa*area*a/denom with
! spacing delta and the stored value c/a in the operator's boundary
! argument reproduces the eliminated flux with the old row's own
! sign, so a boundary strengthens the diagonal exactly as the old
! assembler's lhs does - the suite demonstrates the row. For a = 0
! the constant travels as a source in the balance instead, because
! a zero coefficient carries nothing.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_robin_condition

  use util_precision  , only : dp
  use view_directed, only : directed_graph
  use graph_fractal      , only : graph
  use map_set_store, only : set_store
  use field_stored, only : stored_field
  use view_mesh, only : mesh

  implicit none

  private
  integer, parameter :: COEFFICIENT_LHS = 1
  integer, parameter :: COEFFICIENT_RHS = 2
  integer, parameter :: COEFFICIENT_ADVECTION_LHS = 3
  integer, parameter :: COEFFICIENT_ADVECTION_RHS = 4
  integer, parameter :: COEFFICIENT_OPERATOR = 5

  public :: robin_condition
  public :: robin, dirichlet, neumann

  !===================================================================!
  ! One condition: a tag and three numbers.
  !===================================================================!

  type :: robin_condition

     character(len=:), allocatable :: tag

     real(dp) :: a = 1.0_dp
     real(dp) :: b = 0.0_dp
     real(dp) :: c = 0.0_dp

   contains

     procedure :: faces
     procedure :: lhs_coefficients
     procedure :: rhs_coefficients
     procedure :: advection_lhs_coefficients
     procedure :: advection_rhs_coefficients
     procedure :: operator_coefficients
     procedure :: boundary_values
     procedure :: wall_relation

  end type robin_condition

contains

  !===================================================================!
  ! The three constructors: the general statement and its two
  ! famous specializations.
  !===================================================================!

  pure type(robin_condition) function robin(tag, a, b, c) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: a, b, c

    this % tag = tag
    this % a   = a
    this % b   = b
    this % c   = c

  end function robin

  pure type(robin_condition) function dirichlet(tag, value) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: value

    this = robin(tag, 1.0_dp, 0.0_dp, value)

  end function dirichlet

  pure type(robin_condition) function neumann(tag, flux) result(this)

    character(len=*), intent(in) :: tag
    real(dp)        , intent(in) :: flux

    this = robin(tag, 0.0_dp, 1.0_dp, flux)

  end function neumann

  !===================================================================!
  ! The tag resolved, once: WHICH edges this condition covers.
  !
  ! The answer is a declared subset, so it is a new declared domain, and
  ! the caller's set store says who belongs, what it is called and
  ! which ambient set it sits under. It is an argument because the answer
  ! outlives this call - a set the caller cannot interpret would be no
  ! answer at all.
  !===================================================================!

  subroutine faces(this, m, sets, members)

    class(robin_condition), intent(in)    :: this
    type(mesh)            , intent(in)    :: m
    type(set_store)       , intent(inout) :: sets
    type(graph)       , intent(out)   :: members

    call m % tagged_edges(this % tag, sets, members)

  end subroutine faces

  !===================================================================!
  ! The four coefficient arrays, one entry per tagged face, in the
  ! member set's order. Each reads the mesh's own area and delta at
  ! its face.
  !===================================================================!

  subroutine lhs_coefficients(this, m, kappa, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: kappa
    real(dp), allocatable, intent(out) :: values(:)

    call coefficient_values(this, m, kappa, COEFFICIENT_LHS, values)

  end subroutine lhs_coefficients

  subroutine rhs_coefficients(this, m, kappa, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: kappa
    real(dp), allocatable, intent(out) :: values(:)

    call coefficient_values(this, m, kappa, COEFFICIENT_RHS, values)

  end subroutine rhs_coefficients

  subroutine advection_lhs_coefficients(this, m, vn, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: vn
    real(dp), allocatable, intent(out) :: values(:)

    call coefficient_values(this, m, vn, COEFFICIENT_ADVECTION_LHS, values)

  end subroutine advection_lhs_coefficients

  subroutine advection_rhs_coefficients(this, m, vn, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: vn
    real(dp), allocatable, intent(out) :: values(:)

    call coefficient_values(this, m, vn, COEFFICIENT_ADVECTION_RHS, values)

  end subroutine advection_rhs_coefficients

  !===================================================================!
  ! The operator road, for a > 0. The edge coefficient and the
  ! stored value that make the calculus reproduce the eliminated
  ! flux on a headless edge:
  !
  !      z = c_f*(value - phi_p)/delta  =  lhs*phi_p - rhs
  !
  !      c_f   = kappa*area*a/denom  =  -lhs*delta
  !      value = c/a
  !===================================================================!

  subroutine operator_coefficients(this, m, kappa, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: kappa
    real(dp), allocatable, intent(out) :: values(:)

    call coefficient_values(this, m, kappa, COEFFICIENT_OPERATOR, values)

  end subroutine operator_coefficients

  subroutine coefficient_values(this, m, scale, which, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp)  , intent(in)             :: scale
    integer   , intent(in)             :: which
    real(dp), allocatable, intent(out) :: values(:)

    real(dp), allocatable :: area(:), delta(:)
    integer :: f

    call measures_of(this, m, area, delta)
    allocate(values(size(area)))

    select case (which)
    case (COEFFICIENT_LHS)
       do f = 1, size(values)
          values(f) = -scale * area(f) * this % a &
               & / (delta(f) * denom(this, delta(f)))
       end do
    case (COEFFICIENT_RHS)
       do f = 1, size(values)
          values(f) = -scale * area(f) * this % c &
               & / (delta(f) * denom(this, delta(f)))
       end do
    case (COEFFICIENT_ADVECTION_LHS)
       do f = 1, size(values)
          values(f) = -scale * area(f) * (this % b / delta(f)) &
               & / denom(this, delta(f))
       end do
    case (COEFFICIENT_ADVECTION_RHS)
       do f = 1, size(values)
          values(f) = scale * area(f) * this % c / denom(this, delta(f))
       end do
    case (COEFFICIENT_OPERATOR)
       do f = 1, size(values)
          values(f) = scale * area(f) * this % a / denom(this, delta(f))
       end do
    case default
       error stop 'operation_robin_condition: unknown coefficient projection'
    end select

  end subroutine coefficient_values

  subroutine boundary_values(this, m, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: values(:)

    real(dp), allocatable :: area(:), delta(:)

    call measures_of(this, m, area, delta)
    allocate(values(size(area)))

    if (abs(this % a) > 0.0_dp) then
       values = this % c / this % a
    else
       ! A pure neumann condition has no value to stand in; its
       ! constant travels as a source in the balance instead.
       values = 0.0_dp
    end if

  end subroutine boundary_values

  !===================================================================!
  ! THE WHOLE WALL, IN TWO NUMBERS. Eliminating the face value from
  !
  !      a*phi_b + b*(phi_b - phi_p)/delta = c
  !
  ! gives, with denom = a + b/delta,
  !
  !      phi_b = (1 - w)*phi_p + v        w = a/denom, v = c/denom
  !
  ! and that one line is the entire condition: dirichlet is w = 1,
  ! phi_b = c; neumann is w = 0, phi_b = phi_p + c*delta, the wall
  ! that carries a gradient rather than a value; anything mixed
  ! stands between them. A caller that takes only v can express
  ! dirichlet and nothing else, which is why both numbers travel.
  !===================================================================!

  subroutine wall_relation(this, m, weights, values)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: weights(:)
    real(dp), allocatable, intent(out) :: values(:)

    real(dp), allocatable :: area(:), delta(:)
    integer :: f

    call measures_of(this, m, area, delta)
    allocate(weights(size(area)), values(size(area)))

    do f = 1, size(area)
       weights(f) = this % a / denom(this, delta(f))
       values(f)  = this % c / denom(this, delta(f))
    end do

  end subroutine wall_relation

  !===================================================================!
  ! The mesh's area and delta at this condition's faces, in member
  ! order. The one place the tag is resolved.
  !===================================================================!

  subroutine measures_of(this, m, area, delta)

    class(robin_condition), intent(in) :: this
    type(mesh), intent(in)             :: m
    real(dp), allocatable, intent(out) :: area(:)
    real(dp), allocatable, intent(out) :: delta(:)

    !----------------------------------------------------------------!
    ! The declared set is local to this call, so its interpretation is
    ! local too: this store is not a hidden environment. Nothing that
    ! needs it escapes.
    !----------------------------------------------------------------!

    type(graph)     :: members
    type(set_store) :: sets

    type(stored_field) :: fa, fd
    real(dp), allocatable :: all_areas(:), all_deltas(:)
    integer :: f, e

    call m % tagged_edges(this % tag, sets, members)

    fa = m % face_area()
    call fa % real_vector(all_areas)
    fd = m % face_delta()
    call fd % real_vector(all_deltas)

    allocate(area(sets % num_members_of(members)))
    allocate(delta(sets % num_members_of(members)))

    do f = 1, size(area)
       e = sets % member_of(members, f)
       area(f)  = all_areas(e)
       delta(f) = all_deltas(e)
    end do

  end subroutine measures_of

  !===================================================================!
  ! The shared denominator of every formula: a + b/delta.
  !===================================================================!

  pure real(dp) function denom(this, delta)

    class(robin_condition), intent(in) :: this
    real(dp)              , intent(in) :: delta

    denom = this % a + this % b / delta

  end function denom

end module operation_robin_condition
