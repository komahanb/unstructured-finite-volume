! How much a functional of the trajectory moves when the design does,
! computed three ways that share no arithmetic.
!
! The block is van der Pol at a design of one, so the statement is
! nonlinear and newton has to iterate. The functional is the energy
! summed over the instants with the step that ends at each.
!
!    tangent     one solve in the state, then the functional's
!                gradient read along what it gives
!    adjoint     one solve against the transpose, then the
!                statement's design partial read along what it gives
!    difference  the whole march repeated either side of the design
!                and the functional differenced
!
! The first two read the same three objects through different
! algebra and must agree to round-off. The third repeats everything
! and agrees only if all of it is right, so it is the one that could
! fail on its own.
program sensitivity

  use iso_fortran_env      , only : dp => REAL64
  use view_directed_stored , only : stored_directed_graph
  use field_stored         , only : stored_field
  use operation_family_bdf , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family     , only : family
  use physics_vanderpol    , only : van_der_pol, van_der_pol_energy
  use gti_block            , only : block_residual
  use gti_march            , only : partition, block_of, solved, unknowns_graph
  use gti_sweeps           , only : functional_of, functional_gradient, &
       & design_partial, jacobian_of, by_tangent, by_adjoint

  implicit none

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  integer , parameter :: num_instants = 21
  real(dp), parameter :: duration = 2.0_dp
  real(dp), parameter :: design = 1.0_dp

  call sensitivity_of('bdf 2', bdf_family(2))
  call sensitivity_of('adams-moulton 3', adams_family(3))

contains

  !-------------------------------------------------------------------!
  ! The d-th derivative of the cosine, which is what the carried
  ! instants hold: they are the block's initial data and any
  ! consistent numbers would do.
  !-------------------------------------------------------------------!

  pure real(dp) function exact(d, t) result(q)

    integer , intent(in) :: d
    real(dp), intent(in) :: t

    select case (mod(d, 4))
    case (0)
       q =  cos(t)
    case (1)
       q = -sin(t)
    case (2)
       q = -cos(t)
    case default
       q =  sin(t)
    end select

  end function exact

  function carried_values(scheme, t) result(held)

    class(family), intent(in) :: scheme
    real(dp)     , intent(in) :: t(:)
    real(dp), allocatable :: held(:)

    integer :: k, d

    held = [((exact(d, t(k)), d = 0, degrees - 1), k = 1, scheme % history_depth())]

  end function carried_values

  !-------------------------------------------------------------------!
  ! The functional at one design, which is one whole march.
  !-------------------------------------------------------------------!

  real(dp) function marched(scheme, design_value, q) result(f)

    class(family), intent(in) :: scheme
    real(dp)     , intent(in) :: design_value
    real(dp), allocatable, intent(out) :: q(:)

    type(block_residual) :: rows
    type(stored_directed_graph) :: unknowns, instants
    type(stored_field) :: state, knobs
    real(dp), allocatable :: dt(:), t(:)
    real(dp) :: achieved

    call partition(duration, num_instants, dt, t)
    rows = block_of(scheme, van_der_pol(state_degree), degrees, num_instants, dt, &
         & carried_values(scheme, t))
    call solved(rows, num_instants, degrees, design_value, q, achieved)

    unknowns = unknowns_graph(num_instants, degrees)
    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), size(q))
    knobs    = stored_field('nu', unknowns % vertex_set(), num_instants)
    call state % set_real_vector(q)
    call knobs % set_real_vector(spread(design_value, 1, num_instants))

    f = functional_of(van_der_pol_energy(state_degree), instants, [state, knobs], dt)

  end function marched

  subroutine sensitivity_of(title, scheme)

    character(len=*), intent(in) :: title
    class(family)   , intent(in) :: scheme

    real(dp), parameter :: delta = 1.0e-6_dp

    real(dp), allocatable :: q(:), g(:), rate(:), a(:,:), plus(:), minus(:)
    real(dp) :: f, tangent, adjoint, differenced

    f = marched(scheme, design, q)
    call three_objects(scheme, q, g, rate, a)

    tangent = by_tangent(a, g, rate, 0.0_dp)
    adjoint = by_adjoint(a, g, rate, 0.0_dp)

    differenced = (marched(scheme, design + delta, plus) - &
         &         marched(scheme, design - delta, minus)) / (2.0_dp * delta)

    write(*,'(a)')        ' '
    write(*,'(a)')        ' ' // title // ', van der pol at a design of one'
    write(*,'(a,f16.10)') '   the functional             ', f
    write(*,'(a,f16.10)') '   sensitivity, tangent       ', tangent
    write(*,'(a,f16.10)') '   sensitivity, adjoint       ', adjoint
    write(*,'(a,f16.10)') '   sensitivity, differenced   ', differenced
    write(*,'(a,es16.2)') '   tangent against adjoint    ', abs(tangent - adjoint)
    write(*,'(a,es16.2)') '   tangent against difference ', abs(tangent - differenced)

  end subroutine sensitivity_of

  !-------------------------------------------------------------------!
  ! The three objects the tangent and the adjoint both read: the
  ! functional's gradient in the state, the statement's partial in
  ! the design, and the jacobian.
  !-------------------------------------------------------------------!

  subroutine three_objects(scheme, q, g, rate, a)

    class(family), intent(in) :: scheme
    real(dp)     , intent(in) :: q(:)
    real(dp), allocatable, intent(out) :: g(:), rate(:), a(:,:)

    type(block_residual) :: rows
    type(stored_directed_graph) :: unknowns, instants
    type(stored_field) :: state, knobs
    real(dp), allocatable :: dt(:), t(:)

    call partition(duration, num_instants, dt, t)
    rows = block_of(scheme, van_der_pol(state_degree), degrees, num_instants, dt, &
         & carried_values(scheme, t))

    unknowns = unknowns_graph(num_instants, degrees)
    instants = stored_directed_graph(num_instants, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), size(q))
    knobs    = stored_field('nu', unknowns % vertex_set(), num_instants)
    call state % set_real_vector(q)
    call knobs % set_real_vector(spread(design, 1, num_instants))

    call functional_gradient(van_der_pol_energy(state_degree), instants, &
         & [state, knobs], dt, num_instants, degrees, unknowns % vertex_set(), g)
    call design_partial(rows, unknowns, [state, knobs], num_instants, &
         & unknowns % vertex_set(), rate)
    call jacobian_of(rows, unknowns, [state, knobs], num_instants * degrees, &
         & unknowns % vertex_set(), a)

  end subroutine three_objects

end program sensitivity
