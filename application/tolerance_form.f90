! The analytical form of the weight a block carries, identified from
! the family's own coefficients rather than measured from its matrix.
!
! A family's coefficients are dimensionless; the weight an edge
! finally carries is the coefficient times dt raised to the difference
! of the two degrees it joins. For a difference family every source is
! a value and a row determining the d-th derivative therefore carries
! dt^-d.
!
! The d-th row is the velocity row composed d times, so its
! coefficients are bounded by the d-th power of the velocity row's,
! and this reports whether that bound is attained. If it is, then
!
!      ||A||_inf  =  1  +  (sum_j |a_j|)^d / dt^d
!
! with a the velocity coefficients and the one being the unit a row
! carries on the column it determines. That is a closed form for the
! norm, hence for the floor eps ||A|| ||q|| a march can reach, and
! hence for a tolerance that need not be chosen.
!
! Both sides are computed here: the left from the assembled jacobian,
! the right from the family through the same coupling the scheme is
! built on.
program tolerance_form

  use util_precision  , only : dp
  use operation_coupling    , only : weights_of
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : field
  use field_stored          , only : stored_field
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_grid        , only : uniform_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use gti_march             , only : partition
  use gti_sweeps            , only : jacobian_of
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use gti_chain             , only : one_functional
  use gti_chain             , only : chain_block, march_chain, &
       & chain_system, chain_systems

  implicit none

  write(*,'(a)') ' '
  write(*,'(a)') '  the velocity row, and the d-th row it composes to'
  write(*,'(a)') '  scheme    d    sum|a|    sum|c(d)|   (sum|a|)^d   attained'
  call composed(1, 2)
  call composed(2, 2)
  call composed(3, 2)
  call composed(4, 2)
  call composed(2, 3)
  call composed(3, 3)

  write(*,'(a)') ' '
  write(*,'(a)') '  ||A||_inf, from the family against the assembled jacobian'
  write(*,'(a)') '  scheme    d   instants       dt      from family     assembled   agreement'
  call against(1, 3, 21)
  call against(2, 3, 21)
  call against(3, 3, 21)
  call against(2, 3, 61)
  call against(2, 4, 61)
  call against(3, 3, 41)

  write(*,'(a)') ' '
  write(*,'(a)') '  conditioning, and what the same form removes from it'
  write(*,'(a)') '  the row determining degree d is divided by dt^-d, which is the'
  write(*,'(a)') '  weight the family says it carries'
  write(*,'(a)') ' '
  write(*,'(a)') '  scheme    d   instants       dt      kappa(A)   kappa(DA)     ratio'
  call conditioned(2, 3, 21)
  call conditioned(2, 3, 41)
  call conditioned(2, 3, 61)
  call conditioned(2, 3, 81)
  call conditioned(3, 3, 41)
  call conditioned(2, 4, 41)

contains

  !-------------------------------------------------------------------!
  ! The absolute sum of a family's coefficients on the row determining
  ! the given degree, from the coupling that row reads.
  !-------------------------------------------------------------------!

  real(dp) function row_sum(scheme, order, determines) result(total)

    class(family), intent(in) :: scheme
    integer      , intent(in) :: order, determines

    real(dp), allocatable :: c(:)
    integer :: reach, last, k

    reach = determines * order
    last  = reach + 1

    call weights_of(scheme, last, [(last - k, k = 0, reach)], [(last, k = 0, reach)], &
         & [(1.0_dp, k = 1, last)], [(0, k = 0, reach)], [(determines, k = 0, reach)], c)

    total = sum(abs(c))

  end function row_sum

  subroutine composed(order, determines)

    integer, intent(in) :: order, determines

    real(dp) :: velocity, derived, powered
    character(len=8) :: named

    velocity = row_sum(bdf_family(order), order, 1)
    derived  = row_sum(bdf_family(order), order, determines)
    powered  = velocity ** determines

    write(named,'(a,i0)') 'bdf ', order
    write(*,'(a,a,i5,3f12.4,a)') '  ', named, determines, velocity, derived, powered, &
         & merge('   yes', '    no', abs(derived - powered) <= 1.0e-10_dp * powered)

  end subroutine composed

  !-------------------------------------------------------------------!
  ! The closed form beside the assembled jacobian's infinity norm.
  !-------------------------------------------------------------------!

  subroutine against(order, degrees, instants)

    integer, intent(in) :: order, degrees, instants

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    type(chain_system), allocatable :: systems(:)
    type(bdf_family) :: scheme
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:), a(:,:)
    real(dp) :: achieved, duration, design, predicted, assembled
    character(len=8) :: named
    integer :: k, d, top

    duration = 3.0_dp
    design   = 1.0_dp
    scheme   = bdf_family(order)
    top      = degrees - 1

    allocate(schemes(1))
    allocate(schemes(1) % scheme, source=scheme)
    added = [instants]

    call partition(duration, instants, dt, t)
    held = [((0.0_dp, d = 0, degrees - 1), k = 1, &
         &   scheme % history_depth(degrees - 1))]

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)
    call chain_systems(chain, [one_functional(van_der_pol_energy(degrees - 1))], degrees, &
         & design, systems)

    predicted = 1.0_dp + row_sum(scheme, order, top) / dt(size(dt)) ** top
    call dense_jacobian(chain, degrees, design, a)
    assembled = largest_row(a)

    write(named,'(a,i0)') 'bdf ', order
    write(*,'(a,a,i5,i10,f10.5,2es15.5,f11.3,a)') '  ', named, top, instants, &
         & dt(size(dt)), predicted, assembled, &
         & 100.0_dp * (1.0_dp - abs(predicted - assembled) / assembled), ' %'

  end subroutine against

  !-------------------------------------------------------------------!
  ! kappa in the infinity norm, before and after every row is divided
  ! by the weight its own family gives it.
  !-------------------------------------------------------------------!

  subroutine conditioned(order, degrees, instants)

    integer, intent(in) :: order, degrees, instants

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    type(chain_system), allocatable :: systems(:)
    type(bdf_family) :: scheme
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:), b(:,:), a(:,:)
    real(dp) :: achieved, duration, design, bare, scaled, step
    character(len=8) :: named
    integer :: k, d, i, n

    duration = 3.0_dp
    design   = 1.0_dp
    scheme   = bdf_family(order)

    allocate(schemes(1))
    allocate(schemes(1) % scheme, source=scheme)
    added = [instants]

    call partition(duration, instants, dt, t)
    held = [((0.0_dp, d = 0, degrees - 1), k = 1, &
         &   scheme % history_depth(degrees - 1))]

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)
    call chain_systems(chain, [one_functional(van_der_pol_energy(degrees - 1))], degrees, &
         & design, systems)

    step = dt(size(dt))
    call dense_jacobian(chain, degrees, design, a)
    bare = kappa(a)

    n = size(a, 1)
    allocate(b(n, n))
    do i = 1, n
       d = mod(i - 1, degrees)
       b(i, :) = a(i, :) * step ** d
    end do
    scaled = kappa(b)

    write(named,'(a,i0)') 'bdf ', order
    write(*,'(a,a,i5,i10,f10.5,2es13.4,f10.1)') '  ', named, degrees - 1, &
         & instants, step, bare, scaled, bare / scaled

  end subroutine conditioned

  !-------------------------------------------------------------------!
  ! ||A||_inf ||A^-1||_inf, the inverse by elimination on the identity.
  !-------------------------------------------------------------------!

  real(dp) function kappa(a) result(k)

    real(dp), intent(in) :: a(:,:)

    real(dp), allocatable :: w(:,:), inverse(:,:), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, j, p

    n = size(a, 1)
    allocate(w(n, n), inverse(n, n), row(n))

    w       = a
    inverse = 0.0_dp
    do i = 1, n
       inverse(i, i) = 1.0_dp
    end do

    do j = 1, n
       p = j - 1 + maxloc(abs(w(j:n, j)), dim=1)
       if (p /= j) then
          row          = w(j, :)
          w(j, :)      = w(p, :)
          w(p, :)      = row
          row          = inverse(j, :)
          inverse(j, :) = inverse(p, :)
          inverse(p, :) = row
       end if
       pivot = w(j, j)
       if (abs(pivot) <= tiny(1.0_dp)) then
          k = huge(1.0_dp)
          return
       end if
       w(j, :)       = w(j, :) / pivot
       inverse(j, :) = inverse(j, :) / pivot
       do i = 1, n
          if (i == j) cycle
          factor = w(i, j)
          w(i, :)       = w(i, :) - factor * w(j, :)
          inverse(i, :) = inverse(i, :) - factor * inverse(j, :)
       end do
    end do

    k = largest_row(a) * largest_row(inverse)

  end function kappa

  real(dp) function largest_row(a) result(most)

    real(dp), intent(in) :: a(:,:)

    most = maxval(sum(abs(a), dim=2))

  end function largest_row

  !-------------------------------------------------------------------!
  ! The jacobian of the first block at its solved state, formed from
  ! the compiled tangent, for the measurements below.
  !-------------------------------------------------------------------!

  subroutine dense_jacobian(chain, degrees, design, a)

    type(chain_block), intent(in) :: chain(:)
    integer          , intent(in) :: degrees
    real(dp)         , intent(in) :: design
    real(dp), allocatable, intent(out) :: a(:,:)

    type(stored_directed_graph) :: unknowns
    type(stored_field) :: state, knobs
    integer :: n

    n        = chain(1) % rows % num_unknowns()
    unknowns = stored_directed_graph(n, tails=[integer ::], heads=[integer ::])
    state    = stored_field('state', unknowns % vertex_set(), n)
    knobs    = stored_field('nu', unknowns % vertex_set(), chain(1) % rows % num_points())
    call state % set_real_vector(chain(1) % state)
    call knobs % set_real_vector(spread(design, 1, chain(1) % rows % num_points()))
    call jacobian_of(chain(1) % rows, unknowns, [state, knobs], n, unknowns % vertex_set(), a)

    associate (u1 => degrees); end associate

  end subroutine dense_jacobian

end program tolerance_form
