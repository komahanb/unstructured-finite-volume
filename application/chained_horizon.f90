! A horizon whose blocks are marched by different families, and whose
! layouts therefore differ.
!
! A multistep block holds one set of components per instant; a stage
! block holds, for each step, its stages and then the instant it
! arrives at. So the junction between them is an index map, and the
! only thing a block tells its successor is where among its
! unknowns its k-th instant sits.
!
! Two checks. A chain of two blocks of the same family must give what
! one block of the whole gives, since the same rows are made in the
! same places. Then a chain that changes family part way through is
! expanded, and every order is checked against a difference of the
! order below - which crosses the junction, because the whole chain
! is remarched either side of the design.
program chained_horizon

  use iso_fortran_env       , only : dp => REAL64
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : crouzeix_two_stage
  use operation_grid        , only : uniform_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder
  use gti_march             , only : partition
  use gti_chain             , only : chain_block, march_chain, chain_expansion

  implicit none

  integer , parameter :: state_degree = 2
  integer , parameter :: degrees = state_degree + 1
  integer , parameter :: max_order = 3
  real(dp), parameter :: duration = 2.0_dp
  real(dp), parameter :: delta = 1.0e-4_dp

  call splitting_changes_nothing()
  call across_families('bdf 2 then crouzeix two-stage', bdf_of(2), dirk_of())
  call across_families('crouzeix two-stage then bdf 2', dirk_of(), bdf_of(2))
  call across_families('adams 3 then crouzeix two-stage', adams_of(3), dirk_of())

contains

  pure real(dp) function initial_at(d, t) result(q)

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

  end function initial_at

  function bdf_of(order) result(held)

    integer, intent(in) :: order
    type(family_holder) :: held

    allocate(held % scheme, source=bdf_family(order))

  end function bdf_of

  function adams_of(order) result(held)

    integer, intent(in) :: order
    type(family_holder) :: held

    allocate(held % scheme, source=adams_family(order))

  end function adams_of

  function dirk_of() result(held)

    type(family_holder) :: held

    allocate(held % scheme, source=crouzeix_two_stage())

  end function dirk_of

  !-------------------------------------------------------------------!
  ! One chain marched and expanded.
  !-------------------------------------------------------------------!

  subroutine expanded(schemes, added, design, f)

    type(family_holder), intent(in) :: schemes(:)
    integer            , intent(in) :: added(:)
    real(dp)           , intent(in) :: design
    real(dp), allocatable, intent(out) :: f(:)

    type(chain_block), allocatable :: chain(:)
    real(dp), allocatable :: dt(:), t(:), held(:)
    real(dp) :: achieved
    integer :: k, d, given

    given = schemes(1) % scheme % history_depth(degrees - 1)
    call partition(duration, sum(added), dt, t)
    held = [((initial_at(d, t(k)), d = 0, degrees - 1), k = 1, given)]

    call march_chain(schemes, added, van_der_pol(state_degree), degrees, &
         & uniform_grid(duration), design, held, chain, dt, t, achieved)

    call chain_expansion(chain, van_der_pol(state_degree), &
         & van_der_pol_energy(state_degree), degrees, dt, design, max_order, f)

  end subroutine expanded

  !-------------------------------------------------------------------!
  ! Splitting a horizon among blocks of one family must change
  ! nothing.
  !-------------------------------------------------------------------!

  subroutine splitting_changes_nothing()

    type(family_holder) :: whole(1), split(2)
    real(dp), allocatable :: f_whole(:), f_split(:)

    whole(1) = bdf_of(2)
    split(1) = bdf_of(2)
    split(2) = bdf_of(2)

    call expanded(whole, [20], 1.0_dp, f_whole)
    call expanded(split, [10, 10], 1.0_dp, f_split)

    write(*,'(a)')        ' bdf 2 over twenty instants, in one block and in two'
    write(*,'(a,4f14.8)') '   whole                    ', f_whole
    write(*,'(a,4f14.8)') '   split                    ', f_split
    write(*,'(a,4es14.2)')'   apart                    ', abs(f_whole - f_split)

  end subroutine splitting_changes_nothing

  !-------------------------------------------------------------------!
  ! A chain that changes family, expanded, every order against a
  ! difference of the order below.
  !-------------------------------------------------------------------!

  subroutine across_families(title, one, two)

    character(len=*)   , intent(in) :: title
    type(family_holder), intent(in) :: one, two

    type(family_holder) :: schemes(2)
    real(dp), allocatable :: f(:), plus(:), minus(:)
    real(dp) :: differenced(max_order)
    integer :: m

    schemes(1) = one
    schemes(2) = two

    call expanded(schemes, [10, 10], 1.0_dp, f)
    call expanded(schemes, [10, 10], 1.0_dp + delta, plus)
    call expanded(schemes, [10, 10], 1.0_dp - delta, minus)

    do m = 1, max_order
       differenced(m) = (plus(m - 1) - minus(m - 1)) / (2.0_dp * delta)
    end do

    write(*,'(a)')            ' '
    write(*,'(a)')            ' ' // title
    write(*,'(a,4i14)')       '   derivative order         ', [(m, m = 0, max_order)]
    write(*,'(a,4f14.8)')     '   from the expansion       ', f
    write(*,'(a,14x,3f14.8)') '   differenced              ', differenced
    write(*,'(a,14x,3es14.2)')'   apart                    ', &
         & abs(f(1:max_order) - differenced)

  end subroutine across_families

end program chained_horizon
