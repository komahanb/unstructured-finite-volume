! The shape of a block's jacobian: how far its rows reach either way,
! and how much of the square is empty.
!
! The reach below the diagonal is the history the family reads. The
! reach above it is what a solve cannot march through: a matrix whose
! upper reach is zero is block lower triangular, and needs a forward
! substitution rather than an elimination.
program jacobian_shape

  use util_precision  , only : dp
  use operation_family      , only : family
  use operation_family_bdf  , only : bdf_family
  use operation_family_adams, only : adams_family
  use operation_family_dirk , only : crouzeix_two_stage
  use operation_grid        , only : uniform_grid
  use physics_vanderpol     , only : van_der_pol, van_der_pol_energy
  use gti_expansion         , only : family_holder, expansion
  use gti_driver            , only : cosine
  use gti_march             , only : partition
  use gti_sweeps            , only : jacobian_of
  use view_directed_stored  , only : stored_directed_graph
  use field_stored          , only : stored_field
  use gti_chain             , only : one_functional
  use gti_chain             , only : chain_block, march_chain, &
       & chain_system, chain_systems

  implicit none

  write(*,'(a)') ' '
  write(*,'(a)') '  scheme        unknowns    filled   below   above   per cent full' // &
       & '     largest    on diagonal    largest row'
  call shape_of('bdf 1',   bdf_family(1),        3, 21)
  call shape_of('bdf 2',   bdf_family(2),        3, 21)
  call shape_of('bdf 3',   bdf_family(3),        3, 21)
  call shape_of('adams 2', adams_family(2),      3, 21)
  call shape_of('adams 3', adams_family(3),      3, 21)
  call shape_of('dirk 2',  crouzeix_two_stage(), 3, 21)
  call shape_of('bdf 2',   bdf_family(2),        3, 61)
  call shape_of('bdf 2',   bdf_family(2),        4, 61)

contains

  subroutine shape_of(label, scheme, degrees, instants)

    character(len=*), intent(in) :: label
    class(family)   , intent(in) :: scheme
    integer         , intent(in) :: degrees, instants

    type(family_holder), allocatable :: schemes(:)
    type(chain_block) , allocatable :: chain(:)
    type(expansion), allocatable, target :: tower
    type(chain_system), allocatable :: systems(:)
    integer , allocatable :: added(:)
    real(dp), allocatable :: held(:), dt(:), t(:), a(:,:)
    real(dp) :: achieved, duration, design
    integer :: k, d

    duration = 3.0_dp
    design   = 1.0_dp

    allocate(schemes(1))
    allocate(schemes(1) % scheme, source=scheme)
    added = [instants]

    call partition(duration, instants, dt, t)
    held = [((cosine(d, t(k)), d = 0, degrees - 1), k = 1, &
         &   scheme % history_depth(degrees - 1))]

    call march_chain(schemes, added, van_der_pol(degrees - 1), degrees, &
         & uniform_grid(duration), design, held, chain, tower, dt, t, achieved)
    call chain_systems(chain, tower, [one_functional(van_der_pol_energy(degrees - 1))], degrees, systems)

    call dense_jacobian(chain, degrees, design, a)
    call reported(label, a, degrees, instants)

  end subroutine shape_of

  subroutine reported(label, a, degrees, instants)

    character(len=*), intent(in) :: label
    real(dp)        , intent(in) :: a(:,:)
    integer         , intent(in) :: degrees, instants

    integer  :: n, i, j, filled, below, above
    real(dp) :: biggest, least, on_diagonal, row_most

    n       = size(a, 1)
    biggest = maxval(abs(a))
    least   = 1.0e-12_dp * biggest

    filled = 0
    below  = 0
    above  = 0

    do j = 1, n
       do i = 1, n
          if (abs(a(i, j)) <= least) cycle
          filled = filled + 1
          below  = max(below, i - j)
          above  = max(above, j - i)
       end do
    end do

    on_diagonal = 0.0_dp
    do i = 1, n
       on_diagonal = max(on_diagonal, abs(a(i, i)))
    end do

    row_most = maxval(sum(abs(a), dim=2))

    write(*,'(a,a,i8,i10,i8,i8,f12.2,3es15.4)') '  ', label // repeat(' ', 12 - len(label)), &
         & n, filled, below, above, 100.0_dp * real(filled, dp) / real(n * n, dp), &
         & biggest, on_diagonal, row_most

    associate (u1 => degrees, u2 => instants); end associate

  end subroutine reported

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

end program jacobian_shape
