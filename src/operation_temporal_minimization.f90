!=====================================================================!
! A minimizer over an ordered dependency schedule.
!
! The base minimizer states one residual operator and solves one vector.
! A temporal solve is a different member of the same family: it owns a
! directed schedule of local rules, a data branch paired to that schedule,
! and a traversal that releases data once no later local rule depends on
! it. The local rule may itself contain a nonlinear minimizer; this type
! owns the schedule and the inner solver template, not the application
! object that constructs each local rule.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_temporal_minimization

  use util_precision        , only : dp
  use graph_fractal         , only : graph
  use view_directed         , only : directed_graph
  use operation_action      , only : operation
  use operation_minimization, only : minimizer, state
  use operation_driver      , only : driver, pairing
  use operation_residual    , only : residual_operator
  use operation_multigrid   , only : multigrid
  use operation_newton      , only : newton
  use operation_gmres       , only : gmres
  use operation_elimination , only : elimination
  use view_directed_stored  , only : stored_directed_graph
  use field_calculus        , only : FIELD_REAL
  use field_stored          , only : stored_field

  implicit none

  private
  public :: temporal_minimizer

  type, extends(minimizer) :: temporal_minimizer

     class(minimizer), allocatable :: inner

     integer, allocatable :: member_of(:)
     integer, allocatable :: member_order(:)
     logical :: seed_from_previous = .false.
     ! the seed of a member from its predecessor as a linear transfer
     ! of every tuple, one matrix per member: the Taylor shift of the
     ! stored jet over the step between them; a plain copy when absent
     real(dp), allocatable :: seed_transfer(:,:,:)

     type(driver), private :: runner
     logical     , private :: scheduled = .false.

   contains

     procedure :: name  => temporal_minimizer_name
     procedure :: state => temporal_minimizer_state
     procedure :: pair_with
     procedure :: pairing_of
     procedure :: partition
     procedure :: visits
     procedure :: last_dependent_of
     procedure :: released_after
     procedure :: solve

  end type temporal_minimizer

contains

  pure function temporal_minimizer_name(this) result(name)

    class(temporal_minimizer), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'temporal minimizer'

  end function temporal_minimizer_name

  !===================================================================!
  ! STATE. A temporal minimizer may be stated with an ordinary operation,
  ! in which case it delegates to its inner minimizer, or with a driver,
  ! in which case the driver is the schedule this object solves.
  !===================================================================!

  subroutine temporal_minimizer_state(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(temporal_minimizer), intent(inout)        :: this
    class(operation)          , intent(in)           :: action
    class(directed_graph)     , intent(in)           :: context
    type(graph)               , intent(in)           :: unknown_domain
    integer                   , intent(in)           :: num_unknowns
    integer                   , intent(in), optional :: num_components
    class(directed_graph)     , intent(in), optional :: coupling
    type(stored_field)        , intent(in), optional :: stored_inputs(:)

    select type (scheduled_action => action)
    type is (driver)
       if (allocated(this % action)) deallocate(this % action)
       allocate(this % action, source=action)
       if (allocated(this % graph)) deallocate(this % graph)
       allocate(this % graph, source=context)
       if (allocated(this % coupling)) deallocate(this % coupling)
       if (present(coupling)) allocate(this % coupling, source=coupling)
       if (allocated(this % stored)) deallocate(this % stored)
       if (present(stored_inputs)) allocate(this % stored, source=stored_inputs)
       this % unknown_domain  = unknown_domain
       this % residual_domain = unknown_domain
       this % num_unknowns    = num_unknowns
       this % num_residuals   = num_unknowns
       this % num_components  = 1
       if (present(num_components)) this % num_components = max(num_components, 1)
       call this % declare_arguments(0)
       if (allocated(this % affine)) deallocate(this % affine)
       allocate(this % affine(0))
       this % runner = scheduled_action
       this % scheduled = .true.
       this % diagonal_valid = .false.
    class default
       this % scheduled = .false.
       call state(this, action, context, unknown_domain, num_unknowns, &
            & num_components, coupling, stored_inputs)
    end select

  end subroutine temporal_minimizer_state

  subroutine pair_with(this, connection)

    class(temporal_minimizer), intent(inout) :: this
    type(pairing)             , intent(in)    :: connection

    call require_schedule(this)
    call this % runner % pair_with(connection)

  end subroutine pair_with

  function pairing_of(this) result(connection)

    class(temporal_minimizer), intent(in) :: this
    type(pairing) :: connection

    call require_schedule(this)
    connection = this % runner % pairing_of()

  end function pairing_of

  subroutine partition(this, member_of, member_order, seed_from_previous, seed_transfer)

    class(temporal_minimizer), intent(inout)        :: this
    integer                  , intent(in)           :: member_of(:), member_order(:)
    logical                  , intent(in), optional :: seed_from_previous
    real(dp)                 , intent(in), optional :: seed_transfer(:,:,:)

    integer :: last_member

    if (size(member_of) < 1) then
       error stop 'temporal_minimizer: a partition labels at least one unknown'
    end if
    if (size(member_order) < 1) then
       error stop 'temporal_minimizer: a partition states at least one member'
    end if
    if (any(member_of < 1) .or. any(member_order < 1)) then
       error stop 'temporal_minimizer: partition labels are positive'
    end if

    last_member = maxval(member_of)
    if (maxval(member_order) > last_member) then
       error stop 'temporal_minimizer: the member order names the partition'
    end if

    this % member_of    = member_of
    this % member_order = member_order
    this % seed_from_previous = .false.
    if (present(seed_from_previous)) this % seed_from_previous = seed_from_previous
    if (allocated(this % seed_transfer)) deallocate(this % seed_transfer)
    if (present(seed_transfer)) then
       if (size(seed_transfer, 1) /= size(seed_transfer, 2)) then
          error stop 'temporal_minimizer: a seed transfer is square over the tuple'
       end if
       if (size(seed_transfer, 3) /= last_member) then
          error stop 'temporal_minimizer: one seed transfer per member'
       end if
       this % seed_transfer = seed_transfer
    end if

  end subroutine partition

  function visits(this) result(order)

    class(temporal_minimizer), intent(in) :: this
    integer, allocatable :: order(:)

    call require_schedule(this)
    order = this % runner % visits()

  end function visits

  integer function last_dependent_of(this, datum)

    class(temporal_minimizer), intent(in) :: this
    integer                  , intent(in) :: datum

    call require_schedule(this)
    last_dependent_of = this % runner % last_reader_of(datum)

  end function last_dependent_of

  function released_after(this, step) result(vertices)

    class(temporal_minimizer), intent(in) :: this
    integer                  , intent(in) :: step
    integer, allocatable :: vertices(:)

    call require_schedule(this)
    vertices = this % runner % released_after(step)

  end function released_after

  !===================================================================!
  ! SOLVE. Where a schedule is stated, solving is the traversal of that
  ! schedule. Otherwise this object is a shell around its inner minimizer.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(temporal_minimizer), intent(inout) :: this
    real(dp)                  , intent(in)    :: rhs(:)
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(out)   :: achieved

    if (this % scheduled) then
       if (size(rhs) /= 0 .or. size(x) /= 0) then
          error stop 'temporal_minimizer: a scheduled solve has no flat right side'
       end if
       call require_schedule(this)
       call this % runner % evaluate(this % graph)
       achieved = 0.0_dp
       return
    end if

    if (allocated(this % member_of)) then
       call partitioned_solve(this, rhs, x, achieved)
       return
    end if

    if (.not. allocated(this % inner)) then
       error stop 'temporal_minimizer: an inner minimizer is stated'
    end if

    if (allocated(this % coupling)) then
       if (allocated(this % stored)) then
          call this % inner % state(this % action, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components=this % num_components, &
               & coupling=this % coupling, stored_inputs=this % stored)
       else
          call this % inner % state(this % action, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components=this % num_components, &
               & coupling=this % coupling)
       end if
    else
       if (allocated(this % stored)) then
          call this % inner % state(this % action, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components=this % num_components, &
               & stored_inputs=this % stored)
       else
          call this % inner % state(this % action, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components=this % num_components)
       end if
    end if
    call this % inner % solve(rhs, x, achieved)

  end subroutine solve

  subroutine partitioned_solve(this, rhs, x, achieved)

    class(temporal_minimizer), intent(inout) :: this
    real(dp)                  , intent(in)    :: rhs(:)
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(out)   :: achieved

    type(residual_operator) :: sub
    type(stored_directed_graph) :: unknowns
    type(stored_field), allocatable :: inputs(:)
    class(minimizer), allocatable :: local
    real(dp), allocatable :: y(:), zeros(:), piece(:), before(:)
    integer, allocatable :: member(:), previous(:), all_unknowns(:)
    logical, allocatable :: fixed(:)
    integer :: pass, mm, m, count

    if (.not. allocated(this % inner)) then
       error stop 'temporal_minimizer: an inner minimizer is stated'
    end if
    if (this % num_components /= 1) then
       error stop 'temporal_minimizer: a partitioned residual is scalar in each unknown'
    end if
    if (size(x) /= this % num_unknowns .or. size(rhs) /= this % num_unknowns) then
       error stop 'temporal_minimizer: a partitioned solve is sized by its stated unknowns'
    end if
    if (size(this % member_of) /= this % num_unknowns) then
       error stop 'temporal_minimizer: a partition labels every unknown'
    end if
    if (any(abs(rhs) > 0.0_dp)) then
       error stop 'temporal_minimizer: a partitioned residual solve has zero right hand side'
    end if

    select type (residual => this % action)
    class is (residual_operator)
       fixed = residual % fixed_mask()
       x(residual % fixed_unknowns()) = residual % fixed_values()
       allocate(all_unknowns(this % num_unknowns))
       all_unknowns = [(count, count = 1, this % num_unknowns)]

       call this % begin_imbalance()
       call this % evaluate(x, y)
       achieved = this % norm(y - rhs)
       if (this % halted(achieved, 0)) return

       do pass = 1, this % max_iterations
          before = y
          do mm = 1, size(this % member_order)
             m = this % member_order(mm)
             member = pack(all_unknowns, this % member_of == m)
             if (size(member) < 1) cycle
             if (all(fixed(member))) cycle
             if (this % seed_from_previous .and. pass == 1 .and. mm > 1) then
                previous = pack(all_unknowns, this % member_of == this % member_order(mm - 1))
                if (allocated(this % seed_transfer)) then
                   call seed_from_member(x, member, previous, this % seed_transfer(:, :, m))
                else
                   call seed_from_member(x, member, previous)
                end if
             end if
             sub = residual % constrain(member, x)
             if (residual % version() /= 0) then
                call sub % versioned(abs(residual % version()) * size(this % member_order) + m, &
                     & transposed=residual % transpose_version())
             end if
             unknowns = stored_directed_graph(sub % num_unknowns(), tails=[integer ::], heads=[integer ::])
             call local_stored(this, sub, unknowns, inputs)
             call minimizer_over_member(this, member, local)
             if (allocated(inputs)) then
                call local % state(sub, unknowns, unknowns % vertex_set(), sub % num_unknowns(), &
                     & stored_inputs=inputs)
             else
                call local % state(sub, unknowns, unknowns % vertex_set(), sub % num_unknowns())
             end if
             allocate(zeros(sub % num_unknowns()), source=0.0_dp)
             piece = x(member)
             call local % solve(zeros, piece, achieved)
             x(member) = piece
             deallocate(local)
             if (allocated(zeros)) deallocate(zeros)
             if (allocated(inputs)) deallocate(inputs)
          end do
          call this % evaluate(x, y)
          achieved = this % norm(y - rhs)
          if (this % halted(achieved, pass)) exit
          if (all(y == before)) exit
       end do
    class default
       error stop 'temporal_minimizer: a partitioned solve states a residual operator'
    end select

  end subroutine partitioned_solve

  subroutine minimizer_over_member(this, member, local)

    class(temporal_minimizer), intent(in)  :: this
    integer                  , intent(in)  :: member(:)
    class(minimizer), allocatable, intent(out) :: local

    allocate(local, source=this % inner)
    select type (local)
    type is (multigrid)
       if (allocated(local % aggregates)) then
          local % aggregates = compact_labels(local % aggregates(member))
       end if
    type is (newton)
       ! a Newton solve over the member coarsens by the member's
       ! aggregates as well
       if (allocated(local % inner)) then
          select type (inner => local % inner)
          type is (multigrid)
             if (allocated(inner % aggregates)) then
                inner % aggregates = compact_labels(inner % aggregates(member))
             end if
          type is (gmres)
             ! and so does a multigrid preconditioner under its Krylov solve
             if (allocated(inner % preconditioner)) then
                select type (levels => inner % preconditioner)
                type is (multigrid)
                   if (allocated(levels % aggregates)) then
                      levels % aggregates = compact_labels(levels % aggregates(member))
                   end if
                end select
             end if
          type is (elimination)
             ! the rows eliminated over the member are the member's
             if (allocated(inner % eliminated)) then
                inner % eliminated = inner % eliminated(member)
             end if
          end select
       end if
    end select

  end subroutine minimizer_over_member

  function compact_labels(label) result(mapped)

    integer, intent(in) :: label(:)
    integer, allocatable :: mapped(:)

    integer, allocatable :: representative(:)
    integer :: i, at, n

    allocate(mapped(size(label)), representative(size(label)))
    n = 0
    do i = 1, size(label)
       at = 0
       if (n > 0) at = findloc(representative(1:n), label(i), dim=1)
       if (at == 0) then
          n = n + 1
          representative(n) = label(i)
          at = n
       end if
       mapped(i) = at
    end do

  end function compact_labels

  subroutine local_stored(this, residual, unknowns, inputs)

    class(temporal_minimizer), intent(in)  :: this
    type(residual_operator)  , intent(in)  :: residual
    type(stored_directed_graph), intent(in) :: unknowns
    type(stored_field), allocatable, intent(out) :: inputs(:)

    real(dp), allocatable :: values(:), local(:)
    integer :: i

    if (.not. allocated(this % stored)) return

    allocate(inputs(size(this % stored)))
    do i = 1, size(this % stored)
       if (this % stored(i) % value_kind() /= FIELD_REAL) then
          error stop 'temporal_minimizer: residual stored inputs are real fields'
       end if
       if (this % stored(i) % num_components() /= 1) then
          error stop 'temporal_minimizer: residual stored inputs have one component'
       end if
       call this % stored(i) % real_vector(values)
       if (size(values) == residual % num_points()) then
          local = values
       else if (size(values) > 0) then
          if (all(values == values(1))) then
             allocate(local(residual % num_points()), source=values(1))
          else
             error stop 'temporal_minimizer: partitioned residual inputs are constant on a subproblem'
          end if
       else
          error stop 'temporal_minimizer: partitioned residual inputs are constant on a subproblem'
       end if
       inputs(i) = stored_field(this % stored(i) % name(), unknowns % vertex_set(), &
            & residual % num_points())
       call inputs(i) % set_real_vector(local)
       if (allocated(local)) deallocate(local)
    end do

  end subroutine local_stored

  subroutine seed_from_member(x, member, previous, transfer)

    real(dp), intent(inout) :: x(:)
    integer , intent(in)    :: member(:), previous(:)
    real(dp), intent(in), optional :: transfer(:,:)

    integer :: pieces, i, width

    if (size(previous) < 1) return

    ! every tuple of the member is the transfer of the tuple at the
    ! same position of its predecessor
    if (present(transfer)) then
       width = size(transfer, 1)
       if (size(member) /= size(previous) .or. mod(size(member), width) /= 0) then
          error stop 'temporal_minimizer: a seed transfer maps a predecessor of the same tuples'
       end if
       pieces = size(member) / width
       do i = 1, pieces
          x(member((i - 1) * width + 1:i * width)) = &
               & matmul(transfer, x(previous((i - 1) * width + 1:i * width)))
       end do
       return
    end if

    if (size(member) == size(previous)) then
       x(member) = x(previous)
    else if (mod(size(member), size(previous)) == 0) then
       width  = size(previous)
       pieces = size(member) / width
       do i = 1, pieces
          x(member((i - 1) * width + 1:i * width)) = x(previous)
       end do
    else if (mod(size(previous), size(member)) == 0) then
       width = size(member)
       x(member) = x(previous(size(previous) - width + 1:))
    else
       error stop 'temporal_minimizer: a member is seeded from a compatible predecessor'
    end if

  end subroutine seed_from_member

  subroutine require_schedule(this)

    class(temporal_minimizer), intent(in) :: this

    if (.not. this % scheduled) then
       error stop 'temporal_minimizer: a temporal schedule is stated'
    end if

  end subroutine require_schedule

end module operation_temporal_minimization
