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

  use iso_fortran_env       , only : int64
  use util_precision        , only : dp
  use graph_fractal         , only : graph
  use view_directed         , only : directed_graph
  use operation_action      , only : operation
  use operation_minimization, only : minimizer, state, restrict, compact_labels
  use operation_minimization, only : solve_result, SOLVE_EVALUATED, SOLVE_STAGNATED, SOLVE_INNER_FAILED
  use operation_driver      , only : driver, pairing, vertex_rule
  use operation_residual    , only : residual_operator
  use field_calculus        , only : FIELD_REAL
  use field_stored          , only : stored_field
  use util_tally            , only : tally

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

     type(driver), private :: schedule
     logical     , private :: scheduled = .false.

   contains

     procedure :: name  => temporal_minimizer_name
     procedure :: state => temporal_minimizer_state
     procedure :: restrict => temporal_minimizer_restrict
     procedure :: bind_account => temporal_minimizer_bind_account
     procedure :: storage_entries => temporal_minimizer_storage_entries
     procedure :: pair_with
     procedure :: pairing_of
     procedure :: partition
     procedure :: visits
     procedure :: next_rule
     procedure :: complete
     procedure :: num_completed
     procedure :: advance
     procedure :: advance_with
     procedure :: set_rule
     procedure :: clear_rule
     procedure :: released_after
     procedure :: released_at
     procedure :: expired_at
     procedure :: live_after
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

    call this % initialize_residual_history()
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
       this % schedule = scheduled_action
       this % scheduled = .true.
       this % diagonal_valid = .false.
    class default
       this % scheduled = .false.
       call state(this, action, context, unknown_domain, num_unknowns, &
            & num_components, coupling, stored_inputs)
    end select

  end subroutine temporal_minimizer_state

  subroutine pair_with(this, connection, position)

    class(temporal_minimizer), intent(inout) :: this
    type(pairing)             , intent(in)    :: connection
    integer                  , intent(in), optional :: position

    call require_schedule(this)
    call this % schedule % pair_with(connection, position)
    call this % initialize_residual_history()

  end subroutine pair_with

  function pairing_of(this) result(connection)

    class(temporal_minimizer), intent(in) :: this
    type(pairing) :: connection

    call require_schedule(this)
    connection = this % schedule % pairing_of()

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

  !===================================================================!
  ! The partition of the selected unknowns: their labels relabelled
  ! compactly in order of first appearance, the member order retaining
  ! the members the selection meets in their stated order, and the
  ! seed transfer of each retained member. The inner template is
  ! restricted by the selection itself. A schedule is over rules, not
  ! unknowns, and a scheduled minimizer refuses restriction.
  !===================================================================!

  subroutine temporal_minimizer_restrict(this, selected)

    class(temporal_minimizer), intent(inout) :: this
    integer                  , intent(in)    :: selected(:)

    integer, allocatable :: labels(:), members(:), order(:)
    integer :: k, at, n

    if (this % scheduled) then
       error stop 'temporal_minimizer: a schedule is over rules, and is not restricted to unknowns'
    end if
    call restrict(this, selected)
    if (allocated(this % member_of)) then
       if (any(selected > size(this % member_of))) then
          error stop 'temporal_minimizer: a restriction selects unknowns of the partition'
       end if
       call compact_labels(this % member_of(selected), labels, members)
       this % member_of = labels
       allocate(order(size(this % member_order)))
       n = 0
       do k = 1, size(this % member_order)
          at = findloc(members, this % member_order(k), dim=1)
          if (at == 0) cycle
          n = n + 1
          order(n) = at
       end do
       if (n < 1) then
          error stop 'temporal_minimizer: a restriction meets an ordered member at least'
       end if
       this % member_order = order(1:n)
       if (allocated(this % seed_transfer)) this % seed_transfer = this % seed_transfer(:, :, members)
    end if
    if (allocated(this % inner)) call this % inner % restrict(selected)

  end subroutine temporal_minimizer_restrict

  ! the inner template's entries over the same unknowns
  pure integer(int64) function temporal_minimizer_storage_entries(this, num_unknowns) result(entries)

    class(temporal_minimizer), intent(in) :: this
    integer                  , intent(in) :: num_unknowns

    entries = 0_int64
    if (allocated(this % inner)) entries = this % inner % storage_entries(num_unknowns)

  end function temporal_minimizer_storage_entries

  function visits(this) result(order)

    class(temporal_minimizer), intent(in) :: this
    integer, allocatable :: order(:)

    call require_schedule(this)
    order = this % schedule % visits()

  end function visits

  integer function next_rule(this)
    class(temporal_minimizer), intent(in) :: this
    call require_schedule(this)
    next_rule = this % schedule % next_rule()
  end function next_rule

  logical function complete(this)
    class(temporal_minimizer), intent(in) :: this
    call require_schedule(this)
    complete = this % schedule % complete()
  end function complete

  integer function num_completed(this)
    class(temporal_minimizer), intent(in) :: this
    call require_schedule(this)
    num_completed = this % schedule % num_completed()
  end function num_completed

  subroutine advance(this, executed)
    class(temporal_minimizer), intent(inout) :: this
    integer, intent(out), optional :: executed
    call require_schedule(this)
    call this % schedule % advance(this % graph, executed)
    if (this % schedule % complete()) then
       call this % initialize_residual_history()
       call this % record_result(0.0_dp, this % schedule % num_completed(), SOLVE_EVALUATED)
    end if
  end subroutine advance

  subroutine advance_with(this, rule, executed)
    class(temporal_minimizer), intent(inout) :: this
    class(vertex_rule), intent(inout) :: rule
    integer, intent(out), optional :: executed
    call require_schedule(this)
    call this % schedule % advance_with(this % graph, rule, executed)
    if (this % schedule % complete()) then
       call this % initialize_residual_history()
       call this % record_result(0.0_dp, this % schedule % num_completed(), SOLVE_EVALUATED)
    end if
  end subroutine advance_with

  subroutine set_rule(this, vertex, rule)
    class(temporal_minimizer), intent(inout) :: this
    integer, intent(in) :: vertex
    class(operation), intent(in) :: rule
    call require_schedule(this)
    call this % schedule % set_rule(vertex, rule)
  end subroutine set_rule

  subroutine clear_rule(this, vertex)
    class(temporal_minimizer), intent(inout) :: this
    integer, intent(in) :: vertex
    call require_schedule(this)
    call this % schedule % clear_rule(vertex)
  end subroutine clear_rule

  function released_after(this, step) result(vertices)

    class(temporal_minimizer), intent(in) :: this
    integer                  , intent(in) :: step
    integer, allocatable :: vertices(:)

    call require_schedule(this)
    vertices = this % schedule % released_after(step)

  end function released_after

  function released_at(this, step) result(vertices)
    class(temporal_minimizer), intent(in) :: this
    integer, intent(in) :: step
    integer, allocatable :: vertices(:)
    call require_schedule(this)
    vertices = this % schedule % released_at(step)
  end function released_at

  function expired_at(this, step) result(vertices)
    class(temporal_minimizer), intent(in) :: this
    integer, intent(in) :: step
    integer, allocatable :: vertices(:)
    call require_schedule(this)
    vertices = this % schedule % expired_at(step)
  end function expired_at

  function live_after(this, step) result(vertices)
    class(temporal_minimizer), intent(in) :: this
    integer, intent(in) :: step
    integer, allocatable :: vertices(:)
    call require_schedule(this)
    vertices = this % schedule % live_after(step)
  end function live_after

  !===================================================================!
  ! SOLVE. Where a schedule is stated, solving is the traversal of that
  ! schedule. Otherwise this object is a shell around its inner minimizer.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(temporal_minimizer), intent(inout) :: this
    real(dp)                  , intent(in)    :: rhs(:)
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(out)   :: achieved
    type(solve_result) :: outcome

    if (this % scheduled) then
       if (size(rhs) /= 0 .or. size(x) /= 0) then
          error stop 'temporal_minimizer: a scheduled solve has no flat right side'
       end if
       call require_schedule(this)
       call this % schedule % evaluate(this % graph)
       achieved = 0.0_dp
       call this % initialize_residual_history()
       call this % record_result(achieved, size(this % schedule % visits()), SOLVE_EVALUATED)
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
    call this % initialize_residual_history()
    call this % inner % solve(rhs, x, achieved)
    outcome = this % inner % result()
    call this % record_residual_norm(this % inner % initial_residual_norm())
    call this % record_result(outcome % residual, outcome % iterations, outcome % reason)

  end subroutine solve

  subroutine partitioned_solve(this, rhs, x, achieved)

    class(temporal_minimizer), intent(inout) :: this
    real(dp)                  , intent(in)    :: rhs(:)
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(out)   :: achieved

    type(residual_operator) :: sub
    type(stored_field), allocatable :: inputs(:)
    class(minimizer), allocatable :: local
    real(dp), allocatable :: y(:), zeros(:), member_solution(:), previous_residual(:)
    integer, allocatable :: member(:), previous(:), all_unknowns(:)
    logical, allocatable :: fixed(:)
    integer :: pass, mm, m, count
    type(solve_result) :: outcome

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
       fixed = residual % fixed_indicator()
       x(residual % fixed_unknowns()) = residual % fixed_values()
       allocate(all_unknowns(this % num_unknowns))
       all_unknowns = [(count, count = 1, this % num_unknowns)]

       call this % initialize_residual_history()
       call this % evaluate(x, y)
       achieved = this % norm(y - rhs)
       if (this % terminated(achieved, 0)) return

       do pass = 1, this % max_iterations
          previous_residual = y
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
             call local_stored(this, residual, member, sub, inputs)
             ! the inner minimizer restricted to the member: each
             ! solver maps its own metadata and its children's
             allocate(local, source=this % inner)
             call local % restrict(member)
             if (allocated(inputs)) then
                call local % state(sub, sub % unknown_graph(), sub % unknown_domain(), sub % num_unknowns(), &
                     & stored_inputs=inputs)
             else
                call local % state(sub, sub % unknown_graph(), sub % unknown_domain(), sub % num_unknowns())
             end if
             allocate(zeros(sub % num_unknowns()), source=0.0_dp)
             member_solution = x(member)
             call local % solve(zeros, member_solution, achieved)
             outcome = local % result()
             if (outcome % failed()) then
                call this % evaluate(x, y)
                achieved = this % norm(y - rhs)
                call this % record_result(achieved, pass - 1, SOLVE_INNER_FAILED)
                return
             end if
             x(member) = member_solution
             deallocate(local)
             if (allocated(zeros)) deallocate(zeros)
             if (allocated(inputs)) deallocate(inputs)
          end do
          call this % evaluate(x, y)
          achieved = this % norm(y - rhs)
          if (this % terminated(achieved, pass)) exit
          if (all(y == previous_residual)) then
             call this % record_result(achieved, pass, SOLVE_STAGNATED)
             exit
          end if
       end do
    class default
       error stop 'temporal_minimizer: a partitioned solve states a residual operator'
    end select

  end subroutine partitioned_solve

  !===================================================================!
  ! THE STORED INPUTS RESTRICTED TO A MEMBER. A stored input of a
  ! residual is a field on the residual's point domain P with one
  ! value per point; the constrained residual's points are the
  ! selected points, so the input restricted to the member is the
  ! field on the constrained residual's P' of the values at those
  ! points. Invalid input: a stored input that is not a real field on
  ! P with one value per point.
  !===================================================================!

  subroutine local_stored(this, residual, member, sub, inputs)

    class(temporal_minimizer), intent(in)  :: this
    class(residual_operator) , intent(in)  :: residual
    integer                  , intent(in)  :: member(:)
    class(residual_operator) , intent(in)  :: sub
    type(stored_field), allocatable, intent(out) :: inputs(:)

    real(dp), allocatable :: values(:)
    integer , allocatable :: points(:)
    type(graph) :: whole, selected
    integer :: i

    if (.not. allocated(this % stored)) return

    whole    = residual % design_domain()
    selected = sub % design_domain()
    points   = residual % selected_points(member)
    allocate(inputs(size(this % stored)))
    do i = 1, size(this % stored)
       if (this % stored(i) % value_kind() /= FIELD_REAL) then
          error stop 'temporal_minimizer: residual stored inputs are real fields'
       end if
       if (.not. this % stored(i) % defined_on(whole)) then
          error stop 'temporal_minimizer: residual stored inputs are defined on the residual''s point domain'
       end if
       if (this % stored(i) % num_components() /= 1 .or. &
            & this % stored(i) % num_entries() /= residual % num_points()) then
          error stop 'temporal_minimizer: residual stored inputs have one value per point'
       end if
       call this % stored(i) % real_vector(values)
       inputs(i) = stored_field(this % stored(i) % name(), selected, size(points))
       call inputs(i) % set_real_vector(values(points))
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

  !===================================================================!
  ! Bind the account of this minimizer and of its children.
  !===================================================================!

  subroutine temporal_minimizer_bind_account(this, account)

    class(temporal_minimizer), intent(inout) :: this
    type(tally), pointer, intent(in) :: account

    this % account => account
    if (allocated(this % inner)) call this % inner % bind_account(account)

  end subroutine temporal_minimizer_bind_account

end module operation_temporal_minimization
