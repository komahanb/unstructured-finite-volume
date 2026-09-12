!=====================================================================!
! The elimination of unit-diagonal rows before a linear solve: a
! minimizer like every other, stated on the whole system and solving
! its Schur complement on the unknowns it retains.
!
! The unknowns split into the retained set K and the eliminated set E,
! one flag per unknown. The rows eliminated are the tying rows of a
! derived component - a time derivative read from its family, a
! spatial derivative read from its fit - each linear with a unit
! diagonal once divided by it, reading retained unknowns and
! eliminated unknowns of rows before it: a derivative of the second
! degree reads the first's, so the eliminated block is unit triangular
! in some order,
!
!      [ J_KK    J_KE  ] [ x_K ]   [ b_K ]
!      [ J_EK   I + N  ] [ x_E ] = [ b_E ] ,
!
! N strictly triangular in that order. Then x_E is read by one
! substitution, x_E = (I + N)^-1 (b_E - J_EK x_K), and the retained
! unknowns solve the complement
!
!      (J_KK - J_KE M) x_K = b_K - J_KE (I + N)^-1 b_E ,   M = (I + N)^-1 J_EK ,
!
! M assembled row by row in the order, each row the row of J_EK less
! the rows it reads. The complement is a stencil from the tangent's
! own triples and is passed to the inner minimizer, which never sees
! the eliminated rows: its vectors, its blocks and its Krylov basis
! are over K alone. The solution is the solution of the whole system,
! so a Newton iteration over it takes the same steps whether or not
! the rows are eliminated; only the linear solve changes shape.
!
! The tangent must be explicit: the complement is read from triples,
! and a matrix-free linearisation has none. An eliminated row without
! a diagonal, or whose reads of eliminated unknowns close a cycle, is
! not of the form above, and the solve stops there naming it.
!
! Every coefficient of M and of the complement is summed over its
! dependency paths before a dependant reads it or the inner minimizer
! multiplies it, so the storage counts coefficients, never paths or
! uncombined products. The entries stored are bounded by max_entries,
! the sum of five accounts declared in elimination_storage; a
! construction beyond the limit, or beyond the largest count an index
! array addresses, is refused before its allocation and reported by
! the solve result SOLVE_STORAGE_EXCEEDED.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_elimination

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use iso_fortran_env       , only : int64
  use util_precision        , only : dp
  use graph_fractal         , only : graph
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use operation_action      , only : operation
  use operation_minimization, only : minimizer, state, restrict, saturated_sum
  use operation_minimization, only : solve_result, SOLVE_INNER_FAILED, SOLVE_EXHAUSTED, SOLVE_CONTINUE
  use operation_minimization, only : SOLVE_STORAGE_EXCEEDED
  use operation_stencil     , only : stencil
  use field_stored          , only : stored_field

  implicit none

  private
  public :: elimination, elimination_storage

  ! one entry: a coefficient with its column index
  integer, parameter, public :: entry_bytes = (storage_size(1.0_dp) + storage_size(1)) / 8

  !===================================================================!
  ! The entries an elimination stores, by account, and the limit their
  ! sum is admitted under. Every count is a number of entries; bytes
  ! multiply by entry_bytes. The pattern of the complement stencil and
  ! the minimizer objects themselves are outside these counts.
  !===================================================================!

  type :: elimination_storage

     ! nnz(J) triples read from the stencil and the retained block
     ! J_KK split out, released once the complement is formed
     integer(int64) :: input = 0_int64

     ! J_KE, J_EK, N, the diagonals and M = (I + N)^-1 J_EK, retained
     ! through the solve
     integer(int64) :: substitution = 0_int64

     ! the position of every unknown, one row accumulator over the
     ! retained columns and two permutations of the complement's
     ! entries ordering each row's columns, all released after the
     ! statement; the index permutations grouping the input triples by
     ! row are proportional to the input account and are not counted
     ! again
     integer(int64) :: temporary = 0_int64

     ! nnz(S), the combined complement the inner minimizer is stated on
     integer(int64) :: schur = 0_int64

     ! the entries the inner minimizer declares for the retained unknowns
     integer(int64) :: factorisation = 0_int64

     ! the limit on the sum, in entries
     integer(int64) :: limit = int(huge(1), int64)

   contains

     procedure :: total => storage_total
     procedure :: bytes => storage_bytes
     procedure :: admitted => storage_admitted

  end type elimination_storage

  type, extends(minimizer) :: elimination

     ! one flag per unknown of the stated system: true for a row
     ! eliminated before the inner solve
     logical, allocatable :: eliminated(:)

     class(minimizer), allocatable :: inner

     ! the limit on the entries stored, the sum of the five accounts
     integer :: max_entries = huge(1)

     ! the accounts of the last statement, or at its refusal the
     ! accounts including the increment refused
     type(elimination_storage) :: storage

     ! whether the last statement was refused within its limit
     logical, private :: storage_exceeded = .false.

     ! the unknowns of each set, in the order the inner reads them
     integer, allocatable, private :: retained_at(:), eliminated_at(:)

     ! the coupling blocks as triples in the numbering of each set:
     ! J_KE by row in K and column in E, J_EK by row in E and column in K
     integer , allocatable, private :: ke_row(:), ke_column(:), ek_row(:), ek_column(:)
     real(dp), allocatable, private :: ke_weight(:), ek_weight(:)

     ! N by row: the eliminated unknowns each eliminated row reads,
     ! first_n(e) to first_n(e + 1) - 1, and the substitution order;
     ! every eliminated row is read divided by its diagonal
     integer , allocatable, private :: first_n(:), n_column(:), order(:)
     real(dp), allocatable, private :: n_weight(:), diagonal(:)

   contains

     procedure :: name  => elimination_name
     procedure :: state => elimination_state
     procedure :: restrict => elimination_restrict
     procedure :: storage_entries => elimination_storage_entries
     procedure :: solve => elimination_solve

  end type elimination

contains

  pure function elimination_name(this) result(name)

    class(elimination), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % inner)) then
       name = 'elimination before ' // this % inner % name()
    else
       name = 'elimination'
    end if

  end function elimination_name

  pure integer(int64) function storage_total(this) result(total)

    class(elimination_storage), intent(in) :: this

    total = saturated_sum(this % input, this % substitution)
    total = saturated_sum(total, this % temporary)
    total = saturated_sum(total, this % schur)
    total = saturated_sum(total, this % factorisation)

  end function storage_total

  pure integer(int64) function storage_bytes(this) result(bytes)

    class(elimination_storage), intent(in) :: this

    bytes = this % total()
    if (real(bytes, dp) * real(entry_bytes, dp) > real(huge(bytes), dp)) then
       bytes = huge(bytes)
    else
       bytes = bytes * int(entry_bytes, int64)
    end if

  end function storage_bytes

  !===================================================================!
  ! Whether the accounts are within the limit and within the largest
  ! count a default-integer index array addresses.
  !===================================================================!

  pure logical function storage_admitted(this) result(admitted)

    class(elimination_storage), intent(in) :: this

    admitted = this % total() <= min(this % limit, int(huge(1), int64))

  end function storage_admitted

  !===================================================================!
  ! State the whole system, split its triples by set, divide the
  ! eliminated rows by their diagonals, assemble the complement and
  ! state the inner minimizer on it. Every allocation proportional to
  ! a count of coefficients follows a check of the accounts against
  ! the limit; a refused statement stores no partition.
  !===================================================================!

  subroutine elimination_state(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(elimination)   , intent(inout) :: this
    class(operation)     , intent(in)    :: action
    class(directed_graph), intent(in)    :: context
    type(graph)          , intent(in)    :: unknown_domain
    integer              , intent(in)    :: num_unknowns
    integer              , intent(in), optional :: num_components
    class(directed_graph), intent(in), optional :: coupling
    type(stored_field)   , intent(in), optional :: stored_inputs(:)

    integer , allocatable :: rows(:), columns(:), position(:), kk_row(:), kk_column(:)
    integer , allocatable :: first_ek(:), by_row(:), n_row(:), first_m(:), count_m(:), m_column(:)
    integer , allocatable :: s_row(:), s_column(:), first_s(:)
    real(dp), allocatable :: weights(:), kk_weight(:), m_weight(:), s_weight(:)
    type(stencil) :: complement
    type(stored_directed_graph) :: retained
    integer(int64) :: stored_substitution
    integer :: n, nk, ne, e, r, c, i, nkk, nke, nek, nn, m, total, capacity
    logical :: within

    call state(this, action, context, unknown_domain, num_unknowns, &
         & num_components, coupling, stored_inputs)
    call clear_partition(this)
    this % storage_exceeded = .false.
    if (this % max_entries < 1) then
       error stop 'elimination: the storage limit is one entry at least'
    end if
    this % storage = elimination_storage(limit = int(this % max_entries, int64))

    if (present(num_components)) then
       if (num_components > 1) then
          error stop 'elimination: one value per unknown; a wider right-hand side is not eliminated'
       end if
    end if
    if (.not. allocated(this % inner)) then
       error stop 'elimination: an inner minimizer solves the complement, and none is stored'
    end if
    n = num_unknowns
    if (.not. allocated(this % eliminated)) then
       error stop 'elimination: one flag per unknown states which rows are eliminated'
    end if
    if (size(this % eliminated) /= n) then
       error stop 'elimination: one flag per unknown'
    end if

    select type (action)
    type is (stencil)
       m = action % pattern % num_edges()
    class default
       error stop 'elimination: the complement is read from the explicit tangent, and this &
            &statement is a matrix-vector product without one'
    end select

    nk = count(.not. this % eliminated)
    ne = n - nk
    if (nk < 1) then
       error stop 'elimination: an unknown is retained at least'
    end if

    ! the input triples, the work-space and the inner's requirement
    ! are known before any of them is allocated
    this % storage % input         = int(m, int64)
    this % storage % temporary     = int(n, int64) + int(nk, int64)
    this % storage % factorisation = this % inner % storage_entries(nk)
    if (.not. admitted(this)) return

    select type (action)
    type is (stencil)
       call action % entries(rows, columns, weights)
    end select

    ! the position of every unknown within its set
    allocate(position(n))
    nk = 0
    ne = 0
    do i = 1, n
       if (this % eliminated(i)) then
          ne = ne + 1
          position(i) = ne
       else
          nk = nk + 1
          position(i) = nk
       end if
    end do
    this % retained_at   = pack([(i, i = 1, n)], .not. this % eliminated)
    this % eliminated_at = pack([(i, i = 1, n)],       this % eliminated)

    ! the diagonal of every eliminated row, its own scale
    allocate(this % diagonal(ne), source=0.0_dp)
    do e = 1, m
       if (rows(e) == columns(e) .and. this % eliminated(rows(e))) then
          this % diagonal(position(rows(e))) = this % diagonal(position(rows(e))) + weights(e)
       end if
    end do
    if (any(this % diagonal == 0.0_dp)) then
       error stop 'elimination: an eliminated row has no diagonal: it is not the tying row of &
            &its unknown; under a staged family the states at the stages are the tied &
            &components, so its time derivatives are stated as rows'
    end if

    ! the triples of every block counted before the blocks are stored
    nkk = 0
    nke = 0
    nek = 0
    nn  = 0
    do e = 1, m
       r = rows(e)
       c = columns(e)
       if (this % eliminated(r)) then
          if (this % eliminated(c)) then
             if (r /= c .and. weights(e) /= 0.0_dp) nn = nn + 1
          else
             nek = nek + 1
          end if
       else
          if (this % eliminated(c)) then
             nke = nke + 1
          else
             nkk = nkk + 1
          end if
       end if
    end do
    this % storage % input        = int(m, int64) + int(nkk, int64)
    this % storage % substitution = int(nke, int64) + int(nek, int64) + int(nn, int64) + int(ne, int64)
    if (.not. admitted(this)) return

    allocate(kk_row(nkk), kk_column(nkk), kk_weight(nkk))
    allocate(this % ke_row(nke), this % ke_column(nke), this % ke_weight(nke))
    allocate(this % ek_row(nek), this % ek_column(nek), this % ek_weight(nek))
    allocate(n_row(nn), this % n_column(nn), this % n_weight(nn))
    nkk = 0
    nke = 0
    nek = 0
    nn  = 0
    do e = 1, m
       r = rows(e)
       c = columns(e)
       if (this % eliminated(r)) then
          if (this % eliminated(c)) then
             if (r /= c .and. weights(e) /= 0.0_dp) then
                nn = nn + 1
                n_row(nn)           = position(r)
                this % n_column(nn) = position(c)
                this % n_weight(nn) = weights(e) / this % diagonal(position(r))
             end if
          else
             nek = nek + 1
             this % ek_row(nek)    = position(r)
             this % ek_column(nek) = position(c)
             this % ek_weight(nek) = weights(e) / this % diagonal(position(r))
          end if
       else
          if (this % eliminated(c)) then
             nke = nke + 1
             this % ke_row(nke)    = position(r)
             this % ke_column(nke) = position(c)
             this % ke_weight(nke) = weights(e)
          else
             nkk = nkk + 1
             kk_row(nkk)    = position(r)
             kk_column(nkk) = position(c)
             kk_weight(nkk) = weights(e)
          end if
       end if
    end do
    deallocate(rows, columns, weights)

    ! N and J_EK by row
    call by_rows(ne, n_row, by_row, this % first_n)
    this % n_column = this % n_column(by_row)
    this % n_weight = this % n_weight(by_row)
    call by_rows(ne, this % ek_row, by_row, first_ek)

    ! the substitution order: a row after every eliminated row it reads
    call ordered(ne, this % first_n, this % n_column, this % order)

    ! M = (I + N)^-1 J_EK row by row in the substitution order, then
    ! the complement J_KK - J_KE M row by row over K, each within the
    ! capacity its account admits; a row beyond it is refused with the
    ! account naming the requirement
    stored_substitution = this % storage % substitution
    capacity = admissible(this, stored_substitution) - int(stored_substitution)
    call substitution_rows(this, nk, first_ek, by_row, capacity, first_m, count_m, m_column, m_weight, total, within)
    this % storage % substitution = stored_substitution + int(total, int64)
    if (.not. within) then
       call refused(this)
       return
    end if
    capacity = admissible(this, 0_int64)
    call complement_rows(this, nk, kk_row, kk_column, kk_weight, first_m, count_m, m_column, m_weight, capacity, &
         & s_row, s_column, s_weight, first_s, total, within)
    this % storage % schur = int(total, int64)
    if (.not. within) then
       call refused(this)
       return
    end if

    ! each row's columns ascending, by two permutations of the entries
    this % storage % temporary = this % storage % temporary + 2_int64 * int(total, int64)
    if (.not. admitted(this)) return
    call ordered_columns(nk, first_s, s_row, s_column, s_weight)

    complement = stencil(s_row, s_column, s_weight, spread(0.0_dp, 1, nk), 'eliminated tangent')
    call complement % versioned(action % version(), action % transpose_version())
    retained = stored_directed_graph(nk, tails=[integer ::], heads=[integer ::])
    call this % inner % state(complement, complement % pattern, retained % vertex_set(), nk, &
         & num_components = 1, coupling = complement % pattern)

  end subroutine elimination_state

  !===================================================================!
  ! M = (I + N)^-1 J_EK row by row in the substitution order: the row
  ! of J_EK less the rows of M it reads, equal retained columns summed
  ! before another row reads this one, so storage counts the nonzero
  ! coefficients of M and not dependency paths. total is the entries
  ! stored, or where a row exceeds the capacity the requirement
  ! refused, within then false.
  !===================================================================!

  subroutine substitution_rows(this, nk, first_ek, ek_by_row, capacity, first_m, count_m, m_column, m_weight, &
       & total, within)

    class(elimination), intent(in) :: this
    integer, intent(in) :: nk, first_ek(:), ek_by_row(:), capacity
    integer , allocatable, intent(out) :: first_m(:), count_m(:), m_column(:)
    real(dp), allocatable, intent(out) :: m_weight(:)
    integer, intent(out) :: total
    logical, intent(out) :: within

    integer , allocatable :: column_row(:), active_columns(:)
    real(dp), allocatable :: row_weight(:)
    integer :: ne, i, j, k, e, c, num_active_columns, nonzero

    ne = size(this % order)
    allocate(first_m(ne), count_m(ne))
    allocate(column_row(nk), source=0)
    allocate(active_columns(nk), row_weight(nk), m_column(0), m_weight(0))
    within = .true.
    total = 0
    do i = 1, ne
       e = this % order(i)
       num_active_columns = 0
       do j = first_ek(e), first_ek(e + 1) - 1
          call accumulated(this % ek_column(ek_by_row(j)), this % ek_weight(ek_by_row(j)))
       end do
       do j = this % first_n(e), this % first_n(e + 1) - 1
          c = this % n_column(j)
          do k = first_m(c), first_m(c) + count_m(c) - 1
             call accumulated(m_column(k), -this % n_weight(j) * m_weight(k))
          end do
       end do
       ! the coefficients stored are the nonzero ones: counted exactly
       ! where the active columns would exceed the capacity
       nonzero = num_active_columns
       if (total + nonzero > capacity) then
          nonzero = 0
          do j = 1, num_active_columns
             if (row_weight(active_columns(j)) /= 0.0_dp) nonzero = nonzero + 1
          end do
          if (total + nonzero > capacity) then
             total = total + nonzero
             within = .false.
             return
          end if
       end if
       first_m(e) = total + 1
       call reserve_coefficients(total + nonzero, capacity, total, m_column, m_weight)
       do j = 1, num_active_columns
          c = active_columns(j)
          if (row_weight(c) == 0.0_dp) cycle
          total = total + 1
          m_column(total) = c
          m_weight(total) = row_weight(c)
       end do
       count_m(e) = total - first_m(e) + 1
    end do

  contains

    ! a coefficient of the row being formed: the first of its column
    ! begins the sum at zero
    subroutine accumulated(column, weight)
      integer , intent(in) :: column
      real(dp), intent(in) :: weight
      if (column_row(column) /= e) then
         column_row(column) = e
         num_active_columns = num_active_columns + 1
         active_columns(num_active_columns) = column
         row_weight(column) = 0.0_dp
      end if
      row_weight(column) = row_weight(column) + weight
    end subroutine accumulated

  end subroutine substitution_rows

  !===================================================================!
  ! The complement row by row over K: J_KK less J_KE M, every retained
  ! column summed once per row in the order the triples of that row
  ! arrive - J_KK's, then each J_KE entry's product with its row of M
  ! - so no uncombined product is formed. first_s(r) is the first
  ! entry of row r. total is the entries stored, or where a row
  ! exceeds the capacity the requirement refused, within then false.
  !===================================================================!

  subroutine complement_rows(this, nk, kk_row, kk_column, kk_weight, first_m, count_m, m_column, m_weight, &
       & capacity, s_row, s_column, s_weight, first_s, total, within)

    class(elimination), intent(in) :: this
    integer , intent(in) :: nk, kk_row(:), kk_column(:), first_m(:), count_m(:), m_column(:), capacity
    real(dp), intent(in) :: kk_weight(:), m_weight(:)
    integer , allocatable, intent(out) :: s_row(:), s_column(:), first_s(:)
    real(dp), allocatable, intent(out) :: s_weight(:)
    integer, intent(out) :: total
    logical, intent(out) :: within

    integer , allocatable :: kk_by_row(:), first_kk(:), ke_by_row(:), first_ke(:), column_row(:), active_columns(:)
    real(dp), allocatable :: row_weight(:)
    integer :: r, j, k, e, c, num_active_columns

    call by_rows(nk, kk_row, kk_by_row, first_kk)
    call by_rows(nk, this % ke_row, ke_by_row, first_ke)
    allocate(column_row(nk), source=0)
    allocate(active_columns(nk), row_weight(nk), s_row(0), s_column(0), s_weight(0), first_s(nk + 1))
    within = .true.
    total = 0
    first_s(1) = 1
    do r = 1, nk
       num_active_columns = 0
       do j = first_kk(r), first_kk(r + 1) - 1
          call combined(kk_column(kk_by_row(j)), kk_weight(kk_by_row(j)))
       end do
       do j = first_ke(r), first_ke(r + 1) - 1
          e = ke_by_row(j)
          c = this % ke_column(e)
          do k = first_m(c), first_m(c) + count_m(c) - 1
             call combined(m_column(k), -this % ke_weight(e) * m_weight(k))
          end do
       end do
       if (total + num_active_columns > capacity) then
          total = total + num_active_columns
          within = .false.
          return
       end if
       call reserve_coefficients(total + num_active_columns, capacity, total, s_column, s_weight)
       call reserve_rows(size(s_column), total, s_row)
       do j = 1, num_active_columns
          c = active_columns(j)
          total = total + 1
          s_row(total)    = r
          s_column(total) = c
          s_weight(total) = row_weight(c)
       end do
       first_s(r + 1) = total + 1
    end do

  contains

    ! a coefficient of the row being formed: the first of its column
    ! is the sum's first term
    subroutine combined(column, weight)
      integer , intent(in) :: column
      real(dp), intent(in) :: weight
      if (column_row(column) /= r) then
         column_row(column) = r
         num_active_columns = num_active_columns + 1
         active_columns(num_active_columns) = column
         row_weight(column) = weight
      else
         row_weight(column) = row_weight(column) + weight
      end if
    end subroutine combined

  end subroutine complement_rows

  !===================================================================!
  ! The complement's entries in ascending column order within each
  ! row: grouped by column, then that sequence by row, both groupings
  ! stable, so the triples are in the order of a combined triple list.
  !===================================================================!

  subroutine ordered_columns(nk, first_s, s_row, s_column, s_weight)

    integer, intent(in) :: nk, first_s(:)
    integer , allocatable, intent(inout) :: s_row(:), s_column(:)
    real(dp), allocatable, intent(inout) :: s_weight(:)

    integer, allocatable :: by_column(:), first_column(:), by_row_column(:), next(:)
    integer :: total, j, e, r

    total = first_s(nk + 1) - 1
    call by_rows(nk, s_column(1:total), by_column, first_column)
    allocate(by_row_column(total))
    next = first_s(1:nk)
    do j = 1, total
       e = by_column(j)
       r = s_row(e)
       by_row_column(next(r)) = e
       next(r) = next(r) + 1
    end do
    deallocate(by_column)
    s_row    = s_row(by_row_column)
    s_column = s_column(by_row_column)
    s_weight = s_weight(by_row_column)

  end subroutine ordered_columns

  !===================================================================!
  ! Whether the accounts are admitted; a refusal releases the partial
  ! partition and records the outcome the solve will report.
  !===================================================================!

  logical function admitted(this)

    class(elimination), intent(inout) :: this

    admitted = this % storage % admitted()
    if (.not. admitted) call refused(this)

  end function admitted

  subroutine refused(this)

    class(elimination), intent(inout) :: this

    call clear_partition(this)
    this % storage_exceeded = .true.
    call this % record_result(0.0_dp, 0, SOLVE_STORAGE_EXCEEDED)

  end subroutine refused

  !===================================================================!
  ! The largest capacity an account may grow to within the limit: the
  ! limit less every other account, where the account reads entries.
  !===================================================================!

  pure integer function admissible(this, account) result(capacity)

    class(elimination), intent(in) :: this
    integer(int64)    , intent(in) :: account

    integer(int64) :: others

    others = this % storage % total() - account
    capacity = int(min(this % storage % limit - others, int(huge(1), int64)))

  end function admissible

  !===================================================================!
  ! Capacity for required coefficients, doubling within the admissible
  ! capacity, never below what is required.
  !===================================================================!

  subroutine reserve_coefficients(required, admissible_capacity, used, columns, weights)

    integer, intent(in) :: required, admissible_capacity, used
    integer, allocatable, intent(inout) :: columns(:)
    real(dp), allocatable, intent(inout) :: weights(:)

    integer, allocatable :: larger_columns(:)
    real(dp), allocatable :: larger_weights(:)
    integer :: capacity

    capacity = size(columns)
    if (required <= capacity) return
    if (capacity <= huge(capacity) / 2) then
       capacity = 2 * capacity
    else
       capacity = huge(capacity)
    end if
    capacity = max(required, min(capacity, admissible_capacity))
    allocate(larger_columns(capacity), larger_weights(capacity))
    larger_columns(1:used) = columns(1:used)
    larger_weights(1:used) = weights(1:used)
    call move_alloc(larger_columns, columns)
    call move_alloc(larger_weights, weights)

  end subroutine reserve_coefficients

  ! the row indices at the capacity of the columns beside them
  subroutine reserve_rows(capacity, used, rows)

    integer, intent(in) :: capacity, used
    integer, allocatable, intent(inout) :: rows(:)

    integer, allocatable :: larger(:)

    if (capacity <= size(rows)) return
    allocate(larger(capacity))
    larger(1:used) = rows(1:used)
    call move_alloc(larger, rows)

  end subroutine reserve_rows

  !===================================================================!
  ! The entries the last statement stored, its complement's included:
  ! an enclosing minimizer reads the accounts of the statement made.
  !===================================================================!

  pure integer(int64) function elimination_storage_entries(this, num_unknowns) result(entries)

    class(elimination), intent(in) :: this
    integer           , intent(in) :: num_unknowns

    associate (u1 => num_unknowns); end associate
    entries = this % storage % total()

  end function elimination_storage_entries

  !===================================================================!
  ! The flags of the selected unknowns, in the selection's order. The
  ! inner minimizer is over the retained set K of the whole; the
  ! selection induced on it is the position in K of every retained
  ! selected unknown, in the selection's order. The partition of the
  ! whole is invalid on the selection and is released. Invalid input:
  ! a selection retaining no unknown, or one beyond the stated flags.
  !===================================================================!

  subroutine elimination_restrict(this, selected)

    class(elimination), intent(inout) :: this
    integer           , intent(in)    :: selected(:)

    integer, allocatable :: retained_position(:), retained_of_selection(:)
    integer :: i, n

    call restrict(this, selected)
    call clear_partition(this)
    this % storage_exceeded = .false.
    this % storage = elimination_storage(limit = int(this % max_entries, int64))
    if (.not. allocated(this % eliminated)) then
       error stop 'elimination: one flag per unknown states which rows are eliminated before &
            &the unknowns retained by a restriction are known'
    end if
    if (any(selected > size(this % eliminated))) then
       error stop 'elimination: a restriction selects unknowns of the stated flags'
    end if

    n = size(this % eliminated)
    allocate(retained_position(n))
    retained_position = 0
    do i = 1, n
       if (i > 1) retained_position(i) = retained_position(i - 1)
       if (.not. this % eliminated(i)) retained_position(i) = retained_position(i) + 1
    end do
    retained_of_selection = pack(retained_position(selected), .not. this % eliminated(selected))
    if (size(retained_of_selection) < 1) then
       error stop 'elimination: a restriction retains an unknown at least'
    end if
    this % eliminated = this % eliminated(selected)
    if (allocated(this % inner)) call this % inner % restrict(retained_of_selection)

  end subroutine elimination_restrict

  !===================================================================!
  ! The stored partition released before a statement is read again.
  !===================================================================!

  subroutine clear_partition(this)

    class(elimination), intent(inout) :: this

    if (allocated(this % retained_at))   deallocate(this % retained_at)
    if (allocated(this % eliminated_at)) deallocate(this % eliminated_at)
    if (allocated(this % ke_row))    deallocate(this % ke_row, this % ke_column, this % ke_weight)
    if (allocated(this % ek_row))    deallocate(this % ek_row, this % ek_column, this % ek_weight)
    if (allocated(this % first_n))   deallocate(this % first_n)
    if (allocated(this % n_column))  deallocate(this % n_column, this % n_weight)
    if (allocated(this % order))     deallocate(this % order)
    if (allocated(this % diagonal))  deallocate(this % diagonal)

  end subroutine clear_partition

  !===================================================================!
  ! The permutation gathering triples by row, and the first triple of
  ! every row, first(nrows + 1) one past the last.
  !===================================================================!

  pure subroutine by_rows(nrows, row, by_row, first)

    integer, intent(in) :: nrows, row(:)
    integer, allocatable, intent(out) :: by_row(:), first(:)

    integer, allocatable :: next(:)
    integer :: e, i

    allocate(first(nrows + 1), by_row(size(row)))
    first = 0
    do e = 1, size(row)
       first(row(e) + 1) = first(row(e) + 1) + 1
    end do
    first(1) = 1
    do i = 1, nrows
       first(i + 1) = first(i) + first(i + 1)
    end do
    next = first(1:nrows)
    do e = 1, size(row)
       by_row(next(row(e))) = e
       next(row(e)) = next(row(e)) + 1
    end do

  end subroutine by_rows

  !===================================================================!
  ! An order of the eliminated rows with every row after the rows it
  ! reads: rows reading none first, each row released once all it
  ! reads are ordered. A cycle leaves rows unordered, and is refused.
  !===================================================================!

  pure subroutine ordered(ne, first_n, n_column, order)

    integer, intent(in) :: ne, first_n(:), n_column(:)
    integer, allocatable, intent(out) :: order(:)

    integer, allocatable :: pending(:), readers_first(:), readers(:), reads(:)
    integer :: e, j, placed, at, c

    allocate(order(ne), pending(ne), reads(size(n_column)))
    do e = 1, ne
       pending(e) = first_n(e + 1) - first_n(e)
    end do
    ! the readers of every eliminated row, by the row read
    do e = 1, ne
       reads(first_n(e):first_n(e + 1) - 1) = e
    end do
    call by_rows(ne, n_column, readers, readers_first)
    placed = 0
    do e = 1, ne
       if (pending(e) == 0) then
          placed = placed + 1
          order(placed) = e
       end if
    end do
    at = 0
    do while (at < placed)
       at = at + 1
       c = order(at)
       do j = readers_first(c), readers_first(c + 1) - 1
          e = reads(readers(j))
          pending(e) = pending(e) - 1
          if (pending(e) == 0) then
             placed = placed + 1
             order(placed) = e
          end if
       end do
    end do
    if (placed /= ne) then
       error stop 'elimination: the eliminated rows read one another in a cycle; &
            &state one of their kinds as rows'
    end if

  end subroutine ordered

  !===================================================================!
  ! y <- (I + N)^-1 y by substitution in the stored order.
  !===================================================================!

  pure subroutine substituted(this, y)

    class(elimination), intent(in)    :: this
    real(dp)          , intent(inout) :: y(:)

    integer :: i, e, j

    do i = 1, size(this % order)
       e = this % order(i)
       do j = this % first_n(e), this % first_n(e + 1) - 1
          y(e) = y(e) - this % n_weight(j) * y(this % n_column(j))
       end do
    end do

  end subroutine substituted

  !===================================================================!
  ! Solve the complement for the retained unknowns, then read the
  ! eliminated ones from their rows. The residual reported is that of
  ! the whole system.
  !===================================================================!

  subroutine elimination_solve(this, rhs, x, achieved)

    class(elimination), intent(inout) :: this
    real(dp)          , intent(in)    :: rhs(:)
    real(dp)          , intent(inout) :: x(:)
    real(dp)          , intent(out)   :: achieved

    real(dp), allocatable :: rhs_retained(:), x_retained(:), x_eliminated(:), r(:)
    real(dp) :: inner_achieved
    integer :: e
    type(solve_result) :: outcome

    if (size(x) /= size(rhs)) then
       error stop 'elimination: solution size matches rhs'
    end if
    if (size(rhs) /= size(this % eliminated)) then
       error stop 'elimination: the right-hand side is over the stated unknowns'
    end if

    call this % initialize_residual_history()
    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    call this % record_residual_norm(achieved)

    ! a statement refused within its storage limit stored no
    ! partition: the unknowns are unchanged and the refusal reported
    if (this % storage_exceeded) then
       call this % record_result(achieved, 0, SOLVE_STORAGE_EXCEEDED)
       return
    end if

    ! b_K - J_KE (I + N)^-1 b_E, the eliminated rows divided by their diagonals
    x_eliminated = rhs(this % eliminated_at) / this % diagonal
    call substituted(this, x_eliminated)
    rhs_retained = rhs(this % retained_at)
    do e = 1, size(this % ke_row)
       rhs_retained(this % ke_row(e)) = rhs_retained(this % ke_row(e)) &
            & - this % ke_weight(e) * x_eliminated(this % ke_column(e))
    end do

    x_retained = x(this % retained_at)
    call this % inner % solve(rhs_retained, x_retained, inner_achieved)
    outcome = this % inner % result()

    ! a singular or overflowing inner solve yields no step; its report
    ! is passed on as the residual
    if (.not. ieee_is_finite(inner_achieved) .or. outcome % failed()) then
       call this % record_result(achieved, outcome % iterations, SOLVE_INNER_FAILED)
       return
    end if
    if (inner_achieved > huge(1.0_dp) / 2.0_dp) then
       call this % record_result(achieved, outcome % iterations, SOLVE_INNER_FAILED)
       return
    end if

    ! x_E = (I + N)^-1 (b_E - J_EK x_K)
    x_eliminated = rhs(this % eliminated_at) / this % diagonal
    do e = 1, size(this % ek_row)
       x_eliminated(this % ek_row(e)) = x_eliminated(this % ek_row(e)) &
            & - this % ek_weight(e) * x_retained(this % ek_column(e))
    end do
    call substituted(this, x_eliminated)

    x(this % retained_at)   = x_retained
    x(this % eliminated_at) = x_eliminated

    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    call this % record_result(achieved, outcome % iterations)
    outcome = this % result()
    if (outcome % reason == SOLVE_CONTINUE) then
       call this % record_result(achieved, outcome % iterations, SOLVE_EXHAUSTED)
    end if

  end subroutine elimination_solve

end module operation_elimination
