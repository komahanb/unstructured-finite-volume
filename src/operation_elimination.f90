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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_elimination

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use util_precision        , only : dp
  use graph_fractal         , only : graph
  use view_directed         , only : directed_graph
  use view_directed_stored  , only : stored_directed_graph
  use operation_action      , only : operation
  use operation_minimization, only : minimizer, state, restrict
  use operation_minimization, only : solve_result, SOLVE_INNER_FAILED, SOLVE_EXHAUSTED, SOLVE_CONTINUE
  use operation_stencil     , only : stencil, combine_triples
  use field_stored          , only : stored_field

  implicit none

  private
  public :: elimination

  type, extends(minimizer) :: elimination

     ! one flag per unknown of the stated system: true for a row
     ! eliminated before the inner solve
     logical, allocatable :: eliminated(:)

     class(minimizer), allocatable :: inner

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

  !===================================================================!
  ! State the whole system, split its triples by set, divide the
  ! eliminated rows by their diagonals, assemble the complement and
  ! state the inner minimizer on it.
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
    integer , allocatable :: first_ek(:), by_row(:), rk(:), ck(:)
    integer , allocatable :: n_row(:), first_m(:), count_m(:), m_column(:)
    integer , allocatable :: column_row(:), active_columns(:)
    real(dp), allocatable :: weights(:), kk_weight(:), wk(:), m_weight(:), row_weight(:)
    type(stencil) :: complement
    type(stored_directed_graph) :: retained
    integer :: n, nk, ne, e, r, c, i, j, k, nkk, nke, nek, nn, m, total, num_active_columns

    call state(this, action, context, unknown_domain, num_unknowns, &
         & num_components, coupling, stored_inputs)
    call clear_partition(this)

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
       call action % entries(rows, columns, weights)
    class default
       error stop 'elimination: the complement is read from the explicit tangent, and this &
            &statement is a matrix-vector product without one'
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
    if (nk < 1) then
       error stop 'elimination: an unknown is retained at least'
    end if
    this % retained_at   = pack([(i, i = 1, n)], .not. this % eliminated)
    this % eliminated_at = pack([(i, i = 1, n)],       this % eliminated)

    ! split the triples by set
    m = size(rows)
    allocate(kk_row(m), kk_column(m), kk_weight(m))
    allocate(this % ke_row(m), this % ke_column(m), this % ke_weight(m))
    allocate(this % ek_row(m), this % ek_column(m), this % ek_weight(m))
    allocate(n_row(m), this % n_column(m), this % n_weight(m))
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
    this % ke_row    = this % ke_row(1:nke)
    this % ke_column = this % ke_column(1:nke)
    this % ke_weight = this % ke_weight(1:nke)
    this % ek_row    = this % ek_row(1:nek)
    this % ek_column = this % ek_column(1:nek)
    this % ek_weight = this % ek_weight(1:nek)

    ! N and J_EK by row
    call by_rows(ne, n_row(1:nn), by_row, this % first_n)
    this % n_column = this % n_column(by_row)
    this % n_weight = this % n_weight(by_row)
    call by_rows(ne, this % ek_row, by_row, first_ek)

    ! the substitution order: a row after every eliminated row it reads
    call ordered(ne, this % first_n, this % n_column, this % order)

    ! M = (I + N)^-1 J_EK row by row in that order. Sum equal
    ! retained columns before another row reads this one: storage
    ! then counts the nonzero coefficients of M, not dependency paths.
    allocate(first_m(ne), count_m(ne))
    allocate(column_row(nk), source=0)
    allocate(active_columns(nk), row_weight(nk), m_column(0), m_weight(0))
    total = 0
    do i = 1, ne
       e = this % order(i)
       num_active_columns = 0
       do j = first_ek(e), first_ek(e + 1) - 1
          call accumulated(this % ek_column(by_row(j)), this % ek_weight(by_row(j)))
       end do
       do j = this % first_n(e), this % first_n(e + 1) - 1
          c = this % n_column(j)
          do k = first_m(c), first_m(c) + count_m(c) - 1
             call accumulated(m_column(k), -this % n_weight(j) * m_weight(k))
          end do
       end do
       first_m(e) = total + 1
       call reserve_coefficients(total + num_active_columns, total, m_column, m_weight)
       do j = 1, num_active_columns
          c = active_columns(j)
          if (row_weight(c) == 0.0_dp) cycle
          total = total + 1
          m_column(total) = c
          m_weight(total) = row_weight(c)
       end do
       count_m(e) = total - first_m(e) + 1
    end do

    ! the complement's triples: J_KK, then minus J_KE M
    m = nkk
    do e = 1, nke
       m = m + count_m(this % ke_column(e))
    end do
    allocate(rk(m), ck(m), wk(m))
    rk(1:nkk) = kk_row(1:nkk)
    ck(1:nkk) = kk_column(1:nkk)
    wk(1:nkk) = kk_weight(1:nkk)
    m = nkk
    do e = 1, nke
       c = this % ke_column(e)
       rk(m + 1:m + count_m(c)) = this % ke_row(e)
       ck(m + 1:m + count_m(c)) = m_column(first_m(c):first_m(c) + count_m(c) - 1)
       wk(m + 1:m + count_m(c)) = -this % ke_weight(e) * m_weight(first_m(c):first_m(c) + count_m(c) - 1)
       m = m + count_m(c)
    end do
    call combine_triples(nk, nk, rk, ck, wk, rows, columns, weights)

    complement = stencil(rows, columns, weights, spread(0.0_dp, 1, nk), 'eliminated tangent')
    call complement % versioned(action % version(), action % transpose_version())
    retained = stored_directed_graph(nk, tails=[integer ::], heads=[integer ::])
    call this % inner % state(complement, complement % pattern, retained % vertex_set(), nk, &
         & num_components = 1, coupling = complement % pattern)

  contains

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

  end subroutine elimination_state

  subroutine reserve_coefficients(required, used, columns, weights)
    integer, intent(in) :: required, used
    integer, allocatable, intent(inout) :: columns(:)
    real(dp), allocatable, intent(inout) :: weights(:)
    integer, allocatable :: larger_columns(:)
    real(dp), allocatable :: larger_weights(:)
    integer :: capacity

    capacity = size(columns)
    if (required <= capacity) return
    if (capacity <= huge(capacity) / 2) capacity = 2 * capacity
    capacity = max(required, capacity)
    allocate(larger_columns(capacity), larger_weights(capacity))
    larger_columns(1:used) = columns(1:used)
    larger_weights(1:used) = weights(1:used)
    call move_alloc(larger_columns, columns)
    call move_alloc(larger_weights, weights)
  end subroutine reserve_coefficients

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
