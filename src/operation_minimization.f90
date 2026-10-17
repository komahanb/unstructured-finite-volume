!=====================================================================!
! THE DISCRETE MINIMIZER M_h of a discrete Lagrangian L_h.
!
! L_h is solved block by block along its chain of instants.
! M_h pairs the rows of every block with a minimizer.
! The pairing is like that of a condition with its multiplier in L.
! The minimizer of a block may contain inner minimizers.
! tree(M_h) prints the whole structure:
!
!      (march)            over the blocks, forward for the state,
!         |               in reverse for the multipliers
!      [R_b] - (newton)   the rows of block b
!                 |
!              [J_b] - (elimination)       the jacobian of the rows
!                         |
!                      [S_b] - (dense direct)   the rows that remain
!                                               after the elimination
!
! THE PAIRING. The minimizer stated for M_h is copied to every block.
! Each copy reads from the block's rows what it needs.
! The elimination reads the tying rows: the rows that define a
! derivative from other unknowns, such as a velocity from the
! positions by the time scheme.
! Rows without tying rows are refused, and the reason is named.
! Where no minimizer is stated, the rows of each block choose the
! default (default_minimizer).
!
! THE MINIMAND is the data M_h reads from L_h.
! It is an abstract type, and the discrete Lagrangian extends it.
! It states the blocks, and the earlier blocks each one reads.
! It gives the rows of a block at the state.
! It gives the right sides of every order of the expansion and of the
! adjoint.
! It places the multipliers and assembles the solution.
! M_h reads nothing else of L_h.
!
! THE FACTORS OF A BLOCK. The jacobian of the rows of a block at its
! solution is formed once.
! It serves every order of the expansion along every direction.
! Its transpose is formed once, and serves every order of the adjoint.
! A direct solve factorises each of them once.
! Every later right side is substituted against those factors.
! Over several images the rows of a block are cut into parts.
! Then each right side is solved on the rows linearized at the
! solution, part by part.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_minimization

  use util_precision             , only : dp
  use util_verbosity             , only : verbosity
  use operation_minimizer        , only : minimizer, solve_result
  use operation_newton           , only : newton
  use operation_dense_direct     , only : dense_direct
  use operation_elimination      , only : elimination
  use operation_gmres            , only : gmres
  use operation_gauss_seidel     , only : gauss_seidel
  use operation_stencil          , only : stencil
  use operation_exchange         , only : exchange
  use operation_residual_operator, only : residual_operator, chain_block
  use operation_field            , only : discrete_field
  use field_stored               , only : stored_field
  use operation_driver           , only : driver
  use view_read_write            , only : bipartite_digraph, FIRST_PART, SECOND_PART

  implicit none

  private
  public :: minimand, march_state, block_facts, discrete_minimizer
  public :: default_minimizer, minimizer_of, system_solution, rows_evaluated

  !===================================================================!
  ! THE STATE OF A MARCH: the data the blocks read and write - the
  ! tuples of every column at the instants and the stages, and the
  ! multipliers - as the minimand lays them out. The minimization
  ! passes it and reads nothing in it.
  !===================================================================!

  type, abstract :: march_state
  end type march_state

  !===================================================================!
  ! THE FACTS A BLOCK'S DEFAULT MINIMIZER IS CHOSEN FROM:
  ! whether the manifold has space;
  ! whether a block is the whole line at one moment;
  ! whether the march continues from the instant before each block;
  ! whether the rows have tying rows;
  ! the number of cells, and the components of a point.
  !===================================================================!

  type :: block_facts
     logical :: with_space = .false.
     logical :: whole_line = .false.
     logical :: continued  = .false.
     logical :: tying      = .false.
     integer :: cells      = 1
     integer :: width      = 1
  end type block_facts

  !===================================================================!
  ! THE MINIMAND, the data of a minimization.
  !===================================================================!

  type, abstract :: minimand
   contains
     ! the blocks of the march, the blocks whose instants block b reads
     ! as its history, and the line a block is announced by, empty for
     ! none
     procedure(count_of_blocks)   , deferred :: num_blocks
     procedure(blocks_read)       , deferred :: history_of
     procedure(announcement)      , deferred :: block_label
     ! the facts a default minimizer is chosen from, and the image owning
     ! each unknown of a block's rows
     procedure(facts_of_blocks)   , deferred :: facts
     procedure(owners_of_rows)    , deferred :: unknown_owners
     ! the state from the estimate; block b's rows at the state, its
     ! estimate, its design and its owners - at the solution, with the
     ! jets of every column, where reversed; a column of the block
     ! written into the state
     procedure(state_of_estimate) , deferred :: begun
     procedure(rows_at_state)     , deferred :: block_system
     procedure(column_written)    , deferred :: written
     ! the columns of the expansion, (d, m, column), and the right side
     ! of order m along d at a block, the lower orders read from the
     ! state
     procedure(columns_of_expansion), deferred :: expansion_columns
     procedure(expansion_side)    , deferred :: expansion_right_side
     ! the solution of the march, the state at every point
     procedure(solution_of_state) , deferred :: march_solution
     ! whether the equations are paired with their multiplier, and the
     ! adjoint: its data from the solution, the right side of order o
     ! along d at a block, the multiplier of that order received, the
     ! multipliers of the block placed, the solution with them
     procedure(pairing_stated)    , deferred :: equation_pairing
     procedure(adjoint_of_solution), deferred :: adjoint_begun
     procedure(adjoint_side)      , deferred :: adjoint_right_side
     procedure(multiplier_received), deferred :: adjoint_received
     procedure(multipliers_placed), deferred :: adjoint_placed
     procedure(solution_with_multipliers), deferred :: adjoint_solution
  end type minimand

  abstract interface

     integer function count_of_blocks(this)
       import :: minimand
       class(minimand), intent(in) :: this
     end function count_of_blocks

     function blocks_read(this, b) result(read)
       import :: minimand
       class(minimand), intent(in) :: this
       integer        , intent(in) :: b
       integer, allocatable :: read(:)
     end function blocks_read

     function announcement(this, b) result(line)
       import :: minimand
       class(minimand), intent(in) :: this
       integer        , intent(in) :: b
       character(len=:), allocatable :: line
     end function announcement

     function facts_of_blocks(this) result(facts)
       import :: minimand, block_facts
       class(minimand), intent(in) :: this
       type(block_facts) :: facts
     end function facts_of_blocks

     function owners_of_rows(this, unknowns, stride, npts) result(owner)
       import :: minimand
       class(minimand), intent(in) :: this
       integer        , intent(in) :: unknowns, stride, npts
       integer, allocatable :: owner(:)
     end function owners_of_rows

     subroutine state_of_estimate(this, estimate, state)
       import :: minimand, march_state, discrete_field
       class(minimand)                  , intent(in)  :: this
       type(discrete_field)             , intent(in)  :: estimate
       class(march_state), allocatable  , intent(out) :: state
     end subroutine state_of_estimate

     subroutine rows_at_state(this, b, state, reversed, blk)
       import :: minimand, march_state, chain_block
       class(minimand)   , intent(in)    :: this
       integer           , intent(in)    :: b
       class(march_state), intent(inout) :: state
       logical           , intent(in)    :: reversed
       type(chain_block) , intent(out)   :: blk
     end subroutine rows_at_state

     subroutine column_written(this, blk, x, column, state)
       import :: minimand, march_state, chain_block, dp
       class(minimand)   , intent(in)    :: this
       type(chain_block) , intent(in)    :: blk
       real(dp)          , intent(in)    :: x(:)
       integer           , intent(in)    :: column
       class(march_state), intent(inout) :: state
     end subroutine column_written

     function columns_of_expansion(this) result(columns)
       import :: minimand
       class(minimand), intent(in) :: this
       integer, allocatable :: columns(:,:)
     end function columns_of_expansion

     subroutine expansion_side(this, blk, d, m, state, rhs)
       import :: minimand, march_state, chain_block, dp
       class(minimand)      , intent(in)    :: this
       type(chain_block)    , intent(in)    :: blk
       integer              , intent(in)    :: d, m
       class(march_state)   , intent(inout) :: state
       real(dp), allocatable, intent(out)   :: rhs(:)
     end subroutine expansion_side

     subroutine solution_of_state(this, state, solution)
       import :: minimand, march_state, discrete_field
       class(minimand)     , intent(in)    :: this
       class(march_state)  , intent(inout) :: state
       type(discrete_field), intent(out)   :: solution
     end subroutine solution_of_state

     logical function pairing_stated(this)
       import :: minimand
       class(minimand), intent(in) :: this
     end function pairing_stated

     subroutine adjoint_of_solution(this, solution, state)
       import :: minimand, march_state, discrete_field
       class(minimand)     , intent(in)    :: this
       type(discrete_field), intent(in)    :: solution
       class(march_state)  , intent(inout) :: state
     end subroutine adjoint_of_solution

     subroutine adjoint_side(this, blk, o, d, state, rhs)
       import :: minimand, march_state, chain_block, dp
       class(minimand)      , intent(in)    :: this
       type(chain_block)    , intent(in)    :: blk
       integer              , intent(in)    :: o, d
       class(march_state)   , intent(inout) :: state
       real(dp), allocatable, intent(out)   :: rhs(:)
     end subroutine adjoint_side

     subroutine multiplier_received(this, blk, o, d, x, state)
       import :: minimand, march_state, chain_block, dp
       class(minimand)   , intent(in)    :: this
       type(chain_block) , intent(in)    :: blk
       integer           , intent(in)    :: o, d
       real(dp)          , intent(in)    :: x(:)
       class(march_state), intent(inout) :: state
     end subroutine multiplier_received

     subroutine multipliers_placed(this, blk, state)
       import :: minimand, march_state, chain_block
       class(minimand)   , intent(in)    :: this
       type(chain_block) , intent(inout) :: blk
       class(march_state), intent(inout) :: state
     end subroutine multipliers_placed

     subroutine solution_with_multipliers(this, state, solution)
       import :: minimand, march_state, discrete_field
       class(minimand)     , intent(in)    :: this
       class(march_state)  , intent(inout) :: state
       type(discrete_field), intent(inout) :: solution
     end subroutine solution_with_multipliers

  end interface

  !===================================================================!
  ! THE DISCRETE MINIMIZER: the minimizer stated for the rows of every
  ! block, or none, where the rows of each block choose the default.
  ! The march visits the blocks in the order of the digraph of the
  ! blocks and the data each one writes and reads.
  !===================================================================!

  type :: discrete_minimizer
     class(minimizer) , allocatable :: block
     type(driver)     , allocatable :: march
     type(block_facts) :: facts
     integer :: blocks = 0
   contains
     procedure :: minimize => minimization_minimize
     procedure :: reversed => minimization_reversed
     procedure :: tree     => minimizer_tree
     procedure, private :: block_solved
  end type discrete_minimizer

  interface discrete_minimizer
     module procedure discrete_minimizer_of
  end interface discrete_minimizer

  ! the systems stated for their factors to be retained, numbered from one:
  ! a direct minimizer retains a factorisation for the system of its
  ! own number alone
  integer, save :: systems_stated = 0

contains

  !===================================================================!
  ! THE DISCRETE MINIMIZER OF A MINIMAND, with the minimizer stated for
  ! the rows of every block, or with none for the default.
  ! The march has one rule per block.
  ! Rule b writes the datum of block b.
  ! It reads the data of the blocks in the history of block b.
  !===================================================================!

  function discrete_minimizer_of(data, block) result(this)

    class(minimand) , intent(in)           :: data
    class(minimizer), intent(in), optional :: block
    type(discrete_minimizer) :: this

    integer, allocatable :: from_part(:), from_vertex(:), to_part(:), to_vertex(:), read(:)
    type(bipartite_digraph) :: over
    type(newton) :: plain
    integer :: b

    if (present(block)) allocate(this % block, source=block)
    this % facts  = data % facts()
    this % blocks = data % num_blocks()
    allocate(from_part(0), from_vertex(0), to_part(0), to_vertex(0))
    do b = 1, this % blocks
       read        = data % history_of(b)
       from_part   = [from_part, FIRST_PART, spread(SECOND_PART, 1, size(read))]
       from_vertex = [from_vertex, b, read]
       to_part     = [to_part, SECOND_PART, spread(FIRST_PART, 1, size(read))]
       to_vertex   = [to_vertex, b, spread(b, 1, size(read))]
    end do
    over = bipartite_digraph(this % blocks, this % blocks, from_part, from_vertex, to_part, to_vertex)
    if (allocated(this % block)) then
       this % march = driver(this % block, over)
    else
       this % march = driver(plain, over)
    end if

  end function discrete_minimizer_of

  !===================================================================!
  ! THE MARCH. Forward, each block's rows at the state solved by its
  ! minimizer, the solution written, and each column of the expansion
  ! solved against the factors of the jacobian at the solution; the
  ! solution assembled; with the equations paired, the adjoint in
  ! reverse (reversed) and the solution completed with the
  ! multipliers. converged present, a block whose rows do not converge
  ! returns it false at once; state present, the state of the march.
  !===================================================================!

  subroutine minimization_minimize(this, data, estimate, solution, converged, state)

    class(discrete_minimizer), intent(in)  :: this
    class(minimand)     , intent(in)  :: data
    type(discrete_field), intent(in)  :: estimate
    type(discrete_field), intent(out) :: solution
    logical             , intent(out), optional :: converged
    class(march_state)  , allocatable, intent(out), optional :: state

    class(march_state), allocatable :: march
    type(chain_block) :: blk
    class(minimizer), allocatable :: system
    character(len=:), allocatable :: label
    integer, allocatable :: order(:), columns(:,:)
    real(dp), allocatable :: rhs(:), x(:)
    logical :: met, parts
    real(dp) :: tolerance
    integer :: i, b, j

    if (present(converged)) converged = .true.
    call data % begun(estimate, march)
    columns = data % expansion_columns()
    parts   = rows_in_parts(data)

    order = this % march % visits()
    do i = 1, size(order)
       b = order(i)
       label = data % block_label(b)
       if (verbosity >= 1 .and. len(label) > 0) print '(a)', label
       call data % block_system(b, march, .false., blk)
       call this % block_solved(data, blk, met, size(columns, 2) > 0)
       if (.not. met) then
          if (present(converged)) then
             converged = .false.
             if (present(state)) call move_alloc(march, state)
             return
          end if
          error stop 'operation_minimization: minimize did not converge on the rows of a block'
       end if
       call data % written(blk, blk % q, 0, march)
       if (size(columns, 2) == 0) cycle
       if (.not. parts) call factorised(this, data, blk, .false., system, tolerance)
       do j = 1, size(columns, 2)
          if (verbosity >= 1) print '(a,i0,a,i0)', 'expansion ', columns(1, j), ' along the design, order ', &
               & columns(2, j)
          call data % expansion_right_side(blk, columns(1, j), columns(2, j), march, rhs)
          if (parts) then
             call linearized_solution(this, data, blk, .false., rhs, x)
          else
             call substituted(system, rhs, tolerance, x)
          end if
          call data % written(blk, x, columns(3, j), march)
       end do
    end do
    call data % march_solution(march, solution)

    if (data % equation_pairing()) then
       call this % reversed(data, solution, march)
       call data % adjoint_solution(march, solution)
    end if
    if (present(state)) call move_alloc(march, state)

  end subroutine minimization_minimize

  !===================================================================!
  ! THE ADJOINT of a march, on its state: in reverse, each block's
  ! transposed rows at the solution, every order of the adjoint along
  ! every direction the march expanded - column zero first - solved
  ! against their factors, the multipliers of each order received and
  ! the block's placed.
  !===================================================================!

  subroutine minimization_reversed(this, data, solution, state)

    class(discrete_minimizer), intent(in)    :: this
    class(minimand)     , intent(in)    :: data
    type(discrete_field), intent(in)    :: solution
    class(march_state)  , intent(inout) :: state

    type(chain_block) :: blk
    class(minimizer), allocatable :: system
    character(len=:), allocatable :: label
    integer, allocatable :: order(:), columns(:,:), orders(:,:)
    real(dp), allocatable :: rhs(:), x(:)
    logical :: parts
    real(dp) :: tolerance
    integer :: i, b, j

    parts   = rows_in_parts(data)
    columns = data % expansion_columns()
    allocate(orders(3, 1 + size(columns, 2)))
    orders(:, 1)  = [1, 0, 0]
    orders(:, 2:) = columns
    call data % adjoint_begun(solution, state)
    order = this % march % visits()
    do i = size(order), 1, -1
       b = order(i)
       label = data % block_label(b)
       if (verbosity >= 1 .and. len(label) > 0) print '(a)', 'adjoint  ' // label
       call data % block_system(b, state, .true., blk)
       call inputs_ensured(this, data, blk)
       if (.not. parts) call factorised(this, data, blk, .true., system, tolerance)
       do j = 1, size(orders, 2)
          if (verbosity >= 1 .and. orders(2, j) > 0) print '(a,i0,a,i0)', 'adjoint expansion ', orders(1, j), &
               & ' along the design, order ', orders(2, j)
          call data % adjoint_right_side(blk, orders(2, j), orders(1, j), state, rhs)
          if (parts) then
             call linearized_solution(this, data, blk, .true., rhs, x)
          else
             call substituted(system, rhs, tolerance, x)
          end if
          call data % adjoint_received(blk, orders(2, j), orders(1, j), x, state)
       end do
       call data % adjoint_placed(blk, state)
    end do

  end subroutine minimization_reversed

  !===================================================================!
  ! THE ROWS OF A BLOCK SOLVED by the minimizer of the block, from the
  ! estimate of the block.
  ! The minimizer is the one stated, or the default of the rows.
  ! expanded: the inputs of the rows at the solution are stored in the
  ! block, and the jacobian there is formed from them.
  !===================================================================!

  subroutine block_solved(this, data, blk, converged, expanded)

    class(discrete_minimizer), intent(in)    :: this
    class(minimand)          , intent(in)    :: data
    type(chain_block)        , intent(inout) :: blk
    logical                  , intent(out)   :: converged
    logical                  , intent(in)    :: expanded

    class(minimizer), allocatable :: M_b

    call minimizer_of(blk % rows, data % facts(), M_b, this)
    if (rows_in_parts(data)) then
       call system_solution(M_b, blk % rows, blk % design, blk % q, converged, &
            & owners=data % unknown_owners(size(blk % q), blk % rows % degrees, size(blk % rows % at)))
    else
       call system_solution(M_b, blk % rows, blk % design, blk % q, converged)
    end if
    if (converged .and. expanded) call inputs_ensured(this, data, blk)

  end subroutine block_solved

  ! whether the rows of a block are cut into parts over the images:
  ! several images sharing the cells of a region a block does not read
  ! as one line
  logical function rows_in_parts(data)
    class(minimand), intent(in) :: data
    type(block_facts) :: facts
    facts = data % facts()
    rows_in_parts = num_images() > 1 .and. facts % with_space .and. .not. facts % whole_line
  end function rows_in_parts

  !===================================================================!
  ! THE MINIMIZER OF ROWS: a copy of the minimizer of by, paired with
  ! the rows. Without by, or where by states no minimizer, the default
  ! the rows choose.
  !===================================================================!

  subroutine minimizer_of(rows, facts, M_b, by)

    type(residual_operator)      , intent(in)  :: rows
    type(block_facts)            , intent(in)  :: facts
    class(minimizer), allocatable, intent(out) :: M_b
    type(discrete_minimizer)     , intent(in), optional :: by

    if (present(by)) then
       if (allocated(by % block)) then
          allocate(M_b, source=by % block)
          call paired_with(M_b, rows, facts)
          return
       end if
    end if
    call default_minimizer(facts, M_b, rows)

  end subroutine minimizer_of

  !===================================================================!
  ! THE PAIRING OF A MINIMIZER WITH ROWS: each minimizer in it reads
  ! from the rows what it needs.
  ! The elimination removes the tying rows of the rows.
  ! Rows without tying rows are refused, and the reason is named.
  ! A block of the Gauss-Seidel sweeps is a point of the rows.
  ! Where the rows couple every point of a moment through the line, a
  ! block of the sweeps is the whole line at one moment.
  !===================================================================!

  recursive subroutine paired_with(M_b, rows, facts)

    class(minimizer)       , intent(inout) :: M_b
    type(residual_operator), intent(in)    :: rows
    type(block_facts)      , intent(in)    :: facts

    select type (M_b)
    type is (newton)
       if (allocated(M_b % inner)) call paired_with(M_b % inner, rows, facts)
    type is (elimination)
       if (.not. allocated(rows % tying)) then
          error stop 'operation_minimization: an elimination removes the tying rows of the rows it is paired with, &
               &and these rows state none'
       end if
       M_b % eliminated = rows % tying
       if (allocated(M_b % inner)) call paired_with(M_b % inner, rows, facts)
    type is (gmres)
       if (allocated(M_b % preconditioner)) call paired_with(M_b % preconditioner, rows, facts)
    type is (gauss_seidel)
       M_b % block_width = rows % degrees
       if (facts % whole_line) M_b % block_width = rows % degrees * facts % cells
    end select

  end subroutine paired_with

  !===================================================================!
  ! THE DEFAULT MINIMIZER of the rows of a block.
  !
  ! Without space: Newton over a direct solve.
  ! Where the rows have tying rows, the elimination removes them first.
  ! Newton converges quadratically to round-off.
  ! The sensitivities read the converged state, so the tolerance is
  ! 1e-13.
  !
  ! With space: Newton over GMRES, restarted every 60 iterations.
  ! GMRES stops at a relative residual of 1e-3, the constant forcing
  ! term of Dembo, Eisenstat and Steihaug.
  ! GMRES is preconditioned by four block Gauss-Seidel sweeps.
  ! With two sweeps, GMRES(60) stagnates on the 32 x 32 box.
  ! Four sweeps solve that box in 26 s.
  ! A restart of 120 with two sweeps takes 27 s, at twice the storage.
  ! A block of the sweeps is the components of one point.
  ! It is the whole line at one moment where the rows couple the points
  ! of a moment through the line: the conditions at one end of an open
  ! line, or a law projected on the cells of a line.
  ! A march continued from the instant before begins near the solution,
  ! so its first residual is not the scale of the law.
  ! Its tolerance is 1e-13, the rounding of the residual.
  !
  ! With space and tying rows: the elimination, then a direct solve.
  ! These rows mix positions, rotations and reactions, of different
  ! units.
  ! The direct solve scales the rows, then the columns, to a largest
  ! entry of one before the LU factorisation (equilibrated).
  !
  ! Without rows, the default of the facts alone, not paired.
  !===================================================================!

  subroutine default_minimizer(facts, M_b, rows)

    type(block_facts)            , intent(in)  :: facts
    class(minimizer), allocatable, intent(out) :: M_b
    type(residual_operator)      , intent(in), optional :: rows

    type(block_facts) :: own

    own = facts
    if (present(rows)) then
       own % tying = allocated(rows % tying)
       own % width = rows % degrees
    end if
    if (own % with_space) then
       call require_images()
       if (own % tying) then
          M_b = newton(elimination(dense_direct(equilibrated=.true.)), tolerance=1.0e-13_dp)
       else if (own % continued) then
          M_b = newton(gmres(gauss_seidel(4), restart=60, tolerance=1.0e-3_dp), tolerance=1.0e-13_dp)
       else
          M_b = newton(gmres(gauss_seidel(4), restart=60, tolerance=1.0e-3_dp))
       end if
    else if (own % tying) then
       M_b = newton(elimination(dense_direct()), tolerance=1.0e-13_dp)
    else
       M_b = newton(dense_direct(), tolerance=1.0e-13_dp)
    end if
    if (present(rows)) call paired_with(M_b, rows, own)

  end subroutine default_minimizer

  !===================================================================!
  ! A SYSTEM SOLVED BY ITS MINIMIZER M_b, from q.
  ! M_b is stated on the rows, with their design as the stored input. With owners, over
  ! several images each image solves its part of the rows, the
  ! partitioner's cut of the rows' coupling, and the parts' owned
  ! entries are assembled whole. inputs present, the rows' inputs at
  ! the solution, from which their jacobian there is formed.
  !===================================================================!

  subroutine system_solution(M_b, rows, design, q, converged, owners, inputs)

    class(minimizer)       , intent(in)    :: M_b
    type(residual_operator), intent(in)    :: rows
    type(stored_field)     , intent(in)    :: design
    real(dp)               , intent(inout) :: q(:)
    logical                , intent(out)   :: converged
    integer                , intent(in), optional :: owners(:)
    type(stored_field), allocatable, intent(out), optional :: inputs(:)

    class(minimizer), allocatable :: solver
    type(solve_result) :: outcome
    type(exchange) :: whole, part
    type(residual_operator) :: sub
    real(dp), allocatable :: local(:), y(:)
    real(dp) :: achieved

    allocate(solver, source=M_b)
    if (present(owners)) then
       whole = exchange(owners, width=rows % degrees)
       part  = whole % part(rows % coupling())
       sub   = rows % part_of(part)
       select type (solver)
       type is (newton)
          solver % distribution = part
       end select
       call stated(solver, sub, design)
       local = part % restriction(q)
       call solver % solve(spread(0.0_dp, 1, size(local)), local, achieved)
       q = part % assembled(local)
    else
       call stated(solver, rows, design)
       call solver % solve(spread(0.0_dp, 1, size(q)), q, achieved)
    end if
    outcome   = solver % result()
    converged = outcome % converged()
    if (present(inputs) .and. converged) call rows_evaluated(rows, design, q, y, inputs)

  end subroutine system_solution

  ! the rows at q, y, and their inputs there, from which their jacobian
  ! at q is formed
  subroutine rows_evaluated(rows, design, q, y, inputs)
    type(residual_operator)        , intent(in)  :: rows
    type(stored_field)             , intent(in)  :: design
    real(dp)                       , intent(in)  :: q(:)
    real(dp), allocatable          , intent(out) :: y(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)
    type(newton) :: solver
    call stated(solver, rows, design)
    call solver % evaluate(q, y, inputs)
  end subroutine rows_evaluated

  ! a minimizer stated on rows: their unknowns, their design the stored
  ! input, their coupling the sweeps' colouring - after an elimination
  ! a direct solve has no sweep
  subroutine stated(solver, rows, design)
    class(minimizer)       , intent(inout) :: solver
    type(residual_operator), intent(in)    :: rows
    type(stored_field)     , intent(in)    :: design
    if (allocated(rows % tying)) then
       call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), rows % num_unknowns(), &
            & stored_inputs=[design])
    else
       call solver % state(rows, rows % unknown_graph(), rows % unknown_domain(), rows % num_unknowns(), &
            & stored_inputs=[design], coupling=rows % coupling())
    end if
  end subroutine stated

  ! the rows' inputs at the block's unknowns, where the block has none
  subroutine inputs_ensured(this, data, blk)
    class(discrete_minimizer), intent(in)    :: this
    class(minimand)          , intent(in)    :: data
    type(chain_block)        , intent(inout) :: blk
    real(dp), allocatable :: y(:)
    associate (u1 => this, u2 => data); end associate
    if (allocated(blk % inputs)) return
    call rows_evaluated(blk % rows, blk % design, blk % q, y, blk % inputs)
  end subroutine inputs_ensured

  !===================================================================!
  ! THE FACTORS OF A BLOCK: the jacobian of its rows at the solution -
  ! transposed for the adjoint - stated on the inner minimizer of the
  ! block's newton.
  ! The jacobian is numbered, so that a direct minimizer keeps its
  ! factorisation for every right side substituted after.
  ! tolerance is the one of the block's newton; the solves are refined
  ! to it.
  !===================================================================!

  subroutine factorised(this, data, blk, transposed, system, tolerance)

    class(discrete_minimizer)    , intent(in)    :: this
    class(minimand)              , intent(in)    :: data
    type(chain_block)            , intent(inout) :: blk
    logical                      , intent(in)    :: transposed
    class(minimizer), allocatable, intent(out)   :: system
    real(dp)                     , intent(out)   :: tolerance

    class(minimizer), allocatable :: M_b
    type(stencil) :: jacobian
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:)
    logical :: defined

    call inputs_ensured(this, data, blk)
    call blk % rows % explicit_jacobian(blk % rows % unknown_graph(), blk % rows % bind(blk % inputs), 1, &
         & rows, columns, weights, defined)
    if (.not. defined) then
       error stop 'operation_minimization: the rows of a block state no explicit jacobian to factorise'
    end if
    jacobian = stencil(rows, columns, weights, spread(0.0_dp, 1, size(blk % q)), 'explicit jacobian')
    if (transposed) call jacobian % reverse()
    systems_stated = systems_stated + 1
    call jacobian % versioned(systems_stated)
    call minimizer_of(blk % rows, data % facts(), M_b, this)
    tolerance = M_b % tolerance
    select type (M_b)
    type is (newton)
       call move_alloc(M_b % inner, system)
    class default
       call move_alloc(M_b, system)
    end select
    ! the coupling newton gives its inner minimizer: the rows' own,
    ! the jacobian's pattern after an elimination
    if (allocated(blk % rows % tying)) then
       call system % state(jacobian, jacobian % pattern, blk % rows % unknown_domain(), size(blk % q), &
            & coupling=jacobian % pattern)
    else
       call system % state(jacobian, jacobian % pattern, blk % rows % unknown_domain(), size(blk % q), &
            & coupling=blk % rows % coupling())
    end if

  end subroutine factorised

  !===================================================================!
  ! THE SOLUTION x OF THE SYSTEM STATED, J x = rhs, refined, x <- x +
  ! dx with J dx = rhs - J x by the system's own solve, until the
  ! imbalance is below the tolerance relative to the first, as newton
  ! drives a linear system: a direct system meets it at the first
  ! solve, an iterative one - GMRES to a relative residual of its own -
  ! in a few. A refinement that does not decrease the imbalance is at
  ! its rounding, and ends it.
  !===================================================================!

  subroutine substituted(system, rhs, tolerance, x)

    class(minimizer)     , intent(inout) :: system
    real(dp)             , intent(in)    :: rhs(:), tolerance
    real(dp), allocatable, intent(out)   :: x(:)

    type(solve_result) :: outcome
    real(dp), allocatable :: r(:), dx(:)
    real(dp) :: achieved, initial, imbalance, previous
    integer :: it

    allocate(x(size(rhs)), dx(size(rhs)), source=0.0_dp)
    call system % imbalance(rhs, x, r)
    initial   = system % norm(r)
    imbalance = initial
    do it = 1, system % max_iterations
       if (imbalance <= tolerance * initial) exit
       dx = 0.0_dp
       call system % solve(r, dx, achieved)
       outcome = system % result()
       if (outcome % failed()) then
          error stop 'operation_minimization: a block''s linear system was not solved: ' // outcome % description()
       end if
       x = x + dx
       previous = imbalance
       call system % imbalance(rhs, x, r)
       imbalance = system % norm(r)
       if (imbalance >= previous) exit
    end do

  end subroutine substituted

  ! the solution x of J x = rhs - transposed for the adjoint - by the
  ! minimizer of the block on its rows linearized at the solution, part
  ! by part over the images
  subroutine linearized_solution(this, data, blk, transposed, rhs, x)
    class(discrete_minimizer), intent(in)    :: this
    class(minimand)          , intent(in)    :: data
    type(chain_block)        , intent(inout) :: blk
    logical                  , intent(in)    :: transposed
    real(dp)                 , intent(in)    :: rhs(:)
    real(dp), allocatable    , intent(out)   :: x(:)
    type(residual_operator) :: lin
    class(minimizer), allocatable :: M_b
    logical :: met
    call inputs_ensured(this, data, blk)
    lin = blk % rows % linearize(blk % rows % unknown_graph(), blk % rows % bind(blk % inputs), rhs, &
         & transposed=transposed, version_number=blk % rows % version())
    call minimizer_of(lin, data % facts(), M_b, this)
    allocate(x(size(rhs)), source=0.0_dp)
    call system_solution(M_b, lin, blk % design, x, met, &
         & owners=data % unknown_owners(size(blk % q), blk % rows % degrees, size(blk % rows % at)))
    if (.not. met) error stop 'operation_minimization: a block''s linear system was not solved over its parts'
  end subroutine linearized_solution

  !===================================================================!
  ! THE TREE OF M_h: the march, then the minimizer of the rows of a
  ! block. Each inner minimizer is printed under the one it serves.
  ! The data each minimizer solves are printed in brackets above it.
  !===================================================================!

  function minimizer_tree(this) result(text)

    class(discrete_minimizer), intent(in) :: this
    character(len=:), allocatable :: text

    class(minimizer), allocatable :: M_b
    character(len=24) :: count

    write(count, '(i0)') this % blocks
    text = '(march) over ' // trim(count) // ' blocks, forward for the state, in reverse for the multipliers'
    if (allocated(this % block)) then
       text = text // minimizer_lines(this % block, 1, 'the rows of block b')
    else
       call default_minimizer(this % facts, M_b)
       text = text // minimizer_lines(M_b, 1, 'the rows of block b, the default minimizer')
    end if

  end function minimizer_tree

  ! the lines of M_b and of its inner minimizers, each line begun on a
  ! new line
  recursive function minimizer_lines(M_b, depth, data) result(text)

    class(minimizer), intent(in) :: M_b
    integer         , intent(in) :: depth
    character(len=*), intent(in) :: data
    character(len=:), allocatable :: text

    character(len=:), allocatable :: indent, line
    character(len=24) :: number

    indent = new_line('a') // repeat('   ', depth)
    text   = indent // '[' // data // ']'
    write(number, '(es9.2)') M_b % tolerance
    select type (M_b)
    type is (newton)
       line = '(newton) tolerance ' // trim(adjustl(number))
       if (M_b % chord) line = line // ', the chord method'
       text = text // indent // line
       if (allocated(M_b % inner)) text = text // minimizer_lines(M_b % inner, depth + 1, 'the jacobian of the rows')
    type is (elimination)
       text = text // indent // '(elimination) of the tying rows'
       if (allocated(M_b % inner)) text = text // minimizer_lines(M_b % inner, depth + 1, 'the rows that remain')
    type is (gmres)
       write(number, '(i0)') M_b % restart
       text = text // indent // '(gmres) restarted every ' // trim(number)
       if (allocated(M_b % preconditioner)) then
          text = text // minimizer_lines(M_b % preconditioner, depth + 1, 'the linear system')
       end if
    type is (gauss_seidel)
       write(number, '(i0)') M_b % max_iterations
       text = text // indent // '(gauss-seidel) ' // trim(number) // ' sweeps'
    type is (dense_direct)
       if (M_b % equilibrated) then
          text = text // indent // '(dense direct) rows, then columns, scaled to a largest entry of one'
       else
          text = text // indent // '(dense direct)'
       end if
    class default
       text = text // indent // '(' // M_b % name() // ')'
    end select

  end function minimizer_lines

  !===================================================================!
  ! A build of one image under a launcher of several would run the
  ! whole problem once per process, and is refused.
  !===================================================================!

  subroutine require_images()

    character(len=32) :: value
    integer :: status, ranks

    if (num_images() > 1) return
    ranks = 1
    call get_environment_variable('OMPI_COMM_WORLD_SIZE', value, status=status)
    if (status /= 0) call get_environment_variable('PMI_SIZE', value, status=status)
    if (status == 0) read(value, *, iostat=status) ranks
    if (ranks > 1) then
       error stop 'operation_minimization: this program was built for one image, but the launcher &
            &started several processes; run it without mpirun, or build it with COARRAY=lib'
    end if

  end subroutine require_images

end module operation_minimization
