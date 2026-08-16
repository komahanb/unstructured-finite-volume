!=====================================================================!
! TIME INTEGRATION TOWER . LEVEL 5 . FIELD CALCULUS
!
! The level answers one question: WHERE DO THE NUMERICAL VALUES
! LIVE. Three fields, three domains, and the domains are the answer:
!
!      q0   : Q -> reals  [2, 0]
!      time : T -> reals  [0, 1/2, 1, 3/2, 2]
!      h    : E -> reals  [1/2, 1/2, 1/2, 1/2]
!
!                    THE CENTRAL TRUTH
!
!      STATE VALUES LIVE ON Q, INDEPENDENTLY OF ANY GRAPH.
!
! No graph is constructed at this level. None is needed: production
! field calculus already says a field is a function over one set
! identity, and Q is one. That capability is not
! discovered here - it is EXERCISED here, by a client whose state
! domain is emphatically not anybody's vertex set.
!
! Constructing a two-vertex graph so that q0 could live on its
! vertices would be manufacturing the conflation this tower exists
! to refuse, one rung below the rung that tests it. So the level
! imports the field types and nothing that marches, steps or
! solves.
!
!                    VALUES ARE NOT STRUCTURE
!
! The instant t2 is a MEMBER of T. The real 1.0 is the VALUE of a
! field at t2. Four objects stay apart here where a looser client
! would keep one: the carrier T, the structure Tail/Head/A1/A2 over
! it, the coordinates time : T -> reals, and the sizes h : E -> reals.
!
! And the consistency between values and structure is PROVED rather
! than assumed, using the relations Level 1 earned:
!
!      time(head(e)) - time(tail(e)) = h(e)     for every step e
!
! That is field calculus over established structure. It is not a
! time scheme, it discretizes nothing, and it would hold for any
! monotone coordinate whatever.
!
! Values are read through each domain's LOCAL POSITION, never by
! treating a member's integer as an index.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program time_level_5

  use iso_fortran_env       , only : dp => REAL64
  use time_assert           , only : report, verdict
  use time_assert           , only : NQ, NT, NE, TOL
  use time_assert           , only : C_X, C_Y, T0, T2, T4, E1
  use time_assert           , only : H_STEP, TIME_COORD, Q0
  use fractal_graph        , only : set_graph => graph
  use graph_set_map        , only : set_map
  use graph_binary_relation , only : csr_relation
  use class_graph_field     , only : field
  use time_carriers_fixture , only : time_carriers
  use time_relations_fixture, only : tail_relation, head_relation
  use time_fields_fixture   , only : state_field, instant_coordinates, &
       &                             step_sizes

  implicit none

  type(set_graph)  :: q, t, e
  type(set_map)  :: sets
  type(csr_relation) :: tail, head
  type(field)        :: qf, tf, hf
  integer            :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "time integration tower . level 5 . fields"
  write(*,'(1x,a)') "============================================="

  call time_carriers(sets, q, t, e)
  tail = tail_relation(e, t, sets)
  head = head_relation(e, t, sets)

  qf = state_field(q)
  tf = instant_coordinates(t)
  hf = step_sizes(e)

  call check_domains_by_identity(nfail)
  call check_domains_are_distinct(nfail)
  call check_values_by_local_position(nfail)
  call check_coordinates_agree_with_structure(nfail)

  call verdict(nfail, "level 5")

contains

  !===================================================================!
  ! THE level's central truth: each field answers a DECLARED
  ! carrier, and the state field answers Q without a graph existing
  ! anywhere in this program.
  !===================================================================!

  subroutine check_domains_by_identity(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: d

    d = qf % domain()
    call report(d % same_as(q), &
         & "domain(q0) IS Q - the state field lives on the state " // &
         & "carrier, and no graph was built in this program", nfail)

    d = tf % domain()
    call report(d % same_as(t), &
         & "domain(time) IS T", nfail)

    d = hf % domain()
    call report(d % same_as(e), &
         & "domain(h) IS E - one step size per step", nfail)

    call report(qf % num_entries() .eq. NQ .and. &
         &      tf % num_entries() .eq. NT .and. &
         &      hf % num_entries() .eq. NE, &
         & "and each holds exactly its domain's worth of values: " // &
         & "2, 5 and 4", nfail)

  end subroutine check_domains_by_identity

  !===================================================================!
  ! Three fields, three domains, no two the same. A field carrying
  ! the right number of reals on the wrong carrier would be a
  ! different mathematical object, and only identity says so.
  !===================================================================!

  subroutine check_domains_are_distinct(nfail)

    integer, intent(inout) :: nfail

    type(set_graph) :: dq, dt, de

    dq = qf % domain()
    dt = tf % domain()
    de = hf % domain()

    call report(.not. dq % same_as(dt), &
         & "domain(q0) is NOT domain(time): the state axis is not " // &
         & "the time axis", nfail)
    call report(.not. dq % same_as(de), &
         & "domain(q0) is NOT domain(h): a state coordinate is not " // &
         & "a step", nfail)
    call report(.not. dt % same_as(de), &
         & "and domain(time) is NOT domain(h) either", nfail)

  end subroutine check_domains_are_distinct

  !===================================================================!
  ! Values are read through the DOMAIN'S local position. The
  ! carriers here happen to enumerate 1..n, so member value and
  ! local position coincide - which is exactly when a raw-integer
  ! read would pass while meaning nothing.
  !===================================================================!

  subroutine check_values_by_local_position(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: v(:)
    logical               :: ok
    integer               :: i

    ! Read against the ORACLES in time_assert, never against a
    ! literal retyped here - a test that spells its own expected
    ! value twice has only proved it can copy.
    call qf % get_real_vector(v)
    call report(abs(v(sets % index_in(q, C_X)) - Q0(1)) .lt. TOL .and. &
         &      abs(v(sets % index_in(q, C_Y)) - Q0(2)) .lt. TOL, &
         & "q0(x) = 2 and q0(y) = 0, read at Q's local positions", &
         & nfail)

    call tf % get_real_vector(v)
    call report(abs(v(sets % index_in(t, T0)) - TIME_COORD(1)) .lt. TOL .and. &
         &      abs(v(sets % index_in(t, T2)) - TIME_COORD(3)) .lt. TOL .and. &
         &      abs(v(sets % index_in(t, T4)) - TIME_COORD(5)) .lt. TOL, &
         & "time(t0) = 0, time(t2) = 1, time(t4) = 2", nfail)

    call hf % get_real_vector(v)
    ok = .true.
    do i = 1, sets % size_of(e)
       ok = ok .and. (abs(v(sets % index_in(e, sets % member_of(e, i))) - H_STEP) &
            & .lt. TOL)
    end do
    call report(ok .and. abs(v(sets % index_in(e, E1)) - 0.5_dp) .lt. TOL, &
         & "h(e) = 1/2 at every step - and the uniformity is a " // &
         & "property of the VALUES, not of the type holding them", &
         & nfail)

  end subroutine check_values_by_local_position

  !===================================================================!
  ! THE rung's theorem, and the reason Level 1's relations are still
  ! in scope: the coordinates agree with the structure.
  !
  !      time(head(e)) - time(tail(e)) = h(e)
  !
  ! Neither field knows about the other, and neither knows about
  ! Tail or Head. The agreement is measured, through the relations,
  ! by asking each step which instants it joins.
  !===================================================================!

  subroutine check_coordinates_agree_with_structure(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: tv(:), hv(:)
    integer               :: i, m, from, into
    logical               :: ok

    call tf % get_real_vector(tv)
    call hf % get_real_vector(hv)

    ok = .true.
    do i = 1, sets % size_of(e)
       m    = sets % member_of(e, i)
       from = instant_of(tail, m)
       into = instant_of(head, m)
       ok = ok .and. (from .ne. 0) .and. (into .ne. 0)
       if (from .ne. 0 .and. into .ne. 0) then
          ok = ok .and. &
               & (abs((tv(sets % index_in(t, into)) - &
               &       tv(sets % index_in(t, from))) - &
               &      hv(sets % index_in(e, m))) .lt. TOL)
       end if
    end do

    call report(ok, &
         & "time(head(e)) - time(tail(e)) = h(e) for EVERY step: " // &
         & "the coordinates agree with the structure, and neither " // &
         & "was told about the other", nfail)

    ! The same statement, refused where it should be: the first
    ! instant's coordinate is not a step size, and no arithmetic on
    ! members produced any of this.
    call report(abs(tv(sets % index_in(t, T4)) - tv(sets % index_in(t, T0)) &
         &          - real(NE, dp) * H_STEP) .lt. TOL, &
         & "and the whole span t4 - t0 is four steps of 1/2 - two, " // &
         & "as the specimen says", nfail)

  end subroutine check_coordinates_agree_with_structure

  !===================================================================!
  ! Which instant a step meets under a given relation, found by
  ! asking rather than by arithmetic on the step's number.
  !===================================================================!

  integer function instant_of(r, step)

    type(csr_relation), intent(in) :: r
    integer           , intent(in) :: step

    integer :: j, m

    instant_of = 0
    do j = 1, sets % size_of(t)
       m = sets % member_of(t, j)
       if (r % has([step, m])) instant_of = m
    end do

  end function instant_of

end program time_level_5
