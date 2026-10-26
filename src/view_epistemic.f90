!=====================================================================!
! Module that states the epistemic view: one view of a graph, in
! which its two branches are read as
!
!     branch(1) = Q     the data
!     branch(2) = R     the residual, the operator governing it
!
! Q and R are graphs. The view stores nothing, declares no type and
! adds no state; it classifies the branches' statuses and returns
! references.
!
! The domain of the classification is the four combinations in which
! both branches are UNKNOWN or KNOWN:
!
!                      R UNKNOWN       R KNOWN
!     Q UNKNOWN        void            operator
!     Q KNOWN          data            realized
!
! NULL is outside that domain. The view has no name for a branch that
! is definitely absent, and the 3x3 state space is not forced back
! into the old 2x2:
!
!                  R:  NULL     UNKNOWN     KNOWN
!     Q:  NULL         -        -           -
!         UNKNOWN      -        void        operator
!         KNOWN        -        data        realized
!
!     -  outside the domain: epistemic_defined is false, and
!        epistemic_name refuses
!
! Realized is not solved. A realized graph asserts that both branches
! reference a graph, and nothing more; R(Q) = 0 is a separate property
! with its own predicate, asserted elsewhere.
!
! Unknown is not null. Unrealized is not absent, and neither is an
! empty Q: a Q that references a graph with no members is KNOWN.
!
! No host. A graph is not attached to a graph of another kind; where a
! view needs another graph it references it, or an external map
! associates the two.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_epistemic

  use graph_fractal, only : graph, BRANCH_UNKNOWN, BRANCH_KNOWN

  implicit none

  private

  public :: has_data, has_operator
  public :: epistemic_defined, epistemic_name
  public :: data_of, residual_of

contains

  !===================================================================!
  ! Returns whether the data branch, Q, is KNOWN: the first of the
  ! view's two primitive predicates
  !===================================================================!

  logical function has_data(g)

    type(graph), intent(in) :: g

    has_data = g % branch(1) % status() .eq. BRANCH_KNOWN

  end function has_data

  !===================================================================!
  ! Returns whether the operator branch, R, is KNOWN: the second of the
  ! view's two primitive predicates
  !===================================================================!

  logical function has_operator(g)

    type(graph), intent(in) :: g

    has_operator = g % branch(2) % status() .eq. BRANCH_KNOWN

  end function has_operator

  !===================================================================!
  ! Returns whether the graph is in the domain of the view: whether
  ! both branches are UNKNOWN or KNOWN
  !===================================================================!

  logical function epistemic_defined(g)

    type(graph), intent(in) :: g

    integer :: s

    epistemic_defined = .true.
    do s = 1, 2
       epistemic_defined = epistemic_defined .and. &
            & (g % branch(s) % status() .eq. BRANCH_UNKNOWN .or. &
            &  g % branch(s) % status() .eq. BRANCH_KNOWN)
    end do

  end function epistemic_defined

  !===================================================================!
  ! Returns the canonical name of the graph, for the four combinations
  ! that have one: void, data, operator or realized
  !===================================================================!

  function epistemic_name(g) result(name)

    type(graph)     , intent(in)  :: g
    character(len=:), allocatable :: name

    if (.not. epistemic_defined(g)) then
       error stop 'view_epistemic: epistemic_name requires an epistemically-defined graph, &
            &but epistemic_defined(g) is false'
    end if

    if (has_data(g)) then
       if (has_operator(g)) then
          name = 'realized'
       else
          name = 'data'
       end if
    else
       if (has_operator(g)) then
          name = 'operator'
       else
          name = 'void'
       end if
    end if

  end function epistemic_name

  !===================================================================!
  ! Returns a reference to the data, Q. Refused unless the branch is
  ! KNOWN: neither NULL nor UNKNOWN is a value, and no accessor
  ! manufactures one.
  !===================================================================!

  function data_of(g) result(q)

    type(graph), intent(in) :: g
    type(graph), pointer    :: q

    if (.not. has_data(g)) then
       error stop 'view_epistemic: data_of requires the data branch to be KNOWN, &
            &but has_data(g) is false'
    end if

    q => g % branch(1) % known()

  end function data_of

  !===================================================================!
  ! Returns a reference to the residual, R. Refused unless the branch
  ! is KNOWN.
  !===================================================================!

  function residual_of(g) result(r)

    type(graph), intent(in) :: g
    type(graph), pointer    :: r

    if (.not. has_operator(g)) then
       error stop 'view_epistemic: residual_of requires the operator branch to be KNOWN, &
            &but has_operator(g) is false'
    end if

    r => g % branch(2) % known()

  end function residual_of

end module view_epistemic
