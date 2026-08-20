!=====================================================================!
! EPISTEMIC VIEW
!
! One view of a graph, in which the two branches are interpreted as
!
!     branch(1) = Q     the data
!     branch(2) = R     the residual/operator governing it
!
! Q and R are graphs. The view stores nothing, declares no type and
! adds no state; it classifies branch statuses and returns references.
!
! The domain of the classification is the four combinations in which
! both branches are UNKNOWN or KNOWN:
!
!     (UNKNOWN, UNKNOWN)   void
!     (KNOWN  , UNKNOWN)   data
!     (UNKNOWN, KNOWN  )   operator
!     (KNOWN  , KNOWN  )   realized
!
! NULL is outside that domain. The view has no name for a branch
! that is definitely absent, and the 3x3 state space is not forced
! back into the old 2x2: epistemic_defined answers whether a name
! exists, and epistemic_name refuses when one does not.
!
! REALIZED IS NOT SOLVED. A realized graph asserts that both branches
! reference a graph, and nothing more. R(Q) = 0 is a separate property
! with its own words, asserted elsewhere.
!
! UNKNOWN IS NOT NULL. Unrealized is not absent, and neither is an
! empty Q: a Q that references a graph with no members is KNOWN.
!
! No host. A graph does not ride on a graph of another kind; where a
! view needs another graph it references it, or an external map
! associates the two.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_epistemic

  use graph_fractal, only : graph, GRAPH_UNKNOWN, GRAPH_KNOWN

  implicit none

  private
  public :: has_data, has_operator
  public :: epistemic_defined, epistemic_name
  public :: data_of, residual_of

contains

  !===================================================================!
  ! The two primitive questions of this view.
  !===================================================================!

  logical function has_data(g)

    type(graph), intent(in) :: g

    has_data = g % branch(1) % status() .eq. GRAPH_KNOWN

  end function has_data

  logical function has_operator(g)

    type(graph), intent(in) :: g

    has_operator = g % branch(2) % status() .eq. GRAPH_KNOWN

  end function has_operator

  !===================================================================!
  ! Is this graph in the domain of the view: are both branches
  ! UNKNOWN or KNOWN.
  !===================================================================!

  logical function epistemic_defined(g)

    type(graph), intent(in) :: g

    integer :: s

    epistemic_defined = .true.
    do s = 1, 2
       epistemic_defined = epistemic_defined .and. &
            & (g % branch(s) % status() .eq. GRAPH_UNKNOWN .or. &
            &  g % branch(s) % status() .eq. GRAPH_KNOWN)
    end do

  end function epistemic_defined

  !===================================================================!
  ! The canonical name, for the four combinations that have one.
  !===================================================================!

  function epistemic_name(g) result(said)

    type(graph), intent(in)       :: g
    character(len=:), allocatable :: said

    if (.not. epistemic_defined(g)) then
       error stop 'view_epistemic: NULL has no epistemic name'
    end if

    if (has_data(g)) then
       if (has_operator(g)) then
          said = 'realized'
       else
          said = 'data'
       end if
    else
       if (has_operator(g)) then
          said = 'operator'
       else
          said = 'void'
       end if
    end if

  end function epistemic_name

  !===================================================================!
  ! The two references. Refused unless the branch is KNOWN: neither
  ! NULL nor UNKNOWN is a value, and no accessor manufactures one.
  !===================================================================!

  function data_of(g) result(q)

    type(graph), intent(in) :: g
    type(graph), pointer    :: q

    if (.not. has_data(g)) then
       error stop 'view_epistemic: Q is not KNOWN'
    end if

    q => g % branch(1) % known()

  end function data_of

  function residual_of(g) result(r)

    type(graph), intent(in) :: g
    type(graph), pointer    :: r

    if (.not. has_operator(g)) then
       error stop 'view_epistemic: R is not KNOWN'
    end if

    r => g % branch(2) % known()

  end function residual_of

end module view_epistemic
