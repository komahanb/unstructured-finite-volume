!=====================================================================!
! THE LABEL PROVENANCE SUITE
!
! map_label is the third identity-keyed association, and it signs
! the same storage law as the other two:
!
!     an identity map owns its keys by value;
!     it borrows no graph object merely to recognize it.
!
! What is proved here beyond that law is the one thing a label must
! never become:
!
!     two graphs may carry the SAME label and remain two sets
!     a copied token retrieves the SAME label
!
! Together those say a label names and does not address. The first
! forbids using a label to decide identity; the second says identity,
! once decided, determines the label - so the map is a function of the
! token and of nothing else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program labels

  use graph_fractal  , only : graph
  use map_label, only : label_map

  implicit none

  integer :: failures = 0

  write(*,'(1x,a)') "label provenance suite"

  !===================================================================!
  ! 1 . A LABEL IS NOT AN IDENTITY.
  !
  ! Two independently declared graphs, one string. They answer the
  ! same name and are not the same set - exactly as two sets with
  ! equal extensions are two sets.
  !===================================================================!

  naming_block: block

    type(graph)     :: a, b
    type(label_map) :: names

    call a % declare(); call b % declare()

    call names % bind(a, 'cells')
    call names % bind(b, 'cells')

    call check('1  two graphs may carry one label', &
         & names % label_of(a) .eq. 'cells' .and. &
         &  names % label_of(b) .eq. 'cells')

    call check('1  and remain two sets', .not. a % same_as(b))

    call check('1  both are labelled', &
         & names % labelled(a) .and. names % labelled(b))

  end block naming_block

  !===================================================================!
  ! 2 . A COPIED TOKEN RETRIEVES THE SAME LABEL.
  !
  ! A copy of a graph IS the same declared set, so it must answer the
  ! same name. This is the half that makes the map a function of the
  ! identity rather than of the variable.
  !===================================================================!

  copy_block: block

    type(graph)     :: a, twin
    type(label_map) :: names

    call a % declare()
    call names % bind(a, 'interior_vertices')

    twin = a

    call check('2  a copy is the same set', twin % same_as(a))
    call check('2  and answers the same label', &
         & names % label_of(twin) .eq. 'interior_vertices')

  end block copy_block

  !===================================================================!
  ! 3 . THE UNNAMED ANSWER. Being unnamed is not a contradiction:
  ! metadata is allowed to be absent, and '' is the carriers' own
  ! convention for it.
  !===================================================================!

  unnamed_block: block

    type(graph)     :: a, stranger
    type(label_map) :: names

    call a % declare(); call stranger % declare()
    call names % bind(a, 'faces')

    call check('3  an unbound graph is not labelled', &
         & .not. names % labelled(stranger))

    call check('3  and answers the empty string, not a refusal', &
         & names % label_of(stranger) .eq. '')

    call check('3  while an empty map answers for nobody', &
         & unlabelled_in_empty_map(a))

  end block unnamed_block

  !===================================================================!
  ! 4 . THE STORAGE LAW, over the label map. Bind inside a scope whose
  ! graphs are heap-allocated and destroyed on return, then ask about
  ! copies carrying the same tokens. No TARGET is declared anywhere
  ! below, which is the compile-time half of the proof.
  !===================================================================!

  lifetime_block: block

    type(label_map) :: names
    type(graph)     :: a, b

    call name_and_die(names, a, b)

    call check('4  the map still names a set whose binder is gone', &
         & names % label_of(a) .eq. 'owned_vertices')

    call check('4  and the second, distinctly', &
         & names % label_of(b) .eq. 'borrowed_vertices')

    call check('4  and the two remain two sets', .not. a % same_as(b))

  end block lifetime_block

  !===================================================================!
  ! 5 . A MAP IS A VALUE: intrinsic assignment deep-copies the rows,
  ! keys and labels included.
  !===================================================================!

  value_block: block

    type(label_map) :: names, twin
    type(graph)     :: a, extra

    call a % declare(); call extra % declare()
    call names % bind(a, 'cells')

    twin = names
    call twin % bind(extra, 'faces')

    call check('5  a copied map answers for the keys it copied', &
         & twin % label_of(a) .eq. 'cells')

    call check('5  and growth of the copy does not reach the original', &
         & twin % labelled(extra) .and. .not. names % labelled(extra))

  end block value_block

  !===================================================================!

  if (failures .eq. 0) then
     print *, ''
     print *, ' ALL PROPOSITIONS HOLD'
  else
     print *, ''
     print *, ' FAILURES :', failures
     error stop 'labels: a proposition failed'
  end if

contains

  !===================================================================!
  ! Name two sets and destroy every graph object the map was built
  ! from. No TARGET: after the gate, binding does not ask for one.
  !===================================================================!

  subroutine name_and_die(m, a_key, b_key)

    type(label_map), intent(out) :: m
    type(graph)    , intent(out) :: a_key, b_key

    type(graph), allocatable :: a, b

    allocate(a, b)
    call a % declare(); call b % declare()

    call m % bind(a, 'owned_vertices')
    call m % bind(b, 'borrowed_vertices')

    a_key = a
    b_key = b

    deallocate(a, b)

  end subroutine name_and_die

  logical function unlabelled_in_empty_map(g) result(ok)

    type(graph), intent(in) :: g

    type(label_map) :: empty

    ok = .not. empty % labelled(g) .and. empty % label_of(g) .eq. ''

  end function unlabelled_in_empty_map

  subroutine check(label, ok)

    character(len=*), intent(in) :: label
    logical         , intent(in) :: ok

    if (ok) then
       print *, ' PASS : ', label
    else
       print *, ' FAIL : ', label
       failures = failures + 1
    end if

  end subroutine check

end program labels
