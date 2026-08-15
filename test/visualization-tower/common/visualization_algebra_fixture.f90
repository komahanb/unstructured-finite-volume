!=====================================================================!
! THE SPECIMEN'S DERIVATIONS - earned at LEVEL 2.
!
! Four small operations, each one a naming of something the nucleus
! already does. Nothing here computes a dependency: the relation
! algebra computes it, and this file only says which composition to
! take, in which order, under which name.
!
!                    MATHEMATICS AND SPELLING, SIDE BY SIDE
!
! The repository's composition reads its arguments in the order the
! data flows, and its result in the order mathematics writes it:
!
!      compose_binary(P_AB, P_BC)  =  P_BC o P_AB
!
! so the FIRST argument is applied FIRST. Written out for the three
! derivations this tower needs:
!
!   mathematics                          Fortran
!   -----------------------------------  ---------------------------
!   D_k    = H_k o T_k^T                 compose_binary(T_k^T, H_k)
!   D_2:1  = D_2 o D_1                   compose_binary(D_1, D_2)
!   D_rev  = D_1^T o D_2^T o D_3^T       compose_binary(
!                                          compose_binary(D_3^T, D_2^T),
!                                          D_1^T)
!
! The two directions of reading are genuinely opposite, and this
! comment is the only place in the tower where that is allowed to be
! resolved. Every level above says compose and means the arrow.
!
!                     WHY D_k = H_k o T_k^T
!
! T_k^T sends a member of X_(k-1) to the occurrences that READ it;
! H_k sends an occurrence to the member of X_k it WRITES. Composed,
! they send x to everything x can reach in one operator - which is
! the dependency, derived, with the occurrence carrier forgotten on
! the way through.
!
!                       WHAT `named` IS FOR, AND NOT FOR
!
! compose_binary and transpose_of both name their results after their
! arguments, which after three compositions reads like a receipt
! rather than a relation. `named` re-materializes a binary relation
! over its OWN two domains, from its OWN tuples, under a chosen
! label - so D1 is called D1 in a picture instead of
! 'T1^T then H1'.
!
! It changes the name and nothing else, and Level 2 proves that by
! holding every named copy against the relation it was made from,
! extensionally. A relabelling that quietly altered an extension
! would be caught there, not trusted here.
!
!                        EXTENSION, NOT STAMP
!
! same_extension is the tower's equality. Two relations are equal
! when their domains are the same DECLARED domains, slot for slot,
! and they hold exactly the same tuples - each in the other, both
! ways round. Never by identity: a transpose view and its
! materialized copy are two objects, and the law they satisfy is
! about content. Never by tuple order either: a relation is a set,
! and how its tuples were handed in is no part of what it is.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module visualization_algebra_fixture

  use graph_carrier         , only : member_set
  use graph_relation        , only : relation
  use graph_binary_relation , only : binary_relation, csr_relation
  use graph_binary_relation , only : transposed_view, transpose_of
  use graph_relation_algebra, only : compose_binary

  implicit none

  private
  public :: derive_dependency, derive_composition
  public :: materialized_transpose, named, same_extension

contains

  !===================================================================!
  ! The dependency of one operator, from its two primitive incidence
  ! relations:
  !
  !      D_k  =  H_k o T_k^T  :  X_(k-1) -> X_k
  !
  ! The occurrence carrier E_k is the middle domain of the
  ! composition, and composition is what forgets it. Two occurrences
  ! that read the same member and write the same member therefore
  ! yield ONE tuple - which is the specimen's central law, and it is
  ! the nucleus's constructor that enforces it, not this file.
  !===================================================================!

  type(csr_relation) function derive_dependency(label, tail, head) &
       & result(d)

    character(len=*)              , intent(in) :: label
    class(binary_relation), target, intent(in) :: tail
    class(relation)               , intent(in) :: head

    type(transposed_view) :: reads

    reads = transpose_of(tail)
    d     = named(label, compose_binary(reads, head))

  end function derive_dependency

  !===================================================================!
  ! One dependency after another:
  !
  !      derive_composition(name, first, second)  =  second o first
  !
  ! The argument order is the order of TRAVERSAL, matching
  ! compose_binary; the name says so, so no caller has to remember
  ! which way the circle points.
  !===================================================================!

  type(csr_relation) function derive_composition(label, first, second) &
       & result(d)

    character(len=*), intent(in) :: label
    class(relation) , intent(in) :: first
    class(relation) , intent(in) :: second

    d = named(label, compose_binary(first, second))

  end function derive_composition

  !===================================================================!
  ! The structural reverse of one dependency, as an owned relation.
  !
  ! transpose_of answers a VIEW - O(1), nothing copied, and a
  ! borrower that lives only as long as its base. A view is exactly
  ! right for asking a question and exactly wrong for handing back
  ! out of a function, so this materializes the view's extension into
  ! a relation whole unto itself.
  !
  ! The nucleus still does the transposing. What is added is
  ! ownership, and Level 2 checks that nothing else was.
  !===================================================================!

  type(csr_relation) function materialized_transpose(label, r) result(rt)

    character(len=*)              , intent(in) :: label
    class(binary_relation), target, intent(in) :: r

    type(transposed_view) :: reversed

    reversed = transpose_of(r)
    rt       = named(label, reversed)

  end function materialized_transpose

  !===================================================================!
  ! Re-materialize a binary relation under a chosen name, over its
  ! own domains, from its own tuples.
  !===================================================================!

  type(csr_relation) function named(label, r) result(copy)

    character(len=*), intent(in) :: label
    class(relation) , intent(in) :: r

    class(member_set), allocatable :: first, second
    integer, allocatable           :: table(:,:)

    if (r % arity() .ne. 2) then
       error stop 'visualization_algebra_fixture: only a binary relation is named here'
    end if

    first  = r % domain(1)
    second = r % domain(2)
    call r % tuples(table)

    copy = csr_relation(label, first, second, table)

  end function named

  !===================================================================!
  ! Extensional equality: the same declared domains slot for slot,
  ! the same number of tuples, and every tuple of each held by the
  ! other. Orientation is part of the answer - a relation over
  ! X3 x X0 is not equal to one over X0 x X3, however its tuples
  ! read.
  !===================================================================!

  logical function same_extension(r, s)

    class(relation), intent(in) :: r
    class(relation), intent(in) :: s

    class(member_set), allocatable :: dr, ds
    integer, allocatable           :: tr(:,:), ts(:,:)
    integer                        :: k, j

    same_extension = (r % arity() .eq. s % arity())
    if (.not. same_extension) return

    do k = 1, r % arity()
       dr = r % domain(k)
       ds = s % domain(k)
       same_extension = same_extension .and. dr % same_as(ds)
    end do
    if (.not. same_extension) return

    same_extension = (r % num_tuples() .eq. s % num_tuples())
    if (.not. same_extension) return

    call r % tuples(tr)
    call s % tuples(ts)

    do j = 1, size(tr, 2)
       same_extension = same_extension .and. s % has(tr(:, j))
    end do
    do j = 1, size(ts, 2)
       same_extension = same_extension .and. r % has(ts(:, j))
    end do

  end function same_extension

end module visualization_algebra_fixture
