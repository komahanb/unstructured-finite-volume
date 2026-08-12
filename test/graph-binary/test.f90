!=====================================================================!
! The binary relation suite: the laws of arity two (AGENTS.md,
! level 1, phase 3).
!
! The CSR citizen answers image and preimage as O(degree) slices,
! speaks members in and members out across any carrier concretion -
! the sparse-source check is the inverse-map law at work - and the
! transpose arrives as a borrowing view: O(1) to make, ends
! swapped, its own identity. The involution (R^T)^T = R is tested
! the only way it is promised: by EXTENSION, tuples and domains,
! never by minted identity.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_binary

  use graph_carrier         , only : counted_set, member_set
  use graph_binary_relation , only : csr_relation, transpose_of, &
       &                             transposed_view
  use listed_set_fixture    , only : listed_set

  implicit none

  integer :: nfail

  nfail = 0

  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph binary relation suite (AGENTS phase 3)"
  write(*,'(1x,a)') "============================================="

  call check_csr_contract(nfail)
  call check_sparse_source(nfail)
  call check_transpose_view(nfail)
  call check_involution_is_extensional(nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all binary relation checks passed"
  else
     write(*,'(1x,a,i0,a)') "FAILED: ", nfail, " binary relation check(s)"
     error stop
  end if

contains

  subroutine report(ok, label, nfail)

    logical         , intent(in)    :: ok
    character(len=*), intent(in)    :: label
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! The CSR citizen on cells x faces: fibres both ways, membership
  ! by row scan, set semantics, the named ends.
  !===================================================================!

  subroutine check_csr_contract(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)              :: cells, faces
    type(csr_relation)             :: r
    class(member_set), allocatable :: d
    integer, allocatable           :: idx(:), t(:,:)

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    ! One duplicate handed in, on purpose.
    r = csr_relation('touches', cells, faces, &
         & reshape([1,1,  1,2,  2,2,  3,4,  1,1], [2, 5]))

    call report(r % arity() .eq. 2, &
         & "arity two, answered once for the whole family", nfail)
    call report(r % num_tuples() .eq. 4, &
         & "a tuple handed in twice is in the relation once", nfail)

    call r % image(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "the image is the row, one slice", nfail)
    call r % image(4, idx)
    call report(size(idx) .eq. 0, &
         & "a member relating to nothing has an empty image", nfail)

    call r % preimage(2, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 2]), &
         & "the preimage is the mirrored row, one slice", nfail)
    call r % preimage(5, idx)
    call report(size(idx) .eq. 0, &
         & "and an unrelated face answers empty too", nfail)

    call report(r % has([1, 2]) .and. .not. r % has([2, 1]), &
         & "membership is one row scan, order honoured", nfail)

    call r % tuples(t)
    call report(size(t, 2) .eq. 4, &
         & "the tuple table carries the set, not the handing", nfail)

    d = r % source()
    call report(d % same_as(cells), &
         & "the source is the first slot's domain, by identity", nfail)
    d = r % target()
    call report(d % same_as(faces), &
         & "the target is the second's", nfail)

  end subroutine check_csr_contract

  !===================================================================!
  ! The inverse-map law at work: a sparse listed carrier as source.
  ! Rows live at local indices; questions and answers speak
  ! members; local_index is the only bridge, and outsiders answer
  ! empty.
  !===================================================================!

  subroutine check_sparse_source(nfail)

    integer, intent(inout) :: nfail

    type(listed_set)               :: sensors
    type(counted_set)              :: cells
    type(csr_relation)             :: r
    class(member_set), allocatable :: d
    integer, allocatable           :: idx(:)

    sensors = listed_set('sensors', [10, 20, 30])
    cells   = counted_set('cells', 3)

    r = csr_relation('reads', sensors, cells, &
         & reshape([10,1,  20,1,  20,3], [2, 3]))

    call r % image(20, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [1, 3]), &
         & "a sparse member's image arrives through local_index", nfail)

    call r % image(30, idx)
    call report(size(idx) .eq. 0, &
         & "a member with no tuples answers empty", nfail)

    call r % image(15, idx)
    call report(size(idx) .eq. 0, &
         & "an outsider answers empty: relating to nothing is an answer", &
         & nfail)

    call r % preimage(1, idx)
    call report(size(idx) .eq. 2 .and. all(idx .eq. [10, 20]), &
         & "the preimage speaks members, never rows", nfail)

    call report(r % has([20, 3]) .and. .not. r % has([10, 3]), &
         & "membership crosses the sparse carrier untroubled", nfail)

    d = r % domain(1)
    call report(d % same_as(sensors), &
         & "the sparse slot answers its own declared domain", nfail)

  end subroutine check_sparse_source

  !===================================================================!
  ! The transpose view: O(1) to make, nothing copied, every answer
  ! the base's with the ends swapped - and its own identity, since
  ! a view is a new declared relation over the same extension.
  !===================================================================!

  subroutine check_transpose_view(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)               :: cells, faces
    type(csr_relation), target      :: r
    type(transposed_view)           :: t
    class(member_set), allocatable  :: d
    integer, allocatable            :: fwd(:), rev(:)

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    r = csr_relation('touches', cells, faces, &
         & reshape([1,1,  1,2,  2,2,  3,4], [2, 4]))

    t = transpose_of(r)

    call report(t % num_tuples() .eq. r % num_tuples(), &
         & "the view counts what the base counts", nfail)

    call t % image(2, fwd)
    call r % preimage(2, rev)
    call report(size(fwd) .eq. size(rev) .and. all(fwd .eq. rev), &
         & "the view's image is the base's preimage", nfail)

    call t % preimage(1, fwd)
    call r % image(1, rev)
    call report(size(fwd) .eq. size(rev) .and. all(fwd .eq. rev), &
         & "and its preimage is the base's image", nfail)

    call report(t % has([2, 1]) .and. .not. t % has([1, 2]), &
         & "membership reads backwards through the view", nfail)

    d = t % domain(1)
    call report(d % same_as(faces), &
         & "the view's source is the base's target", nfail)
    d = t % domain(2)
    call report(d % same_as(cells), &
         & "and its target the base's source", nfail)

    call report(.not. t % same_as(r), &
         & "a view signs its own token: identity is not extension", nfail)
    call report(t % name() == 'touches^T', &
         & "and wears its derivation in its name", nfail)

  end subroutine check_transpose_view

  !===================================================================!
  ! The involution (R^T)^T = R, tested the only way it is promised:
  ! by extension. The double view carries R's tuples over R's
  ! domains - and still is not R by identity, because same_as
  ! answers stamps and no canonicalization is promised.
  !===================================================================!

  subroutine check_involution_is_extensional(nfail)

    integer, intent(inout) :: nfail

    type(counted_set)               :: cells, faces
    type(csr_relation) , target    :: r
    type(transposed_view), target  :: t
    type(transposed_view)          :: tt
    class(member_set), allocatable :: da, db
    integer, allocatable           :: rt(:,:), ttt(:,:)
    integer                        :: j
    logical                        :: ok

    cells = counted_set('cells', 4)
    faces = counted_set('faces', 5)

    r = csr_relation('touches', cells, faces, &
         & reshape([1,1,  1,2,  2,2,  3,4], [2, 4]))

    t  = transpose_of(r)
    tt = transpose_of(t)

    call r % tuples(rt)
    call tt % tuples(ttt)

    ok = size(ttt, 2) .eq. size(rt, 2)
    if (ok) then
       do j = 1, size(rt, 2)
          ok = ok .and. all(ttt(:, j) .eq. rt(:, j))
       end do
    end if
    call report(ok, &
         & "the double transpose carries the base's tuples", nfail)

    da = tt % domain(1)
    db = r % domain(1)
    call report(da % same_as(db), &
         & "over the base's own domains, slot for slot", nfail)

    call report(.not. tt % same_as(r), &
         & "extension returns; identity, unpromised, does not", nfail)

  end subroutine check_involution_is_extensional

end program test_graph_binary
