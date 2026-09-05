!=====================================================================!
! THE VALUED INTERPRETER - earned at LEVEL 5.
!
! It draws the same grid the Level-4 renderer draws, and writes a
! coefficient where Level 4 writes a hash:
!
!      #   structurally present
!      0   structurally present, carrying the value zero
!      .   structurally ABSENT
!
! so that within this representation
!
!      0  is not  .
!
! and that inequality is executable rather than prose.
!
!                    TWO QUESTIONS, ASKED SEPARATELY
!
! Every cell is the answer to two independent questions, in this order
! and never the other:
!
!      IS IT THERE      glyph_at(D, x, y), which is relation % has -
!                       the Level-4 renderer's own structural
!                       primitive, called here rather than
!                       reimplemented
!
!      WHAT IS IT       the unique occurrence e with T(e)=x and
!                       H(e)=y, and then w(e)
!
! THE SECOND QUESTION IS NEVER ALLOWED TO ANSWER THE FIRST. Nothing
! here reads a coefficient to decide whether a cell is filled; there
! is no
!
!      if (coefficient /= 0) draw('#')
!
! anywhere, and there must never be. A cell is absent because the
! relation does not hold the tuple, and for no other reason. If value
! ever decided presence, a zero coefficient would silently delete a
! dependency - which is precisely the confusion this level exists to
! forbid.
!
!                       THE DEPENDENCY DIRECTION
!
!      L4 structural renderer
!              ^
!      L5 valued renderer
!
! and not the other way. Level 4 knows nothing about fields and must
! stay able to render mathematics that has no numbers at all; this
! module imports Level 4's renderer, and Level 4 imports nothing of
! this.
!
! Nothing in the Level-4 renderer was modified to make this work. The
! valued picture is drawn on the same page as the structural one and
! takes its width and its stub POSITION from a structural picture of
! the same relation - so the two line up row for row when printed
! together, and no measurement is duplicated as a literal.
!
!                    WHAT IS NOT RENDERED, AND WHY
!
! Only D1, D2 and D3 get coefficients. D2:1, D3:1 and every transpose
! stay structural, because a coefficient for a COMPOSED dependency
! would require choosing an algebra for numerical composition - sums
! over intermediate members, which is operator mathematics and not
! field calculus. Level 5 asks only whether values can inhabit
! already-defined structural occurrences.
!
! There is no w^T either. Transposing a coefficient is a numerical
! act, and the structural reverse Gate A established needs no help
! from one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module valued_renderer_fixture

  use iso_fortran_env               , only : dp => REAL64
  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use map_label      , only : label_map
  use relation_finitary                , only : relation
  use field_calculus          , only : FIELD_REAL
  use field_stored             , only : stored_field
  use visualization_carriers_fixture, only : label_for
  use structural_renderer_fixture   , only : picture, sparsity_picture
  use structural_renderer_fixture   , only : glyph_at

  implicit none

  private
  public :: valued_sparsity_picture
  public :: coefficients_fit, occurrence_joining, occurrences_joining
  public :: value_at, value_token
  public :: ABSENT

  !-------------------------------------------------------------------!
  ! What an absent cell is written as. It is a glyph and not a number,
  ! which is the whole point: there is no value there to print,
  ! because there is no dependency there to carry one.
  !-------------------------------------------------------------------!

  character(len=*), parameter :: ABSENT = '.'

contains

  !===================================================================!
  ! Do these coefficients belong to these occurrences.
  !
  ! Asked by IDENTITY - the field's declared domain must BE the
  ! occurrence carrier, not merely count the same. |X0| = |E2| = 4 in
  ! this specimen, and a four-valued field on X0 is not a coefficient
  ! field for E2.
  !
  ! This invents no validation mechanism. It asks the field for its
  ! domain and the domain for its identity, which is what the nucleus
  ! has always offered; the level above exposes those semantics rather
  ! than adding to them.
  !===================================================================!

  logical function coefficients_fit(w, occurrences)

    class(stored_field)     , intent(in) :: w
    type(graph), intent(in) :: occurrences

    type(graph) :: on

    on = w % domain()

    coefficients_fit = on % same_as(occurrences)
    coefficients_fit = coefficients_fit .and. (w % num_components() .eq. 1)
    coefficients_fit = coefficients_fit .and. &
         &             (w % value_kind() .eq. FIELD_REAL)

  end function coefficients_fit

  !===================================================================!
  ! The occurrence that joins x to y, or zero if none does.
  !
  ! Found by walking the occurrence carrier in ITS declaration order
  ! and asking the two primitive relations - never by arithmetic on
  ! raw member integers, and never by assuming an occurrence's number
  ! is a position in a grid.
  !===================================================================!

  integer function occurrence_joining(tail, head, occurrences, from, to, sets)

    class(relation)  , intent(in) :: tail, head
    type(set_map)  , intent(in) :: sets
    type(graph), intent(in) :: occurrences
    integer          , intent(in) :: from, to

    integer :: k, e

    occurrence_joining = 0
    do k = 1, sets % num_members_of(occurrences)
       e = sets % member_of(occurrences, k)
       if (tail % has([e, from]) .and. head % has([e, to])) then
          occurrence_joining = e
          return
       end if
    end do

  end function occurrence_joining

  !===================================================================!
  ! How many occurrences join x to y. One, for every tuple of a direct
  ! dependency in this specimen - and the level above proves it,
  ! because a coefficient picture is only well defined where the seat
  ! is unique.
  !===================================================================!

  integer function occurrences_joining(tail, head, occurrences, from, to, sets)

    class(relation)  , intent(in) :: tail, head
    type(graph), intent(in) :: occurrences
    integer          , intent(in) :: from, to
    type(set_map)  , intent(in) :: sets

    integer :: k, e

    occurrences_joining = 0
    do k = 1, sets % num_members_of(occurrences)
       e = sets % member_of(occurrences, k)
       if (tail % has([e, from]) .and. head % has([e, to])) then
          occurrences_joining = occurrences_joining + 1
       end if
    end do

  end function occurrences_joining

  !===================================================================!
  ! The value at one occurrence. The convenience form: it fetches the
  ! whole vector and indexes it by the domain's own local_index, which
  ! is where the field says that member's value sits. The renderer
  ! below fetches once and indexes many times, as the field contract
  ! intends.
  !===================================================================!

  real(dp) function value_at(w, occurrences, member, sets)

    class(stored_field)     , intent(in) :: w
    type(graph), intent(in) :: occurrences
    integer          , intent(in) :: member
    type(set_map)  , intent(in) :: sets

    real(dp), allocatable :: values(:)
    integer               :: seat

    call w % real_vector(values)

    seat = sets % index_in(occurrences, member)
    if (seat .lt. 1 .or. seat .gt. size(values)) then
       error stop 'valued_renderer_fixture: that member has no seat in this field'
    end if

    value_at = values(seat)

  end function value_at

  !===================================================================!
  ! How a value is written. DISPLAY ONLY: a whole number loses its
  ! empty fractional part so the picture reads like arithmetic rather
  ! than like output.
  !
  ! Note carefully that this is the one place a coefficient's VALUE is
  ! inspected, and it decides only how the number looks. It never
  ! decides whether a cell is filled - a zero prints as 0 and stays in
  ! the picture.
  !===================================================================!

  function value_token(v) result(text)

    real(dp), intent(in) :: v

    character(len=:), allocatable :: text
    character(len=32)             :: buf

    if (abs(v) .lt. 1.0e9_dp .and. v .eq. real(nint(v), dp)) then
       write(buf, '(i0)') nint(v)
    else
       write(buf, '(f0.3)') v
    end if

    text = trim(buf)

  end function value_token

  !===================================================================!
  ! THE COEFFICIENT PICTURE.
  !
  !      rows    = codomain, in ITS declaration order
  !      columns = domain,   in ITS declaration order
  !      cell    = w(e) where the dependency is present
  !                '.' where it is not
  !
  ! Drawn on the page a structural picture of the same relation
  ! occupies, with the same stub, so the two align row for row.
  !===================================================================!

  type(picture) function valued_sparsity_picture(d, tail, head, &
       & occurrences, w, sets, labels) result(pic)

    class(relation)  , intent(in) :: d
    class(relation)  , intent(in) :: tail, head
    type(graph), intent(in) :: occurrences
    class(stored_field)     , intent(in) :: w
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    type(picture)                  :: page
    type(graph) :: cols, rows
    character(len=:) , allocatable :: cell
    real(dp)         , allocatable :: values(:)
    integer :: stub, wide, i, j, at, e, seat, width

    if (d % arity() .ne. 2) then
       error stop 'valued_renderer_fixture: a coefficient picture reads a binary relation'
    end if
    if (.not. coefficients_fit(w, occurrences)) then
       error stop 'valued_renderer_fixture: those coefficients do not live on those occurrences'
    end if

    cols = d % domain(1)
    rows = d % domain(2)

    ! The page, and where its first column stands - taken from a
    ! STRUCTURAL picture of this same relation rather than from a
    ! measurement copied out of Level 4.
    page  = sparsity_picture(d, sets, labels)
    width = len(page % line)
    stub  = first_nonblank(page % line(2))

    call w % real_vector(values)

    wide = max(widest_label(cols, sets, labels), widest_value(values)) + 3

    allocate(character(len=width) :: pic % line(2 + sets % num_members_of(rows)))
    pic % line = repeat(' ', width)

    call put(pic % line(1), 1, d % name() // ' VALUES')

    do j = 1, sets % num_members_of(cols)
       at = stub + (j - 1) * wide
       call put(pic % line(2), at, &
            &   right(label_for(cols, sets % member_of(cols, j), labels), wide - 1))
    end do

    do i = 1, sets % num_members_of(rows)
       call put(pic % line(2 + i), 1, label_for(rows, sets % member_of(rows, i), labels))
       do j = 1, sets % num_members_of(cols)
          at = stub + (j - 1) * wide

          ! FIRST the structural question, asked of the relation.
          if (glyph_at(d, sets % member_of(cols, j), sets % member_of(rows, i)) .eq. '#') then

             ! ONLY THEN the value question.
             e = occurrence_joining(tail, head, occurrences, &
                  &                 sets % member_of(cols, j), sets % member_of(rows, i), sets)
             if (e .eq. 0) then
                error stop 'valued_renderer_fixture: a present dependency with no occurrence to seat it'
             end if
             seat = sets % index_in(occurrences, e)
             cell = value_token(values(seat))
          else
             cell = ABSENT
          end if

          call put(pic % line(2 + i), at, right(cell, wide - 1))
       end do
    end do

  end function valued_sparsity_picture

  !===================================================================!
  ! Small mechanics.
  !===================================================================!

  integer function first_nonblank(line)

    character(len=*), intent(in) :: line

    integer :: k

    first_nonblank = 1
    do k = 1, len(line)
       if (line(k:k) .ne. ' ') then
          first_nonblank = k
          return
       end if
    end do

  end function first_nonblank

  function right(text, wide) result(padded)

    character(len=*), intent(in) :: text
    integer         , intent(in) :: wide

    character(len=:), allocatable :: padded

    if (len(text) .ge. wide) then
       padded = text
    else
       padded = repeat(' ', wide - len(text)) // text
    end if

  end function right

  integer function widest_label(carrier, sets, labels)

    type(graph), intent(in) :: carrier
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    integer :: k

    widest_label = 1
    do k = 1, sets % num_members_of(carrier)
       widest_label = max(widest_label, &
            &             len(label_for(carrier, sets % member_of(carrier, k), labels)))
    end do

  end function widest_label

  integer function widest_value(values)

    real(dp), intent(in) :: values(:)

    integer :: k

    widest_value = len(ABSENT)
    do k = 1, size(values)
       widest_value = max(widest_value, len(value_token(values(k))))
    end do

  end function widest_value

  subroutine put(line, at, text)

    character(len=*), intent(inout) :: line
    integer         , intent(in)    :: at
    character(len=*), intent(in)    :: text

    if (at + len(text) - 1 .gt. len(line)) then
       error stop 'valued_renderer_fixture: the picture is wider than its page'
    end if

    line(at : at + len(text) - 1) = text

  end subroutine put

end module valued_renderer_fixture
