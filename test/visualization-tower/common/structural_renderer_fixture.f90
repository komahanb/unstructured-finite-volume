!=====================================================================!
! THE STRUCTURAL INTERPRETER - earned at LEVEL 4, and at no level
! below it.
!
! The smallest thing that can turn a typed relation into something a
! person can look at. It is a TEST FIXTURE and not production: the
! tower is finding out whether the relational nucleus already carries
! everything a picture needs, and answering that question by writing
! a renderer that may use nothing but the nucleus.
!
!                    WHAT THIS RENDERER IS NOT ALLOWED TO KNOW
!
! It does not know the specimen. It has never heard of a, p, u or m;
! it does not know that D1 has four columns, or that its first column
! is full, or that X0 comes before X1. Every one of those facts is
! obtained, at the moment of drawing, by asking:
!
!      domain(1), domain(2)          which carriers this relation
!                                    relates
!      sets % num_members_of(carrier)              how many rows, how many columns
!      sets % member_of(carrier, i)           WHICH member stands at position i
!      relation % has([col, row])    whether this cell is filled
!      relation % name()             what to call the picture
!      label_for(carrier, member, labels)    what the reader calls a member
!
! and nothing else. Point it at a relation it has never seen and it
! draws that one instead.
!
!                    ORDER COMES FROM DECLARATION, ALWAYS
!
! Every axis is walked as member(1), member(2), ..., member(n). Never
! sorted, never assumed to be 1..n, never read off the tuple table -
! which is a set and has no order to read. A carrier that declares
! its members { 30, 10, 20 } is drawn 30, 10, 20, and Level 4 checks
! exactly that with a carrier built to be hostile.
!
! The same walk serves both representations, so the sparsity picture
! and the dependency listing cannot disagree: they ask the identical
! question, cell by cell, in the identical order.
!
!                       THE PICTURE IS DOWNSTREAM
!
! The flow is one-directional and this file sits at its end:
!
!      primitive incidence -> relation algebra -> dependency
!                                                     -> representation
!
! Nothing here computes a dependency, and nothing here may be
! consulted about one. A representation that could change what is
! true would not be a representation.
!
!                        THE CHAIN VALIDATES ITSELF
!
! chain_picture refuses to draw a chain that does not compose: stage
! k's codomain must be stage k+1's domain, by structural identity. A
! picture that could show X0 -> X1 -> X2 for relations that do not
! meet would be a lie with arrows in it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module structural_renderer_fixture

  use graph_fractal        , only : graph
  use map_set        , only : set_map
  use map_label      , only : label_map
  use relation_finitary                , only : relation
  use visualization_carriers_fixture, only : label_for

  implicit none

  private
  public :: picture, stage
  public :: sparsity_picture, chain_picture, dependency_listing
  public :: glyph_at

  integer, parameter :: LINE = 132

  !-------------------------------------------------------------------!
  ! The two glyphs, and the two field widths. A cell is filled or it
  ! is not; there is no third answer, because there is no coefficient
  ! to be small.
  !-------------------------------------------------------------------!

  character, parameter :: FILLED = '#'
  character, parameter :: EMPTY  = '.'

  !-------------------------------------------------------------------!
  ! The narrowest the row-label stub may be. A wider label widens it;
  ! nothing narrows it, so a one-letter carrier and a three-letter one
  ! line their grids up in the same column.
  !
  ! Named MIN_STUB rather than STUB because Fortran does not
  ! distinguish case, and a local named `stub` would silently shadow
  ! it - which is exactly the bug this comment was written after.
  !-------------------------------------------------------------------!

  integer  , parameter :: MIN_STUB = 7

  !===================================================================!
  ! A finished representation: lines of text, and the two things a
  ! test wants to do with them - look at one, and print them all.
  !===================================================================!

  type :: picture

     character(len=LINE), allocatable :: line(:)

   contains

     procedure :: emit
     procedure :: rows => picture_rows
     procedure :: at   => picture_at

  end type picture

  !===================================================================!
  ! One leg of a chain. A Fortran array carries one dynamic type, so
  ! a chain of relations of assorted concretions sits in wrappers -
  ! the same device the nucleus uses for signature slots and graph
  ! seats.
  !===================================================================!

  type :: stage

     class(relation), allocatable :: leg

  end type stage

  interface stage
     module procedure make_stage
  end interface stage

contains

  type(stage) function make_stage(leg) result(this)

    class(relation), intent(in) :: leg

    allocate(this % leg, source=leg)

  end function make_stage

  !===================================================================!
  ! REPRESENTATION C/D/E - the sparsity picture of one dependency.
  !
  !      rows    = codomain, in ITS declaration order
  !      columns = domain,   in ITS declaration order
  !      glyph   = has([column member, row member])
  !
  ! The name line comes from the relation itself, so a picture always
  ! says which relation it is a picture of.
  !===================================================================!

  type(picture) function sparsity_picture(r, sets, labels) result(pic)

    class(relation), intent(in) :: r
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    type(graph) :: cols, rows
    integer                        :: stub, wide, i, j, at

    if (r % arity() .ne. 2) then
       error stop 'structural_renderer_fixture: a sparsity picture reads a binary relation'
    end if

    cols = r % domain(1)
    rows = r % domain(2)

    stub = max(widest(rows, sets, labels), MIN_STUB) + 2
    wide = widest(cols, sets, labels) + 1

    allocate(pic % line(2 + sets % num_members_of(rows)))
    pic % line = repeat(' ', LINE)

    call put(pic % line(1), 1, r % name())

    ! The header: the domain's members, at their own positions.
    do j = 1, sets % num_members_of(cols)
       at = stub + (j - 1) * wide
       call put(pic % line(2), at, label_for(cols, sets % member_of(cols, j), labels))
    end do

    ! One line per member of the codomain, in its own order.
    do i = 1, sets % num_members_of(rows)
       call put(pic % line(2 + i), 1, label_for(rows, sets % member_of(rows, i), labels))
       do j = 1, sets % num_members_of(cols)
          at = stub + (j - 1) * wide
          call put(pic % line(2 + i), at, &
               &   glyph_at(r, sets % member_of(cols, j), sets % member_of(rows, i)))
       end do
    end do

  end function sparsity_picture

  !===================================================================!
  ! Is this cell filled. The one question a picture asks about
  ! content, and it is the relation's own membership test.
  !===================================================================!

  character function glyph_at(r, from, to)

    class(relation), intent(in) :: r
    integer        , intent(in) :: from, to

    if (r % has([from, to])) then
       glyph_at = FILLED
    else
       glyph_at = EMPTY
    end if

  end function glyph_at

  !===================================================================!
  ! REPRESENTATION A/B - the chain, forward or reversed.
  !
  !      X0 --D1--> X1 --D2--> X2 --D3--> X3
  !
  ! Every carrier name and every arrow label is read off the stages;
  ! the reverse chain is not this line with its words swapped, it is
  ! this same routine handed the transposed relations.
  !
  ! REFUSED: a chain whose legs do not meet. Stage k's codomain must
  ! be stage k+1's domain, by structural identity - never by size,
  ! which would happily let X1 hand over to X2.
  !===================================================================!

  type(picture) function chain_picture(title, legs, labels) result(pic)

    character(len=*), intent(in) :: title
    type(stage)     , intent(in) :: legs(:)
    type(label_map), intent(in) :: labels

    type(graph) :: here, there
    character(len=:) , allocatable :: text
    integer                        :: k

    if (size(legs) .lt. 1) then
       error stop 'structural_renderer_fixture: a chain has at least one leg'
    end if

    do k = 1, size(legs)
       if (legs(k) % leg % arity() .ne. 2) then
          error stop 'structural_renderer_fixture: a chain is drawn from binary relations'
       end if
    end do

    do k = 1, size(legs) - 1
       here  = legs(k)     % leg % domain(2)
       there = legs(k + 1) % leg % domain(1)
       if (.not. here % same_as(there)) then
          error stop 'structural_renderer_fixture: this chain does not compose'
       end if
    end do

    here = legs(1) % leg % domain(1)
    text = carrier_name(here, labels)
    do k = 1, size(legs)
       there = legs(k) % leg % domain(2)
       text  = text // ' --' // legs(k) % leg % name() // '--> ' // &
            &  carrier_name(there, labels)
    end do

    allocate(pic % line(2))
    pic % line = repeat(' ', LINE)
    call put(pic % line(1), 1, title)
    call put(pic % line(2), 1, text)

  end function chain_picture

  !===================================================================!
  ! REPRESENTATION A, the deterministic form. Where ASCII geometry
  ! gets fragile, the semantic content does not:
  !
  !      a -> p
  !      b -> p q
  !
  ! one line per member of the domain, in declaration order, listing
  ! what it reaches in the codomain's declaration order. Same walk,
  ! same question, same order as the sparsity picture - so the two
  ! representations cannot tell different stories.
  !===================================================================!

  type(picture) function dependency_listing(r, sets, labels) result(pic)

    class(relation), intent(in) :: r
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    type(graph) :: cols, rows
    character(len=:) , allocatable :: text
    integer                        :: i, j
    logical                        :: any_reached

    if (r % arity() .ne. 2) then
       error stop 'structural_renderer_fixture: a listing reads a binary relation'
    end if

    cols = r % domain(1)
    rows = r % domain(2)

    allocate(pic % line(1 + sets % num_members_of(cols)))
    pic % line = repeat(' ', LINE)

    call put(pic % line(1), 1, r % name())

    do j = 1, sets % num_members_of(cols)
       text        = label_for(cols, sets % member_of(cols, j), labels) // ' ->'
       any_reached = .false.
       do i = 1, sets % num_members_of(rows)
          if (r % has([sets % member_of(cols, j), sets % member_of(rows, i)])) then
             text        = text // ' ' // label_for(rows, sets % member_of(rows, i), labels)
             any_reached = .true.
          end if
       end do
       if (.not. any_reached) text = text // ' (nothing)'
       call put(pic % line(1 + j), 1, text)
    end do

  end function dependency_listing

  !===================================================================!
  ! Looking at a picture, and printing one.
  !===================================================================!

  pure integer function picture_rows(this)

    class(picture), intent(in) :: this

    picture_rows = size(this % line)

  end function picture_rows

  function picture_at(this, k) result(text)

    class(picture), intent(in) :: this
    integer       , intent(in) :: k

    character(len=:), allocatable :: text

    text = trim(this % line(k))

  end function picture_at

  subroutine emit(this)

    class(picture), intent(in) :: this

    integer :: k

    do k = 1, size(this % line)
       write(*,'(4x,a)') trim(this % line(k))
    end do

  end subroutine emit

  !===================================================================!
  ! Small mechanics. The widest label a carrier will produce, found
  ! by asking about every member it has - not by guessing from its
  ! size.
  !===================================================================!

  integer function widest(carrier, sets, labels)

    type(graph), intent(in) :: carrier
    type(set_map)  , intent(in) :: sets
    type(label_map), intent(in) :: labels

    integer :: k

    widest = 1
    do k = 1, sets % num_members_of(carrier)
       widest = max(widest, len(label_for(carrier, sets % member_of(carrier, k), labels)))
    end do

  end function widest

  function carrier_name(carrier, labels) result(text)

    type(graph), intent(in) :: carrier
    type(label_map), intent(in) :: labels

    character(len=:), allocatable :: text

    text = labels % label_of(carrier)
    if (len(text) .eq. 0) text = '?'

  end function carrier_name

  subroutine put(line, at, text)

    character(len=*), intent(inout) :: line
    integer         , intent(in)    :: at
    character(len=*), intent(in)    :: text

    if (at + len(text) - 1 .gt. len(line)) then
       error stop 'structural_renderer_fixture: the picture is wider than its page'
    end if

    line(at : at + len(text) - 1) = text

  end subroutine put

end module structural_renderer_fixture
