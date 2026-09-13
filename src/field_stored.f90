!=====================================================================!
! The concrete field: values over a domain.
!
! One concrete type serves every field in the tower. Its domain is a
! set GRAPH, and the domain's identity is the only datum that ever
! distinguishes a cell field from a face field; the field stores no
! side flag. Because there is exactly one concrete field, a plain
! Fortran array can store a collection of them.
!
!                  WHAT THE FIELD STORES OF ITS DOMAIN
!
!      type(graph) :: graph       which set        O(1)
!      integer     :: num_entries    how many         O(1)
!
! and nothing else. An earlier version stored a COPY of the whole
! domain object, which for a listed domain meant a copy of the member
! list: 40 fields on a 200 000-member domain stored 28.7 MB of
! duplicated extension, measured, against 30.5 MB predicted for
! exactly that duplication. The extension is now stored once, in
! whatever set map the caller owns.
!
! No information was lost. A field only ever returned two values
! about its domain - WHICH and HOW MANY - and both are retained by
! value. The copy already froze the count, so freezing it explicitly
! changes no behaviour; it only stops the copy from costing
! O(N_extent).
!
!=====================================================================!
!
!                        THE VALUE-KIND RULE
!
! A field stores one kind of value at a time, in the one store every
! field inherits from field_calculus, where the ten adapters are
! written once. From that, three rules that hold for all of them:
!
!      check first   a caller checks value_kind() before reading a
!                    vector
!
!      wrong getter  returns a zero-length array. No conversion and
!                    no inference happens, and a pure procedure has
!                    no error path. The zero-length result is the
!                    indicator
!
!      any setter    replaces both the values and the kind. Setting
!                    reals onto a field that stored integers makes it
!                    a real field
!
! No conversion happens anywhere. A field that stores boundary names
! does not return them as numbers.
!
!=====================================================================!
!
!                       WHERE A VALUE IS STORED
!
! A field stores its values in the order the domain lists its
! members, and stores the components of one member contiguously:
!
!      member          7        7        3        3
!      component       1        2        1        2
!                   +--------+--------+--------+--------+
!      values       |  v(1)  |  v(2)  |  v(3)  |  v(4)  |
!                   +--------+--------+--------+--------+
!
!      position = (entry_position - 1) * num_components + component
!
! Everything that reads a flat vector out of a field depends on this -
! a linear solver, a file writer, a matrix adapter. It is the reason
! this library needs no degree-of-freedom index map of its own.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module field_stored

  use util_precision  , only : dp
  use field_calculus, only : field
  use graph_fractal , only : graph
  use view_directed , only : directed_graph
  use token_identity, only : token

  implicit none

  private
  public :: stored_field, typed_field_domain

  !===================================================================!
  ! One field: the description and the store the abstract field
  ! contains, stated by the constructor and nothing more.
  !===================================================================!

  type, extends(field) :: stored_field
   contains
     procedure :: assign_in => stored_field_assign_in
  end type stored_field

  !===================================================================!
  ! The domain a family of fields shares. A primal state, a tangent,
  ! a costate, a residual, and a functional's state differ in meaning,
  ! not in the rule for choosing the graph they live on. This value
  ! states that graph identity and shape once.
  !===================================================================!

  type :: typed_field_domain

     type(graph), private :: graph
     integer    , private :: ne = 0
     integer    , private :: nc = 1

   contains

     procedure :: domain         => typed_field_domain_graph
     procedure :: num_entries    => typed_field_domain_num_entries
     procedure :: num_components => typed_field_domain_num_components

     procedure :: real_field       => typed_field_domain_real_field
     procedure :: state            => typed_field_domain_state
     procedure :: functional_state => typed_field_domain_functional_state
     procedure :: design           => typed_field_domain_design
     procedure :: direction        => typed_field_domain_direction
     procedure :: tangent          => typed_field_domain_tangent
     procedure :: costate          => typed_field_domain_costate
     procedure :: residual         => typed_field_domain_residual
     procedure :: forcing          => typed_field_domain_forcing
     procedure :: solution         => typed_field_domain_solution

  end type typed_field_domain

  !===================================================================!
  ! Constructor. Name the field, state its domain, state how many
  ! components each entry has. The values are set afterwards through
  ! a setter, which is also what fixes the kind.
  !===================================================================!

  interface stored_field
     module procedure create
  end interface stored_field

  interface typed_field_domain
     module procedure create_domain
     module procedure create_placed
  end interface typed_field_domain

contains

  !===================================================================!
  ! Build an empty field on a domain. The domain's identity states
  ! whether this is a cell field or a face field; the field does not
  ! store the fact a second time.
  !===================================================================!

  type(stored_field) function create(label, domain, num_entries, num_components, unit_name) &
       & result(this)

    character(len=*), intent(in)           :: label
    type(graph) , intent(in)           :: domain
    integer         , intent(in)           :: num_entries
    integer         , intent(in), optional :: num_components
    character(len=*), intent(in), optional :: unit_name

    call this % describe(label, domain, num_entries, num_components, unit_name)

  end function create

  !===================================================================!
  ! Build the common domain of a family of fields on the vertex set
  ! of a directed graph: the identity is the graph's vertex set and
  ! the extent is its vertex count, so no caller states an extent
  ! read from an array. The stored_field setter refuses a value
  ! vector that does not fill entries times components.
  !===================================================================!

  type(typed_field_domain) function create_placed(points, num_components) result(this)

    class(directed_graph), intent(in)           :: points
    integer              , intent(in), optional :: num_components

    if (present(num_components)) then
       this = create_domain(points % vertex_set(), points % num_vertices(), num_components)
    else
       this = create_domain(points % vertex_set(), points % num_vertices())
    end if

  end function create_placed

  !===================================================================!
  ! Build the common domain of a family of fields from an identity
  ! and its extent. A graph value is an identity without a
  ! cardinality, so this form is for an owner that stores the declared
  ! cardinality of that same identity (a discrete domain, a residual's
  ! unknown and point sets, a linearization's image read from the
  ! statement's result); a caller with the directed graph places
  ! through it instead.
  !===================================================================!

  type(typed_field_domain) function create_domain(domain, num_entries, num_components) result(this)

    type(graph), intent(in)           :: domain
    integer    , intent(in)           :: num_entries
    integer    , intent(in), optional :: num_components
    type(token) :: identity
    character(len=150) :: message

    identity = domain % id()
    if (.not. identity % declared()) then
       error stop 'field_stored: create_domain requires a declared graph, but this domain''s &
            &identity is undeclared'
    end if
    if (num_entries < 0) then
       write(message,'(a,i0)') 'field_stored: create_domain requires a nonnegative extent; &
            &num_entries = ', num_entries
       error stop trim(message)
    end if

    this % graph = domain
    this % ne    = num_entries
    this % nc    = 1
    if (present(num_components)) this % nc = num_components
    if (this % nc < 1) then
       write(message,'(a,i0)') 'field_stored: create_domain requires at least one component; &
            &num_components = ', this % nc
       error stop trim(message)
    end if

  end function create_domain

  type(graph) function typed_field_domain_graph(this) result(domain)

    class(typed_field_domain), intent(in) :: this

    domain = this % graph

  end function typed_field_domain_graph

  pure integer function typed_field_domain_num_entries(this) result(num_entries)

    class(typed_field_domain), intent(in) :: this

    num_entries = this % ne

  end function typed_field_domain_num_entries

  pure integer function typed_field_domain_num_components(this) result(num_components)

    class(typed_field_domain), intent(in) :: this

    num_components = this % nc

  end function typed_field_domain_num_components

  type(stored_field) function typed_field_domain_real_field(this, label, values, unit_name) &
       & result(data)

    class(typed_field_domain), intent(in)           :: this
    character(len=*)    , intent(in)           :: label
    real(dp)            , intent(in)           :: values(:)
    character(len=*)    , intent(in), optional :: unit_name

    if (present(unit_name)) then
       data = stored_field(label, this % graph, this % ne, &
            & num_components=this % nc, unit_name=unit_name)
    else
       data = stored_field(label, this % graph, this % ne, &
            & num_components=this % nc)
    end if
    call data % set_real_vector(values)

  end function typed_field_domain_real_field

  type(stored_field) function typed_field_domain_state(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('state', values)

  end function typed_field_domain_state

  type(stored_field) function typed_field_domain_functional_state(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('functional state', values)

  end function typed_field_domain_functional_state

  type(stored_field) function typed_field_domain_design(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('design', values)

  end function typed_field_domain_design

  type(stored_field) function typed_field_domain_direction(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('direction', values)

  end function typed_field_domain_direction

  type(stored_field) function typed_field_domain_tangent(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('tangent', values)

  end function typed_field_domain_tangent

  type(stored_field) function typed_field_domain_costate(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('costate', values)

  end function typed_field_domain_costate

  type(stored_field) function typed_field_domain_residual(this, values, label) result(data)

    class(typed_field_domain), intent(in)           :: this
    real(dp)            , intent(in)           :: values(:)
    character(len=*)    , intent(in), optional :: label

    if (present(label)) then
       data = this % real_field(label, values)
    else
       data = this % real_field('residual', values)
    end if

  end function typed_field_domain_residual

  type(stored_field) function typed_field_domain_forcing(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('forcing', values)

  end function typed_field_domain_forcing

  type(stored_field) function typed_field_domain_solution(this, values) result(data)

    class(typed_field_domain), intent(in) :: this
    real(dp)            , intent(in) :: values(:)

    data = this % real_field('solution', values)

  end function typed_field_domain_solution

  !===================================================================!
  ! Place this value at a location that is a stored field. A location
  ! of any other type is an error: the caller requested a copy the
  ! location cannot store.
  !===================================================================!

  subroutine stored_field_assign_in(this, location)

    class(stored_field), intent(in)    :: this
    class(field)       , intent(inout) :: location

    select type (location)
    type is (stored_field)
       location = this
    class default
       error stop 'field_stored: stored_field_assign_in requires the location to be a &
            &stored_field, but the actual location has a different dynamic type'
    end select

  end subroutine stored_field_assign_in

end module field_stored
