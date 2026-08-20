!=====================================================================!
! The mesh: a stored graph that lives in space.
!
! LEVEL 1 OF THE STRATIFICATION - for measurements are calculus
! content: values riding on structure, no goal anywhere in them.
! This is the tower's one inheritance crossing: the mesh IS a graph
! - cells are the vertices, interior faces the edges, boundary
! faces the edges without heads - so it extends the stored graph
! rather than holding one.
!
!     +-----+-----+
!     |  1  |  2  |            (1)---(2)
!     +-----+-----+             |     |        the same mesh, seen
!     |  3  |  4  |            (3)---(4)       as its graph
!     +-----+-----+
!
! Everything else is measurement hung on that graph, carried as typed
! fields with compiled names:
!
!      cell_volume()    one number per cell
!      cell_centre()    three per cell
!      face_area()      one per face
!      face_delta()     the centre-to-centre distance, one per face
!      face_normal()    three per face
!      face_centre()    three per face
!      face_weights()   the interpolation weight, one per face
!
! No string names any of these. The dictionary in
! geometry-to-operator-mapping.md says which operator argument each
! one feeds; an operator receives those numbers at construction and
! never reads the mesh.
!
! THE GATE AT LOAD. Every geometry array must match the structure it
! measures - one volume per cell, one area per face - and the
! constructor stops the program on a mismatch rather than storing a
! lie. A mesh that loads is a mesh whose measurements fit.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh

  use iso_fortran_env    , only : dp => REAL64, error_unit
  use field_stored  , only : stored_field
  use view_directed_stored        , only : stored_directed_graph
  use graph_fractal      , only : graph

  implicit none

  private
  public :: mesh

  !===================================================================!
  ! One mesh: the inherited structure, plus seven measurements.
  !===================================================================!

  type, extends(stored_directed_graph) :: mesh

     type(stored_field) :: volumes
     type(stored_field) :: cell_centres
     type(stored_field) :: areas
     type(stored_field) :: deltas
     type(stored_field) :: normals
     type(stored_field) :: face_centers_
     type(stored_field) :: weights

   contains

     procedure :: cell_volume
     procedure :: cell_centre
     procedure :: face_area
     procedure :: face_delta
     procedure :: face_normal
     procedure :: face_centre
     procedure :: face_weights

  end type mesh

  interface mesh
     module procedure create
  end interface mesh

contains

  !===================================================================!
  ! Build a mesh from its structure and its measurements. The
  ! structure arguments are the stored graph's own; the geometry
  ! arrives as plain arrays, one entry per cell or per face, vector
  ! quantities three wide in entry order. The gate: every array must
  ! fit the structure, or the constructor stops.
  !===================================================================!

  impure type(mesh) function create(nv, tails, heads, volumes, &
       & cell_centres, areas, deltas, normals, face_centres, weights, &
       & vtags, etags, number) result(this)

    integer         , intent(in)           :: nv
    integer         , intent(in)           :: tails(:)
    integer         , intent(in)           :: heads(:)
    real(dp)        , intent(in)           :: volumes(:)
    real(dp)        , intent(in)           :: cell_centres(:)
    real(dp)        , intent(in)           :: areas(:)
    real(dp)        , intent(in)           :: deltas(:)
    real(dp)        , intent(in)           :: normals(:)
    real(dp)        , intent(in)           :: face_centres(:)
    real(dp)        , intent(in)           :: weights(:)
    character(len=*), intent(in), optional :: vtags(:)
    character(len=*), intent(in), optional :: etags(:)
    integer         , intent(in), optional :: number

    type(graph) :: cells, faces
    integer :: ne

    ! The structure first, through the parent's own constructor.
    this % stored_directed_graph = stored_directed_graph(nv, tails=tails, heads=heads, &
         & vtags=vtags, etags=etags, number=number)

    ne = this % num_edges()

    ! The gate: measurements must fit the structure they measure.
    call gate(size(volumes)      == nv    , 'one volume per cell')
    call gate(size(cell_centres) == 3 * nv, 'three centre parts per cell')
    call gate(size(areas)        == ne    , 'one area per face')
    call gate(size(deltas)       == ne    , 'one delta per face')
    call gate(size(normals)      == 3 * ne, 'three normal parts per face')
    call gate(size(face_centres) == 3 * ne, 'three centre parts per face')
    call gate(size(weights)      == ne    , 'one weight per face')

    ! Geometry rides the graph's OWN carriers, so a field's domain
    ! answers the mesh identity every consumer will ask about.
    cells = this % vertex_set()
    faces = this % edge_set()

    this % volumes       = stored_field('cell_volume' , cells, nv, unit_name='m3')
    this % cell_centres  = stored_field('cell_centre' , cells, nv, num_components=3, unit_name='m')
    this % areas         = stored_field('face_area'   , faces, ne, unit_name='m2')
    this % deltas        = stored_field('face_delta'  , faces, ne, unit_name='m')
    this % normals       = stored_field('face_normal' , faces, ne, num_components=3, unit_name='-')
    this % face_centers_ = stored_field('face_centre' , faces, ne, num_components=3, unit_name='m')
    this % weights       = stored_field('face_weights', faces, ne, unit_name='-')

    call this % volumes       % set_real_vector(volumes)
    call this % cell_centres  % set_real_vector(cell_centres)
    call this % areas         % set_real_vector(areas)
    call this % deltas        % set_real_vector(deltas)
    call this % normals       % set_real_vector(normals)
    call this % face_centers_ % set_real_vector(face_centres)
    call this % weights       % set_real_vector(weights)

  end function create

  !===================================================================!
  ! The gate itself: state what failed, then stop. A mesh with wrong
  ! measurements must not exist.
  !===================================================================!

  subroutine gate(fits, what)

    logical         , intent(in) :: fits
    character(len=*), intent(in) :: what

    if (fits) return

    write(error_unit, *) 'mesh gate: expected ', what
    error stop 'mesh: a measurement does not fit the structure'

  end subroutine gate

  !===================================================================!
  ! The seven answers. Each is a copy of the stored field, so a
  ! caller may read it, hand its values to an operator, and never
  ! reach back into the mesh.
  !===================================================================!

  type(stored_field) function cell_volume(this)

    class(mesh), intent(in) :: this

    cell_volume = this % volumes

  end function cell_volume

  type(stored_field) function cell_centre(this)

    class(mesh), intent(in) :: this

    cell_centre = this % cell_centres

  end function cell_centre

  type(stored_field) function face_area(this)

    class(mesh), intent(in) :: this

    face_area = this % areas

  end function face_area

  type(stored_field) function face_delta(this)

    class(mesh), intent(in) :: this

    face_delta = this % deltas

  end function face_delta

  type(stored_field) function face_normal(this)

    class(mesh), intent(in) :: this

    face_normal = this % normals

  end function face_normal

  type(stored_field) function face_centre(this)

    class(mesh), intent(in) :: this

    face_centre = this % face_centers_

  end function face_centre

  type(stored_field) function face_weights(this)

    class(mesh), intent(in) :: this

    face_weights = this % weights

  end function face_weights

end module view_mesh
