!=====================================================================!
! The field suite: the laws of the one field on the one domain kind
! (level 5). A field lives on a set graph - ambient, carved,
! nested, or empty - stores by domain enumeration, and carries its
! domain's identity wherever it is copied.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_field

  use iso_fortran_env  , only : dp => REAL64
  use graph_fractal           , only : graph
  use map_set_representation, only : counted_set_representation, &
       & listed_set_representation
  use map_set           , only : set_map
  use map_inclusion     , only : inclusion_map, declared_subobject
  use field_stored, only : field

  implicit none

  type(graph)         :: cells, walls, hot, none
  type(set_map)       :: sets
  type(inclusion_map) :: inclusions
  type(field)         :: q, w, h, z, copy
  type(graph)         :: dom
  real(dp), allocatable :: v(:)
  integer :: kk
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph field suite (level 5)"
  write(*,'(1x,a)') "============================================="

  call cells % declare()
  call sets % bind(cells, counted_set_representation(6))

  call walls % declare()
  call sets       % bind(walls, listed_set_representation([5, 2, 6]))
  call inclusions % include_in(walls, cells)

  call hot % declare()
  call sets       % bind(hot, listed_set_representation([6, 2]))
  call inclusions % include_in(hot, walls)

  call none % declare()
  call sets       % bind(none, listed_set_representation([integer ::]))
  call inclusions % include_in(none, cells)

  q = field('q', cells, 6, ncomp=2)
  call q % set_real_vector([(1.0_dp * kk, kk = 1, 12)])
  dom = q % domain()
  call report(dom % same_as(cells), &
       & "an ambient field's domain is the carrier, by identity", nfail)
  call report(q % num_entries() .eq. 6, &
       & "entries count the domain", nfail)
  call q % get_real_vector(v)
  call report(abs(v((4 - 1) * 2 + 2) - 8.0_dp) < 1.0d-14, &
       & "values interleave at (entry-1)*ncomp + component", nfail)

  w = field('w', walls, 3)
  call w % set_real_vector([50.0_dp, 20.0_dp, 60.0_dp])
  call w % get_real_vector(v)
  call report(abs(v(sets % index_in(walls, 2)) - 20.0_dp) < 1.0d-14 .and. &
       &      abs(v(1) - 50.0_dp) < 1.0d-14, &
       & "a subset field stores by declaration, not by member value", nfail)

  h = field('h', hot, 2)
  call h % set_real_vector([600.0_dp, 200.0_dp])
  dom = h % domain()
  call report(declared_subobject(dom, cells, inclusions), &
       & "a nested-subset field still descends from the ground", nfail)

  z = field('z', none, 0)
  call z % set_real_vector([real(dp) ::])
  call report(z % num_entries() .eq. 0, &
       & "the empty-domain field is lawful, zero entries", nfail)

  copy = w
  dom = copy % domain()
  call report(dom % same_as(walls), &
       & "a copy carries its domain's identity along", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all field checks passed"
  else
     error stop
  end if

contains

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

  ! silence unused warnings for k in implied do
end program test_graph_field
