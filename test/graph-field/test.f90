!=====================================================================!
! The field suite: the laws of the one field on the one domain kind
! (level 5). A field lives on a member_set - ambient, subset,
! nested, or empty - stores by domain enumeration, and carries its
! domain's identity wherever it is copied.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_graph_field

  use iso_fortran_env  , only : dp => REAL64
  use graph_carrier    , only : counted_set, subset_set, member_set
  use class_graph_field, only : field

  implicit none

  type(counted_set) :: cells
  type(subset_set)  :: walls, hot, none
  type(field)       :: q, w, h, z, copy
  class(member_set), allocatable :: dom
  real(dp), allocatable :: v(:)
  integer :: kk
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "graph field suite (level 5)"
  write(*,'(1x,a)') "============================================="

  cells = counted_set('cells', 6)
  walls = subset_set('walls', cells, [5, 2, 6])
  hot   = subset_set('hot'  , walls, [6, 2])
  none  = subset_set('none' , cells, [integer ::])

  q = field('q', cells, ncomp=2)
  call q % set_real_vector([(1.0_dp * kk, kk = 1, 12)])
  call q % domain(dom)
  call report(dom % same_as(cells), &
       & "an ambient field's domain is the carrier, by identity", nfail)
  call report(q % num_entries() .eq. 6, &
       & "entries count the domain", nfail)
  call q % get_real_vector(v)
  call report(abs(v((4 - 1) * 2 + 2) - 8.0_dp) < 1.0d-14, &
       & "values interleave at (entry-1)*ncomp + component", nfail)

  w = field('w', walls)
  call w % set_real_vector([50.0_dp, 20.0_dp, 60.0_dp])
  call w % get_real_vector(v)
  call report(abs(v(walls % local_index(2)) - 20.0_dp) < 1.0d-14 .and. &
       &      abs(v(1) - 50.0_dp) < 1.0d-14, &
       & "a subset field stores by declaration, not by member value", nfail)

  h = field('h', hot)
  call h % set_real_vector([600.0_dp, 200.0_dp])
  call h % domain(dom)
  call report(dom % is_subobject_of(cells), &
       & "a nested-subset field still descends from the ground", nfail)

  z = field('z', none)
  call z % set_real_vector([real(dp) ::])
  call report(z % num_entries() .eq. 0, &
       & "the empty-domain field is lawful, zero entries", nfail)

  copy = w
  call copy % domain(dom)
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
