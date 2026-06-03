!=====================================================================!
! Read a run from a plain text file - the library stays problem free,
! the config decides the problem. One keyword per line:
!
!   mesh      test/box-3.msh
!   equation  diffusion
!   kappa     1.0
!   source    0.0
!   boundary  front  dirichlet 5.0
!   boundary  wall   neumann   0.0
!   boundary  skin   robin     1.0 1.0 2.0     ! a*phi + b*dphi/dn = c
!   solver    cg                               ! cg | sor | gs | gj
!   omega     1.5
!   max_tol   1.0e-12
!   max_it    200
!   output    box.vtu
!
! Blank lines and lines starting with '#' are ignored.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_config

  use iso_fortran_env, only : dp => REAL64
  use class_file     , only : file
  use class_string   , only : string

  implicit none

  private
  public :: config

  type :: config

     type(string) :: meshfile
     type(string) :: equation         ! "diffusion"
     real(dp)     :: kappa = 1.0_dp
     real(dp)     :: source = 0.0_dp

     ! Boundary specs (parallel arrays, nbc entries)
     integer                   :: nbc = 0
     type(string), allocatable :: bc_name(:)
     type(string), allocatable :: bc_kind(:)   ! dirichlet|neumann|robin
     real(dp)    , allocatable :: bc_a(:), bc_b(:), bc_c(:)

     type(string) :: solver           ! cg | sor | gs | gj
     real(dp)     :: omega   = 1.5_dp
     real(dp)     :: max_tol = 1.0d-12
     integer      :: max_it  = 200

     ! Transient marching - dt > 0 turns it on (backward euler)
     real(dp)     :: dt     = 0.0_dp
     real(dp)     :: tinit  = 0.0_dp
     real(dp)     :: tfinal = 0.0_dp

     type(string) :: output

   contains

     procedure :: print

  end type config

  interface config
     module procedure create
  end interface config

contains

  !===================================================================!
  ! Read a config from a text file
  !===================================================================!

  type(config) function create(filename) result(this)

    character(len=*), intent(in) :: filename

    type(file)                :: cfile
    type(string), allocatable :: lines(:)
    type(string), allocatable :: tok(:)
    integer :: iline, ntok, ib

    ! Defaults
    this % equation = string("diffusion")
    this % solver   = string("cg")
    this % output   = string("output.vtu")

    cfile = file(filename, 256)
    call cfile % read_lines(lines)

    ! First pass: count boundary lines so the bc arrays can be sized
    this % nbc = 0
    do iline = 1, size(lines)
       call split_words(lines(iline) % str, ntok, tok)
       if (ntok .lt. 1) cycle
       if (trim(tok(1) % str) .eq. "boundary") this % nbc = this % nbc + 1
    end do

    allocate(this % bc_name(this % nbc), this % bc_kind(this % nbc))
    allocate(this % bc_a(this % nbc), this % bc_b(this % nbc), this % bc_c(this % nbc))

    ! Second pass: fill
    ib = 0
    do iline = 1, size(lines)

       call split_words(lines(iline) % str, ntok, tok)
       if (ntok .lt. 1) cycle
       if (trim(tok(1) % str(1:1)) .eq. "#") cycle

       select case (trim(tok(1) % str))

       case ("mesh");     this % meshfile = string(trim(tok(2) % str))
       case ("equation"); this % equation = string(trim(tok(2) % str))
       case ("kappa");    this % kappa  = tok(2) % asreal()
       case ("source");   this % source = tok(2) % asreal()
       case ("solver");   this % solver = string(trim(tok(2) % str))
       case ("omega");    this % omega   = tok(2) % asreal()
       case ("max_tol");  this % max_tol = tok(2) % asreal()
       case ("max_it");   this % max_it  = tok(2) % asinteger()
       case ("dt");       this % dt      = tok(2) % asreal()
       case ("t_init");   this % tinit   = tok(2) % asreal()
       case ("t_final");  this % tfinal  = tok(2) % asreal()
       case ("output");   this % output  = string(trim(tok(2) % str))

       case ("boundary")
          ib = ib + 1
          this % bc_name(ib) = string(trim(tok(2) % str))
          this % bc_kind(ib) = string(trim(tok(3) % str))
          this % bc_a(ib) = 0.0_dp; this % bc_b(ib) = 0.0_dp; this % bc_c(ib) = 0.0_dp
          select case (trim(tok(3) % str))
          case ("dirichlet")   ! a*phi = c
             this % bc_a(ib) = 1.0_dp
             this % bc_c(ib) = tok(4) % asreal()
          case ("neumann")     ! b*dphi/dn = c
             this % bc_b(ib) = 1.0_dp
             this % bc_c(ib) = tok(4) % asreal()
          case ("robin")       ! a*phi + b*dphi/dn = c
             this % bc_a(ib) = tok(4) % asreal()
             this % bc_b(ib) = tok(5) % asreal()
             this % bc_c(ib) = tok(6) % asreal()
          case default
             print *, "config: unknown bc kind '", trim(tok(3) % str), "'"
             error stop
          end select

       case default
          ! ignore unknown keywords (comments handled above)
       end select

    end do

    ! Minimal validation
    if (.not. allocated(this % meshfile % str)) then
       print *, "config: no 'mesh <file>' line in ", trim(filename)
       error stop
    end if
    if (this % nbc .eq. 0) then
       print *, "config: warning - no boundary conditions specified (all default to phi=0)"
    end if

  end function create

  !===================================================================!
  ! Split a line into whitespace-delimited words (collapses runs of
  ! spaces/tabs, so alignment in the config file is harmless).
  !===================================================================!

  subroutine split_words(str, nw, words)
    character(len=*)          , intent(in)  :: str
    integer                   , intent(out) :: nw
    type(string), allocatable , intent(out) :: words(:)
    type(string) :: tmp(64)
    integer      :: i, n, start
    logical      :: inword, inquote
    character(len=1), parameter :: tab = char(9)
    character(len=1), parameter :: quote = '"'
    n = len_trim(str)
    nw = 0; inword = .false.; inquote = .false.; start = 1
    do i = 1, n
       if (str(i:i) .eq. quote) then
          ! A double-quoted span is one word (quotes stripped), so names
          ! with spaces work, e.g. boundary "cylindrical wall" dirichlet 5
          if (.not. inquote) then
             inquote = .true.; inword = .true.; start = i + 1
          else
             inquote = .false.; inword = .false.
             nw = nw + 1; tmp(nw) = string(str(start:i-1))
          end if
       else if (inquote) then
          cycle   ! inside quotes: keep accumulating
       else if (str(i:i) .ne. ' ' .and. str(i:i) .ne. tab) then
          if (.not. inword) then
             inword = .true.; start = i
          end if
       else
          if (inword) then
             inword = .false.; nw = nw + 1; tmp(nw) = string(str(start:i-1))
          end if
       end if
    end do
    if (inword .and. .not. inquote) then
       nw = nw + 1; tmp(nw) = string(str(start:n))
    end if
    allocate(words(max(nw,1)))
    if (nw .gt. 0) words(1:nw) = tmp(1:nw)
  end subroutine split_words

  !===================================================================!
  ! Print the parsed configuration
  !===================================================================!

  subroutine print(this)
    class(config), intent(in) :: this
    integer :: i
    write(*,'(1x,a,a)')    "mesh     : ", this % meshfile % str
    write(*,'(1x,a,a)')    "equation : ", this % equation % str
    write(*,'(1x,a,es12.4,a,es12.4)') "kappa    : ", this % kappa, "   source : ", this % source
    write(*,'(1x,a,a,a,es12.4,a,i0)') "solver   : ", this % solver % str, &
         & "   max_tol : ", this % max_tol, "   max_it : ", this % max_it
    do i = 1, this % nbc
       write(*,'(1x,a,a,1x,a,3es12.4)') "boundary : ", trim(this % bc_name(i) % str), &
            & trim(this % bc_kind(i) % str), this % bc_a(i), this % bc_b(i), this % bc_c(i)
    end do
    write(*,'(1x,a,a)')    "output   : ", this % output % str
  end subroutine print

end module class_config
