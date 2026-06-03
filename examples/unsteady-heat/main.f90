!=====================================================================!
! Transient heat conduction driver - the problem lives entirely in a
! config file, nothing is hardcoded here. Backward-euler marches the
! diffusion equation and writes a vtu snapshot every few steps so the
! evolving field can be animated.
!
!   ./run <config-file>
!
! The march is done with class_time_integrator (backward euler); calling
! it over successive windows is what lets each frame be written out.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program unsteady_heat

  use iso_fortran_env       , only : dp => real64
  use class_config          , only : config
  use class_gmsh_loader     , only : gmsh_loader
  use class_mesh            , only : mesh
  use class_assembler       , only : assembler
  use class_diffusion       , only : diffusion
  use class_time_integrator , only : time_integrator
  use class_paraview_writer , only : paraview_writer
  use class_string          , only : string
  use module_verbosity      , only : set_verbosity

  implicit none

  type(config)                        :: cfg
  character(len=256)                  :: cfgfile
  class(gmsh_loader)    , allocatable :: gmsh
  class(mesh)           , allocatable :: grid
  class(assembler)      , allocatable :: fvm
  type(time_integrator)               :: ti
  class(paraview_writer), allocatable :: pw
  real(dp)              , allocatable :: phi(:), phinew(:)
  type(string)                        :: labels(1)
  real(dp)                            :: t, out_dt
  integer                             :: i, frame, nframes, nout, nsteps

  ! Read the run
  call get_command_argument(1, cfgfile)
  if (len_trim(cfgfile) .eq. 0) then
     print *, "usage: run <config-file>"
     error stop
  end if

  cfg = config(trim(cfgfile))
  call cfg % print()
  call set_verbosity(cfg % verbosity)

  if (cfg % dt .le. 0.0_dp) then
     print *, "unsteady-heat is transient - set dt / t_init / t_final in the config"
     error stop
  end if

  ! Geometry
  allocate(gmsh, source = gmsh_loader(trim(cfg % meshfile % str)))
  allocate(grid, source = mesh(gmsh))

  ! Assembler + physics (set the equation before the boundaries)
  allocate(fvm, source = assembler(grid))
  call fvm % set_equation(diffusion(cfg % kappa, source = cfg % source))

  ! Boundary conditions, by physical group name
  do i = 1, cfg % nbc
     select case (trim(cfg % bc_kind(i) % str))
     case ("dirichlet")
        call fvm % set_dirichlet(trim(cfg % bc_name(i) % str), cfg % bc_c(i))
     case ("neumann")
        call fvm % set_neumann(trim(cfg % bc_name(i) % str), cfg % bc_c(i))
     case ("robin")
        call fvm % set_robin(trim(cfg % bc_name(i) % str), cfg % bc_a(i), cfg % bc_b(i), cfg % bc_c(i))
     end select
  end do

  ! Frame schedule - aim for ~10 vtu snapshots across the run
  nsteps  = nint((cfg % tfinal - cfg % tinit)/cfg % dt)
  nout    = max(1, nsteps/10)
  nframes = nsteps/nout
  out_dt  = cfg % dt*real(nout, dp)

  ! One integrator advances the state by nout steps per call
  ti = time_integrator(fvm, cfg % tinit, cfg % tinit + out_dt, cfg % dt, cfg % max_it, cfg % max_tol)

  ! Output
  labels(1) = string("phi")
  allocate(pw, source = paraview_writer(fvm % grid))

  ! Initial state phi = 0
  allocate(phi(fvm % num_state_vars))
  phi = 0.0_dp
  t   = cfg % tinit

  call pw % write(frame_name(cfg % output % str, 0), reshape(phi, [size(phi), 1]), labels)
  call report(0, t, phi)

  ! March, one window (nout steps) per frame
  do frame = 1, nframes
     call ti % integrate(phi, phinew)
     phi = phinew
     t   = t + out_dt
     call pw % write(frame_name(cfg % output % str, frame), reshape(phi, [size(phi), 1]), labels)
     call report(frame, t, phi)
  end do

  print *, "wrote", nframes + 1, "snapshots over t =", cfg % tinit, "..", t

contains

  !===================================================================!
  ! Snapshot file name "<stem>_NNNN.vtu" from an "<stem>.vtu" base
  !===================================================================!

  function frame_name(base, k) result(name)

    character(len=*), intent(in) :: base
    integer         , intent(in) :: k

    character(len=:), allocatable :: name
    character(len=4)              :: kk
    integer                       :: dotpos

    write(kk, '(i4.4)') k
    dotpos = index(base, ".vtu", back = .true.)

    if (dotpos .gt. 0) then
       name = base(1:dotpos-1)//"_"//kk//".vtu"
    else
       name = trim(base)//"_"//kk//".vtu"
    end if

  end function frame_name

  !===================================================================!
  ! Print min / max / mean of one frame
  !===================================================================!

  subroutine report(k, time, u)

    integer , intent(in) :: k
    real(dp), intent(in) :: time
    real(dp), intent(in) :: u(:)

    write(*,'(1x,a,i4,a,f9.2,a,3es13.5)') &
         & "frame ", k, "  t =", time, "  min/max/mean =", &
         & minval(u), maxval(u), sum(u)/size(u)

  end subroutine report

end program unsteady_heat
