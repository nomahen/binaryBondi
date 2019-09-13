!!****if* source/Simulation/SimulationMain/hlaBoundary/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_data, ONLY : gr_globalMe 
  implicit none

#include "constants.h"
#include "Flash.h"

  ! local variables
  integer :: i,j,k, io, AllocateStat
  real    :: soundspeed, G, pi

  ! gas  params
  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('sim_rbound', sim_rbound)

  ! sink params
  call RuntimeParameters_get('sim_radius', sim_radius)
  call RuntimeParameters_get('sim_max_ref_rad', sim_max_ref_rad)
  call RuntimeParameters_get('sim_presfactor', sim_presfactor)
  call RuntimeParameters_get("sim_tRamp",sim_tRamp)

  ! grf params
  call RuntimeParameters_get("sim_useGrf", sim_useGrf)
  call RuntimeParameters_get("sim_grf", sim_grf)
  call RuntimeParameters_get("sim_numGrf", sim_numGrf)

  ! amom profile params
  call RuntimeParameters_get("sim_amProf", sim_amProf)
  call RuntimeParameters_get("sim_kepFact", sim_kepFact)

  ! Copies of variables from gravity Config, to be used in Simulation_initBlock
  call RuntimeParameters_get("totmass", sim_totmass)
  call RuntimeParameters_get("sma", sim_sma)

  call PhysicalConstants_get("newton", G)
  call PhysicalConstants_get("pi", pi)

  ! sim_totmass is normalized to bondi units
  sim_totmass = 1.0/G
  sim_period = 2.0*pi*sqrt((sim_sma**3.0)/(G*sim_totmass))
  sim_w      = 2.0*pi/sim_period

  ! NO: mdot to 0 to start
  mdot(:) = 0.0
  sink_mass_new(:) = 0.0
  sink_mass_old(:) = 0.0

  ! NO: set sound speed to unity
  soundspeed = 1.0

  ! NO: Set pressure accordingly
  sim_pambient = soundspeed*soundspeed*sim_rhoAmbient/sim_gamma


  !! AA: make a file to store accreted mass data
  if (gr_globalMe ==  MASTER_PE) then
  
     open (30, file="mdot.dat", iostat=io,  &
          status='unknown', position='append')
     write(30,*)                &
     'time                       ', &
     'mdot1                      ', &
     'mdot2                      '
     close(30)
  endif

  ! NO: Parse our (g)aussian (r)andom (f)ield data
  if (sim_useGrf) then
    allocate ( sim_dataGrf(sim_numGrf,sim_numGrf,sim_numGrf))

    open(2, file = sim_grf, status = 'old')

    ! fill grid values from GRF generated in python
    do i=1,sim_numGrf
      do j=1,sim_numGrf
        do k=1,sim_numGrf
        read(2,*) sim_dataGrf(i,j,k)
        enddo
      enddo
    enddo

  
    close(2)
  endif  

  !! Making sure I got the units right!
  if (gr_globalMe  == MASTER_PE) then
  
  
     print*,"######## Checking Simulation parameters  ##################"
     print*,"sim_rhoAmbient = ", sim_rhoAmbient
     print*,"sim_pAmbient = ",sim_pAmbient
     print*,"sim_gamma = ",sim_gamma
     print*,"Sound Speed = ", soundspeed
     print*,"Using GRF = ", sim_useGrf
     print*,"GRF filename = ", sim_grf
     print*,"GRF grid points per dim = ", sim_numGrf
     print*,"##############################################"

  endif

end subroutine Simulation_init
