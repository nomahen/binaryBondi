!!****if* source/Simulation/SimulationMain/hlaBoundary/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!!
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***

subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)

  use Driver_data, ONLY: dr_simTime, dr_globalMe, dr_globalComm
  use Simulation_data
  use Grid_interface
  use Eos_interface, ONLY : Eos_wrapped, Eos
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Gravity_data, only: grv_factor, grv_ptxpos, grv_ptypos, grv_ptzpos, grv_ptmass, grv_totmass

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Flash_mpi.h"

  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList
  integer, OPTIONAL, intent(IN):: pass

  real :: my_simTime
  real, dimension(EOS_NUM) :: eosData   !to recalc eos qtys inside sink

  integer :: error, iii, size

  integer :: block_no

  integer :: i,j,k,blk
  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC,eosRange,bcs

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  logical :: gcell = .true.

  integer :: istat, AllocateStatus
  real, pointer :: solnData(:,:,:,:)

  INTEGER :: xCoordSize, yCoordSize, zCoordSize
  REAL, ALLOCATABLE, DIMENSION(:) :: x, y, z

  !AA: Need this stuff for accretion
  real, dimension(MDIM) :: delta
  real, dimension(2) ::  xx, yy, zz
  real, dimension(2) ::  mnew_loc, mold_loc  !local accreted mass on each processor
  real :: dist(2)
  real ::  dx, dy, dz
  real ::  dvol, Pi, G

  call PhysicalConstants_get("pi",Pi)
  call PhysicalConstants_get("newton",G)

  yy(:) = 0.
  zz(:) = 0. 

  dy = 1.0
  dz = 1.0
!--------------------------------------------------

  !! AA: loop over all the blocks and apply the sink conditions

  mnew_loc(:) = 0.d0
  mold_loc(:) = 0.d0

  ! turn on pt mass over time
  if(dr_simTime .lt. sim_tRamp) then
        grv_factor(1) = - G * grv_ptmass(1) * dr_simTime / sim_tRamp
        grv_factor(2) = - G * grv_ptmass(2) * dr_simTime / sim_tRamp
  else
        grv_factor(1) = - G * grv_ptmass(1)
        grv_factor(2) = - G * grv_ptmass(2)
  endif

  do block_no = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(block_no), blkLimits, blkLimitsGC)
     call Grid_getDeltas(blockList(block_no),delta)

     xCoordSize = blkLimitsGC(HIGH, IAXIS)
     yCoordSize = blkLimitsGC(HIGH, JAXIS)
     zCoordSize = blkLimitsGC(HIGH, KAXIS)

     allocate(x(xCoordSize), stat=istat)
     allocate(y(yCoordSize), stat=istat)
     allocate(z(zCoordSize), stat=istat)

     ! get coordinates of current cell, store in x, y, z
     call Grid_getCellCoords(IAXIS, blockList(block_no), CENTER, .TRUE., x, xCoordSize)
     call Grid_getCellCoords(JAXIS, blockList(block_no), CENTER, .TRUE., y, yCoordSize)
     call Grid_getCellCoords(KAXIS, blockList(block_no), CENTER, .TRUE., z, zCoordSize)

     ! if current block is a leaf
     call Grid_getBlkPtr(blockList(block_no), solnData)

     !! AA: loop over all the cells in the current block
     ! Loop z
     do k=blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        if (NDIM > 2) then
           zz(1) = z(k) - grv_ptzpos(1)
           zz(2) = z(k) - grv_ptzpos(2)
           dz = delta(3)
        end if

        ! loop y
        do j=blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
            if (NDIM > 1) then
               yy(1) = y(j) - grv_ptypos(1)
               yy(2) = y(j) - grv_ptypos(2)
               dy = delta(2)
            end if

           ! loop x
           do i=blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              dx = delta(1)

              xx(1) = x(i) - grv_ptxpos(1)
              xx(2) = x(i) - grv_ptxpos(2)

              dist(1) = sqrt(xx(1)**2. + yy(1)**2. + zz(1)**2.)  !dist from cell to sink
              dist(2) = sqrt(xx(2)**2. + yy(2)**2. + zz(2)**2.)  !dist from cell to sink
              dvol = dx*dy*dz                        !cell volume
              
              ! NO: accretion loop (scale accretion radius to mass) (primary)
              if (abs(dist(1)) .le. (sim_radius*grv_ptmass(1)/grv_totmass) ) then !if cell in sink, apply bc
               !add old and new cell masses to running totals in mold_loc and mnew_loc
               mnew_loc(1) = mnew_loc(1) + solnData(DENS_VAR,i,j,k)*dvol

	       !flor DENS and PRES in the soluData array
               solnData(DENS_VAR,i,j,k) = sim_rhoAmbient*sim_presfactor
               solnData(PRES_VAR,i,j,k) = sim_pAmbient*sim_presfactor
               solnData(VELX_VAR,i,j,k) = 0.0
               solnData(VELY_VAR,i,j,k) = 0.0
               solnData(VELZ_VAR,i,j,k) = 0.0

               mold_loc(1) = mold_loc(1) + solnData(DENS_VAR,i,j,k)*dvol

	       !pass new values to eos unit
               eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
               eosData(EOS_PRES) = solnData(PRES_VAR,i,j,k)

	       !ask eos to recalc all gas vars using new PRES and DENS
               call Eos(MODE_DENS_PRES,1,eosData)

               !update solnData for those blocks with the new gas values
               solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
               solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)

              ! NO: accretion loop (scale accretion radius to mass) (secondary)
              elseif (abs(dist(2)) .le. (sim_radius*grv_ptmass(2)/grv_totmass) ) then !if cell in sink, apply bc

               !add old and new cell masses to running totals in mold_loc and mnew_loc
               mnew_loc(2) = mnew_loc(2) + solnData(DENS_VAR,i,j,k)*dvol

	       !floor DENS and PRES in the soluData array
               solnData(DENS_VAR,i,j,k) = sim_rhoAmbient*sim_presfactor
               solnData(PRES_VAR,i,j,k) = sim_pAmbient*sim_presfactor
               solnData(VELX_VAR,i,j,k) = 0.0
               solnData(VELY_VAR,i,j,k) = 0.0
               solnData(VELZ_VAR,i,j,k) = 0.0

               mold_loc(2) = mold_loc(2) + solnData(DENS_VAR,i,j,k)*dvol

	       !pass new values to eos unit
               eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
               eosData(EOS_PRES) = solnData(PRES_VAR,i,j,k)

	       !ask eos to recalc all gas vars using new PRES and DENS
               call Eos(MODE_DENS_PRES,1,eosData)

               !update solnData for those blocks with the new gas values
               solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
               solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
              endif ! NO: end accreiton loop

              ! NO: boundary loop (distance evaluated from origin)
              if (abs(sqrt(x(i)**2.0 + y(j)**2.0 + z(k)**2.0)) .ge. sim_rbound) then
               ! NO: set density, pressure to ambient values; velocity to zero
               solnData(DENS_VAR,i,j,k) = sim_rhoAmbient
               solnData(PRES_VAR,i,j,k) = sim_pAmbient
               ! NKO: Input vx,vy based on angular momentum profiles
               select case (sim_amProf)
                 case ("none")
                   solnData(VELX_VAR,i,j,k) = 0.
                   solnData(VELY_VAR,i,j,k) = 0.
                   solnData(VELZ_VAR,i,j,k) = 0.
                 case ("rigid")
                   ! NKO: The min statements are for the radii > sim_rbound. If
                   ! rigidly rotating, these would be super-Keplerian; min statement
                   ! limits the rotation profiles so they're only ever Keplerian
                   ! This should not affect the solution because the BC is reset
                   ! at sim_rbound every timestep. 
                   solnData(VELX_VAR,i,j,k) = -y(j)*sqrt(G*sim_totmass/sim_rbound**3.0)*sim_kepFact
                   solnData(VELX_VAR,i,j,k)  = min(solnData(VELX_VAR,i,j,k),-y(j)*&
                                               sqrt(G*sim_totmass/(x(i)*x(i)+y(j)*y(j)+z(k)*z(k))**1.5))
                   solnData(VELY_VAR,i,j,k) = x(i)*sqrt(G*sim_totmass/sim_rbound**3.0)*sim_kepFact
                   solnData(VELY_VAR,i,j,k)  = min(solnData(VELY_VAR,i,j,k),x(i)*&
                                               sqrt(G*sim_totmass/(x(i)*x(i)+y(j)*y(j)+z(k)*z(k))**1.5)) 
                   solnData(VELZ_VAR,i,j,k) = 0.
                 case ("shell")
                   solnData(VELX_VAR,i,j,k) = -y(j)*sqrt(G*sim_totmass/(x(i)*x(i)+y(j)*y(j)+z(k)*z(k))**1.5)*sim_kepFact
                   solnData(VELY_VAR,i,j,k) = x(i)*sqrt(G*sim_totmass/(x(i)*x(i)+y(j)*y(j)+z(k)*z(k))**1.5)*sim_kepFact
                   solnData(VELZ_VAR,i,j,k) = 0.
                 case ("corotating")
                   solnData(VELX_VAR,i,j,k) = -1.0*y(j)*sim_w
                   solnData(VELY_VAR,i,j,k) = x(i)*sim_w
                   solnData(VELZ_VAR,i,j,k) =  0.0
               end select
 
               ! corrections here for rotating frame
               solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + y(j)*sim_w
               solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) - x(i)*sim_w


	       !pass new values to eos unit
               eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
               eosData(EOS_PRES) = solnData(PRES_VAR,i,j,k)

	       !ask eos to recalc all gas vars using new PRES and DENS
               call Eos(MODE_DENS_PRES,1,eosData)

               !update solnData for those blocks with the new gas values
               solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
               solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT)
               solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
              endif ! NO: end boundary loop

           end do !i
        end do !j
     end do !k

     call Grid_releaseBlkPtr(blockList(block_no), solnData)

     call Eos_wrapped(MODE_DENS_EI, blkLimitsGC, blockList(block_no))
     deallocate(x)
     deallocate(y)
     deallocate(z)
  end do !block_no


  !! sum mnew_loc and mold_loc across all procs. put the totals in sink_mass_new and sink_mass_old
  CALL MPI_ALLREDUCE(mnew_loc(1), sink_mass_new(1), 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, error)
  CALL MPI_ALLREDUCE(mold_loc(1), sink_mass_old(1), 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, error)
  CALL MPI_ALLREDUCE(mnew_loc(2), sink_mass_new(2), 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, error)
  CALL MPI_ALLREDUCE(mold_loc(2), sink_mass_old(2), 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, error)

  CALL MPI_ALLREDUCE(dr_simTime, my_simTime, 1, FLASH_REAL, MPI_MIN, MPI_COMM_WORLD, error)

  !! Calculate instantaneous accretion rate
  mdot(1) = (sink_mass_new(1) - sink_mass_old(1)) / dt
  mdot(2) = (sink_mass_new(2) - sink_mass_old(2)) / dt

  !! Open the data file and write the calculated mdot
  if (dr_globalMe .EQ. MASTER_PE) then
     open(30, file='mdot.dat', position='APPEND')
     write(30,'(F40.20)',ADVANCE="No") my_simTime
     write(30,'(A)',ADVANCE="No") '    '
     write(30,'(F40.20)',ADVANCE="No") mdot(1)/Pi !! BH Mdot units
     write(30,'(A)',ADVANCE="No") '    '
     write(30,'(F40.20)',ADVANCE="No") mdot(2)/Pi !! BH Mdot units
     write(30,*) ! Next line
     close (30)
  end if

  return
end subroutine Driver_sourceTerms
