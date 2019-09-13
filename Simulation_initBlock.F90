!!****if* source/Simulation/SimulationMain/hlaBoundary/Simulation_initBlock
!!
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId,
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Grid_data, ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) ::  blockId


  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  xDist, yDist, zDist
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek
  real     ::  dist
  logical  ::  validGeom
  integer  :: istat

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData

  logical :: gcell = .true.

  real :: G

  call PhysicalConstants_get("newton", G)

  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1)
        else
           dzz = zCoord(k) - zCoord(k-1)
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)

     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1)
           else
              dyy = yCoord(j) - yCoord(j-1)
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1)
           else
              dxx = xCoord(i) - xCoord(i-1)
           endif
           
           ! NKO: If using gaussian random field in initial conditions, input it
           ! & here.
           if (sim_useGrf) then 
             ! NKO: Monster if loop to check if the coordinates are outside the
             ! & computational domain. If so, this is a guard cell of an exterior
             ! & block, and our GRF wont interpolate here; this wont influence and
             ! & flow and can be set to the ambient density. 
             if ((xx .lt. gr_imin) .or. (xx .gt. gr_imax) .or. (yy .lt. gr_jmin) .or. &
                 (yy .gt. gr_jmax) .or. (zz .lt. gr_kmin) .or. (zz .gt. gr_kmax)) then
               rho = sim_rhoAmbient
             else
               call trilinear_interpolation(xx,yy,zz,gr_imin,gr_imax,gr_jmin,gr_jmax,gr_jmin,gr_jmax,sim_numGrf,sim_dataGrf,rho)
               rho = rho*sim_rhoAmbient
             endif
           else
             rho = sim_rhoAmbient
           endif
           p   = sim_pAmbient
           
           ! NKO: Input vx,vy based on angular momentum profiles
           select case (sim_amProf)
             case ("none")
               vx  = 0.
               vy  = 0.
               vz  = 0.
             case ("rigid")
               ! NKO: The min statements are for the radii > sim_rbound. If
               ! rigidly rotating, these would be super-Keplerian; min statement
               ! limits the rotation profiles so they're only ever Keplerian
               ! This should not affect the solution because the BC is reset
               ! at sim_rbound every timestep. 
               vx  = -yy*sqrt(G*sim_totmass/sim_rbound**3.0)*sim_kepFact
               vx  = min(vx,-yy*sqrt(G*sim_totmass/(xx*xx+yy*yy+zz*zz)**1.5))
               vy  = xx*sqrt(G*sim_totmass/sim_rbound**3.0)*sim_kepFact
               vy  = min(vy,xx*sqrt(G*sim_totmass/(xx*xx+yy*yy+zz*zz)**1.5)) 
               vz  = 0.
             case ("shell")
               vx  = -yy*sqrt(G*sim_totmass/(xx*xx+yy*yy+zz*zz)**1.5)*sim_kepFact
               vy  = xx*sqrt(G*sim_totmass/(xx*xx+yy*yy+zz*zz)**1.5)*sim_kepFact
               vz  = 0.
             case ("corotating")
               vx  = -1.0*yy*sim_w
               vy  = xx*sim_w
               vz  = 0.0
           end select

           ! make corrections for rotating reference frame
           ! w is in positive z direction
           vx = vx + yy*sim_w
           vy = vy - xx*sim_w

           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !  assume gamma-law equation of state
           e   = p/(sim_gamma-1.)
           e   = e/rho + ek
           e   = max (e, sim_smallP)

           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
              solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
           end if
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz

#else
           if (NSPECIES > 0) then
              ! putting in the value of the default species
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)

              !if there is only one species, this loop will not execute
              do n = SPECIES_BEGIN+1, SPECIES_END

                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)

#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock


!******************************************************************************

!  Routine:     trilinear_interpolation()

!  Description: Given coordinates x,y,z and domain boundaries xmin,xmax, ymin,
!               ymax, zmin,zmax, total number of grid points per side N, grid
!               of Gaussian Random Field values in dimension (N,N,N), return
!               trilinearly interpolated GRF value p. 

subroutine trilinear_interpolation(x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,N,grid,p)
  implicit none
  
  real, intent(IN) :: x,y,z
  real, intent(IN) :: xmin,ymin,zmin
  real, intent(IN) :: xmax,ymax,zmax
  integer, intent(IN) :: N
  real, dimension(N,N,N), intent(IN) :: grid
  real, intent(out) :: p

  real :: dx,dy,dz
  real, dimension(N) :: xgrid, ygrid, zgrid
  integer :: i, xx, yy, zz
  real :: x0,y0,z0,x1,y1,z1,xd,yd,zd
  real, dimension(2,2,2) :: c
  real, dimension(2,2) :: b
  real, dimension(2) :: a
  !!
  dx = (xmax - xmin) / (N-1)
  dy = (ymax - ymin) / (N-1)
  dz = (zmax - zmin) / (N-1)

  do i = 1,N
    xgrid(i) = xmin + dx*(i-1)
    ygrid(i) = ymin + dy*(i-1)
    zgrid(i) = zmin + dz*(i-1)
  end do

  call find(xgrid, N, x, xx)
  call find(ygrid, N, y, yy)
  call find(zgrid, N, z, zz)

  x0 = xgrid(xx)
  x1 = xgrid(xx+1)
  xd = (x-x0)/(x1-x0)
  y0 = ygrid(yy)
  y1 = ygrid(yy+1)
  yd = (y-y0)/(y1-y0)
  z0 = zgrid(zz)
  z1 = zgrid(zz+1)
  zd = (z-z0)/(z1-z0)

  c(1,1,1) = grid(xx,yy,zz)
  c(1,1,2) = grid(xx,yy,zz+1)
  c(1,2,1) = grid(xx,yy+1,zz)
  c(1,2,2) = grid(xx,yy+1,zz+1)
  c(2,1,1) = grid(xx+1,yy,zz)
  c(2,1,2) = grid(xx+1,yy,zz+1)
  c(2,2,1) = grid(xx+1,yy+1,zz)
  c(2,2,2) = grid(xx+1,yy+1,zz+1)

  b(1,1) = c(1,1,1)*(1.0-xd) + c(2,1,1)*xd
  b(1,2) = c(1,1,2)*(1.0-xd) + c(2,1,2)*xd
  b(2,1) = c(1,2,1)*(1.0-xd) + c(2,2,1)*xd
  b(2,2) = c(1,2,2)*(1.0-xd) + c(2,2,2)*xd

  a(1) = b(1,1)*(1.0-yd) + b(2,1)*yd
  a(2) = b(1,2)*(1.0-yd) + b(2,2)*yd

  p = a(1)*(1.0-zd) + a(2)*zd

  return
end subroutine trilinear_interpolation


!******************************************************************************

!  Routine:     find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine find(x, N, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: N
  integer, intent(OUT):: i
  real, intent(IN)    :: x(N), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(N)) then

     i = N

  else

     il = 1
     ir = N
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif
  return

end subroutine find
