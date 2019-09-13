!!****if* source/Grid/GridMain/paramesh/gr_unmarkRefineByLogRadius
!!
!! NAME
!!  gr_unmarkRefineByLogRadius
!!
!!  
!! SYNOPSIS 
!!  call gr_unmarkRefineByLogRadius(real(in) :: xc, 
!!                                  real(in) :: yc,
!!                                  real(in) :: zc)
!!  
!! DESCRIPTION
!!  Cancel refinement flags for all blocks that are too far away from a
!!  center given by (xc,yc,zc).
!!  The determination whether a block is 'too far away' depends on the
!!  current block size as well as the distance, and is made using the
!!  runtime parameter gr_lrefineMaxRedRadiusFact.
!!  
!! ARGUMENTS 
!!  xc -   Center of the interval/circle/sphere : IAXIS
!!  yc -                                          JAXIS
!!  zc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  
!! SIDE EFFECTS
!!
!!  Elements in the PARAMESH logical array refine(:) may be modified.
!!
!! NOTES
!! 
!!  This routine has not been tested well. It probably should be viewed only as a
!!  guideline for a user's implementation.
!!  
!!  If the geometry is SPHERICAL or POLAR, the distance is measured in the radial
!!  direction (X-direction) alone, and is taken as the radial distance from a
!!  sphere of radius given by xc.  In particular, the distance is the distance
!!  from the coordinate center if xc = 0.0.
!!
!!  3D cylindrical geometry is not supported.
!!***

subroutine gr_unmarkRefineByLogRadius(xc, yc, zc)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry, gr_lrefineMaxRedRadiusSq
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  real, dimension(2), intent(IN)      :: xc, yc, zc

! Local data

  real, dimension(MDIM) :: blockCenter
  real                  :: blockSmallestSide
  real,dimension(2)     :: dist2, radiusSq
  integer               :: b
  real :: bxl, bxr 

  if((NDIM == 3).and.(gr_geometry == CYLINDRICAL)) then
     call Driver_abortFlash("gr_unmarkRefineByLogRadius: 3D Cylindrical not supported")
  end if

  if((gr_geometry == CARTESIAN).or.(gr_geometry == CYLINDRICAL)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSmallestSide = minval(bsize(1:NDIM,b))
           radiusSq(1) = (BlockCenter(1) - xc(1))**2
           radiusSq(2) = (BlockCenter(1) - xc(2))**2
           if (NDIM > 1) then
              radiusSq(1) = radiusSq(1) + (BlockCenter(2) - yc(1))**2
              radiusSq(2) = radiusSq(2) + (BlockCenter(2) - yc(2))**2
           endif
           if ((NDIM == 3).and.(gr_geometry==CARTESIAN)) then
              radiusSq(1) = radiusSq(1) + (BlockCenter(3) - zc(1))**2
              radiusSq(2) = radiusSq(2) + (BlockCenter(3) - zc(2))**2
           endif

           ! Now compare the ratio of block's smallest side to distance from center
           ! to a threshold.  If the ratio is less, then cancel a pending refinement.
           if (all(blockSmallestSide**2 < gr_lrefineMaxRedRadiusSq * radiusSq)) then
              refine(b) = .false.
           end if
        end if
     end do

  else
     call Driver_abortFlash("MarkRefine: geometry spec is wrong")
     !-------------------------------------------------------------------------------
  end if
  return
end subroutine gr_unmarkRefineByLogRadius
