!!****if* source/Simulation/SimulationMain/hlaCluster/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation
!!
!! PARAMETERS
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  ! boundary
  real, save    :: sim_pAmbient, sim_rhoAmbient, sim_gamma
  real, save    :: sim_smallX, sim_smallRho, sim_smallP, sim_mass
  real, save    :: sim_rbound

  ! sinks
  real, save    :: sim_max_ref_rad, sim_presfactor
  real, save    :: sim_radius, sim_tRamp

  ! binary parameters (copies of runtime parameters set in gravity Config)
  ! must copy because Gravity_init gets called after the grid is initialized,
  ! and we need these variables in Simulation_initBlock
  real, save    :: sim_totmass, sim_sma, sim_period, sim_w

  ! integral variables
  real, dimension(2), save    :: sink_mass_old, sink_mass_new, mdot

  ! general
  integer, save :: sim_meshMe

  ! for reading gaussian random field
  logical, save                       :: sim_useGrf
  character(len=40), save             :: sim_grf
  integer, save                       :: sim_numGrf
  real, allocatable, dimension(:,:,:) :: sim_dataGrf

  ! angular momentum profiles
  character(len=40), save :: sim_amProf
  real, save              :: sim_kepFact

end module Simulation_data
