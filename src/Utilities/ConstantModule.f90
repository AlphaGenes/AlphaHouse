!###############################################################################
!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     ConstantModule.f90
!
! DESCRIPTION:
!> @brief    Holds various constants
!
!> @details  Holds constants that are used in different programs.
!
!> @author   Various
!
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
!
!-------------------------------------------------------------------------------
module ConstantModule

    integer, parameter :: MissingPhaseCode = 9 !< value for missing phase values
    integer, parameter :: ErrorPhaseCode = -1 !< Error code for phase
    integer, parameter :: MissingGenotypeCode = 3 !< Value for missing parents
    character, parameter :: EMPTY_PARENT = '0' !< Value for unknown parents
    character, parameter :: EMPTYID = '0' !< value for unknown ids
    integer, parameter :: IDLENGTH = 32 
    integer, parameter :: IDINTLENGTH = 8
    integer, parameter :: generationThreshold = 1000
    integer, parameter :: OFFSPRINGTHRESHOLD = 150
    integer, parameter :: NOGENERATIONVALUE = -9999
    integer, parameter :: FILELENGTH = 256
    integer, parameter :: SPECOPTIONLENGTH = 300
    integer, parameter :: MISSINGGENDERCODE = -9
    integer , parameter ::  DEFAULTDICTSIZE = 4993 !< Default size for the hashtable
    integer, parameter :: DICT_NULL = -214748364 !<default value for when hash table doesn't match
end Module ConstantModule
