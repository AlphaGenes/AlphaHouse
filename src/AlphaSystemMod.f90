
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaSystemMod.f90
!
! DESCRIPTION:
!> @brief    Module cotaining system procedures to be used to interface with an OS.
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 DWilson - Initial Version
!
!-------------------------------------------------------------------------------

module AlphaSystemMod

    contains

        integer function checkSystemCommand(command)
            use ifport
            
            character(len=*) :: command
            checkSystemCommand = SYSTEM(command)

        end function checkSystemCommand











end module AlphaSystemMod