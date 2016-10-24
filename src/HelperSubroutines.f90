


!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     HelperSubroutines.f90
!
! DESCRIPTION:
!> @brief   Simple subroutines and functions that will be used throughout the AlphaSuite
!
!
!> @author   David Wilson david.wilson@roslin.ed.ac.uk
!
!> @date     October 22, 2016
!
!> @version  0.0.1 (alpha)
!
!------------------------------------------------------------------------------
module HelperSubroutines

    contains

    function getIntegerSeedFromFile(seedFile) result(res)
        character(len=*),intent(in), optional :: seedFile
        integer :: res, unit

        if (present(seedFile)) then
            open(newunit=unit, file=seedFile, status="old")
        else
            open(newunit=unit, file="Seed.txt", status="old")
        endif

        read(unit, * ) res


        return

    end function getIntegerSeedFromFile


end module HelperSubroutines