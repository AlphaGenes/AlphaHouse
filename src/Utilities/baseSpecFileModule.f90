
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     baseSpecFileModule.f90
!
! DESCRIPTION:
!> @brief    Module cotaining abstract class for base spec files
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
module baseSpecFileModule
	use iso_fortran_env
	type, abstract :: baseSpecFile

	character(len=512) :: resultFolderPath !< Path where results should go
	character(len=512) :: plinkinputfile !< prepend to plink file
	logical :: plinkBinary
	integer(kind=int32) :: nsnp
	integer(kind=1) :: SexOpt,HetGameticStatus, HomGameticStatus

	end type




end module baseSpecFileModule