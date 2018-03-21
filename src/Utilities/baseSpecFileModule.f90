
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     baseSpecFileModule.f90
!
! DESCRIPTION:
!> @brief    Module cotaining abstract class for base spec files
!> Designed to work with programs that work on a single chromosome
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
	logical :: plinkBinary !< are the plink files binary
	integer(kind=int32) :: nsnp !< number of snp for this chromosme
	integer(kind=1) :: SexOpt,HetGameticStatus, HomGameticStatus
	logical :: plinkOutput !< if true - output in plink format
	logical :: stopAfterPlink !< if true -- will call exit after plink output has been created - 
	integer(kind=int32), dimension(:), allocatable :: useChroms !< Array containing chromosomes to do
	contains
		! procedure :: validateBase
		procedure :: validate => validateBase
	end type baseSpecFile

	contains

		function validateBase(params) result(res)

			class(baseSpecFile) :: params
			LOGICAL :: res
			res = .true.


		end function validateBase


		subroutine writeOutSpecOptions(params)
		class(baseSpecFile), intent(in) :: params
		integer :: unit,i
			character(len=:), allocatable :: tmpString
		end subroutine writeOutSpecOptions

end module baseSpecFileModule
