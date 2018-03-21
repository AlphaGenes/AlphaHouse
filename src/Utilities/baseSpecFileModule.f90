
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
	type :: baseSpecFile

	character(len=512) :: resultFolderPath !< Path where results should go
	character(len=512) :: plinkinputfile !< prepend to plink file
	character(len=300):: PedigreeFile = "NoPedigree",GenotypeFile="Genotypes.txt" !< Pedigree and genotype file - again from the plink input
	logical :: plinkBinary !< are the plink files binary
	integer(kind=int32) :: nsnp !< number of snp for this chromosme
	integer(kind=1) :: SexOpt,HetGameticStatus, HomGameticStatus
	logical :: plinkOutput !< if true - output in plink format
	logical :: stopAfterPlink = .false. !< if true -- will call exit after plink output has been created -
	integer(kind=int32), dimension(:), allocatable :: useChroms !< Array containing chromosomes to do
	character(len=20) :: programName = "UNSPECFIED"
	character(len=32) :: version = "UNSPECFIED"
	contains
		! procedure :: validateBase
		procedure :: validate => validateBase
		procedure :: writeSpec => writeOutBase
		procedure :: copy => copyBase
	end type baseSpecFile

	contains


		subroutine copyBase(old, new)

			class(baseSpecFile),intent(in) ::  old
			class(baseSpecFile),allocatable, intent(out) ::  new


			allocate(baseSpecFile::new)
			new%resultFolderPath = old%resultFolderPath
			new%plinkinputfile = old%plinkinputfile
			new%PedigreeFile = old%PedigreeFile
			new%GenotypeFile = old%GenotypeFile
			new%plinkBinary = old%plinkBinary
			new%nsnp = old%nsnp
			new%SexOpt = old%SexOpt
			new%HetGameticStatus = old%HetGameticStatus
			new%HomGameticStatus = old%HomGameticStatus
			new%plinkOutput = old%plinkOutput
			new%stopAfterPlink = old%stopAfterPlink
			new%useChroms = old%useChroms
			new%programName = old%programName
			new%version = old%version

		end subroutine copyBase
		function validateBase(params) result(res)

			class(baseSpecFile) :: params
			LOGICAL :: res
			res = .true.


		end function validateBase


		subroutine writeOutBase(params,path)
			class(baseSpecFile), intent(in) :: params
			integer :: unit,i
			character(len=*), optional, intent(in):: path
			character(len=:), allocatable :: tmpString


			if (present(path)) then
				open(newunit=unit, file=path,status='unknown')
			else
				open(newunit=unit, file=trim(params%programName)//"Spec.txt",status='unknown')
			endif

			write(unit,* ) "ERROR - writing of spec files not implemented for this program, please contact developers"
		end subroutine writeOutBase

end module baseSpecFileModule


