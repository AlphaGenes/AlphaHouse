
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     alphafullchrom.f90
!
! DESCRIPTION:
!> @brief    Module containing definition of pedigree.
!
!> @details Functions fro running alphaprograms designed to work on single chromosomes across multiple chromosomes
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2017-09-26 Dwilson - Initial Version

!-------------------------------------------------------------------------------
module alphaFullChromModule

	use baseSpecFileModule

	interface
	subroutine runProgram (specFileInput, ped)
		use baseSpecFileModule
		use pedigreeModule

		class(baseSpecFile),target :: specFileInput
		type(pedigreeHolder), optional :: ped
	end subroutine runProgram
end interface
contains

#ifdef MPIACTIVE

	!---------------------------------------------------------------------------
	!< @brief Runs plink based on how many chromosomes are on the file. Parallelizes this across clusters
	!< @details takes in a BASESPECFILE object, and then a function pointer for the main function of an alphaproram
	!< @author  David Wilson david.wilson@roslin.ed.ac.uk
	!---------------------------------------------------------------------------
	subroutine runPlink(plinkPre, specfile, funPointer)
		use CompatibilityModule
		use pedigreeModule
		use baseSpecFileModule
		use MPIUtilities

		type(pedigreeHolder) :: ped
		character(len=*) :: plinkPre
		procedure(runProgram), pointer, intent(in):: funPointer
		character(len=128), dimension(:), allocatable :: chromPaths
		integer :: i
		integer, dimension(:), allocatable :: nsnps
		logical :: sexChroms
		integer :: totalToDo, curChrom
		class(baseSpecFile) :: specfile

		call initialiseMPI


		if (specfile%plinkBinary) then
			call readPlink(plinkPre, ped, chromPaths,nsnps, sexChroms)
		else
			call readPlinkNoneBinary(plinkPre, ped, chromPaths,nsnps, sexChroms)
		endif

		totalToDo = (size(chromPaths))/mpiSize
		print *,"TOTALTODO", totalTOdo
		if (totalToDo <1) then
			if (mpiRank+1 < nCoreLengths) return
		endif

		do i=1, totalToDo
			curChrom = (mpiRank+1)+((i-1) * size(chromPaths) )


			specFile%resultFolderPath = chromPaths(i)
			specFile%nsnp = nsnps(i)
			! write(chromPath,'(a,i0)') "chr",i
			! result=makedirqq(prepend//trim(chromPath))

			if (i > size(chromPaths)-2 .and. sexChroms) then
				if (i == size(chromPaths)-1) then !< x chrom
					specFile%SexOpt = 1
					specFile%HetGameticStatus=1
					specFile%HomGameticStatus=2
				else !< y chrom
					specFile%SexOpt = 1
					specFile%HetGameticStatus=2
					specFile%HomGameticStatus=1
				endif
			endif

			call ped%wipeGenotypeAndPhaseInfo


			call ped%addGenotypeInformationFromFile(chromPaths(i)//"genotypes.txt",nsnps(i))


			call funPointer(specFile,ped)

		enddo
		call endMPI

	end subroutine runPlink

#else


	!---------------------------------------------------------------------------
	!< @brief Runs plink based on how many chromosomes are on the file.
	!< @details takes in a BASESPECFILE object, and then a function pointer for the main function of an alphaproram
	!< @author  David Wilson david.wilson@roslin.ed.ac.uk
	!---------------------------------------------------------------------------
	subroutine runPlink(plinkPre, specfile, funPointer)
		use CompatibilityModule
		use pedigreeModule
		use baseSpecFileModule

		type(pedigreeHolder) :: ped
		character(len=*) :: plinkPre
		procedure(runProgram), pointer, intent(in):: funPointer
		character(len=128), dimension(:), allocatable :: chromPaths
		integer :: i
		integer, dimension(:), allocatable :: nsnps
		logical :: sexChroms
		class(baseSpecFile) :: specfile

		if (specfile%plinkBinary) then
			call readPlink(plinkPre, ped, chromPaths,nsnps, sexChroms)
		else
			call readPlinkNoneBinary(plinkPre, ped, chromPaths,nsnps, sexChroms)
		endif

		call ped%printPedigreeOriginalFormat("PLINKPED.txt")
		do i=1, size(chromPaths)

			specFile%resultFolderPath = chromPaths(i)
			specFile%nsnp = nsnps(i)
			! write(chromPath,'(a,i0)') "chr",i
			! result=makedirqq(prepend//trim(chromPath))

			if (i > size(chromPaths)-2 .and. sexChroms) then
				if (i == size(chromPaths)-1) then !< x chrom
					specFile%SexOpt = 1
					specFile%HetGameticStatus=1
					specFile%HomGameticStatus=2
				else !< y chrom
					specFile%SexOpt = 1
					specFile%HetGameticStatus=2
					specFile%HomGameticStatus=1
				endif
			endif



			print *,"path:",chromPaths(i)
			call ped%wipeGenotypeAndPhaseInfo
					! first chrom should already be read in
			! if (i /= 1) then
			print *,"using ",trim(chromPaths(i))//"genotypes.txt"
				call ped%addGenotypeInformationFromFile(trim(chromPaths(i))//"genotypes.txt",nsnps(i),initAll=1)
			! endif


			print *,"starting function run"
			call funPointer(specFile,ped)
			print *,"Finished function run"
			


		enddo


	end subroutine runPlink
#endif
end module alphaFullChromModule


