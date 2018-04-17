#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"


#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)



#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""


#endif
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
		type(pedigreeHolder), target, optional :: ped
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
		use ifport

		type(pedigreeHolder) :: ped
		character(len=*) :: plinkPre
		procedure(runProgram), pointer, intent(in):: funPointer
		character(len=128), dimension(:), allocatable :: chromPaths
		integer :: i
		integer, dimension(:), allocatable :: nsnps
		integer :: sexChroms
		integer :: totalToDo, curChrom
		class(baseSpecFile) :: specfile
		real(kind=real64), dimension(:) ,allocatable :: lengths
		integer, dimension(:) ,allocatable :: basepairs

		call initialiseMPI


		if (specfile%plinkBinary) then
			call readPlink(plinkPre, ped, chromPaths,nsnps, sexChroms)
		else
			call readPlinkNoneBinary(plinkPre, ped, chromPaths,nsnps, sexChroms)
		endif

		totalToDo = (size(chromPaths))/mpiSize
		if (totalToDo <1) then
			if (mpiRank+1 < nCoreLengths) return
		endif

		do i=1, totalToDo

			if (allocated(specFile%useChroms)) then
				! Check if we are only doing a subset of chromsomes
				if (.not. any(specFile%useChroms == i)) cycle
			endif
			curChrom = (mpiRank+1)+((i-1) * size(chromPaths) )
			result=makedirqq("MultiChromResults")
			path = "MultiChromResults/" // curChrom
			result=makedirqq(path)
			CALL chdir(path)

			specFile%resultFolderPath = chromPaths(i)
			specFile%nsnp = nsnps(i)
			specFile%CurrChrom = i
			! write(chromPath,'(a,i0)') "chr",i
			! result=makedirqq(prepend//trim(chromPath))

			if (i > size(chromPaths)-2 .and. sexChroms /= 0) then
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

			call ped%addGenotypeInformationFromFile(chromPaths(i)//"genotypes.txt",nsnps(i),initAll=1)
			call ped%setSnpBasePairs(trim(chromPaths(i))//"snpBasepairs.txt",nsnps(i))
			call ped%setSnpLengths(trim(chromPaths(i))//"snplengths.txt",nsnps(i))

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
		type(plinkInfoType) :: plinkInfo
		class(baseSpecFile) :: specfile

		if (specfile%plinkBinary) then
			call readPlink(plinkPre, ped, chromPaths,plinkInfo,specFile%useChroms)
		else
			call readPlinkNoneBinary(plinkPre, ped, chromPaths,plinkInfo,specFile%useChroms)
		endif

		call ped%printPedigreeOriginalFormat("PLINKPED.txt")
		do i=1, size(chromPaths)

			if (allocated(specFile%useChroms)) then
				if (.not. any(specFile%useChroms == i)) cycle
			endif
			! Check if we are only doing a subset of chromsomes
			specFile%CurrChrom = i

			specFile%resultFolderPath = trim(chromPaths(i))//"results"
			specFile%nsnp = plinkInfo%nsnpsPerChromosome(i)
			print *,"doing chrom ", i
			if (i > size(chromPaths)-2 .and. plinkInfo%sexChrom /=0) then
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
			print *,"using ",trim(chromPaths(i))//"genotypes.txt",plinkInfo%nsnpsPerChromosome(i)

			call ped%addGenotypeInformationFromFile(trim(chromPaths(i))//"genotypes.txt",plinkInfo%nsnpsPerChromosome(i),initAll=1)
			call ped%setSnpBasePairs(trim(chromPaths(i))//"snpBasepairs.txt",plinkInfo%nsnpsPerChromosome(i))
			call ped%setSnpLengths(trim(chromPaths(i))//"snplengths.txt",plinkInfo%nsnpsPerChromosome(i))


			print *,"starting function run"
			if (.not. specFile%validate()) then
				write(error_unit, *) "ERROR - Spec file validation has failed"
			endif

			if (specfile%stopAfterPlink) then
				call addNeccessaryOutputForProgram(specfile,chromPaths, ped)
				print *, "Stopping program after plink output - neccessary files have been copied."
				call exit(0)
			endif
			call funPointer(specFile,ped)
			print *,"Finished function run"

		enddo

		if (specfile%plinkOutput) then
			call writePedFile(ped,plinkInfo,specfile,chromPaths)
			call writeMapFile(plinkInfo)
			call writeRefFile(plinkInfo)
		endif

	end subroutine runPlink


	subroutine addNeccessaryOutputForProgram(specFile,chromPaths, ped)

		use ifport
		use baseSpecFileModule
		use pedigreeModule
		use alphahousemod

		class(basespecfile), intent(in) :: specfile
		character(len=128), dimension(:), allocatable,intent(in) :: chromPaths
		type(PedigreeHolder), intent(in) :: ped
		class(basespecfile), allocatable :: specFileTemp
		integer :: status,i
		character(len=:),allocatable :: exePath

		call specFile%copy(specFileTemp)
		specFileTemp%plinkinputfile = ""
		specFileTemp%PlinkOutput = .false.
		specFileTemp%resultFolderPath = "Results"
		specFileTemp%PedigreeFile = "Pedigree.txt"

		do i=1, size(chromPaths)
			if (allocated(specFile%useChroms)) then
				! Check if we are only doing a subset of chromsomes
				if (.not. any(specFile%useChroms == i)) cycle
			endif
			print *, "PATH: ", trim(chromPaths(i))
			call specFileTemp%writeSpec(trim(chromPaths(i)) //trim(specFile%programName)//"Spec.txt")
			! copy program executable
			call getExecutablePath(exePath)
#ifdef _WIN32
			status = SYSTEMQQ("mklink " // trim(chrompaths(i))//trim(specFile%programName) // " " // trim(exePath))
#else
			status = SYSTEMQQ("ln -sf " //  trim(exePath)// " " // trim(chrompaths(i)) // "/.")
#endif
			! status = SYSTEMQQ(COPY // " " //  trim(exePath)// " " // chrompaths(i))
			call ped%printPedigreeOriginalFormat(trim(chrompaths(i))//"pedigree.txt")
		end do

	end subroutine addNeccessaryOutputForProgram

#endif
end module alphaFullChromModule






