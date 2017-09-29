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

		

		! subroutine fullChromSplit(pedigreeFile, specfile, funPointer)
		! 	use pedigreeModule
			
		! 	implicit none
		! 	character(len=*) :: pedigreeFile
		! 	integer :: nChroms, result
		! 	type(specFilebase) :: specfile
		! 	type(pedigreeHolder) :: ped
		! 	character(len=300) :: chromPath

		! 	procedure(runProgram), pointer, intent(in):: funPointer

		! 	! TODO currently only opperate on a single pedigree, can be chnaged to take in more using array
		! 	! Will use more memory
		! 	! This means this can be done in parallel, even mpi
		! 	ped = PedigreeHolder(pedigreeFile)

		! 	do i=1, nChroms

		! 		specFile%resultFolderPath = chromPath
		! 		write(chromPath,'(a,i0)') "chr",i
		! 		result=makedirqq(prepend//trim(chromPath))
		! 		if (result /= 0) then
		! 			write(error_unit,*) "ERROR - creating directories has failed with error code:", result
		! 		endif

		! 		call ped%wipeGenotypeAndPhaseInfo


		! 		call ped%addGenotypeInformationFromFile("genoFile")


		! 		call funPointer(specFile,ped)
		! 		! cwd for program should be settable
		! 		! Can be done using CHDIR command

		! 	enddo



		! end subroutine fullChromSplit


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

				call ped%wipeGenotypeAndPhaseInfo


				call ped%addGenotypeInformationFromFile(chromPaths(i)//"genotypes.txt",nsnps(i))


				call funPointer(specFile,ped)
				

				


				! cwd for program should be settable
				! Can be done using CHDIR command

			enddo


		end subroutine runPlink

end module alphaFullChromModule
