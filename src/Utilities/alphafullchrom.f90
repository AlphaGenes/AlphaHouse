module alphaFullChromModule

    use AlphaImputeModule


      interface
     function runProgram (specFileInput, ped)
        type(specFilebase) :: specFileInput
        type(pedigreeHolder) :: ped
     end function runProgram

contains 


    subroutine fullChromSplit(pedigreeFile, specfile, funPointer)
        implicit none
        character(len=*) :: pedigreeFile
        integer :: nChroms, result

        character(len=300) :: chromPath

        procedure(runProgram), pointer, intent(in):: funPointer

        ! TODO currently only opperate on a single pedigree, can be chnaged to take in more using array
        ! Will use more memory
        ! This means this can be done in parallel, even mpi
        ped = PedigreeHolder(pedigreeFile)

        do i=1, nChroms

            specFile%resultFolderPath = chromPath
            write(chromPath,'(a,i0)') "chr",i
            result=makedirqq(prepend//trim(chromPath))
            if (result /= 0) then
                write(error_unit,*) "ERROR - creating directories has failed with error code:", result
            endif

            call ped%wipeGenotypeAndPhaseInfo


            call ped%addGenotypeInformationFromFile("genoFile")


            call funPointer(specFile,ped)
            ! cwd for program should be settable
            ! Can be done using CHDIR command
            
        enddo 

        

    end subroutine fullChromSplit




subroutine parseMapFile(filename)

    character(len=*) :: filename

    


end subroutine parseMapFile


end module alphaFullChromModule