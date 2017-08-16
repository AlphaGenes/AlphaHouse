module alphaFullChrom

    use AlphaImputeModule


      interface
     function specFileInput (filename)
        type(specFilebase) :: specFileInput
        character(len=*), intent (in) :: file\
     end function func

contains 


    subroutine fullChromSplit(pedigreeFile, specfile, funPointer)
        implicit none
        character(len=*) :: pedigreeFile
        integer :: nChroms, result

        character(len=300) :: chromPath

        procedure(forced), pointer, intent(in):: funPointer

        ! TODO currently only opperate on a single pedigree, can be chnaged to take in more using array
        ! Will use more memory
        ! This means this can be done in parallel, even mpi
        ped = PedigreeHolder(pedigreeFile)

        do i=1, nChroms

            write(chromPath,'(a,i0)') "chr",i
            result=makedirqq(prepend//trim(chromPath))
            if (result /= 0) then
                write(error_unit,*) "ERROR - creating directories has failed with error code:", result
            endif

            call ped%wipeGenotypeAndPhaseInfo


            call ped%addGenotypeInformationFromFile("genoFile")


            call funPointer
            ! cwd for program should be settable
            ! Can be done using CHDIR command
            
        enddo 

        

    end subroutine fullChromSplit


end module alphaFullChrom