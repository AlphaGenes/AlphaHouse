module alphaFullChrom


contains 
    subroutine fullChromSplit(prepend,args)
        implicit none
        real :: args
        integer :: nChroms, result

        character(len=300) :: chromPath



        ! TODO currently only opperate on a single pedigree, can be chnaged to take in more using array
        ! Will use more memory
        ped = PedigreeHolder( )

        do i=1, nChroms

            write(chromPath,'(a,i0)') "chr",i
            result=makedirqq(prepend//trim(chromPath))
            if (result /= 0) then
                write(error_unit,*) "ERROR - creating directories has failed with error code:", result
            endif

            call ped%wipeGenotypeAndPhaseInfo


            call ped%addGenotypeInformationFromFile("genoFile")

            ! cwd for program should be settable
            ! Can be done using CHDIR command
            
        enddo 

        

    end subroutine fullChromSplit


end module alphaFullChrom