module alphaFullChrom


contains 
    subroutine fullChromSplit(prepend,args)
        implicit none
        real :: args
        integer :: nChroms, result

        character(len=300) :: chromPath



        do i=1, nChroms

            write(chromPath,'(a,i0)') "chr",i
            result=makedirqq(prepend//trim(chromPath))

            if (result /= 0) then
                write(error_unit,*) "ERROR - creating directories has failed with error code:", result
        enddo 

        

    end subroutine fullChromSplit


end module alphaFullChrom