module informationModule


    use AlphaStatMod
    use PedigreeModule
    use AlphaHouseMod

    implicit none
    contains
    function checkYield(ped) result(res)

        type(PedigreeHolder) :: ped
        integer(kind=1) ,dimension(:,:), allocatable  :: array
        real :: res
        array = ped%getGenotypesAsArrayWitHMissing()
        res = 1 - real(count(array == 9))/real(size(array))



    end function checkYield



    function calculateAccuracyPerAnimal(ped, trueFile, outputPerAnimalPath, snpErrorPath) result(meanAccuracy)
        
        type(PedigreeHolder) :: ped
        character(len=*),intent(in) :: trueFile
        character(len=*),intent(in),optional :: outputPerAnimalPath, snpErrorPath

        integer(kind=1), allocatable, dimension(:) :: tmpArray
        integer :: animal,i, unit, lines, snpErrorUnit

        character(len=IDLENGTH), dimension(:), allocatable :: tmpIDs 
        real,dimension(:), allocatable :: accuracies
        real :: meanAccuracy
        ! if nothing is genotyped, we can't do anything of worth here
        if (ped%nGenotyped <1) then 
            return
        endif

        ! set accuracies as 1, so not pedigreed animals won't affect
        accuracies = 1
        
        lines = countLines(trueFile)
        allocate(accuracies(lines))
        allocate(tmpIds(lines))
        ! allocate the
        allocate(tmpArray(ped%pedigree(ped%genotypeMap(1))%individualgenotype%length))
        open(newunit=unit, file=trueFile, status="old")

        if (present(snpErrorPath)) then
            open(newunit=snpErrorUnit, file=trim(snpErrorPath), status='unknown')
        endif
        do i=1,lines

            read(unit,*) tmpIDs(i), tmpArray

            animal = ped%dictionary%getValue(tmpIDs(i))
            if (animal /= DICT_NULL) Then
                ! print *,ped%pedigree(animal)%individualGenotype%getgenotype(232)
                if (present(snpErrorPath)) then
                    accuracies(i) = corAccuracies(tmpArray, ped%pedigree(animal)%individualGenotype%toIntegerArray(), i, snpErrorUnit)
                else
                    accuracies(i) = corAccuracies(tmpArray, ped%pedigree(animal)%individualGenotype%toIntegerArray(), i)
                endif
            endif
        enddo
        close(unit)

        if (present(outputPerAnimalPath)) then
            open(newunit=unit, file=outputPerAnimalPath, status="unknown" ) 
            do i=1,lines
                write(unit,*) tmpIds(i),accuracies(i)
            enddo
            close(unit)
        endif

        if (present(snpErrorPath)) then
            close(snpErrorUnit)
        endif
        meanAccuracy = 0
        meanAccuracy = sum(accuracies) / lines
            

    end function calculateAccuracyPerAnimal




    function calculateAccuracyPerSnp(ped, trueFile, outperSnpPath) result(meanAccuracy)
        
        type(PedigreeHolder) :: ped
        character(len=*),intent(in) :: trueFile
        character(len=*),intent(in), optional :: outperSnpPath 

        integer(kind=1), allocatable, dimension(:,:) :: tmpArray, pedArray
        integer :: animal,i

        character(len=IDLENGTH), dimension(:), allocatable :: tmpIDs 
        real,dimension(:), allocatable :: accuracies
        real :: meanAccuracy
        integer :: nsnps,unit,lines

        ! function assumes order of input true file is the same as genotype file

        nsnps = ped%pedigree(ped%genotypeMap(1))%individualgenotype%length

        ! if nothing is genotyped, we can't do anything of worth here
        if (ped%nGenotyped <1) then 
            return
        endif

        ! set accuracies as 1, so not pedigreed animals won't affect
        
        
        lines = countLines(trueFile)
        
        allocate(tmpIds(lines))
        ! allocate the
        allocate(tmpArray(lines,nsnps))
        allocate(pedArray(lines,nsnps))

        open(newunit=unit, file=trueFile, status="old")
        do i=1,lines
            read(unit,*) tmpIDs(i), tmpArray(i,:)


            animal = ped%dictionary%getValue(tmpIds(i))

            if (animal /= DICT_NULL) then
                pedArray(i,:) = ped%pedigree(animal)%individualGenotype%toIntegerArray() 
            endif
        enddo
        close(unit)

        allocate(accuracies(nsnps))
        accuracies = 1
        do i = 1, nsnps
        ! compare snps for genotyped animals
            accuracies(i) = corAccuracies(pedArray(:,i), tmpArray(:,i))
        enddo 

        if (present(outperSnpPath)) then
            open(newunit=unit, file=outperSnpPath, status="unknown" ) 
            do i=1,nsnps
                write(unit,*) i, accuracies(i)
            enddo
            close(unit)
        endif


        meanAccuracy = sum(accuracies) / nsnps
            

    end function calculateAccuracyPerSnp


    function correlation(one, two) result(res)
        real(kind=real64), dimension(:), intent(in) :: one, two
        real(kind=real64) :: res

        res = (sum(one*two)-sum(one)*sum(two))/ (sqrt(sum(one**2) - sum(one)**2) * sqrt(sum(two**2)-sum(two)**2))


    end function

    function corAccuracies(one, two, index, unit) result(res)

        integer(kind=1), dimension(:),intent(in) :: one,two
        real(kind=real64) :: res
        integer, optional :: index
        integer, optional :: unit

        integer :: incorrect, count, i
        incorrect = 0
        count = 0

    
        do i =1,size(one)

            if (one(i) == 9 .or. two(i) == 9) cycle
            count = count +1
            if (one(i)/=two(i)) then
                if (present(unit)) then
                    write(unit,*) "error at animal at line: ",index," at position ",i," with vals:",one(i),two(i)
                endif

                incorrect = incorrect + 1   
            endif
        enddo 

        res = 1 - real(incorrect) / real(count)

    end function corAccuracies
    ! subroutine checkPhaseYield

end module informationModule