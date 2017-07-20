module informationModule


    use AlphaStatMod
    use PedigreeModule
    use AlphaHouseMod

    contains
    function checkYield(ped) result(res)

        type(PedigreeHolder) :: ped
        integer(kind=1) ,dimension(:,:), allocatable  :: array
        real :: res
        array = ped%getGenotypesAsArrayWitHMissing()
        res = 1 - real(count(array == 9))/real(size(array))



    end function checkYield



    function calculateAccuracyPerAnimal(ped, trueFile, outputPerAnimalPath) result(meanAccuracy)
        
        type(PedigreeHolder) :: ped
        character(len=*),intent(in) :: trueFile
        character(len=*),intent(in),optional :: outputPerAnimalPath   

        integer(kind=1), allocatable, dimension(:) :: tmpArray
        integer :: animal,i, unit

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
        do i=1,lines

            read(unit,*) tmpIDs(i), tmpArray

            animal = ped%dictionary%getValue(tmpIDs(i))
            if (animal /= DICT_NULL) Then
                accuracies(i) = corAccuracies(tmpArray, ped%pedigree(animal)%individualGenotype%toIntegerArray())
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
        integer :: nsnps,unit

        ! function assumes order of input true file is the same as genotype file

        nsnps = ped%pedigree(ped%genotypeMap(1))%individualgenotype%length

        ! if nothing is genotyped, we can't do anything of worth here
        if (ped%nGenotyped <1) then 
            return
        endif

        ! set accuracies as 1, so not pedigreed animals won't affect
        accuracies = 1
        
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




    function corAccuracies(one, two) result(res)

        integer(kind=1), dimension(:),intent(in) :: one,two
        real(kind=real32) :: res

        integer :: incorrect
        incorrect = 0


        do i =1,size(one)

            if (one(i) == 9 .or. two(i) == 9) cycle

            if (one(i)/=two(i)) incorrect = incorrect + 1   
        enddo 

        res = 1 - real(incorrect) / real(size(one))
    end function corAccuracies
    ! subroutine checkPhaseYield

end module informationModule