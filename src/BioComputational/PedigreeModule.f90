
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     PedigreeModule.f90
!
! DESCRIPTION:
!> @brief    Module containing logic of Pedigree
!> @details  Module contains functions to read in and set up pedgiree data structure
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     January 4, 2017
!
!> @version  1.0.0
!
!
!-------------------------------------------------------------------------------




module PedigreeModule
    use IndividualModule
    use IndividualLinkedListModule
 use HashModule
 use constantModule
 use AlphaHouseMod, only : Int2Char

 private addOffspringsAfterReadIn

type PedigreeHolder

    type(Individual), pointer, dimension(:) :: Pedigree !have to use pointer here as otherwise won't let me point to it
    type(IndividualLinkedList) :: Founders !linked List holding all founders
    type(IndividualLinkedList),allocatable, dimension(:) :: generations !linked List holding each generation
    type(DictStructure) :: dictionary
    integer(kind=int32) :: pedigreeSize, nDummys !pedigree size cannot be bigger than 2 billion animals
    integer(kind=int32) :: maxPedigreeSize
    integer(kind=int32), dimension(:), allocatable :: sortedIndexList

    integer, dimension(:) , allocatable :: genotypeMap ! map going from genotypeMap(1:nAnisG) = newID 

    integer :: maxGeneration
    contains
        procedure :: destroyPedigree
        procedure :: setPedigreeGenerationsAndBuildArrays
        procedure :: outputSortedPedigree
        procedure :: setOffspringGeneration
        procedure :: addGenotypeInformationFromFile
        procedure :: addGenotypeInformationFromArray
        generic :: addGenotypeInformation => addGenotypeInformationFromArray, addGenotypeInformationFromFile
        procedure :: outputSortedPedigreeInAlphaImputeFormat
        procedure :: isDummy
        procedure :: sortPedigreeAndOverwrite
        procedure :: sortPedigreeAndOverwriteWithDummyAtTheTop
        procedure :: makeRecodedPedigreeArray
        procedure :: printPedigree
        procedure :: getMatePairsAndOffspring
        procedure :: getAllGenotypesAtPosition

end type PedigreeHolder

type RecodedPedigreeArray
  integer(kind=int32) :: nInd
  character(len=IDLENGTH), allocatable, dimension(:) :: originalId
  integer(kind=int32), allocatable, dimension(:) :: generation
  integer(kind=int32), allocatable, dimension(:,:) :: id
  contains
      procedure :: init    => initRecodedPedigreeArray
      procedure :: destroy => destroyRecodedPedigreeArray
      procedure :: write   => writeRecodedPedigreeArray
end type

    interface PedigreeHolder
        module procedure initPedigree
        module procedure initPedigreeArrays
        module procedure initEmptyPedigree
        module procedure initPedigreeGenotypeFiles
    end interface PedigreeHolder

    interface Sort !Sorts into generation list
        module procedure :: setPedigreeGenerationsAndBuildArrays
    end interface Sort
contains


    function initEmptyPedigree() result(pedStructure)
        use iso_fortran_env
        type(PedigreeHolder) :: pedStructure

        pedStructure%dictionary = DictStructure()
        pedStructure%pedigreeSize = 0
        pedStructure%nDummys = 0
        pedStructure%maxPedigreeSize = DEFAULTDICTSIZE
    end function initEmptyPedigree

    !---------------------------------------------------------------------------
    !< @brief Constructor for pedigree class
    !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    function initPedigree(fileIn, numberInFile, genderFile) result(pedStructure)
        use AlphaHouseMod, only : countLines
        use iso_fortran_env
        type(PedigreeHolder) :: pedStructure
        character(len=*),intent(in) :: fileIn !< path of pedigree file
        character(len=*), intent(in),optional :: genderFile !< path to gender file
        integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
        
        character(len=IDLENGTH) :: tmpId,tmpSire,tmpDam
        integer(kind=int32) :: stat, fileUnit,tmpSireNum, tmpDamNum, tmpGender,tmpIdNum
        integer(kind=int64) :: nIndividuals
        integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
        integer :: tmpAnimalArrayCount
        integer :: i
        integer(kind=int64) :: sizeDict
        logical :: sireFound, damFound

        pedStructure%nDummys = 0
        tmpAnimalArrayCount = 0
        
        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        sizeDict = nIndividuals
        pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
        allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
        pedStructure%pedigreeSize = nIndividuals
        pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
        pedStructure%maxGeneration = 0
        open(newUnit=fileUnit, file=fileIn, status="old")

        do i=1,nIndividuals

            sireFound = .false.
            damFound = .false.

            read(fileUnit,*) tmpId,tmpSire,tmpDam
            call pedStructure%dictionary%addKey(tmpId, i)

            pedStructure%Pedigree(i) =  Individual(trim(tmpId),trim(tmpSire),trim(tmpDam), i) !Make a new individual based on info from ped

            if (tmpSire /= EMPTY_PARENT) then !check sire is defined in pedigree
                tmpSireNum = pedStructure%dictionary%getValue(tmpSire)
                if (tmpSireNum /= DICT_NULL) then
                    sireFound = .true.
                endif
            endif

            if (tmpDam /= EMPTY_PARENT) then
                tmpDamNum = pedStructure%dictionary%getValue(tmpDam)
                if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
                    damFound = .true.
                endif
            endif
            if (tmpSire == EMPTY_PARENT .and. tmpDam == EMPTY_PARENT) then !if animal is a founder
                pedStructure%Pedigree(i)%founder = .true.
                call pedStructure%Founders%list_add(pedStructure%Pedigree(i))

            else if (sireFound == .false. .or. damFound == .false. ) then
                tmpAnimalArrayCount = tmpAnimalArrayCount +1
                tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
            else ! if sire and dam are both found
                pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))
                call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male

                pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))
                call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
            endif
        enddo

        close(fileUnit)



        ! if we want gender info read in rather than calculated on the fly, lets do it here
        if (present(genderFile)) then !read in gender here
            open(newUnit=fileUnit, file=genderFile, status="old")
            do
                read (fileUnit,*, IOSTAT=stat) tmpId,tmpGender
                if (stat /=0) exit

                tmpIdNum = pedStructure%dictionary%getValue(tmpId)
                if (tmpIdNum /= DICT_NULL) then
                    pedStructure%Pedigree(i)%gender = int(tmpGender)
                else
                    write(error_unit, *) "ERROR: Gender  defined for an animal that does not exist in Pedigree!"
                    write(error_unit, *) "Amimal:",tmpId
                endif
            end do


        endif


    call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)

    deallocate(tmpAnimalArray)
    write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys
    end function initPedigree




      !---------------------------------------------------------------------------
    !< @brief Constructor for pedigree class using Genotype File format
    !< @details Constructor builds pedigree, without any sorting being done. 
    !< If no pedigree file is supplied, all animals are founders
    !< If an animal is in the pedigree, but not in the genotypeFile, this animal is not created (no dummy either) and the animals parents are set to empty.
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    function initPedigreeGenotypeFiles(fileIn, numberInFile, nSnp, pedFile) result(pedStructure)
        use AlphaHouseMod, only : countLines
        use iso_fortran_env
        type(PedigreeHolder) :: pedStructure
        integer, intent(in) :: nSnp !< number of snps to read
        character(len=*),intent(in) :: fileIn !< path of Genotype file
        integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
        character(len=*),intent(in), optional :: pedFile !< path of pedigreeFile file

        character(len=IDLENGTH) :: tmpId
        integer(kind=int32) :: stat, fileUnit,tmpIdNum
        integer(kind=int64) :: nIndividuals
        integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
        integer :: tmpAnimalArrayCount
        integer :: i
        integer, dimension(:), allocatable :: tmpGeno
        integer(kind=int64) :: sizeDict


        allocate(tmpGeno(nsnp))
        pedStructure%nDummys = 0
        tmpAnimalArrayCount = 0
        
        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        sizeDict = nIndividuals
        pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
        allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
        pedStructure%pedigreeSize = nIndividuals
        pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
        pedStructure%maxGeneration = 0
        open(newUnit=fileUnit, file=fileIn, status="old")

        do i=1,nIndividuals

            read(fileUnit,*) tmpId,tmpGeno(:)
            call pedStructure%dictionary%addKey(tmpId, i)
            pedStructure%Pedigree(i) =  Individual(trim(tmpId),"0","0", i) !Make a new individual based on info from ped
        enddo

        close(fileUnit)



        ! if we want gender info read in rather than calculated on the fly, lets do it here
        if (present(pedFile)) then !read in gender here

            block
                character(len=IDLENGTH) :: tmpSire, tmpDam
                integer :: tmpSireNum, tmpDamNum
                
                open(newUnit=fileUnit, file=pedFile, status="old")
                do
                    read (fileUnit,*, IOSTAT=stat) tmpId,tmpSire,tmpDam
                    if (stat /=0) exit

                    
                        
                    tmpIdNum = pedStructure%dictionary%getValue(tmpId)
                    if (tmpIdNum /= DICT_NULL) then
                        if (trim(tmpSire) == "0" .and. trim(tmpDam) == "0") then
                            call pedStructure%Founders%list_add(pedStructure%pedigree(tmpIdNum))
                        else
                            tmpSireNum = pedStructure%dictionary%getValue(tmpSire)
                            tmpDamNum = pedStructure%dictionary%getValue(tmpDam)
                            if (tmpSireNum /= DICT_NULL) then
                                pedStructure%pedigree(tmpIdNum)%sirePointer => pedStructure%pedigree(tmpSireNum)
                                call pedStructure%pedigree(tmpSireNum)%addOffspring(pedStructure%pedigree(tmpIdNum))
                            ! else
                                ! TODO ask daniel what to do about dummys
                            endif
                            if (tmpDamNum /= DICT_NULL) then
                                pedStructure%pedigree(tmpIdNum)%damPointer => pedStructure%pedigree(tmpDamNum)
                                call pedStructure%pedigree(tmpDamNum)%addOffspring(pedStructure%pedigree(tmpIdNum))
                            ! else
                            !     ! TODO ask daniel what to do about dummys
                            endif
                        endif
                    else
                        ! We only care about animals that are in genotype file
                        continue
                    endif
                end do
            
            end block
        else
        ! Otherwise, with no animals, every animal is a founder 
            do i=1, pedStructure%pedigreeSize
               call pedStructure%Founders%list_add(pedStructure%pedigree(i)) 
            enddo
        endif

    end function initPedigreeGenotypeFiles

    !---------------------------------------------------------------------------
    !< @brief Constructor for pedigree class
    !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals.
    !< Takes Two arrays as inputs rather than files
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    function initPedigreeArrays(pedArray, genderArray) result(pedStructure)
        use iso_fortran_env
        type(PedigreeHolder) :: pedStructure

        character(len=IDLENGTH), dimension(:,:), intent(in) :: pedArray !< array detailing pedigree of format ped([id, sireId, damId], index)
        integer, dimension(:) ,optional , intent(in):: genderArray !< gender array corresponding to index in pedArray
        integer(kind=int32) :: tmpSireNum, tmpDamNum
        integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
        integer :: tmpAnimalArrayCount
        integer :: i
        integer(kind=int64) :: sizeDict
        logical :: sireFound, damFound

        pedStructure%nDummys = 0
        tmpAnimalArrayCount = 0

        sizeDict = size(pedArray)
        pedStructure%maxPedigreeSize = size(pedArray) + (size(pedArray) * 4)
        allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
        pedStructure%pedigreeSize = size(pedArray)
        pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(size(pedArray))) !allocate to nIndividuals in case all animals are in incorrect order of generations
        pedStructure%maxGeneration = 0

        do i=1,size(pedArray)

            sireFound = .false.
            damFound = .false.

            call pedStructure%dictionary%addKey(pedArray(1,i), i)

            pedStructure%Pedigree(i) =  Individual(pedArray(1,i),pedArray(2,i),pedArray(3,i), i) !Make a new individual based on info from ped

            if (pedArray(i,2) /= EMPTY_PARENT) then !check sire is defined in pedigree
                tmpSireNum = pedStructure%dictionary%getValue(pedArray(2,i))
                if (tmpSireNum /= DICT_NULL) then
                    sireFound = .true.
                endif
            endif

            if (pedArray(i,3) /= EMPTY_PARENT) then
                tmpDamNum = pedStructure%dictionary%getValue(pedArray(3,i))
                if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
                    damFound = .true.
                endif
            endif
            if (pedArray(2,i) == EMPTY_PARENT .and. pedArray(3,i) == EMPTY_PARENT) then !if animal is a founder
                pedStructure%Pedigree(i)%founder = .true.
                call pedStructure%Founders%list_add(pedStructure%Pedigree(i))

            else if (sireFound == .false. .or. damFound == .false. ) then
                tmpAnimalArrayCount = tmpAnimalArrayCount +1
                tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
            else ! if sire and dam are both found
                pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))
                call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male

                pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))
                call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
            endif
        enddo


        if (present(genderArray)) then

            do i=1, size(genderArray)

                if (genderArray(i) /=MISSINGGENDERCODE) then
                    pedStructure%pedigree(i)%gender = genderArray(i)
                endif
            enddo
        endif

    call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)

    deallocate(tmpAnimalArray)
    write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys
    end function initPedigreeArrays



    subroutine addOffspringsAfterReadIn(pedStructure, tmpAnimalArray, tmpAnimalArrayCount)
        use ConstantModule, only : IDLENGTH
        class(PedigreeHolder) :: pedStructure
        integer, dimension(:), intent(in) :: tmpAnimalArray !< array containing indexes of tmp animals
        integer, intent(in) :: tmpAnimalArrayCount !< number of animals actually in tmpAnimalArray
        logical :: sireFound, damFound


        integer(kind=int32) :: tmpSireNum, tmpDamNum
        integer(kind=int32) :: i, tmpCounter
        character(len=IDLENGTH) :: tmpSire,tmpDam,tmpCounterStr
        
        tmpCounter = 0
                !check animals that didn't have parental information initially
        ! this is done to avoid duplication when a pedigree is sorted
        do i=1,tmpAnimalArrayCount
            sireFound = .false.
            damFound = .false.
            tmpSire = pedStructure%Pedigree(tmpAnimalArray(i))%getSireId()
            tmpSireNum = pedStructure%dictionary%getValue(tmpSire)
            tmpDam = pedStructure%Pedigree(tmpAnimalArray(i))%getDamId()
            tmpDamNum = pedStructure%dictionary%getValue(tmpDam)

            if (tmpSire /= EMPTY_PARENT) then
                if (tmpSireNum /= DICT_NULL) then !if sire has been found in hashtable
                    pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                    call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                    call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male

                else !if sire is defined but not in the pedigree, create him
                    ! check if the tmp animal has already been created
                    tmpSireNum = pedStructure%dictionary%getValue("dum"//trim(tmpSire))
                    if (tmpSireNum == DICT_NULL) then
                        pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                        pedStructure%nDummys = pedStructure%nDummys + 1
                        if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                            write(error_unit,*) "ERROR: too many undefined animals"
                            stop
                            ! TODO do a move alloc here to avoid this hack
                        endif
                        pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual("dum"//trim(tmpSire),'0','0', pedStructure%pedigreeSize)
                        call pedStructure%dictionary%addKey("dum"//trim(tmpSire), pedStructure%pedigreeSize)
                        pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                        pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                        call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                        call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                        pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                    else
                        pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                        call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                        call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                    endif
                endif
                sireFound = .true.
            endif

            if (tmpDam /= EMPTY_PARENT) then
                if (tmpDamNum /= DICT_NULL) then !if dam has been found
                        pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                        call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                else
                    ! Check for defined animals that have nit been set in pedigree
                    tmpDamNum = pedStructure%dictionary%getValue("dum"//trim(tmpDam))
                    if (tmpDamNum == DICT_NULL) then !If dummy animal has not already been set in pedigree
                        pedStructure%nDummys = pedStructure%nDummys + 1
                        pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                        if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                            write(error_unit,*) "ERROR: too many undefined animals"
                            stop
                        endif
                        pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual("dum"//trim(tmpDam),'0','0', pedStructure%pedigreeSize)
                        call pedStructure%dictionary%addKey("dum"//trim(tmpDam), pedStructure%pedigreeSize)
                        pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                        pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                        call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))

                        call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                        pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                    else
                        pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                        call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a sire, it should be male, dam female
                    endif

                endif
                damFound = .true.
            endif


            if (.not. damFound .XOR. .not. sireFound) then
                tmpCounter =  tmpCounter + 1
                write(tmpCounterStr, '(I3.3)') tmpCounter
                pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                pedStructure%nDummys = pedStructure%nDummys + 1
                    if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                        write(error_unit,*) "ERROR: too many undefined animals"
                        stop
                    endif
                    pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual("dum"//trim(tmpCounterStr),'0','0', pedStructure%pedigreeSize)
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                    if (.not. damFound) then
                        pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                    else if (.not. sireFound) then
                        pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                    endif
                    call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                    call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
            endif
        enddo

    end subroutine addOffspringsAfterReadIn



    !---------------------------------------------------------------------------
    !< @brief distructor for pedigree class
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine destroyPedigree(this)
        class(PedigreeHolder) :: this
        integer :: i

        do i=1,this%pedigreeSize
            call this%Pedigree(i)%destroyIndividual
            if (allocated(this%generations)) then
                deallocate(this%generations)
            endif
        enddo
        call this%Founders%destroyLinkedList
        call this%dictionary%destroy !destroy dictionary as we no longer need it
            deallocate(this%pedigree)

    end subroutine destroyPedigree



    subroutine addGenotypeInformationFromArray(this, array)

        use AlphaHouseMod, only : countLines
        implicit none
        class(PedigreeHolder) :: this
        integer(kind=1),allocatable,dimension (:,:), intent(in) :: array !< array should be dimensions nanimals, nsnp
        integer :: i 
        do i=1,this%pedigreeSize-this%nDummys
            call this%pedigree(i)%setGenotypeArray(array(i,:))
        enddo

    end subroutine addGenotypeInformationFromArray

    subroutine addGenotypeInformationFromFile(this, genotypeFile, nsnps, nAnnisG)

        use AlphaHouseMod, only : countLines
        implicit none
        class(PedigreeHolder) :: this
        character(len=*) :: genotypeFile
        character(len=IDLENGTH) :: tmpID
        integer,intent(in) :: nsnps
        integer,intent(in),optional :: nAnnisG
        integer(kind=1), allocatable, dimension(:) :: tmpSnpArray
        integer :: i, j,fileUnit, nAnnis,tmpIdNum


            if (present(nAnnisG)) then
                nAnnis = nAnnisG
            else
                nAnnis = countLines(genotypeFile)
            endif

            allocate(tmpSnpArray(nsnps))
            open(newUnit=fileUnit, file=genotypeFile, status="old")
            do i=1, nAnnis
                read (fileUnit,*) tmpId,tmpSnpArray(:)
                do j=1,nsnps
                  if ((tmpSnpArray(j)<0).or.(tmpSnpArray(j)>2)) tmpSnpArray(j)=9
                enddo
                tmpIdNum = this%dictionary%getValue(tmpId)
                if (tmpIdNum == DICT_NULL) then
                    write(error_unit, *) "ERROR: Genotype info for non existing animal"
                else
                    call this%pedigree(tmpIdNum)%setGenotypeArray(tmpSnpArray)
                endif
            enddo



    end subroutine addGenotypeInformationFromFile

    !---------------------------------------------------------------------------
    !< @brief builds correct generation information by looking at founders
    !< This is effectively a sort function for the pedigree
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine setPedigreeGenerationsAndBuildArrays(this)

        implicit none
        class(PedigreeHolder) :: this
        integer :: i
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        tmpIndNode => this%Founders%first
        allocate(this%generations(0:generationThreshold))
        do i=1, this%Founders%length
            call this%setOffspringGeneration(0, tmpIndNode%item)

            tmpIndNode => tmpIndNode%next
        end do

    end subroutine setPedigreeGenerationsAndBuildArrays


    !---------------------------------------------------------------------------
    !< @brief returns true if individual at given index is a isDummy
    !< if 0 is given, return false
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !< @param[in] file path (string)
    !---------------------------------------------------------------------------
    logical function isDummy(this, id)
        implicit none
        class(PedigreeHolder) :: this
        integer, intent(in) :: id !< ID to check if animal is dummy

        if (id == 0) then
            isDummy = .false.
        else if (this%pedigree(id)%isDummy) then
            isDummy = .true.
        else
            isDummy = .false.
        endif
    end function isDummy


    !---------------------------------------------------------------------------
    !< @brief writes sorted pedigree information to either a file or stdout
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !< @param[in] file path (string)
    !---------------------------------------------------------------------------
    subroutine outputSortedPedigree(this,file)

        use iso_fortran_env, only : output_unit
        class(PedigreeHolder) :: this
        character(len=*), intent(in), optional :: file !< output path for sorted pedigree
        integer :: unit, i,h,sortCounter
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        character(len=:), allocatable :: fmt

        sortCounter = 0
         if (.not. allocated(this%sortedIndexList)) then
            allocate(this%sortedIndexList(this%pedigreeSize))
        endif

        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays
        endif

        if (present(file)) then
            open(newUnit=unit, file=file, status="unknown")
        else
            unit = output_unit
        endif

        fmt = "(3a"//Int2Char(IDLENGTH)//", i"//Int2Char(IDINTLENGTH)//")"
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                sortCounter = sortCounter + 1
                    this%sortedIndexList(sortCounter) = tmpIndNode%item%id
                    write(unit, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId, tmpIndNode%item%generation
                    ! write(*, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
                    tmpIndNode => tmpIndNode%next
            end do
        enddo
        if (present(file)) then
            close(unit)
        endif
    end subroutine outputSortedPedigree

    subroutine outputSortedPedigreeInAlphaImputeFormat(this, file)
        use iso_fortran_env, only : output_unit
        class(PedigreeHolder) :: this
        character(len=*), intent(in), optional :: file !< output path for sorted pedigree
        integer :: unit, i,h, sortCounter
        type(IndividualLinkedListNode), pointer :: tmpIndNode


        if (.not. allocated(this%sortedIndexList)) then
            allocate(this%sortedIndexList(this%pedigreeSize))
        endif
        sortCounter = 0
        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays
        endif
        if (present(file)) then
            open(newUnit=unit, file=file, status="unknown")
        else
            unit = output_unit
        endif

        block
            integer :: sireId, damId
            character(len=:), allocatable :: fmt
            fmt = "(3i"//Int2Char(IDINTLENGTH)//", a"//Int2Char(IDLENGTH)//")"
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                    if (associated(tmpIndNode%item%damPointer)) then
                        damId = tmpIndNode%item%damPointer%id
                    else
                        damId = 0
                    endif
                    if (associated(tmpIndNode%item%sirePointer)) then
                        sireId = tmpIndNode%item%sirePointer%id
                    else
                        sireId = 0
                    endif
                    sortCounter = sortCounter +1
                    this%sortedIndexList(sortCounter) = tmpIndNode%item%id
                     write (unit, fmt) tmpIndNode%item%id,sireId,damId, tmpIndNode%item%originalID
                    ! write(*,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
                    tmpIndNode => tmpIndNode%next
            end do
        enddo
        endblock
        if (present(file)) then !avoids closing stdout
            close(unit)
        endif
    end subroutine outputSortedPedigreeInAlphaImputeFormat



    subroutine sortPedigreeAndOverwrite(this)
        use iso_fortran_env, only : output_unit, int64
        class(PedigreeHolder) :: this
        integer :: i,h, pedCounter, tmpId
        integer(kind=int64) :: sizeDict
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        type(Individual), pointer, dimension(:) :: newPed
        type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList
        type(IndividualLinkedList) :: dummyList
        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays
        endif
        pedCounter = 0
        call this%dictionary%destroy()
        call this%founders%destroyLinkedList()
        sizeDict  =this%pedigreeSize
        allocate(newPed(this%maxPedigreeSize))
        this%dictionary = DictStructure(sizeDict)
        allocate(newGenerationList(0:this%maxGeneration))
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                 if (tmpIndNode%item%isDummy) then
                    call dummyList%list_add(tmpIndNode%item)
                    cycle
                 endif
                 pedCounter = pedCounter +1
                 call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)
                 newPed(pedCounter) = tmpIndNode%item
                 newPed(pedCounter)%id = pedCounter
                 call newPed(pedCounter)%resetOffspringInformation ! reset offsprings
                 if (associated(newPed(pedCounter)%sirePointer)) then
                    tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
                    if (tmpID /= DICT_NULL) then
                        call newPed(tmpId)%addOffspring(newPed(pedCounter))
                        newPed(pedCounter)%sirePointer=> newPed(tmpId)
                    endif
                endif
                if (associated(newPed(pedCounter)%damPointer)) then
                    tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
                    if (tmpID /= DICT_NULL) then
                        call newPed(tmpId)%addOffspring(newPed(pedCounter))
                        newPed(pedCounter)%damPointer=> newPed(tmpId)
                    endif
                endif
                if (i ==0 ) then !if object is afounder add to founder array
                   call  this%founders%list_add(newPed(pedCounter))
                endif

                call newGenerationList(i)%list_add(newPed(pedCounter))
                tmpIndNode => tmpIndNode%next
            end do
        enddo
        ! call move_alloc(newPed, this%pedigree)
        ! print *, "size:", size(this%pedigree),":",size(newPed)

        tmpIndNode => dummyList%first

        ! add dummys to end of pedigree
        do i=1, dummyList%length
            pedCounter = pedCounter +1
            call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)
            newPed(pedCounter) = tmpIndNode%item
            newPed(pedCounter)%id = pedCounter
            call newPed(pedCounter)%resetOffspringInformation ! reset offsprings

            do h=1, tmpIndNode%item%nOffs
                tmpId =  this%dictionary%getValue(tmpIndNode%item%offsprings(h)%p%originalID)
                call newPed(pedCounter)%addOffspring(newPed(tmpID))
                if(trim(tmpIndNode%item%offsprings(h)%p%sireId) == trim(newPed(pedCounter)%originalID)) then
                    newPed(tmpId)%sirePointer => newPed(pedCounter)
                else
                    newPed(tmpId)%damPointer => newPed(pedCounter)
                endif
            enddo
            call  this%founders%list_add(newPed(pedCounter))

            call newGenerationList(0)%list_add(newPed(pedCounter))
            tmpIndNode => tmpIndNode%next

        enddo


        this%pedigree = newPed
        this%generations = newGenerationList
    end subroutine sortPedigreeAndOverwrite


        !---------------------------------------------------------------------------
    !< @brief Output pedigree to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine sortPedigreeAndOverwriteWithDummyAtTheTop(this)
        use iso_fortran_env, only : output_unit, int64
        class(PedigreeHolder) :: this
        integer :: i,h, pedCounter, tmpId
        integer(kind=int64) :: sizeDict
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        type(Individual), pointer, dimension(:) :: newPed
        type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList
        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays
        endif
        pedCounter = 0
        call this%dictionary%destroy()
        call this%founders%destroyLinkedList()
        sizeDict  =this%pedigreeSize
        allocate(newPed(this%maxPedigreeSize))
        this%dictionary = DictStructure(sizeDict)
        allocate(newGenerationList(0:this%maxGeneration))
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                 pedCounter = pedCounter +1
                 call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)
                 newPed(pedCounter) = tmpIndNode%item
                 newPed(pedCounter)%id = pedCounter
                 call newPed(pedCounter)%resetOffspringInformation ! reset offsprings
                 if (associated(newPed(pedCounter)%sirePointer)) then
                    tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
                    if (tmpID /= DICT_NULL) then
                        call newPed(tmpId)%addOffspring(newPed(pedCounter))
                        newPed(pedCounter)%sirePointer=> newPed(tmpId)
                    endif
                endif
                if (associated(newPed(pedCounter)%damPointer)) then
                    tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
                    if (tmpID /= DICT_NULL) then
                        call newPed(tmpId)%addOffspring(newPed(pedCounter))
                        newPed(pedCounter)%damPointer=> newPed(tmpId)
                    endif
                endif
                if (i ==0 ) then !if object is afounder add to founder array
                   call  this%founders%list_add(newPed(pedCounter))
                endif

                call newGenerationList(i)%list_add(newPed(pedCounter))
                tmpIndNode => tmpIndNode%next
            end do
        enddo


        this%pedigree = newPed
        this%generations = newGenerationList
    end subroutine sortPedigreeAndOverwriteWithDummyAtTheTop

    !---------------------------------------------------------------------------
    !< @brief Output pedigree to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine printPedigree(this)
          class(PedigreeHolder) :: this
          integer ::i
          do i= 1, this%pedigreeSize
            print *, this%pedigree(i)%id, this%pedigree(i)%getIntegerVectorOfRecodedIds()
          enddo
    end subroutine printPedigree

    !---------------------------------------------------------------------------
    !< @brief Sets generation of an individual and his children recursively
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !< @param[in] generation (integer)
    !< @param[in] pointer to an individual
    !---------------------------------------------------------------------------
    recursive subroutine setOffspringGeneration(this,generation, indiv)
        type(Individual),pointer, intent(inout) :: indiv
        integer, intent(in) :: generation
        class(pedigreeHolder):: this


        integer :: i

        if (generation > indiv%generation) then
            if (indiv%generation /= 0 .and. indiv%generation > 0) then !remove from other list, as has already been set
                call this%generations(indiv%generation)%list_remove(indiv)
            endif

            indiv%generation = generation

            call this%generations(indiv%generation)%list_add(indiv)
            if(indiv%generation > this%maxGeneration) then
                this%maxGeneration = indiv%generation
            endif
        endif
        if ( indiv%nOffs /= 0) then
            do i=1,indiv%nOffs

                call this%setOffspringGeneration(indiv%generation+1,indiv%OffSprings(i)%p)
            enddo
        endif
    end subroutine setOffspringGeneration

    !---------------------------------------------------------------------------
    !< @brief Constructor for recodedPedigreeArray
    !< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
    !< @date   December 22, 2016
    !---------------------------------------------------------------------------
    pure subroutine initRecodedPedigreeArray(this, n)
        implicit none
        class(recodedPedigreeArray), intent(inout) :: This !< @return initialized recoded pedigree array
        integer(int32), intent(in) :: n                    !< number of individuals in pedigree

        this%nInd = n

        if (allocated(this%originalId)) then
            deallocate(this%originalId)
        end if
        allocate(this%originalId(0:n))
        this%originalId = EMPTYID

        if (allocated(this%generation)) then
            deallocate(this%generation)
        end if
        allocate(this%generation(0:n))
        this%generation = 0

        if (allocated(this%id)) then
            deallocate(this%id)
        end if
        allocate(this%id(3, 0:n))
        this%id = 0
    end subroutine

    !---------------------------------------------------------------------------
    !< @brief Destructor for recodedPedigreeArray
    !< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
    !< @date   December 22, 2016
    !---------------------------------------------------------------------------
    pure subroutine destroyRecodedPedigreeArray(this)
        implicit none
        class(recodedPedigreeArray), intent(inout) :: this !< @return recodedPedigreeArray that will be destructed

        if (allocated(this%originalId)) then
            deallocate(this%originalId)
        end if

        if (allocated(this%generation)) then
            deallocate(this%generation)
        end if

        if (allocated(this%id)) then
            deallocate(this%id)
        end if
    end subroutine

    !---------------------------------------------------------------------------
    !< @brief Write recodedPedigreeArray to a file or stdout
    !< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
    !< @date   December 22, 2016
    !---------------------------------------------------------------------------
    subroutine writeRecodedPedigreeArray(this, file)
        use iso_fortran_env, only : output_unit, int32
        implicit none
        class(recodedPedigreeArray), intent(in) :: this !< recodedPedigreeArray that will be written
        character(len=*), intent(in), optional :: file  !< If present File name, else stdout

        integer(int32) :: unit, ind
        character(len=:), allocatable :: fmt

        fmt = "(4i"//Int2Char(IDINTLENGTH)//", a1, 3a"//Int2Char(IDLENGTH)//")"
        if (present(File)) then
          open(newunit=unit, file=trim(file), status="unknown")
        else
          unit = output_unit
        end if
        do ind = 1, this%nInd
          write(unit, fmt) this%id(1:3, ind), this%generation(ind), "", this%originalId(this%id(1:3, ind))
        end do
        if (present(File)) then
          close(unit)
        end if
    end subroutine

    !---------------------------------------------------------------------------
    !< @brief Sorts and recodes pedigree
    !< @details Sorts pedigree such that parents preceede children and recodes ID to 1:n
    !< @author David Wilson david.wilson@roslin.ed.ac.uk & Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
    !< @date   December 20, 2016
    !---------------------------------------------------------------------------
    subroutine makeRecodedPedigreeArray(this, recPed)
        implicit none
        class(pedigreeHolder), intent(inout) :: this      !< object to operate on
        type(recodedPedigreeArray), intent(out) :: recPed !< @return recoded pedigree array
        integer :: counter,i,h
        type(IndividualLinkedListNode), pointer :: tmpIndNode

        call RecPed%init(n=this%pedigreeSize)

        if (.not. allocated(this%generations)) then
            call this%sortPedigreeAndOverwriteWithDummyAtTheTop
        end if

        counter = 0
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                counter = counter + 1
                recPed%originalId(counter) = tmpIndNode%item%originalId
                recPed%generation(counter) = tmpIndNode%item%generation
                recPed%id(1:3,counter) = tmpIndNode%item%getIntegerVectorOfRecodedIds()
                tmpIndNode => tmpIndNode%next
            end do
        end do
    end subroutine makeRecodedPedigreeArray



! TODO get (unique) list of mates. Sorted by Generation.

! function should return offspringlist (=recID) and listOfParents(2, nMatingPairs))

    !---------------------------------------------------------------------------
    !< @brief returns list of mates and offspring for those mate pairs for given pedigree
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    function getAllGenotypesAtPosition(this, position) result(res)

        class(pedigreeHolder) :: this
        integer, intent(in) :: position
        integer(KIND=1), allocatable, dimension(:) :: res
        integer :: counter, i
        allocate(res(this%pedigreeSize))
        res = 9
            counter = 0
        ! TODO can do in parallel
            do i=1, this%pedigreeSize

                if (this%pedigree(i)%isGenotyped()) then
                    counter = counter +1
                    res(counter) = this%pedigree(i)%individualGenotype%getGenotype(position)
                endif

            enddo

    end function getAllGenotypesAtPosition



    !---------------------------------------------------------------------------
    !< @brief returns list of mates and offspring for those mate pairs for given pedigree
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine getMatePairsAndOffspring(this, offSpringList, listOfParents, nMatingPairs)

        use AlphaHouseMod, only : generatePairing

        class(pedigreeHolder), intent(inout) :: this      !< Pedigree object
        integer, dimension(:, :), allocatable, intent(out) :: listOfParents !< indexed by (sire/dam, mateID) = recodedId
        integer, intent(out) :: nMatingPairs
        type(IndividualLinkedList),allocatable, dimension(:) :: offspringList !< list off spring based on index of parents mateID
        
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        type(DictStructure) :: dictionary
        integer(kind=int64) :: tmpPairingKey 
        character(len=20) :: tmpPairingKeyStr ! 19 is largest characters for int64 so 20 just to be safe [and that its round]
        integer :: i,h,j

        dictionary = DictStructure()
        nMatingPairs = 0
        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays 
        endif
        if (allocated(listOfParents)) then
            deallocate(listOfParents)
        endif
        allocate(listOfParents(2,this%pedigreeSize))
        
        if (allocated(offspringList)) then
            deallocate(offspringList)
        endif
        allocate(offspringList(this%pedigreeSize))

        do i=0,this%maxGeneration
        tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                
                do j=1, tmpIndNode%item%nOffs
                    if(associated(tmpIndNode%item,tmpIndNode%item%offsprings(j)%p%sirePointer)) then


                        tmpPairingKey = generatePairing(tmpIndNode%item%offsprings(j)%p%sirePointer%id, tmpIndNode%item%offsprings(j)%p%damPointer%id)
                        write(tmpPairingKeyStr, '(i0)') tmpPairingKey
                        if (.not. dictionary%hasKey(tmpPairingKeyStr)) then
                            nMatingPairs = nMatingPairs + 1
                            listOfParents(1,nMatingPairs) = tmpIndNode%item%id
                            listOfParents(2,nMatingPairs) = tmpIndNode%item%offsprings(j)%p%damPointer%id
                            call offspringList(nMatingPairs)%list_add(tmpIndNode%item%offsprings(j)%p)
                            call dictionary%addKey(tmpPairingKeyStr, nMatingPairs)
                        else
                            call offspringList(dictionary%getValue(tmpPairingKeyStr))%list_add(tmpIndNode%item%offsprings(j)%p)
                        ! TODO make sure using list rather than array here is not terrible
                        ! TODO make sure if only checking sire is good enough
                        endif
                    endif
                enddo
                tmpIndNode => tmpIndNode%next
                
            enddo
        enddo

        call dictionary%destroy()

    end subroutine getMatePairsAndOffspring
        
end module PedigreeModule
