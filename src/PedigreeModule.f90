module PedigreeModule
    use IndividualModule
    use IndividualLinkedListModule

    character, parameter :: EMPTY_PARENT = '0'
    integer, parameter :: IDLENGTH = 32

type PedigreeHolder

    type(Individual), pointer, dimension(:) :: Pedigree !have to use pointer here as otherwise won't let me point to it
    type(IndividualLinkedList) :: Founders !linked List holding all founders
    contains 
        procedure  :: destroyPedigree

end type PedigreeHolder


    interface PedigreeHolder
        module procedure initPedigree
    end interface PedigreeHolder


contains
    

    function initPedigree(fileIn, numberInFile) result(pedStructure)
        use AlphaHouseMod, only : countLines
        use HashModule
        type(PedigreeHolder) :: pedStructure
        character(len=*) :: fileIn
        character(len=IDLENGTH) :: tmpId,tmpSire,tmpDam
        integer(kind=int32),optional :: numberInFile
        integer(kind=int32) :: nIndividuals, fileUnit,tmpSireNum, tmpDamNum
        type(DictStructure) :: dictionary 
        integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
        integer :: tmpAnimalArrayCount = 0
        logical :: sireFound, damFound

        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        allocate(pedStructure%Pedigree(nIndividuals))
        dictionary = DictStructure(nIndividuals) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
        open(newUnit=fileUnit, file=fileIn, status="old")
       
        do i=1,nIndividuals
            
            sireFound = .false.
            damFound = .false.
            
            read(fileUnit,*) tmpId,tmpSire,tmpDam
            call dictionary%addKey(tmpId, i)
            
            pedStructure%Pedigree(i) =  Individual(trim(tmpId),trim(tmpSire),trim(tmpDam), i) !Make a new individual based on info from ped

            if (tmpSire /= EMPTY_PARENT) then !check sire is defined in pedigree
                tmpSireNum = dictionary%getValue(tmpSire)
                if (tmpSireNum /= DICT_NULL) then
                    pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                    call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))
                    call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                    sireFound = .true.
                endif
            endif

            if (tmpDam /= EMPTY_PARENT) then !check dam is defined in pedigree
                tmpDamNum = dictionary%getValue(tmpDam)
                if (tmpDamNum /= DICT_NULL) then !if 
                    pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                    call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))
                    call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                    damFound = .true.
                endif                
            endif
            if (tmpSire == EMPTY_PARENT .and. tmpDam == EMPTY_PARENT) then !if animal is a founder
                call pedStructure%Pedigree(i)%SetObjectAsFounder()
                call pedStructure%Founders%list_add(pedStructure%Pedigree(i))         
            else if (.not. sireFound .and. .not. damFound ) then
                tmpAnimalArrayCount = tmpAnimalArrayCount +1
                tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
            endif 
        enddo


        do i=1,tmpAnimalArrayCount !check animals that didn't have parental information initially
            tmpSire = pedStructure%Pedigree(tmpAnimalArray(i))%getSireId()
            tmpSireNum = dictionary%getValue(tmpSire)
            tmpDam = pedStructure%Pedigree(tmpAnimalArray(i))%getDamId()
            tmpDamNum = dictionary%getValue(tmpDam)
            if (tmpSireNum /= DICT_NULL) then
                pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
            endif
            if (tmpDamNum /= DICT_NULL) then !if 
                    pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                    call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                    call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                    damFound = .true.
            
            else if (tmpSireNum == DICT_NULL .and. tmpDamNum == DICT_NULL) then
                call pedStructure%Pedigree(tmpAnimalArray(i))%SetObjectAsFounder()
                call pedStructure%Founders%list_add(pedStructure%Pedigree(tmpAnimalArray(i)))     
            endif
        enddo
    deallocate(tmpAnimalArray)
    call dictionary%destroy !destroy dictionary as we no longer need it
    end function initPedigree


    subroutine destroyPedigree(this)
        class(PedigreeHolder) :: this
        integer :: i
        
        do i=1,size(this%Pedigree)
            call this%Pedigree(i)%destroyIndividual
            ! call this%Pedigree(i)%Founders
            ! TODO do founders need dsitruction
        enddo
    end subroutine destroyPedigree





end module PedigreeModule