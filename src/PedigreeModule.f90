module PedigreeModule
    use IndividualModule
    use IndividualLinkedListModule

    character, parameter :: EMPTY_PARENT = '0'
    integer, parameter :: IDLENGTH = 32
    integer, parameter :: generationThreshold = 1000



type PedigreeHolder

    type(Individual), pointer, dimension(:) :: Pedigree !have to use pointer here as otherwise won't let me point to it
    type(IndividualLinkedList) :: Founders !linked List holding all founders
    type(IndividualLinkedList),allocatable, dimension(:) :: generations !linked List holding all founders
    integer :: maxGeneration
    contains 
        procedure :: destroyPedigree
        procedure :: setPedigreeGenerationsAndBuildArrays
        procedure :: outputSortedPedigree
        procedure :: setOffspringGeneration

end type PedigreeHolder

    interface PedigreeHolder
        module procedure initPedigree
    end interface PedigreeHolder

contains
    

        !---------------------------------------------------------------------------
    !> @brief distructor for pedigree class 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] file path (string)
    !---------------------------------------------------------------------------
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
        integer :: i
        logical :: sireFound, damFound

        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        allocate(pedStructure%Pedigree(nIndividuals))
        
        dictionary = DictStructure(nIndividuals) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
        pedStructure%maxGeneration = 0
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

        !check animals that didn't have parental information initially
        ! this is done to avoid duplication when a pedigree is sorted
        do i=1,tmpAnimalArrayCount 
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


    !---------------------------------------------------------------------------
    !> @brief distructor for pedigree class 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] file path (string)
    !---------------------------------------------------------------------------
    subroutine destroyPedigree(this)
        class(PedigreeHolder) :: this
        integer :: i
        
        do i=1,size(this%Pedigree)
            call this%Pedigree(i)%destroyIndividual
            if (allocated(this%generations)) then
                deallocate(this%generations)
            endif
            ! call this%Pedigree(i)%Founders
            ! TODO do founders need dsitruction
        enddo
    end subroutine destroyPedigree


    !---------------------------------------------------------------------------
    !> @brief builds correct generation information by looking at founders 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
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
    !> @brief writes sorted pedigree information to either a file or stdout
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] file path (string)
    !---------------------------------------------------------------------------
    subroutine outputSortedPedigree(this,file)

        use iso_fortran_env, only : output_unit
        class(PedigreeHolder) :: this
        character(len=*), intent(in), optional :: file
        integer :: unit, i,h
        type(IndividualLinkedListNode), pointer :: tmpIndNode
        
        if (.not. allocated(this%generations)) then
            call this%setPedigreeGenerationsAndBuildArrays
        endif
        print *,"finished building generation arrays"

        if (present(file)) then
            open(newUnit=unit, file=file, status="unknown")
        else
            unit = output_unit
        endif
        
        do i=1, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length-1
                ! write(unit,*) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId, tmpIndNode%item%generation
                print *,"generation:",this%generations(i)%length,h
                print *,tmpIndNode%item%nOffs
                print *,tmpIndNode%item%originalID
                tmpIndNode => tmpIndNode%next
            end do
        enddo
        if (present(file)) then
            close(unit)
        endif
    end subroutine outputSortedPedigree


    !---------------------------------------------------------------------------
    !> @brief Sets generation of an individual and his children recursively
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] generation (integer) 
    !> @param[in] pointer to an individual
    !---------------------------------------------------------------------------
    recursive subroutine setOffspringGeneration(this,generation, indiv)
        type(Individual),pointer, intent(inout) :: indiv
        integer, intent(in) :: generation
        class(pedigreeHolder):: this

        integer :: i
        
        if (indiv%generation /= 0) then !remove from other list, as has already been set
            call this%generations(indiv%generation)%list_remove(indiv)
        endif
        indiv%generation = generation
        call this%generations(indiv%generation)%list_add(indiv)
        if(indiv%generation > this%maxGeneration) then
            this%maxGeneration = indiv%generation
        endif
        if ( indiv%nOffs /= 0) then
            do i=1,indiv%nOffs
                
                call this%setOffspringGeneration(indiv%generation+1,indiv%OffSprings(i)%p)
            enddo
        endif
    end subroutine setOffspringGeneration





end module PedigreeModule