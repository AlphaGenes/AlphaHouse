module PedigreeModule
    use IndividualModule
    use IndividualLinkedListModule
 use HashModule
    character, parameter :: EMPTY_PARENT = '0'
    integer, parameter :: IDLENGTH = 32
    integer, parameter :: generationThreshold = 1000


type PedigreeHolder
    
    type(Individual), pointer, dimension(:) :: Pedigree !have to use pointer here as otherwise won't let me point to it
    type(IndividualLinkedList) :: Founders !linked List holding all founders
    type(IndividualLinkedList),allocatable, dimension(:) :: generations !linked List holding all founders
    type(DictStructure) :: dictionary 
    integer(kind=int32) :: pedigreeSize !pedigree size cannot be bigger than 2 billion animals
    integer(kind=int32) :: maxPedigreeSize

    integer :: maxGeneration
    contains 
        procedure :: destroyPedigree
        procedure :: setPedigreeGenerationsAndBuildArrays
        procedure :: outputSortedPedigree
        procedure :: setOffspringGeneration
        procedure :: addGenotypeInformation
        procedure :: outputSortedPedigreeInAlphaImputeFormat


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
    function initPedigree(fileIn, numberInFile, genderFile) result(pedStructure)
        use AlphaHouseMod, only : countLines
        use HashModule
        use iso_fortran_env
        type(PedigreeHolder) :: pedStructure
        character(len=*),intent(in) :: fileIn
        character(len=*), intent(in),optional :: genderFile
        character(len=IDLENGTH) :: tmpId,tmpSire,tmpDam,tmpCounterStr
        integer(kind=int32),optional,intent(in) :: numberInFile
        integer(kind=int32) :: stat, fileUnit,tmpSireNum, tmpDamNum, tmpGender,tmpIdNum
        integer(kind=int64) :: nIndividuals
        integer :: tmpCounter = 0
        integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
        integer :: tmpAnimalArrayCount = 0
        integer :: i
        integer(kind=int64) :: sizeDict
        logical :: sireFound, damFound

        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        sizeDict = nIndividuals
        pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals)
        allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))   
        pedStructure%pedigreeSize = nIndividuals
        pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
        allocate(tmpAnimalArray(nIndividuals+10)) !allocate to nIndividuals in case all animals are in incorrect order of generations
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
                    pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                    if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                        write(error_unit,*) "ERROR: too many undefined animals"
                        stop
                        ! TODO do a move alloc here to avoid this hack
                    endif
                    pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual("dum"//trim(tmpSire),'0','0', pedStructure%pedigreeSize)
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                    pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                    call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                    call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.  
                endif
                sireFound = .true.
            endif

            if (tmpDam /= EMPTY_PARENT) then
                if (tmpDamNum /= DICT_NULL) then !if dam has been found            
                        pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                        call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female       
                else
                    pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                    if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                        write(error_unit,*) "ERROR: too many undefined animals"
                        stop
                    endif
                    pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual("dum"//trim(tmpDam),'0','0', pedStructure%pedigreeSize)
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                    pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                    call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))

                    call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                    pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true. 
                endif
                damFound = .true.
            endif


            if (.not. damFound .XOR. .not. sireFound) then
                tmpCounter =  tmpCounter + 1
                write(tmpCounterStr, '(I3.3)') tmpCounter
                pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
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



    deallocate(tmpAnimalArray)
    
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

    subroutine addGenotypeInformation(this, genotypeFile, nsnps, nAnnisG)
         
        use AlphaHouseMod, only : countLines
        implicit none
        class(PedigreeHolder) :: this
        character(len=*) :: genotypeFile
        character(len=IDLENGTH) :: tmpID
        integer,intent(in) :: nsnps
        integer,intent(in),optional :: nAnnisG
        integer, allocatable, dimension(:) :: tmpSnpArray
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
                    allocate(this%pedigree(tmpIdNum)%genotype(nsnps))
                    this%pedigree(tmpIdNum)%genotype = tmpSnpArray
                endif
            enddo



    end subroutine addGenotypeInformation
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
        
        do i=0, this%maxGeneration
            tmpIndNode => this%generations(i)%first
            do h=1, this%generations(i)%length
                    write(unit,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId, tmpIndNode%item%generation
                    ! write(*,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
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
        
        block 
            integer :: sireId, damId
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

                     write (unit,'(3i20,a20)') tmpIndNode%item%id,sireId,damId, tmpIndNode%item%originalID
                    ! write(*,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
                    tmpIndNode => tmpIndNode%next
            end do
        enddo
        endblock
        if (present(file)) then !avoids closing stdout
            close(unit)
        endif
    end subroutine outputSortedPedigreeInAlphaImputeFormat
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





end module PedigreeModule