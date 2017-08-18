
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!<@file     PedigreeModule.f90
!
! DESCRIPTION:
!<@brief    Module containing logic of Pedigree
!<@details  Module contains functions to read in and set up pedgiree data structure
!
!<@author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!<@date     January 4, 2017
!
!<@version  1.0.0
!
!
!-------------------------------------------------------------------------------




module PedigreeModule
    use IndividualModule
    use IndividualLinkedListModule
    use HashModule
    use constantModule
    use AlphaHouseMod, only : Int2Char
    implicit none


    private addOffspringsAfterReadIn

    type PedigreeHolder

        type(Individual), pointer, dimension(:) :: Pedigree !have to use pointer here as otherwise won't let me point to it
        type(IndividualLinkedList) :: Founders !linked List holding all founders
        type(IndividualLinkedList),allocatable, dimension(:) :: generations !linked List holding each generation
        type(DictStructure) :: dictionary ! hashmap of animal ids to index in pedigree
        integer(kind=int32) :: pedigreeSize, nDummys !pedigree size cannot be bigger than 2 billion animals
        integer(kind=int32) :: unknownDummys !< dummys that have been set by having one unknown parent
        integer(kind=int32) :: maxPedigreeSize ! maximum size pedigree can be

        integer, dimension(:) , allocatable :: inputMap ! map going from inputMap(1:pedsize) = recID -- this is to maintain the order of the original pedigree

        integer, dimension(:) , allocatable :: genotypeMap ! map going from genotypeMap(1:nAnisG) = recID
        type(DictStructure) :: genotypeDictionary ! maps id to location in genotype map
        integer(kind=int32) :: nGenotyped ! number of animals that are genotyped

        integer, dimension(:) , allocatable :: hdMap ! map going from genotypeMap(1:nHd) = recID
        type(DictStructure) :: hdDictionary ! maps id to location in genotype map
        integer(kind=int32) :: nHd ! number of animals that are genotyped hd
        integer :: maxGeneration ! largest generation

        integer :: nsnpsPopulation ! number of snps for 
        integer(kind=1) :: isSorted !integer saying how pedigree is sorted. 0 is not sorted. 1 is sorted with all dummys at end, 2 is sorted with unknown dummys at end, 3 is sorted with dummys at top 

        type(IndividualLinkedList) :: sireList, damList !< lists containing all sires and dams

    contains
        procedure :: calculatePedigreeCorrelationNoInbreeding
        procedure :: copyPedigree
        procedure :: destroyPedigree
        procedure :: setPedigreeGenerationsAndBuildArrays
        procedure :: outputSortedPedigree
        procedure :: setOffspringGeneration
        procedure :: addGenotypeInformationFromFile
        procedure :: addPhaseInformationFromFile
        procedure :: addGenotypeInformationFromArray
        generic :: addGenotypeInformation => addGenotypeInformationFromArray, addGenotypeInformationFromFile
        procedure :: outputSortedPedigreeInAlphaImputeFormat
        procedure :: isDummy
        procedure :: sortPedigreeAndOverwrite
        procedure :: sortPedigreeAndOverwriteWithDummyAtTheTop
        procedure :: makeRecodedPedigreeArray
        procedure :: printPedigree
        procedure :: printPedigreeOriginalFormat
        procedure :: getMatePairsAndOffspring
        procedure :: getAllGenotypesAtPosition
        procedure :: getAllGenotypesAtPositionWithUngenotypedAnimals
        procedure :: getPhaseAtPosition
        procedure :: getPhaseAtPositionUngenotypedAnimals
        procedure :: setAnimalAsGenotyped
        procedure :: getGenotypesAsArray
        procedure :: getPhaseAsArray
        procedure :: getGenotypesAsArrayWitHMissing
        procedure :: countMissingGenotypesNoDummys  
        procedure :: setGenotypeFromArray
        procedure :: setPhaseFromArray
        procedure :: getNumGenotypesMissing
        procedure :: getGenotypedFounders
        procedure :: getSireDamGenotypeIDByIndex
        procedure :: setAnimalAsHD
        procedure :: getSireDamHDIDByIndex
        procedure :: getGenotypePercentage
        procedure :: writeOutGenotypes
        procedure :: WriteoutPhase
        procedure :: createDummyAnimalAtEndOfPedigree
        procedure :: addAnimalAtEndOfPedigree
        procedure :: addSequenceFromFile
        procedure :: setAnimalAsGenotypedSequence
        procedure :: convertSequenceDataToArray
        procedure :: getSequenceAsArrayWithMissing
        procedure :: writeOutGenders
        procedure :: readInGenders
        procedure :: MakeGenotype
        procedure :: PhaseComplement
        procedure :: homozygoticFillIn
        procedure :: wipeGenotypeAndPhaseInfo
    end type PedigreeHolder

    ! TODO - overload == and = functions

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
                        module procedure initPedigreeIntArrays
                            ! initPedigreeFromOutputFileFolder initialiser from file
                        end interface PedigreeHolder

                        interface Sort !Sorts into generation list
                            module procedure :: setPedigreeGenerationsAndBuildArrays
                            end interface Sort
                        contains


                            !---------------------------------------------------------------------------
                            !< @brief constructor for creating an empty pedgiree
                            !< @details Constructs an empty pedigree 
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initEmptyPedigree(nsnps) result(pedStructure)
                                use iso_fortran_env
                                type(PedigreeHolder) :: pedStructure
                                integer, optional :: nsnps

                                pedStructure%dictionary = DictStructure()
                                pedStructure%pedigreeSize = 0
                                pedStructure%nDummys = 0
                                pedStructure%nGenotyped = 0
                                pedStructure%nHd = 0
                                pedStructure%maxPedigreeSize = DEFAULTDICTSIZE
                                pedStructure%nsnpsPopulation = 0

                                if (present(nsnps)) then
                                    pedStructure%nsnpsPopulation = nsnps
                                endif
                            end function initEmptyPedigree


                            !---------------------------------------------------------------------------
                            !< @brief Sorts pedigree, and overwrites all fields to new values
                            !< @details effectively, does a deep copy to sort pedigree based on generation, but puts dummys at bottom
                            !< If value is given for unknownDummysAtEnd, then only unknown dummys will be put at the end
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function copyPedigree(this) result(res)
                                use iso_fortran_env, only : output_unit, int64
                                class(PedigreeHolder),intent(in) :: this
                                type(pedigreeHolder) :: res
                                integer :: i, tmpId,tmpGenotypeMapIndex
                                integer(kind=int64) :: sizeDict
                                type(Individual), pointer, dimension(:) :: newPed
                                type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList

                                if (.not. allocated(this%generations)) then
                                    call this%setPedigreeGenerationsAndBuildArrays
                                endif
                                i = 0


                                res%pedigreeSize = this%pedigreeSize
                                res%nDummys = this%nDummys
                                res%unknownDummys = this%unknownDummys
                                res%maxPedigreeSize = this%maxPedigreeSize
                                res%inputMap = this%inputMap
                                res%genotypeMap = this%genotypeMap
                                call copyHashTable(res%genotypeDictionary, this%genotypeDictionary)
                                res%nGenotyped =this%nGenotyped
                                res%hdMap = this%hdMap
                                call copyHashTable(res%hdDictionary, this%hdDictionary)
                                res%nHd = this%nHd
                                res%maxGeneration = this%maxGeneration
                                res%nsnpsPopulation = this%nsnpsPopulation
                                res%isSorted =this%isSorted

                                ! call res%dictionary%destroy()
                                ! call res%founders%destroyLinkedList()
                                ! call res%sireList%destroyLinkedList()
                                ! call res%damList%destroyLinkedList()


                                sizeDict  =this%pedigreeSize
                                allocate(newPed(this%maxPedigreeSize))

                                res%dictionary = DictStructure(sizeDict)
                                allocate(newGenerationList(0:this%maxGeneration))

                                do i=1, this%pedigreeSize
                                    ! update dictionary index
                                    call res%dictionary%addKey(this%pedigree(i)%originalID,i)

                                    !  update genotype map

                                    if (res%nGenotyped > 0) then
                                        tmpGenotypeMapIndex = res%genotypeDictionary%getValue(this%pedigree(i)%originalID)
                                        if (tmpGenotypeMapIndex /= DICT_NULL) then
                                            res%genotypeMap(tmpGenotypeMapIndex) = i
                                        endif
                                    endif
                                    ! Update hd map
                                    if (res%nHd > 0) then
                                        tmpGenotypeMapIndex = res%hdDictionary%getValue(this%pedigree(i)%originalID)
                                        if (tmpGenotypeMapIndex /= DICT_NULL) then
                                            res%hdMap(tmpGenotypeMapIndex) = i
                                        endif
                                    endif
                                    ! Set the location to the pedigree to the new value
                                    newPed(i) = this%pedigree(i)

                                    ! take the original id, and update it
                                    if(.not. newPed(i)%isDummy) then
                                        res%inputMap(newPed(i)%originalPosition) = i
                                    endif
                                    newPed(i)%id = i
                                    call newPed(i)%resetOffspringInformation ! reset offsprings
                                    if (associated(newPed(i)%sirePointer)) then
                                        tmpId =  res%dictionary%getValue(newPed(i)%sirePointer%originalID)
                                        if (tmpID /= DICT_NULL) then
                                            call newPed(tmpId)%addOffspring(newPed(i))
                                            newPed(i)%sirePointer=> newPed(tmpId)
                                            if (newPed(tmpId)%nOffs == 1) then
                                                call res%sireList%list_add(newPed(tmpId))
                                            endif
                                        endif
                                    endif

                                    if (associated(newPed(i)%damPointer)) then
                                        tmpId =  res%dictionary%getValue(newPed(i)%damPointer%originalID)
                                        if (tmpID /= DICT_NULL) then
                                            call newPed(tmpId)%addOffspring(newPed(i))
                                            newPed(i)%damPointer=> newPed(tmpId)
                                            if (newPed(tmpId)%nOffs == 1) then
                                                call res%damList%list_add(newPed(tmpId))
                                            endif
                                        endif
                                    endif

                                    if (i ==0 ) then !if object is afounder add to founder array
                                        call  res%founders%list_add(newPed(i))
                                    endif

                                    if (res%isSorted /= 0) then
                                        call newGenerationList(newPed(i)%generation)%list_add(newPed(i))
                                    endif
                                enddo

                                res%pedigree => newPed
                                res%generations = newGenerationList

                            end function copyPedigree



                            subroutine findMendelianInconsistencies(ped, threshold)
                                implicit none
                                type(pedigreeHolder) :: ped
                                real, intent(in) :: threshold !< threshold for removing animals
                                integer :: sireG, damG
                                type(individual),pointer :: sire, dam
                                integer :: i,g,j
                                integer :: GenConflict(3,ped%pedigreesize)   ! score for 1) paternal link, 2) maternal link, and 3) sum of both
                                integer :: PedConflict(4,ped%pedigreesize)   ! count of father-individual and mother-individual genotype combinations available (1,2) and conflicts thereof (3,4)

                                do j=1,ped%pedigree(i)%individualGenotype%length

                                    do i=1,ped%pedigreeSize


                                        if (associated(ped%pedigree(i)%sirePointer)) then

                                            sire => ped%pedigree(i)%sirePointer
                                            dam => ped%pedigree(i)%damPointer
                                            g = ped%pedigree(i)%individualGenotype%getGenotype(j)
                                            sireG = sire%individualGenotype%getGenotype(j)
                                            damG = dam%individualGenotype%getGenotype(j)
                                            select case(g)


                                            case(0)
                                                if (sireG == 2) then
                                                    if (damG == 2) then
                                                        write(error_unit, *) "not valid - both"
                                                    else
                                                        write(error_unit, *) "not valid paternal"
                                                        PedConflict(1:2,i)=PedConflict(1:2,i)+1; GenConflict(1:2,i)=(/+1,-1/);
                                                    endif
                                                else if (sireG == 9 .and. damG == 2) then
                                                    GenConflict(1:2,i)=(/-1,+1/)
                                                    write(error_unit, *) "not valid maternal"
                                                else
                                                    if (.not. (damG == 0 .or. damG == 1 .or. damG == 9)) then
                                                        ! TODO add removal
                                                        GenConflict(1:2,i)=(/-1,+1/)
                                                        write(error_unit, *) "not valid maternal"
                                                    endif
                                                endif


                                            case(1)

                                                if (sireG == 1 .or. sireG == 9) then 
                                                    ! valid                                                    ! 
                                                else if (sireG == 0) then
                                                    if ( .not. (damG == 1 .or. damG == 2 .or. damG == 9)) then
                                                        ! TODO not valid
                                                    endif

                                                else if (sireG == 2) then
                                                    if (.not. (damG == 0 .or. damG == 1 .or. damG == 9)) then
                                                        ! TODO not valid

                                                    endif
                                                endif

                                            case(2)

                                                if (sireG == 9 .or. sireG == 1 .or. sireG == 2) then

                                                    if (damG == 0) then
                                                        ! maternal conflict

                                                    endif
                                                else if (sireG == 0) then

                                                    if (damG == 0) then
                                                        ! both paternal and maternal conflict
                                                    else
                                                        ! paternal conflict
                                                    endif

                                                else 
                                                    write(error_unit, *) "ERROR - invalid genotpe configuration. ABORT"
                                                endif

                                            end select

                                            GenConflict(3,i)=sum(GenConflict(1:2,i))
                                            GenConflict(3,sire%id)=GenConflict(3,sire%id)+GenConflict(1,i)
                                            GenConflict(3,dam%id)=GenConflict(3,dam%id)+GenConflict(2,i)
                                        else 
                                            write(error_unit,*) "Animal is founder"
                                            ! cycle
                                        endif
                                    enddo !< pedigree loop
                                    ! call removeGenotype(ped,position, genConflict)
                                enddo

                            end subroutine findMendelianInconsistencies


                            ! subroutine removeGenotypeMendellian(ped,position, genConflict)
                          
                            !     ! integer, dimension(:,:), allocatable :: pedConflicts !< format pedConflicts(pedsize, 1,2)




                            !     do i=1,ped%pedigreeSize

                            !         if (ped%pedigree(i)%founder) cycle

                            !         damId = ped%pedigree(i)%damId%id
                            !         sireId = ped%pedigree(i)%sirePointer%id

                            !         if (GenConflict(3,sireId) < GenConflict(3,i) then
                            !             ped%pedigree(i)%setGenotype(position,9)

                            !         else if (GenConflict(3,sireId) == GenConflict(3,i) then
                            !             ! set both to missing
                            !             ped%pedigree(sireId)%setGenotype(position,9)
                            !             ped%pedigree(i)%setGenotype(position,9)
                            !         else
                            !             ped%pedigree(sireId)%setGenotype(position,9)
                            !         endif

                            ! end subroutine       

                            !---------------------------------------------------------------------------
                            !< @brief Sets phase information from an array
                            !< @details sets phase objects for animal when given an array
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            pure subroutine setPhaseFromArray(this, array)

                                use HaplotypeModule

                                class(PedigreeHolder), intent(inout) :: this
                                integer(kind=1),dimension(:,:,:),intent(in) :: array !< array should be of format (recodedindId, snp, allele )
                                integer :: i

                                do i=1, size(array,1)
                                    this%pedigree(i)%individualPhase(1) = newHaplotypeInt(array(i,:,1))
                                    this%pedigree(i)%individualPhase(2) = newHaplotypeInt(array(i,:,2))
                                enddo

                            end subroutine setPhaseFromArray


                            subroutine homozygoticFillIn(this)

                                class(PedigreeHolder), intent(inout) :: this
                                ! Impute phase information for homozygous cases

                                integer :: i,j
                                !$OMP PARALLEL DO &
                                !$OMP PRIVATE(i,j)
                                do i=1, this%pedigreeSize
                                    do j=1, this%pedigree(this%genotypeMap(1))%individualGenotype%length
                                        if (this%pedigree(i)%individualGenotype%getGenotype(j) == 2) then
                                            call this%pedigree(i)%individualPhase(1)%setPhase(j,1)
                                            call this%pedigree(i)%individualPhase(2)%setPhase(j,1)
                                        else if (this%pedigree(i)%individualGenotype%getGenotype(j) == 0) then
                                            call this%pedigree(i)%individualPhase(1)%setPhase(j,0)
                                            call this%pedigree(i)%individualPhase(2)%setPhase(j,0)
                                        endif
                                    enddo
                                enddo
                                !$OMP END PARALLEL DO

                            end subroutine homozygoticFillIn

                            !---------------------------------------------------------------------------
                            !<@brief Calculate correlation between pedigree 
                            !<Calculates the correlation between individuals in the pedigree, assuming no inbreeding.   Will also handle groups - just give
                            !<the group parents a unique ID.   It assumes that the pedigree has been sorted such that the unknown dummies are at the end
                            !<(sortPedigreeAndOverwrite(1)).
                            !<@author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
                            !---------------------------------------------------------------------------
                            function calculatePedigreeCorrelationNoInbreeding(pedigreeIn, additiveVarianceIn) result(valuesOut)
                                class(PedigreeHolder), intent(in):: pedigreeIn
                                real(real64), dimension(:,:), allocatable:: valuesOut
                                real(real64), intent(in), optional:: additiveVarianceIn
                                integer:: damKnown, sireKnown
                                integer::numLevels, ID, sireID, damID
                                integer:: i
                                real(real64):: D

                                numLevels = pedigreeIn%pedigreeSize-pedigreeIn%UnknownDummys

                                allocate(valuesOut(numLevels, numLevels))
                                valuesOut = 0_real64

                                do i = 1,numLevels
                                    damKnown = 0
                                    sireKnown = 0
                                    damID=0
                                    sireID=0
                                    ID = pedigreeIn%pedigree(i)%ID
                                    if (associated(pedigreeIn%pedigree(i)%sirePointer)) then
                                        if (.not. pedigreeIn%pedigree(i)%sirePointer%isUnknownDummy) then
                                            sireKnown = 1
                                            sireID = pedigreeIn%pedigree(i)%sirePointer%ID
                                        end if
                                    end if

                                    if (associated(pedigreeIn%pedigree(i)%damPointer)) then
                                        if (.not. pedigreeIn%pedigree(i)%damPointer%isUnknownDummy) then
                                            damKnown = 1
                                            damID = pedigreeIn%pedigree(i)%damPointer%ID
                                        end if
                                    end if

                                    D = 4.0_real64/(2.0_real64+damKnown+sireKnown) !If both are known, then this is set to be 2,if only 1 is known this is 4.0/3.0

                                    valuesOut(ID, ID) = valuesOut(ID,ID)+D

                                    if (sireKnown .ne. 0) then
                                        valuesOut(sireID, ID) = valuesOut(sireID, ID) -0.5*D
                                        valuesOut(ID, sireID) = valuesOut(ID, sireID) - 0.5*D
                                        valuesOut(sireID, sireID) = valuesOut(sireId, sireID) +0.25*D
                                    end if

                                    if (damKnown .ne. 0) then
                                        valuesOut(damID, ID) = valuesOut(damID, ID) -0.5*D
                                        valuesOut(ID, damID) = valuesOut(ID, damID) -0.5*D
                                        valuesOut(damID, damID) = valuesOut(damId, damID) +0.25*D
                                    end if

                                    if ((sireKnown .ne. 0) .and. (damKnown.ne.0)) then
                                        valuesOut(sireID, damID) = valuesOut(sireID, damID) + 0.25*D
                                        valuesOut(damID, sireID) = valuesOut(damID, sireID) +0.25*D
                                    end if
                                end do

                                if (present(additiveVarianceIn)) then
                                    valuesOut = valuesOut*additiveVarianceIn
                                end if

                            end function calculatePedigreeCorrelationNoInbreeding



                            !---------------------------------------------------------------------------
                            !< @brief Sets pedigree genotype from integer array 
                            !< @details Sets animals genotype from integer array
                            !< NOTE: does not set animals as genotyped
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------  
                            subroutine setGenotypeFromArray(this, array)
                                use GenotypeModule

                                class(PedigreeHolder), intent(inout)  :: this
                                integer(kind=1),dimension(:,:) :: array !< array should be of format (recodedindId, snp)
                                integer :: i

                                do i=1, size(array,1)
                                    this%pedigree(i)%individualGenotype = newGenotypeInt(array(i,:))
                                enddo

                            end subroutine setGenotypeFromArray


                            !---------------------------------------------------------------------------
                            !< @brief Constructor for pedigree class
                            !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigree(fileIn, numberInFile, genderFile, nsnps) result(pedStructure)
                                use AlphaHouseMod, only : countLines
                                use iso_fortran_env
                                type(PedigreeHolder) :: pedStructure
                                character(len=*),intent(in) :: fileIn !< path of pedigree file
                                character(len=*), intent(in),optional :: genderFile !< path to gender file
                                integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
                                integer, optional, intent(in) :: nsnps !< number of snps for the population

                                character(len=IDLENGTH) :: tmpId,tmpSire,tmpDam
                                integer(kind=int32) :: stat, fileUnit,tmpSireNum, tmpDamNum, tmpGender,tmpIdNum
                                integer(kind=int64) :: nIndividuals
                                integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
                                integer :: tmpAnimalArrayCount,i
                                integer(kind=int64) :: sizeDict
                                logical :: sireFound, damFound

                                pedStructure%isSorted = 0
                                pedStructure%nDummys = 0
                                pedStructure%nHd = 0
                                tmpAnimalArrayCount = 0
                                pedStructure%nGenotyped = 0
                                pedStructure%nsnpsPopulation = 0

                                if (present(nsnps)) then
                                    pedStructure%nsnpsPopulation = nsnps
                                endif
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
                                allocate(pedStructure%inputMap(nIndividuals))
                                pedStructure%maxGeneration = 0
                                open(newUnit=fileUnit, file=fileIn, status="old")

                                do i=1,nIndividuals
                                    sireFound = .false.
                                    damFound = .false.
                                    read(fileUnit,*) tmpId,tmpSire,tmpDam

                                    if (trim(tmpId) == trim(tmpsire) .or. trim(tmpDam) == trim(tmpId)) then

                                        write(error_unit,*) "Error: Animal ", trim(tmpId), " has been given itself as a parent. please fix this. Error on line:",i
                                        stop
                                    endif
                                    call pedStructure%dictionary%addKey(tmpId, i)

                                    pedStructure%Pedigree(i) = Individual(trim(tmpId),trim(tmpSire),trim(tmpDam), i, nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
                                    pedStructure%Pedigree(i)%originalPosition = i
                                    pedStructure%inputMap(i) = i
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

                                        if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                                            call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum))
                                        endif
                                        pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))
                                        if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                                            call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum))
                                        endif
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
                                            write(error_unit, *) "WARNING: Gender  defined for an animal that does not exist in Pedigree!"
                                            write(error_unit, *) "Amimal:",tmpId
                                        endif
                                    end do
                                endif
                                call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)
                                deallocate(tmpAnimalArray)
                                !  write(output_unit, *) "Number of animals in Pedigree:",pedStructure%pedigreeSize-pedStructure%nDummys
                                !  write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys

                            end function initPedigree


                            !---------------------------------------------------------------------------
                            !< @brief Constructor for pedigree class taking in files that have been written out
                            !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigreeFromOutputFile(pedigreeFile, genderFile, pedGenotypeFile, phaseFile, nsnps) result(pedStructure)

                                character(len=*), intent(in) :: pedigreeFile
                                character(len=*), intent(in), optional :: genderFile, pedGenotypeFile, phaseFile
                                integer, intent(in) :: nsnps
                                type(PedigreeHolder) :: pedStructure


                                if (present(genderFile)) then
                                    pedStructure = initPedigree(pedigreeFile, nsnps=nsnps, genderFile=genderFile)
                                else 
                                    pedStructure = initPedigree(pedigreeFile, nsnps=nsnps)

                                endif

                                if (present(pedGenotypeFile)) then
                                    call pedStructure%addGenotypeInformationFromFile(pedGenotypeFile, nsnps)
                                endif

                                if (present(phaseFile)) then
                                    call pedStructure%addPhaseInformationFromFile(phaseFile, nsnps)
                                endif


                            end function initPedigreeFromOutputFile


                            !---------------------------------------------------------------------------
                            !< @brief creates a pedigree object from input folder
                            !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigreeFromOutputFileFolder(folder, nsnps) result(pedStructure)

                                character(len=*), intent(in) :: folder !< input folder
                                integer, intent(in) :: nsnps
                                type(PedigreeHolder) :: pedStructure
                                character(len=:), allocatable :: pedigreeFile, genotypeFile, phaseFile,genderFile


                                pedigreeFile = "pedigree.txt"
                                genotypeFile = "genotypes.txt"
                                phaseFile = "phase.txt"
                                genderFile = "gender.txt"
                                pedStructure = initPedigree(folder//pedigreeFile, nsnps=nsnps, genderFile=folder//genderFile)
                                call pedStructure%addGenotypeInformationFromFile(folder//genotypeFile, nsnps)
                                call pedStructure%addPhaseInformationFromFile(folder//phaseFile, nsnps)


                            end function initPedigreeFromOutputFileFolder


                            !---------------------------------------------------------------------------
                            !< @brief writes out the pedigree to predefined filess
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !---------------------------------------------------------------------------
                            subroutine writeOutPedigree(this,pedigreeFolder)
                                class(PedigreeHolder) :: this
                                character(len=*), intent(in), optional :: pedigreeFolder
                                character(len=:), allocatable :: pedigreeFile, genotypeFile, phaseFile,genderFile


                                pedigreeFile = "pedigree.txt"
                                genotypeFile = "genotypes.txt"
                                phaseFile = "phase.txt"
                                genderFile = "gender.txt"


                                call this%printPedigreeOriginalFormat(pedigreeFolder//pedigreeFile)
                                call this%WriteoutPhase(pedigreeFolder//phaseFile)
                                call this%writeOutGenotypes(pedigreeFolder//genotypeFile) !< just output genotype information for animals that are set to genotpyyed
                                call this%writeOutGenders(pedigreeFolder//genderFile)

                            end subroutine writeOutPedigree




                            !---------------------------------------------------------------------------
                            !< @brief Constructor for pedigree class using Genotype File format
                            !< @details Constructor builds pedigree, without any sorting being done.
                            !< If no pedigree file is supplied, all animals are founders
                            !< If an animal is in the pedigree, but not in the genotypeFile, this animal is still created as a dummy!
                            !< If the animal is in the genotype file, but not in the pedigree, it is added!
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigreeGenotypeFiles(fileIn, numberInFile, nSnp,GenotypeFileFormatIn, pedFile, genderfile) result(pedStructure)
                                use AlphaHouseMod, only : countLines
                                use iso_fortran_env
                                type(PedigreeHolder) :: pedStructure
                                integer, intent(in) :: nSnp !< number of snps to read
                                integer , intent(in), optional :: GenotypeFileFormatIn
                                character(len=*),intent(in) :: fileIn !< path of Genotype file
                                integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
                                character(len=*),intent(in), optional :: pedFile !< path of pedigree file
                                character(len=*),intent(in), optional :: genderfile !< path of gender file
                                character(len=IDLENGTH) :: tmpId
                                integer(kind=int32) :: fileUnit
                                integer(kind=int64) :: nIndividuals
                                integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
                                integer :: tmpAnimalArrayCount
                                integer :: i,j,GenotypeFileFormat
                                integer(kind=1), dimension(:), allocatable :: tmpGeno
                                integer(kind=int64) :: sizeDict
                                integer(kind=1), dimension(nSnp * 2) :: WorkVec

                                allocate(tmpGeno(nsnp))
                                pedStructure%isSorted = 0
                                pedStructure%nHd = 0
                                pedStructure%nGenotyped = 0
                                pedStructure%nsnpsPopulation = nsnp

                                if (present(numberInFile)) then
                                    nIndividuals = numberInFile
                                else
                                    nIndividuals = countLines(fileIn)
                                endif

                                if  (present(pedFile)) then
                                    if (present(genderFile)) then
                                        pedStructure = PedigreeHolder(pedFile, genderFile=genderFile, nsnps=nSnp)
                                    else
                                        pedStructure = PedigreeHolder(pedFile, nsnps=nSnp)
                                    endif
                                else
                                    pedStructure%nDummys = 0
                                    tmpAnimalArrayCount = 0
                                    sizeDict = nIndividuals
                                    pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
                                    allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
                                    pedStructure%pedigreeSize = nIndividuals
                                    pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
                                    allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
                                    allocate(pedStructure%inputMap(nIndividuals))
                                    pedStructure%maxGeneration = 0
                                endif

                                if (present(GenotypeFileFormatIn)) then
                                    GenotypeFileFormat = GenotypeFileFormatIn
                                else
                                    GenotypeFileFormat = 1
                                endif

                                open(newUnit=fileUnit, file=fileIn, status="old")

                                do i=1,nIndividuals
                                    select case(GenotypeFileFormat)
                                    case(1)
                                        read(fileUnit,*) tmpId,tmpGeno(:)
                                    case(2)
                                        read (fileUnit, *) tmpId, tmpGeno(:)
                                        read (fileUnit, *) tmpId, tmpGeno(:)
                                    case(3)
                                        read(fileUnit,*) tmpId,WorkVec(:)
                                    end select

                                    if (GenotypeFileFormat == 3) then
                                        do j=1, nsnp
                                            tmpGeno(j) = MissingGenotypeCode
                                            if ((WorkVec(j*2 - 1) == 1).and.(WorkVec(j*2) == 1)) tmpGeno(j) = 0
                                            if ((WorkVec(j*2 - 1) == 1).and.(WorkVec(j*2) == 2)) tmpGeno(j) = 1
                                            if ((WorkVec(j*2 - 1) == 2).and.(WorkVec(j*2) == 1)) tmpGeno(j) = 1
                                            if ((WorkVec(j*2 - 1) == 2).and.(WorkVec(j*2) == 2)) tmpGeno(j) = 2
                                        enddo
                                    endif
                                    if (present(pedFile)) then
                                        j = pedStructure%dictionary%getValue(tmpID)

                                        if ( j == DICT_NULL) then
                                            call pedStructure%addAnimalAtEndOfPedigree(tmpID,tmpGeno)
                                        else
                                            call pedStructure%setAnimalAsGenotyped(j,tmpGeno)
                                            call pedStructure%setAnimalAsHd(j)
                                        endif
                                    else
                                        call pedStructure%dictionary%addKey(tmpId, i)
                                        pedStructure%Pedigree(i) =  Individual(trim(tmpId),"0","0", i, nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
                                        pedStructure%Pedigree(i)%originalPosition = i
                                        pedStructure%inputMap(i) = i
                                        call pedStructure%setAnimalAsGenotyped(i,tmpGeno)
                                        call pedStructure%setAnimalAsHd(i)
                                    endif
                                enddo

                                close(fileUnit)
                            end function initPedigreeGenotypeFiles

                            !---------------------------------------------------------------------------
                            !< @brief Constructor for pedigree class
                            !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals.
                            !< Takes Two arrays as inputs rather than files
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigreeArrays(pedArray, genderArray, nsnps) result(pedStructure)
                                use iso_fortran_env
                                type(PedigreeHolder) :: pedStructure

                                character(len=IDLENGTH), dimension(:,:), intent(in) :: pedArray !< array detailing pedigree of format ped([id, sireId, damId], index)
                                integer, dimension(:) ,optional , intent(in):: genderArray !< gender array corresponding to index in pedArray
                                integer ,optional , intent(in):: nsnps !< number of snps to initialse to 
                                integer(kind=int32) :: tmpSireNum, tmpDamNum
                                integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
                                integer :: tmpAnimalArrayCount
                                integer :: i
                                integer(kind=int64) :: sizeDict
                                logical :: sireFound, damFound

                                pedStructure%nHd = 0
                                pedStructure%nGenotyped = 0
                                pedStructure%nDummys = 0
                                tmpAnimalArrayCount = 0
                                pedStructure%nsnpsPopulation = 0

                                if (present(nsnps)) then
                                    pedStructure%nsnpsPopulation = nsnps
                                endif

                                pedStructure%isSorted = 0
                                sizeDict = size(pedArray)
                                pedStructure%maxPedigreeSize = size(pedArray) + (size(pedArray) * 4)
                                allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
                                pedStructure%pedigreeSize = size(pedArray)
                                pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
                                allocate(tmpAnimalArray(size(pedArray))) !allocate to nIndividuals in case all animals are in incorrect order of generations
                                allocate(pedStructure%inputMap(size(pedArray)))
                                pedStructure%maxGeneration = 0

                                do i=1,size(pedArray)

                                    sireFound = .false.
                                    damFound = .false.

                                    call pedStructure%dictionary%addKey(pedArray(1,i), i)

                                    pedStructure%Pedigree(i) =  Individual(pedArray(1,i),pedArray(2,i),pedArray(3,i), i,nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
                                    pedStructure%Pedigree(i)%originalPosition = i 
                                    pedStructure%inputMap(i) = i
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

                                        if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                                            call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
                                        endif
                                        pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))

                                        if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                                            call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum)) ! add animal to dam list
                                        endif
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


                            !---------------------------------------------------------------------------
                            !< @brief Constructor for pedigree class
                            !< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals.
                            !< Takes Two arrays as inputs rather than files
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function initPedigreeIntArrays(pedArray, genderArray) result(pedStructure)
                                use iso_fortran_env
                                type(PedigreeHolder) :: pedStructure

                                integer, dimension(:,:), intent(in) :: pedArray !< array detailing pedigree of format ped([id, sireId, damId], index)
                                integer, dimension(:) ,optional , intent(in):: genderArray !< gender array corresponding to index in pedArray
                                integer(kind=int32) :: tmpSireNum, tmpDamNum
                                integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
                                integer :: tmpAnimalArrayCount
                                integer :: i
                                integer(kind=int64) :: sizeDict
                                character(len=IDLENGTH) :: tmpID,tmpSireId, tmpDamID

                                logical :: sireFound, damFound

                                pedStructure%nHd = 0
                                pedStructure%nGenotyped = 0
                                pedStructure%nDummys = 0
                                tmpAnimalArrayCount = 0

                                pedStructure%isSorted = 0
                                sizeDict = size(pedArray)
                                pedStructure%maxPedigreeSize = size(pedArray) + (size(pedArray) * 4)
                                allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
                                pedStructure%pedigreeSize = size(pedArray)
                                pedStructure%dictionary = DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
                                allocate(tmpAnimalArray(size(pedArray))) !allocate to nIndividuals in case all animals are in incorrect order of generations
                                allocate(pedStructure%inputMap(size(pedArray)))
                                pedStructure%maxGeneration = 0

                                do i=1,size(pedArray)

                                    sireFound = .false.
                                    damFound = .false.
                                    write(tmpId,*) pedArray(1,i)
                                    write(tmpSireId,*) pedArray(2,i)
                                    write(tmpDamID,*) pedArray(3,i)
                                    call pedStructure%dictionary%addKey(tmpId, i)
                                    pedStructure%Pedigree(i) = Individual(tmpId,tmpSireId,tmpDamID, i,nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
                                    pedStructure%Pedigree(i)%originalPosition = i
                                    pedStructure%inputMap(i) = i
                                    if (pedArray(i,2) /= EMPTY_PARENT) then !check sire is defined in pedigree
                                        tmpSireNum = pedStructure%dictionary%getValue(tmpSireId)
                                        if (tmpSireNum /= DICT_NULL) then
                                            sireFound = .true.
                                        endif
                                    endif

                                    if (pedArray(i,3) /= EMPTY_PARENT) then
                                        tmpDamNum = pedStructure%dictionary%getValue(tmpSireId)
                                        if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
                                            damFound = .true.
                                        endif
                                    endif
                                    if (tmpSireId == EMPTY_PARENT .and. tmpDamID == EMPTY_PARENT) then !if animal is a founder
                                        pedStructure%Pedigree(i)%founder = .true.
                                        call pedStructure%Founders%list_add(pedStructure%Pedigree(i))

                                    else if (sireFound == .false. .or. damFound == .false. ) then
                                        tmpAnimalArrayCount = tmpAnimalArrayCount +1
                                        tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
                                    else ! if sire and dam are both found
                                        pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                                        call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))

                                        if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                                            call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
                                        endif
                                        pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                                        call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))

                                        if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
                                            call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                                            call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum)) ! add animal to sire list
                                        endif
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
                            end function initPedigreeIntArrays

                            !---------------------------------------------------------------------------
                            !< @brief Helper function to avoid code duplication
                            !< required by constructors to determine animals that need offspring info added
                            !< due to it not being initially available in the pedigree
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine addOffspringsAfterReadIn(pedStructure, tmpAnimalArray, tmpAnimalArrayCount)
                                use ConstantModule, only : IDLENGTH,EMPTY_PARENT
                                class(PedigreeHolder), intent(inout) :: pedStructure
                                integer, dimension(:), intent(in) :: tmpAnimalArray !< array containing indexes of tmp animals
                                integer, intent(in) :: tmpAnimalArrayCount !< number of animals actually in tmpAnimalArray
                                logical :: sireFound, damFound
                                integer(kind=int32) :: tmpSireNum, tmpDamNum
                                integer(kind=int32) :: i, tmpCounter
                                character(len=IDLENGTH) :: tmpSire,tmpDam,tmpCounterStr

                                tmpCounter = 0 !< counter for dummy animals

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

                                        ! check that we've not already defined the parent above
                                        if (tmpSireNum /= DICT_NULL .and. .not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer)) then !if sire has been found in hashtable
                                            pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                                            call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))

                                            if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
                                                call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                                                call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
                                            endif
                                            ! check that we've not already defined the parent above
                                        else if (.not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer)) then!if sire is defined but not in the pedigree, create him
                                            ! check if the tmp animal has already been created
                                            tmpSireNum = pedStructure%dictionary%getValue(trim(tmpSire))
                                            if (tmpSireNum == DICT_NULL) then
                                                pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                                                pedStructure%nDummys = pedStructure%nDummys + 1
                                                if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                                                    write(error_unit,*) "ERROR: too many undefined animals"
                                                    stop
                                                endif
                                                pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual(trim(tmpSire),'0','0', pedStructure%pedigreeSize,nsnps=pedStructure%nsnpsPopulation)
                                                call pedStructure%dictionary%addKey(trim(tmpSire), pedStructure%pedigreeSize)
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                                                pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                                                call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))

                                                if (pedStructure%Pedigree(pedStructure%pedigreeSize)%nOffs == 1) then
                                                    call pedStructure%Pedigree(pedStructure%pedigreeSize)%setGender(1) !if its a sire, it should be male
                                                    call pedStructure%sireList%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize)) ! add animal to sire list
                                                endif
                                                call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                                            else
                                                pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
                                                call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                                if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
                                                    call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
                                                    call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
                                                endif
                                            endif
                                        endif
                                        sireFound = .true.
                                    endif

                                    if (tmpDam /= EMPTY_PARENT) then
                                        ! check that we've not already defined the parent above
                                        if (tmpDamNum /= DICT_NULL .and. .not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%damPointer)) then !if dam has been found
                                            pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                                            call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                            call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
                                            if (pedStructure%Pedigree(tmpDamnum)%nOffs == 1) then
                                                call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamnum)) ! add animal to dam list
                                            endif
                                            ! check that we've not already defined the parent above
                                        else if (.not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%damPointer)) then
                                            ! Check for defined animals that have nit been set in pedigree
                                            tmpDamNum = pedStructure%dictionary%getValue(trim(tmpDam))
                                            if (tmpDamNum == DICT_NULL) then !If dummy animal has not already been set in pedigree
                                                pedStructure%nDummys = pedStructure%nDummys + 1
                                                pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                                                if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                                                    write(error_unit,*) "ERROR: too many undefined animals"
                                                    stop
                                                endif
                                                pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual(trim(tmpDam),'0','0', pedStructure%pedigreeSize,nsnps=pedStructure%nsnpsPopulation)
                                                call pedStructure%dictionary%addKey(trim(tmpDam), pedStructure%pedigreeSize)
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                                                pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)
                                                call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                                if (pedStructure%Pedigree(pedStructure%pedigreeSize)%nOffs == 1) then
                                                    call pedStructure%Pedigree(pedStructure%pedigreeSize)%setGender(2) !if its a dam, it should be female
                                                    call pedStructure%damList%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize)) ! add animal to sire list
                                                endif
                                                call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                                            else
                                                pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
                                                call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                                if (pedStructure%Pedigree(pedStructure%pedigreeSize)%nOffs == 1) then
                                                    call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a sire, it should be male, dam female
                                                    call pedStructure%damList%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize)) ! add animal to sire list
                                                endif

                                                call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamnum)) ! add animal to dam list
                                            endif

                                        endif
                                        damFound = .true.
                                    endif


                                    if (.not. damFound .OR. .not. sireFound) then

                                        if (.not. damFound) then
                                            tmpCounter =  tmpCounter + 1
                                            write(tmpCounterStr, '(a,I3.3)') dummyAnimalPrepre,tmpCounter
                                            pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                                            pedStructure%nDummys = pedStructure%nDummys + 1
                                            if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                                                write(error_unit,*) "ERROR: too many undefined animals"
                                                stop
                                            endif
                                            pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual(trim(tmpCounterStr),'0','0', pedStructure%pedigreeSize,nsnps=pedStructure%nsnpsPopulation)
                                            pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.
                                            if (tmpDam == EMPTY_PARENT) then
                                                pedStructure%unknownDummys = pedStructure%unknownDummys+1
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%isUnknownDummy = .true.
                                            endif


                                            pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)

                                            call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                            if (pedStructure%Pedigree(pedStructure%pedigreeSize)%nOffs == 1) then
                                                call pedStructure%Pedigree(pedStructure%pedigreeSize)%setGender(2)
                                                call pedStructure%damList%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize)) ! add animal to sire list
                                            endif

                                            ! add dummy animal to correct lists and dictionaries
                                            call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                                            call pedStructure%dictionary%addKey(tmpCounterStr, pedStructure%pedigreeSize)
                                            pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                                        endif
                                        if (.not. sireFound) then
                                            tmpCounter =  tmpCounter + 1
                                            write(tmpCounterStr, '(a,I3.3)')  dummyAnimalPrepre,tmpCounter
                                            pedStructure%pedigreeSize = pedStructure%pedigreeSize + 1
                                            pedStructure%nDummys = pedStructure%nDummys + 1
                                            if (pedStructure%pedigreeSize > pedStructure%maxPedigreeSize) then
                                                write(error_unit,*) "ERROR: too many undefined animals"
                                                stop
                                            endif
                                            pedStructure%Pedigree(pedStructure%pedigreeSize) =  Individual(trim(tmpCounterStr),'0','0', pedStructure%pedigreeSize,nsnps=pedStructure%nsnpsPopulation)
                                            if (tmpSire == EMPTY_PARENT) then
                                                pedStructure%unknownDummys = pedStructure%unknownDummys+1
                                                pedStructure%Pedigree(pedStructure%pedigreeSize)%isUnknownDummy = .true.
                                            endif
                                            pedStructure%Pedigree(pedStructure%pedigreeSize)%isDummy = .true.

                                            pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(pedStructure%pedigreeSize)

                                            call pedStructure%Pedigree(pedStructure%pedigreeSize)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
                                            if (pedStructure%Pedigree(pedStructure%pedigreeSize)%nOffs == 1) then
                                                call pedStructure%Pedigree(pedStructure%pedigreeSize)%setGender(1)
                                                call pedStructure%sireList%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize)) ! add animal to sire list
                                            endif

                                            ! add animals to correct lists and dictinoaries
                                            call pedStructure%Founders%list_add(pedStructure%Pedigree(pedStructure%pedigreeSize))
                                            pedStructure%Pedigree(pedStructure%pedigreeSize)%founder = .true.
                                            call pedStructure%dictionary%addKey(tmpCounterStr, pedStructure%pedigreeSize)
                                        endif

                                    endif
                                enddo

                            end subroutine addOffspringsAfterReadIn





                            !---------------------------------------------------------------------------
                            !< @brief returns a list of animals that are genotyped, and are classed as founders
                            !< Animals are classed as founders if they have no ancestors that are genotyped in a given number of generations
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getGenotypedFounders(this, numberOfGenerations) result(genotypedFounders)


                                class(pedigreeHolder) :: this
                                integer, intent(in) :: numberOfGenerations
                                type(IndividualLinkedList) :: genotypedFounders
                                integer :: i
                                do i=1, this%pedigreeSize
                                    if (this%pedigree(i)%isGenotypedNonMissing()) then

                                        if (this%pedigree(i)%founder) then
                                            call genotypedFounders%list_add(this%pedigree(i))
                                        else if (.not. this%pedigree(i)%hasGenotypedAnsestors(numberOfGenerations)) then !< this checks if ancestors are not genotyped given a number
                                            call genotypedFounders%list_add(this%pedigree(i))
                                        endif

                                    endif
                                enddo



                            end function getGenotypedFounders



                            subroutine wipeGenotypeAndPhaseInfo(this)
                                class(pedigreeHolder) :: this

                                integer :: i
                                
                                do i=1,this%pedigreeSize

                                    deallocate(this%pedigree(i)%individualPhase)
                                    deallocate(this%pedigree(i)%individualGenotype)
                                enddo

                            end subroutine wipeGenotypeAndPhaseInfo

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

                                enddo
                                if (allocated(this%generations)) then

                                    do i=0, this%maxGeneration
                                        call this%generations(i)%destroyLinkedList
                                    enddo
                                    deallocate(this%generations)
                                endif
                                call this%sireList%destroyLinkedList
                                call this%damList%destroyLinkedList
                                call this%Founders%destroyLinkedListFinal

                                if (allocated(this%inputMap)) then
                                    deallocate(this%inputMap)
                                endif
                                call this%dictionary%destroy !destroy dictionary as we no longer need it
                                if (this%nGenotyped > 0) then
                                    call this%genotypeDictionary%destroy
                                    deallocate(this%genotypeMap)
                                endif

                                if (this%nHd > 0) then
                                    call this%hdDictionary%destroy
                                    deallocate(this%hdMap)
                                endif

                            end subroutine destroyPedigree


                            !---------------------------------------------------------------------------
                            ! DESCRIPTION:
                            !<@brief     Adds genotype information to pedigree from a 2 dimensional array
                            !
                            !<@author     David Wilson, david.wilson@roslin.ed.ac.uk
                            !
                            !<@date       October 25, 2016
                            !---------------------------------------------------------------------------
                            subroutine addGenotypeInformationFromArray(this, array)

                                use AlphaHouseMod, only : countLines
                                implicit none
                                class(PedigreeHolder) :: this
                                integer(kind=1),allocatable,dimension (:,:), intent(in) :: array !< array should be dimensions nanimals, nsnp
                                integer :: i
                                do i=1,this%pedigreeSize-this%nDummys !< assumes dummys are at end, otherwise this will NOT work
                                    call this%setAnimalAsGenotyped(i, array(i,:))
                                enddo

                            end subroutine addGenotypeInformationFromArray


                            !---------------------------------------------------------------------------
                            ! DESCRIPTION:
                            !<@brief     Adds genotype information to pedigree from a file
                            !
                            !<@author     David Wilson, david.wilson@roslin.ed.ac.uk
                            !
                            !<@date       October 25, 2016
                            !--------------------------------------------------------------------------
                            subroutine addGenotypeInformationFromFile(this, genotypeFile, nsnps, nAnnisG, startSnp, endSnp, lockIn)

                                use AlphaHouseMod, only : countLines
                                implicit none
                                class(PedigreeHolder) :: this
                                character(len=*) :: genotypeFile
                                character(len=IDLENGTH) :: tmpID
                                integer,intent(in) :: nsnps
                                integer,intent(in),optional :: nAnnisG
                                integer, intent(in),optional :: startSnp, endSnp
                                logical, intent(in), optional :: lockIn
                                integer(kind=1), allocatable, dimension(:) :: tmpSnpArray 
                                integer :: i, j,fileUnit, nAnnis,tmpIdNum
                                integer :: count, end
                                logical :: lock

                                if (present(lockIn)) then
                                    lock = lockIn
                                else 
                                    lock = .false.

                                endif
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
                                        if ((tmpSnpArray(j)<0).or.(tmpSnpArray(j)>2)) then
                                            tmpSnpArray(j)=9
                                        endif    
                                    enddo

                                    tmpIdNum = this%dictionary%getValue(tmpId)
                                    if (tmpIdNum == DICT_NULL) then
                                        write(error_unit, *) "WARNING: Genotype info for non existing animal here:",trim(tmpId), " file:", trim(genotypeFile), " line:",i
                                        write(error_unit, *) "Animal will be added as a founder to pedigree"

                                        call this%addAnimalAtEndOfPedigree(trim(tmpID), tmpSnpArray)

                                    else

                                        if (present(startSnp)) then
                                            count = 0

                                            if (present(endSnp)) then
                                                end =endSnp
                                            else 
                                                end = nsnps
                                            endif
                                            call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray(startSnp:endSnp),lock)
                                        else if (present(endSnp)) then
                                            count = 0
                                            call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray(1:endsnp),lock)
                                        else 
                                            call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray,lock)
                                        endif
                                    endif
                                enddo
                                write(output_unit,*) "NOTE: Number of Genotyped animals: ",this%nGenotyped

                            end subroutine addGenotypeInformationFromFile

                            !---------------------------------------------------------------------------
                            ! DESCRIPTION:
                            !<@brief     Adds phase information to pedigree from a file
                            !
                            !<@author    Daniel Money, daniel.money@roslin.ed.ac.uk
                            !
                            !<@date       June 19, 2017
                            !--------------------------------------------------------------------------
                            subroutine addPhaseInformationFromFile(this, phaseFile, nsnps, nAnnisG)

                                use AlphaHouseMod, only : countLines
                                implicit none
                                class(PedigreeHolder) :: this
                                character(len=*) :: phaseFile
                                character(len=IDLENGTH) :: tmpID
                                integer,intent(in) :: nsnps
                                integer,intent(in),optional :: nAnnisG
                                integer(kind=1), allocatable, dimension(:) :: tmpSnpArray
                                integer :: i, j, h, fileUnit, nAnnis,tmpIdNum

                                if (present(nAnnisG)) then
                                    nAnnis = nAnnisG
                                else
                                    nAnnis = countLines(phaseFile)
                                endif

                                allocate(tmpSnpArray(nsnps))
                                open(newUnit=fileUnit, file=phaseFile, status="old")
                                do i=1, nAnnis
                                    do h = 1, 2
                                        read (fileUnit,*) tmpId,tmpSnpArray(:)
                                        do j=1,nsnps
                                            if ((tmpSnpArray(j)<0).or.(tmpSnpArray(j)>1)) tmpSnpArray(j)=9
                                        enddo
                                        tmpIdNum = this%dictionary%getValue(tmpId)
                                        if (tmpIdNum == DICT_NULL) then
                                            write(error_unit, *) "WARNING: Phase info for non existing animal here:",trim(tmpId), " file:", trim(phaseFile), " line:",i
                                        else
                                            call this%pedigree(tmpIdNum)%setPhaseArray(h, tmpSnpArray)
                                        endif
                                    end do
                                enddo

                            end subroutine addPhaseInformationFromFile


                            subroutine addSequenceFromFile(this, seqFile, nsnps, nAnisGIn,maximumReads, startSnp, endSnp)

                                use AlphaHouseMod, only : countLines
                                use ConstantModule, only : IDLENGTH,DICT_NULL
                                implicit none
                                class(PedigreeHolder) :: this
                                character(len=*) :: seqFile
                                integer,intent(in) :: nsnps
                                integer, intent(in), optional :: maximumReads
                                integer,intent(in),optional :: nAnisGIn
                                integer,intent(in),optional :: startSnp, endSnp 
                                integer :: nanisG, end
                                ! type(Pedigreeholder), intent(inout) :: genotype
                                integer(KIND=1), allocatable, dimension(:) :: tmp
                                integer,allocatable, dimension(:) :: ref, alt
                                integer :: unit, tmpID,i,j
                                character(len=IDLENGTH) :: seqid !placeholder variables

                                if (.not. Present(nAnisGIn)) then
                                    NanisG = countLines(seqFile)/2
                                else 
                                    nanisG = nAnisGIn
                                endif

                                open(newunit=unit,FILE=trim(seqFile),STATUS="old") !INPUT FILE
                                allocate(ref(nsnps))
                                allocate(alt(nsnps))
                                ! tmp = 9
                                print *, "Number of animals in seq file", NanisG
                                do i=1,nAnisG
                                    read (unit,*) seqid, ref(:)
                                    read (unit,*) seqid, alt(:)

                                    tmpID = this%dictionary%getValue(seqid)
                                    if (present(maximumReads)) then
                                        do j=1,nsnps
                                            if (ref(j)>=maximumReads) ref(j)=maximumReads-1
                                            if (alt(j)>=maximumReads) alt(j)=maximumReads-1
                                        enddo
                                    endif

                                    if (tmpID /= DICT_NULL) then

                                        if (present(startSnp)) then
                                            if (present(endSnp)) then
                                                end = endsnp
                                            else
                                                end = nsnps
                                            endif
                                            call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref(startSnp:end),alt(startSnp:end))
                                        else if (present(endSnp)) then
                                            call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref(1:endsnp),alt(1:endSnp))
                                        else
                                            call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref,alt)
                                        endif
                                    endif
                                end do

                                close(unit)
                            end subroutine addSequenceFromFile



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
                                    call this%setOffspringGeneration(tmpIndNode%item)
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
                                integer :: unit, i,h
                                type(IndividualLinkedListNode), pointer :: tmpIndNode
                                character(len=:), allocatable :: fmt

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
                                        write(unit, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId, tmpIndNode%item%generation
                                        ! write(*, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
                                        tmpIndNode => tmpIndNode%next
                                    end do
                                enddo
                                if (present(file)) then
                                    close(unit)
                                endif
                            end subroutine outputSortedPedigree


                            !---------------------------------------------------------------------------
                            !< @brief Output pedigree to stdout in the format recodedID, recodedSireId, recodedDamId, originalId
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine outputSortedPedigreeInAlphaImputeFormat(this, file)
                                use iso_fortran_env, only : output_unit
                                class(PedigreeHolder) :: this
                                character(len=*), intent(in), optional :: file !< output path for sorted pedigree
                                integer :: unit, i,h, sortCounter
                                type(IndividualLinkedListNode), pointer :: tmpIndNode
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




                            !---------------------------------------------------------------------------
                            !< @brief Sorts pedigree, and overwrites all fields to new values
                            !< @details effectively, does a deep copy to sort pedigree based on generation, but puts dummys at bottom
                            !< If value is given for unknownDummysAtEnd, then only unknown dummys will be put at the end
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine sortPedigreeAndOverwrite(this, unknownDummysAtEnd)
                                use iso_fortran_env, only : output_unit, int64
                                class(PedigreeHolder) :: this
                                integer :: i,h, pedCounter, tmpId,tmpGenotypeMapIndex
                                integer(kind=int64) :: sizeDict
                                type(IndividualLinkedListNode), pointer :: tmpIndNode
                                type(Individual), pointer, dimension(:) :: newPed
                                type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList
                                type(IndividualLinkedList) :: dummyList
                                integer, intent(in) , optional :: unknownDummysAtEnd !< if this option is specified, then only unknown dummies are put at end

                                if (.not. allocated(this%generations)) then
                                    call this%setPedigreeGenerationsAndBuildArrays
                                endif
                                pedCounter = 0

                                call this%dictionary%destroy()
                                call this%founders%destroyLinkedList()
                                call this%sireList%destroyLinkedList()
                                call this%damList%destroyLinkedList()
                                sizeDict  =this%pedigreeSize
                                allocate(newPed(this%maxPedigreeSize))
                                this%dictionary = DictStructure(sizeDict)
                                allocate(newGenerationList(0:this%maxGeneration))
                                do i=0, this%maxGeneration
                                    tmpIndNode => this%generations(i)%first
                                    do h=1, this%generations(i)%length

                                        if (present(unknownDummysAtEnd)) then
                                            if (tmpIndNode%item%isUnknownDummy) then
                                                call dummyList%list_add(tmpIndNode%item)
                                                tmpIndNode => tmpIndNode%next
                                                cycle
                                            endif
                                        else  
                                            if (tmpIndNode%item%isDummy) then
                                                call dummyList%list_add(tmpIndNode%item)
                                                tmpIndNode => tmpIndNode%next
                                                cycle
                                            endif
                                        endif

                                        pedCounter = pedCounter +1
                                        ! update dictionary index
                                        call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)

                                        !  update genotype map
                                        if (this%nGenotyped > 0) then
                                            tmpGenotypeMapIndex = this%genotypeDictionary%getValue(tmpIndNode%item%originalID)
                                            if (tmpGenotypeMapIndex /= DICT_NULL) then
                                                this%genotypeMap(tmpGenotypeMapIndex) = pedCounter
                                            endif
                                        endif
                                        ! Update hd map
                                        if (this%nHd > 0) then
                                            tmpGenotypeMapIndex = this%hdDictionary%getValue(tmpIndNode%item%originalID)
                                            if (tmpGenotypeMapIndex /= DICT_NULL) then
                                                this%hdMap(tmpGenotypeMapIndex) = pedCounter
                                            endif
                                        endif
                                        ! Set the location to the pedigree to the new value
                                        newPed(pedCounter) = tmpIndNode%item

                                        ! take the original id, and update it
                                        if(.not. newPed(pedCounter)%isDummy) then
                                            this%inputMap(newPed(pedCounter)%originalPosition) = pedCounter
                                        endif
                                        newPed(pedCounter)%id = pedCounter
                                        call newPed(pedCounter)%resetOffspringInformation ! reset offsprings
                                        if (associated(newPed(pedCounter)%sirePointer)) then
                                            tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
                                            if (tmpID /= DICT_NULL) then
                                                call newPed(tmpId)%addOffspring(newPed(pedCounter))
                                                newPed(pedCounter)%sirePointer=> newPed(tmpId)
                                                if (newPed(tmpId)%nOffs == 1) then
                                                    call this%sireList%list_add(newPed(tmpId))
                                                endif
                                            endif
                                        endif

                                        if (associated(newPed(pedCounter)%damPointer)) then
                                            tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
                                            if (tmpID /= DICT_NULL) then
                                                call newPed(tmpId)%addOffspring(newPed(pedCounter))
                                                newPed(pedCounter)%damPointer=> newPed(tmpId)
                                                if (newPed(tmpId)%nOffs == 1) then
                                                    call this%damList%list_add(newPed(tmpId))
                                                endif
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
                                    call newPed(pedCounter)%resetOffspringInformation() ! reset offsprings

                                    do h=1, tmpIndNode%item%nOffs
                                        tmpId =  this%dictionary%getValue(tmpIndNode%item%offsprings(h)%p%originalID)
                                        call newPed(pedCounter)%addOffspring(newPed(tmpID))
                                        if(associated(tmpIndNode%item%offsprings(h)%p%sirePointer, tmpIndNode%item)) then
                                            newPed(tmpId)%sirePointer => newPed(pedCounter)
                                            if (newPed(pedCounter)%nOffs == 1) then
                                                call this%sireList%list_add(newPed(pedCounter))
                                            endif
                                        else
                                            newPed(tmpId)%damPointer => newPed(pedCounter)
                                            if (newPed(pedCounter)%nOffs == 1) then
                                                call this%damList%list_add(newPed(pedCounter))
                                            endif
                                        endif
                                    enddo
                                    call  this%founders%list_add(newPed(pedCounter))

                                    call newGenerationList(0)%list_add(newPed(pedCounter))
                                    tmpIndNode => tmpIndNode%next

                                enddo

                                call dummyList%destroyLinkedList()

                                do i=1,this%pedigreeSize
                                    call this%Pedigree(i)%destroyIndividual
                                enddo
                                deallocate(this%pedigree)

                                this%pedigree => newPed
                                do i = 0, this%maxGeneration
                                    call this%generations(i)%destroyLinkedList
                                enddo
                                this%generations = newGenerationList

                                !
                                if (present(unknownDummysAtEnd)) then
                                    this%isSorted = 2
                                else 
                                    this%isSorted = 1
                                endif
                            end subroutine sortPedigreeAndOverwrite


                            !---------------------------------------------------------------------------
                            !< @brief Output pedigree to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine sortPedigreeAndOverwriteWithDummyAtTheTop(this)
                                use iso_fortran_env, only : output_unit, int64
                                class(PedigreeHolder) :: this
                                integer :: i,h, pedCounter, tmpId,tmpGenotypeMapIndex
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
                                call this%sireList%destroyLinkedList()
                                call this%damList%destroyLinkedList()
                                sizeDict  =this%pedigreeSize
                                allocate(newPed(this%maxPedigreeSize))
                                this%dictionary = DictStructure(sizeDict)
                                allocate(newGenerationList(0:this%maxGeneration))
                                do i=0, this%maxGeneration
                                    tmpIndNode => this%generations(i)%first
                                    do h=1, this%generations(i)%length
                                        pedCounter = pedCounter +1
                                        call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)

                                        !  update genotype map
                                        if (this%nGenotyped > 0) then
                                            tmpGenotypeMapIndex = this%genotypeDictionary%getValue(tmpIndNode%item%originalID)
                                            if (tmpGenotypeMapIndex /= DICT_NULL) then
                                                this%genotypeMap(tmpGenotypeMapIndex) = pedCounter
                                            endif
                                        endif
                                        ! Update hd map
                                        if (this%nHd > 0) then
                                            tmpGenotypeMapIndex = this%hdDictionary%getValue(tmpIndNode%item%originalID)
                                            if (tmpGenotypeMapIndex /= DICT_NULL) then
                                                this%hdMap(tmpGenotypeMapIndex) = pedCounter
                                            endif
                                        endif

                                        newPed(pedCounter) = tmpIndNode%item
                                        if (.not. newPed(pedCounter)%isDummy) then
                                            ! take the original id, and update it - we don't want dummies in this list
                                            this%inputMap(newPed(pedCounter)%originalPosition) = pedCounter
                                        endif
                                        newPed(pedCounter)%id = pedCounter
                                        call newPed(pedCounter)%resetOffspringInformation ! reset offsprings

                                        if (associated(newPed(pedCounter)%sirePointer)) then
                                            tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
                                            if (tmpID /= DICT_NULL) then
                                                call newPed(tmpId)%addOffspring(newPed(pedCounter))
                                                newPed(pedCounter)%sirePointer=> newPed(tmpId)
                                                if (newPed(tmpId)%nOffs == 1) then
                                                    call this%sireList%list_add(newPed(tmpId))
                                                endif

                                            endif
                                        endif

                                        if (associated(newPed(pedCounter)%damPointer)) then
                                            tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
                                            if (tmpID /= DICT_NULL) then
                                                call newPed(tmpId)%addOffspring(newPed(pedCounter))
                                                newPed(pedCounter)%damPointer=> newPed(tmpId)
                                                if (newPed(tmpId)%nOffs == 1) then
                                                    call this%damList%list_add(newPed(tmpId))
                                                endif
                                            endif
                                        endif
                                        if (i ==0 ) then !if object is afounder add to founder array
                                            call  this%founders%list_add(newPed(pedCounter))
                                        endif

                                        call newGenerationList(i)%list_add(newPed(pedCounter))
                                        tmpIndNode => tmpIndNode%next
                                    end do
                                enddo
                                do i = 0, this%maxGeneration
                                    call this%generations(i)%destroyLinkedList
                                enddo

                                do i=1,this%pedigreeSize
                                    call this%Pedigree(i)%destroyIndividual
                                enddo
                                deallocate(this%pedigree)

                                this%pedigree => newPed
                                this%generations = newGenerationList
                                this%isSorted = 3
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
                                    print *, this%pedigree(i)%originalId, this%pedigree(i)%getIntegerVectorOfRecodedIds()
                                enddo
                            end subroutine printPedigree


                            !---------------------------------------------------------------------------
                            !< @brief Output pedigree to stdout in the format originalID,sireId,damId
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine printPedigreeOriginalFormat(this, filePath)
                                class(PedigreeHolder) :: this
                                character(len=*), optional :: filePath
                                integer ::i,unit

                                if(present(filepath)) then
                                    open(newunit=unit, file=filePath, status="unknown")
                                else
                                    unit = output_unit
                                endif 
                                do i= 1, this%pedigreeSize
                                    write(unit,*)  this%pedigree(i)%originalId, this%pedigree(i)%sireId, this%pedigree(i)%damId
                                enddo
                            end subroutine printPedigreeOriginalFormat


                            !---------------------------------------------------------------------------
                            !< @brief Output genders to file
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine writeOutGenders(this, filepath)
                                class(PedigreeHolder) :: this
                                character(len=*), intent(in) :: filepath
                                integer ::i, unit

                                open(newunit= unit, file= filepath, status="new")
                                do i= 1, this%pedigreeSize
                                    write(unit,*) this%pedigree(i)%originalId, this%pedigree(i)%gender
                                enddo

                                close(unit)
                            end subroutine writeOutGenders


                            !---------------------------------------------------------------------------
                            !< @brief read in gender information
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine readInGenders(this, filepath)
                                class(PedigreeHolder) :: this
                                character(len=*), intent(in) :: filepath

                                integer ::i, tmpGender,tmp,unit
                                character(len=IDLENGTH) :: tmpId

                                open(newunit= unit, file= filepath, status="old")
                                do i= 1, this%pedigreeSize
                                    read(unit,*) tmpId, tmpGender
                                    tmp = this%dictionary%getValue(tmpId)
                                    if (tmp/=DICT_NULL) then
                                        this%pedigree(tmp)%gender = tmpGender
                                    else
                                        write(error_unit,*) "WARNING: Gender info exists for animal not in the pedigree" 
                                    endif
                                enddo
                            end subroutine readInGenders



                            !---------------------------------------------------------------------------
                            !< @brief Output  of animals that are genotyped
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine writeOutGenotypes(this, filename)
                                class(PedigreeHolder) :: this
                                character(*), intent(in) :: filename
                                integer ::i, fileUnit

                                open(newUnit=fileUnit,file=filename,status="unknown")
                                do i= 1, this%nGenotyped
                                    write(fileUnit,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)')  this%pedigree(this%genotypeMap(i))%originalId, this%pedigree(this%genotypeMap(i))%individualGenotype%toIntegerArray()
                                enddo
                                close(fileUnit)
                            end subroutine writeOutGenotypes


                            !---------------------------------------------------------------------------
                            !< @brief Output genotypes to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
                            !< for all animals not just the ones that are genotyped\
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine writeOutGenotypesAll(this, filename)
                                class(PedigreeHolder) :: this
                                character(*), intent(in) :: filename
                                integer ::i, fileUnit

                                open(newUnit=fileUnit,file=filename,status="unknown")
                                do i= 1, this%pedigreeSize
                                    write(fileUnit,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)')  this%pedigree(i)%originalId, this%pedigree(i)%individualGenotype%toIntegerArray()
                                enddo
                                close(fileUnit)
                            end subroutine writeOutGenotypesAll



                            !---------------------------------------------------------------------------
                            !< @brief Outputs phase to file
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine WriteoutPhase(this, filename)
                                class(PedigreeHolder) :: this
                                character(*), intent(in) :: filename
                                integer ::i, fileUnit

                                open(newUnit=fileUnit,file=filename,status="unknown")
                                do i= 1, this%pedigreeSize
                                    write(fileUnit,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)')  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(1)%toIntegerArray()
                                    write(fileUnit,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)')  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(2)%toIntegerArray()
                                enddo
                                close(fileUnit)
                            end subroutine WriteoutPhase



                            !---------------------------------------------------------------------------
                            !< @brief Sets generation of an individual and his children recursively
                            !< @details makes assumption that both parents also exist, and that is how generation is got
                            !< both parents generation has to be set for this to work
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    Febuary 17, 2016
                            !< @param[in] generation (integer)
                            !< @param[in] pointer to an individual
                            !---------------------------------------------------------------------------
                            recursive subroutine setOffspringGeneration(this, indiv)
                                type(Individual),pointer, intent(inout) :: indiv
                                class(pedigreeHolder):: this


                                integer :: i
                                if (indiv%generation /= NOGENERATIONVALUE) then !< animal has already been set so return
                                    return
                                endif
                                if (.not. indiv%founder) then
                                    ! if the generation of both parents has been set, add one to the greater one
                                    if (indiv%sirePointer%generation /= NOGENERATIONVALUE .and. indiv%damPointer%generation /= NOGENERATIONVALUE) then

                                      indiv%generation = max(indiv%sirePointer%generation, indiv%damPointer%generation)+1
                                     !   if (indiv%sirePointer%generation > indiv%damPointer%generation) then
                                     !       indiv%generation = indiv%sirePointer%generation + 1
                                     !   else
                                     !       indiv%generation = indiv%damPointer%generation + 1
                                     !   endif

                                    else !otherwise, both parents have not been set so return, as animal will get checked later
                                        return
                                    endif

                                else
                                    indiv%generation = 0
                                endif

                                call this%generations(indiv%generation)%list_add(indiv)
                                if(indiv%generation > this%maxGeneration) then
                                    this%maxGeneration = indiv%generation
                                endif

                                if ( indiv%nOffs /= 0) then
                                    do i=1,indiv%nOffs

                                        call this%setOffspringGeneration(indiv%OffSprings(i)%p)
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
                                class(pedigreeHolder), intent(in) :: this      !< object to operate on
                                type(recodedPedigreeArray), intent(out) :: recPed !< @return recoded pedigree array
                                integer :: counter,i,h
                                type(IndividualLinkedListNode), pointer :: tmpIndNode

                                call RecPed%init(n=this%pedigreeSize)

                                call this%sortPedigreeAndOverwriteWithDummyAtTheTop

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


                            !---------------------------------------------------------------------------
                            !< @brief returns array of genotypes for all animals 
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getAllGenotypesAtPosition(this, position) result(res)
                                use constantModule, only : MISSINGPHASECODE
                                class(pedigreeHolder) :: this
                                integer, intent(in) :: position
                                integer(KIND=1), allocatable, dimension(:) :: res
                                integer :: counter, i
                                allocate(res(this%nGenotyped))
                                res = MISSINGPHASECODE
                                counter = 0

                                do i=1, this%nGenotyped

                                    counter = counter +1
                                    res(counter) = this%pedigree(this%genotypeMap(i))%individualGenotype%getGenotype(position)
                                    if (res(counter) /= 0 .and. res(counter) /= 1 .and. res(counter) /= 2 .and. res(counter) /= MISSINGGENOTYPECODE) then
                                        res(counter) = MISSINGGENOTYPECODE
                                    endif

                                enddo

                            end function getAllGenotypesAtPosition




                            !---------------------------------------------------------------------------
                            !< @brief returns list of mates and offspring for those mate pairs for given pedigree
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getPhaseAtPosition(this, position, allele) result(res)
                                use constantModule, only : MISSINGPHASECODE
                                class(pedigreeHolder) :: this
                                integer, intent(in) :: position
                                integer, intent(in) :: allele
                                integer(KIND=1), allocatable, dimension(:) :: res
                                integer :: counter, i
                                allocate(res(this%nGenotyped))
                                res = MISSINGPHASECODE
                                counter = 0

                                do i=1, this%nGenotyped

                                    counter = counter +1
                                    res(counter) = this%pedigree(this%genotypeMap(i))%individualPhase(allele)%getPhase(position)
                                    if (res(counter) /= 0 .and. res(counter) /= 1 .and. res(counter) /= 2 .and. res(counter) /= MISSINGPHASECODE) then
                                        res(counter) = MISSINGPHASECODE
                                    endif

                                enddo

                            end function getPhaseAtPosition


                            function getPhaseAtPositionUngenotypedAnimals(this, position, allele) result(res)
                                use constantModule, only : MISSINGPHASECODE
                                class(pedigreeHolder) :: this
                                integer, intent(in) :: position
                                integer, intent(in) :: allele
                                integer(KIND=1), allocatable, dimension(:) :: res
                                integer :: counter, i
                                allocate(res(this%pedigreeSize))
                                res = MISSINGPHASECODE
                                counter = 0

                                do i=1, this%nGenotyped

                                    counter = counter +1
                                    res(this%genotypeMap(i)) = this%pedigree(this%genotypeMap(i))%individualPhase(allele)%getPhase(position)
                                    if (res(this%genotypeMap(i)) /= 0 .and. res(this%genotypeMap(i)) /= 1 .and. res(this%genotypeMap(i)) /= 2 .and. res(this%genotypeMap(i)) /= MISSINGPHASECODE) then
                                        res(this%genotypeMap(i)) = MISSINGPHASECODE
                                    endif

                                enddo

                            end function getPhaseAtPositionUngenotypedAnimals



                            !---------------------------------------------------------------------------
                            !< @brief returns list of mates and offspring for those mate pairs for given pedigree
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getAllGenotypesAtPositionWithUngenotypedAnimals(this, position) result(res)
                                use constantModule, only : MISSINGPHASECODE
                                class(pedigreeHolder) :: this
                                integer, intent(in) :: position
                                integer(KIND=1), allocatable, dimension(:) :: res
                                integer :: i
                                allocate(res(this%pedigreeSize))
                                res = MISSINGPHASECODE

                                do i=1, this%pedigreeSize

                                    res(i) = this%pedigree(i)%individualGenotype%getGenotype(position)
                                    if (res(i) /= 0 .and. res(i) /= 1 .and. res(i) /= 2 .and. res(i) /= MISSINGPHASECODE) then
                                        res(i) = MISSINGPHASECODE
                                    endif

                                enddo

                            end function getAllGenotypesAtPositionWithUngenotypedAnimals

                            !---------------------------------------------------------------------------
                            !< @brief returns array of what percentages an animal has been genotyped
                            !<
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getGenotypePercentage(this) result(res)
                                use constantModule, only : MISSINGPHASECODE
                                class(pedigreeHolder) :: this
                                real(KIND=real64), allocatable, dimension(:) :: res
                                integer(kind=1), allocatable, dimension(:) :: indGenotypeArray
                                logical, dimension(:), allocatable :: genotypedAtMarker
                                integer :: i

                                allocate(res(this%pedigreeSize))
                                res = 0

                                do i=1, this%pedigreeSize

                                    if (this%pedigree(i)%isGenotyped()) then

                                        ! print *, "genotyped", indGenotypeArray
                                        indGenotypeArray = this%pedigree(i)%individualGenotype%toIntegerArray()
                                        genotypedAtMarker = ((indGenotypeArray == 0 .or. indGenotypeArray == 1) .or. indGenotypeArray == 2)
                                        res(i) = count(genotypedAtMarker)*1d0/size(genotypedAtMarker)
                                    endif
                                enddo

                            end function getGenotypePercentage



                            !---------------------------------------------------------------------------
                            !< @brief returns array of genotype information as is used by alphaimpute in format (0:nGenotyped, nSnp)
                            !<
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getGenotypesAsArray(this) result(res)

                                class(pedigreeHolder) :: this
                                integer(kind=1) ,dimension(:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
                                integer :: i


                                allocate(res(0:this%pedigreeSize, this%pedigree(this%genotypeMap(1))%individualGenotype%length))
                                res = 9
                                do i=1, this%nGenotyped
                                    res(this%genotypeMap(i),:) = this%pedigree(this%genotypeMap(i))%individualGenotype%toIntegerArray()
                                enddo

                            end function getGenotypesAsArray


                            function getPhaseAsArray(this) result(res)

                                class(pedigreeHolder) :: this
                                integer(kind=1) ,dimension(:,:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
                                integer :: i


                                allocate(res(0:this%pedigreeSize, this%pedigree(this%genotypeMap(1))%individualGenotype%length,2))
                                res = 9
                                do i=1, this%nGenotyped
                                    res(this%genotypeMap(i),:,1) = this%pedigree(this%genotypeMap(i))%individualPhase(1)%toIntegerArray()
                                    res(this%genotypeMap(i),:,2) = this%pedigree(this%genotypeMap(i))%individualPhase(2)%toIntegerArray()
                                enddo

                            end function getPhaseAsArray


                            !---------------------------------------------------------------------------
                            !< @brief returns array of genotype information as is used by alphaimpute in format (0:pedSized, nSnp)
                            !< This takes the genotype info even if an animal is not genotyped
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getGenotypesAsArrayWitHMissing(this) result(res)

                                class(pedigreeHolder) :: this
                                integer(kind=1) ,dimension(:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
                                integer :: i


                                allocate(res(0:this%pedigreeSize, this%pedigree(this%genotypeMap(1))%individualGenotype%length))
                                do i=1, this%pedigreeSize
                                    res(i,:) = this%pedigree(i)%individualGenotype%toIntegerArray()
                                enddo

                            end function getGenotypesAsArrayWitHMissing

                            !---------------------------------------------------------------------------
                            !< @brief returns integer value of number of missing genotypes accross all genotypes
                            !< @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !< @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getNumGenotypesMissing(this) result(count)

                                class(pedigreeHolder) :: this
                                integer :: count,i

                                count = 0
                                !$omp parallel do reduction(+:count)
                                do i=1, this%nGenotyped
                                    count = count + this%pedigree(this%genotypeMap(i))%individualGenotype%numMissing()
                                enddo
                                !$omp end parallel do
                            end function getNumGenotypesMissing



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
                                                endif
                                            endif
                                        enddo
                                        tmpIndNode => tmpIndNode%next

                                    enddo
                                enddo

                                call dictionary%destroy()

                            end subroutine getMatePairsAndOffspring


                            !---------------------------------------------------------------------------
                            !<@brief Sets the individual to be genotyped.
                            !<If geno array is not given, animal will still be set to genotyped. It is up to the callee
                            !<if the animal has enough snps set to actually genotyped
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine setAnimalAsGenotyped(this, individualIndex, geno, lockIn)

                                class(pedigreeHolder) :: this
                                integer, intent(in) :: individualIndex !< index of animal to get genotyped
                                integer(KIND=1), dimension(:),optional, intent(in) :: geno !< One dimensional array of genotype information
                                logical, intent(in), optional :: lockIn
                                logical :: lock

                                if (present(lockIn)) then
                                    lock = lockIn
                                else 
                                    lock = .false.
                                endif



                                if (this%nGenotyped == 0) then
                                    this%genotypeDictionary = DictStructure()
                                    allocate(this%genotypeMap(this%pedigreeSize))
                                    this%genotypeMap = 0

                                else if (this%nGenotyped > this%pedigreeSize) then
                                    ! Following error should never appear
                                    write(error_unit,*) "Error: animals being genotyped that are bigger than ped structure size!"
                                else if (this%genotypeDictionary%getValue(this%pedigree(individualIndex)%originalID) /= DICT_NULL) then
                                    ! if animal has already been genotyped, overwrite array, but don't increment
                                    if (present(geno)) then
                                        call this%pedigree(individualIndex)%setGenotypeArray(geno,lock)
                                    endif
                                    return
                                endif

                                this%nGenotyped = this%nGenotyped+1
                                call this%genotypeDictionary%addKey(this%pedigree(individualIndex)%originalID, this%nGenotyped)
                                if (present(geno)) then
                                    call this%pedigree(individualIndex)%setGenotypeArray(geno,lock)
                                endif
                                this%genotypeMap(this%nGenotyped) = individualIndex

                            end subroutine setAnimalAsGenotyped



                            !---------------------------------------------------------------------------
                            !<@brief Sets the individual to be genotyped.
                            !<If geno array is not given, animal will still be set to genotyped. It is up to the callee
                            !<if the animal has enough snps set to actually genotyped
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine setAnimalAsGenotypedSequence(this, individualIndex, geno, referAllele, alterAllele)

                                class(pedigreeHolder) :: this
                                integer, intent(in) :: individualIndex !< index of animal to get genotyped
                                integer(KIND=1), dimension(:),optional, intent(in) :: geno !< One dimensional array of genotype information
                                integer, dimension(:), intent(in) :: referAllele, alterAllele
                                if (this%nGenotyped == 0) then
                                    this%genotypeDictionary = DictStructure()
                                    allocate(this%genotypeMap(this%pedigreeSize))

                                else if (this%nGenotyped > this%pedigreeSize) then
                                    ! Following error should never appear
                                    write(error_unit,*) "Error: animals being genotyped that are bigger than ped structure size!"
                                else if (this%genotypeDictionary%getValue(this%pedigree(individualIndex)%originalID) /= DICT_NULL) then
                                    ! if animal has already been genotyped, overwrite array, but don't increment
                                    if (present(geno)) then
                                        call this%pedigree(individualIndex)%setGenotypeArray(geno)
                                    endif
                                    return
                                endif

                                this%nGenotyped = this%nGenotyped+1
                                call this%genotypeDictionary%addKey(this%pedigree(individualIndex)%originalID, this%nGenotyped)
                                if (present(geno)) then
                                    call this%pedigree(individualIndex)%setGenotypeArray(geno)
                                endif

                                this%pedigree(individualIndex)%referAllele = referAllele
                                this%pedigree(individualIndex)%alterAllele = alterAllele

                                this%genotypeMap(this%nGenotyped) = individualIndex

                            end subroutine setAnimalAsGenotypedSequence


                            !---------------------------------------------------------------------------
                            !<@brief  Converts the sequence data counts to a 3D array (similar to phase) of size of pedigree
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function convertSequenceDataToArray(this) result(res)

                                class(PedigreeHolder), intent(in) :: this
                                integer, dimension(:,:,:),allocatable :: res
                                integer :: i

                                do i =1, this%pedigreeSize

                                    if (.not. allocated(this%pedigree(i)%referAllele)) cycle

                                    if(.not. allocated(res)) then
                                        allocate(res(this%pedigreesize,size(this%pedigree(i)%referAllele),2))
                                        res = 0
                                    endif

                                    res(i,:,1) = this%pedigree(i)%referAllele
                                    res(i,:,2) = this%pedigree(i)%alterAllele
                                enddo

                            end function convertSequenceDataToArray

                            function getSequenceAsArrayWithMissing(this, index) result(res)

                                class(PedigreeHolder), intent(in) :: this
                                integer, dimension(:,:),allocatable :: res
                                integer :: i, index

                                if(.not. allocated(res)) then
                                    allocate(res(this%pedigreesize,2))
                                endif

                                res = 0
                                do i =1, this%pedigreeSize

                                    if (.not. allocated(this%pedigree(i)%referAllele)) cycle

                                    res(i,1) = this%pedigree(i)%referAllele(index)
                                    res(i,2) = this%pedigree(i)%alterAllele(index)
                                enddo

                            end function getSequenceAsArrayWithMissing

                            !---------------------------------------------------------------------------
                            !<@brief Returns either the individuals id, the sires id or dams id based on
                            !<which index is passed.

                            !<THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function getSireDamGenotypeIDByIndex(this,ind, index) result(v)
                                use iso_fortran_env, only : ERROR_UNIT
                                class(PedigreeHolder), intent(in) :: this
                                type(Individual), intent(in) :: ind
                                character(len=IDLENGTH) :: tmp
                                integer, intent(in) :: index !< index of geno index to return (1 for this, 2 for sire, 3 for dam)
                                integer:: v

                                v = 0
                                select case (index)
                                case(1)
                                    tmp = ind%originalId
                                    v = this%genotypeDictionary%getValue(tmp)
                                    if (v == DICT_NULL) then
                                        v = 0
                                    endif
                                case(2)
                                    if (associated(ind%sirePointer)) then
                                        tmp = ind%sirePointer%originalId
                                        v = this%genotypeDictionary%getValue(tmp)
                                        if (v == DICT_NULL) then
                                            v = 0
                                        endif
                                    endif
                                case(3)
                                    if (associated(ind%damPointer)) then
                                        tmp = ind%damPointer%originalId
                                        v = this%genotypeDictionary%getValue(tmp)
                                        if (v == DICT_NULL) then
                                            v = 0
                                        endif
                                    endif
                                case default
                                    write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
                                end select
                                return
                            end function getSireDamGenotypeIDByIndex


                            !---------------------------------------------------------------------------
                            !<@brief Sets the individual to be genotyped at high density.
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine setAnimalAsHD(this, indId)
                                use iso_fortran_env
                                class(PedigreeHolder) :: this
                                integer, intent(in) :: indId

                                ! if index not in pedigree return.
                                if (indId > this%pedigreeSize) then
                                    write(error_unit, *) "warning - setAnimalAsHD was given an index that was out of range"
                                    return
                                endif
                                if (.not. this%pedigree(indId)%genotyped) then
                                    write(error_unit, *) "warning - setAnimalAsHD was given an index of animal that was not genotyped"
                                    write(error_unit, *) "animal has ID:", trim(this%pedigree(indId)%originalID), " and recoded ID:", indid
                                endif
                                if (this%nHd == 0) then
                                    this%hdDictionary = DictStructure()
                                    allocate(this%hdMap(this%pedigreeSize))
                                    this%hdMap = 0
                                endif

                                if (this%hdDictionary%getValue(this%pedigree(indId)%originalId) ==DICT_NULL) then
                                    this%nHd = this%nHd + 1
                                    this%pedigree(indId)%hd = .true.

                                    this%hdMap(this%nHd) = indId
                                    call this%hdDictionary%addKey(this%pedigree(indId)%originalId, this%nHd)
                                else
                                    this%pedigree(indId)%hd = .true.
                                endif

                            end subroutine setAnimalAsHD



                            !---------------------------------------------------------------------------
                            !<@brief Returns either the individuals id in hd index, the sires id or dams id based on
                            !<which index is passed.
                            !<THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            ! PARAMETERS:
                            !<@param[in] index - the index
                            !<@return hdIndex of animal based on index
                            !---------------------------------------------------------------------------
                            function getSireDamHDIDByIndex(this,ind, index) result(v)
                                use iso_fortran_env, only : ERROR_UNIT
                                class(PedigreeHolder), intent(in) :: this
                                type(Individual), intent(in) :: ind
                                character(len=IDLENGTH) :: tmp
                                integer, intent(in) :: index !< index of hd index to return (1 for this, 2 for sire, 3 for dam)
                                integer:: v

                                v = 0
                                select case (index)
                                case(1)
                                    tmp = ind%originalId
                                    v = this%hdDictionary%getValue(tmp)
                                    if (v == DICT_NULL) then
                                        v = 0
                                    endif
                                case(2)
                                    if (associated(ind%sirePointer)) then
                                        tmp = ind%sirePointer%originalId
                                        v = this%hdDictionary%getValue(tmp)
                                        if (v == DICT_NULL) then
                                            v = 0
                                        endif
                                    endif
                                case(3)
                                    if (associated(ind%damPointer)) then
                                        tmp = ind%damPointer%originalId
                                        v = this%hdDictionary%getValue(tmp)
                                        if (v == DICT_NULL) then
                                            v = 0
                                        endif
                                    endif
                                case default
                                    write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
                                end select
                                return
                            end function getSireDamHDIDByIndex


                            !---------------------------------------------------------------------------
                            !<@brief creates a new dummy animal at end of pedigree
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            ! PARAMETERS:
                            !---------------------------------------------------------------------------
                            subroutine createDummyAnimalAtEndOfPedigree(this,dummyId, offspringId)
                                class(PedigreeHolder) :: this
                                integer, optional :: offspringId !< offspring recoded id canbe given here
                                integer, intent(out) :: dummyId
                                character(len=IDLENGTH) :: tmpCounterStr

                                this%pedigreeSize = this%pedigreeSize+1
                                this%nDummys = this%nDummys + 1

                                this%isSorted = 0

                                write(tmpCounterStr, '(a,I3.3)') dummyAnimalPrepre,this%nDummys
                                this%Pedigree(this%pedigreeSize) =  Individual(tmpCounterStr ,'0','0', this%pedigreeSize,nsnps=this%nsnpsPopulation)
                                call this%dictionary%addKey(tmpCounterStr, this%pedigreeSize)
                                this%Pedigree(this%pedigreeSize)%isDummy = .true.
                                call this%Founders%list_add(this%Pedigree(this%pedigreeSize))
                                this%Pedigree(this%pedigreeSize)%founder = .true.

                                if (present(offspringId)) then
                                    if (offspringId > this%pedigreeSize) then
                                        write(error_unit,*) "ERROR - dummy list given index larger than pedigree"
                                    endif

                                    call this%Pedigree(this%pedigreeSize)%AddOffspring(this%pedigree(offspringId))

                                    if (.not. associated(this%pedigree(offspringId)%sirePointer)) then
                                        this%pedigree(offspringId)%sirePointer => this%Pedigree(this%pedigreeSize)
                                    else if (.not. associated(this%pedigree(offspringId)%damPointer)) then
                                        this%pedigree(offspringId)%damPointer => this%Pedigree(this%pedigreeSize)
                                    else
                                        write(error_unit,*) "ERROR - dummy animal given offspring that already has both parents!"
                                    end if

                                endif
                                dummyId = this%pedigreeSize
                            end subroutine createDummyAnimalAtEndOfPedigree




                            !---------------------------------------------------------------------------
                            !<@brief creates a new animal at end of pedigree
                            !<If genotype is supplied, animal is set to hd
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine addAnimalAtEndOfPedigree(this, originalID, geno)
                                class(PedigreeHolder) :: this
                                character(len=IDLENGTH) ,intent(in):: OriginalId
                                integer(kind=1), dimension(:), intent(in), optional :: geno

                                ! change pedigree to no longer be sorted

                                this%issorted = 0

                                this%pedigreeSize = this%pedigreeSize+1
                                this%Pedigree(this%pedigreeSize) =  Individual(OriginalId ,'0','0', this%pedigreeSize,nsnps=this%nsnpsPopulation)
                                call this%dictionary%addKey(OriginalId, this%pedigreeSize)
                                this%Pedigree(this%pedigreeSize)%isDummy = .false.
                                call this%Founders%list_add(this%Pedigree(this%pedigreeSize))
                                this%Pedigree(this%pedigreeSize)%founder = .true.

                                if (present(geno)) then
                                    call this%setAnimalAsGenotyped(this%pedigreeSize, geno)
                                    !   TODO make sure animal is actually hd
                                    call this%setAnimalAsHD(this%pedigreeSize)
                                endif
                            end subroutine addAnimalAtEndOfPedigree



                            !---------------------------------------------------------------------------
                            !> @brief Sets all individuals genotype from the Haplotype
                            !> @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !---------------------------------------------------------------------------
                            SUBROUTINE MakeGenotype(this)
                                class(PedigreeHolder) :: this
                                ! Any individual that has a missing genotype information but has both alleles
                                ! known, has its genotype filled in as the sum of the two alleles
                                integer :: i

                                !$OMP PARALLEL DO &
                                !$OMP PRIVATE(i)
                                do i=1,this%pedigreeSize
                                    call this%pedigree(i)%makeIndividualGenotypeFromPhase()    
                                enddo 
                                !$OMP END PARALLEL DO


                            END SUBROUTINE MakeGenotype

                            !#############################################################################################################################################################################################################################



                            !---------------------------------------------------------------------------
                            !> @brief Sets the individual haplotypes from the compilement if animal is genotyped for all animals
                            !> @author  David Wilson david.wilson@roslin.ed.ac.uk
                            !> @date    October 26, 2016
                            !---------------------------------------------------------------------------
                            subroutine PhaseComplement(this)
                                class(PedigreeHolder) :: this
                                ! If the genotype at a locus for an individual is known and one of its alleles has been determined
                                ! then impute the missing allele as the complement of the genotype and the known phased allele
                                integer :: i

                                !$OMP PARALLEL DO &
                                !$OMP PRIVATE(i)
                                do i=1,this%pedigreeSize
                                    call this%pedigree(i)%makeIndividualPhaseCompliment()
                                enddo
                                !$OMP END PARALLEL DO

                            end subroutine PhaseComplement

                            !---------------------------------------------------------------------------
                            !<@brief  counts missing snps across all animals at every snp
                            !<Returns a count of missing snps across all animals at every snp
                            !<@author  David Wilson david.wilson@roslin.ed.ac.uk
                            !<@date    October 26, 2016
                            !---------------------------------------------------------------------------
                            function countMissingGenotypesNoDummys(this) result(res)
                                integer :: res
                                integer :: i
                                class(PedigreeHolder) :: this

                                res = 0
                                do i=1, this%pedigreeSize

                                    if (this%pedigree(i)%isDummy) cycle

                                    res = res + this%pedigree(i)%individualGenotype%numMissing()
                                end do
                            end function countMissingGenotypesNoDummys

                        end module PedigreeModule
