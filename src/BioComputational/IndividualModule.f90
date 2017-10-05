
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     IndividualModule.f90
!
! DESCRIPTION:
!> @brief    Module holding logic of Individual class. This represents the pedigree.
!
!> @details  Module holds the class to represent pedigree, and contains class subroutines.
!
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 Dwilson - Initial Version

!-------------------------------------------------------------------------------

module IndividualModule
    use constantModule, only : OFFSPRINGTHRESHOLD, NOGENERATIONVALUE
    use genotypeModule
    use HaplotypeModule
    use IntegerLinkedListModule

    use iso_fortran_env
    
    implicit none

    public :: Individual,individualPointerContainer,operator ( == ),compareIndividual, initIndividual
    
    private

    ! This type is required to have an array of pointers
    type IndividualPointerContainer
        type(Individual), pointer :: p

        contains
          final:: deallocateIndividualPointer
    end type IndividualPointerContainer

    type Individual
        character(len=:), allocatable :: originalID
        character(len=:), allocatable :: sireID
        character(len=:), allocatable :: damID
        integer :: generation
        integer :: id
        integer :: originalPosition
        integer(kind=1) :: gender 
        type(individual), pointer :: sirePointer
        type(individual), pointer :: damPointer
        type(individualPointerContainer), allocatable :: OffSprings(:) !holds array of given size
        integer :: nOffs  = 0 !number of offspring
        logical(kind=1) :: Founder     = .false.
        logical(kind=1) :: Genotyped   = .false.
        logical(kind=1) :: Sequenced   = .false.
        logical(kind=1) :: isPhased    = .false.
        logical(kind=1) :: HD          = .false.
        logical(kind=1) :: isDummy     = .false.  ! if this animal is not in the pedigree, this will be true
        logical(kind=1) :: isUnknownDummy = .false.
        type(genotype),allocatable :: individualGenotype, individualGenotypeSubset
        type(Haplotype),allocatable,dimension(:) :: individualPhase, individualPhaseSubset
        integer,dimension(:), allocatable :: referAllele, AlterAllele

        real(kind=real32), allocatable, dimension(:) :: genotypeProbabilities
        real(kind=real32), allocatable, dimension(:,:) :: phaseProbabilities 
        integer(kind=1), dimension(:,:), allocatable :: seg !< should be dimension nsnps,2
        integer, allocatable :: nHighDensityOffspring
        ! plant stuff
        integer :: MyFamilyAsOffspring = -99
        real(kind=real64) :: MendelianInconsistenciesPercFamily = 0
        logical(kind=1):: isInbred
        logical(kind=1)::isImputed
        logical(kind=1) :: IgnoreMe =.false.

        integer, allocatable, dimension(:) :: StrandBreakArray
        integer(kind=1), allocatable, dimension(:,:) :: PartialInformativeHaplotype, PartialInformativeHaplotypeFull
        integer(kind=1), allocatable, dimension(:,:) :: ParentalHaplotypes, ParentalHaplotypesFull
        ! plant stuff testing 
        integer :: NumberOfSwitches = 0
        integer, allocatable, dimension(:) :: BestSwitchPosArray, ParentInferredSwitchCount
        double precision :: MyGenoYield, MyGenoAccuracy, MyGenoAccuracyInf, MyGenoYieldInf, MyGenoCorrect, MyGenoCorrectInf
        double precision :: MyPhase1YieldInf, MyPhase1AccInf,  MyPhase1CorrectInf, MyPhase2YieldInf, MyPhase2AccInf, MyPhase2CorrectInf
        double precision, allocatable, dimension( :,:) :: MyPhaseYield, MyPhaseAcc, MyPhaseCorrect

        type(IntegerLinkedList) :: families

        integer, allocatable, dimension(:) :: inconsistencies !< number of consistencies an individual has overall, so each offsprings inconsistencies will add to this.


        contains
            procedure :: initIndividual
            procedure :: getSireDamByIndex
            procedure :: isGenotyped
            procedure :: isGenotypedNonMissing
            procedure :: setGenotypeArray
            procedure :: setPhaseArray
            procedure :: setGenotypeArraySubset
            procedure :: setPhaseArraySubset
            procedure :: isHD
            procedure :: SetHD
            procedure :: getSireId
            procedure :: getDamID
            procedure :: setSequenceArray
            procedure :: GetNumberOffsprings
            procedure :: GetOffsprings
            procedure :: AddOffspring
            procedure :: setGender
            final :: destroyIndividual
            procedure :: setGeneration
            procedure :: getSireDamObjectByIndex
            procedure :: getSireDamNewIDByIndex
            procedure :: getSireDamNewIDByIndexNoDummy
            procedure :: getIntegerVectorOfRecodedIds
            procedure :: getCharacterVectorOfRecodedIds
            procedure :: getPaternalGrandSireRecodedIndex
            procedure :: getMaternalGrandSireRecodedIndex
            procedure :: getPaternalGrandDamRecodedIndex
            procedure :: getMaternalGrandDamRecodedIndex
            procedure :: getParentGenderBasedOnIndex
            procedure :: hasDummyParent
            procedure :: hasDummyParentsOrGranparents
            procedure :: isDummyBasedOnIndex
            procedure :: isUnknownDummyBasedOnIndex
            procedure :: getPaternalGrandSireRecodedIndexNoDummy
            procedure :: getMaternalGrandSireRecodedIndexNoDummy
            procedure :: getPaternalGrandDamRecodedIndexNoDummy
            procedure :: getMaternalGrandDamRecodedIndexNoDummy
            procedure :: getIntegerVectorOfRecodedIdsNoDummy
            procedure :: resetOffspringInformation
            procedure :: removeOffspring
            procedure :: writeIndividual
            procedure :: initPhaseArrays
            procedure :: hasGenotypedAnsestors
            procedure :: getSeg
            procedure :: setSeg
            procedure :: setSegToMissing
            procedure :: makeIndividualPhaseCompliment
            procedure :: makeIndividualGenotypeFromPhase
            procedure :: countHighDensityOffspring
            procedure :: addFamily
            procedure :: initialisePlantArrays
            procedure :: getProbabilitiesFromOwnGenotypeAndPhase
            procedure :: initPhaseAndGenotypes
            procedure :: initGenotype
            generic:: write(formatted)=> writeIndividual

            
    end type Individual

    interface Individual
        module procedure initIndividual
    end interface Individual
    
    interface operator ( == )
        module procedure compareIndividual
    end interface operator ( == )

contains

  subroutine deallocateIndividualPointer(this)
    type(IndividualPointerContainer), intent(inout):: this

    this%p => null()

  end subroutine deallocateIndividualPointer


        !---------------------------------------------------------------------------
    !> @brief Constructor for individual class.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine initIndividual(this, originalID,sireIDIn,damIDIn, id, generation,gender, nsnps, probabilites)
        class(Individual) :: this
        character(*), intent(in) :: originalID,sireIDIn,damIDIn
        integer, intent(in), Optional :: generation
        integer, intent(in) :: id
        integer, intent(in), Optional :: nsnps !< number of snps to initialise default genotype class
        integer(kind=1), intent(in), Optional :: gender
        logical, optional :: probabilites !< if present, allocate probabilites


        allocate(character(len=len(originalID)) ::this%originalID)
        allocate(character(len=len(sireIDIn)) ::this%sireID)
        allocate(character(len=len(sireIDIn)) ::this%damID)
        allocate(this%OffSprings(OFFSPRINGTHRESHOLD))
        this%originalID = originalID
        this%id = id
        this%gender = -9
        this%sireId = sireIDIn
        this%damId = damIDIn
        if (present(generation)) then
            this%generation = generation
        else
            this%generation = NOGENERATIONVALUE
        endif
        if (present(gender)) then
            this%gender = gender
        endif

        if (present(nsnps)) then
            if (nsnps /= 0) then
                allocate(this%individualPhase(2))
                allocate(this%individualGenotype)
                allocate(this%inconsistencies(nsnps))
                call this%individualGenotype%newGenotypeMissing(nsnps)
                call this%individualPhase(1)%newHaplotypeMissing(nsnps)
                call this%individualPhase(2)%newHaplotypeMissing(nsnps)
                this%inconsistencies =  0
            endif

            if (present(probabilites)) then
                allocate(this%genotypeProbabilities(nsnps))
                allocate(this%phaseProbabilities(2,nsnps))
            endif
        endif



    end subroutine initIndividual

     !---------------------------------------------------------------------------
    !> @brief Deallocates individual object
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine destroyIndividual(this)
        type(Individual) :: this

        if (allocated(this%offsprings)) then
            deallocate(this%offsprings)
        endif


        if (allocated(this%originalId)) then
            deallocate(this%originalID)
            deallocate(this%sireID)
            deallocate(this%damID)
        endif

        if (allocated(this%referAllele)) then
            deallocate(this%referAllele)
            deallocate(this%alterAllele)
        endif

        if (allocated(this%seg)) then
            deallocate(this%seg)
        endif

        if (allocated(this%individualPhase)) then
            deallocate(this%individualPhase)
        endif

        if (allocated(this%individualGenotype)) then
            deallocate(this%individualGenotype)
        endif

        if (allocated(this%genotypeProbabilities)) then
            deallocate(this%genotypeProbabilities)
        endif

         if (allocated(this%phaseProbabilities)) then
            deallocate(this%phaseProbabilities)
        endif

    
        if (allocated(this%BestSwitchPosArray)) then

            deallocate(this%BestSwitchPosArray)
        endif
        
        if (allocated(this%ParentInferredSwitchCount)) then

            deallocate(this%ParentInferredSwitchCount)
        endif
    if (allocated(this%StrandBreakArray)) then

            deallocate(this%StrandBreakArray)
        endif

        if (allocated(this%PartialInformativeHaplotype)) then

            deallocate(this%PartialInformativeHaplotype)
        endif

          if (allocated(this%PartialInformativeHaplotypeFull)) then

            deallocate(this%PartialInformativeHaplotypeFull)
        endif

        if (allocated(this%ParentInferredSwitchCount)) then

            deallocate(this%ParentInferredSwitchCount)
        endif

        if (allocated(this%MyPhaseYield)) then

            deallocate(this%MyPhaseYield)
        endif

         if (allocated(this%MyPhaseAcc)) then

            deallocate(this%MyPhaseAcc)
        endif

        if (allocated(this%MyPhaseCorrect)) then

            deallocate(this%MyPhaseCorrect)
        endif
        if (allocated(this%inconsistencies)) then
            deallocate(this%inconsistencies)
        endif

    end subroutine destroyIndividual



    subroutine setSeg(this, location, parent, value) 

        class(individual) :: this
        integer, intent(in) :: location, parent, value

        if (.not. allocated(this%seg)) then
            allocate(this%seg(this%individualGenotype%length,2 ))
        endif


        this%seg(location, parent) = value

    end subroutine setSeg


        subroutine setSegToMissing(this) 

        class(individual) :: this

        if (.not. allocated(this%seg)) then
            allocate(this%seg(this%individualGenotype%length,2 ))
        endif


        this%seg = MISSINGGENOTYPECODE

    end subroutine setSegToMissing


    function getSeg(this, location,parent) result(res)
        implicit none

        class(individual) :: this
        integer, intent(in) :: location, parent
        integer(kind =1) :: res
        res = this%seg(location, parent)

        
    end function getSeg
     !---------------------------------------------------------------------------
    !> @brief Returns true if individuals are equal, false otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    logical function compareIndividual(l1,l2)
        class(Individual), intent(in) :: l1
        type (Individual) , intent(in) :: l2 !< individuals to compare

        if (l1%id == l2%id .and. l1%sireID == l2%sireID .and. l1%damID == l2%damID) then
            compareIndividual=.true.
        else
            compareIndividual=.false.
        endif

        return
    end function compareIndividual

        !---------------------------------------------------------------------------
    !< @brief returns true if individual has genotyped ancestors up to certrain gen, false otherwise
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !< @date    October 26, 2016
    !---------------------------------------------------------------------------
    recursive function hasGenotypedAnsestors(ind,count) result(res)
    class(individual) ,intent(in) :: ind !< individual to check ancestors of
    integer, intent(in) :: count !< how many generations should we look for
    logical :: res

    if (count == 0) then
        res = .false.
        return
    endif

    if (associated(ind%damPointer)) then

        if (ind%damPointer%isGenotypedNonMissing()) then
            res= .true.
            return
        endif
    endif
    if (associated(ind%sirePointer)) then

        if (ind%sirePointer%isGenotypedNonMissing()) then
            res= .true.
            return
        endif
        res = hasGenotypedAnsestors(ind%sirePointer, count -1)
        if (res) then
            return
        endif
    endif

    ! This is done to insure no recursion is done if it is not neccessary, as recursion is more expensive than the branch, which can be trivially optimised
    if (associated(ind%damPointer)) then
        res = hasGenotypedAnsestors(ind%damPointer, count -1)
        if (res) then
            return
        endif
    endif

    end function hasGenotypedAnsestors

         !---------------------------------------------------------------------------
    !> @brief output for individual
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine writeIndividual(dtv, unit, iotype, v_list, iostat, iomsg)
        class(Individual), intent(in) :: dtv         !< Object to write.
        integer, intent(in) :: unit         !< Internal unit to write to.
        character(*), intent(in) :: iotype  !< LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)    !< parameters from fmt spec.
        integer, intent(out) :: iostat      !< non zero on error, etc.
        character(*), intent(inout) :: iomsg  !< define if iostat non zero.

        write(unit, *, iostat = iostat, iomsg = iomsg) dtv%originalID,",", dtv%sireId,",", dtv%damID 

    end subroutine writeIndividual


     !---------------------------------------------------------------------------
    !> @brief Returns either the individuals id, the sires id or dams id based on
    !> which index is passed.

    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    function getSireDamByIndex(this, index) result(v)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual), intent(in) :: this
        integer, intent(in) :: index !< index of value to return (1 for thisID, 2 for sireID, 3 for damID)
        character(:),allocatable :: v
          
        select case (index)
            case(1)
                v = this%originalID
            case(2)
                if (associated(this%sirePointer)) then
                  v = this%sirePointer%originalId
                else
                  v = "0"
                endif

            case(3)
                if (associated(this%damPointer)) then
              v = this%damPointer%originalId
              else
                  v = "0"
            endif
            case default
                write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
        end select
        return
    end function getSireDamByIndex

    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of paternal grand sire, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getPaternalGrandSireRecodedIndex(this)
        class(individual) :: this

        if (associated(this%sirePointer)) then
            if (associated(this%sirePointer%sirePointer)) then
                getPaternalGrandSireRecodedIndex = this%sirePointer%sirePointer%id
                return
            endif
        endif
        getPaternalGrandSireRecodedIndex = 0
    end function getPaternalGrandSireRecodedIndex


      !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of paternal grand sire, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getPaternalGrandSireRecodedIndexNoDummy(this)
        class(individual) :: this

        if (associated(this%sirePointer)) then
            if (associated(this%sirePointer%sirePointer)) then
                if (.not. this%sirePointer%sirePointer%isDummy) then
                    getPaternalGrandSireRecodedIndexNoDummy = this%sirePointer%sirePointer%id
                    return
                endif   
            endif
        endif
        getPaternalGrandSireRecodedIndexNoDummy = 0
    end function getPaternalGrandSireRecodedIndexNoDummy

    !---------------------------------------------------------------------------
    !> @brief Returns the gender of individual if index of 1 is given.
    !> returns the gender of sire if index of 1 is given (if available), 
    !> returns the gender of dam if index of 2 is given (if available).
    !> if specified parent gender is not available, 0 will be returned
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    function getParentGenderBasedOnIndex(this, index) result(v)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual), intent(in) :: this
        integer, intent(in) :: index !< index of animal to get gender of (1 for this, 2 for sireGender, 3 for damGender)
        integer :: v
        select case (index)
            case(1)
                v = this%gender
            case(2)
                if (associated(this%sirePointer)) then
                    v = this%sirePointer%gender
                else
                    v = 0
                endif
            case(3)
                if (associated(this%damPointer)) then
                    v = this%damPointer%gender
                else
                    v = 0
                endif
            case default
                write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
        end select
    end function getParentGenderBasedOnIndex

    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of maternal grand sire, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getMaternalGrandSireRecodedIndex(this)
        class(individual) :: this

        if (associated(this%damPointer)) then
            if (associated(this%damPointer%sirePointer)) then
                getMaternalGrandSireRecodedIndex = this%damPointer%sirePointer%id
                return
            endif
        endif
        getMaternalGrandSireRecodedIndex = 0
    end function getMaternalGrandSireRecodedIndex


    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of maternal grand sire, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getMaternalGrandSireRecodedIndexNoDummy(this)
        class(individual) :: this

        if (associated(this%damPointer)) then
            if (associated(this%damPointer%sirePointer)) then
                if (.not. this%damPointer%sirePointer%isDummy) then
                    getMaternalGrandSireRecodedIndexNoDummy = this%damPointer%sirePointer%id
                    return
                endif
            endif
        endif
        getMaternalGrandSireRecodedIndexNoDummy = 0
    end function getMaternalGrandSireRecodedIndexNoDummy

    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of paternal grand dam, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getPaternalGrandDamRecodedIndex(this)
        class(individual) :: this

        if (associated(this%sirePointer)) then
            if (associated(this%sirePointer%damPointer)) then
                getPaternalGrandDamRecodedIndex = this%sirePointer%damPointer%id
                return
            endif
        endif
        getPaternalGrandDamRecodedIndex = 0
    end function getPaternalGrandDamRecodedIndex

     !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of paternal grand dam, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getPaternalGrandDamRecodedIndexNoDummy(this)
        class(individual) :: this

        if (associated(this%sirePointer)) then
            if (associated(this%sirePointer%damPointer)) then
                if (.not. this%sirePointer%damPointer%isDummy) then
                    getPaternalGrandDamRecodedIndexNoDummy = this%sirePointer%damPointer%id
                    return
                endif
            endif
        endif
        getPaternalGrandDamRecodedIndexNoDummy = 0
    end function getPaternalGrandDamRecodedIndexNoDummy

    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of maternal grand dam, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getMaternalGrandDamRecodedIndex(this)
        class(individual) :: this

        if (associated(this%damPointer)) then
            if (associated(this%damPointer%damPointer)) then
                getMaternalGrandDamRecodedIndex = this%damPointer%damPointer%id
                return
            endif
        endif
        getMaternalGrandDamRecodedIndex = 0
    end function getMaternalGrandDamRecodedIndex



    !---------------------------------------------------------------------------
    !> @brief Returns the index in the pedigree of maternal grand dam, or 0 otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    integer function getMaternalGrandDamRecodedIndexNoDummy(this)
        class(individual) :: this

        if (associated(this%damPointer)) then
            if (associated(this%damPointer%damPointer)) then
                if (.not. this%damPointer%damPointer%isDummy) then
                    getMaternalGrandDamRecodedIndexNoDummy = this%damPointer%damPointer%id
                    return
                endif
            endif
        endif
        getMaternalGrandDamRecodedIndexNoDummy = 0
    end function getMaternalGrandDamRecodedIndexNoDummy




    !---------------------------------------------------------------------------
    !> @brief Returns an array of recoded id's where index 1 is individuals id,
    !> index 2 is sire's recoded ID (0 if not available),
    !> index 3 is dam's recoded ID (0 if not available)
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    function getIntegerVectorOfRecodedIds(this) result(res)
        class(Individual) :: this!< 
        integer :: res(3)

        res = 0
        res(1) = this%id
        if (associated(this%sirePointer)) then
            res(2) = this%sirePointer%id
        endif

        if (associated(this%damPointer)) then
            res(3) = this%damPointer%id
        endif

    end function getIntegerVectorOfRecodedIds



        !---------------------------------------------------------------------------
    !> @brief Returns an array of recoded id's where index 1 is individuals id,
    !> index 2 is sire's recoded ID (0 if not available),
    !> index 3 is dam's recoded ID (0 if not available)
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    function getCharacterVectorOfRecodedIds(this) result(res)
        use constantModule, only : IDLENGTH
        class(Individual) :: this
        character(IDLENGTH) :: res(3)

        res = '0'
        res(1) = this%originalID
        if (associated(this%sirePointer)) then
            res(2) = this%sirePointer%originalID
        endif

        if (associated(this%damPointer)) then
            res(3) = this%damPointer%originalID
        endif

    end function getCharacterVectorOfRecodedIds

        !---------------------------------------------------------------------------
    !> @brief Returns an array of recoded id's where index 1 is individuals id,
    !> index 2 is sire's recoded ID (0 if not available or if animal is a dummy),
    !> index 3 is dam's recoded ID (0 if not available or if animal is a dummy)
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    function getIntegerVectorOfRecodedIdsNoDummy(this) result(res)
        class(Individual) :: this
        integer :: res(3)

        res = 0
        res(1) = this%id
        if (associated(this%sirePointer)) then
            if (this%sirePointer%isDummy) then
                res(2) = 0
            else
                res(2) = this%sirePointer%id
            endif
        endif

        if (associated(this%damPointer)) then
            if (this%damPointer%isDummy) then
                res(3) = 0
            else
                res(3) = this%damPointer%id
            endif
        endif

    end function getIntegerVectorOfRecodedIdsNoDummy
        
        


!---------------------------------------------------------------------------
    !> @brief Returns either the individual object, the sires object or dams object based on
    !> which index is passed.

    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    function getSireDamObjectByIndex(this, index) result(v)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual),target, intent(in) :: this
        integer, intent(in) :: index !< index of object to return (1 for this, 2 for sire, 3 for dam)
        type(individual), pointer :: v
        select case (index)
            case(1)
                v => this
            case(2)
                v => this%sirePointer
            case(3)
                v => this%damPointer
            case default
                write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
        end select
        return
    end function getSireDamObjectByIndex


!---------------------------------------------------------------------------
    !> @brief returns true if index of corresponding parent is dummy
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    logical function isDummyBasedOnIndex(this, index)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual),target, intent(in) :: this
        integer, intent(in) :: index !< index of object to return (1 for this, 2 for sire, 3 for dam)
        select case (index)
            case(1)
                isDummyBasedOnIndex = this%isDummy
            case(2)
                if (associated(this%sirePointer)) then
                    isDummyBasedOnIndex = this%sirePointer%isDummy
                    return
                endif
            case(3)
                if (associated(this%damPointer)) then
                    isDummyBasedOnIndex = this%damPointer%isDummy
                    return
                endif
            case default
                write(error_unit, *) "error: getSireDamObjectByIndex has been given an out of range value"
        end select
        isDummyBasedOnIndex = .false.
        
    end function isDummyBasedOnIndex


!---------------------------------------------------------------------------
    !> @brief returns true if index of corresponding parent is dummy
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
    !---------------------------------------------------------------------------
    logical function isUnknownDummyBasedOnIndex(this, index)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual),target, intent(in) :: this
        integer, intent(in) :: index !< index of object to return (1 for this, 2 for sire, 3 for dam)
        select case (index)
            case(1)
                isUnknownDummyBasedOnIndex = this%isUnknownDummy
            case(2)
                if (associated(this%sirePointer)) then
                    isUnknownDummyBasedOnIndex = this%sirePointer%isUnknownDummy
                    return
                endif
            case(3)
                if (associated(this%damPointer)) then
                    isUnknownDummyBasedOnIndex = this%damPointer%isUnknownDummy
                    return
                endif
            case default
                write(error_unit, *) "error: getSireDamObjectByIndex has been given an out of range value"
        end select
        isUnknownDummyBasedOnIndex = .false.
        
    end function isUnknownDummyBasedOnIndex


         !---------------------------------------------------------------------------
    !> @brief Returns either the individuals id, the sires id or dams id based on
    !> which index is passed.

    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !---------------------------------------------------------------------------
    function getSireDamNewIDByIndex(this, index) result(v)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual), intent(in) :: this
        integer, intent(in) :: index !< index of object to return (1 for this, 2 for sire, 3 for dam)
        integer:: v
        select case (index)
            case(1)
                v = this%id
            case(2)
                if (associated(this%sirePointer)) then
                    v = this%sirePointer%id
                else
                    v = 0
                endif
            case(3)
                if (associated(this%damPointer)) then
                    v = this%damPointer%id
                else
                    v = 0
                endif
            case default
                write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
        end select
        return
    end function getSireDamNewIDByIndex


            !---------------------------------------------------------------------------
        !> @brief Returns either the individuals id, the sires id or dams id based on
        !> which index is passed. 0 is returned if no parent or if parent is a dummy

        !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
        !> @author  David Wilson david.wilson@roslin.ed.ac.uk
        !> @date    October 26, 2016
        ! PARAMETERS:
        !> @param[in] index - the index
        !---------------------------------------------------------------------------
        function getSireDamNewIDByIndexNoDummy(this, index) result(v)
        use iso_fortran_env, only : ERROR_UNIT
        class(Individual), intent(in) :: this
        integer, intent(in) :: index !< index of object to return (1 for this, 2 for sire, 3 for dam)
        integer:: v
        v = 0
        select case (index)
            case(1)
                v = this%id
            case(2)
                if (associated(this%sirePointer)) then
                    if (.not. this%sirePointer%isDummy) then
                        v = this%sirePointer%id
                    endif
                endif
            case(3)
                if (associated(this%damPointer)) then
                    if (.not. this%damPointer%isDummy) then
                        v = this%damPointer%id
                    endif
                endif
            case default
                write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
        end select
        return
    end function getSireDamNewIDByIndexNoDummy

    !---------------------------------------------------------------------------
    !> @brief returns true if either paretns are a dummy animal
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    logical function hasDummyParent(this)
        class(Individual), intent(in) :: this
        if (associated(this%sirePointer)) then
            if (this%sirePointer%isDummy) then
                hasDummyParent = .true.
                return
            endif
        endif

        if (associated(this%damPointer)) then
            if (this%damPointer%isDummy) then
                hasDummyParent = .true.
                return
            endif
        endif


        hasDummyParent = .false.

    end function hasDummyParent



       !---------------------------------------------------------------------------
    !> @brief returns true if either parents or grandparents are a dummy animal
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    logical function hasDummyParentsOrGranparents(this)
        class(Individual), intent(in) :: this
        if (associated(this%sirePointer)) then
            if (this%sirePointer%isDummy) then
                hasDummyParentsOrGranparents = .true.
                return
            else
                if (associated(this%sirePointer%sirePointer)) then
                   if (this%sirePointer%sirePointer%isDummy) then
                        hasDummyParentsOrGranparents = .true.
                        return
                    endif
                endif
                if (associated(this%sirePointer%damPointer)) then
                   if (this%sirePointer%damPointer%isDummy) then
                        hasDummyParentsOrGranparents = .true.
                        return
                    endif
                endif
            endif
        endif

        if (associated(this%damPointer)) then
            if (this%damPointer%isDummy) then
                hasDummyParentsOrGranparents = .true.
                return
            else
                if (associated(this%damPointer%sirePointer)) then
                   if (this%damPointer%sirePointer%isDummy) then
                        hasDummyParentsOrGranparents = .true.
                        return
                    endif
                endif
                if (associated(this%damPointer%damPointer)) then
                   if (this%damPointer%damPointer%isDummy) then
                        hasDummyParentsOrGranparents = .true.
                        return
                    endif
                endif
            endif

        endif


        hasDummyParentsOrGranparents = .false.

    end function hasDummyParentsOrGranparents


      
    !---------------------------------------------------------------------------
    !> @brief Resets the offspring information for a given animal
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine resetOffspringInformation(this)
        class(Individual) :: this
    
        deallocate(this%offsprings)
        allocate(this%OffSprings(OFFSPRINGTHRESHOLD))
    
        this%noffs = 0
    
    
    end subroutine resetOffspringInformation 

  

    !---------------------------------------------------------------------------
    !> @brief boolean function returning true if object is a founder of a generation
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if object is a founder of a generation
    !---------------------------------------------------------------------------
    elemental function isFounder(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%Founder
    end function isFounder

    !---------------------------------------------------------------------------
    !> @brief sets object to be a founder
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if object is a founder of a generation
    !---------------------------------------------------------------------------
    subroutine SetObjectAsFounder(this)
        class(Individual), intent(inout) :: this
        this%Founder = .true.
    end subroutine SetObjectAsFounder

    !---------------------------------------------------------------------------
    !> @brief returns true if object is genotyped
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if object is genotyped
    !---------------------------------------------------------------------------
    elemental function isGenotyped(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%Genotyped
    end function isGenotyped
!---------------------------------------------------------------------------
    !> @brief returns true if object is genotyped
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if object is genotyped
    !---------------------------------------------------------------------------
    elemental function isGenotypedNonMissing(this) result(ans)
        class(Individual), intent(in) :: this
        integer(kind=1), dimension(:), allocatable :: geno
        logical :: ans
        ans = this%Genotyped

        if(ans) then
            geno = this%individualGenotype%toIntegerArray()
            ans = any(geno == 0 .or. geno == 1 .or. geno==2)
        endif
    end function isGenotypedNonMissing


    !---------------------------------------------------------------------------
    !> @brief Sets the individual genotype from the Haplotype
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine makeIndividualGenotypeFromPhase(this)
        class(Individual), intent(inout) :: this
        
        call this%IndividualGenotype%setFromHaplotypesIfMissing(this%individualPhase(1),this%individualPhase(2)) 
    end subroutine makeIndividualGenotypeFromPhase



    !---------------------------------------------------------------------------
    !> @brief Sets the individual haplotypes from the compilement if animal is genotyped
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine makeIndividualPhaseCompliment(this)
        class(Individual), intent(inout) :: this
        type (haplotype) :: comp1, comp2
        
        call this%individualPhase(1)%setErrorToMissing()
        call this%individualPhase(2)%setErrorToMissing()

        ! allocate(comp1)
        ! allocate(comp2)
        comp2 = this%individualGenotype%complement(this%individualPhase(1))
        comp1 = this%individualGenotype%complement(this%individualPhase(2))

        call comp1%setErrorToMissing()
        call comp2%setErrorToMissing()

        call this%individualPhase(1)%setFromOtherIfMissing(comp1)
        call this%individualPhase(2)%setFromOtherIfMissing(comp2)
        
        ! deallocate(comp1)
        ! deallocate(comp2)

    end subroutine makeIndividualPhaseCompliment

    !---------------------------------------------------------------------------
    !> @brief Sets the individual to be genotyped.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine setGenotypeArray(this, geno, lockIn)
        use constantModule
        
        class(Individual), intent(inout) :: this
        integer(KIND=1), dimension(:), intent(in) :: geno !< One dimensional array of genotype information
        logical, intent(in), optional :: lockIn

        this%Genotyped = .true.
        !TODO this%Genotyped = any(geno == 1 .or. geno == 2 .or. geno == 0)
        ! this%Genotyped = any(geno == 1 .or. geno == 2 .or. geno == 0)

        if (present(lockIn)) then

            if (lockIn) then
                call this%individualGenotype%Genotype(Geno,lock=1)
                return
            endif
        endif

        if (allocated(this%individualGenotype)) then
            deallocate(this%individualGenotype)
            allocate(this%individualGenotype)
        endif

        call this%individualGenotype%Genotype(Geno)
    end subroutine setGenotypeArray


    !---------------------------------------------------------------------------
    !> @brief Sets the individual to be genotyped.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine setSequenceArray(this, referAllele, alterAllele)
        use constantModule
        
        class(Individual), intent(inout) :: this
        integer, dimension(:), intent(in) :: referAllele, alterAllele

        this%Genotyped = .true.
        this%Sequenced = .true.
        
        this%referAllele = referAllele
        this%alterAllele  = alterAllele
    end subroutine setSequenceArray


    !---------------------------------------------------------------------------
    !> @brief For plants, sets subset genotypes to use (for imputation)
    !> @author  Serap Gonen serap.gonen@roslin.ed.ac.uk
    !> @date    August 30, 2017
    !---------------------------------------------------------------------------
    subroutine setGenotypeArraySubset(this, geno)
        use constantModule

        class(Individual), intent(inout) :: this
        integer(KIND=1), dimension(:), intent(in) :: geno !< One dimensional array of genotype information

        call this%individualGenotypeSubset%Genotype(geno)

    end subroutine setGenotypeArraySubset


    !---------------------------------------------------------------------------
    !> @brief Sets the individual's phase.
    !> @author  Daniel Money, daniel.money@roslin.ed.ac.uk
    !> @date    June 19, 2017
    !---------------------------------------------------------------------------
    subroutine setPhaseArray(this, hap, phase)
        use constantModule
        
        class(Individual), intent(inout) :: this
        integer, intent(in) :: hap
        integer(KIND=1), dimension(:), intent(in) :: phase !< One dimensional array of phase information

        !TODO: Should we be checking that genotype is set and if not set it to missing?
        call this%individualPhase(hap)%Haplotype(phase)
    end subroutine setPhaseArray


    !---------------------------------------------------------------------------
    !> @brief For plants, sets subset phase to use (for imputation)
    !> @author  Serap Gonen serap.gonen@roslin.ed.ac.uk
    !> @date    August 30, 2017
    !---------------------------------------------------------------------------
    subroutine setPhaseArraySubset(this, hap, phase)
        use constantModule

        class(Individual), intent(inout) :: this
        integer, intent(in) :: hap
        integer(KIND=1), dimension(:), intent(in) :: phase !< One dimensional array of phase information
                

        call this%individualPhaseSubset(hap)%Haplotype(phase)

    end subroutine setPhaseArraySubset

    !---------------------------------------------------------------------------
    !> @brief initialises an individual phases arrays given the number of snps
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine initPhaseArrays(this, nsnp)
        use constantModule
        class(Individual) :: this
        integer, intent(in) :: nsnp

        if (allocated(this%individualPhase)) then
            deallocate(this%individualPhase)
        endif
        if (.not. allocated(this%individualPhase)) then
            allocate(this%individualPhase(2))
        endif
        call this%individualPhase(1)%Haplotype(nSnp)
        call this%individualPhase(2)%Haplotype(nSnp)
    end subroutine initPhaseArrays


     !---------------------------------------------------------------------------
    !> @brief initialises an individual genotype given the number of snps
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine initGenotype(this, nsnp)
        use constantModule
        class(Individual) :: this
        integer, intent(in) :: nsnp

        call this%individualGenotype%Genotype(nSnp)
    end subroutine initGenotype


    !---------------------------------------------------------------------------
    !> @brief initialises an individuals phase and genotype
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine initPhaseAndGenotypes(this, nsnp)
        use constantModule
        class(Individual) :: this
        integer, intent(in) :: nsnp

        call this%initGenotype(nsnp)
        call this%initPhaseArrays(nsnp)
    end subroutine initPhaseAndGenotypes

    !---------------------------------------------------------------------------
    !> @brief returns true if the individual is genotyped at high density.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if individual is HD
    !---------------------------------------------------------------------------
    elemental function isHD(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%HD
    end function isHD

    !---------------------------------------------------------------------------
    !> @brief Sets the individual to be genotyped at high density.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine SetHD(this)
        class(Individual), intent(inout) :: this
        this%HD = .true.
    end subroutine SetHD

    !------------------------\---------------------------------------------------
    !> @brief Adds an individual as offspring
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine AddOffspring(this, offspringToAdd)
        class(Individual), intent(inout) :: this
        class(Individual),target, intent(in) :: offspringToAdd
        type(individualPointerContainer), allocatable :: tmp(:)
        integer :: motherId, fatherId
        this%nOffs = this%nOffs + 1
        
        motherId = this%getSireDamNewIDByIndexNoDummy(3)
        fatherId = this%getSireDamNewIDByIndexNoDummy(2)
        if (offspringTOAdd%id == motherId .or. offspringToAdd%id == fatherID) then

            write(error_unit,*) "ERROR: Animal ", this%originalID ," has been given animal ", offspringToAdd, " as both parent and offspring"
            stop
        
        endif
        if (this%nOffs > OFFSPRINGTHRESHOLD) then
            allocate(tmp(this%nOffs))
            tmp(1:size(this%Offsprings)) = this%Offsprings
            call move_alloc(tmp,this%OffSprings)
        endif
        this%OffSprings(this%nOffs)%p => offspringToAdd
    end subroutine AddOffspring

    !---------------------------------------------------------------------------
    !> @brief gets number of offspring of individual
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return integer of number of individuals
    !---------------------------------------------------------------------------
    elemental function GetNumberOffsprings(this) result(ans)
        class(Individual), intent(in) :: this
        integer :: ans

        ans = this%nOffs

    end function GetNumberOffsprings

    !---------------------------------------------------------------------------
    !> @brief gets array of *individual* objects that are this individuals offspring
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return array of pointers of individuals which are offspring of this parent
    !---------------------------------------------------------------------------
    subroutine GetOffsprings(this, Offsprings)

        type(individualPointerContainer), allocatable :: Offsprings(:)
        class(Individual), intent(in) :: this

        Offsprings = this%Offsprings

    end subroutine GetOffsprings


  

    !---------------------------------------------------------------------------
    !> @brief Sets gender info of individual. 1 signifies male, 2 female. 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine setGender(this,gender) 
        class(individual) :: this
        integer(kind=1),intent(in) :: gender !< gender (1 or 2) to set animal
        this%gender = gender
    end subroutine setGender

    !---------------------------------------------------------------------------
    !> @brief returns gender info of individual. 1 signifies male, 2 female. 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return gender (integer kind(1)) 
    !---------------------------------------------------------------------------
    function getGender(this) result(gender)
        class(individual) :: this
        integer(kind=1) :: gender
        gender = this%gender
    end function getGender

    !---------------------------------------------------------------------------
    !> @brief returns gender info of individual. 1 signifies male, 2 female. 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return String of SireID
    !---------------------------------------------------------------------------
    function getSireID(this) result(v)
        class(Individual), intent(in) :: this
        character(:),allocatable :: v
        v = this%sireID
    end function getSireID


    function getDamID(this) result(v)
        class(Individual), intent(in) :: this
        character(:),allocatable :: v
        v = this%damID
    end function getDamID

    !---------------------------------------------------------------------------
    !> @brief Sets generation info of individual
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] generation (integer) 
    !---------------------------------------------------------------------------
    subroutine setGeneration(this,generation) 
        class(individual) :: this
        integer,intent(in) :: generation
        this%generation = generation
    end subroutine setGeneration



    !---------------------------------------------------------------------------
    !> @brief removes offspring information and disassociatespointers for given animal
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @param[in] generation (integer) 
    !---------------------------------------------------------------------------
    subroutine removeOffspring(this, offspring)
        use iso_fortran_env
        class(individual ) :: this
        type(individual) :: offspring
        integer :: i,h, old
        logical :: found
        found =.false.

        do i=1, this%nOffs
            
            if (compareIndividual(this%offsprings(i)%p, offspring)) then

                if (compareIndividual(this%offsprings(i)%p%sirePointer, this)) then
                    this%offsprings(i)%p%sirePointer => null() 
                    old = this%nOffs-1
                    do h=i, old 
                        this%offsprings(h)%p => this%offsprings(h+1)%p
                    enddo
                    this%nOffs = this%nOffs - 1
                    found = .true.
                    exit
                else if (compareIndividual(this%offsprings(i)%p%damPointer, this)) then
                    this%offsprings(i)%p%damPointer => null() 
                    old = this%nOffs-1
                    do h=i, old
                        this%offsprings(h)%p => this%offsprings(h+1)%p
                    enddo
                    this%nOffs = this%nOffs - 1
                    found = .true.
                    exit
                else
                    write(error_unit,*) "WARNING: parent isn't present in offspring to remove"
                endif
                return
            endif
        enddo


        if (.not. found) then
            write(error_unit,*) "WARNING: unknown offspring trying to be removed from parent"
        endif
    end subroutine removeOffspring




  !---------------------------------------------------------------------------
    !> @brief counts the number of offspring flagged as hd
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !---------------------------------------------------------------------------
    function countHighDensityOffspring(this) result(count)
        class(individual) :: this 
        integer :: count !< the number of high density offspring this individual has
        integer :: i

        if (allocated(this%nHighDensityOffspring)) then
            count = this%nHighDensityOffspring
        else 
            count = 0
            allocate(this%nHighDensityOffspring)
            do i=1,this%nOffs
                if (this%OffSprings(i)%p%hd) then
                    count = count +1
                endif
            enddo

            this%nHighDensityOffspring = count
        endif

    end function countHighDensityOffspring


    !---------------------------------------------------------------------------
    !> @brief adds index for family to the individual's families linked list
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine addFamily(this, int)
        class(individual) :: this
        integer :: int !< the input integer to add to the linked list

        call this%families%list_add(int)


    end subroutine addFamily

    !---------------------------------------------------------------------------
    !> @brief initialises probabilities arrays for both the genotype and the phase
    !> @details this is done  from the genoptype and phase structure
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine getProbabilitiesFromOwnGenotypeAndPhase(this)
        use iso_fortran_env
        class(individual) :: this
        

        this%genotypeProbabilities = this%individualGenotype%toIntegerArray()
        this%phaseProbabilities(1,:) = this%IndividualPhase(1)%toIntegerArrayWithErrors()
        this%phaseProbabilities(2,:) = this%IndividualPhase(2)%toIntegerArrayWithErrors()

    end subroutine getProbabilitiesFromOwnGenotypeAndPhase

    !---------------------------------------------------------------------------
    !> @brief arrays for plant imputation
    !> @author  serap gonen serap.gonen@roslin.ed.ac.uk
    !> @date    September 08, 2017
    !---------------------------------------------------------------------------

    subroutine initialisePlantArrays(this, nSnpSubset, nSnpAll)
        use constantModule, only: MissingPlantArrayCode

        class(individual) :: this
        integer, intent(in) :: nSnpSubset, nSnpAll
        
        ! allocate(this%individualGenotypeSubset)
        allocate(this%individualPhaseSubset(2))
        allocate(this%ParentalHaplotypes(2,nSnpSubset))
        allocate(this%StrandBreakArray(nSnpSubset))
        allocate(this%PartialInformativeHaplotype(2,nSnpSubset))
        allocate(this%ParentInferredSwitchCount(2))
        allocate(this%ParentalHaplotypesFull(2,nSnpAll))
        allocate(this%PartialInformativeHaplotypeFull(2,nSnpAll))
        call this%individualGenotypeSubset%newGenotypeMissing(nSnpSubset)
        call this%individualPhaseSubset(1)%newHaplotypeMissing(nSnpSubset)
        call this%individualPhaseSubset(2)%newHaplotypeMissing(nSnpSubset)
        this%ParentalHaplotypes = MissingPlantArrayCode
        this%StrandBreakArray = 0
        this%ParentInferredSwitchCount = 0
        this%PartialInformativeHaplotype = MissingPlantArrayCode
        this%ParentalHaplotypesFull = MissingPlantArrayCode
        this%PartialInformativeHaplotypeFull =MissingPlantArrayCode
    end subroutine initialisePlantArrays

end module IndividualModule

