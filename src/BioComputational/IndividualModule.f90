
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
    implicit none

    public :: Individual,individualPointerContainer,operator ( == ),compareIndividual
    
    private

    ! This type is required to have an array of pointers
    type IndividualPointerContainer
        type(Individual), pointer :: p
    end type IndividualPointerContainer

    type Individual

        character(len=:), allocatable :: originalID
        character(len=:), allocatable :: sireID
        character(len=:), allocatable :: damID
        integer :: generation
        integer :: id
        integer(kind=1) :: gender 
        type(individual), pointer :: sirePointer
        type(individual), pointer :: damPointer
        type(individualPointerContainer), allocatable :: OffSprings(:) !holds array of given size
        integer :: nOffs  = 0 !number of offspring
        logical :: Founder     = .false.
        logical :: Genotyped   = .false.
        logical :: HD          = .false.
        logical :: isDummy     = .false.  ! if this animal is not in the pedigree, this will be true

        type(genotype) :: individualGenotype
        integer(kind=1), allocatable, dimension(:,:) :: phase !where size is the number of sn
        contains
            procedure :: getSireDamByIndex
            procedure :: isGenotyped
            procedure :: isGenotypedNonMissing
            procedure :: setGenotypeArray
            procedure :: isHD
            procedure :: SetHD
            procedure :: getSireId
            procedure :: getDamID
            procedure :: GetNumberOffsprings
            procedure :: GetOffsprings
            procedure :: AddOffspring
            procedure :: setGender
            procedure :: destroyIndividual
            procedure :: setGeneration
            procedure :: getSireDamObjectByIndex
            procedure :: getSireDamNewIDByIndex
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
            procedure :: getPaternalGrandSireRecodedIndexNoDummy
            procedure :: getMaternalGrandSireRecodedIndexNoDummy
            procedure :: getPaternalGrandDamRecodedIndexNoDummy
            procedure :: getMaternalGrandDamRecodedIndexNoDummy
            procedure :: getIntegerVectorOfRecodedIdsNoDummy
            procedure :: resetOffspringInformation
            procedure :: writeIndividual
            generic:: write(formatted) => writeIndividual
    end type Individual

    interface Individual
        module procedure initIndividual
    end interface Individual
    
    interface operator ( == )
        module procedure compareIndividual
    end interface operator ( == )

contains

        !---------------------------------------------------------------------------
    !> @brief Constructor for individual class.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    function initIndividual(originalID,sireIDIn,damIDIn, id, generation,gender) result (this)
        type(Individual) :: this
        character(*), intent(in) :: originalID,sireIDIn,damIDIn
        integer, intent(in), Optional :: generation
        integer, intent(in) :: id
        integer(kind=1), intent(in), Optional :: gender


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

    end function initIndividual

     !---------------------------------------------------------------------------
    !> @brief Deallocates individual object
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine destroyIndividual(this)
        class(Individual) :: this
        deallocate(this%offsprings)
        deallocate(this%originalID)
        deallocate(this%sireID)
        deallocate(this%damID)
    end subroutine destroyIndividual
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
                v = this%sireID
            case(3)
                v = this%damID
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
                isDummyBasedOnIndex = .false.
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
    !> @brief Returns either the individuals id, the sires id or dams id based on
    !> which index is passed.

    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    ! PARAMETERS:
    !> @param[in] index - the index
    !> @return .True. if file exists, otherwise .false.
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
    !> @brief Sets the individual to be genotyped.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine setGenotypeArray(this, geno)
        class(Individual), intent(inout) :: this
        integer(KIND=1), dimension(:), intent(in) :: geno !< One dimensional array of genotype information
        this%Genotyped = .true.
        !TODO this%Genotyped = any(geno == 1 .or. geno == 2 .or. geno == 0)
        ! this%Genotyped = any(geno == 1 .or. geno == 2 .or. geno == 0)
        this%individualGenotype = Genotype(Geno)
    end subroutine setGenotypeArray

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

    !---------------------------------------------------------------------------
    !> @brief Adds an individual as offspring
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine AddOffspring(this, offspringToAdd)
        class(Individual), intent(inout) :: this
        class(Individual),target, intent(in) :: offspringToAdd
        type(individualPointerContainer), allocatable :: tmp(:)
        this%nOffs = this%nOffs + 1
        
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






end module IndividualModule

