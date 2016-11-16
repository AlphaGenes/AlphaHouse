
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
    implicit none

    integer, parameter :: OFFSPRINGTHRESHOLD = 100
    public :: BuildOffspringInfortmation,Individual,individualPointerContainer,operator ( == )
    
    private

    ! This type is required to have an array of pointers
    type IndividualPointerContainer
        type(Individual), pointer :: p
    end type IndividualPointerContainer

    type Individual
        character(len=:), allocatable :: originalID
        integer :: generation
        integer :: OldGlobalID
        integer :: id
        integer :: sireID
        integer :: damID
        integer(kind=1) :: gender 
        type(Individual), pointer :: sirePointer
        type(Individual), pointer :: damPointer
        type(individualPointerContainer), allocatable :: OffSprings(:)
        integer :: nOffs  = 0
        logical :: Founder     = .false.
        logical :: Genotyped   = .false.
        logical :: HD          = .false.
        logical :: initialised = .false.
        character(len=:), allocatable :: path
        contains
            procedure :: getIdSireDamArrayFormat
            procedure :: getSireDamByIndex
            procedure :: init => initIndividual
            procedure :: isInitialised
            procedure :: isGenotyped
            procedure :: SetAsGenotyped
            procedure :: isHD
            procedure :: SetHD
            procedure :: isFounder
            procedure :: SetObjectAsFounder
            procedure :: GetNumberOffsprings
            procedure :: GetOffsprings
            procedure :: AddOffspring

      ! TODO contains writeIndividualFUNCTION
    end type Individual

    interface operator ( == )
        module procedure compareIndividual
    end interface operator ( == )

contains


     !---------------------------------------------------------------------------
    !> @brief   returns an array where index 1 is the individuals id,
    !> index 2 is the sire id, and index 3 is the dam ID. 
    !> THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    pure function getIdSireDamArrayFormat(this) result(r)
        class(Individual), intent(in) :: this
        integer :: r(3)
        r(1) = this%id
        r(2) = this%sireID
        r(3) = this%damID
        return
    end function getIdSireDamArrayFormat



     !---------------------------------------------------------------------------
    !> @brief Returns true if individuals are equal, false otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    logical function compareIndividual(l1,l2)
        class(Individual), intent(in) :: l1,l2

        if (l1%id == l2%id .and. l1%sireID == l2%sireID .and. l1%damID == l2%damID) then
            compareIndividual=.true.
        else
            compareIndividual=.false.
        endif

        return
    end function compareIndividual

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
        integer, intent(in) :: index
        integer :: v
        select case (index)
            case(1)
                v = this%id
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
    !> @brief Constructor for siredam class.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine initIndividual(this, originalID, OldGlobalID, id, sireID, damID, generation,gender, path)
        class(Individual), intent(inout) :: this
        character(*), intent(in) :: originalID
        integer, intent(in) :: OldGlobalID
        integer, intent(in) :: id, sireID, damID
        integer, intent(in), Optional :: generation
        integer(kind=1), intent(in), Optional :: gender
        character(*),intent(in), Optional :: path
        character(len=512) :: tempPath

        if (.not. this%initialised) then
            allocate(character(len=len(originalID)) ::this%originalID)
            allocate(this%OffSprings(OFFSPRINGTHRESHOLD))
            this%originalID = originalID
            this%id = id
            this%OldGlobalID = OldGlobalID
            this%sireID = sireID
            this%damID = damID
            if (present(generation)) then
                this%generation = generation
            endif
            if (present(gender)) then
                this%gender = gender
            endif
            if (present(Path)) then
                allocate(character(len=len(path)) ::this%path)
                this%path = path
            else
                call getcwd(tempPath)
                this%path = tempPath
            endif
            if (sireID==0 .and. damID==0) then
                this%Founder = .true.
            end if
            this%initialised = .true.
        endif
    end subroutine initIndividual


    !---------------------------------------------------------------------------
    !> @brief boolean function returning true if object is initialised
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !> @return .True. if object is initialised
    !---------------------------------------------------------------------------
    elemental function isInitialised(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%initialised
    end function isInitialised

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
    !> @brief Sets the individual to be genotyped.
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine SetAsGenotyped(this)
        class(Individual), intent(inout) :: this
        this%Genotyped = .true.
    end subroutine SetAsGenotyped

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
            call move_alloc(this%OffSprings,tmp)
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
    !> @brief builds offspring information given a vector of individual objects that are sorted by generation (although this does not really matter)
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine BuildOffspringInfortmation(individuals)
        class(Individual),target, dimension(:), allocatable, intent(inout) :: individuals
        integer :: i, tmpSire, tmpDam

        do i=size(individuals),1,-1  ! start at the end of sorted array and build backwards 
            tmpSire = individuals(i)%sireID
            tmpDam = individuals(i)%damID

            if (tmpSire /= 0) then
                individuals(i)%sirePointer => individuals(tmpSire)
                call individuals(tmpSire)%AddOffspring(individuals(i))

            endif

            if (tmpDam /= 0) then
                individuals(i)%damPointer => individuals(tmpDam)
                call individuals(tmpDam)%AddOffspring(individuals(i))
            endif

        enddo
    end subroutine BuildOffspringInfortmation


end module IndividualModule
