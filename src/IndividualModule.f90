
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
    integer, parameter :: OFFSPRINGTHRESHOLD = 150
    public :: Individual,individualPointerContainer,operator ( == )
    
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
        type(individualPointerContainer), allocatable :: OffSprings(:)
        integer :: nOffs  = 0
        logical :: Founder     = .false.
        logical :: Genotyped   = .false.
        logical :: HD          = .false.
        
        contains
            procedure :: getSireDamByIndex
            procedure :: isGenotyped
            procedure :: SetAsGenotyped
            procedure :: isHD
            procedure :: SetHD
            procedure :: isFounder
            procedure :: getSireId
            procedure :: getDamID
            procedure :: SetObjectAsFounder
            procedure :: GetNumberOffsprings
            procedure :: GetOffsprings
            procedure :: AddOffspring
            procedure :: setGender
            procedure :: destroyIndividual
            procedure :: setGeneration
      ! TODO contains writeIndividualFUNCTION
    end type Individual

    interface Individual
        module procedure initIndividual
    end interface Individual
    
    interface operator ( == )
        module procedure compareIndividual
    end interface operator ( == )

contains



     !---------------------------------------------------------------------------
    !> @brief Deallocates individual object
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine destroyIndividual(this)
        class(Individual) :: this
        deallocate(this%OffSprings)
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
    !> @brief Constructor for siredam class.
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
        this%sireId = sireIDIn
        this%damId = damIDIn
        if (present(generation)) then
            this%generation = generation
        endif
        if (present(gender)) then
            this%gender = gender
        endif

    end function initIndividual

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
    !> @param[in] gender (integer kind(1)) 
    !---------------------------------------------------------------------------
    subroutine setGender(this,gender) 
        class(individual) :: this
        integer(kind=1),intent(in) :: gender
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
