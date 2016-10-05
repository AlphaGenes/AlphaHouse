
module PedigreeTable
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
        type(Individual), pointer :: sirePointer
        type(Individual), pointer :: damPointer
        ! integer, allocatable :: OffSprings(:)
        type(individualPointerContainer), allocatable :: OffSprings(:)
        integer :: nOffs  = 0
        logical :: Founder     = .false.
        logical :: Genotyped   = .false.
        logical :: HD          = .false.
        logical :: initialised = .false.
        character(len=:), allocatable :: path
        contains
            procedure :: getIdSireDam => getIdSireDamFromLine
            procedure :: getSireDamByIndex
            procedure :: init => initLine
            procedure :: isInitialised => isInitialisedLine
            procedure :: isGenotyped => isGenotypedLine
            procedure :: SetGenotyped => SetGenotypedLine
            procedure :: isHD => isHDLine
            procedure :: SetHD => SetHDLine
            procedure :: isFounder => isFounderLine
            procedure :: SetFounder => SetFounderLine
            procedure :: GetNumberOffsprings
            procedure :: GetOffsprings
            procedure :: AddOffspring

      ! TODO contains writeIndividualFUNCTION
    end type Individual


    interface operator ( == )
        module procedure compareIndividual
    end interface operator ( == )

contains

    pure function getIdSireDamFromLine(this) result(r)
        class(Individual), intent(in) :: this
        integer :: r(3)
        r(1) = this%id
        r(2) = this%sireID
        r(3) = this%damID
        return
    end function getIdSireDamFromLine

    logical function compareIndividual(l1,l2)
        class(Individual), intent(in) :: l1,l2

        if (l1%id == l2%id .and. l1%sireID == l2%sireID .and. l1%damID == l2%damID) then
            compareIndividual=.true.
        else
            compareIndividual=.false.
        endif

        return
    end function compareIndividual

    pure function getSireDamByIndex(this, index) result(v)
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
        end select
        return
    end function getSireDamByIndex

    subroutine initLine(this, originalID, OldGlobalID, id, sireID, damID, generation, path)
        class(Individual), intent(inout) :: this
        character(*), intent(in) :: originalID
        integer, intent(in) :: OldGlobalID
        integer, intent(in) :: id, sireID, damID
        integer, intent(in) :: generation
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
            this%generation = generation
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
    end subroutine initLine

    function isInitialisedLine(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%initialised
    end function isInitialisedLine

    elemental function isFounderLine(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%Founder
    end function isFounderLine

    subroutine SetFounderLine(this)
        class(Individual), intent(inout) :: this
        this%Founder = .true.
    end subroutine SetFounderLine

    elemental function isGenotypedLine(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%Genotyped
    end function isGenotypedLine


    subroutine SetGenotypedLine(this)
        class(Individual), intent(inout) :: this
        this%Genotyped = .true.
    end subroutine SetGenotypedLine

    elemental function isHDLine(this) result(ans)
        class(Individual), intent(in) :: this
        logical :: ans
        ans = this%HD
    end function isHDLine

    subroutine SetHDLine(this)
        class(Individual), intent(inout) :: this
        this%HD = .true.
    end subroutine SetHDLine

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

    elemental function GetNumberOffsprings(this) result(ans)
        class(Individual), intent(in) :: this
        integer :: ans

        ans = this%nOffs

    end function GetNumberOffsprings


    subroutine GetOffsprings(this, Offsprings)

        type(individualPointerContainer), allocatable :: Offsprings(:)
        class(Individual), intent(in) :: this

        Offsprings = this%Offsprings

    end subroutine GetOffsprings

! subroutine builds offspring information given a vector of individual objects that are sorted by generation (although this does not really matter)
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


end module PedigreeTable
