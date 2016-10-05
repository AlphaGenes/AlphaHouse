
module PedigreeTable
    implicit none
    public :: pedigreeLine, getIdSireDamFromLine, initLine
    public :: GetNumberOffsprings
    type PedigreeLine
        character(len=:), allocatable :: ID
        integer :: generation
        integer :: OldGlobalID
        integer :: NewId
        integer :: NewSire
        integer :: NewDam
        integer, allocatable :: OffSprings(:)
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

      ! contains writePedigreeLINEFUNCTION
    end type PedigreeLine

    interface operator ( == )
        module procedure comparepedigreeline
    end interface operator ( == )

contains

    pure function getIdSireDamFromLine(this) result(r)
        class(PedigreeLine), intent(in) :: this
        integer :: r(3)
        r(1) = this%NewId
        r(2) = this%NewSire
        r(3) = this%NewDam
        return
    end function getIdSireDamFromLine

    logical function comparePedigreeLine(l1,l2)
        class(pedigreeLine), intent(in) :: l1,l2

        if (l1%newID == l2%newId .and. l1%newSire == l2%newSire .and. l1%newDam == l2%newDam) then
            comparePedigreeLine=.true.
        else
            comparePedigreeLine=.false.
        endif

        return
    end function comparePedigreeLine

    pure function getSireDamByIndex(this, index) result(v)
        class(pedigreeLine), intent(in) :: this
        integer, intent(in) :: index
        integer :: v
        select case (index)
            case(1)
                v = this%NewId
            case(2)
                v = this%NewSire
            case(3)
                v = this%newDam
        end select
        return
    end function getSireDamByIndex

    subroutine initLine(this, ID, OldGlobalID, NewId, NewSire, NewDam, generation, path)
        class(pedigreeLine), intent(inout) :: this
        character(*), intent(in) :: ID
        integer, intent(in) :: OldGlobalID
        integer, intent(in) :: NewId, NewSire, NewDam
        integer, intent(in) :: generation
        character(*),intent(in), Optional :: path
        character(len=512) :: tempPath

        if (.not. this%initialised) then
            allocate(character(len=len(ID)) ::this%ID)
            this%ID = ID
            this%NewId = NewId
            this%OldGlobalID = OldGlobalID
            this%NewSire = NewSire
            this%NewDam = NewDam
            this%generation = generation
            if (present(Path)) then
                allocate(character(len=len(path)) ::this%path)
                this%path = path
            else
                call getcwd(tempPath)
                this%path = tempPath
            endif
            if (NewSire==0 .and. NewDam==0) then
                this%Founder = .true.
            end if
            this%initialised = .true.
        endif
    end subroutine initLine

    function isInitialisedLine(this) result(ans)
        class(pedigreeLine), intent(in) :: this
        logical :: ans
        ans = this%initialised
    end function isInitialisedLine

    elemental function isFounderLine(this) result(ans)
        class(pedigreeLine), intent(in) :: this
        logical :: ans
        ans = this%Founder
    end function isFounderLine

    subroutine SetFounderLine(this)
        class(pedigreeLine), intent(inout) :: this
        this%Founder = .true.
    end subroutine SetFounderLine

    elemental function isGenotypedLine(this) result(ans)
        class(pedigreeLine), intent(in) :: this
        logical :: ans
        ans = this%Genotyped
    end function isGenotypedLine


    subroutine SetGenotypedLine(this)
        class(pedigreeLine), intent(inout) :: this
        this%Genotyped = .true.
    end subroutine SetGenotypedLine

    elemental function isHDLine(this) result(ans)
        class(pedigreeLine), intent(in) :: this
        logical :: ans
        ans = this%HD
    end function isHDLine

    subroutine SetHDLine(this)
        class(pedigreeLine), intent(inout) :: this
        this%HD = .true.
    end subroutine SetHDLine

    subroutine IncreaseOffsprings(this)
        class(pedigreeLine), intent(inout) :: this
        this%nOffs = this%nOffs + 1
    end subroutine IncreaseOffsprings

    subroutine AddOffspring(this, OffIndex)
        class(pedigreeLine), intent(inout) :: this
        integer, intent(in) ::OffIndex

        integer, allocatable :: tmp(:)

        allocate(tmp(this%nOffs))
        tmp(1:size(this%Offsprings)) = this%Offsprings !Save previous array
        call move_alloc(tmp,this%OffSprings)
        this%Offsprings(this%nOffs) = OffIndex
    end subroutine AddOffspring

    subroutine GetNumberOffsprings(LineList, nAnis)
        type(pedigreeLine), intent(inout) :: LineList(0:nAnis)
        integer, intent(in) :: nAnis

        integer :: i

        do i = nAnis, 0, -1
            if (LineList(i)%NewSire /= 0) then
                call IncreaseOffsprings(LineList(LineList(i)%NewSire))
            end if
            if (LineList(i)%NewDam /= 0) then
                call IncreaseOffsprings(LineList(LineList(i)%NewDam))
            end if
        end do
    end subroutine GetNumberOffsprings

    subroutine GetOffsprings(LineList, nAnis)
        type(pedigreeLine), intent(inout) :: LineList(0:nAnis)
        integer, intent(in) :: nAnis

        integer :: i

        do i = nAnis, 0, -1
            if (LineList(i)%NewSire /= 0) then
                call IncreaseOffsprings(LineList(LineList(i)%NewSire))
                call AddOffspring(LineList(LineList(i)%NewSire),LineList(i)%NewId)
            end if
            if (LineList(i)%NewDam /= 0) then
                call IncreaseOffsprings(LineList(LineList(i)%NewDam))
                call AddOffspring(LineList(LineList(i)%NewDam),LineList(i)%NewId)
            end if
        end do
    end subroutine GetOffsprings

end module PedigreeTable
