
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     IndividualHelperModule.f90
!
! DESCRIPTION:
!> @brief    Module providing aditional functionality of IndividualModule.f90
!> seperate module due to fortrans circular dependency issues
!
!> @details  Module holds functions to get full, half, and all sibs
!
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!
!> @date     November 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-11-26 Dwilson - Initial Version

!-------------------------------------------------------------------------------

module IndividualHelperModule

    use IndividualModule

    ! procedures
    public :: getFullSibs, getSibs, getOnlyHalfSibs

    contains


    !---------------------------------------------------------------------------
    !> @brief Returns linked list of Full sibs of animal passed in
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @return linked list of full sibs
    !---------------------------------------------------------------------------
    function getFullSibs(indiv) result(res)
        use IndividualLinkedListModule
        class(individual) :: indiv
        type(IndividualLinkedList) :: res
        integer :: i

        if (associated(indiv%sirePointer) .and. associated(indiv%damPointer)) then
            ! loop through sire offsprings

            do i=1, indiv%sirePointer%nOffs
                if (indiv%sirePointer%OffSprings(i)%p%id == indiv%id) cycle
                if(indiv%sirePointer%OffSprings(i)%p%damId == indiv%damID) then
                    call res%list_add(indiv%sirePointer%OffSprings(i)%p)
                endif
            enddo

        endif


    end function getFullSibs


   !---------------------------------------------------------------------------
    !> @brief Returns linked list of Full and half sibs of animal passed in
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @return linked list of full and half sibs
    !---------------------------------------------------------------------------
    function getSibs(indiv) result(res)
        use IndividualLinkedListModule
        use HashModule
        class(individual) :: indiv
        type(IndividualLinkedList) :: res
        type(DictStructure) :: dict
        integer :: i, sireNum, damNum
        integer(kind=int64) ::tmpSize

        sireNum = 0
        damNum = 0

        if (associated(indiv%sirePointer) .and. associated(indiv%damPointer)) then
            ! loop through sire offsprings

            tmpSize = indiv%sirePointer%nOffs + indiv%damPointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = indiv%sirePointer%nOffs
            damNum = indiv%damPointer%nOffs

        else if (associated(indiv%sirePointer)) then
            tmpSize = indiv%sirePointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = indiv%sirePointer%nOffs
        else if (associated(indiv%damPointer)) then
            tmpSize = indiv%damPointer%nOffs
            dict = DictStructure(tmpSize)
            damNum = indiv%damPointer%nOffs
        else
            return !NO parents so return empty list
        endif

        do i=1, sireNum
            if (indiv%sirePointer%OffSprings(i)%p%id == indiv%id) cycle
            call res%list_add(indiv%sirePointer%OffSprings(i)%p)
            call dict%addKey(indiv%sirePointer%OffSprings(i)%p%originalID,1)
        enddo

        do i=1, damNum
            if (indiv%damPointer%OffSprings(i)%p%id == indiv%id) cycle 
            if(.not. dict%hasKey(indiv%damPointer%OffSprings(i)%p%originalID)) then
                call res%list_add(indiv%damPointer%OffSprings(i)%p)
            endif  
        enddo
    

    end function getSibs


    !---------------------------------------------------------------------------
    !> @brief Returns linked list of half sibs of animal passed in
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @return linked list of half sibs
    !---------------------------------------------------------------------------
    function getOnlyHalfSibs(indiv) result(res)

        use IndividualLinkedListModule
        use HashModule
        implicit none
        class(individual) :: indiv
        type(IndividualLinkedList) :: res
        type(DictStructure) :: dict
        integer :: i, sireNum, damNum
        integer(kind=int64) ::tmpSize

        sireNum = 0
        damNum = 0

        if (associated(indiv%sirePointer) .and. associated(indiv%damPointer)) then
            ! loop through sire offsprings

            tmpSize = indiv%sirePointer%nOffs + indiv%damPointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = indiv%sirePointer%nOffs
            damNum = indiv%damPointer%nOffs

        else if (associated(indiv%sirePointer)) then
            tmpSize = indiv%sirePointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = indiv%sirePointer%nOffs
        else if (associated(indiv%damPointer)) then
            tmpSize = indiv%damPointer%nOffs
            dict = DictStructure(tmpSize)
            damNum = indiv%damPointer%nOffs
        else
            return !NO parents so return empty list
        endif

        do i=1, sireNum
            if (indiv%sirePointer%OffSprings(i)%p%id == indiv%id) cycle
            call res%list_add(indiv%sirePointer%OffSprings(i)%p)
            call dict%addKey(indiv%sirePointer%OffSprings(i)%p%originalID,1)
        enddo

        do i=1, damNum
            if (indiv%damPointer%OffSprings(i)%p%id == indiv%id) cycle 
            if(.not. dict%hasKey(indiv%damPointer%OffSprings(i)%p%originalID)) then
                call res%list_add(indiv%sirePointer%OffSprings(i)%p)
            else !as we only want half sibs, remove it
                call res%list_remove(indiv%sirePointer%OffSprings(i)%p)
            endif  
        enddo
    

    end function getOnlyHalfSibs
end module IndividualHelperModule

