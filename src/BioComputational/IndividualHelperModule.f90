
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
    public :: getFullSibs, getSibs, getOnlyHalfSibs, getMates, getAncestors

    contains


    !---------------------------------------------------------------------------
    !> @brief Returns linked list of Full sibs of animal passed in
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @return linked list of full sibs
    !---------------------------------------------------------------------------
    function getFullSibs(indiv) result(res)
        use IndividualLinkedListModule
        class(individual) :: indiv !< individual to get full sibs of
        type(IndividualLinkedList) :: res !< list of full sibs
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
        class(individual) :: indiv !< individual to get full and half sibs on
        type(IndividualLinkedList) :: res !< linked list of full and half sibs
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
        class(individual) :: indiv !< individual to get half sibs of
        type(IndividualLinkedList) :: res !< linked list of half sibs
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



      !---------------------------------------------------------------------------
    !> @brief Returns linked list of individuals that are shared between two parents
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @param[in] type(individual) parent1, parent2
    !> @return linked list of shared children between input parents
    !---------------------------------------------------------------------------
    function getSharedKids(individualOne, IndividualTwo) result(res)
        use IndividualModule
        use IndividualLinkedListModule
        type(Individual),target, intent(in) :: individualOne, IndividualTwo !< individuals to compare kids
        type(IndividualLinkedList) :: res !< linked list of shared kids between both parents
        integer ::  i
        do i=1, individualOne%nOffs

            if (associated(individualOne%offsprings(i)%p%sirePointer, IndividualOne)) then
                if (associated(individualOne%offsprings(i)%p%damPointer, IndividualTwo)) then
                    call res%list_add(individualOne%offsprings(i)%p)
                endif
            else if (associated(individualOne%offsprings(i)%p%damPointer, IndividualOne)) then
                if (associated(individualOne%offsprings(i)%p%sirePointer, IndividualTwo)) then
                    call res%list_add(individualOne%offsprings(i)%p)
                endif
            endif  
        enddo


    end function getSharedKids


    !---------------------------------------------------------------------------
    !> @brief Returns linked list of mates that an individual has had 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !> @param[in] type(individual) parent1
    !> @return linked list of mates of parent
    !---------------------------------------------------------------------------
    function getMates(individualOne) result(res)
        use IndividualModule
        use IndividualLinkedListModule
        type(Individual),target, intent(in) :: individualOne !< individual to get mates of
        type(IndividualLinkedList) :: res !< linked list of individual's mates
        integer ::  i
        do i=1, individualOne%nOffs

            if (individualOne%offsprings(i)%p%sirePointer == IndividualOne) then
                if (associated(individualOne%offsprings(i)%p%damPointer)) then
                    if (.not. res%contains(individualOne%offsprings(i)%p%damPointer)) then
                        call res%list_add(individualOne%offsprings(i)%p%damPointer)
                    endif
                endif
            else if (individualOne%offsprings(i)%p%damPointer == IndividualOne) then
                if (associated(individualOne%offsprings(i)%p%sirePointer)) then
                    if (.not. res%contains(individualOne%offsprings(i)%p%sirePointer)) then
                        call res%list_add(individualOne%offsprings(i)%p%sirePointer)
                    endif
                endif
            endif  
        enddo

        print *,res%length

    end function getMates

    !---------------------------------------------------------------------------
    !> @brief Returns true if animals are mates, false otherwise 
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !---------------------------------------------------------------------------
    function areMates(ind1, ind2) result(res)

        use IndividualModule
        use IndividualLinkedListModule

        type(Individual), intent(in) :: ind1, ind2 !< animals to determine if mates

        type(IndividualLinkedList) :: mates1
        logical :: res !< true if animals are mates

        res = .false.
        mates1 = getMates(ind1)
        res = mates1%contains(ind2)

    end function areMates

        !---------------------------------------------------------------------------
    !> @brief Returns true if animals share mates, .false otherwise
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !---------------------------------------------------------------------------
    function doShareMates(ind1, ind2) result(res)
        use IndividualModule
        use IndividualLinkedListModule

        type(Individual), intent(in) :: ind1, ind2
        type(IndividualLinkedListNode),pointer :: tmpInd

        type(IndividualLinkedList) :: mates1, mates2
        logical :: res
        integer :: i
        res = .false.
        mates1 = getMates(ind1)
        mates2 = getMates(ind2)

        if (.not. associated(mates1%first) .or. .not. associated(mates2%first)) then
            res = .false.
            return
        endif
        tmpInd => mates1%first
        do i=1, mates1%length
            if (mates2%contains(tmpInd%item)) then
                res = .true.
                return
            endif
            tmpInd => tmpInd%next
        enddo

    end function doShareMates

    !---------------------------------------------------------------------------
    !> @brief creates list of ancestors of given animal
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    November 26, 2016
    !---------------------------------------------------------------------------
    recursive subroutine getAncestors(ind, disCount, res, distList, cap)
    ! TODO write tests
        use IndividualModule
        use IndividualLinkedListModule
        use IntegerLinkedListModule
        type(Individual), intent(in) :: ind !< individual to get ancestors of
        type(integerLinkedList),intent(inout) ::distList !< list of genetic distance to indivs in res
        type(IndividualLinkedList),intent(inout) :: res
        integer, intent(in) :: disCount !< count of current distance
        integer,intent(in) :: cap !< should be set to max int if no cap is desired
        if (disCount >= cap) then
            Return
        endif

        if (associated(ind%sirePointer)) then
            if (.not. res%contains(ind%sirePointer)) then
                call res%list_add(ind%sirePointer)
                call distList%list_add(disCount+1)
            endif
            call getAncestors(ind%sirePointer,disCount+1, res, distList, cap)
        endif
        if (associated(ind%damPointer)) then
            if (.not. res%contains(ind%damPointer)) then
                call res%list_add(ind%damPointer)
                call distList%list_add(disCount+1)
            endif
            call getAncestors(ind%damPointer,disCount+1,res,distList, cap)
        endif
end subroutine getAncestors




!---------------------------------------------------------------------------
!> @brief calculates the minimum genetic distance between two animals
!> If no relationship is found a distance of MAXINT32 is returned.
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    November 26, 2016
!---------------------------------------------------------------------------
integer function calcGenDistance(ind1, ind2)
    use IndividualModule
    use IndividualLinkedListModule
    use IntegerLinkedListModule
    use ConstantModule

    type(individual), intent(in) :: ind1,ind2 !< individuals to calculate genetic distance of
    type(IndividualLinkedListNode) :: tmp1, tmp2
    type(integerLinkedList)::dist1,dist2
    type(IndividualLinkedList):: an1,an2
    integer :: val,i,h
    calcGenDistance = MAXINT32 !< init to a large num
    
    call getAncestors(ind1, 0, an1, dist1, MAXINT32) !< take ancestors of both animals

    if (associated(an1%first)) then

        tmp1 = an1%first

        ! loop through ancestors of ind1 and see if ind2 is an ancestor
        do i=1, an1%length
            
            if ( tmp1%item == ind2) then

                calcGenDistance = dist1%list_get_nth(i)
                return
            endif
            tmp1 = tmp1%next
        enddo
    endif


    call getAncestors(ind2, 0, an2, dist2, MAXINT32)

    if (associated(an2%first)) then
        tmp2 = an2%first
    ! loop through ancestors of ind1 and see if ind2 is an ancestor
        do i=1, an2%length
            
            if (tmp2%item == ind1) then

                calcGenDistance = dist2%list_get_nth(i)
                return
            endif
            tmp2 = tmp2%next
        enddo
    endif 

    ! if they are not ancestors and either list is empty then animals aren't connected
    if (.not. associated(an1%first) .or. .not. associated(an2%first)) then
        calcGenDistance = MAXINT32
        return
    end if
    tmp1 = an1%first

    ! loop through animals of both animals and see where they interest
    do i=1, an1%length
        tmp2 = an2%first
        do h=1, an2%length

            if (tmp1%item == tmp2%item) then
                val = dist1%list_get_nth(i) + dist2%list_get_nth(h)
                if (calcGenDistance > val ) then
                    calcGenDistance = val
                endif
            endif
            tmp2 = tmp2%next
        enddo
        tmp1 = tmp1%next
    enddo


end function calcGenDistance

end module IndividualHelperModule

