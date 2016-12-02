module IndividualHelperModule

    use IndividualModule

    contains

    function getFullSibs(this) result(res)
        use IndividualLinkedListModule
        class(individual) :: this
        type(IndividualLinkedList) :: res
        integer :: i

        if (associated(this%sirePointer) .and. associated(this%damPointer)) then
            ! loop through sire offsprings

            do i=1, this%sirePointer%nOffs
                if (this%sirePointer%OffSprings(i)%p%id == this%id) cycle
                if(this%sirePointer%OffSprings(i)%p%damId == this%damID) then
                    call res%list_add(this%sirePointer%OffSprings(i)%p)
                endif
            enddo

        endif


    end function getFullSibs


! half sibs and fullSibs
    function getSibs(this) result(res)
        use IndividualLinkedListModule
        use HashModule
        class(individual) :: this
        type(IndividualLinkedList) :: res
        type(DictStructure) :: dict
        integer :: i, sireNum, damNum
        integer(kind=int64) ::tmpSize

        sireNum = 0
        damNum = 0

        if (associated(this%sirePointer) .and. associated(this%damPointer)) then
            ! loop through sire offsprings

            tmpSize = this%sirePointer%nOffs + this%damPointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = this%sirePointer%nOffs
            damNum = this%damPointer%nOffs

        else if (associated(this%sirePointer)) then
            tmpSize = this%sirePointer%nOffs
            dict = DictStructure(tmpSize)
            sireNum = this%sirePointer%nOffs
        else if (associated(this%damPointer)) then
            tmpSize = this%damPointer%nOffs
            dict = DictStructure(tmpSize)
            damNum = this%damPointer%nOffs
        else
            return !NO parents so return empty list
        endif

        do i=1, sireNum
            if (this%sirePointer%OffSprings(i)%p%id == this%id) cycle
            call res%list_add(this%sirePointer%OffSprings(i)%p)
            call dict%addKey(this%sirePointer%OffSprings(i)%p%originalID,1)
        enddo

        do i=1, damNum
            if (this%damPointer%OffSprings(i)%p%id == this%id) cycle 
            if(.not. dict%hasKey(this%damPointer%OffSprings(i)%p%originalID)) then
                call res%list_add(this%sirePointer%OffSprings(i)%p)
            endif  
        enddo
    

    end function getSibs


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

