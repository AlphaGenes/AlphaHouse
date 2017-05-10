module SortingModule
    implicit none

    character(:), dimension(:), pointer :: sortArray


    type group
	    integer :: order    
	    real :: value       
	end type group
    contains
	 
	recursive subroutine QSortNew(a,na)
	 
		integer, intent(in) :: nA
		type (group), dimension(nA), intent(inout) :: A
		 
		integer :: left, right, marker
		real :: random, pivot
		type(group) :: temp
	 
	    if (nA > 1) then
	 
	        call random_number(random)
	        pivot = A(int(random*real(nA-1))+1)%value   
	        left = 0
	        right = nA + 1
	 
	        do while (left < right)
	            right = right - 1
	            do while (A(right)%value > pivot)
	                right = right - 1
	            enddo
	            left = left + 1
	            do while (A(left)%value < pivot)
	                left = left + 1
	            enddo
	            if (left < right) then
	                temp = A(left)
	                A(left) = A(right)
	                A(right) = temp
	            endif
	        enddo
	 
	        if (left == right) then
	            marker = left + 1
	        else
	            marker = left
	        endif
	 
	        call QSortNew(A(:marker-1),marker-1)
	        call QSortNew(A(marker:),nA-marker+1)
	 
	    endif
	 
	end subroutine QSortNew

    subroutine sortWithIndex(array, indexes)
    character(*), dimension(:), intent(inout), target :: array
    integer, dimension(size(array)), intent(out) :: indexes

    character(len(array)), dimension(size(array)) :: tempArray
    integer :: i
    integer(8) :: as, es

    sortArray => array
    as = size(indexes)
    es = sizeof(indexes(1))

    do i = 1, size(array)
        indexes(i) = i
    end do

    call qsort(indexes, as, es, sortCompare)

    do i = 1, size(array)
        tempArray(i) = array(indexes(i))
    end do

    array = tempArray
    end subroutine sortWithIndex

    function compare(a, b) result(i)
    character(*), intent(in) :: a, b
    integer(2) :: i

    if (a < b) then
        i = -1
    end if
    if (a > b) then
        i = 1
    end if
    if (a == b) then
        i = 0
    end if
    !    print *, a, b, i
    end function compare

    function sortCompare(a, b) result(i)
    integer, intent(in) :: a, b
    integer(2) :: i

    i = compare(sortArray(a), sortArray(b))
    end function sortCompare
    end module SortingModule

