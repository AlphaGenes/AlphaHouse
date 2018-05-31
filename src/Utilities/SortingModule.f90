

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     SortingModule.f90
!
! DESCRIPTION:
!> @brief    Module containing sort functions
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 DWilson - Initial Version
!
!-------------------------------------------------------------------------------
module SortingModule
	implicit none

	character(:), dimension(:), pointer :: sortArray

	type group
	integer :: order
	real :: value
end type group

    interface heapsort
        module procedure heapsortR, heapsortI
    end interface heapsort
contains


	subroutine heapsortR(a)
		real, intent(in out) :: a(0:)
		integer :: start, n, bottom
		real :: temp

		n = size(a)
		do start = (n - 2) / 2, 0, -1
			call siftdownR(a, start, n);
		end do

		do bottom = n - 1, 1, -1
			temp = a(0)
			a(0) = a(bottom)
			a(bottom) = temp;
			call siftdownR(a, 0, bottom)
		end do
	end subroutine heapsortR

	subroutine siftdownR(a, start, bottom)

		real, intent(in out) :: a(0:)
		integer, intent(in) :: start, bottom
		integer :: child, root
		real :: temp

		root = start
		do while(root*2 + 1 < bottom)
			child = root * 2 + 1

			if (child + 1 < bottom) then
				if (a(child) < a(child+1)) child = child + 1
			end if

			if (a(root) < a(child)) then
				temp = a(child)
				a(child) = a (root)
				a(root) = temp
				root = child
			else
				return
			end if
		end do

	end subroutine siftdownR

    subroutine heapsortI(a)
        integer, intent(in out) :: a(0:)
        integer :: start, n, bottom
        real :: temp

        n = size(a)
        do start = (n - 2) / 2, 0, -1
            call siftdownI(a, start, n);
        end do

        do bottom = n - 1, 1, -1
            temp = a(0)
            a(0) = a(bottom)
            a(bottom) = temp;
            call siftdownI(a, 0, bottom)
        end do
    end subroutine heapsortI

    subroutine siftdownI(a, start, bottom)

        integer, intent(in out) :: a(0:)
        integer, intent(in) :: start, bottom
        integer :: child, root
        real :: temp

        root = start
        do while(root*2 + 1 < bottom)
            child = root * 2 + 1

            if (child + 1 < bottom) then
                if (a(child) < a(child+1)) child = child + 1
            end if

            if (a(root) < a(child)) then
                temp = a(child)
                a(child) = a (root)
                a(root) = temp
                root = child
            else
                return
            end if
        end do

    end subroutine siftdownI

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

	!---------------------------------------------------------------------------
	!> @brief   Takes a string array and returns a list of sorted indicies
	!> @detail  As an example [c, a, d, b] would return [2, 4, 1, 3].  Item 2
	!>          is the smallest, item 4 the next smallest etc
	!>          Should probably be a function rather than a subroutine.  Not
	!>          sure why it is not.
	!> @author  Daniel Money, daniel.money@roslin.ed.ac.uk
	!> @date    Aug 21, 2017
	!> @return  Array of sorted indexes
	!---------------------------------------------------------------------------
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

	!---------------------------------------------------------------------------
	!> @brief   Compares two strings
	!> @detail  Uses standard Fortran sort order
	!> @author  Daniel Money, daniel.money@roslin.ed.ac.uk
	!> @date    Aug 21, 2017
	!> @return  -1 if first character is smaller, 0 if the same, 1 if second
	!>          integer is smaller.
	!---------------------------------------------------------------------------
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
	end function compare

	!---------------------------------------------------------------------------
	!> @brief   Compares two elements of the sort array
	!> @author  Daniel Money, daniel.money@roslin.ed.ac.uk
	!> @date    Aug 21, 2017
	!> @return  -1 if array(a) is smaller than array(b) , 0 if the same, 1
	!>          otherwiae.
	!---------------------------------------------------------------------------
	function sortCompare(a, b) result(i)
		integer, intent(in) :: a, b
		integer(2) :: i

		i = compare(sortArray(a), sortArray(b))
	end function sortCompare
end module SortingModule


