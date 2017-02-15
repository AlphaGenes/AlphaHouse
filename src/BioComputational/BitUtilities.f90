!###############################################################################
!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     BitUtilities.f90
!
! DESCRIPTION:
!> @brief    Algorthms for doing bit operations
!
!> @details  
!
!> @author   Daniel Money (I believe)
!
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
!-------------------------------------------------------------------------------
module BitUtilities
  use ISO_Fortran_Env, only: int64
  
  contains
  
    !> @brief Counts the number of positive bits in an array 
    !> @details Takes in an array of 64 bit integers and counts the number of bits in the array which are set to one.
  function bitCount(bits) result(c)
    integer(kind=int64), dimension(:), intent(in) :: bits
    
    integer :: c
    
    integer :: i
    
    c = 0
    do i = 1, size(bits)
      c = c + POPCNT(bits(i))
    end do
  end function bitCount
  
  !>@brief Converts an array of integers (x64) to bits
  !> @details What a confusing name.   But this takes in an array of 64 bit integers and converts it to a bit array.
  function bitToIntegerArray(b) result(array)
    integer(kind=int64), dimension(:), intent(in) :: b
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(size(b) * 64))
    
    cursection = 1
    curpos = 1
    do i = 1, size(b) * 64
      if (btest(b(cursection),curpos)) then
	array(i) = 1
      else
	array(i) = 0
      end if
      
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function bitToIntegerArray
  
end module BitUtilities
