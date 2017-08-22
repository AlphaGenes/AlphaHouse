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
!> @author   Daniel Money
!
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
!-------------------------------------------------------------------------------
module BitUtilities
  use ISO_Fortran_Env, only: int64
  
  contains
  
  !-------------------------------------------------------------------------------
    !> @brief Counts the number of bits set across all integers in the array 
    !> @details Takes in an array of 64 bit integers and counts the number of bits in the array which are set to one.
    !-------------------------------------------------------------------------------
  function bitCount(bits) result(c)
    integer(kind=int64), dimension(:), intent(in) :: bits
    
    integer :: c
    
    integer :: i
    
    c = 0
    do i = 1, size(bits)
      c = c + POPCNT(bits(i))
    end do
  end function bitCount




  !-------------------------------------------------------------------------------
  !>@brief returns true if bit is set at position
  !-------------------------------------------------------------------------------
  function testBit(bits, pos) result (res)
    integer(kind=int64), dimension(:), intent(in) :: bits !< input array of bits 

    integer, intent(in) :: pos !< position 

    logical :: res !< true if position is set

    integer :: cursection, curpos


    if (pos == 0) then 
      res = .false.
      return
    endif
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1

    if (btest(bits(cursection),curpos)) then
        res = .true.

    else
        res = .false.
    end if
end function testBit
  
  !-------------------------------------------------------------------------------
  !>@brief Converts bits to an integer array
  !> @details Converts each bit in an integer to a zero or one in a integer array.  Can do this for multiple integers (hence the
  !> array input).
  !-------------------------------------------------------------------------------
  function bitToIntegerArray(b) result(array)
    integer(kind=int64), dimension(:), intent(in) :: b
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(size(b) * 64))
    
    cursection = 1
    curpos = 0
    do i = 1, size(b) * 64
      if (btest(b(cursection),curpos)) then
	array(i) = 1
      else
	array(i) = 0
      end if
      
      curpos = curpos + 1
      if (curpos == 64) then
	curpos = 0
	cursection = cursection + 1
      end if
    end do
  end function bitToIntegerArray
  
end module BitUtilities
