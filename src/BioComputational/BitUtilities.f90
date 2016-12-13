module BitUtilities
  
  contains
  
  function bitCount(bits) result(c)
    integer(kind=8), dimension(:), intent(in) :: bits
    
    integer :: c
    
    integer :: i
    
    c = 0
    do i = 1, size(bits)
      c = c + POPCNT(bits(i))
    end do
  end function bitCount
  
  function bitToIntegerArray(b) result(array)
    integer(kind=8), dimension(:), intent(in) :: b
    
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