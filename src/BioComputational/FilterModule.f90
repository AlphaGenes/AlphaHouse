module FilterModule
  use constantModule, only : MissingGenotypeCode
  
  contains
  
  function filterSNPmissing(genos, threshold) result(newGenos)
    integer(kind=1), dimension(:,:), intent(in) :: genos
    double precision :: threshold    
    
    integer(kind=1), pointer, dimension(:,:) :: newGenos
    
    logical, dimension(size(genos,2)) :: pass
    integer :: numPass
    integer :: i, j, newi
    integer :: present
    
    numPass = 0    
    do i = 1, size(genos,2)
      present = 0
      do j = 1, size(genos,1)
	if (genos(j, i) /= MissingGenotypeCode) then
	  present = present + 1
	end if
      end do
      pass(i) = ((float(present) / size(genos,1)) >= threshold)
      if (pass(i)) then
	numPass = numPass + 1
      end if
    end do
    
    allocate(newGenos(size(genos,1),numPass))
    
    newi = 0
    
    do i = 1, size(genos,2)
      if (pass(i)) then
	newi = newi + 1
	newGenos(:,newi) = genos(:,i)
      end if
    end do
    
  end function filterSNPmissing

end module FilterModule