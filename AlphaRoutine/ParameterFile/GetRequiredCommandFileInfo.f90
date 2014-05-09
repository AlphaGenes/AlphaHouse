subroutine GetRequiredCommandFileInfo(LFileNameC,LNumOfLineI,LNumOfColumnI,LReturnWordC)

implicit none

integer :: i
integer :: LNumOfLineI,LNumOfColumnI
character(len=2000) :: LFileNameC,LDumC,LReturnWordC
character(len=2000),allocatable,dimension(:) :: LWordVecC
logical :: LCommandPresentL

!open(unit=1,file=trim(LFileNameC),status='old')
!
!allocate(LWordVecC(LNumOfColumnI))
!
!do i=1,LNumOfLineI
!	read (1,*) LDumC
!enddo
!
!read (1,*) LWordVecC(:)
!
!close (1)
!
!LReturnWordC=LWordVecC(LNumOfColumnI)	
!

end subroutine GetRequiredCommandFileInfo
