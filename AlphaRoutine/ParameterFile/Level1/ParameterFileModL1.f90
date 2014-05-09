module ParameterFileModL1

implicit none

contains


subroutine CheckKeyWordPresence(LKeyWordC,LFileNameC,LSearchStartI,LSeachEndI,LNumOfLineI,LWordPresentL,LSubModeI)

implicit none

integer :: i
integer :: LnLineI,LNumOfLineI,LTmpI,LUnitNumberI,LSearchStartI,LSeachEndI,LSubModeI
character(len=2000) :: LDumC,LKeyWordC,LFileNameC,LMatchWordC
logical :: LWordPresentL

open(unit=1,file=trim(LFileNameC),status='old')

LWordPresentL=.false.

if (LSubModeI==1) then

	LnLineI=0
	do

		read (1,*,iostat=LTmpI) LMatchWordC
		LnLineI=LnLineI+1
		if (LTmpI/=0) then
			LnLineI=LnLineI-1
			exit
		else
			if (trim(LMatchWordC)==trim(LKeyWordC)) then	
				LWordPresentL=.true.
				LNumOfLineI=LnLineI
			endif
		endif

	enddo
	
endif

if (LSubModeI==2) then

	do i=1,LSearchStartI-1
		read (1,*) LDumC
	enddo

	do i=LSearchStartI,LSeachEndI
		read (1,*,iostat=LTmpI) LMatchWordC
		LnLineI=i
		if (trim(LMatchWordC)==trim(LKeyWordC)) then	
			LWordPresentL=.true.
			LNumOfLineI=LnLineI
			exit
		endif

	enddo

endif	

close(1)

end subroutine CheckKeyWordPresence

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


end module ParameterFileModL1