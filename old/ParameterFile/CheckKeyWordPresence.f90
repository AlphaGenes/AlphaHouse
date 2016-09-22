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
