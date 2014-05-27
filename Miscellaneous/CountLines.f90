subroutine CountLines(LFileNameC,LnLineI)
implicit none

character(len=*) :: LFileNameC
character(len=3000) :: LDumC
integer :: LnLineI,f

open (unit=1,file=trim(LFileNameC),status="old")

LnLineI=0
do
	read (1,*,iostat=f) LDumC
	LnLineI=LnLineI+1
	if (f/=0) then
		LnLineI=LnLineI-1
		exit
	endif
enddo

close(1)

end subroutine CountLines
