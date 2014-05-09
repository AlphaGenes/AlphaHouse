subroutine ConvertCharacterToLowerCase(LWordC)

implicit none

integer :: i,Del
character(len=2000) :: LWordC

Del = iachar('a')-iachar('A')
do i = 1, len_trim(LWordC)                                                   
    if (lge(LWordC(i:i),'A') .and. lle(LWordC(i:i),'Z')) then                  
        LWordC(i:i) = achar(iachar(LWordC(i:i)) + Del)                         
    endif
enddo


end subroutine ConvertCharacterToLowerCase



