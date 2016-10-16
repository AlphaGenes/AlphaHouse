program main
use inputfile
implicit none

 character(len=:), allocatable:: fileName

 fileName = getInputFile()

 write(*,*) fileName, len(fileName)

!read(*, *) file
!if (trim(file) =="d") then
!Call initialise()
!else
!call initialise(file)
!end if
!write(*,*) inputFileName, len(inputFileName)

end Program

