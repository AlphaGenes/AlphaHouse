program main
use inputfile
implicit none

 character(len=:), allocatable:: fileName

 call getInputFile(filename)

 write(*,*) fileName, len(fileName)

end Program

