!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file    Input.f90 
!
! DESCRIPTION:
!> @brief    Holds a subroutine to get spec file
!
!> @details   This holds two subroutines to get the input spec file.   The input spec file is by default set to be the name of the
!>program (gotten from get_command_arg(0)) appended by "Spec.txt", i.e. AlphaAnalyseSpec.txt.   You can change the name of the input
!>file by passing the flag "-f <fileNameWanted>" into the program from the command line.
!> 
!> It is also possible to set what the input filename should be in your program by using setInputFile.   This is not recommended
!though, as it means that trying to change the input file name via -f will stop working.
!!           
!
!> @author   Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
!
!> @date     January 24, 2017
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
!
!-------------------------------------------------------------------------------
module inputfile
  use iso_fortran_env
  implicit none

  private
  public:: getInputFile, setInputFile

  ! inputFileName, initialise

  integer,parameter:: fileTypeLength=2, fileLength=1000
  character(len=:), allocatable :: inputFileName
  character(len=fileTypeLength), parameter :: inputType="-f"
  character(len=2)::test

contains

!>@brief Sets the input file name
!> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
subroutine setInputFile(fileNameIn)
  character(len=*):: fileNameIn !< the filename that you want to set

  inputFileName = fileNameIn

end subroutine setInputFile

!> @brief A function to return the input file name
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
function getInputFile() result (inputFile)
  character(len=:), allocatable:: inputFile !< the name of the file send out
  if (.not. allocated(inputFileName)) then
    call initialise()
  end if
  inputFile=inputFileName
end function getInputFile

!> @brief A subroutine that will set up the input file name
!> @details This subroutine will either set the input file name to be the name of the program followed by "Spec.txt", unless the -f
!> flag has been used to send in a filename via the command line.   If the -f flag was used, then that file name is used instead.
subroutine initialise(filePathIn)
  character(len=*), intent(in), optional::filePathIn
  character(len=fileLength):: currentWorkingDirectory
  integer(kind=int32):: commandArgument
  if (present(filePathIn)) then
    inputFileName=trim(filePathIn)
  else
    commandArgument = getFileCommandArgumentNumber()
    call get_command_argument(commandArgument,currentWorkingDirectory)
    if (commandArgument==0) then
    write(currentWorkingDirectory,*) trim(currentWorkingDirectory(1:len(currentWorkingDirectory))),"Spec.txt"
    end if
    call checkInputFileExists(currentWorkingDirectory)
  end if


end subroutine initialise

function getFileCommandArgumentNumber() result(argumentNumber)
  integer(int32)::commandLineArguments
  integer(int32):: argumentNumber

  integer(int32)::i
  argumentNumber=0

  commandLineArguments = command_argument_count()
  do i = 1, commandLineArguments
    if (checkIfFileInput(i)) then
      argumentNumber = i+1
    end if
  end do
end function getFileCommandArgumentNumber

function checkIfFileInput(i) result(isInputFileTagAndCommandExists)
  logical isInputFileTagAndCommandExists, isInputFileTag
  integer(int32) :: i, nextCommandArgumentExists
  integer(int32) :: ArgumentLength
  character(len=fileTypeLength)::checkIfInputFile

  isInputFileTagAndCommandExists = .false.
  call get_command_argument(i,checkIfInputFile, length=ArgumentLength)
  isInputFileTag = (checkIfInputFile==inputType) .and. (ArgumentLength==2)
  if (isInputFileTag) then
    call get_command_argument(i+1,status=nextCommandArgumentExists)
    if (nextCommandArgumentExists==0) then
      isInputFileTagAndCommandExists =.true.
    else if ((nextCommandArgumentExists>0)) then
      write(*,*) "Input file retrieval failed"
      call writeOutInstructions()
      stop 100
    else if (nextCommandArgumentExists==-1) then
      write(*,*) "Input file length to long. Length is limited to 1000 characters"
      stop 101
    else
      stop 102
    end if
  end if
end function checkIfFileInput


subroutine checkInputFileExists(defaultInputFilePath) 
  character(len=*)::defaultInputFilePath
  logical::defaultInputFilePathExists

  Inquire(file=trim(defaultInputFilePath), exist=defaultInputFilePathExists)

  if (defaultInputFilePathExists) then
    inputFileName=trim(defaultInputFilePath)
  else
    write(*,"(A,A,A)") "InputFile " ,trim(defaultInputFilePath), " doesn't exist"
    stop 1
  end if
end subroutine checkInputFileExists

subroutine writeOutInstructions()
  character(len = 1000):: executableName
  character(len=:), allocatable:: defaultPostfix

  defaultPostfix = "Spec.txt"

  call get_command_argument(0, executableName)
  write(*,*) "Tell me the input file by using ", inputType, " YourFileName."
  write(*,*) 
  write(*,*) "Alternatively don't pass me ", inputType, "and I will use the default inputFile. & 
    &This is the name of the exacutable followed by ", defaultPostfix, "The name of  &
    &this exacutable is ", trim(executableName(3:len(executableName))) , "and"
  write(*,*) "the default filename is ", trim(executableName(3:len(executableName))),defaultPostfix
end subroutine writeOutInstructions

end module inputfile
