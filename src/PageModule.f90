module pageModule
  use iso_fortran_env
  use lineModule
  use stringModule
  implicit none

  private
  public:: assignment(=), Page

  type:: Page
    private
    type(Line), dimension(:), allocatable:: lines
    character(len=:), allocatable:: pageName
    contains
      procedure, private :: initPageWithArray
      procedure, private :: initChar
      procedure, private :: initPageWithPage
      procedure, private :: readInInputFile
      procedure, public  :: getLine
      procedure, public  :: getWord => getPageWord
      procedure, public  :: set => readInInputFile
      procedure, public  :: getWithStackTrace
      procedure, public  :: setNumLines
      procedure, public  :: getNumLines
      procedure, public  :: getName
      procedure, public  :: setName
      procedure, public  :: getNumWords => getPageTypeNumWords

  end type

!  interface set
!    module procedure setWithStackTraceChar, setWithStackTraceArbLine
!  end interface set

  interface assignment (=)
    module procedure initChar, initPageWithArray, initPageWithPage
  end interface 
  contains

    function getLine(this, numLine) result (res)
      class(Page), intent(in):: this
      integer, intent(in):: numLine

      character(len=:), allocatable:: res

      integer:: i

      res = ""
      do i = 1, this%getNumWords(numLine)
        res = res // this%getWord(numLine, i)
      end do
    end function getLine

    function getPageTypeNumWords(this, i) result(numWords)
      class(Page), intent(in):: this
      integer, intent(in):: i
      integer:: numWords

      numWords = this%lines(i)%getNumWords()

    end function getPageTypeNumWords

    function getNumLines(this) result(numLines)
      integer(int32):: numLines
      class(Page), intent(in):: this

      numLines = size(this%lines)
    end function getNumLines

    subroutine setNumLines(this, i)
      integer, intent(in):: i
      class(Page), intent(inout):: this

      if (allocated(this%lines)) then
        deallocate(this%lines)
      end if
      allocate(this%lines(i))
    end subroutine setNumLines

    function getName(this) result(nameOut)
      class(Page), intent(in):: this
      character(len=:), allocatable:: nameOut

      if (allocated(this%pageName)) then
        nameOut = this%pageName
      else
        nameOut = "Undefined"
      end if
    end function

    subroutine setName(this, nameIn)
      class(Page), intent(inout):: this
      character(len=*):: nameIn

      this%pageName = nameIn
    end subroutine setName

    function getPageWord(this, lineNum, wordNum) result(charOut)
      class(Page), intent(in):: this
      integer, intent(in) :: lineNum, wordNum
      character(len=:), allocatable:: charOut

      if (lineNum>size(this%lines)) then
        write(*,*) "Attempting to read element ", lineNum, " from file ", this%getName(), ". This line doesn't exist"
        stop
      end if

      if (wordNum>size(this%lines(lineNum)%words)) then
        write(*,*) "Attempting to read word ", wordNum, "from line ", lineNum, " from file ", this%getName(), ". This word doesn't exist"
!        stop
      end if

      charOut =  getWithStackTrace(this, lineNum, wordNum)
      end function getPageWord

      subroutine setWithStackTraceChar(this, i, charIn)
        type(Page), intent(inout):: this
        integer, intent(in):: i
        character(len=*), dimension(:):: charIn

        type(Line):: line

        line = charIn

        this%lines(i) = line
      end subroutine setWithStackTraceChar

      subroutine setWithStackTraceArbLine(this, arbLineIn, i)
        type(Page), intent(inout):: this
        integer, intent(in):: i
        type(Line):: arbLineIn

        this%lines(i) = arbLineIn

      end subroutine setWithStackTraceArbLine

      function getWithStackTrace(this, lineNum, wordNum) result (charOut)
        class(Page), intent(in):: this
        character(len=:), allocatable:: charOut
        integer(int32), intent(in):: lineNum, wordNum
        
        charOut = this%lines(lineNum)%words(wordNum)%line
      end function getWithStackTrace

    subroutine initChar(this, charArray)
      class(Page), intent(inout):: this
      character(len=*), dimension(:,:):: charArray

      integer:: numLines, numwords
      integer::i,j

      numLines = size(charArray(1,:))
      numWords = size(charArray(:,1))

      allocate(this%lines(numLines))
      do i = 1, numLines
        allocate(this%lines(i)%words(numWords))
        do j = 1, numWords
          this%lines(i)%words(j) = charArray(i,j)
        end do
      end do
      end subroutine initChar

      subroutine initPageWithPage(this, pageIn)
        class(Page), intent(inout):: this
        type(Page), intent(in):: pageIn

        this%lines = pageIn%lines
      end subroutine initPageWithPage
    subroutine initPageWithArray(this, array)
      class(Page), intent(inout):: this
      type(Line), dimension(:):: array

      this%lines = array

      end subroutine initPageWithArray
      
!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief      Reads in a file into a page 
!
!> @details    Given an filename, this reads the lines of that file into a page type - taking away all of the comments (anything
!following a # by default, though this can be changed by passing in your own comment character).   If no file is specified, it
!attempts to read the default file name, which it gets from inputFile. I suspect that I should take this out and force you to pass
!in a file name, thus removing the dependency on inputFile.   This seems better actually.
!
!> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
!
!> @date       October 25, 2016
!
! PARAMETERS:
!> @param[optional] fileNameIn Name of the file to be used as input file
!> @return inputFile Allocatable character array where each element of the array corresponds to a single line from the input
!---------------------------------------------------------------------------
  subroutine readInInputFile(this, inputFileName, commentCharacterIn, delimiters)
    use AlphaHouseMod, only: CountLines
    class(Page), intent(inout):: this
    type(String), dimension(:), allocatable:: tempArray 
    type(Line), dimension(:), allocatable:: tempLine
    character(len=1), intent(in), optional:: commentCharacterIn
    character(len=1), dimension(:), intent(in), optional:: delimiters
    character(len=1)::commentCharacter
    character(len=*), intent(in):: inputFileName

    character(len=10000):: temp
    integer(int32):: numLines, fileId, commentPos
    integer(int32):: i

    if (present(commentCharacterIn)) then
      write(commentCharacter, "(A)")  commentCharacterIn
    else
      write(commentCharacter, *)  "#"
    end if

    numLines =  CountLines(inputFileName)

    allocate(tempArray(numLines))
    allocate(tempLine(numLines))
    open(newunit=fileId, file = inputFileName, action="read")
    do i = 1, numLines
      read(fileId, "(A)") temp
      commentPos = index(temp, commentCharacter)-1
      if (commentPos<=0) then
        commentPos = len(trim(temp))
      end if
      tempArray(i) = trim(temp(:commentPos))
      if (present(delimiters)) then
        tempLine(i) = tempArray(i)%split(delimiters)
      else
        write(*,*) tempArray(i)%getSize()
        tempLine(i) = tempArray(i)%split()
      end if
    end do
    this = tempLine
    this%pageName = inputFileName
    deallocate(tempArray)
    deallocate(tempLine)

  end subroutine readInInputFile
  end module pageModule
