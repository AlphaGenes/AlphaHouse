module LineModule
  use iso_fortran_env
  use stringModule
  implicit none

  private

  public:: assignment(=), operator(==), Line
  type :: Line
   type(String), allocatable, dimension(:):: words
    contains
      procedure:: add => addAWord
      procedure, private:: removeByIndex
      procedure, private:: removeFirstName
      generic:: remove => removeFirstName, removeByIndex
      procedure:: getWordAsString
      procedure:: getWord
      procedure, private:: setArbitaryLengthLine
      procedure, private:: setArbitaryLengthLineChar
      procedure:: getNumWords
      procedure, private:: writeFormattedLineType
      procedure, private:: writeUnformattedLineType
      procedure, private:: readLineType
      procedure, private:: readUnformattedLineType
      generic:: write(formatted) => writeFormattedLineType
      generic:: write(unformatted) => writeUnformattedLineType
      generic:: read(formatted) => readLineType
      generic:: read(unformatted) => readUnformattedLineType
  end type

  interface assignment (=)
    module procedure setArbitaryLengthLine
    module procedure setArbitaryLengthLineChar
    module procedure setSingleChar
  end interface 

  interface operator (==)
    module procedure compareLine
  end interface 
  contains

    !>@brief Removes the first occurance of a word in a line
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine removeFirstName(self, charIn)
      class(Line), intent(inout):: self
      character(len=*)::charIn !< word that you want to remove from the Line

      integer:: i

      do i = 1, self%getNumWords()
        if (self%words(i) == charIn) then
          call self%remove(i)
          return
        end if
      end do
    end subroutine removeFirstName

    !>@brief Removes a word based on the index position
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine removeByIndex(self, intIn)
      class(Line), intent(inout):: self 
      integer(int32), intent(in):: intIn !< position you want to remove

      type(String), allocatable, dimension(:):: temp

      allocate(temp(self%getNumWords()-1))
      temp(1:intIn-1) = self%words(1:intIn-1)
      write(*,*) size(temp), intIn
      if (intIn< size(temp)) then
        temp(intIn:) = self%words(intIn+1:)
      end if
      deallocate(self%words)
      allocate(self%words(size(temp)))

      self%words(:) = temp(:)

    end subroutine removeByIndex

    subroutine addAWord(self, charIn)
      class(Line), intent(inout):: self
      character(len=*), intent(in):: charIn

      integer, dimension(1):: newSize
      type(String), dimension(1):: newString
      type(String), dimension(:), allocatable:: testString

      if (allocated(self%words)) then
        newString(1) = charIn
        newSize(1) = self%getNumWords()+1
        allocate(testString(newSize(1)))
        testString(1:self%getNumWords()) = self%words
        testString(newSize(1)) = newString(1)
!        testString =  reshape(self%words, newSize, newString)
        deallocate(self%words)
        allocate(self%words(newSize(1)))
        self%words(:) = testString(:)
      else
        allocate(self%words(1))
        self%words(1) = charIn
      end if

    end subroutine addAWord

    function getWord(this, i) result (wordOut)
      class(Line), intent(in):: this
      integer, intent(in):: i
      character(len=:), allocatable:: wordOut

      wordOut = this%words(i)%line
    end function getWord

    function getWordAsString(this, i) result (stringOut)
      class(Line), intent(in):: this
      integer, intent(in):: i
      type(String):: stringOut

      stringOut = this%getWord(i)
    end function getWordAsString


    function compareLine(this, lineIn) result (same)
      logical:: same
      logical:: sameStrings
      class(Line), intent(in):: this
      type(Line), intent(in):: lineIn
      integer(int32):: i

      sameStrings = .true.
      if (this%getNumWords() == lineIn%getNumWords()) then
        do i = 1, this%getNumWords()
          sameStrings = sameStrings .and. this%getWord(i)==lineIn%getWord(i)
        end do
        same = sameStrings
      else
        same = .false.
      end if
      end function compareLine

      

    function getNumWords(this) result(numWords)
      class(Line), intent(in):: this
      integer(int32)::numWords

      if (allocated(this%words)) then
        numWords = size(this%words)
      else
        numWords = 0
      end if
    end function 

    subroutine setSingleChar(this, lineIn)
      class(Line), intent(inout):: this
      character(len=*) :: lineIn

      if (allocated(this%words)) then
        deallocate(this%words)
      end if
      allocate(this%words(1))
      this%words(1) = lineIn
      end subroutine setSingleChar

    subroutine setArbitaryLengthLine(this, lineIn)
      class(Line), intent(inout)::this
      type(String), intent(in), dimension(:) ::lineIn

      
       this%words = lineIn
    end subroutine setArbitaryLengthLine

    subroutine setArbitaryLengthLineChar(this, lineIn)
    class(Line), intent(inout):: this
    character(len=*), dimension(:):: lineIn

    integer:: i
    
    if (allocated(this%words)) then
      deallocate(this%words)
    end if
    allocate(this%words(size(lineIn)))
    do i = 1, size(lineIn)
    this%words(i) = lineIn(i)
    end do
    end subroutine setArbitaryLengthLineChar


    subroutine readUnformattedLineType(dtv, unit, iostat, iomsg)
      class(Line), intent(inout) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.
    
        read(unit,*,iostat=iostat, iomsg=iomsg) dtv%words
    end subroutine readUnformattedLineType

    subroutine readLineType(dtv, unit, iotype, vlist, iostat, iomsg)
      class(Line), intent(inout) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: vlist(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.
    
      read(unit,iotype,iostat=iostat, iomsg=iomsg) dtv%words
    end subroutine readLineType

    subroutine writeFormattedLineType(dtv, unit, iotype, v_list, iostat, iomsg)
      class(Line), intent(in) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: v_list(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      character(len=:), allocatable:: writeOutChar
      integer:: i
      
      writeOutChar = ""
      do i = 1, dtv%getNumWords()
        writeOutChar = writeOutChar // " " // dtv%getWord(i)
    end do
      write(unit, "(A)", iostat = iostat, iomsg = iomsg) writeOutChar
    
      end subroutine writeFormattedLineType

      subroutine writeUnformattedLineType(dtv, unit, iostat, iomsg)
        class(Line), intent(in)::dtv
        integer, intent(in):: unit
        integer, intent(out) :: iostat      ! non zero on error, etc.
        character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      character(len=:), allocatable:: writeOutChar
      integer:: i
      
      writeOutChar = ""
      do i = 1, dtv%getNumWords()
        writeOutChar = writeOutChar // " " //dtv%getWord(i)
    end do
        write(unit, "(A)", iostat = iostat, iomsg = iomsg) writeOutChar
        end subroutine writeUnformattedLineType

end module LineModule
