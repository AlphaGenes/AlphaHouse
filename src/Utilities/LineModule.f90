!###############################################################################
!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     LineModule.f90
!
! DESCRIPTION:
!> @brief    Holds a line of strings
!
!> @details  Allows you to modify strings, where strings are arbitary length characters. This is a list of those stings, like a line
!> of a page.
!
!> @author   Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
!
!-------------------------------------------------------------------------------
module LineModule
  use iso_fortran_env
  use stringModule
  implicit none

  private

  public:: assignment(=), operator(==), Line
  type :: Line
   type(String), allocatable, dimension(:):: words
    contains
      procedure, private:: addAWordWithChar
      procedure, private:: addAWordWithString
      generic:: add => addAWordWithChar, addAWordWithString
      procedure, private:: removeByIndex
      procedure, private:: removeFirstName
      generic:: remove => removeFirstName, removeByIndex
      procedure:: has => hasWithin
      procedure:: getWordAsString
      procedure:: getWord
      procedure, private:: setArbitaryLengthLine
      procedure, private:: setArbitaryLengthLineChar
      procedure:: getNumWords
      procedure, private:: setWordWithChar
      procedure, private:: setWordWithString
      generic:: setWord => setWordWithChar, setWordWithString
      procedure, private:: writeFormattedLineType
      procedure, private:: writeUnformattedLineType
      procedure, private:: readLineType
      procedure, private:: readUnformattedLineType
      procedure:: removeAll
      generic:: write(formatted) => writeFormattedLineType
      generic:: write(unformatted) => writeUnformattedLineType
      generic:: read(formatted) => readLineType
      generic:: read(unformatted) => readUnformattedLineType
      final:: deallocateLine
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
    !> @brief Checks to see if charecter is contained in Line
    !> @details Checks each string to see if it is the same as the character passed in.   If it is, then it returns the number of
    !>the string holding that character.   If it doesn't have the character it returns 0.   It has an optional logical parameter.
    !>Setting this to true will cause it to disregard case.   By default it will be false (i.e. case sensitive).
    pure integer function hasWithin(self, charIn, isCaseSensitive) result (indexOut)
      class(Line), intent(in):: self
      character(len=*), intent(in):: charIn
      logical, intent(in), optional:: isCaseSensitive

      logical:: caseSensitiveUsed
      integer:: i

      if (present(isCaseSensitive)) then
        caseSensitiveUsed = isCaseSensitive
      else
        caseSensitiveUsed = .false.
      end if

      indexOut = 0

      if (caseSensitiveUsed) then
        do i = 1, size(self%words)
          call self%words(i)%toLowerCase()
        end do
      end if
        do i =1, size(self%words)
          if (self%words(i) == charIn) then
            indexOut = i
          end if
        end do
    end function hasWithin


    !> @brief Sets a single word
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine setWordWithString(self, i, stringIn)
      class(Line), intent(inout):: self
      integer, intent(in):: i !< The index of the word you want to set
      type(String), intent(in):: stringIn !<What you want to set the word to
      self%words(i) = stringIn
    end subroutine setWordWithString
    !> @brief Sets a single word
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine setWordWithChar(self, i, charIn)
      class(Line), intent(inout):: self
      integer, intent(in):: i !< The index of the word you want to set
      character(len=*), intent(in):: charIn !<What you want to set the word to
      self%words(i) = charIn
    end subroutine setWordWithChar

    !> @brief Removes all words
    !> @details Deallocates the array of strings that are being used to hold the words
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine removeAll(self)
      class(Line), intent(inout):: self

      call deallocateLine(self)
    end subroutine removeAll

    !> @brief Final sub to deallocate Line module
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    subroutine deallocateLine(self)
      type(Line), intent(inout):: self

      if (allocated(self%words)) then
        deallocate(self%words)
      end if
    end subroutine deallocateLine
    
    

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
      if (intIn<= size(temp)) then
        temp(intIn:) = self%words(intIn+1:)
      end if
      deallocate(self%words)
      allocate(self%words(size(temp)))

      self%words(:) = temp(:)

    end subroutine removeByIndex

    subroutine addAWordWithChar(self, charIn)
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

    end subroutine addAWordWithChar

    subroutine addAWordWithString(self, stringIn)
      type(String), intent(in):: stringIn
      class(Line), intent(inout):: self

      call self%add(stringIn%line)
    end subroutine addAWordWithString

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
