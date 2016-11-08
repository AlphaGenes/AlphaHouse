module LineModule
  use iso_fortran_env
  use stringModule
  implicit none

  private

  public:: assignment(=), operator(==), Line
  type :: Line
   type(String), allocatable, dimension(:):: words
!    character(100)::words
    contains
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
    module procedure setArbitaryLengthLine, setArbitaryLengthLineChar, setSingleChar
  end interface 

  interface operator (==)
    module procedure compareLine
  end interface 
  contains

    function getWord(this, i) result (wordOut)
      class(Line), intent(in):: this
      integer, intent(in):: i
      character(len=:), allocatable:: wordOut

      wordOut = this%words(i)%line
      end function getWord


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

      numWords = size(this%words)
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
        writeOutChar = writeOutChar // dtv%getWord(i)
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
        writeOutChar = writeOutChar // dtv%getWord(i)
    end do
        write(unit, "(A)", iostat = iostat, iomsg = iomsg) writeOutChar
        end subroutine writeUnformattedLineType

end module LineModule
