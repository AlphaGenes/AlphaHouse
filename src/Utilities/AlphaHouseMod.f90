
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaHouseMod.f90
!
! DESCRIPTION:
!> @brief    Alpha basic (house) subroutines and functions
!
!> @details  Stuff that is commonly used in Alpha suite of software and does not
!!           fit into any other module
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 GGorjanc - Initial Version
! 2016-10-28 DdeBurca - added FileCheck, changed Int2char/char2Int to support 64 bit integers
! 2016-11-12 DdeBurca - Added Char2Real and Char2Double
! 2016-11-22 DdeBurca - Added isDelim
!-------------------------------------------------------------------------------
module AlphaHouseMod

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit

  implicit none

  private
  ! Methods
  public :: Append, CountLines, int2Char, Real2Char, RandomOrder, ToLower, FindLoc, Match
  public :: removeWhitespace, parseToFirstWhitespace, splitLineIntoTwoParts
  public :: checkFileExists, char2Int, char2Real, char2Double, Log2Char
  public :: isDelim, PrintElapsedTime, intToChar, SetSeed

  public :: generatePairing, unPair
  public :: CountLinesWithBlankLines
  public:: countColumns
  !> @brief List of characters for case conversion in ToLower
  CHARACTER(*),PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*),PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface Real2Char
    module procedure Real322Char, Real642Char
  end interface

  interface int2Char
    module procedure Int82Char, Int322Char, Int642Char
  end interface

  interface intToChar
    module procedure int82CharArray, int322CharArray, int642CharArray
  end interface

  interface char2Int
    module procedure char2Int32
  end interface

  interface FindLoc
    module procedure FindLocC, FindLocI, FindLocS, FindLocD
  end interface

  interface Match
    module procedure MatchC, MatchI, MatchS, MatchD
  end interface

  interface Append
    module procedure AppendReal64, AppendChar
  end interface

  interface countColumns
    module procedure countColumnsMultiDelim, countColumnsSingleDelim
  end interface

  contains
  !> @brief A function that counts how many colums in a file
  !> Counts the number of times that a columns in a file.   Assumes that repeated delimiters are all a single delimiter, (i.e. a,,,a
  !!> would have two colums, with the delimiter being ,).   If the first character of the line is a delimiter, then that delimiter
  !> isn't counted as having a column associated with it, with the same thing being true for the final character (i.e. all of 
  !> ",A, A" and "A, A,", ",A,A," and "A, A" have two columns).  This combines with the repeated delimiters being the same (i.e.
  !> ",,,A,A,,,," has two columns).
  !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  integer function countColumnsMultiDelim(fileNameIn, delimiterIn) result (numColumnsOut)
    character(len=*), intent(in):: fileNameIn !< File to count columns 
    character(len=1), dimension(:):: delimiterIn !< List of delimiters to use

    character(len=:), allocatable:: tempChar

    integer:: fileSize, fileUnit, filePosition, fileSizeLeft
    integer:: IOStatus, i, finalLetter
    logical:: fileExists, lineEnd, previousIsDelim

    numColumnsOut = 0
    Inquire(file=fileNameIn, size=fileSize, exist=fileExists)

    if (fileSize<1000) then
      allocate(character(len=fileSize):: tempChar)
    else
      allocate(character(len=1000):: tempChar)
    end if

    ioStatus = 0
    filePosition = 1

    if (fileExists) then
      open(newunit=fileUnit, file=fileNameIn, action="read", status="old", access="stream")
      read(fileUnit, pos=1) tempChar
      finalLetter = len(tempChar)

      if (any(delimiterIn==tempChar(1:1))) then
        previousIsDelim = .true.
      else
        previousIsDelim = .false.
      end if
        
      readLoop: do while (ioStatus==0)
        filePosition = filePosition+len(tempChar)
        do i = 1, len(tempChar)
          if (any(delimiterIn == tempChar(i:i))) then
            if (.not. previousIsDelim) then
              numColumnsOut = numColumnsOut+1
              previousIsDelim =.true.
            end if
          else
            previousIsDelim = .false.
          end if
          if (tempChar(i:i) == new_line("a")) then
            finalLetter = i-1
            exit readLoop
          end if
        end do
        fileSizeLeft = fileSize-filePosition

        if (fileSizeLeft< len(tempChar)) then
          deallocate(tempChar)
          allocate(character(len=fileSizeLeft):: tempChar)
          finalLetter = fileSizeLeft
        end if
        read(fileUnit, pos=filePosition, iostat = ioStatus) tempChar
      end do readLoop

      if (.not. any(delimiterIn==tempChar(finalLetter:finalLetter))) then
        numColumnsOut = numColumnsOut+1
      end if
      close(fileUnit)
    else
      numColumnsOut = -1
    end if

  end function countColumnsMultiDelim


  !> @brief get the number of columns for a single delimiter
  !> @details A wraparound for countColumnsMultiDelim when using a single delimiter
  !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk

  integer function countColumnsSingleDelim(fileNameIn, delimiterIn) result (numColumnsOut)
    character(len=*), intent(in):: fileNameIn
    character(len=1), intent(in):: delimiterIn

    character(len=len(delimiterIn)), dimension(1):: delimiterUsed

    delimiterUsed(1) = delimiterIn

    numColumnsOut = countColumns(fileNameIn, delimiterUsed)

  end function countColumnsSingleDelim
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Checks to see if the character passed in is the same as the delimiters
  !
  !> @details     Checks to see if the character passed in is the same as the delimiters
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[in] charIn (character(len=1)) character to be checked
  !> @param[in] delimiters(character(len=1), dimension(:)) delimiters to be checked agains
  !> @param[out] logical - true if same as a delimiter, otherwise false
  !---------------------------------------------------------------------------
  function isDelim(charIn, delimiters)
    character(len=*):: charIn
    character(len = *), dimension(:):: delimiters
    logical:: isDelim
    integer:: i

    isDelim = .false.
    do i= 1, size(delimiters)
      isDelim = isDelim .or. charIn==delimiters(i)
    end do
  end function isDelim

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief      Check if a fileName exists
   !
   !> @details    given a filepath, this function checks if that file exists.   If so, it returns true, otherwise, it returns false
   !
   !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
   !
   !> @date       October 25, 2016
   !
   ! PARAMETERS:
   !> @param[in] fileName Name of the file to be checked
   !> @return .True. if file exists, otherwise .false.
   !---------------------------------------------------------------------------
   function checkFileExists(filename) result(fileExists)
     character(len=*), intent(in):: filename
     logical:: fileExists

     Inquire(file = filename, exist = fileExists)
   end function checkFileExists

    !---------------------------------------------------------------------------
    !> @brief   Count number of lines in a file including blank lines
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function CountLinesWithBlankLines(FileName) result(nLines)
      implicit none

      ! Arguments
      character(len=*),intent(in) :: FileName !< file
      integer(int32)              :: nLines   !< @return number of lines in a file

      ! Other
      integer(int32) :: f,Unit


      nLines=0
      f=0
      open(newunit=Unit,file=trim(FileName),status="old")
      do
        read(Unit,*,iostat=f)
        nLines=nLines+1
        if (f /= 0) then
          nLines=nLines-1
          exit
        end if
      end do
      close(Unit)
      return
    end function

    !---------------------------------------------------------------------------
    !> @brief   Count number of lines in a file excluding blank lines
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function CountLines(FileName) result(nLines)
      implicit none

      ! Arguments
      character(len=*),intent(in) :: FileName !< file
      integer(int32)              :: nLines   !< @return number of lines in a file

      ! Other
      integer(int32) :: f,Unit

      character(len=300) :: DumC

      nLines=0
      f=0
      open(newunit=Unit,file=trim(FileName),status="old")
      do
        read(Unit,*,iostat=f) DumC
        nLines=nLines+1
        if (f /= 0) then
          nLines=nLines-1
          exit
        end if
      end do
      close(Unit)
      return
    end function

    !###########################################################################
    !---------------------------------------------------------------------------
    !> @brief   returns string that prececdes first whitespace (used to be parsestringwindows)
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    October 18, 2016
    !---------------------------------------------------------------------------
    function parseToFirstWhitespace(str) result(newstr)
        implicit none

        integer :: k
        character(len=*), intent(in) :: str
        character(len=:), allocatable :: newstr

        character(len=512) :: dummyStr

        dummyStr = ToLower(str)
        k=1
        newstr=""
        do
            if (dummyStr(k:k) /= " ") then
                newstr = dummyStr(1:k)
                k=k+1
            else
                exit
            endif
        enddo

    end function parseToFirstWhitespace

    !###########################################################################
    !---------------------------------------------------------------------------
    !> @brief   returns string without whitespace
    !> @author  David Wilson, david.wilson@roslin.ed.ac.uk, Diarmaid de Burca, diarmaid.deBurca@ed.ac.uk
    !> @details http://stackoverflow.com/questions/27179549/removing-whitespace-in-string
    !           DDB: Modified it to also work for empty strings (just returns an empty string)
    !> @date    October 18, 2016
    !           October 28m 2016
    !---------------------------------------------------------------------------

    recursive function removewhitespace(string) result(res)
        character(len=*), intent(in) :: string
        character, parameter:: ch = ' '
        character(:), allocatable :: res

        if (len(string)==0) then
          res = ''
        else if (len(string)==1) then
           if (string==ch) then
              res = ''
           else
              res = string
           end if
        else
           if (string(1:1)==ch) then
              res = removewhitespace(string(2:))
           else
              res = string(1:1)//removewhitespace(string(2:))
           end if
        end if
    end function removewhitespace

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert character to int32
    !> @details See http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90. Returns a 32 bit
    !integer.
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Char2Int32(c) result(Res)
      implicit none

      character(*), intent(in) :: c   !< character
      integer(int32)           :: Res !< @return integer

      read(c, *) Res
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int8 to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 9, 2017
    !---------------------------------------------------------------------------
    pure function Int82Char(i,fmt) result(Res)
      implicit none

      integer(int8),intent(in)         :: i   !< integer
      character(*),intent(in),optional :: fmt !< format
      character(:),allocatable         :: Res !< @return character

      character(range(i)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) i
      else
        write(Tmp,"(i0)") i
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int32 to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function Int322Char(i,fmt) result(Res)
      implicit none

      integer(int32),intent(in)        :: i   !< integer
      character(*),intent(in),optional :: fmt !< format
      character(:),allocatable         :: Res !< @return character

      character(range(i)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) i
      else
        write(Tmp,"(i0)") i
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int64 to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran. Converts 64 bit integer.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    October 28, 2016
    !---------------------------------------------------------------------------
    pure function Int642Char(i,fmt) result(Res)
      implicit none

      integer(int64),intent(in)        :: i   !< integer
      character(*),intent(in),optional :: fmt !< format
      character(:),allocatable         :: Res !< @return character

      character(range(i)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) i
      else
        write(Tmp,"(i0)") i
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int8 to character
    !> @details Converts an integer to a character.   Character out size is given by a parameter.   Usable with arrays.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    December 13th, 2016
    !---------------------------------------------------------------------------
    function int82CharArray(i, sizeIn, fmt) result (res)
      integer(int8), intent(in), dimension(:):: i
      integer(int32), intent(in):: sizeIn
      character(len=*), intent(in), optional:: fmt
      character(len=sizeIn), dimension(size(i)):: res
      integer::j

     do j = 1, size(i)
      if (present(fmt)) then
        write(res(j), fmt) i(j)
      else
        write(res(j), "(i0)") i(j)
      end if
    end do
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int32 to character
    !> @details Converts an integer to a character.   Character out size is given by a parameter.   Usable with arrays.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    December 13th, 2016
    !---------------------------------------------------------------------------
    function int322CharArray(i, sizeIn, fmt) result (res)
      integer(int32), intent(in), dimension(:):: i
      integer(int32), intent(in):: sizeIn
      character(len=*), intent(in), optional:: fmt
      character(len=sizeIn), dimension(size(i)):: res
      integer::j

      do j = 1, size(i)
        if (present(fmt)) then
          write(res(j), fmt) i(j)
        else
          write(res(j), "(i0)") i(j)
        end if
      end do
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert int64 to character
    !> @details Converts an integer to a character.   Character out size is given by a parameter.   Usable with arrays.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    December 13th, 2016
    !---------------------------------------------------------------------------
    function int642CharArray(i, sizeIn, fmt) result (res)
      integer(int64), intent(in), dimension(:):: i
      integer(int64), intent(in):: sizeIn
      character(len=*), intent(in), optional:: fmt
      character(len=sizeIn), dimension(size(i)):: res
      integer::j

      do j = 1, size(i)
        if (present(fmt)) then
          write(res(j), fmt) i(j)
        else
          write(res(j), "(i0)") i(j)
        end if
      end do
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Convert logical to character (0/1)
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   December 23, 2016
    !---------------------------------------------------------------------------
    function Log2Char(l) result(Res)
      implicit none

      logical, intent(in) :: l   !< logical
      character(len=1)    :: Res !< @return character

      if (l) then
        Res = "1"
      else
        Res = "0"
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert character to single precision real
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    November 15, 2016
    !---------------------------------------------------------------------------

    function char2Double(charIn, format) result (realOut)
      implicit none

      real(real64):: realOut
      character(*), intent(in), optional:: format
      character(*), intent(in):: charIn

      if (present(format)) then
        read(charIn, format) realOut
      else
        read(charIn, *) realOut
      end if
    end function char2Double

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert character to single precision real
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    November 15, 2016
    !---------------------------------------------------------------------------

    function char2Real(charIn, format) result (realOut)
      implicit none

      real(real32):: realOut
      character(*), intent(in), optional:: format
      character(*), intent(in):: charIn

      if (present(format)) then
        read(charIn, format) realOut
      else
        read(charIn, *) realOut
      end if
    end function char2Real

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert real32 to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Real322Char(r,fmt) result(Res)
      implicit none

      real(real32),intent(in)          :: r   !< real
      character(*),intent(in),optional :: fmt !< format
      character(:),allocatable         :: Res !< @return character

      character(range(r)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) r
      else
        write(Tmp,"(f)") r
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Convert real64 to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Real642Char(r,fmt) result(Res)
      implicit none

      real(real64),intent(in)          :: r   !< real
      character(*),intent(in),optional :: fmt !< format
      character(:),allocatable         :: Res !< @return character

      character(range(r)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) r
      else
        write(Tmp,"(f)") r
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Generate a random ordering of the integers 1, 2, ..., n
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function RandomOrder(n) result(Order)
      implicit none

      ! Arguments
      integer(int32),intent(in)  :: n        !< number of values to shuffle
      integer(int32)             :: Order(n) !< @return randomly ordered integers

      ! Other
      integer(int32) :: i,j,k

      real(real64) :: wk(n)

      do i=1,n
        Order(i)=i
      end do

      ! Starting at the end, swap the current last indicator with one
      ! randomly chosen from those preceeding it.
      call random_number(wk)
      do i=n,2,-1
        j=1 + i * wk(i)
        if (j < i) then
          k=Order(i)
          Order(i)=Order(j)
          Order(j)=k
        end if
      end do

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Chnage case to lower
    !> @details See https://groups.google.com/forum/#!topic/comp.lang.fortran/CKx1L2Ahkxg
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function ToLower(StringIn) result(StringOut)
      implicit none

      character(len=*),intent(in) :: StringIn  !< input string
      character(len(StringIn))    :: StringOut !< @return output string in lower case
      integer :: i,n

      ! Copy input string
      StringOut=StringIn

      ! Convert case character by character
      do i=1,len(StringOut)
        n=index(UPPER_CASE,StringOut(i:i))
        if (n /= 0) StringOut(i:i)=LOWER_CASE(n:n)
      end do
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Find position of the value in an integer vector
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function FindLocI(Val,Vec) result(i)
      implicit none
      integer(int32),intent(in) :: Val    !< value
      integer(int32),intent(in) :: Vec(:) !< vector
      integer(int32)            :: i      !< @return position, 0 for no match

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Find position of the value in a character vector
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function FindLocC(Val,Vec) result(i)
      implicit none
      character(len=*),intent(in) :: Val    !< value
      character(len=*),intent(in) :: Vec(:) !< vector
      integer(int32)              :: i      !< @return position, 0 for no match

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        if (trim(Val) == trim(Vec(j))) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Find position of the value in a real (single precision) vector
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function FindLocS(Val,Vec) result(i)
      implicit none
      real(real32),intent(in) :: Val    !< value
      real(real32),intent(in) :: Vec(:) !< vector
      integer(int32)          :: i      !< @return position, 0 for no match

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        if (.not. (Val .lt. Vec(j) .or. Val .gt. Vec(j))) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Find position of the value in a real (double precision) vector
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure function FindLocD(Val,Vec) result(i)
      implicit none
      real(real64),intent(in) :: Val    !< value
      real(real64),intent(in) :: Vec(:) !< vector
      integer(int32)          :: i      !< @return position, 0 for no match

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        if (.not. (Val .lt. Vec(j) .or. Val .gt. Vec(j))) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Match one set of values onto another - character
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 4, 2017
    !-------------------------------------------------------------------------
    pure function MatchC(Set, TargetSet) result(Result)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: Set(:)              !< A set
      character(len=*), intent(in) :: TargetSet(:)        !< Target set
      integer(int32), allocatable, dimension(:) :: Result !< @return Locations of set values in the other set, 0 for no match

      ! Other
      integer(int32) :: i, n

      n = size(Set)
      allocate(Result(n))
      do i = 1, n
        Result(i) = FindLoc(Val=Set(i), Vec=TargetSet)
      end do
    end function

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Match one set of values onto another - integer
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 4, 2017
    !-------------------------------------------------------------------------
    pure function MatchI(Set, TargetSet) result(Result)
      implicit none

      ! Arguments
      integer(int32), intent(in) :: Set(:)                !< A set
      integer(int32), intent(in) :: TargetSet(:)          !< Target set
      integer(int32), allocatable, dimension(:) :: Result !< @return Locations of set values in the other set, 0 for no match

      ! Other
      integer(int32) :: i, n

      n = size(Set)
      allocate(Result(n))
      do i = 1, n
        Result(i) = FindLoc(Val=Set(i), Vec=TargetSet)
      end do
    end function

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Match one set of values onto another - single real
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 4, 2017
    !-------------------------------------------------------------------------
    pure function MatchS(Set, TargetSet) result(Result)
      implicit none

      ! Arguments
      real(real32), intent(in) :: Set(:)                  !< A set
      real(real32), intent(in) :: TargetSet(:)            !< Target set
      integer(int32), allocatable, dimension(:) :: Result !< @return Locations of set values in the other set, 0 for no match

      ! Other
      integer(int32) :: i, n

      n = size(Set)
      allocate(Result(n))
      do i = 1, n
        Result(i) = FindLoc(Val=Set(i), Vec=TargetSet)
      end do
    end function

    !###########################################################################

    !-------------------------------------------------------------------------
    !> @brief  Match one set of values onto another - double real
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 4, 2017
    !-------------------------------------------------------------------------
    pure function MatchD(Set, TargetSet) result(Result)
      implicit none

      ! Arguments
      real(real64), intent(in) :: Set(:)                  !< A set
      real(real64), intent(in) :: TargetSet(:)            !< Target set
      integer(int64), allocatable, dimension(:) :: Result !< @return Locations of set values in the other set, 0 for no match

      ! Other
      integer(int32) :: i, n

      n = size(Set)
      allocate(Result(n))
      do i = 1, n
        Result(i) = FindLoc(Val=Set(i), Vec=TargetSet)
      end do
    end function

    !###########################################################################
    !---------------------------------------------------------------------------
    !> @brief   splits string into initial and second part, where second part is another array
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    October 18, 2016
    !---------------------------------------------------------------------------
    subroutine splitLineIntoTwoParts(line,first, second)

     implicit none

    integer :: lenin,i
    integer :: sCount1, sCount2, fCount ! SCount1 starts from 0 to avoid extra allocation
    integer, parameter :: numberAfterComma = 10000 ! Can be max 10000 params after first comma
    character(len=*), intent(in) :: line
    character :: c
    character(len=300) :: first
    character (len=300),dimension(:), allocatable,intent(out) :: second
    character (len=300), dimension(:), allocatable :: tmp
    logical :: useSecond
    first = " "
    useSecond = .false.
    sCount1 = 0
    sCount2 = 1
    fCount = 1
    if (allocated(second)) deallocate(second)
    allocate(second(numberAfterComma)) ! Allocate the second array
    lenin=len(line)
    lenin=len_trim(line(1:lenin))  ! Trim string as no point dealing with trailing whitespace
    do i=1,lenin
        c=line(i:i)
        if(.not. (ichar(c) == 9 .or. c == " " .or. c=="")) then  !if c is not tab or whitespace
            if (c == ',') then
                sCount1 = sCount1 + 1
                scount2 = 1 ! reset scount 2
                useSecond = .true.
                second(sCount1) = " "
            else if (useSecond) then
                second(sCount1)(sCount2:sCount2) = c
                scount2 = sCount2 + 1
            else
                first(fCount:fCount) = c
                fCount = fCount + 1
            endif
        else
            !We know that it is either a tab or whitespace, so why add them
            continue
        endif
    enddo
    if (sCount1 > 0) then
        allocate(tmp(sCount1)) !Allocate temp to count after
        tmp(1:sCount1) = second(1:sCount1) ! make tmp contain values of second
        call move_alloc(tmp, second)
    else
        deallocate(second) ! Deallocate second if nothing has come afterwards
    endif
    end subroutine splitLineIntoTwoParts

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Set seed
    !> @details Standard Fortran seed approach
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Set seed, potentially created file (SeedFile), and potentially
    !!          returned seed value (Out)
    !---------------------------------------------------------------------------
    subroutine SetSeed(Seed,SeedFile,Out)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional    :: Seed     !< A number to initialize RNG
      character(len=*),intent(in), optional :: SeedFile !< File to save the seed in
      integer(int32),intent(out),optional   :: Out      !< Make the seed value available outside

      ! Other
      integer(int32) :: Size,Unit
      integer(int32),allocatable :: SeedList(:)

      ! Get the size of seed array
      call random_seed(size=Size)
      allocate(SeedList(Size))

      ! Set seed
      if (present(Seed)) then ! using the given value
        SeedList(1)=Seed
        SeedList(2:Size)=1
        call random_seed(put=SeedList)
      else                    ! using system/compiler value
        call random_seed
        call random_seed(get=SeedList)
      end if

      ! Save to a file
      if (present(SeedFile)) then
        open(newunit=Unit,file=trim(SeedFile),status="unknown")
        write(Unit,*) SeedList(1)
        close(Unit)
      end if

      ! Output
      if (present(Out)) then
          Out=SeedList(1)
      end if
      deallocate(SeedList)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print elapsed time in a nice way
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    December 22, 2016
    !> @return  Print on standard output
    !---------------------------------------------------------------------------
    subroutine PrintElapsedTime(Start, End)
      implicit none
      real(real32) :: Start !< Start time from cpu_time()
      real(real32) :: End   !< End time from cpu_time()

      real(real32) :: Total
      integer(int32) :: Hours, Minutes, Seconds

      Total = 0.0
      Hours = 0
      Minutes = 0
      Seconds = 0

      Total = End - Start
      Minutes = int(Total / 60)
      Seconds = int(Total - (Minutes * 60))
      Hours = int(Minutes / 60)
      Minutes = Minutes - (Hours * 60)

      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") " Process ended"
      write(STDOUT, "(a,f20.2,a,3(i4,a))") " Time elapsed: ", Total, " seconds => ",&
                                           Hours,   " hours ",&
                                           Minutes, " minutes ",&
                                           Seconds, " seconds"
      write(STDOUT, "(a)") ""
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   szudzik pairing function in fortran
    !> @details Generates a unique pairing based on two integers
    !> If input is (N,M) space, output will be (N*M) space.
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !---------------------------------------------------------------------------
    elemental function generatePairing(xin,yin) result(res)

      integer(int32), intent(in) :: xin, yin
      integer(int32) :: x, y
      integer(int64) :: res

      ! ensures that order (e.g. [1,2] and [2,1]) doesn't matter
       if (xin < yin) then
        x = yin
        y = xin
      else
        x = xin
        y = yin
      endif

      if (x >= y) then
        res = x * x + x + y

      else
        res =y * y + x
      endif

    end function generatePairing

    !---------------------------------------------------------------------------
    !> @brief   szudzik unpairing function in fortran
    !> @details returns two integers that generated unique number based on pair
    !> If using tuples, will overflow very quickly.
    !< @author  David Wilson david.wilson@roslin.ed.ac.uk
    !---------------------------------------------------------------------------
    subroutine unPair(num, xout, yout)
      integer(int64), intent(in) :: num !< number to unPair
      integer(int32), intent(out) :: xout , yout !< numbers used to get pairing function

      integer :: x,y
      real(kind=real32) :: sqrtz, sqrz
      sqrtz = floor(SQRT(real(num)))
      sqrz = sqrtz * sqrtz

      if ((num-sqrz) >= sqrtz) then
        x = sqrtz
        y = num - sqrz - sqrtz
      else
        x = num - sqrz
        y = sqrtz
      endif

      ! ensures that order (e.g. [1,2] and [2,1]) doesn't matter
      if (y > x) then
        xout = y
        yout = x
      else
        xout = x
        yout = y
      endif
    end subroutine unPair

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Append value y at the end of a vector x - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   May 1, 2017
    !---------------------------------------------------------------------------
    pure subroutine AppendReal64(x, y)
      implicit none
      real(real64), intent(inout), allocatable :: x(:) !< @return Appended vector
      real(real64), intent(in)                 :: y    !< Value
      integer(int32) :: n
      real(real64), allocatable :: Tmp(:)
      if (allocated(x)) then
        n = size(x)
        allocate(Tmp(n))
        Tmp = x
        deallocate(x)
        allocate(x(n + 1))
        x(1:n) = Tmp
        n = n + 1
        x(n) = y
      else
        n = 1
        allocate(x(n))
        x(n) = y
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Append value y at the end of a vector x - char
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   May 1, 2017
    !> @todo   Can we get rid of Len?
    !---------------------------------------------------------------------------
    pure subroutine AppendChar(x, y, len)
      implicit none
      character(len=len), intent(inout), allocatable :: x(:) !< @return Appended vector
      character(len=*),   intent(in)                 :: y    !< Value
      integer(int32), intent(in)                     :: Len  !< Length of character value
      integer(int32) :: n
      character(len=len), allocatable :: Tmp(:)
      if (allocated(x)) then
        n = size(x)
        allocate(Tmp(n))
        Tmp = x
        deallocate(x)
        allocate(x(n + 1))
        x(1:n) = Tmp
        n = n + 1
        x(n) = y
      else
        n = 1
        allocate(x(n))
        x(n) = y
      end if
    end subroutine

    !###########################################################################

end module

!###############################################################################
