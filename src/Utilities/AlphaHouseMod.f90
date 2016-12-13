
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

  use ISO_Fortran_Env

  implicit none

  private
  ! Methods
  public :: CountLines,int2Char,Real2Char,RandomOrder,ToLower,FindLoc,SetSeed,removeWhitespace,parseToFirstWhitespace,splitLineIntoTwoParts
  public :: checkFileExists, char2Int, char2Int64, char2Real, char2Double
  public:: isDelim

  !> @brief List of characters for case conversion in ToLower
  CHARACTER(*),PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*),PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  !> @brief Real2Char interface
  interface Real2Char
    module procedure RealS2Char,RealD2Char
  end interface

  !>@brief char2Int interface
  interface char2Int
    module procedure char2Int32
  end interface

  !> @brief Integer to character interface
  interface int2Char
    module procedure Int2Char32, Int2Char64
  end interface 
  !> @brief FindLoc interface
  interface FindLoc
    module procedure FindLocC, FindLocI, FindLocS, FindLocD
  end interface

  contains
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
    !> @brief   Count number of lines in a file
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function CountLines(FileName) result(nLines)
      implicit none

      ! Arguments
      character(len=*),intent(in) :: FileName !< file
      integer(int32)              :: nLines   !@result number of lines in a file

      ! Other
      integer(int32) :: f,Unit

      character(len=300) :: DumC

      nLines=0
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
    !> @brief   Convert character to integer
    !> @details See http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90. Returns a 32 bit
    !integer.
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Char2Int32(c) result(Res)
      implicit none

      character(*), intent(in) :: c   !< character
      integer(int32)           :: Res !@result integer

      read(c, *) Res
      return
    end function

    !---------------------------------------------------------------------------
    !> @brief   Convert character to integer
    !> @details See http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90. Returns a 64 bit
    !integer.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Char2Int64(c) result(Res)
      implicit none

      character(*), intent(in) :: c   !< character
      integer(int64)           :: Res !@result integer

      read(c, *) Res
      return
    end function

    !###########################################################################
    !---------------------------------------------------------------------------
    !> @brief   Convert integer to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran. Converts 64 bit integer.
    !> @author  Diarmaid de Burca, diarmaid.deburca@ed.ac.uk 
    !> @date    October 28, 2016
    !---------------------------------------------------------------------------
    function Int2Char64(i,fmt) result(Res)
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

    !---------------------------------------------------------------------------
    !> @brief   Convert integer to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function Int2Char32(i,fmt) result(Res)
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
    !> @brief   Convert real (single precision) to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function RealS2Char(r,fmt) result(Res)
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
    !> @brief   Convert character to single precision real
    !> @details Convert character to single precision real
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
    !> @details Convert character to single precision real
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
    !> @brief   Convert real (single precision) to character
    !> @details See http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function RealD2Char(r,fmt) result(Res)
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
    !> @details TODO
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
    function ToLower(StringIn) result(StringOut)
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
    function FindLocI(Val,Vec) result(i)
      implicit none
      integer(int32),intent(in) :: Val    !< value
      integer(int32),intent(in) :: Vec(:) !< vector
      integer(int32)            :: i      !< @return position

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
    function FindLocC(Val,Vec) result(i)
      implicit none
      character(len=*),intent(in) :: Val    !< value
      character(len=*),intent(in) :: Vec(:) !< vector
      integer(int32)              :: i      !< @return position

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
    !> @brief   Find position of the value in a real (single precision) vector
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function FindLocS(Val,Vec) result(i)
      implicit none
      real(real32),intent(in) :: Val    !< value
      real(real32),intent(in) :: Vec(:) !< vector
      integer(int32)          :: i      !< @return position

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        !> @todo handle floating point representation
        if (Val == Vec(j)) then
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
    function FindLocD(Val,Vec) result(i)
      implicit none
      real(real64),intent(in) :: Val    !< value
      real(real64),intent(in) :: Vec(:) !< vector
      integer(int32)          :: i      !< @return position

      integer(int32) :: j
      i=0
      do j=1,size(Vec)
        !> @todo handle floating point representation
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

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
      integer(int32),intent(in),optional  :: Seed     !< A number to initialize RNG
      character(len=*),optional           :: SeedFile !< File to save the seed in
      integer(int32),intent(out),optional :: Out      !< Make the seed value available outside

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
end module

!###############################################################################
