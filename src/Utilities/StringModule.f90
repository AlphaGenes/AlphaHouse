Module stringModule
  use AlphaHouseMod, only: isDelim
  use iso_fortran_env
  implicit none


  private

  public:: assignment(=), operator(==)
  public:: String

  type :: String
    character(len=:), allocatable:: line
  contains
    procedure, private:: getSubStringStartAndEnd
    procedure, private:: getSubstringEnd
    procedure:: getPosition
    procedure:: toLowerCase => convertToLowerCaseString
    procedure:: toLowerFun
    procedure:: writeType
    procedure:: readType
    procedure:: getSize
    procedure:: split
    procedure:: getNumOccurances
    generic:: getSubString => getSubStringStartAndEnd, getSubstringEnd
    generic:: write(formatted) => writeType
    generic:: read(formatted) => readType
  end type 


  interface assignment (=)
    procedure setString
  end interface 

  interface operator (==)
    module procedure compareString, compareCharacter
  end interface 

contains

  pure function getNumOccurances(self, characterIn) result(numOccurances)
    class(String), intent(in):: self 
    character(len=*), intent(in):: characterIn
    integer:: numOccurances
    logical:: same


    integer:: stringCounter, characterCounter, characterMax, stringMax, tempStringCounter

    characterMax = len(characterIn)
    stringMax = len(self%line)

    numOccurances = 0
    if (characterMax>stringMax) then
      return
    end if

    if (characterIn == "") then
      return
    end if

    stringCounter = 1
    characterCounter = 1

    do while (stringCounter< stringMax)
      same = self%line(stringCounter:stringCounter) == characterIn(characterCounter:characterCounter)
      if (same) then
        tempStringCounter = stringCounter
        do while (same .and. characterCounter .lt. characterMax)
        characterCounter = characterCounter+1
        tempStringCounter = stringCounter +1
        same = self%line(tempStringCounter:tempStringCounter) == characterIn(characterCounter:characterCounter)
        end do
        if (same) then
          numOccurances = numOccurances+1
          stringCounter = stringCounter+characterMax-1
        end if
        characterCounter = 1
      end if
      stringCounter = stringCounter+1
    end do
  end function getNumOccurances

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Converts the input data into lower-case
  !
  !> @details    Given a type String, this converts it from upper case to lower case"
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[inout] input fileInput to be converted to lowercase 
  !---------------------------------------------------------------------------
  pure subroutine convertToLowerCaseString(this)
    use AlphaHouseMod
    class(String), intent(inout):: this

    this%line = ToLower(this%line)
  end subroutine convertToLowerCaseString

  pure function toLowerFun(self) result (lowerCaseOut)
    use AlphaHouseMod, only:toLower
    character(len=:), allocatable:: lowerCaseOut
    class(String), intent(in):: self

    lowerCaseOut = ToLower(self%line)
  end function toLowerFun

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Compares to see if a string is the same as another string
  !
  !> @details    Given a type String, used in comparison to see if the characters that a string holds are the same as the characters
  !that another string holds.   If they are, returns true. Otherwise, returns false.
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[in] stringIn (type: String)- string to compare to
  !> @param[out] logical
  !---------------------------------------------------------------------------
  pure function compareString(this, stringIn) result (same)
    logical:: same
    class(String), intent(in):: this
    type(String), intent(in) :: stringIn

    same = this%line .eq. stringIn%line
  end function compareString

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Compares to see if a string is the same as a character
  !
  !> @details    Given a type String, used in comparison to see if the characters that a string holds are the same as a characters.
  !If they are, returns true. Otherwise, returns false. Note that it doesn't do any trimming, it expects you to do that yourself
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[in] stringIn (type: character(len=*))- character to compare to. 
  !> @param[out] logical
  !---------------------------------------------------------------------------
  pure function compareCharacter(this, stringIn) result (same)
    logical:: same
    class(String), intent(in):: this
    character(len=*), intent(in) :: stringIn

    same = stringIn .eq. this%line
  end function compareCharacter



  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Gives number of characters in string
  !
  !> @details    Given a type String, this tells you how many characters are in the string
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[out] sizeOut. Integer( int32).  
  !---------------------------------------------------------------------------

  function getSize(this) result(sizeOut)
    class(String), intent(in):: this
    integer(int32):: sizeOut

    sizeOut = len( this%line)
  end function getSize

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief     Gives a substring on the initial string,
  !
  !> @details    given a start position and an end position, this will return the substring between the start and end positions.
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       December 1st, 2016
  !
  ! PARAMETERS:
  !---------------------------------------------------------------------------

  function getSubStringStartAndEnd(self, SubStringStart, subStringEnd) result (subString)
    class(String), intent(in):: self
    integer, intent(in):: SubstringStart, subStringEnd
    type(String):: subString

    subString = self%line(subStringStart:subStringEnd)
  end function getSubStringStartAndEnd

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief     Gives a substring on the initial string,
  !
  !> @details    given an end position, this will return the substring between the start of the string and end positions given.
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       December 1st, 2016
  !
  ! PARAMETERS:
  !---------------------------------------------------------------------------
  function getSubstringEnd(self, subStringPosition) result (subString)
    class(String), intent(in):: self
    integer, intent(in):: subStringPosition
    type(String):: subString

    subString = self%line(:subStringPosition)
  end function getSubStringEnd

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Gives position of a character in the string
  !
  !> @details    Given a character (of arbitary length) as input, this function will give you the first occurance of that character in the string.
  !Should the character not exist in the string, it will return a zero instead.
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       December 1st, 2016
  !
  !---------------------------------------------------------------------------
  function getPosition(self, charIn) result (position)
    class(String), intent(in):: self
    character(len=*):: charIn
    integer(int32):: position
    position = index(self%line, charIn)
  end function getPosition



  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Splits a string up
  !
  !> @details    Given a type String, this will return an array of Strings, where each element in the string is the substring
  !between the delimiters.   Default delimiter is a comma.   Any number of delimiters can be used, but they all need to be of size
  !one (i.e., splitting on a single character).   Also, the substrings will NOT include the delimiters.   This function will also
  !not split anything inside brackets "()".   Given an empty string it will return an array of size 1 holding an empty string.  This
  !function will also ignore consecutive delimiters (so if you have ,,, it will treat that as a single delimiter).
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[in] optional : delimitersIn array of characters of length one.   Gives the delimiters which should be used to split the
  !string up.
  !---------------------------------------------------------------------------
  function split(this, delimitersIn) result( components)
    character(len=1), dimension(:), optional :: delimitersIn
    class(String), intent(in):: this
    type(String), allocatable, dimension(:):: components
    character(len = :), allocatable:: line
    character(len=1), dimension(:), allocatable:: delimiters 

    integer(int32):: i
    integer, dimension(:,:), allocatable::numSplit

    if (present(delimitersIn)) then
      !delimiters = delimitersIn
      allocate(delimiters(size(delimitersIn)))
      do i = 1, size(delimitersIn)
         delimiters(i) = delimitersIn(i)
      end do
    else
      allocate(delimiters(1))
      delimiters(1) = ","
    end if

    line = this%line
    if (len(line)==0) then
      allocate(components(1))
      components(1) = this%line
      return
    end if


    !First find out how many times we want to split it up
    call getSplitPositions(line, numSplit, delimiters)

    allocate(components(size(numSplit(:,1))))

    if (len(line)==0) then
      components(1)=""
    else
      do i = 1, size(numSplit(:,1))
        components(i) = line(numSplit(i,1):numSplit(i,2))
      end do
    end if

  end function split

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Gets the positions by which to split the string
  !
  !> @details    Given a type String, this will return a two dimensional array of size numSplits+1, where numSplits is equal to the
  !number of times the delimiters appear in non-consecutive order.   This array holds the start and end position of the split.
  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  ! PARAMETERS:
  !> @param[in] delimiters - array of characters of length one.   Gives the character on which to split the string
  !> @param[out] SplitPositions - 2d integer array.   Gives the start and end positions of the split 
  !---------------------------------------------------------------------------

  subroutine getSplitPositions(this, SplitPositions, delimiters)
    !        type(String), intent(in)::this
    character(len=*):: this
    character(len=1), dimension(:):: delimiters
    integer(int32), dimension(:, :), allocatable:: SplitPositions

    integer(int32):: numSplits, currentSplit


    integer(int32)::characterPosition

    numSplits=1
    characterPosition=1
    if (len(this) .gt. 0) then
      do while (isDelim(this(characterPosition:characterPosition), delimiters))
        characterPosition = characterPosition+1
        if (characterPosition .ge. len(this)) then
          exit
        end if

      end do
    end if
    do while (characterPosition<len(this))
      !Everything in parenthesis is kept
      if (this(characterPosition:characterPosition) == "(") then
        do while (this(characterPosition:characterPosition) .ne. ")")
          characterPosition = characterPosition+1
          if (characterPosition>len(this)) then
            exit
          end if
        end do
        !if it is a deliminator increase the number of splits
      else if (isDelim(this(characterPosition:characterPosition), delimiters )) then
        numSplits = numSplits+1
        !increase once no matter how many delimitators follow each other
        do while (.true.)
          characterPosition = characterPosition+1
          if (characterPosition .ge. len(this) .or. .not. isDelim(this(characterPosition:characterPosition), delimiters )) then
            exit
          end if
        end do
      else
        characterPosition = characterPosition+1
      end if
    end do


    !Now allocate array to hold positions
    !reset the character position to first character
    characterPosition=1
    allocate(SplitPositions(numSplits, 2))
    SplitPositions(1,1) = 1
    !Now get the SNP positions
    if (len(this) .gt. 0) then
      do while (isDelim(this(characterPosition:characterPosition), delimiters))
        characterPosition = characterPosition+1
        SplitPositions(1,1)=characterPosition
        if (characterPosition .ge. len(this)) then
          exit
        end if

      end do
    end if
    currentSplit =1

    do while (characterPosition<len(this))
      if (this(characterPosition:characterPosition) == "(") then
        do while (this(characterPosition:characterPosition) .ne. ")")
          characterPosition = characterPosition+1
          if (characterPosition>len(this)) then
            exit
          end if
        end do
      else if (isDelim(this(characterPosition:characterPosition), delimiters)) then
        SplitPositions(currentSplit,2)=characterPosition-1

        currentSplit = currentSplit+1
        do while (.true.)
          characterPosition = characterPosition+1
          if (characterPosition .ge. len(this) .or. .not. isDelim(this(characterPosition:characterPosition), delimiters)) then
            SplitPositions(currentSplit,1) = characterPosition
            exit
          end if
        end do
      else
        characterPosition = characterPosition+1
      end if
    end do
    SplitPositions(numSplits, 2) = len(this)
    characterPosition = len(this)
    do while (isDelim(this(characterPosition:characterPosition), delimiters))
      characterPosition = characterPosition-1
      if (characterPosition ==0) then
        SplitPositions(numSplits,2) = 1
        exit
      end if
      SplitPositions(numSplits, 2) = characterPosition
    end do
  end subroutine getSplitPositions


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Sets the formatted read for a String
  !
  !> @details    Sets the formatted read for a String.   4 inputs are required by fortran.

  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  !---------------------------------------------------------------------------

  subroutine readType(dtv, unit, iotype, vlist, iostat, iomsg)
    class(String), intent(inout) :: dtv         ! Object to write.
    integer, intent(in) :: unit         ! Internal unit to write to.
    character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
    integer, intent(in) :: vlist(:)    ! parameters from fmt spec.
    integer, intent(out) :: iostat      ! non zero on error, etc.
    character(*), intent(inout) :: iomsg  ! define if iostat non zero.

    character(len=100000):: temp

    read(unit,*,iostat=iostat, iomsg=iomsg) temp
    dtv%line = trim(temp)
  end subroutine readType

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief      Sets the formatted write for a String
  !
  !> @details    Sets the formatted write for a String.   4 inputs are required by fortran.

  !
  !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
  !
  !> @date       October 25, 2016
  !
  !---------------------------------------------------------------------------
  subroutine writeType(dtv, unit, iotype, v_list, iostat, iomsg)
    class(String), intent(in) :: dtv         ! Object to write.
    integer, intent(in) :: unit         ! Internal unit to write to.
    character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
    integer, intent(in) :: v_list(:)    ! parameters from fmt spec.
    integer, intent(out) :: iostat      ! non zero on error, etc.
    character(*), intent(inout) :: iomsg  ! define if iostat non zero.

    write(unit, "(A)", iostat = iostat, iomsg = iomsg) dtv%line
  end subroutine writeType

  subroutine setString(this,  lineIn) 
    type(String), intent(inout) :: this
    character(len=*), intent(in):: lineIn

    this%line = lineIn
  end subroutine setString


end Module stringModule





