Module stringModule
  use iso_fortran_env
  implicit none


  private

  public:: assignment(=), operator(==)
  public:: String

  type :: String
    character(len=:), allocatable:: line
    contains
      procedure:: toLowerCase => convertToLowerCaseString
      procedure:: writeType
      procedure:: readType
      procedure:: getSize
      procedure:: split
!      generic:: assignment(=) => setLine
      generic:: write(formatted) => writeType
      generic:: read(formatted) => readType
  end type 


  interface assignment (=)
    procedure setLine
  end interface 
  
  interface operator (==)
    module procedure compareString, compareCharacter
  end interface 

!  interface operator (.eq.)
!    module pro
  contains

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
  subroutine convertToLowerCaseString(this)
    use AlphaHouseMod
    class(String), intent(inout):: this

    this%line = ToLower(this%line)
  end subroutine convertToLowerCaseString

    function compareString(this, stringIn) result (same)
      logical:: same
      class(String), intent(in):: this
      type(String), intent(in) :: stringIn

      same = this%line .eq. stringIn%line
      end function compareString

    function compareCharacter(this, stringIn) result (same)
      logical:: same
      class(String), intent(in):: this
      character(len=*), intent(in) :: stringIn

      same = stringIn .eq. this%line
    end function compareCharacter


    function getSize(this) result(sizeOut)
      class(String), intent(in):: this
      integer:: sizeOut

      sizeOut = len( this%line)
      end function getSize
    

          !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Splits the string into substrings based on provided delimiters
    !
    !
    !> @author     Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    !
    !> @date       October 25, 2016
    !
    ! PARAMETERS:
    !> @delimitersIn[in], optional : characters to split the string at
    !> @return Array of strings that are split at any occurance of delimiters.
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
          delimiters = delimitersIn
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

        !replace all tabs with spaces
!        do i =1 , len(line)
!          if (line(i:i) == char(9)) then
!            line(i:i) = " "
!          end if
!        end do

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
          do while (isSplitPosition(this(characterPosition:characterPosition), delimiters))
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
          else if (isSplitPosition(this(characterPosition:characterPosition), delimiters )) then
            numSplits = numSplits+1
            !increase once no matter how many delimitators follow each other
            do while (.true.)
              characterPosition = characterPosition+1
              if (characterPosition .ge. len(this) .or. .not. isSplitPosition(this(characterPosition:characterPosition), delimiters )) then
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
          do while (isSplitPosition(this(characterPosition:characterPosition), delimiters))
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
          else if (isSplitPosition(this(characterPosition:characterPosition), delimiters)) then
            SplitPositions(currentSplit,2)=characterPosition-1
            
            currentSplit = currentSplit+1
            do while (.true.)
              characterPosition = characterPosition+1
              if (characterPosition .ge. len(this) .or. .not. isSplitPosition(this(characterPosition:characterPosition), delimiters)) then
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
        do while (isSplitPosition(this(characterPosition:characterPosition), delimiters))
          characterPosition = characterPosition-1
          if (characterPosition ==0) then
            SplitPositions(numSplits,2) = 1
            exit
          end if
          SplitPositions(numSplits, 2) = characterPosition
        end do
      end subroutine getSplitPositions
      
      function isSplitPosition(charIn, delimiters)
        character(len=*):: charIn
        character(len = *), dimension(:):: delimiters
        logical:: isSplitPosition
        integer:: i

        isSplitPosition = .false.
        do i= 1, size(delimiters)
          isSplitPosition = isSplitPosition .or. charIn==delimiters(i)
        end do
      end function isSplitPosition

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

    subroutine writeType(dtv, unit, iotype, v_list, iostat, iomsg)
      class(String), intent(in) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: v_list(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      write(unit, "(A)", iostat = iostat, iomsg = iomsg, advance="no") dtv%line
      end subroutine writeType

    subroutine setLine(this,  lineIn) 
      class(String), intent(inout) :: this
      character(len=*), intent(in):: lineIn

      this%line = lineIn
    end subroutine setLine


end Module stringModule



    

