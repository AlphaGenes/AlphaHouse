
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     HaplotypeModule.f90
!
! DESCRIPTION:
!> @brief    Module cotaining a type representing a Haplotype and associated
!>	     procedures
!
!> @details  currently only contains integer and real heap sort procedures 
!
!> @author   Daniel Money, daniel.money@roslin.ed.ac.uk
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 dmoney - initial version
! 2017-05-25 dwilso18 - revision with comments
!
!-------------------------------------------------------------------------------



module HaplotypeModule
  
    use constantModule, only: MissingPhaseCode,ErrorPhaseCode
    implicit none
    
    
    type :: Haplotype
      !Reminder as to how the data is stored...
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Phase       Phase   Missing !
      ! 0           0       0       !
      ! 1           1       0       !
      ! Missing     0       1       !
      ! Error       1       1       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(kind=8), dimension(:), allocatable :: phase
      integer(kind=8), dimension(:), allocatable :: missing
      integer :: sections
      integer :: overhang
      integer :: length
    contains
    procedure :: toIntegerArray => haplotypeToIntegerArray
    procedure :: toIntegerArrayWithErrors
    procedure :: getPhase
    procedure :: setPhase
    procedure :: overlap
    procedure :: mismatches
    procedure :: noMismatches
    procedure :: compatible
    procedure :: merge
    procedure :: numberMissing
    procedure :: numberMissingOrError
    procedure :: numberNotMissing
    procedure :: numberNotMissingOrError
    procedure :: numberError
    procedure :: compareHaplotype
    procedure :: numberSame
    procedure :: fullyPhased
    procedure :: setUnphased
    procedure :: getLength
    procedure :: isMissing
    procedure :: numberBothNotMissing
    procedure :: setFromOther
    procedure :: setFromOtherIfMissing
    procedure :: setZeroBits
    procedure :: setOneBits
    procedure :: setZero
    procedure :: setOne
    procedure :: subset
    procedure :: setSubset
    procedure :: equalHap
    procedure :: isSubset
    procedure :: setErrorToMissing
    procedure :: readUnformattedHaplotype
    procedure :: readFormattedHaplotype
    procedure :: writeFormattedHaplotype
    procedure :: writeunFormattedHaplotype

    generic:: write(formatted)=> writeFormattedHaplotype
    generic:: write(unformatted)=> writeunFormattedHaplotype
    generic:: read(formatted) => readFormattedHaplotype
    generic:: read(unformatted) => readunFormattedHaplotype
  end type Haplotype
  
  interface Haplotype
    module procedure newHaplotypeInt
    module procedure newHaplotypeBits
    module procedure newHaplotypeMissing
    module procedure newHaplotypeHaplotype
  end interface Haplotype
  
contains

    !---------------------------------------------------------------------------
    !> @brief	Constructs a new Haplotype from a integer array
    !> @date    May 25, 2017
    !> @return	New haplotype object 
    !---------------------------------------------------------------------------
  pure function newHaplotypeInt(hap) result(h)
    integer(kind=1), dimension(:), intent(in) :: hap
    
    type(Haplotype) :: h
    
    integer :: i, cursection, curpos
    
    h%length = size(hap,1)
    
    h%sections = h%length / 64 + 1
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    cursection = 1
    curpos = 0
    h%phase = 0
    h%missing = 0
    do i = 1, h%length
        select case (hap(i))
        case (0)
        ! Nothing to do due to defaults
        case (1)
            h%phase(cursection) = ibset(h%phase(cursection), curpos)
        case (MissingPhaseCode)
            h%missing(cursection) = ibset(h%missing(cursection), curpos)    
        case default
            h%phase(cursection) = ibset(h%phase(cursection), curpos)
            h%missing(cursection) = ibset(h%missing(cursection), curpos)
        end select
        curpos = curpos + 1
        if (curpos == 64) then
            curpos = 0
            cursection = cursection + 1
        end if
    end do
  end function newHaplotypeInt

    !---------------------------------------------------------------------------
    !> @brief	Constructs a new Haplotype from bit arrays
    !> @date    May 25, 2017
    !> @return	New haplotype object 
    !---------------------------------------------------------------------------
  function newHaplotypeBits(phase, missing, length) result(h)
    integer(kind=8), dimension(:), allocatable, intent(in) :: phase, missing
    integer :: length
    
    type(Haplotype) :: h
    
    integer :: i
    
    h%length = length
    h%sections = size(phase,1)
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    h%phase = phase
    h%missing = missing
    
    do i = 64 - h%overhang, 63
        h%phase(h%sections) = ibclr(h%phase(h%sections), i)
        h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
    
  end function newHaplotypeBits

    !---------------------------------------------------------------------------
    !> @brief	Constructs a new Haplotype as  a copy of another Haplotype
    !> @date    May 25, 2017
    !> @return	New haplotype object 
    !---------------------------------------------------------------------------  
  function newHaplotypeHaplotype(oh) result(h)
    class(Haplotype) :: oh
    
    type(Haplotype) :: h
    
    h = Haplotype(oh%phase, oh%missing, oh%length)
  end function newHaplotypeHaplotype
  
  
    !---------------------------------------------------------------------------
    !> @brief	Constructs a new Haplotype with all positions set to missing
    !> @date    May 25, 2017
    !> @return	New haplotype object 
    !---------------------------------------------------------------------------      
  function newHaplotypeMissing(length) result(h)
    integer, intent(in) :: length
    
    type(Haplotype) :: h
    
    integer :: i
    
    h%length = length
    h%sections = h%length / 64 + 1
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    h%phase = 0
    h%missing = 0
    h%missing = NOT(h%missing)

    
    do i = 64 - h%overhang, 63
        h%phase(h%sections) = ibclr(h%phase(h%sections), i)
        h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
  end function newHaplotypeMissing
  
  
    !---------------------------------------------------------------------------
    !> @brief	Converts a Haplotype to an integer array
    !> @detail	Error states are coded as missing (i.e. 9)
    !> @date    May 25, 2017
    !> @return	Integer array representing the haplotype
    !---------------------------------------------------------------------------  
  function haplotypeToIntegerArray(h) result(array)
    class(Haplotype), intent(in) :: h
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(h%length))
    
    cursection = 1
    curpos = 0
    do i = 1, h%length
        if (btest(h%missing(cursection),curpos)) then
            if (btest(h%phase(cursection),curpos)) then
                array(i) = MissingPhaseCode
            else
                array(i) = MissingPhaseCode
            end if
        else
            if (btest(h%phase(cursection),curpos)) then
                array(i) = 1
            else
                array(i) = 0
            end if
        end if
      
      curpos = curpos + 1
      if (curpos == 64) then
        curpos = 0
        cursection = cursection + 1
      end if
    end do
  end function haplotypeToIntegerArray  
  
    !---------------------------------------------------------------------------
    !> @brief	Converts a Haplotype to an integer array
    !> @detail	Error states are coded diffeerent to missing (i.e. as -1)
    !> @date    May 25, 2017
    !> @return	Integer array representing the haplotype
    !---------------------------------------------------------------------------   
  function toIntegerArrayWithErrors(h) result(array)
    class(Haplotype), intent(in) :: h
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(h%length))
    
    cursection = 1
    curpos = 0
    do i = 1, h%length
        if (btest(h%missing(cursection),curpos)) then
            if (btest(h%phase(cursection),curpos)) then
                array(i) = ErrorPhaseCode
            else
                array(i) = MissingPhaseCode
            end if
        else
            if (btest(h%phase(cursection),curpos)) then
                array(i) = 1
            else
                array(i) = 0
            end if
        end if
      
      curpos = curpos + 1
      if (curpos == 64) then
        curpos = 0
        cursection = cursection + 1
      end if
    end do
  end function toIntegerArrayWithErrors

    !---------------------------------------------------------------------------
    !> @brief	Compares two haplotypes
    !> @date    May 25, 2017
    !> @return	Whether the two haplotypes are identical
    !---------------------------------------------------------------------------   
  function compareHaplotype(h1, h2) result(same)
    class(Haplotype), intent(in) :: h1, h2
    
    logical :: same
    
    integer :: i
    
    if (h1%length /= h2%length) then
        same = .false.
    else
        same = .true.
        do i = 1, h1%sections
            same = same .and. (h1%phase(i) == h2%phase(i)) .and. (h1%missing(i) == h2%missing(i))
        end do
    end if
  end function compareHaplotype
  
  
    !---------------------------------------------------------------------------
    !> @brief	Gets the phase of a snp at the given position
    !> @date    May 25, 2017
    !> @return	The phase
    !---------------------------------------------------------------------------   
  function getPhase(h, pos) result (phase)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: pos
    
    integer(kind=1) :: phase
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1
  
    if (btest(h%missing(cursection),curpos)) then
        if (btest(h%phase(cursection),curpos)) then
            phase = ErrorPhaseCode
        else
            phase = MissingPhaseCode
        end if
    else
        if (btest(h%phase(cursection),curpos)) then
            phase = 1
        else
            phase = 0
        end if
    end if
  end function getPhase
  
    !---------------------------------------------------------------------------
    !> @brief	Sets the phase of a snp at the given position
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------   
  subroutine setPhase(h, pos, phase)
    class(Haplotype) :: h
    integer, intent(in) :: pos
    integer(kind=1) :: phase

    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1 
    
    select case (phase)
      case (0)
    h%phase(cursection) = ibclr(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
      case (1)
    h%phase(cursection) = ibset(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
      case (MissingPhaseCode)
    h%phase(cursection) = ibclr(h%phase(cursection), curpos)
    h%missing(cursection) = ibset(h%missing(cursection), curpos)
      case default
    h%phase(cursection) = ibset(h%phase(cursection), curpos)
    h%missing(cursection) = ibset(h%missing(cursection), curpos)
    end select
  end subroutine setPhase
  
    !---------------------------------------------------------------------------
    !> @brief	Counts the overlap between two haplotypes
    !> @detail  Counts the number of snps where both haplotypes are present
    !>		  (i.e. not missing and not error)
    !> @date    May 25, 2017
    !> @return	The overlap between the two haplotypes
    !---------------------------------------------------------------------------   
  function overlap(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND(NOT(h1%missing(i)), NOT(h2%missing(i))))
    end do
    
    num = num - h1%overhang
  end function overlap

    !---------------------------------------------------------------------------
    !> @brief	Counts the mismatches (i.e. opposing homozygotes) between two haplotypes
    !> @date    May 25, 2017
    !> @return	The number of mismatches between the two haplotypes
    !---------------------------------------------------------------------------  
  function mismatches(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
    
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND( IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
        IXOR(h1%phase(i), h2%phase(i)) ))
    end do
  end function mismatches
  
    !---------------------------------------------------------------------------
    !> @brief	Tests whether there's any mismatches (i.e. opposing homozygotes) between two haplotypes
    !> @date    May 25, 2017
    !> @return	Whether there are any mismathes (true if there are none)
    !---------------------------------------------------------------------------  
  function noMismatches(h1, h2) result (noMiss)
    class(Haplotype), intent(in) :: h1, h2
    
    logical :: noMiss
    
    integer :: i
    
    noMiss = .true.
    do i = 1, h1%sections
      if ((IAND( IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
        IXOR(h1%phase(i), h2%phase(i)) )) /= 0) then
	noMiss = .false.
	exit
      end if
    end do
  end function noMismatches
  
    !---------------------------------------------------------------------------
    !> @brief	Returns whether two haplotypes are compatible
    !> @detail  Returns true if the number of mismatches is less than the given
    !>		  threshold and the overlap is greater than or equal to the given
    !>		  threshold.
    !> @date    May 25, 2017
    !> @return	The overlap between the two haplotypes
    !---------------------------------------------------------------------------  
  function compatible(h1, h2, allowedMismatches, minOverlap) result(c)
    class(Haplotype), intent(in) :: h1, h2
    integer, intent(in) :: allowedMismatches, minOverlap
    
    logical :: c
    
    c = ((mismatches(h1,h2) <= allowedMismatches) .and. (overlap(h1,h2) >= minOverlap))
  end function compatible
  
    !---------------------------------------------------------------------------
    !> @brief	Merges two haplotypes
    !> @detail  Merges two haplotypes.  Uses the following rules per snp:
    !>		  1) Either haplotype is error then error
    !>		  2) Both haplotypes are missing then missing
    !>		  3) One haplotype is missing the other phased then that phase
    !>		  4) Both haplotypes phased the same then that phase
    !>		  5) Both haplotypes phased but opposing then error
    !> @date    May 25, 2017
    !> @return	The merged haplotype
    !---------------------------------------------------------------------------  
  function merge(h1,h2) result(h)
    class(Haplotype), intent(in) :: h1, h2
    
    type(Haplotype) :: h
    
    integer :: i
    
    h%sections = h1%sections
    h%overhang = h1%overhang
    h%length = h1%length
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    
    do i = 1, h1%sections
      h%missing(i) = IOR( IOR( &
	! Either is error
	IOR( IAND(h1%missing(i), h1%phase(i)), IAND(h2%missing(i), h2%phase(i))), &
	! Both are missing (or error)
	IAND(h1%missing(i), h2%missing(i))), &
	! Both are present both opposed
	IAND(IAND(NOT(h1%missing(i)),NOT(h2%missing(i))), IXOR(h1%phase(i), h2%phase(i))))
      h%phase(i) = IOR(h1%phase(i), h2%phase(i))
    end do
  end function merge
  
    !---------------------------------------------------------------------------
    !> @brief	Tests whether the haplotypes is fully phased
    !> @date    May 25, 2017
    !> @return	Whether the haplotype is fully phased
    !---------------------------------------------------------------------------    
  function fullyPhased(h) result(fully)
    class(Haplotype) :: h
    
    logical :: fully
    
    integer :: i
    
    fully = .true.
    
    do i = 1, h%sections
        fully = fully .and. (h%missing(i) == 0)
    end do
  end function fullyPhased
  
    !---------------------------------------------------------------------------
    !> @brief	Returns the number of missing snps in the haplotype
    !> @detail	Error snps are NOT missing for this function
    !> @date    May 25, 2017
    !> @return	The number of missing snps
    !--------------------------------------------------------------------------- 
  function numberMissing(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
        num = num + POPCNT(IAND(h%missing(i), NOT(h%phase(i))))
    end do
    
  end function numberMissing
  
    !---------------------------------------------------------------------------
    !> @brief	Returns the number of missing or error snps in the haplotype
    !> @date    May 25, 2017
    !> @return	The number of missing or error snps
    !--------------------------------------------------------------------------- 
  function numberMissingOrError(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
        num = num + POPCNT(h%missing(i))
    end do
    
  end function numberMissingOrError

    !---------------------------------------------------------------------------
    !> @brief	Returns the number of error snps in the haplotype
    !> @date    May 25, 2017
    !> @return	The number of error snps
    !---------------------------------------------------------------------------   
  function numberError(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
        num = num + POPCNT(IAND(h%missing(i), h%phase(i)))
    end do
    
  end function numberError
  
    !---------------------------------------------------------------------------
    !> @brief	Returns the number of not missing snps
    !> @detail	Error snps are NOT missing for this function
    !> @date    May 25, 2017
    !> @return	The number of not missing snps
    !---------------------------------------------------------------------------     
  function numberNotMissing(h) result(num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    num = h%length - h%numberMissing()
  end function numberNotMissing

    !---------------------------------------------------------------------------
    !> @brief	Returns the number of snps not missing or error (i.e. phased)
    !> @date    May 25, 2017
    !> @return	The number of phased snps
    !---------------------------------------------------------------------------      
  function numberNotMissingOrError(h) result(num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    num = h%length - h%numberMissingOrError()
  end function numberNotMissingOrError
  
    !---------------------------------------------------------------------------
    !> @brief	Returns the number of snps that are the same in both haplotypes
    !> @date    May 25, 2017
    !> @return	The number of similar snps
    !---------------------------------------------------------------------------    
  function numberSame(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
    
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT( IAND( &
          IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
          NOT(IEOR(h1%phase(i), h2%phase(i))) ))
    end do
    num = num - h1%overhang
  end function numberSame
  
    !---------------------------------------------------------------------------
    !> @brief	Sets the whole haplotype to be unphased
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------
  subroutine setUnphased(h)
    class(Haplotype) :: h
    
    integer :: i
    
    h%phase = 0
    h%missing = 0
    h%missing = NOT(h%missing)
    
    do i = 64 - h%overhang, 63
        h%phase(h%sections) = ibclr(h%phase(h%sections), i)
        h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
  end subroutine setUnphased
  
    !---------------------------------------------------------------------------
    !> @brief	Gets the length of the haplotype
    !> @date    May 25, 2017
    !> @return	The length of the haplotype
    !---------------------------------------------------------------------------  
  function getLength(h) result(l)
    class(Haplotype), intent(in) :: h
    integer :: l
    
    l = h%length
  end function getLength

    !---------------------------------------------------------------------------
    !> @brief	Tests whether the given snp is missing
    !> @date    May 25, 2017
    !> @return	Whether the given snp is missing
    !---------------------------------------------------------------------------    
  function isMissing(h, pos) result (missing)
   class(Haplotype), intent(in) :: h
   integer, intent(in) :: pos

   logical :: missing

   integer :: cursection, curpos

   cursection = (pos-1) / 64 + 1
   curpos = pos - (cursection - 1) * 64 - 1


   missing = BTEST(h%missing(cursection), curpos)
  end function isMissing
  
    !---------------------------------------------------------------------------
    !> @brief	Counts the number of snps that are present in both haplotypes
    !> @detail	Error snps are counted as missing for this function. (This
    !>		  may need to change
    !> @date    May 25, 2017
    !> @return	The number of snps present in both haplotypes
    !---------------------------------------------------------------------------    
  function numberBothNotMissing(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
        num = num + POPCNT(NOT(IOR(h1%missing(i), h2%missing(i))))
    end do
    
    num = num - h1%overhang
  end function numberBothNotMissing

    !---------------------------------------------------------------------------
    !> @brief	Sets one haplotypes from another
    !> @detail  Uses the following rules per snp:
    !>		  1) If haplotype 2 is phased then return that phase
    !>		  2) If haplotype 1 is phased but haplotype 2 is missing then
    !>			return that phase
    !>		  3) If both are unphased and either error then error
    !>		  4) If both are missing then missing
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setFromOther(h, oh)
    class(Haplotype) :: h
    type(Haplotype) :: oh
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(h%missing(i), oh%missing(i))
      h%phase(i) = IOR(IAND(oh%missing(i),h%phase(i)),IAND(NOT(oh%missing(i)),oh%phase(i)))
    end do
  end subroutine setFromOther
  
    !---------------------------------------------------------------------------
    !> @brief	Sets one haplotypes from another if the first is missing
    !> @detail  Uses the following rules per snp:
    !>		  1) If haplotype 1 is phased then return that phase
    !>		  2) If haplotype 1 is missing but haplotype 2 is phased then
    !>			return that phase
    !>		  3) If both are unphased and both error then error
    !>		  TODO point 3 oddness
    !>		  4) If both are missing then missing
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setFromOtherIfMissing(h, oh)
    class(Haplotype) :: h
    type(Haplotype) :: oh
    
    integer :: i
    
    do i = 1, h%sections
      h%phase(i) = IOR(IAND(NOT(h%missing(i)),h%phase(i)), IAND(h%missing(i), oh%phase(i)))
      h%missing(i) = IAND(h%missing(i), oh%missing(i))      
    end do
  end subroutine setFromOtherIfMissing

    !---------------------------------------------------------------------------
    !> @brief	Sets the snps in a haplotype to one
    !> @detail  Sets a snp to be one if the same position in the bit array is set
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setOneBits(h, array)
    class(Haplotype) :: h
    integer(kind=8), dimension(:), intent(in) :: array
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(NOT(array(i)), h%missing(i))
      h%phase(i) = IOR(array(i), h%phase(i))
    end do
  end subroutine setOneBits
  
    !---------------------------------------------------------------------------
    !> @brief	Sets the snps in a haplotype to zero
    !> @detail  Sets a snp to be zero if the same position in the bit array is set
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setZeroBits(h, array)
    class(Haplotype) :: h
    integer(kind=8), dimension(:), intent(in) :: array
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(NOT(array(i)), h%missing(i))
      h%phase(i) = IAND(NOT(array(i)), h%phase(i))
    end do
  end subroutine setZeroBits
  
    !---------------------------------------------------------------------------
    !> @brief	Sets the snps in the given position to zero
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setZero(h, pos)
    class(Haplotype) :: h
    integer, intent(in) :: pos
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1
    
    h%phase(cursection) = ibclr(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
  end subroutine setZero

    !---------------------------------------------------------------------------
    !> @brief	Sets the snps in the given position to one
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------  
  subroutine setOne(h, pos)
    class(Haplotype) :: h
    integer, intent(in) :: pos
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1
    
    h%phase(cursection) = ibset(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
  end subroutine setOne
  
    !---------------------------------------------------------------------------
    !> @brief	Returns a new haplotype that is a subset of the given one
    !> @date    May 25, 2017
    !> @return  The subset haplotype
    !---------------------------------------------------------------------------  
  function subset(h,start,finish) result (sub)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: start, finish
    
    type(Haplotype) :: sub
    
    integer :: offset, starti, i
    
    sub%length = finish - start + 1

    sub%sections = sub%length / 64 + 1
    sub%overhang = 64 - (sub%length - (sub%sections - 1) * 64)
    
    allocate(sub%phase(sub%sections))
    allocate(sub%missing(sub%sections))
    
    starti = (start - 1) / 64  + 1
    offset = mod(start - 1, 64)
    
    if (offset == 0) then
      do i = 1, sub%sections
	sub%phase(i) = h%phase(i + starti - 1)
	sub%missing(i) = h%missing(i + starti - 1)
      end do
    else
      do i = 1, sub%sections - 1
	sub%phase(i) = IOR(ISHFT(h%phase(i + starti - 1),-offset), ISHFT(h%phase(i + starti), 64 - offset))
	sub%missing(i) = IOR(ISHFT(h%missing(i + starti - 1),-offset), ISHFT(h%missing(i + starti), 64 - offset))
      end do
      
      if (sub%sections + starti - 1 == h%sections) then
	sub%phase(sub%sections) = ISHFT(h%phase(sub%sections + starti - 1),offset)
	sub%missing(sub%sections) = ISHFT(h%missing(sub%sections + starti - 1),offset)
      else
	sub%phase(sub%sections) = IOR(ISHFT(h%phase(sub%sections + starti - 1),-offset), ISHFT(h%phase(sub%sections + starti), 64 - offset))
	sub%missing(sub%sections) = IOR(ISHFT(h%missing(sub%sections + starti - 1),-offset), ISHFT(h%missing(sub%sections + starti), 64 - offset))
      end if
    end if    
          
    do i = 64 - sub%overhang, 63
      sub%phase(sub%sections) = ibclr(sub%phase(sub%sections), i)
      sub%missing(sub%sections) = ibclr(sub%missing(sub%sections), i)
    end do
  end function subset

    !---------------------------------------------------------------------------
    !> @brief	Sets a subset of the  current haplotype 
    !> to be the same as another haplotype
    !> @date    May 25, 2017
    !> @return  The subset haplotype
    !---------------------------------------------------------------------------   
  subroutine setSubset(h, sub, start)
    class(Haplotype) :: h
    class(Haplotype), intent(in) :: sub
    integer, intent(in) :: start
    
    integer(kind=8) :: mask, startmask, endmask, shifted   
    integer :: starti, endi, offset
    integer :: i
    
    starti = (start - 1) / 64  + 1
    endi = (start + sub%length - 1) / 64 + 1
    offset = mod(start - 1, 64)
    
    if (offset == 0) then
      do i = 1, sub%sections - 1
	h%phase(i + starti - 1) = sub%phase(i)
	h%missing(i + starti - 1) = sub%missing(i)
      end do
      mask = 0
      do i = 1, 64 - sub%overhang
	mask = ibset(mask, i - 1)
      end do
      h%phase(sub%sections) = IOR(IAND(mask, sub%phase(sub%sections)), IAND(NOT(mask), h%phase(sub%sections + starti - 1)))
      h%missing(sub%sections) = IOR(IAND(mask, sub%missing(sub%sections)), IAND(NOT(mask), h%missing(sub%sections + starti - 1)))
    else
      startmask = 0
      do i = offset+1, MIN(64, offset + sub%length)
	startmask = ibset(startmask, i - 1)
      end do      
      shifted = ISHFT(sub%phase(1), offset)
      h%phase(starti) = IOR(IAND(startmask, shifted), IAND(NOT(startmask), h%phase(starti)))
      shifted = ISHFT(sub%missing(1), offset)
      h%missing(starti) = IOR(IAND(startmask, shifted), IAND(NOT(startmask), h%missing(starti)))

      do i = starti + 1, endi - 1
	h%phase(i) = IOR(ISHFT(sub%phase(i-starti), offset - 64), ISHFT(sub%phase(i-starti+1), offset))
	h%missing(i) = IOR(ISHFT(sub%missing(i-starti), offset - 64), ISHFT(sub%missing(i-starti+1), offset))
      end do
      
      if (endi > 1) then
	endmask = 0
	if (offset - sub%overhang > 0) then
	  do i = 1, offset - sub%overhang
	    endmask = ibset(endmask, i - 1)
	  enddo
	  shifted = ISHFT(sub%phase(sub%sections), offset - 64)
	  h%phase(endi) = IOR(IAND(endmask, shifted), IAND(NOT(endmask), h%phase(endi)))
	  shifted = ISHFT(sub%missing(sub%sections), offset - 64)
	  h%missing(endi) = IOR(IAND(endmask, shifted), IAND(NOT(endmask), h%missing(endi)))
	else
	  do i = 1, 64 + offset - sub%overhang
	    endmask = ibset(endmask, i - 1)
	  end do
	  shifted = IOR(ISHFT(sub%phase(sub%sections-1), offset - 64), ISHFT(sub%phase(sub%sections), offset))
	  h%phase(endi) = IOR(IAND(endmask, shifted), IAND(NOT(endmask), h%phase(endi)))
	  shifted = IOR(ISHFT(sub%missing(sub%sections-1), offset - 64), ISHFT(sub%missing(sub%sections), offset))
	  h%missing(endi) = IOR(IAND(endmask, shifted), IAND(NOT(endmask), h%missing(endi)))
	endif
      end if
    end if
    
  end subroutine setSubset
  
    !---------------------------------------------------------------------------
    !> @brief	Tests whether two haplotypes are the same. 
    !> @date    May 25, 2017
    !> @return  Whether the two haplotypes are the same
    !---------------------------------------------------------------------------  
  function equalHap(h1, h2) result (equal)
    class(Haplotype) :: h1, h2
    
    logical :: equal
    
    integer :: i
    
    equal = h1%compareHaplotype(h2)
  end function equalHap

    !---------------------------------------------------------------------------
    !> @brief	Tests whether one haplotype is a "subset" of another.
    !> @detail	Tests whether all non-missing snps in h2 are also non-missing
    !>		    in h1
    !> @date    May 25, 2017
    !> @return  Whether h2 is a "subset" of h1
    !---------------------------------------------------------------------------  
  function isSubset(h1, h2) result(is)
    class(Haplotype) :: h1, h2
    
    logical :: is
    
    integer :: i
    
    is = .true.
    do i = 1, h1%sections
      if (h1%missing(i) .and. .not.(h2%missing(i))) then
	is = .false.
	exit
      end if
    end do
  end function isSubset

    !---------------------------------------------------------------------------
    !> @brief	Sets any errors in the haplotype to be missing 
    !> @date    May 25, 2017
    !---------------------------------------------------------------------------     
  subroutine setErrorToMissing(h)
    class(Haplotype) :: h
    
    integer :: i
    
    do i = 1, h%sections
      h%phase(i) = IAND(NOT(h%missing(i)), h%phase(i))
    end do
  end subroutine setErrorToMissing


  subroutine writeFormattedHaplotype(dtv, unit, iotype, v_list, iostat, iomsg)
      class(Haplotype), intent(in) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: v_list(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      integer(kind=1), dimension(:), allocatable :: array


      array = dtv%toIntegerArray()
      write(unit, "(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)", iostat = iostat, iomsg = iomsg) array
    
  end subroutine writeFormattedHaplotype


      subroutine writeunFormattedHaplotype(dtv, unit, iostat, iomsg)
        class(Haplotype), intent(in)::dtv
        integer, intent(in):: unit
        integer, intent(out) :: iostat      ! non zero on error, etc.
        character(*), intent(inout) :: iomsg  ! define if iostat non zero.
  
        integer(kind=1), dimension(:), allocatable :: array

        array = dtv%toIntegerArray()
        write(unit, "(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)", iostat = iostat, iomsg = iomsg) array

      end subroutine writeunFormattedHaplotype

         subroutine readUnformattedHaplotype(dtv, unit, iostat, iomsg)
        class(haplotype), intent(inout) :: dtv         ! Object to write.
        integer, intent(in) :: unit         ! Internal unit to write to.
        integer, intent(out) :: iostat      ! non zero on error, etc.
        character(*), intent(inout) :: iomsg  ! define if iostat non zero.
        integer(kind=1), dimension(:), allocatable :: array


        read(unit,"(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)",iostat=iostat, iomsg=iomsg) array

        call wrapper(dtv, array)
    end subroutine readUnformattedHaplotype



    pure subroutine wrapper(hap, array) 

      type(haplotype), intent(out) :: hap
      integer(kind=1), dimension(:),intent(in), allocatable :: array

      hap = newHaplotypeInt(array)
    
    end subroutine wrapper


    subroutine readFormattedHaplotype(dtv, unit, iotype, vlist, iostat, iomsg)
      class(Haplotype), intent(inout) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: vlist(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.
      integer(kind=1), dimension(:), allocatable :: array
      read(unit,iotype, iostat=iostat, iomsg=iomsg) array

      call wrapper(dtv, array)
    end subroutine readFormattedHaplotype

    
end module HaplotypeModule
  