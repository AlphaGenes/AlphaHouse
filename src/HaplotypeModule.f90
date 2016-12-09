module HaplotypeModule
  implicit none
  private
  !! This should go in a constants module but for now
  integer, parameter :: MissingPhaseCode = 9
  integer, parameter :: ErrorPhaseCode = -1
  
  type, public :: Haplotype
    private
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Phase       Phase   Missing !
    ! 0           0       0       !
    ! 1           1       0       !
    ! Missing     0       1       !
    ! Error       1       1       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), dimension(:), pointer, public :: phase
    integer(kind=8), dimension(:), pointer, public :: missing
    integer, public :: sections
    integer :: overhang
    integer :: length
  contains
    procedure :: toIntegerArray
    procedure :: getPhaseMod
    procedure :: setPhaseMod
    procedure :: overlapMod
    procedure :: mismatchesMod
    procedure :: compatibleMod
    procedure :: mergeMod
    procedure :: numberMissing
    procedure :: numberNotMissing
    procedure :: numberError
    procedure :: compareHaplotype
    procedure :: numberSame
    procedure :: fullyPhased
    procedure :: setUnphased
    procedure :: getLength
    procedure :: isMissing
    procedure :: numberBothNotMissing
  end type Haplotype
  
  interface Haplotype
    module procedure newHaplotypeInt
    module procedure newHaplotypeBits
    module procedure newHaplotypeMissing
  end interface Haplotype
  
contains
  
  function newHaplotypeInt(hap) result(h)
    integer(kind=1), dimension(:), intent(in) :: hap
    
    type(Haplotype) :: h
    
    integer :: i, cursection, curpos
    
    h%length = size(hap,1)
    
    h%sections = h%length / 64 + 1
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    cursection = 1
    curpos = 1
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
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function newHaplotypeInt
  
  function newHaplotypeBits(phase, missing, length) result(h)
    integer(kind=8), dimension(:), pointer, intent(in) :: phase, missing
    integer :: length
    
    type(Haplotype) :: h
    
    integer :: i
    
    h%phase => phase
    h%missing => missing
    h%length = length
    h%sections = size(phase,1)
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    
    do i = 64 - h%overhang + 1, 64
      h%phase(h%sections) = ibclr(h%phase(h%sections), i)
      h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
    
  end function newHaplotypeBits
  
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

    
    do i = 64 - h%overhang + 1, 64
      h%phase(h%sections) = ibclr(h%phase(h%sections), i)
      h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
  end function newHaplotypeMissing
  
  function toIntegerArray(h) result(array)
    class(Haplotype), intent(in) :: h
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(h%length))
    
    cursection = 1
    curpos = 1
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
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function toIntegerArray
  
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
  
  function getPhaseMod(h, pos) result (phase)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: pos
    
    integer(kind=1) :: phase
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64
  
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
  end function getPhaseMod
  
  subroutine setPhaseMod(h, pos, phase)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: pos
    integer(kind=1) :: phase

    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64    
    
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
  end subroutine setPhaseMod
  
  function overlapMod(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND(NOT(h1%missing(i)), NOT(h2%missing(i))))
    end do
    
    num = num - h1%overhang
  end function overlapMod
  
  function mismatchesMod(h1, h2) result (num)
    class(Haplotype), intent(in) :: h1, h2
    
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND( IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
	    IXOR(h1%phase(i), h2%phase(i)) ))
    end do
  end function mismatchesMod
  
  function compatibleMod(h1, h2, allowedMismatches, minOverlap) result(c)
    class(Haplotype), intent(in) :: h1, h2
    integer, intent(in) :: allowedMismatches, minOverlap
    
    logical :: c
    
    c = ((mismatchesMod(h1,h2) <= allowedMismatches) .and. (overlapMod(h1,h2) >= minOverlap))
  end function compatibleMod
  
  function mergeMod(h1,h2) result(h)
    class(Haplotype), intent(in) :: h1, h2
    
    type(Haplotype) :: h
    
    integer :: i
    
    h%sections = h1%sections
    h%overhang = h1%overhang
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    
    do i = 1, h1%sections
      h%missing(i) = IOR( &
	!Both not missing but opposed
	IAND(IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
	IXOR(h1%phase(i), h2%phase(i))), &
	!Both missing
	IAND(h1%missing(i),h2%missing(i)) &
	)
      h%phase(i) = IAND( &
        ! Not missing (phase should always be zero if missing)
	NOT(h%missing(i)), &
	! One of the phases is 1 (no need to test for missing as phase is only one if not missing)
	IOR(h1%phase(i), h2%phase(i)) &
	)
    end do
  end function mergeMod
  
  function fullyPhased(h) result(fully)
    class(Haplotype) :: h
    
    logical :: fully
    
    integer :: i
    
    fully = .true.
    
    do i = 1, h%sections
      fully = fully .and. (h%missing(i) == 0)
    end do
  end function fullyPhased
  
  function numberMissing(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
      num = num + POPCNT(IAND(h%missing(i), NOT(h%phase(i))))
    end do
    
  end function numberMissing
  
  function numberError(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
      num = num + POPCNT(IAND(h%missing(i), h%phase(i)))
    end do
    
  end function numberError
  
  function numberNotMissing(h) result(num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    num = h%length - h%numberMissing()
  end function numberNotMissing
  
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
  end function numberSame
  
  subroutine setUnphased(h)
    class(Haplotype) :: h
    
    integer :: i
    
    h%phase = 0
    h%missing = 0
    h%missing = NOT(h%missing)
    
    do i = 64 - h%overhang + 1, 64
      h%phase(h%sections) = ibclr(h%phase(h%sections), i)
      h%missing(h%sections) = ibclr(h%missing(h%sections), i)
    end do
  end subroutine setUnphased
  
  function getLength(h) result(l)
    class(Haplotype), intent(in) :: h
    integer :: l
    
    l = h%length
  end function getLength
  
  function isMissing(h, pos) result (missing)
   class(Haplotype), intent(in) :: h
   integer, intent(in) :: pos

   logical :: missing

   integer :: cursection, curpos

   cursection = (pos-1) / 64 + 1
   curpos = pos - (cursection - 1) * 64


   missing = BTEST(h%missing(cursection), curpos)
  end function isMissing
  
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
    
end module HaplotypeModule
  