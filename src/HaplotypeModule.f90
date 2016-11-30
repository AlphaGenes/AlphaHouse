module HapMod
  implicit none
!  private
 public 
  !! This should go in a constants module but for now
  integer, parameter :: MissingPhaseCode = 9
  
  type HaplotypeType
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Phase       Phase   Missing !
    ! 0           0       0       !
    ! 1           1       0       !
    ! Missing     0       1       !
    ! Not Allowed 1       1       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), dimension(:), pointer :: phase
    integer(kind=8), dimension(:), pointer :: missing
    integer :: sections
    integer :: overhang
    integer :: length
  contains
    procedure :: toIntegerArray
    procedure :: getPhaseMod
    procedure :: overlapMod
    procedure :: mismatchesMod
    procedure :: compatibleMod
    procedure :: mergeMod
    procedure :: numberMissing
    procedure :: compareHaplotype
  end type HaplotypeType
  
!  interface Haplotype
!    module procedure newHaplotypeInt
!    module procedure newHaplotypeBits
!  end interface Haplotype
  
!  interface operator ( == )
!    module procedure compareHaplotype
!  end interface operator ( == )
  
contains
  
  subroutine newHaplotypeN(hap, h)
    integer(kind=1), dimension(:), intent(in) :: hap
    
    type(HaplotypeType) :: h
    
    integer :: nSnps
    integer :: i, cursection, curpos
    
    h%length = size(hap,1)
    
    h%sections = nSnps / 64 + 1
    h%overhang = 64 - (nSnps - (h%sections - 1) * 64)
    
    allocate(h%phase(h%sections))
    allocate(h%missing(h%sections))
    cursection = 1
    curpos = 1
    h%phase = 0
    h%missing = 0
    do i = 1, h%length
      select case (hap(i))
      case (1)
	h%phase(cursection) = ibset(h%phase(cursection), curpos)
      case default
	h%missing(cursection) = ibset(h%missing(cursection), curpos)	
      end select
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end subroutine newHaplotypeN
  
  function newHaplotypeBits(phase, missing, length) result(h)
    integer(kind=8), dimension(:), pointer, intent(in) :: phase, missing
    integer :: length
    
    type(HaplotypeType) :: h
    
    integer :: i
    
    h%phase => phase
    h%missing => missing
    h%length = length
    h%sections = size(phase,1)
    h%overhang = 64 - (h%length - (h%sections - 1) * 64)
    
    do i = 64 - h%overhang + 1, 64
      h%phase(h%sections) = ibclr(h%phase(h%sections), i)
      h%missing(h%sections) = ibclr(h%phase(h%sections), i)
    end do
  end function newHaplotypeBits
  
  function toIntegerArray(h) result(array)
    class(HaplotypeType), intent(in) :: h
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(h%length))
    
    cursection = 1
    curpos = 1
    do i = 1, h%length
      if (btest(h%missing(cursection),curpos)) then
	array(i) = MissingPhaseCode
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
    class(HaplotypeType), intent(in) :: h1, h2
    
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
    class(HaplotypeType), intent(in) :: h
    integer, intent(in) :: pos
    
    integer :: phase
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64
  
    if (btest(h%missing(cursection),curpos)) then
	phase = MissingPhaseCode
    else
      if (btest(h%phase(cursection),curpos)) then
	phase = 1
      else
	phase = 0
      end if
    end if
  end function getPhaseMod
  
  function overlapMod(h1, h2) result (num)
    class(HaplotypeType), intent(in) :: h1, h2
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND(NOT(h1%missing(i)), NOT(h2%missing(i))))
    end do
    
    num = num - h1%overhang
  end function overlapMod
  
  function mismatchesMod(h1, h2) result (num)
    class(HaplotypeType), intent(in) :: h1, h2
    
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h1%sections
      num = num + POPCNT(IAND( IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), &
	    IXOR(h1%phase(i), h2%phase(i)) ))
    end do
  end function mismatchesMod
  
  function compatibleMod(h1, h2, allowedMismatches, minOverlap) result(c)
    class(HaplotypeType), intent(in) :: h1, h2
    integer, intent(in) :: allowedMismatches, minOverlap
    
    logical :: c
    
    c = ((mismatchesMod(h1,h2) <= allowedMismatches) .and. (overlapMod(h1,h2) >= minOverlap))
  end function compatibleMod
  
  function mergeMod(h1,h2) result(h)
    class(HaplotypeType), intent(in) :: h1, h2
    
    type(HaplotypeType) :: h
    
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
    class(HaplotypeType) :: h
    
    logical :: fully
    
    integer :: i
    
    fully = .true.
    
    do i = 1, h%sections
      fully = fully .and. (h%missing(i) == 0)
    end do
  end function fullyPhased
  
  function numberMissing(h) result (num)
    class(HaplotypeType), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
      num = num + POPCNT(h%missing(i))
    end do
    
  end function numberMissing
    
end module HapMod
  