module HaplotypeModule
  
    use constantModule, only: MissingPhaseCode,ErrorPhaseCode
    implicit none
    
    
    type :: Haplotype
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Phase       Phase   Missing !
      ! 0           0       0       !
      ! 1           1       0       !
      ! Missing     0       1       !
      ! Error       1       1       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer(kind=8), dimension(:), pointer :: phase
      integer(kind=8), dimension(:), pointer :: missing
      integer :: sections
      integer :: overhang
      integer :: length
    contains
    procedure :: toIntegerArray => haplotypeToIntegerArray
    procedure :: toIntegerArrayWithErrors
    procedure :: getPhaseMod
    procedure :: setPhaseMod
    procedure :: overlapMod
    procedure :: mismatchesMod
    procedure :: compatibleMod
    procedure :: mergeMod
    procedure :: numberMissing
    procedure :: numberMissingOrError
    procedure :: numberNotMissing
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
  end type Haplotype
  
  interface Haplotype
    module procedure newHaplotypeInt
    module procedure newHaplotypeBits
    module procedure newHaplotypeMissing
    module procedure newHaplotypeHaplotype
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
  
  function newHaplotypeBits(phase, missing, length) result(h)
    integer(kind=8), dimension(:), pointer, intent(in) :: phase, missing
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
    
    do i = 64 - h%overhang + 1, 64
        h%phase(h%sections) = ibclr(h%phase(h%sections), i - 1)
        h%missing(h%sections) = ibclr(h%missing(h%sections), i - 1)
    end do
    
  end function newHaplotypeBits
  
  function newHaplotypeHaplotype(oh) result(h)
    class(Haplotype) :: oh
    
    type(Haplotype) :: h
    
    h = Haplotype(oh%phase, oh%missing, oh%length)
  end function newHaplotypeHaplotype
  
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
        h%phase(h%sections) = ibclr(h%phase(h%sections), i - 1)
        h%missing(h%sections) = ibclr(h%missing(h%sections), i - 1)
    end do
  end function newHaplotypeMissing
  
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
  end function getPhaseMod
  
  subroutine setPhaseMod(h, pos, phase)
    class(Haplotype), intent(in) :: h
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
  
  function numberMissingOrError(h) result (num)
    class(Haplotype), intent(in) :: h
        
    integer :: num
    
    integer :: i
    
    num = 0
    
    do i = 1, h%sections
        num = num + POPCNT(h%missing(i))
    end do
    
  end function numberMissingOrError
  
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
    num = num - h1%overhang
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
   curpos = pos - (cursection - 1) * 64 - 1


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
  
  subroutine setFromOther(h, oh)
    class(Haplotype) :: h
    type(Haplotype) :: oh
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(h%missing(i), oh%missing(i))
      h%phase(i) = IOR(IAND(oh%missing(i),h%phase(i)),IAND(NOT(oh%missing(i)),oh%phase(i)))
    end do
  end subroutine setFromOther
  
  subroutine setFromOtherIfMissing(h, oh)
    class(Haplotype) :: h
    type(Haplotype) :: oh
    
    integer :: i
    
    do i = 1, h%sections
      h%phase(i) = IOR(IAND(NOT(h%missing(i)),h%phase(i)), IAND(h%missing(i), oh%phase(i)))
      h%missing(i) = IAND(h%missing(i), oh%missing(i))      
    end do
  end subroutine setFromOtherIfMissing
  
  subroutine setOneBits(h, array)
    class(Haplotype) :: h
    integer(kind=8), dimension(:), intent(in) :: array
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(NOT(array(i)), h%missing(i))
      h%phase(i) = IOR(array(i), h%phase(i))
    end do
  end subroutine setOneBits
  
  subroutine setZeroBits(h, array)
    class(Haplotype) :: h
    integer(kind=8), dimension(:), intent(in) :: array
    
    integer :: i
    
    do i = 1, h%sections
      h%missing(i) = IAND(NOT(array(i)), h%missing(i))
      h%phase(i) = IAND(NOT(array(i)), h%phase(i))
    end do
  end subroutine setZeroBits
  
  subroutine setZero(h, pos)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: pos
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1
    
    h%phase(cursection) = ibclr(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
  end subroutine setZero
  
  subroutine setOne(h, pos)
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: pos
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1
    
    h%phase(cursection) = ibset(h%phase(cursection), curpos)
    h%missing(cursection) = ibclr(h%missing(cursection), curpos)
  end subroutine setOne
  
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
          
    do i = 64 - sub%overhang + 1, 64
      sub%phase(sub%sections) = ibclr(sub%phase(sub%sections), i - 1)
      sub%missing(sub%sections) = ibclr(sub%missing(sub%sections), i - 1)
    end do
  end function subset
  
  subroutine setSubset(h, sub, start)
    class(Haplotype), intent(in) :: h, sub
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
  
  function equalHap(h1, h2) result (equal)
    class(Haplotype) :: h1, h2
    
    logical :: equal
    
    integer :: i
    
    equal = .true.
    do i = 1, h1%sections
      equal = equal .and. (h1%phase(i) == h2%phase(i)) .and. (h1%missing(i) == h2%missing(i))
    end do
  end function equalHap
    
end module HaplotypeModule
  