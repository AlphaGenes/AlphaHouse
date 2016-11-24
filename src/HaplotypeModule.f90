module HaplotypeModule
  implicit none
  
  !! This should go in a constants module but for now
  integer, parameter :: MissingPhaseCode = 9
  
  type Haplotype
    private
    integer(kind=8), dimension(:), pointer :: phase
    integer(kind=8), dimension(:), pointer :: missing
    integer :: sections
    integer :: overhang
    integer :: length
  contains
    procedure :: toIntegerArray
  end type Haplotype
  
  interface Haplotype
    module procedure newHaplotypeInt
  end interface Haplotype
  
  interface operator ( == )
    module procedure compareHaplotype
  end interface operator ( == )
  
contains
  
  function newHaplotypeInt(hap) result (h)
    integer(kind=1), dimension(:), intent(in) :: hap
    
    type(Haplotype) :: h
    
    integer :: nSnps
    integer :: i, cursection, curpos
    
    h%length = size(hap,1)
    
    h%sections = nSnps / 64 + 1
    h%overhang = 64 - (nSnps - (h%sections - 1) * 64)
    
    allocate(h%phase(h%sections))
    cursection = 1
    curpos = 1
    h%phase = 0
    h%missing = 0
    do i = 1, h%length
      select case (hap(i))
!      case (0)
!	bits(cursection) = ibclr(bits(cursection), curpos)
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
  end function newHaplotypeInt
  
  function toIntegerArray(h) result(array)
    class(Haplotype), intent(in) :: h
    
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
  
end module HaplotypeModule
  