module GenotypeModule
  implicit none
  private
  
  !! This should go in a constants module but for now
  integer, parameter :: MissingGenotypeCode = 3
  
  type, public :: Genotype
    private
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Genotype    Homo    Additional !
    ! 0           1       0          !
    ! 1           0       0          !
    ! 2           1       1          !
    ! Missing     0       1          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=8), dimension(:), pointer :: homo
    integer(kind=8), dimension(:), pointer :: additional
    integer :: sections
    integer :: overhang
    integer :: length
  contains
    procedure :: toIntegerArray
    procedure :: getGenotype
    procedure :: compatibleHaplotypes
    procedure :: compatibleHaplotype
    procedure :: numIncommon
    procedure :: mismatches
    procedure :: compareGenotype
    procedure :: getLength
    procedure :: complement
  end type Genotype
  
  interface Genotype
    module procedure newGenotypeInt
  end interface Genotype
  
contains
  
  function newGenotypeInt(geno) result (g)
    integer(kind=1), dimension(:), intent(in) :: geno
    
    type(Genotype) :: g
    
    integer :: i, cursection, curpos
    
    g%length = size(geno,1)
    
    g%sections = g%length / 64 + 1
    g%overhang = 64 - (g%length - (g%sections - 1) * 64)
    
    allocate(g%homo(g%sections))
    allocate(g%additional(g%sections))
    cursection = 1
    curpos = 1
    g%homo = 0
    g%additional = 0
    do i = 1, g%length
      select case (geno(i))
      case (0)
	g%homo(cursection) = ibset(g%homo(cursection), curpos)
      case (1)
	! Nothing to do due to defaults
      case (2)
	g%homo(cursection) = ibset(g%homo(cursection), curpos)
	g%additional(cursection) = ibset(g%additional(cursection), curpos)
      case default
	g%additional(cursection) = ibset(g%additional(cursection), curpos)	
      end select
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function newGenotypeInt
  
  function toIntegerArray(g) result(array)
    class(Genotype), intent(in) :: g
    
    integer(kind=1), dimension(:), allocatable :: array
    
    integer :: i, cursection, curpos
    
    allocate(array(g%length))
    
    cursection = 1
    curpos = 1
    do i = 1, g%length
    if (btest(g%homo(cursection),curpos)) then
      if (btest(g%additional(cursection),curpos)) then
	array(i) = 2
      else
	array(i) = 0
      end if
    else
      if (btest(g%additional(cursection),curpos)) then
	array(i) = MissingGenotypeCode
      else
	array(i) = 1
      end if
    end if
      
      curpos = curpos + 1
      if (curpos == 65) then
	curpos = 1
	cursection = cursection + 1
      end if
    end do
  end function toIntegerArray
  
  function compareGenotype(g1, g2) result(same)
    class(Genotype), intent(in) :: g1, g2
    
    logical :: same
    
    integer :: i
    
    if (g1%length /= g2%length) then
      same = .false.
    else
      same = .true.
      do i = 1, g1%sections
	same = same .and. (g1%homo(i) == g2%homo(i)) .and. (g1%additional(i) == g2%additional(i))
      end do
    end if
  end function compareGenotype
  
  function getGenotype(g, pos) result (genotype)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos
    
    integer :: genotype
    
    integer :: cursection, curpos
    
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64
  
    if (btest(g%homo(cursection),curpos)) then
      if (btest(g%additional(cursection),curpos)) then
	genotype = 2
      else
	genotype = 0
      end if
    else
      if (btest(g%additional(cursection),curpos)) then
	genotype = MissingGenotypeCode
      else
	genotype = 1
      end if
    end if
  end function getGenotype
  
  function numOppose(g1, g2) result(num)
    class(Genotype), intent(in) :: g1, g2
    
    integer :: num
    
    integer :: i
    
    do i = 1, g1%sections
	num = num + POPCNT(IAND(IAND(g1%homo(i), g2%homo(i)), &
	IEOR(g1%additional(i), g1%additional(i))))
    end do
  end function numOppose
  
  function numIncommon(g1, g2) result(num)
    class(Genotype), intent(in) :: g1, g2
    
    integer :: num
    
    integer :: i
    
    num = 0
    do i = 1, g1%sections
	num = num + POPCNT(IAND(IOR(g1%homo(i), NOT(g1%additional(i))), &
	IOR(g2%homo(i), NOT(g2%additional(i)))))
    end do
    
    num = num - g1%overhang
  end function numIncommon
  
  function complement(g,h) result(c)
    use HaplotypeModule
    class(Genotype), intent(in) :: g
    class(Haplotype), intent(in) :: h
    
    type(Haplotype) :: c
    
    integer :: i
    integer(kind=8), dimension(:), pointer :: phase, missing
    
    allocate(phase(g%sections))
    allocate(missing(g%sections))
    do i = 1, g%sections
      !phase = IOR( IAND(h%phase(i), h%missing(i)), &
      phase(i) = IOR( h%missing(i), &
		   IOR( IAND(IAND(h%phase(i),NOT(h%missing(i))),g%homo(i)), &
		        IAND(NOT(IOR(h%phase(i),h%missing(i))), NOT(IEOR(g%homo(i), g%additional(i))))))
      missing(i) = IOR (h%missing(i), &
		    IOR( IAND(NOT(h%phase(i)),g%additional(i)), &
			 IAND(h%phase(i),IEOR(g%homo(i),g%additional(i))) ) )
    end do
    
    c = Haplotype(phase,missing,g%length)
  end function complement
    
  function compatibleHaplotypes(g, h1, h2, threshold) result(c)
    use HaplotypeModule
    class(Genotype), intent(in) :: g
    class(Haplotype), intent(in) :: h1, h2
    integer, intent(in) :: threshold
    
    logical :: c
    
    integer :: i, num
    
    num = 0
    
    do i = 1, g%sections
      
      num = num + POPCNT( &
	IOR(IOR( &
	! Genotype is zero !
	IAND( &
	  ! Genotype is 0
	  IAND(g%homo(i), NOT(g%additional(i))), &
	  ! One of the haplotypes is 1
	  IOR( h1%phase(i), h2%phase(i)) ), &

	! Genotype is one !
	IAND( &
	  ! Genotype is 1
	  IAND(NOT(g%homo(i)), NOT(g%additional(i))), &
	  ! Both of the haplotypes are not missing but are not opposed
	  IAND(IAND(NOT(h1%missing(i)), NOT(h2%missing(i))), NOT(IEOR(h1%phase(i), h2%phase(i)))))), &

	! Genotype is two !
	IAND( &
	  ! Genotype is 2
	  IAND(g%homo(i), g%additional(i)), &
	  ! One of the haplotypes is 0
	  IOR(IAND(NOT(h1%phase(i)), NOT(h1%missing(i))), IAND(NOT(h2%phase(i)), NOT(h2%missing(i)))))))
    end do
    
    num = num - g%overhang
	
    c = num <= threshold
  end function compatibleHaplotypes
  
  !! THIS WILL NEED MIN OVERLAP AT SOME POINT !!
  function compatibleHaplotype(g, h, threshold) result(c)
    use HaplotypeModule
    class(Genotype), intent(in) :: g
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: threshold
    
    logical :: c
    
    integer :: i, num
    
    num = 0
    
    do i = 1, g%sections
      num = num + POPCNT( IAND( &
	!! Genotype is zero !!
      	IAND( &
	  ! Genotype is 0
	  IAND(g%homo(i), NOT(g%additional(i))), &
	  ! The haplotypes is 1
	  h%phase(i) ), &
	!! Genotype is two
	IAND( &
	  ! Genotype is 2
	  IAND(g%homo(i), g%additional(i)), &
	  ! One of the haplotypes is 0
	  IAND(NOT(h%phase(i)), NOT(h%missing(i))))))
    end do
    
    c = num <= threshold
  end function compatibleHaplotype
  
  function mismatches(g1, g2) result(c)
    class(Genotype), intent(in) :: g1, g2
    integer :: c
    
    integer :: i

    c = 0
    do i = 1, g1%sections
	c = c + POPCNT(IAND(IAND(g1%homo(i), g2%homo(i)), &
	IEOR(g1%additional(i), g2%additional(i))))
    end do
  end function mismatches
  
  function getLength(g) result(l)
    class(Genotype), intent(in) :: g
    integer :: l
    
    l = g%length
  end function getLength
  
end module GenotypeModule
  

