module GenotypeModule
  implicit none
  
  !! This should go in a constants module but for now
  integer, parameter :: MissingGenotypeCode = 3
  
  type Genotype
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
  end type Genotype
  
  interface Genotype
    module procedure newGenotypeInt
  end interface Genotype
  
  interface operator ( == )
    module procedure compareGenotype
  end interface operator ( == )
  
contains
  
  function newGenotypeInt(geno) result (g)
    integer(kind=1), dimension(:), intent(in) :: geno
    
    type(Genotype) :: g
    
    integer :: nSnps
    integer :: i, cursection, curpos
    
    g%length = size(geno,1)
    
    g%sections = nSnps / 64 + 1
    g%overhang = 64 - (nSnps - (g%sections - 1) * 64)
    
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
      phase = IAND (NOT(h%missing(i)), &
		    IOR( IAND(h%phase(i),AND(g%homo(i),g%additional(i))), &
			 NOT(IOR(h%phase(i),IOR(g%homo(i),g%additional(i))))))
      missing = IOR (h%missing(i), &
		    IOR( IAND(g%additional(i),NOT(h%phase(i))), &
			 IAND(h%phase(i),IEOR(g%homo(i),g%additional(i)))))
    end do
    
    c = newHaplotypeBits(phase,missing,g%length)
  end function complement
    
  
end module GenotypeModule
  

