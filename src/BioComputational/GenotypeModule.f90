
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     GenotypeModule.f90
!
! DESCRIPTION:
!> @brief    Module cotaining a type representing a Genotype and associated
!>	     procedures
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!> @author   Daniel Money, daniel.money@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 dmoney - initial version
! 2017-04-28 dwilso18 - revision with comments
!
!-------------------------------------------------------------------------------




module GenotypeModule

    use constantModule, only : MissingGenotypeCode
    use iso_fortran_env
  implicit none

    !---------------------------------------------------------------------------
    !> @brief   Represents genotypes.  
    !> @detail	Stored as bit arrays for computational efficiency.
    !> @date    May 27, 2017
    !---------------------------------------------------------------------------
  type :: Genotype
  !Reminder as to how the data is stored...
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Genotype    Homo    Additional !
  ! 0           1       0          !
  ! 1           0       0          !
  ! 2           1       1          !
  ! Missing     0       1          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=int64), allocatable, dimension(:) :: homo
  integer(kind=int64), allocatable, dimension(:) :: additional
  integer :: sections
  integer :: overhang
  integer :: length
  contains
  procedure :: toIntegerArray => genotypeToIntegerArray
  procedure :: getGenotype
  procedure :: setGenotype
  procedure :: compatibleHaplotypes
  procedure :: compatibleHaplotype
  procedure :: numIncommon
  procedure :: numMissing
  procedure :: mismatches
  procedure :: compareGenotype
  procedure :: getLength
  procedure :: complement
  procedure :: numNotMissing
  procedure :: setHaplotypeFromGenotype
  procedure :: isZero
  procedure :: isTwo
  procedure :: isMissing
  procedure :: isHomo
  procedure :: numberErrors
  procedure :: getErrors
  procedure :: setHaplotypeFromGenotypeIfError
  procedure :: setHaplotypeFromGenotypeIfMissing
  procedure :: numberErrorsSingle
  procedure :: getErrorsSingle
  procedure :: subset
  procedure :: setFromHaplotypesIfMissing
  procedure :: setFromOtherIfMissing
  procedure :: readFormattedGenotype
  procedure :: readunFormattedGenotype
  procedure :: writeFormattedGenotype
  procedure :: writeunFormattedGenotype
  final :: destroyGenotype
  generic:: write(formatted)=> writeFormattedGenotype
  generic:: write(unformatted)=> writeunFormattedGenotype
  generic:: read(formatted) => readFormattedGenotype
  generic:: read(unformatted) => readunFormattedGenotype
  end type Genotype

    !---------------------------------------------------------------------------
    !> @brief   Represents mendelian consistencies and inconsistencies for a trio.
    !> @detail  Stores 6 bit arrays representing consistent / inconsistent snps
    !>          for each of paternal, maternal and individual.
    !>
    !>          Consistent is where we have the information to know a snp is
    !>          consistent.  Similarly for inconsistent.  If neither is set for
    !>          a snp that means there is missing information and so consistentcy
    !>          can not be determined.
    !>
    !>          Paternal represents where the paternal and individual snps are
    !>          consistent or inconsistent.  This includes the case where the
    !>          individual is a het and the parents are inconsistent.  Maternal
    !>          is defined similarly.  Individual is whether an individual is
    !>          consistent or inconsistent with both it's parents.
    !> @date    May 27, 2017
    !---------------------------------------------------------------------------
    type :: Mendelian
        integer(kind=int64), dimension(:), allocatable :: paternalInconsistent, maternalInconsistent, individualInconsistent
        integer(kind=int64), dimension(:), allocatable :: paternalConsistent, maternalConsistent, individualConsistent
    end type Mendelian

  interface Genotype
    module procedure newGenotypeInt
    module procedure newGenotypeHap
    module procedure newGenotypeMissing
  end interface Genotype

      contains


      subroutine destroyGenotype(g)
        type(Genotype) :: g
        if (allocated(g%homo)) then
          deallocate(g%homo)
          deallocate(g%additional)
        endif
    end subroutine destroyGenotype
    !---------------------------------------------------------------------------
    !> @brief	Constructs a new Genotype from a integer array
    !> @date    November 26, 2016
    !> @return	New genotype object 
    !---------------------------------------------------------------------------
      pure function newGenotypeInt(geno) result (g)
        integer(kind=1), dimension(:), intent(in) :: geno

        type(Genotype) :: g

        integer :: i, cursection, curpos

        g%length = size(geno,1)

        g%sections = (g%length - 1) / 64 + 1
        g%overhang = 64 - (g%length - (g%sections - 1) * 64)

        allocate(g%homo(g%sections))
        allocate(g%additional(g%sections))
        cursection = 1
        curpos = 0
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
            if (curpos == 64) then
                curpos = 0
                cursection = cursection + 1
            end if
        end do
    end function newGenotypeInt
    
    !---------------------------------------------------------------------------
    !> @brief	Constructs a new genotype from two haplotype objects
    !> @date    November 26, 2016
    !> @return	New genotype object 
    !---------------------------------------------------------------------------
    function newGenotypeHap(h1, h2) result (g)
      use HaplotypeModule
      type(Haplotype) :: h1, h2
      
      type(Genotype) :: g
      integer :: i
      
      g%length = h1%length
      g%sections = h1%sections
      g%overhang = h1%overhang
      
      allocate(g%homo(g%sections))
      allocate(g%additional(g%sections))
      
      do i = 1, g%sections
	g%homo(i) = IOR( &
	  ! BOTH HOMO MAJOR
			      IAND( IAND(NOT(h1%phase(i)), NOT(h1%missing(i))), IAND(NOT(h2%phase(i)), NOT(h2%missing(i)))), &
	  ! BOTH HOMO MINOR
			      IAND( IAND(h1%phase(i), NOT(h1%missing(i))), IAND(h2%phase(i), NOT(h2%missing(i)))) )
	
	g%additional(i) = IOR ( &
	  ! BOTH HOMO MINOR
			      IAND( IAND(h1%phase(i), NOT(h1%missing(i))), IAND(h2%phase(i), NOT(h2%missing(i)))), &
	  ! EITHER MISSING / ERROR
			      IOR(h1%missing(i), h2%missing(i)) )
      end do
      
      do i = 64 - g%overhang, 63
	g%homo(g%sections) = ibclr(g%homo(g%sections), i)
      end do 
    end function newGenotypeHap

    
  !---------------------------------------------------------------------------
  !> @brief   Constructs a new Genotype with all positions set to missing
  !> @date    May 25, 2017
  !> @return  New haplotype object 
  !--------------------------------------------------------------------------- 
  function newGenotypeMissing(length) result(g)
    integer, intent(in) :: length
    
    type(Genotype) :: g
    
    integer :: i
    
    g%length = length
    g%sections = (g%length - 1) / 64 + 1
    g%overhang = 64 - (g%length - (g%sections - 1) * 64)
    allocate(g%homo(g%sections))
    allocate(g%additional(g%sections))
    g%homo = 0
    g%additional = 0
    g%additional = NOT(g%additional)

    
    do i = 64 - g%overhang, 63
        g%homo(g%sections) = ibclr(g%homo(g%sections), i)
        g%additional(g%sections) = ibclr(g%additional(g%sections), i)
    end do
  end function newGenotypeMissing

  
    !---------------------------------------------------------------------------
    !> @brief	Converts a Genotype to an integer array
    !> @detail	Error states are coded as missing (i.e. 9)
    !> @date    November 26, 2016
    !> @return	Integer array representing the genotype
    !---------------------------------------------------------------------------  
    function genotypeToIntegerArray(g, nsnp) result(array)
        class(Genotype), intent(in) :: g
        integer, intent(in),optional :: nsnp
        integer(kind=1), dimension(:), allocatable :: array

        integer :: i, cursection, curpos
        integer :: iterator !< true number of snps


      
	
        if (present(nsnp)) then

          if (nsnp < g%length) then
            iterator =  nsnp
          else
            iterator =  g%length
          endif
          allocate(array(nsnp))
          
        else
          allocate(array(g%length))
          iterator = g%length
        endif

        array = MissingGenotypeCode
        cursection = 1
        curpos = 0
        do i = 1,iterator
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
	  if (curpos == 64) then
	      curpos = 0
	      cursection = cursection + 1
	  end if
	end do
    end function genotypeToIntegerArray

  !---------------------------------------------------------------------------
  !> @brief	Compares two genotypes
  !> @date	November 26, 2016
  !> @return	Whether the two genotypes are identical
  !---------------------------------------------------------------------------   
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


  !---------------------------------------------------------------------------
  !> @brief	Gets the genotype of a snp at the given position
  !> @date	November 26, 2016
  !> @return	The genotype
  !---------------------------------------------------------------------------  
function getGenotype(g, pos) result (genotype)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos !< snp 

    integer :: genotype

    integer :: cursection, curpos


    if (pos == 0) then 
      genotype = 9
      return
    endif
    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1

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

  !---------------------------------------------------------------------------
  !> @brief	Sets the genotype of a snp at the given position
  !> @date	November 26, 2016
  !---------------------------------------------------------------------------
subroutine setGenotype(g, pos, val)
    class(Genotype), intent(inout) :: g
    integer, intent(in) :: val
    integer, intent(in) :: pos

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1

    select case (val)
      case (0)
          g%homo(cursection) = ibset(g%homo(cursection), curpos)
          g%additional(cursection) = ibclr(g%additional(cursection), curpos)
      case (1)
          ! Nothing to do due to default
          g%homo(cursection) = ibclr(g%homo(cursection), curpos)
          g%additional(cursection) = ibclr(g%additional(cursection), curpos)
      case (2)
          g%homo(cursection) = ibset(g%homo(cursection), curpos)
          g%additional(cursection) = ibset(g%additional(cursection), curpos)
      case default
          g%additional(cursection) = ibset(g%additional(cursection), curpos)
          g%homo(cursection) = ibclr(g%homo(cursection), curpos)
    end select 

end subroutine setGenotype


  !---------------------------------------------------------------------------
  !> @brief   Returns the number of opposing homozygotes between two genotypes
  !> @date    November 26, 2016
  !> @return  The number of opposing homozygotes
  !---------------------------------------------------------------------------
function numOppose(g1, g2) result(num)
    class(Genotype), intent(in) :: g1, g2

    integer :: num

    integer :: i

    do i = 1, g1%sections
        num = num + POPCNT(IAND(IAND(g1%homo(i), g2%homo(i)), &
            IEOR(g1%additional(i), g1%additional(i))))
    end do
end function numOppose

  !---------------------------------------------------------------------------
  !> @brief   Returns the number of snps that are non-misisng in both genotypes
  !> @date    November 26, 2016
  !> @return  The number of snps in common
  !---------------------------------------------------------------------------
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

  !---------------------------------------------------------------------------
  !> @brief   Returns the complement of a haplotype given a genotype
  !> @detail  Calculates what the other haplotype must be given the genotype
  !>		and one haplotype.  Returns error if the genotype and haplotype
  !>		are incompatible or if the haplotype is error. Returns missing 
  !>            if either the genotype or haplotype are missing.  Else returns
  !>		the appropriate phase
  !> @date    May 26, 2017
  !> @return  The complement haplotype
  !---------------------------------------------------------------------------
function complement(g,h) result(c)
    use HaplotypeModule
    class(Genotype), intent(in) :: g
    class(Haplotype), intent(in) :: h

    type(Haplotype) :: c

    integer :: i
    integer(kind=int64), dimension(:), allocatable :: phase, missing

    allocate(phase(g%sections))
    allocate(missing(g%sections))
    do i = 1, g%sections
      phase(i) = IOR( IAND(h%phase(i), h%missing(i)), &
!      phase(i) = IOR( h%missing(i), &
         IOR( IAND(IAND(h%phase(i),NOT(h%missing(i))),g%homo(i)), &
            IAND(NOT(IOR(h%phase(i),h%missing(i))), NOT(IEOR(g%homo(i), g%additional(i))))))
      missing(i) = IOR (h%missing(i), &
        IOR( IAND(NOT(h%phase(i)),g%additional(i)), &
           IAND(h%phase(i),IEOR(g%homo(i),g%additional(i))) ) )
    end do

    c = Haplotype(phase,missing,g%length)
    deallocate(phase)
    deallocate(missing)
end function complement

  !---------------------------------------------------------------------------
  !> @brief   Tests whether two haplotypes are compatible with a genotype
  !> @detail  Allows threhold amount of error.  An error is where genotype is 0,
  !>		either phase is 1 or genotype is 2, either phase is 0.  Genotype
  !>		is one but both haplotypes are the same phase is also an error.
  !>		If the total number of errors is less than or equal to threshold
  !>		return true.
  !> @date    May 26, 2017
  !> @return  Whether the haplotypes are compatible with the genotype
  !---------------------------------------------------------------------------
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

  !---------------------------------------------------------------------------
  !> @brief   Tests whether a haplotype is compatible with a genotype
  !> @detail  Allows threhold amount of error.  An error is where genotype is 0,
  !>		phase is 1 or genotype is 2, phase is 0.  If the total number of
  !>		errors is less than or equal to threshold return true.
  !> @date    May 26, 2017
  !> @return  Whether the haplotype is compatible with the genotype
  !---------------------------------------------------------------------------
function compatibleHaplotype(g, h, threshold) result(c)
!! THIS WILL NEED MIN OVERLAP AT SOME POINT !!
    use HaplotypeModule
    class(Genotype), intent(in) :: g
    class(Haplotype), intent(in) :: h
    integer, intent(in) :: threshold

    logical :: c

    integer :: i, num

  
    num = 0

    do i = 1, g%sections
      num = num + POPCNT( IOR( &
        !! Genotype is zero !!
        IAND( &
          ! Genotype is 0
          IAND(g%homo(i), NOT(g%additional(i))), &
          ! The haplotypes is 1
          IAND(h%phase(i), NOT(h%missing(i))) ), &
        !! Genotype is two
        IAND( &
          ! Genotype is 2
          IAND(g%homo(i), g%additional(i)), &
          ! One of the haplotypes is 0
          IAND(NOT(h%phase(i)), NOT(h%missing(i))))))
    end do
    
  c = num <= threshold
end function compatibleHaplotype


  !---------------------------------------------------------------------------
  !> @brief   Counts the number of mismatches (opposing homozygotes) between
  !>		two genotypes.
  !> @date    May 26, 2017
  !> @return  The number of opposing homozygotes
  !---------------------------------------------------------------------------
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

  !---------------------------------------------------------------------------
  !> @brief   Returns the length of the genotype
  !> @date    May 26, 2017
  !> @return  The length of the genotype
  !---------------------------------------------------------------------------
function getLength(g) result(l)
  ! THIS SHOULD PROBABLY DISAPPEAR AT SOME POINT AS WE DON'T NORMALLY USE GETTERS
  ! BUT IT'S BEING USED...
    class(Genotype), intent(in) :: g
    integer :: l

    l = g%length
end function getLength

  !---------------------------------------------------------------------------
  !> @brief   Returns the number of not missing snps in the genotype
  !> @date    May 26, 2017
  !> @return  The number of not missing snps
  !---------------------------------------------------------------------------
function numNotMissing(g) result(c)
    class(Genotype), intent(in) :: g

    integer :: c

    integer :: i

    c = 0
    do i = 1, g%sections
      c = c + POPCNT(NOT(IAND(NOT(g%homo(i)), g%additional(i))))
  end do
end function numNotMissing

  !---------------------------------------------------------------------------
  !> @brief   Returns the number of missing snps in the genotype
  !> @date    May 26, 2017
  !> @return  The number of missing snps
  !---------------------------------------------------------------------------
function numMissing(g) result(c)
    class(Genotype), intent(in) :: g

    integer :: c

    c = g%length - g%numNotMissing()
end function numMissing

  !---------------------------------------------------------------------------
  !> @brief   Sets the Haplotype from the Genotype
  !> @detail  If the genotypes means the phase is known set the phase in the
  !>		haplotype, else set it to missing
  !> @date    May 26, 2017
  !---------------------------------------------------------------------------
subroutine setHaplotypeFromGenotype(g, h)
    use HaplotypeModule

    class(Haplotype) :: h
    class(Genotype), intent(in) :: g

    integer :: i

    do i = 1, g%sections
      h%missing(i) = NOT(g%homo(i))

      h%phase(i) = IAND(g%homo(i), g%additional(i))
  end do
  do i = 64 - g%overhang, 63
      h%missing(g%sections) = ibclr(h%missing(g%sections), i)
    end do
end subroutine setHaplotypeFromGenotype

  !---------------------------------------------------------------------------
  !> @brief   Sets the Haplotype from the Genotype if error is set true for that snp.
  !> @detail  If error is true and the genotypes means the phase is known
  !>		set the phase in the haplotype, else set it to missing
  !> @date    May 26, 2017
  !---------------------------------------------------------------------------
  subroutine setHaplotypeFromGenotypeIfError(g, h, errors)
    use HaplotypeModule

    type(Haplotype), intent(in), pointer :: h
    class(Genotype), intent(in) :: g
    integer(kind=int64), intent(in), dimension(:) :: errors

    integer :: i

    do i = 1, g%sections
      h%missing(i) = IOR(IAND(NOT(errors(i)), h%missing(i)), IAND(errors(i), NOT(g%homo(i))))

      h%phase(i) = IOR(IAND(NOT(errors(i)), h%phase(i)), IAND(errors(i), IAND(g%homo(i), g%additional(i))))
    end do
    do i = 64 - g%overhang, 63
      h%missing(g%sections) = ibclr(h%missing(g%sections), i)
    end do
  end subroutine setHaplotypeFromGenotypeIfError

  !---------------------------------------------------------------------------
  !> @brief   Sets the Haplotype from the Genotype if the Haplotype is missing
  !>		or error
  !> @detail  If the phase is misisng or error and the genotypes means the 
  !>		phase is known set the phase in the haplotype, else set it to
  !>		missing
  !> @date    May 26, 2017
  !---------------------------------------------------------------------------
  subroutine setHaplotypeFromGenotypeIfMissing(g, h)
    ! NOT SURE THIS DEALS WITH ERRORS SENSIBLY
    use HaplotypeModule

    type(Haplotype) :: h
    class(Genotype), intent(in) :: g

    integer :: i

    do i = 1, g%sections
      h%phase(i) = IOR(IAND(NOT(h%missing(i)), h%phase(i)), IAND(h%missing(i), IAND(g%homo(i), g%additional(i))))
      h%missing(i) = IAND(h%missing(i), NOT(g%homo(i)))
    end do
    do i = 64 - g%overhang, 63
      h%missing(g%sections) = ibclr(h%missing(g%sections), i)
    end do
  end subroutine setHaplotypeFromGenotypeIfMissing

  !---------------------------------------------------------------------------
  !> @brief   Tests whether the given snp has a zero genotype
  !> @date    May 26, 2017
  !> @return  Whether the given snp has a zero genotype
  !---------------------------------------------------------------------------	
  function isZero(g, pos) result (zero)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: zero

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    zero = BTEST(IAND(g%homo(cursection),NOT(g%additional(cursection))), curpos)
end function isZero

  !---------------------------------------------------------------------------
  !> @brief   Tests whether the given snp has a two genotype
  !> @date    May 26, 2017
  !> @return  Whether the given snp has a two genotype
  !---------------------------------------------------------------------------	
function isTwo(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(IAND(g%homo(cursection),g%additional(cursection)), curpos)
end function isTwo

  !---------------------------------------------------------------------------
  !> @brief   Tests whether the given snp is set missing
  !> @date    May 26, 2017
  !> @return  Whether the given snp is set missing
  !---------------------------------------------------------------------------	
function isMissing(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(IAND(NOT(g%homo(cursection)),g%additional(cursection)), curpos)
end function isMissing

  !---------------------------------------------------------------------------
  !> @brief   Tests whether the given snp is a homozygote
  !> @date    May 26, 2017
  !> @return  Whether the given snp is a homozygote
  !---------------------------------------------------------------------------	
function isHomo(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(g%homo(cursection), curpos)
  end function isHomo
  
  !---------------------------------------------------------------------------
  !> @brief   Returns the location of errors between a genotype and a
  !>		haplotypes.
  !> @detail  An error is where the haplotypes and the genotype are
  !>		inconsistent.  Only snps which are non-missing in the 
  !>		haplotypes and the genotype are tested.
  !> @date    May 26, 2017
  !> @return  The location of errors as a bit array
  !--------------------------------------------------------------------------- 
  function getErrorsSingle(g,h) result(errors)
    use HaplotypeModule

    class(Genotype) :: g
    class(Haplotype) :: h

    integer(kind=int64), dimension(:), pointer :: errors

    integer :: i

    integer(kind=int64) :: allnotmissing, zeroerror, twoerror

    allocate(errors(g%sections))

    do i = 1, g%sections
      allnotmissing = IAND(IOR(g%homo(i), NOT(g%additional(i))), NOT(h%missing(i)))

      zeroerror = IAND(IAND(g%homo(i), NOT(g%additional(i))), h%phase(i))
      twoerror = IAND(IAND(g%homo(i), g%additional(i)), NOT(h%phase(i)))

      errors(i) = IAND(allnotmissing, IOR(zeroerror, twoerror))
    end do

    do i = 64 - g%overhang, 63
      errors(g%sections) = ibclr(errors(g%sections), i)
    end do
  end function getErrorsSingle

  !---------------------------------------------------------------------------
  !> @brief   Returns the number of errors between a genotype and a
  !>		haplotype.
  !> @detail  An error is where the haplotype and the genotype are
  !>		inconsistent.  Only snps which are non-missing in both 
  !>		haplotype and the genotype are tested.
  !> @date    May 26, 2017
  !> @return  The number of errors
  !--------------------------------------------------------------------------- 
  function numberErrorsSingle(g,h) result(c)
    use HaplotypeModule
    use BitUtilities

    class(Genotype) :: g
    class(Haplotype) :: h

    integer :: c

    c = bitCount(g%getErrorsSingle(h))

  end function numberErrorsSingle

  !---------------------------------------------------------------------------
  !> @brief   Returns the location of errors between a genotype and two
  !>		haplotypes.
  !> @detail  An error is where the two haplotypes and the genotype are
  !>		inconsistent.  Only snps which are non-missing in both 
  !>		haplotypes and the genotype are tested.
  !> @date    May 26, 2017
  !> @return  The location of errors as a bit array
  !--------------------------------------------------------------------------- 
  function getErrors(g,h1,h2) result(errors)
    use HaplotypeModule

    class(Genotype) :: g
    class(Haplotype) :: h1, h2

    integer(kind=int64), dimension(:), allocatable :: errors

    integer :: i

    integer(kind=int64) :: allnotmissing, zeroerror, oneerror, twoerror

    allocate(errors(g%sections))

    do i = 1, g%sections
      allnotmissing = IAND(IOR(g%homo(i), NOT(g%additional(i))), IAND(NOT(h1%missing(i)), NOT(h2%missing(i))))

      zeroerror = IAND(IAND(g%homo(i), NOT(g%additional(i))), IOR(h1%phase(i), h2%phase(i)))
      twoerror = IAND(IAND(g%homo(i), g%additional(i)), IOR(NOT(h1%phase(i)), NOT(h2%phase(i))))
      oneerror = IAND(IAND(NOT(g%homo(i)), NOT(g%additional(i))), NOT(IEOR(h1%phase(i), h2%phase(i))))

      errors(i) = IAND(allnotmissing, IOR(  IOR(zeroerror, twoerror), oneerror))
    end do

    do i = 64 - g%overhang, 63
      errors(g%sections) = ibclr(errors(g%sections), i)
    end do

  end function getErrors

  !---------------------------------------------------------------------------
  !> @brief   Returns the number of errors between a genotype and two
  !>		haplotypes.
  !> @detail  An error is where the two haplotypes and the genotype are
  !>		inconsistent.  Only snps which are non-missing in both 
  !>		haplotypes and the genotype are tested.
  !> @date    May 26, 2017
  !> @return  The number of errors
  !--------------------------------------------------------------------------- 
  function numberErrors(g,h1,h2) result(c)
    use HaplotypeModule
    use BitUtilities

    class(Genotype) :: g
    class(Haplotype) :: h1, h2

    integer :: c

    c = bitCount(g%getErrors(h1,h2))

  end function numberErrors

    !---------------------------------------------------------------------------
    !> @brief	Returns a new genotype that is a subset of the given one
    !> @date    May 26, 2017
    !> @return  The subset genootype
    !---------------------------------------------------------------------------  
  function subset(g,start,finish) result (sub)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: start, finish
    
    type(Genotype) :: sub
    
    integer :: offset, starti, i
    
    sub%length = finish - start + 1

    sub%sections = sub%length / 64 + 1
    sub%overhang = 64 - (sub%length - (sub%sections - 1) * 64)
    
    allocate(sub%homo(sub%sections))
    allocate(sub%additional(sub%sections))
    
    starti = (start - 1) / 64  + 1
    offset = mod(start - 1, 64)
    
    if (offset == 0) then
      do i = 1, sub%sections
	sub%homo(i) = g%homo(i + starti - 1)
	sub%additional(i) = g%additional(i + starti - 1)
      end do
    else
      do i = 1, sub%sections - 1
	sub%homo(i) = IOR(ISHFT(g%homo(i + starti - 1),-offset), ISHFT(g%homo(i + starti), 64 - offset))
	sub%additional(i) = IOR(ISHFT(g%additional(i + starti - 1),-offset), ISHFT(g%additional(i + starti), 64 - offset))
      end do
      
      if (sub%sections + starti - 1 == g%sections) then
	sub%homo(sub%sections) = ISHFT(g%homo(sub%sections + starti - 1),-offset)
	sub%additional(sub%sections) = ISHFT(g%additional(sub%sections + starti - 1),-offset)
      else
	sub%homo(sub%sections) = IOR(ISHFT(g%homo(sub%sections + starti - 1),-offset), ISHFT(g%homo(sub%sections + starti), 64 - offset))
	sub%additional(sub%sections) = IOR(ISHFT(g%additional(sub%sections + starti - 1),-offset), ISHFT(g%additional(sub%sections + starti), 64 - offset))
      end if
    end if    
          
    do i = 64 - sub%overhang, 63
      sub%homo(sub%sections) = ibclr(sub%homo(sub%sections), i)
      sub%additional(sub%sections) = ibclr(sub%additional(sub%sections), i)
    end do
  end function subset
  
  !---------------------------------------------------------------------------
  !> @brief   Sets the genotype from the two haplotypes if the genotype is missing
  !> @date    May 26, 2017
  !--------------------------------------------------------------------------- 
  subroutine setFromHaplotypesIfMissing(g,h1,h2)
    use HaplotypeModule
    class(Genotype) :: g
    type(Haplotype), intent(in) :: h1, h2
    
    type(Genotype) :: gFromH
    
    gFromH = Genotype(h1,h2)
    call g%setFromOtherIfMissing(gFromH)
    
  end subroutine setFromHaplotypesIfMissing

  !---------------------------------------------------------------------------
  !> @brief   Sets the genotype from another genotype if the genotype is missing
  !> @date    May 26, 2017
  !---------------------------------------------------------------------------   
  subroutine setFromOtherIfMissing(g, o)
    class(Genotype) :: g
    type(Genotype), intent(in) :: o
    
    integer :: i
    integer(kind=int64) :: origHomo
    
    do i = 1, g%sections
    origHomo = g%homo(i)
      g%homo(i) = IOR( &
		    ! g is not missing, use g
		    IAND(IOR(g%homo(i), NOT(g%additional(i))), g%homo(i)), &
		    ! g is missing use o
		    IAND(IAND(NOT(g%homo(i)), g%additional(i)), o%homo(i)) )


      g%additional(i) = IOR( &
		    ! g is not missing, use g
		    IAND(IOR(origHomo, NOT(g%additional(i))), g%additional(i)), &
		    ! g is missing use o
		    IAND(IAND(NOT(origHomo), g%additional(i)), o%additional(i)) )
    end do
  end subroutine setFromOtherIfMissing


    !---------------------------------------------------------------------------
    !> @brief   Return a type reprsenting Mendelian consisentincies /
    !>          inconsistencies between an individual and its parents.
    !> @date    July 24, 2017
    !> @return  New mendelian object 
    !---------------------------------------------------------------------------
    function mendelianInconsistencies(g, pg, mg) result (mend)
        class(Genotype), intent(in) :: g
        type(Genotype), intent(in) :: pg, mg

        type(Mendelian) :: mend
        integer(kind=int64) :: het, patpresent, matpresent, indpresent

        integer :: i

        allocate(mend%paternalInconsistent(g%sections))
        allocate(mend%maternalInconsistent(g%sections))
        allocate(mend%individualInconsistent(g%sections))

        allocate(mend%paternalConsistent(g%sections))
        allocate(mend%maternalConsistent(g%sections))
        allocate(mend%individualConsistent(g%sections))


        do i = 1, g%sections
            het = IAND( &
                ! Individual is het
                NOT(IOR(g%homo(i),g%additional(i))), &
                ! Parents are same homo
                IAND( &
                    ! Parents are both homo
                    IAND(pg%homo(i),mg%homo(i)), &
                    ! Parents are same homo
                    NOT(IEOR(pg%additional(i), mg%additional(i))) &
                    ) &
                )

            mend%paternalInconsistent(i) = IOR( &
                ! Individual is het and bad
                het, &
                ! Individual is homo and paternal is opposite homo
                IAND(IAND(g%homo(i), pg%homo(i)), IEOR(g%additional(i), pg%additional(i))) &
                )

            mend%paternalInconsistent(i) = IOR( &
                ! Individual is het and bad
                het, &
                ! Individual is homo and maternal is opposite homo)
                IAND(IAND(g%homo(i), mg%homo(i)), IEOR(g%additional(i), mg%additional(i)))&
                )

            mend%individualInconsistent(i) = IOR(mend%paternalInconsistent(i), mend%maternalInconsistent(i))

            patpresent = IOR(pg%homo(i), NOT(pg%additional(i)))
            matpresent = IOR(mg%homo(i), NOT(mg%additional(i)))
            indpresent = IOR(g%homo(i), NOT(g%additional(i)))

            mend%paternalConsistent(i) = IAND( &
                ! Individual and paternal are both not missing
                IAND(indpresent, patpresent), &
                ! Not paternal inconsistent
                NOT(mend%paternalInconsistent(i)) &
                )

            mend%maternalConsistent(i) = IAND( &
                ! Individual and paternal are both not missing
                IAND(indpresent, matpresent), &
                ! Not paternal inconsistent
                NOT(mend%maternalInconsistent(i)) &
                )

            mend%individualConsistent(i) = IAND( &
                ! Individual and both parents are not missing
                IAND(indpresent, IAND(patpresent, matpresent)), &
                ! Not individual inconsistent
                NOT(mend%individualInconsistent(i)) &
                )
        end do

    end function mendelianInconsistencies

    !---------------------------------------------------------------------------
    !> @brief   Returns the position of the first het or -1 if none
    !> @date    August 8, 2017
    !> @return  Position of the first het 
    !---------------------------------------------------------------------------
    function firstHet(g) result(pos)
        class(Genotype), intent(in) :: g

        integer :: pos

        integer :: i, j, endj
        integer(kind=8) :: isHet

        do i = 1, g%sections
            isHet = IAND(NOT(g%homo(i)),NOT(g%additional(i)))

            if (i == g%sections) then
                endj = 63 - g%overhang
            else
                endj = 63
            end if

            do j = 0, endj
                if (BTEST(isHet,j)) then
                    pos = 64 * (i - 1) + j + 1
                    return
                end if
            end do
        end do

        pos = -1
    end function firstHet


    subroutine writeFormattedGenotype(dtv, unit, iotype, v_list, iostat, iomsg)
      class(Genotype), intent(in) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: v_list(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      integer(kind=1), dimension(:), allocatable :: array


      array = dtv%toIntegerArray()
      write(unit, "(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)", iostat = iostat, iomsg = iomsg) array
    
      end subroutine writeFormattedGenotype


      subroutine writeunFormattedGenotype(dtv, unit, iostat, iomsg)
        class(Genotype), intent(in)::dtv
        integer, intent(in):: unit
        integer, intent(out) :: iostat      ! non zero on error, etc.
        character(*), intent(inout) :: iomsg  ! define if iostat non zero.
  
        integer(kind=1), dimension(:), allocatable :: array

        array = dtv%toIntegerArray()
        write(unit, "(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)", iostat = iostat, iomsg = iomsg) array

      end subroutine writeunFormattedGenotype



      subroutine readUnformattedGenotype(dtv, unit, iostat, iomsg)
        class(genotype), intent(inout) :: dtv         ! Object to write.
        integer, intent(in) :: unit         ! Internal unit to write to.
        integer, intent(out) :: iostat      ! non zero on error, etc.
        character(*), intent(inout) :: iomsg  ! define if iostat non zero.
        integer(kind=1), dimension(:), allocatable :: array


        read(unit,"(20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1,20000i1)",iostat=iostat, iomsg=iomsg) array

        call wrapper(dtv, array)
    end subroutine readUnformattedGenotype



    pure subroutine wrapper(geno, array) 

      type(genotype), intent(out) :: geno
      integer(kind=1), dimension(:),intent(in), allocatable :: array

      geno = newGenotypeInt(array)
    
    end subroutine wrapper


    subroutine readFormattedGenotype(dtv, unit, iotype, vlist, iostat, iomsg)
      class(Genotype), intent(inout) :: dtv         ! Object to write.
      integer, intent(in) :: unit         ! Internal unit to write to.
      character(*), intent(in) :: iotype  ! LISTDIRECTED or DTxxx
      integer, intent(in) :: vlist(:)    ! parameters from fmt spec.
      integer, intent(out) :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.
      integer(kind=1), dimension(:), allocatable :: array
      read(unit,iotype, iostat=iostat, iomsg=iomsg) array

      call wrapper(dtv, array)
    end subroutine readFormattedGenotype

end module GenotypeModule


