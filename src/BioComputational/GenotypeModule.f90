module GenotypeModule

    use constantModule, only : MissingGenotypeCode
  implicit none
  private

  type, public :: Genotype
  private
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Genotype    Homo    Additional !
  ! 0           1       0          !
  ! 1           0       0          !
  ! 2           1       1          !
  ! Missing     0       1          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=8), allocatable, dimension(:), public :: homo
  integer(kind=8), allocatable, dimension(:), public :: additional
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
  procedure :: numNotMissing
  procedure :: setHaplotype
  procedure :: isZero
  procedure :: isTwo
  procedure :: isMissing
  procedure :: isHomo
  procedure :: numberErrors
  procedure :: getErrors
  procedure :: setHaplotypeIfError
  procedure :: setHaplotypeIfMissing
  procedure :: numberErrorsSingle
  procedure :: getErrorsSingle
  procedure :: subset
  end type Genotype

  interface Genotype
    module procedure newGenotypeInt
  end interface Genotype

      contains

      pure function newGenotypeInt(geno) result (g)
        integer(kind=1), dimension(:), intent(in) :: geno

        type(Genotype) :: g

        integer :: i, cursection, curpos

        g%length = size(geno,1)

        g%sections = g%length / 64 + 1
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

    function toIntegerArray(g) result(array)
        class(Genotype), intent(in) :: g

        integer(kind=1), dimension(:), allocatable :: array

        integer :: i, cursection, curpos
	
        allocate(array(g%length))

        cursection = 1
        curpos = 0
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
	  if (curpos == 64) then
	      curpos = 0
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
      phase(i) = IOR( IAND(h%phase(i), h%missing(i)), &
!      phase(i) = IOR( h%missing(i), &
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

function numNotMissing(g) result(c)
    class(Genotype), intent(in) :: g

    integer :: c

    integer :: i

    c = 0
    do i = 1, g%sections
      c = c + POPCNT(NOT(IAND(NOT(g%homo(i)), g%additional(i))))
  end do
end function numNotMissing


subroutine setHaplotype(g, h)
    use HaplotypeModule

    class(Haplotype), intent(in) :: h
    class(Genotype), intent(in) :: g

    integer :: i

    do i = 1, g%sections
      h%missing(i) = NOT(g%homo(i))

      h%phase(i) = IAND(g%homo(i), g%additional(i))
  end do
  do i = 64 - g%overhang + 1, 64
      h%missing(g%sections) = ibclr(h%missing(g%sections), i - 1)
    end do
  end subroutine setHaplotype

  subroutine setHaplotypeIfError(g, h, errors)
    use HaplotypeModule

    type(Haplotype), intent(in), pointer :: h
    class(Genotype), intent(in) :: g
    integer(kind=8), intent(in), dimension(:), pointer :: errors

    integer :: i

    do i = 1, g%sections
      h%missing(i) = IOR(IAND(NOT(errors(i)), h%missing(i)), IAND(errors(i), NOT(g%homo(i))))

      h%phase(i) = IOR(IAND(NOT(errors(i)), h%phase(i)), IAND(errors(i), IAND(g%homo(i), g%additional(i))))
    end do
    do i = 64 - g%overhang + 1, 64
      h%missing(g%sections) = ibclr(h%missing(g%sections), i - 1)
    end do
  end subroutine setHaplotypeIfError

  subroutine setHaplotypeIfMissing(g, h)
    use HaplotypeModule

    type(Haplotype), intent(in) :: h
    class(Genotype), intent(in) :: g

    integer :: i

    do i = 1, g%sections
      h%phase(i) = IOR(IAND(NOT(h%missing(i)), h%phase(i)), IAND(h%missing(i), IAND(g%homo(i), g%additional(i))))
      h%missing(i) = IAND(h%missing(i), NOT(g%homo(i)))
    end do
    do i = 64 - g%overhang + 1, 64
      h%missing(g%sections) = ibclr(h%missing(g%sections), i)
    end do
  end subroutine setHaplotypeIfMissing

  function isZero(g, pos) result (zero)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: zero

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    zero = BTEST(IAND(g%homo(cursection),NOT(g%additional(cursection))), curpos)
end function isZero

function isTwo(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(IAND(g%homo(cursection),g%additional(cursection)), curpos)
end function isTwo

function isMissing(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(IAND(NOT(g%homo(cursection)),g%additional(cursection)), curpos)
end function isMissing

function isHomo(g, pos) result (two)
    class(Genotype), intent(in) :: g
    integer, intent(in) :: pos

    logical :: two

    integer :: cursection, curpos

    cursection = (pos-1) / 64 + 1
    curpos = pos - (cursection - 1) * 64 - 1


    two = BTEST(g%homo(cursection), curpos)
  end function isHomo

  function getErrorsSingle(g,h) result(errors)
    use HaplotypeModule

    class(Genotype) :: g
    class(Haplotype) :: h

    integer(kind=8), dimension(:), pointer :: errors

    integer :: i

    integer(kind=8) :: allnotmissing, zeroerror, twoerror

    allocate(errors(g%sections))

    do i = 1, g%sections
      allnotmissing = IAND(IOR(g%homo(i), NOT(g%additional(i))), NOT(h%missing(i)))

      zeroerror = IAND(IAND(g%homo(i), NOT(g%additional(i))), h%phase(i))
      twoerror = IAND(IAND(g%homo(i), g%additional(i)), NOT(h%phase(i)))

      errors(i) = IAND(allnotmissing, IOR(zeroerror, twoerror))
    end do

    do i = 64 - g%overhang + 1, 64
      errors(g%sections) = ibclr(errors(g%sections), i - 1)
    end do
  end function getErrorsSingle

  function numberErrorsSingle(g,h) result(c)
    use HaplotypeModule
    use BitUtilities

    class(Genotype) :: g
    class(Haplotype) :: h

    integer :: c

    c = bitCount(g%getErrorsSingle(h))

  end function numberErrorsSingle


  function getErrors(g,h1,h2) result(errors)
    use HaplotypeModule

    class(Genotype) :: g
    class(Haplotype) :: h1, h2

    integer(kind=8), dimension(:), pointer :: errors

    integer :: i

    integer(kind=8) :: allnotmissing, zeroerror, oneerror, twoerror

    allocate(errors(g%sections))

    do i = 1, g%sections
      allnotmissing = IAND(IOR(g%homo(i), NOT(g%additional(i))), IAND(NOT(h1%missing(i)), NOT(h2%missing(i))))

      zeroerror = IAND(IAND(g%homo(i), NOT(g%additional(i))), IOR(h1%phase(i), h2%phase(i)))
      twoerror = IAND(IAND(g%homo(i), g%additional(i)), IOR(NOT(h1%phase(i)), NOT(h2%phase(i))))
      oneerror = IAND(IAND(NOT(g%homo(i)), NOT(g%additional(i))), NOT(IEOR(h1%phase(i), h2%phase(i))))

      errors(i) = IAND(allnotmissing, IOR(  IOR(zeroerror, twoerror), oneerror))
    end do

    do i = 64 - g%overhang + 1, 64
      errors(g%sections) = ibclr(errors(g%sections), i - 1)
    end do

  end function getErrors

  function numberErrors(g,h1,h2) result(c)
    use HaplotypeModule
    use BitUtilities

    class(Genotype) :: g
    class(Haplotype) :: h1, h2

    integer :: c

    c = bitCount(g%getErrors(h1,h2))

  end function numberErrors
  
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
	sub%homo(sub%sections) = ISHFT(g%homo(sub%sections + starti - 1),offset)
	sub%additional(sub%sections) = ISHFT(g%additional(sub%sections + starti - 1),offset)
      else
	sub%homo(sub%sections) = IOR(ISHFT(g%homo(sub%sections + starti - 1),-offset), ISHFT(g%homo(sub%sections + starti), 64 - offset))
	sub%additional(sub%sections) = IOR(ISHFT(g%additional(sub%sections + starti - 1),-offset), ISHFT(g%additional(sub%sections + starti), 64 - offset))
      end if
    end if    
          
    do i = 64 - sub%overhang + 1, 64
      sub%homo(sub%sections) = ibclr(sub%homo(sub%sections), i)
      sub%additional(sub%sections) = ibclr(sub%additional(sub%sections), i)
    end do
  end function subset

end module GenotypeModule


