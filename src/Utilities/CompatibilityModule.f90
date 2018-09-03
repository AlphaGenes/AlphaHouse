#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"


#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""


#endif
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     CompatibilityModule.f90
!
! DESCRIPTION:
!> @brief    Module cotaining subroutines to deal with text PLINK format
!> @details  currently only contains integer and real heap sort procedures
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     January 4, 2017
!
!> @version  1.0.0
!
!
!-------------------------------------------------------------------------------
module CompatibilityModule
	use integerLinkedListModule
	use HashModule
	use PedigreeModule

	implicit none

	type :: bimHolder
	character(len=2) :: ref,alt
	character(len=IDLENGTH) :: id
	character(len=2) :: chrom !<either an integer, or 'X'/'Y'/'XY'/'MT'
	integer(kind=int64) :: pos, chrompos

end type bimHolder

type Chromosome

type(integerLinkedList) :: snps
end type Chromosome


type plinkInfoType

type(chromosome), dimension(:), allocatable :: chromosomes

integer :: nChroms !< total number of chromosomes
integer :: totalSnps !< total numbers of snps across all chromosomes
character(len=2), dimension(:),allocatable :: referenceAllelePerSnps, alternateAllelePerSnps
integer(kind=1), dimension(:,:), allocatable ::  genotypes	 !<  genotypes for all chroms (nanimals, totalsnps)
integer(kind=1), dimension(:,:,:), allocatable ::  phase
type(DictStructure ) :: dict
real(real64), dimension(:,:) ,allocatable :: lengths !< array of snp lengths in morgans, in format (chrom, snps)
integer, dimension(:,:) ,allocatable :: basepairs !< array of basepair positions, in format (chrom, snps)
integer :: sexChrom = 0 !<if not 0, then this is the sex chrom
integer, dimension(:), allocatable :: nsnpsPerChromosome !< nsnps per chromosome
character(len=128), dimension(:), allocatable :: snpName

contains
procedure :: initPlinkInfoType
final :: destroyPlinkInfoType

end type plinkInfoType
contains


function createBimInfoFromGenotypes(genotypes) result(bimOut)
	type(bimHolder),allocatable, dimension(:) :: bimOut
	integer(kind=1), dimension(:,:) :: genotypes
	integer :: nsnps,i,nanimals
	character(len=16) :: snpNumber 

	nsnps = size(genotypes,2)
	nanimals = size(genotypes,1)
	allocate(bimOut(nsnps))
	
	do i=1,nsnps
		write (snpNumber, '(a,I13.13)') "SNP",i
		bimOut(i)%id = snpNumber
		bimOut(i)%chrom = "1"
		bimOut(i)%pos = 0
		bimOut(i)%chromPos = 1
		bimOut(i)%ref = "1"
		bimOut(i)%alt = "2"
	enddo

end function createBimInfoFromGenotypes


subroutine initPlinkInfoType(plinkInfo, ped)

	class(plinkInfoType) :: plinkInfo

	type(PedigreeHolder) :: ped
	integer :: i


	plinkInfo%nChroms = 1
	plinkInfo%totalSnps = ped%nsnpsPopulation
	allocate(plinkInfo%referenceAllelePerSnps(ped%nsnpsPopulation))
	allocate(plinkInfo%alternateAllelePerSnps(ped%nsnpsPopulation))
	plinkInfo%referenceAllelePerSnps = "1"
	plinkInfo%alternateAllelePerSnps = "2"
	plinkInfo%genotypes = ped%getGenotypesAsArrayWitHMissing()
	plinkInfo%phase = ped%getPhaseAsArrayWithMissing()
	allocate(plinkInfo%lengths(1,ped%nsnpsPopulation))
	plinkInfo%lengths = 0
	allocate(plinkInfo%basepairs(1,ped%nsnpsPopulation))
	plinkInfo%basepairs = 0
	allocate(plinkInfo%nsnpsPerChromosome(1))
	plinkInfo%nsnpsPerChromosome = ped%nsnpsPopulation
	allocate(plinkInfo%snpName(ped%nsnpsPopulation))
	do i=1,ped%nsnpsPopulation
		write(plinkInfo%snpName(i),'(a,i2)') "SNP",i
	enddo
end subroutine initPlinkInfoType

!---------------------------------------------------------------------------
!< @brief destructor for plinkInfoType
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine destroyPlinkInfoType(plinkInfo)
	type(plinkInfoType),intent(inout)  :: plinkInfo

	if (allocated(plinkInfo%referenceAllelePerSnps)) deallocate(plinkInfo%referenceAllelePerSnps)
	if (allocated(plinkInfo%alternateAllelePerSnps)) deallocate(plinkInfo%alternateAllelePerSnps)
	if (allocated(plinkInfo%genotypes)) deallocate(plinkInfo%genotypes)
	if (allocated(plinkInfo%phase)) deallocate(plinkInfo%phase)
	if (allocated(plinkInfo%lengths)) deallocate(plinkInfo%lengths)
	if (allocated(plinkInfo%basepairs)) deallocate(plinkInfo%basepairs)
	if (allocated(plinkInfo%nsnpsPerChromosome	)) deallocate(plinkInfo%nsnpsPerChromosome)
	if (allocated(plinkInfo%snpName)) deallocate(plinkInfo%snpName)

end subroutine destroyPlinkInfoType



!---------------------------------------------------------------------------
!< @brief reads fam file to pedigree object
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readFamFile(ped,pedFile)
	use ConstantModule, only : IDLENGTH
	use AlphaHouseMod, only : countLines
	use PedigreeModule

	type(pedigreeHolder), intent(out) :: ped !< Pedigree object that is returned
	character(len=*), intent(in) :: pedFile !< .ped file generated by plink
	character(len=IDLENGTH) :: familyID,gender,phenotype
	integer :: fileUnit,stat, i,lines
	character(len=IDLENGTH),dimension(:,:), allocatable :: pedArray
	integer, allocatable, dimension(:) :: genderArray, phenotypeArray

	lines=  countLines(pedFile)

	allocate(pedArray(3,lines))
	allocate(genderArray(lines))
	allocate(phenotypeArray(lines))

	open(newUnit=fileUnit, file=pedFile, status="old")
	print *,"TOTAL ANS",lines
	do i=1, lines
		read(fileUnit,*) familyID,pedArray(1,i),pedArray(2,i),pedArray(3,i),gender,phenotype
		read(gender,*,iostat=stat)  genderArray(i)
		read(phenotype,*,iostat=stat)  phenotypeArray(i)
	enddo
	call initPedigreeArrays(ped,pedArray, genderArray)

	print *, "Animals in fam file:",ped%pedigreeSize," without dummies:",ped%addedRealAnimals

	call ped%printPedigreeOriginalFormat("pedigreeOutput.txt")

end subroutine readFamFile

!---------------------------------------------------------------------------
!< @brief Reads PLINK binary format into ped datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readPlink(binaryFilePre, ped, outputPaths,plinkInfo, useChroms)
	use HashModule
	use PedigreeModule

	character(len=*),intent(in) :: binaryFilePre !< part before file extension
	type(pedigreeholder), intent(out) :: ped !< pedigree object returned
	character(len=128), dimension(:), allocatable, intent(out) :: outputPaths !< output paths for each chromosome
	type(bimHolder) , allocatable, dimension(:) :: bimInfo
	type(plinkInfoType) :: plinkInfo
	integer, dimension(:), intent(in),allocatable, optional :: useChroms

	! TODO change to pointer rather than copy
	call readFamFile(ped,trim(binaryFilePre)//".fam")
	call readBim(trim(binaryFilePre)//".bim",bimInfo,plinkInfo)
	print *,"READ BIM"
	call readBED(trim(binaryFilePre)//".bed",ped,4,plinkInfo)
	print *,"READ BED"

	if (present(useChroms)) then
		if(allocated(useChroms)) then
			call createOutputFiles(ped, outputPaths,plinkInfo,useChroms)
			return
		endif
	endif

	call createOutputFiles(ped, outputPaths,plinkInfo)


end subroutine readPlink

!---------------------------------------------------------------------------
!< @brief writes output files
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine createOutputFiles(ped, outputPaths,plinkInfo, useChroms)
	use ifport
	use ALphaHouseMod, only : countLines

	type(pedigreeholder), intent(in) :: ped !< pedigree object taken in
	type(plinkInfoType), intent(in) :: plinkInfo
	character(len=128), dimension(:), allocatable,intent(out) :: outputPaths !< returns the output path for each chromosome

	integer :: outChrF, outChrP,bpFile,lengthFile
	character(len=128) :: path, outChrFile,fmt
	logical, dimension(:), allocatable :: maskedLogi
	integer, dimension(:), allocatable :: masked
	integer :: result, i, h,p,refAlleleUnit,count
	integer, dimension(:), intent(in), optional :: useChroms
	path = "fullGenome/"
	result=makedirqq(path)

	print *, "Total number of Chromosomes:",plinkInfo%nChroms

	allocate(outputPaths(plinkInfo%nChroms))
	allocate(maskedLogi(size(plinkInfo%genotypes,2))) !< alloc to total number of snps
	do i =1, plinkInfo%nChroms


		if (present(useChroms)) then
			if (.not. any(useChroms == i)) cycle
		endif
		write(outChrFile, '(a,a,i0.2,a)') trim(path),trim("chr"),i,DASH

		result=makedirqq(outChrFile)
		outputPaths(i) = outChrFile

		open(newunit=outChrF, file=trim(outChrFile)//"genotypes.txt", status="unknown")
		open(newunit=outChrP, file=trim(outChrFile)//"phase.txt", status="unknown")
		open(newunit=lengthFile, file=trim(outChrFile)//"snplengths.txt", status="unknown")
		open(newunit=bpFile, file=trim(outChrFile)//"snpBasepairs.txt", status="unknown")

		!< uses genotype mask to determine which snps to use
		masked = plinkInfo%chromosomes(i)%snps%convertToArray()
		maskedLogi = .false.
		do h =1, size(masked)
			maskedLogi(masked(h)) = .true.
		enddo


		open(newunit=refAlleleUnit, file=trim(outChrFile)//"refAlleles", status="unknown")

		! alleles = pack(referenceAllelePerSnps, maskedLogi)

		write(refAlleleUnit, '(3a)') "snp in Chrom", "snp pos","Ref Allele"
		count = 0
		do p=1, size(plinkInfo%referenceAllelePerSnps)
			! do p=1,size(alleles)
			! write(refAlleleUnit, '(1i5,a5)') p, referenceAllelePerSnps(p)
			if (maskedLogi(p)) then
				count = count +1
				write(refAlleleUnit, '(2i5,a5)') count,p, plinkInfo%referenceAllelePerSnps(p)
			endif
		enddo
		close(refAlleleUnit)





		! set up the pedigree to avoid read in
		! array = 9
		! do p=1,ped%addedRealAnimals
		! 	array(p,:) = pack(plinkInfo%genotypes(p,:), maskedLogi)
		! end do

		write(fmt, '(a,i10,a)')   "(",plinkInfo%nsnpsPerChromosome(i), "F7.2)"
		write(lengthFile, fmt) pack(plinkInfo%lengths(i,:), maskedLogi)
		write(fmt, '(a,i10,a)')   "(",plinkInfo%nsnpsPerChromosome(i), "i10)"
		write(bpFile, fmt) pack(plinkInfo%basepairs(i,:), maskedLogi)


		write(fmt, '(a,i10,a)')  "(a20,", plinkInfo%nsnpsPerChromosome(i), "i3)"

		do p=1,ped%addedRealAnimals
			write(outChrF,fmt) ped%pedigree(p)%originalId,pack(plinkInfo%genotypes(p,:), maskedLogi)
			write(outChrp,fmt) ped%pedigree(p)%originalId,pack(plinkInfo%phase(p,:,1), maskedLogi)
			write(outChrp,fmt) ped%pedigree(p)%originalId,pack(plinkInfo%phase(p,:,2), maskedLogi)
		end do
		print *, "Finished writeout of chromosome ", i

		close(outChrF)
		close(outChrP)
		close(lengthFile)
		close(bpFile)
	enddo


end subroutine createOutputFiles



!---------------------------------------------------------------------------
!< @brief Reads bim file to datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readBim(bimFile, bimInfo,plinkInfo)
	use HashModule
	use AlphaHouseMod
	use ConstantModule

	character(len=*), intent(in) :: bimFile
	type(bimHolder) , allocatable, dimension(:), intent(out) :: bimInfo !< extra info provided by BIM file
	type(plinkInfoType), intent(out) :: plinkInfo
	character :: ref,alt
	character(len=IDLENGTH) :: id

	integer(kind=int64) :: pos
	real(kind=real64) :: chrompos

	real(kind=real64), dimension(:,:) ,allocatable :: tmpLengths
	integer, dimension(:,:) ,allocatable :: tmpbasePairs
	integer, dimension(:), allocatable :: temparray

	integer :: i, unit,chromCount,curChromSnpCount
	character(len=2) :: chrom,prevChrom

	plinkInfo%nChroms = 0
	curChromSnpCount = 0
	allocate(plinkInfo%nsnpsPerChromosome(LARGECHROMNUMBER))
	plinkInfo%nsnpsPerChromosome = 0
	call plinkInfo%dict%DictStructure()
	plinkInfo%totalSnps = countLines(bimFile)

	open(newUnit=unit, file=bimFile, status='old')
	allocate(plinkInfo%chromosomes(LARGECHROMNUMBER))
	allocate(plinkInfo%snpName(plinkInfo%totalSnps))
	allocate(bimInfo(plinkInfo%totalSnps))
	allocate(plinkInfo%lengths(LARGECHROMNUMBER,plinkInfo%totalSnps))
	allocate(plinkInfo%basepairs(LARGECHROMNUMBER,plinkInfo%totalSnps))
	chromCount = 1
	plinkInfo%lengths = 0
	plinkInfo%basepairs = 0
	do i =1, plinkInfo%totalSnps

		read(unit, *) chrom, id,chrompos, pos ,ref, alt

		if (i == 1) then
			prevChrom = chrom
		endif

		! if we've moved on to the next chromsome


		! if we 've moved on to the next chromosome

		if (chrom /=prevChrom) then

			plinkInfo%nsnpsPerChromosome(chromCount) = curChromSnpCount
			curChromSnpCount = 0
			chromCount = chromCount + 1
			prevChrom = chrom

			if (chrom == "X" .or. chrom == 'Y') then
				plinkInfo%sexChrom = chromCount
			endif
			if (chrom == 'XY' .or. chrom == 'MT') then
				write(error_unit,*) "WARNING - No support currently for XY or MT chromosomes"
			endif

			if (chromCount > plinkInfo%nChroms) then
				plinkInfo%nChroms = chromCount
			endif
		endif
		! last chrom
		if (i == plinkInfo%totalSnps) then
			plinkInfo%nsnpsPerChromosome(chromCount) = curChromSnpCount + 1
			if (chromCount > plinkInfo%nChroms) then
				plinkInfo%nChroms = chromCount
			endif
		end if

		curChromSnpCount = curChromSnpCount + 1

		plinkInfo%lengths(chromCount,curChromSnpCount) = chromPos
		plinkInfo%basepairs(chromCount,curChromSnpCount) = pos

		call plinkInfo%dict%addKey(id, i)
		plinkInfo%snpName(i) = id
		bimInfo(i)%chrom = chrom
		!  TODO clean this up

		bimInfo(i)%id = id


		bimInfo(i)%chrompos = chrompos
		bimInfo(i)%pos = pos
		bimInfo(i)%ref = ref
		bimInfo(i)%alt = alt

		call plinkInfo%chromosomes(chromCount)%snps%list_add(i)
	end do

	if (chromCount /= LARGECHROMNUMBER) then
		allocate(temparray(chromCount))
		temparray(1:chromCount) = plinkInfo%nsnpsPerChromosome(1:chromCount)
		call move_alloc(temparray,plinkInfo%nsnpsPerChromosome)

		allocate(tmpLengths(chromCount,plinkInfo%totalSnps))
		tmpLengths(1:chromCount,:) = plinkInfo%lengths(1:chromCount,:)
		call move_alloc(tmpLengths,plinkInfo%lengths)


		allocate(tmpbasePairs(chromCount,plinkInfo%totalSnps))
		tmpbasePairs(1:chromCount,:) = plinkInfo%basepairs(1:chromCount,:)
		call move_alloc(tmpbasePairs,plinkInfo%basepairs)
	endif
	close (unit)


end subroutine readBim

!---------------------------------------------------------------------------
!< @brief Reads BED Files
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine readBED(bed, ped, minor, plinkInfo)

	use PedigreeModule
	use genotypeModule
	implicit none

	! Arguments
	character(*), intent(in) :: bed !< bed file name
	type(PedigreeHolder), intent(in) :: ped !< pedigree read in from .fam
	integer, intent(in) ::  minor !< if the first allele is minor or major
	type(plinkInfoType), intent(inout) :: plinkInfo
	integer :: status

	!! Types
	INTEGER, PARAMETER :: Byte = SELECTED_INT_KIND(1) ! Byte

	!! Local arguments
	integer(Byte) :: readplinkmode, element, plinkmode
	integer(Byte), dimension(2) :: readmagicnumber, magicnumber
	!logical :: checkmaf
	integer :: stat, i, j, k, snpcount, majorcount
	integer, dimension(4) :: codes, phasecodes
	!integer, dimension(:), allocatable :: domasksnps
	real :: allelefreq
	integer :: bedInUnit
	! Supported formats as per plink 1.9.
	!data magicnumber/X"6C",X'0000001B' /,  plinkmode/X'01'/
	data magicnumber /108,27/, plinkmode /1/


	print *,"start BED read"
	allocate(plinkInfo%genotypes(ped%addedRealAnimals,plinkInfo%totalSnps))
	allocate(plinkInfo%phase(ped%addedRealAnimals,plinkInfo%totalSnps,2))
	plinkInfo%genotypes(:,:) = MISSINGGENOTYPECODE
	plinkInfo%phase(:,:,:) = MISSINGPHASECODE

	if (minor == 1) then
		codes = (/ 0, 1, 2, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 0, 1, 1, MISSINGGENOTYPECODE /)
	else
		codes = (/ 2, 1, 0, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 1, 1, 0, MISSINGGENOTYPECODE /)
	endif

	open(newunit=bedInUnit, file=bed, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED')
	print *, "start read"
	read(bedInUnit) readmagicnumber, readplinkmode
	if (all(readmagicnumber /= magicnumber) ) then
		status=-1
		close(bedInUnit)
		return
	endif
	if (readplinkmode /= plinkmode) then
		status=-2
		close(bedInUnit)
		return
	endif


	j=0  ! Sample-index
	k=1  ! SNP-index
	snpcount = 0
	majorcount = 0
	outer: do


		read(bedInUnit, iostat=stat) element
		if (stat /= 0) exit
		inner: do i=0,6,2
			j = j + 1
			snpcount = snpcount + 1
			select case(IBITS(element, i, 2))
				case (0) ! homozygote
				plinkInfo%genotypes(j,k) = codes(1)
				plinkInfo%phase(j,k,:) = phasecodes(1)
				case (1) ! missing
				plinkInfo%genotypes(j,k) = codes(4)
				plinkInfo%phase(j,k,:) = phaseCodes(4)
				snpcount = snpcount - 1
				case (2) ! heterozygote
				plinkInfo%genotypes(j,k) = codes(2)
				! Set to missing - as we don't know which snp is which
				plinkInfo%phase(j,k,1) = MISSINGPHASECODE
				plinkInfo%phase(j,k,2) = MISSINGPHASECODE
				majorcount = majorcount + 1
				case (3) ! homozygote, minor
				plinkInfo%genotypes(j,k) = codes(3)
				plinkInfo%phase(j,k,:) = phasecodes(3)
				majorcount = majorcount + 2
			endselect
			if (j == ped%addedRealAnimals) then
				if (snpcount /= 0) then
					allelefreq = majorcount / (snpcount*2.)
				endif
				j = 0
				snpcount = 0
				majorcount = 0
				k = k + 1
				cycle outer
			endif
		enddo inner
		! print *, "in outer"
	enddo outer
	close(bedInUnit)

	if (stat == -1) stat=0

	print *, "finished"

end subroutine readBED

!---------------------------------------------------------------------------
!< @brief Reads PLINK non binary version .map, .ref(optional) and .ped files
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine readPlinkNoneBinary(filePre,ped,outputPaths,plinkInfo,useChroms)
	use HashModule
	use PedigreeModule
	use AlphaHouseMod
	use ifport

	character(len=*),intent(in) :: filePre
	character(len=128), dimension(:), allocatable,intent(out) :: outputPaths
	integer, dimension(:), intent(in),allocatable, optional :: useChroms
	type(plinkInfoType) :: plinkInfo
	type(pedigreeHolder) :: ped

	call readMap(trim(filePre)//".map",plinkInfo)
	call readRef(trim(filePre)//".ref", plinkInfo)
	call readPedFile(trim(filePre)//".ped",ped, plinkInfo)
	if (present(useChroms)) then
		if(allocated(useChroms)) then
			call createOutputFiles(ped, outputPaths,plinkInfo,useChroms)
			return
		endif
	endif

	call createOutputFiles(ped, outputPaths,plinkInfo)

end subroutine readPlinkNoneBinary


!---------------------------------------------------------------------------
!< @brief Reads plink map file into datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readMap(filename,plinkInfo)
	use HashModule
	use AlphaHouseMod

	character(len=*),intent(in) :: filename

	type(plinkInfoType), intent(out) :: plinkInfo

	integer, dimension(:) ,allocatable :: tempArray
	real(real64), dimension(:,:) ,allocatable ::tmpLengths
	integer, dimension(:,:) ,allocatable :: tmpbasePairs
	integer :: unit,i,chromCount,basepair

	real(kind=real64) :: length
	character(len=2) :: chrom,prevChrom
	character(len=128) :: id

	plinkInfo%totalSnps = countLines(fileName)
	print *,"TOTAL SNPS IN MAP:", plinkInfo%totalSnps
	open(newunit=unit, file=filename, status='OLD')

	allocate(plinkInfo%chromosomes(LARGECHROMNUMBER))
	allocate(plinkInfo%lengths(LARGECHROMNUMBER, plinkInfo%totalSnps))
	allocate(plinkInfo%basepairs(LARGECHROMNUMBER, plinkInfo%totalSnps))
	allocate(plinkInfo%snpName(plinkInfo%totalSnps))
	allocate(plinkInfo%nsnpsPerChromosome(LARGECHROMNUMBER))

	call plinkInfo%dict%DictStructure()
	plinkInfo%nChroms = 0
	plinkInfo%nsnpsPerChromosome = 0
	plinkInfo%lengths = 0
	plinkInfo%basepairs = 0
	chromCount = 0
	prevChrom = 'MT'

	write(*,*) "Start reading map file"
	do i=1,plinkInfo%totalSnps

		read(unit, *) chrom, id,length, basepair

		if (chrom /=prevChrom) then
			chromCount = chromCount + 1
			prevChrom = chrom
			if (chrom == "X" .or. chrom == 'Y') then
				plinkInfo%sexChrom = chromCount
			endif
			if (chrom == 'XY' .or. chrom == 'MT') then
				write(error_unit,*) "WARNING - No support currently for XY or MT chromosomes"
			endif
		endif

		plinkInfo%nsnpsPerChromosome(chromCount) = plinkInfo%nsnpsPerChromosome(chromCount) + 1
		plinkInfo%lengths(chromCount,plinkInfo%nsnpsPerChromosome(chromCount)) = length
		plinkInfo%basepairs(chromCount,plinkInfo%nsnpsPerChromosome(chromCount) ) = basepair
		call plinkInfo%dict%addKey(id, i)
		plinkInfo%snpName(i) = id
		call plinkInfo%chromosomes(chromCount)%snps%list_add(i)
	enddo

	if (chromCount /= LARGECHROMNUMBER) then
		allocate(temparray(chromCount))
		temparray(1:chromCount) = plinkInfo%nsnpsPerChromosome(1:chromCount)
		call move_alloc(temparray,plinkInfo%nsnpsPerChromosome)

		allocate(tmpLengths(chromCount,plinkInfo%totalSnps))
		tmpLengths(1:chromCount,:) = plinkInfo%lengths(1:chromCount,:)
		call move_alloc(tmpLengths,plinkInfo%lengths)


		allocate(tmpbasePairs(chromCount,plinkInfo%totalSnps))
		tmpbasePairs(1:chromCount,:) = plinkInfo%basepairs(1:chromCount,:)
		call move_alloc(tmpbasePairs,plinkInfo%basepairs)
	endif

	plinkInfo%nChroms = chromCount
	write(*,*) "Finished reading map file"
end subroutine readMap


!---------------------------------------------------------------------------
!< @brief Sets the individual to be genotyped at high density.
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readPedFile(filename,ped, plinkInfo)
	use PedigreeModule
	use ConstantModule
	use AlphaHouseMod

	character(len=*), intent(in) :: filename
	type(pedigreeHolder), intent(out) :: ped

	character(len=IDLENGTH), dimension(:,:), allocatable :: pedArray
	type(plinkInfoType), intent(inout) :: plinkInfo
	character(len=2),dimension(:,:), allocatable :: alleles !<(nanimals. nsnp*2) size nsnp x2 (for each allele,)
	integer, allocatable, dimension(:) :: genderArray, phenotypeArray
	integer :: size,cursnp,i,j,stat, fileUnit,gender,phenotype
	character(len=1) :: all1, all2
	character(len=IDLENGTH),dimension(:),allocatable :: FAMILYID
	integer, dimension(4) :: codes,phaseCodes
	integer :: minor

	size = countlines(filename)

	allocate(familyID(size))
	allocate(pedArray(3,size))
	allocate(genderArray(size))
	allocate(phenotypeArray(size))



	allocate(alleles(size, plinkInfo%totalSnps*2))
	allocate(plinkInfo%phase(size,plinkInfo%totalSnps,2))
	allocate(plinkInfo%genotypes(size,plinkInfo%totalSnps))

	plinkInfo%genotypes(:,:) = MISSINGGENOTYPECODE
	plinkInfo%phase(:,:,:) = MISSINGPHASECODE
	open(newunit=fileUnit, file=fileName, status="unknown")

	! minor is always false atm - this is because we always want the same ref allele to be used here - potentially used for future
	minor = 4
	if (minor == 1) then
		codes = (/ 0, 1, 2, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 0, 1, 1, MISSINGPHASECODE /)
	else
		codes = (/ 2, 1, 0, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 1, 1, 0, MISSINGPHASECODE /)
	endif
	write(*,*) "Start Reading Ped File, number of ans:",size
	do i=1,size
		read(fileUnit,*) familyID(i),pedArray(1,i),pedArray(2,i),pedArray(3,i),genderArray(i),phenotypeArray(i), alleles(i,:)
		print *, "finished formatting for animal:",i
	enddo

	close(fileUnit)


	write(*,*) "Finished reading - now processing"

	! check if reference alleles have been passed in
	if (.not. allocated(plinkInfo%referenceAllelePerSnps)) then
		allocate(plinkInfo%referenceAllelePerSnps(plinkInfo%totalSnps))
		allocate(plinkInfo%alternateAllelePerSnps(plinkInfo%totalSnps))

		do j=1,plinkInfo%totalSnps*2,2
			cursnp = (j/2) + 1
			block
				character(len=2) :: one,two
				integer :: onec,twoc
				one = getFirstNonMissing(alleles, j)
				onec = 0
				twoc = 0
				do i=1,size
					all1 = alleles(i,j)
					all2 = alleles(i,j+1)
					! check for first allele
					if (all1 == one ) then
						onec = onec + 1
					else
						if (all2 == '0') cycle
						if (twoc == 0) then
							two = all1
						endif

						twoc = twoc + 1
					endif
					! check for second allele
					if (all2 == one ) then
						onec = onec + 1
					else
						if (all2 == '0') cycle
						if (twoc == 0) then
							two = all1
						endif
						twoc = twoc + 1
					endif
				enddo
				if (onec < twoc) then
					plinkInfo%referenceAllelePerSnps(cursnp) = one
					plinkInfo%alternateAllelePerSnps(cursnp) = two
				else
					plinkInfo%referenceAllelePerSnps(cursnp) = two
					plinkInfo%alternateAllelePerSnps(cursnp) = one
				endif
			end block
		enddo
	endif

	do j=1,plinkInfo%totalSnps*2,2
		cursnp = (j/2) + 1
		! ref is always the first - allele in snp -
		! referenceAllelePerSnps(cursnp) = alleles(i,j+1)
		! if (.not. allocated(referenceAllelePerSnps)) then
		! 	referenceAllelePerSnps(cursnp) = getFirstHetPosition(alleles, cursnp)
		! endif
		do i=1,size !< loop through animals
			all1 = alleles(i,j)
			all2 = alleles(i,j+1)
			if (all1 == '0' .or. all2 == '0') then
				plinkInfo%genotypes(i,curSnp) = MISSINGGENOTYPECODE

			else if (all1 == all2) then

				if (all1 == plinkInfo%referenceAllelePerSnps(curSnp) ) then
					plinkInfo%genotypes(i,curSnp) = codes(1)
					plinkInfo%phase(i,cursnp,:) = phasecodes(1)
				else
					plinkInfo%genotypes(i,curSnp) = codes(3)
					plinkInfo%phase(i,cursnp,:) = phasecodes(3)

				endif
			else !< means they are different
				plinkInfo%genotypes(i,curSnp) = codes(2)

				if (all1 == plinkInfo%referenceAllelePerSnps(cursnp)) then
					plinkInfo%phase(i,cursnp,1) = phasecodes(3)
					plinkInfo%phase(i,cursnp,2) =  phasecodes(1)
				else
					plinkInfo%phase(i,cursnp,1) = phasecodes(1)
					plinkInfo%phase(i,cursnp,2) = phasecodes(3)
				endif

			endif
		enddo
	enddo
	write(*,*) "Finished Reading Ped File"
	call initPedigreeArrays(ped,pedArray, genderArray)

end subroutine readPedFile


!---------------------------------------------------------------------------
!< @brief Returns the first non missing character in allele array
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
character function getFirstNonMissing(alleles, snp)
	character(len=2),dimension(:,:), allocatable, intent(in) :: alleles !< alleles, in format (nanimals, nsnps)
	integer, intent(in) :: snp
	integer :: i, h
	do i=1, size(alleles,1) !< loops through number of animals
		do h=snp,snp+1
			if (alleles(i,h) /= '0') then
				getFirstNonMissing = alleles(i,h)
				return
			endif
		enddo
	enddo

	write(error_unit,*) "WARNING - no reference alleles for snp!"
	getFirstNonMissing = '0'
end function getFirstNonMissing




!---------------------------------------------------------------------------
!< @brief returns the first snp that is a het at position pos
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
character function getFirstHetPosition(alleles, snp)
	character(len=2),dimension(:,:), allocatable, intent(in) :: alleles
	character(len=2) :: all1, all2
	integer, intent(in) :: snp !< snp position
	integer :: i

	do i=1, size(alleles,1) !< loops through number of animals
		all1 = alleles(i,snp)
		all2 = alleles(i,snp+1)

		! adds support for number format
		if (all1 == "1" .or. all2 == "1") then
			getFirstHetPosition = "1"
			return
		endif
		if (all1 == all2) then
			cycle
		else
			getFirstHetPosition = all1
			return
		endif
	enddo

	write(error_unit,*) "WARNING - no reference alleles for snp!"
	getFirstHetPosition = all1
end function getFirstHetPosition


!---------------------------------------------------------------------------
!< @brief gets reference alleles from ref file if it exists there
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine readRef(file, plinkInfo)
	use AlphaHouseMod, only :countLines
	character(len=*),intent(in) :: file
	type(plinkInfoType),intent(inout) :: plinkInfo

	character(len=2) :: refAllele, minor,id
	logical :: fileExists
	integer :: lines,i, unit,snpId

	inquire( file=file, exist=fileExists )


	if (fileExists) then

		lines = countLines(file)
		open(newunit=unit, file=file, status="old")
		allocate(plinkInfo%referenceAllelePerSnps(lines))
		allocate(plinkInfo%alternateAllelePerSnps(lines))
		do i=1,lines
			read(unit, *) id, refAllele, minor
			snpId = plinkInfo%dict%getValue(id)
			if (snpId /= DICT_NULL) then
				plinkInfo%referenceAllelePerSnps(snpId) = refAllele
			else
				write(error_unit, *) "WARNING: snp " , id, " not found in map file"
			endif
		enddo
		close(unit)
	else
		write(error_unit,*) "Warning: no .ref file found for this dataset"
	endif

end subroutine readRef

!---------------------------------------------------------------------------
!< @brief writes out Ped file of results
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    November 26, 2017
!---------------------------------------------------------------------------
subroutine writePedFile(ped,plinkInfo,outputPath, paths)
	use AlphaHouseMod, only : countColumns
	use basespecfileModule

	type(plinkInfoType), intent(inout) :: plinkInfo
	type(PedigreeHolder) :: ped
	character(len=128), optional,dimension(:), allocatable, intent(in) :: paths !< output paths for each chromosome
	character(len=2), dimension(:,:), allocatable :: outputAlleles
	character(len=128) :: fmt
	character(len=*),intent(in) :: outputPath
	integer(kind=1) :: phase1,phase2
	integer :: snpCounts = 0,outcounts=0, pedUnit,i,nsnps,p,j


	print *, "merging plink output"
	if (.not. present(paths)) then

		call plinkInfo%initPlinkInfoType(ped)
	endif

	allocate(outputAlleles(ped%pedigreeSize, plinkInfo%totalSnps*2)) !outputphase

	outputAlleles = "0"

	! TODO talk with mr Gottardo about this
	if (.not. allocated(plinkInfo%referenceAllelePerSnps)) then
		allocate(plinkInfo%referenceAllelePerSnps(plinkInfo%totalSnps))
		allocate(plinkInfo%alternateAllelePerSnps(plinkInfo%totalSnps))
		plinkInfo%referenceAllelePerSnps = "1"
		plinkInfo%alternateAllelePerSnps = "2"
	endif

	do i =1,plinkInfo%nChroms

		nsnps = plinkInfo%nsnpsPerChromosome(i)

		if (present(paths)) then !< if multiple chromosomes - load in that file
			call ped%addPhaseInformationFromFile(trim(paths(i))//"phase.txt",nsnps)
		endif
		do j=1,nsnps
			snpCounts = snpCounts + 1
			outcounts = outcounts + 1
			do p=1,ped%pedigreeSize
				phase1 = ped%pedigree(p)%individualPhase(1)%getPhase(j)
				phase2 = ped%pedigree(p)%individualPhase(2)%getPhase(j)
				if (phase1 == 1) then
					outputAlleles(p,outcounts) = plinkInfo%referenceAllelePerSnps(snpCounts)
				else if (phase1 == 0) then
					outputAlleles(p,outcounts) = plinkInfo%alternateAllelePerSnps(snpCounts)
				else
					outputAlleles(p,outcounts) = '0'
				endif

				if (phase2 == 1) then
					outputAlleles(p,outcounts+1) = plinkInfo%referenceAllelePerSnps(snpCounts)
				else if (phase2 == 0) then
					outputAlleles(p,outcounts+1) = plinkInfo%alternateAllelePerSnps(snpCounts) !< -1 because snp are for output
				else
					outputAlleles(p,outcounts+1) = '0'
				endif

			enddo
			outcounts = outcounts + 1
		enddo
	enddo

	print *, "Writing plink output to file"
	open(newunit=pedUnit,file=trim(outputPath)//".ped", status='unknown')
	write(fmt, '(a,i10,a)') '(3a20,i2,a20,',2*plinkInfo%totalSnps, 'a2)'
	print *, "Writing alleles per animal"
	do p=1, ped%pedigreeSize
		write(pedUnit,  fmt) ped%pedigree(p)%originalId,ped%pedigree(p)%sireId,ped%pedigree(p)%damId,ped%pedigree(p)%gender,'0 ', outputAlleles(p,:)
	enddo

	close(pedunit)


end subroutine writePedFile


subroutine writeMapFile(plinkInfo, path)
	type(plinkInfoType), intent(in) :: plinkInfo
	character(len=*), intent(in) :: path
	integer :: unit, snpCount,i,h
	open(newunit=unit, file=trim(path)//".map", status='unknown')
	snpCount = 0


	do i=1, plinkInfo%nChroms
		do h=1, plinkInfo%nsnpsPerChromosome(i)
			snpCount = snpCount + 1
			write(unit,'(I4,a30,F10.5, I10)') i,trim(plinkInfo%snpName(snpCount)),plinkInfo%lengths(i,h),plinkInfo%basepairs(i,h)
		enddo
	enddo

	close(unit)


end subroutine writeMapFile


subroutine writeRefFile(plinkInfo, path)
	type(plinkInfoType), intent(in) :: plinkInfo
	character(len=*), intent(in) :: path
	integer :: unit, snpCount,i,h
	open(newunit=unit, file=trim(path)//".ref", status='unknown')
	snpCount = 0
	do i=1, plinkInfo%nChroms

		do h=1, plinkInfo%nsnpsPerChromosome(i)
			snpCount = snpCount + 1
			write(unit,'(a30,a8,a8)') trim(plinkInfo%snpName(snpCount)),trim(plinkInfo%referenceAllelePerSnps(snpCount)),trim(plinkInfo%alternateAllelePerSnps(snpCount))
		enddo
	enddo

	close(unit)


end subroutine writeRefFile


subroutine WriteBedFile(bed, minor, genotypes)

	use PedigreeModule
	use genotypeModule
	implicit none

	! Arguments
	character(len=*), intent(in) :: bed !< bed file name
	integer, intent(in) ::  minor !< if the first allele is minor or major
	integer(kind=1),dimension(:,:), intent(in) :: genotypes !< Genotypes should be in format (nanimals, nsnps)

	!! Types
	INTEGER, PARAMETER :: Byte = SELECTED_INT_KIND(1) ! Byte

	!! Local arguments
	integer(Byte) :: element, plinkmode
	integer(Byte), dimension(2) ::  magicnumber
	integer :: stat, i, animals, snps
	integer, dimension(4) :: codes, phasecodes
	!integer, dimension(:), allocatable :: domasksnps
	integer :: bedUnit
	! Supported formats as per plink 1.9.
	!data magicnumber/X"6C",X'0000001B' /,  plinkmode/X'01'/
	data magicnumber /108,27/, plinkmode /1/


	! TODO have to write function to convert ped to to plinkInfo.

	if (minor == 1) then
		codes = (/ 0, 1, 2, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 0, 1, 1, MISSINGGENOTYPECODE /)
	else
		codes = (/ 2, 1, 0, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 1, 1, 0, MISSINGGENOTYPECODE /)
	endif

	open(newunit=bedUnit, file=bed, status='REPLACE', ACCESS='STREAM', FORM='UNFORMATTED')
	write(bedUnit) magicnumber, plinkmode

	element = 0
	animals = 0
	snps = 1	
	outer: do
		inner: do i=0,6,2
			if (snps > size(genotypes,2)) exit outer
			animals = animals + 1
			if (genotypes(animals,snps) == codes(1)) then
				element = ibclr(element,i)
				element = ibclr(element,i+1)
			else if (genotypes(animals,snps) ==codes(4)) then
				element = ibset(element,i)
				element = ibclr(element,i+1)
			else if (genotypes(animals,snps) == codes(2)) then
				element = ibclr(element,i)
				element = ibset(element,i+1)
			else if (genotypes(animals,snps) ==codes(3)) then
				element = ibset(element,i)
				element = ibset(element,i+1)
			end if
			if (animals == size(genotypes,1)) then
				animals = 0
				snps = snps +1				
				exit
				
			endif
		enddo inner
		write(bedUnit, iostat=stat) element
	enddo outer
	close(bedUnit)


	print *, "finished"

end subroutine WriteBedFile


subroutine writeFamFile(ped,famFile)
	use PedigreeModule

	type(pedigreeHolder), intent(in) :: ped !< Pedigree object that is returned
	character(len=*), intent(in) :: famFile !< .ped file generated by plink
	character(len=2) :: phenotype
	integer :: fileUnit, i
	
	phenotype = "-9"

	open(newUnit=fileUnit, file=famFile, status="REPLACE")
	do i=1, ped%addedRealAnimals
		write(fileUnit,'(4a32,a1, i3,a1, a3)') ped%pedigree(ped%inputMap(i))%familyID,ped%pedigree(ped%inputMap(i))%originalId,ped%pedigree(ped%inputMap(i))%sireId,ped%pedigree(ped%inputMap(i))%damId," ",ped%pedigree(ped%inputMap(i))%gender," ",phenotype
	enddo


	close(fileUnit)
end subroutine writeFamFile



subroutine writeBimFile(bimFile, bimInfo)
	use HashModule
	use AlphaHouseMod
	use ConstantModule

	character(len=*), intent(in) :: bimFile
	type(bimHolder) , allocatable, dimension(:), intent(in) :: bimInfo !< extra info provided by BIM file
	integer :: i, unit

	open(newUnit=unit, file=bimFile, status="REPLACE")

	do i =1, size(bimInfo)
		write(unit, '(a2,a10,I6,I6,a1,a2,a2)') bimInfo(i)%chrom, bimInfo(i)%id,bimInfo(i)%chrompos, bimInfo(i)%pos ," ",bimInfo(i)%ref, bimInfo(i)%alt
	end do

	close (unit)

end subroutine writeBimFile


subroutine writeOutPlinkBinary(ped,path,bimInfo)

	type(PedigreeHolder) :: ped
	character(len=*), intent(in) :: path !< file prepend (before the .)
	type(bimHolder),dimension(:), allocatable, optional  :: bimInfo

	if (present(bimInfo)) then
		call writeBimFile(trim(path)//".bim", bimInfo)
	else
		call writeBimFile(trim(path)//".bim", createBimInfoFromGenotypes(ped%getGenotypesAsArrayWitHMissing()))
	endif 

	call writeFamFile(ped,trim(path)//".fam")

	call WriteBedFile(trim(path)//".bed",2, ped%getGenotypesAsArrayRealAnimals())


end subroutine writeOutPlinkBinary


subroutine writeOutPlinkNonBinary(ped,path,plinkInfoIn)

	type(PedigreeHolder) :: ped
	character(len=*), intent(in) :: path !< file prepend (before the .)
	type(plinkInfoType), optional :: plinkInfoIn
	type(plinkInfoType) :: plinkInfoObj
	
	if (present(plinkInfoIn)) then
		plinkInfoObj = plinkInfoIn
	else
		call plinkInfoObj%initPlinkInfoType(ped)
	endif
	
	
	call writePedFile(ped,plinkInfoObj,path)
	call writeMapFile(plinkInfoObj,path)
	call writeRefFile(plinkInfoObj,path)

end subroutine writeOutPlinkNonBinary



end module CompatibilityModule





















