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
	character(len=1) :: ref,alt
	character(len=IDLENGTH) :: id
	character(len=2) :: chrom !<either an integer, or 'X'/'Y'/'XY'/'MT'
	integer(kind=int64) :: pos, chrompos
end type bimHolder

type Chromosome

type(integerLinkedList) :: snps
end type Chromosome
contains






!---------------------------------------------------------------------------
!< @brief reads fam file to pedigree object
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
function readFamFile(pedFile) result(ped)
	use ConstantModule, only : IDLENGTH
	use AlphaHouseMod, only : countLines
	use PedigreeModule

	type(pedigreeHolder) :: ped !< Pedigree object that is returned
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



		! TODO potentially add family ID here to animals
		! pedArray(1,i) = familyID // ":" // pedArray(1,i)
		! pedArray(2,i) = familyID // ":" // pedArray(2,i)
		! pedArray(3,i) = familyID // ":" // pedArray(3,i)

		! write(*,'(3a20)') pedArray(1,i),pedArray(2,i),pedArray(3,i)
		read(gender,*,iostat=stat)  genderArray(i)
		read(phenotype,*,iostat=stat)  phenotypeArray(i)
	enddo


	call initPedigreeArrays(ped,pedArray, genderArray)

	print *, "ANS in ped",ped%pedigreeSize," without dum:",ped%pedigreeSize-ped%nDummys

	call ped%printPedigreeOriginalFormat("pedigreeOutput.txt")

end function readFamFile

!---------------------------------------------------------------------------
!< @brief Reads PLINK binary format into ped datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readPlink(binaryFilePre, ped, outputPaths,nsnps, sexChrom)
	use HashModule
	use PedigreeModule

	character(len=*),intent(in) :: binaryFilePre !< part before file extension
	type(pedigreeholder), intent(out) :: ped !< pedigree object returned
	logical, intent(out) :: sexChrom !< true if a sex chromosome is present
	integer,dimension(:), allocatable, intent(out) :: nsnps !< number of snps per chromosome
	character(len=128), dimension(:), allocatable, intent(out) :: outputPaths !< output paths for each chromosome

	type(DictStructure) :: dict
	integer:: maxChroms
	type(bimHolder) , allocatable, dimension(:) :: bimInfo
	type(Chromosome), dimension(:), allocatable :: chroms
	integer(kind=1), dimension(:,:), allocatable ::  genotypes
	integer(kind=1), dimension(:,:,:), allocatable ::  phase
	integer :: totalSnps
	real(kind=real64), dimension(:,:) ,allocatable:: lengths
	integer, dimension(:,:) ,allocatable :: basepairs

	! TODO change to pointer rather than copy
	ped = readFamFile(trim(binaryFilePre)//".fam")
	call readBim(trim(binaryFilePre)//".bim",dict,bimInfo,nsnps,totalSnps,chroms,maxChroms, sexChrom,lengths,basepairs)
	print *,"READ BIM"
	call readBED(trim(binaryFilePre)//".bed",totalSnps,ped,4, genotypes, phase)
	print *,"READ BED"

	call createOutputFiles(genotypes,chroms, phase,lengths,basepairs,maxChroms,nsnps,ped, outputPaths)
end subroutine readPlink

!---------------------------------------------------------------------------
!< @brief writes output files
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine createOutputFiles(genotypes,chroms, phase,lengths, basepairs,maxChroms,nsnps,ped, outputPaths,referenceAllelePerSnps)
	use ifport
	use ALphaHouseMod, only : countLines

	integer(kind=1), dimension(:,:), allocatable,intent(in) ::  genotypes	 !<  genotypes for all chroms (nanimals, totalsnps)
	integer(kind=1), dimension(:,:,:), allocatable,intent(in) ::  phase !< array of phase (nanimals, totalsnps,2 )
	real(kind=real64), dimension(:,:) ,allocatable,intent(in) :: lengths
	integer, dimension(:,:) ,allocatable,intent(in) :: basepairs
	integer, intent(in) :: maxChroms
	integer,dimension(:), allocatable, intent(in) :: nsnps
	type(Chromosome), dimension(:), allocatable, intent(in) :: chroms
	character(len=1),dimension(:), allocatable,optional, intent(in) :: referenceAllelePerSnps !<array saying for which snp the reference allele is
	! character(len=1),dimension(:), allocatable :: alleles
	type(pedigreeholder), intent(in) :: ped !< pedigree object taken in
	character(len=128), dimension(:), allocatable,intent(out) :: outputPaths !< returns the output path for each chromosome

	integer :: outChrF, outChrP,bpFile,lengthFile
	character(len=128) :: path, outChrFile,fmt
	logical, dimension(:), allocatable :: maskedLogi
	integer, dimension(:), allocatable :: masked
	integer :: result, i, h,p,refAlleleUnit,count
	integer(kind=1),allocatable,dimension (:,:) :: array
	path = "fullGenome/"
	result=makedirqq(path)

	print *, "Total number of Chromosomes:",maxChroms

	allocate(outputPaths(maxChroms))
	allocate(maskedLogi(size(genotypes,2))) !< alloc to total number of snps
	do i =1, maxChroms

		write(outChrFile, '(a,a,i0.2,a)') trim(path),trim("chr"),i,DASH

		result=makedirqq(outChrFile)
		outputPaths(i) = outChrFile

		open(newunit=outChrF, file=trim(outChrFile)//"genotypes.txt", status="unknown")
		open(newunit=outChrP, file=trim(outChrFile)//"phase.txt", status="unknown")
		open(newunit=lengthFile, file=trim(outChrFile)//"snplengths.txt", status="unknown")
		open(newunit=bpFile, file=trim(outChrFile)//"snpBasepairs.txt", status="unknown")

		!< uses genotype mask to determine which snps to use
		masked = chroms(i)%snps%convertToArray()
		maskedLogi = .false.
		do h =1, size(masked)
			maskedLogi(masked(h)) = .true.
		enddo

		if (present(referenceAllelePerSnps)) then
			open(newunit=refAlleleUnit, file=trim(outChrFile)//"refAlleles", status="unknown")

			! alleles = pack(referenceAllelePerSnps, maskedLogi)

			write(refAlleleUnit, '(3a5)') "snp in Chrom", "snp pos","Ref Allele"
			count = 0
			do p=1, size(referenceAllelePerSnps)
				! do p=1,size(alleles)
				! write(refAlleleUnit, '(1i5,a5)') p, referenceAllelePerSnps(p)
				if (maskedLogi(p)) then
					count = count +1
					write(refAlleleUnit, '(2i5,a5)') count,p, referenceAllelePerSnps(p)
				endif
			enddo
			close(refAlleleUnit)
		endif

		allocate(array(ped%pedigreeSize-ped%nDummys, nsnps(i)))
		write(fmt, '(a,i10,a)')  "(a20,", nsnps(i), "i3)"
		! set up the pedigree to avoid read in
		array = 9
		do p=1,ped%pedigreeSize-ped%nDummys
			array(p,:) = pack(genotypes(p,:), maskedLogi)
		end do

		write(lengthFile, *) pack(lengths(i,:), maskedLogi)
		write(bpFile, *) pack(basepairs(i,:), maskedLogi)

		do p=1,ped%pedigreeSize-ped%nDummys
			write(outChrF,fmt) ped%pedigree(p)%originalId,array(p,:)
			write(outChrp,fmt) ped%pedigree(p)%originalId,pack(phase(p,:,1), maskedLogi)
			write(outChrp,fmt) ped%pedigree(p)%originalId,pack(phase(p,:,2), maskedLogi)
		end do
		print *, "Finished writeout of chromosome ", i

		close(outChrF)
		close(outChrF)
		close(outChrP)
		close(lengthFile)
		close(bpFile)
		deallocate(array)
	enddo


end subroutine createOutputFiles



!---------------------------------------------------------------------------
!< @brief Reads bim file to datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readBim(bimFile, dict, bimInfo,nsnps,totalSnps,chroms, maxChroms, hasSexChrom,lengths,basepairs)
	use HashModule
	use AlphaHouseMod
	use ConstantModule

	character(len=*), intent(in) :: bimFile
	type(DictStructure),intent(out) :: dict !< dicitionary tying identifier to snp
	type(bimHolder) , allocatable, dimension(:), intent(out) :: bimInfo !< extra info provided by BIM file
	type(Chromosome),dimension(:), allocatable, intent(out) :: chroms !< chromosome info
	integer,intent(out), dimension(:), allocatable :: nsnps !< nsnps per chromosome
	logical, intent(out) :: hasSexChrom !< true if there is a sex chrom
	integer,intent(out) :: maxChroms
	real(kind=real64), dimension(:,:) ,allocatable, intent(out) :: lengths !< lengths of nsnp in morgan format is (chrom, totalsnp)
	integer, dimension(:,:) ,allocatable,intent(out) :: basepairs !< basepairs format is (chrom, totalsnp)

	character :: ref,alt
	character(len=IDLENGTH) :: id

	integer(kind=int64) :: pos
	real(kind=real64) :: chrompos

	integer,intent(out) :: totalSnps
	real(kind=real64), dimension(:,:) ,allocatable :: tmpLengths
	integer, dimension(:,:) ,allocatable :: tmpbasePairs
	integer, dimension(:), allocatable :: temparray

	integer :: i, unit,chromCount,curChromSnpCount
	character(len=2) :: chrom,prevChrom

	maxChroms = 0
	hasSexChrom = .false.
	curChromSnpCount = 0
	allocate(nsnps(LARGECHROMNUMBER))
	nsnps = 0
	call dict%DictStructure()
	totalSnps = countLines(bimFile)

	open(newUnit=unit, file=bimFile, status='old')
	allocate(chroms(LARGECHROMNUMBER))
	allocate(bimInfo(totalSnps))
	allocate(lengths(LARGECHROMNUMBER,totalSnps))
	allocate(basepairs(LARGECHROMNUMBER,totalSnps))
	chromCount = 1
	do i =1, totalSnps

		read(unit, *) chrom, id,chrompos, pos ,ref, alt

		if (i == 1) then
			prevChrom = chrom
		endif

		! if we've moved on to the next chromsome
		if (chrom == "X" .or. chrom == 'Y') then
			hasSexChrom = .true.
		endif

		if (chrom == 'XY' .or. chrom == 'MT') then
			write(error_unit,*) "WARNING - No support currently for XY or MT chromosomes"
		endif

		! if we 've moved on to the next chromosome
		if (chrom /=prevChrom .or. i == totalSnps) then

			! set the count to he current numbers
			if (i == totalSnps) then
				curChromSnpCount = curChromSnpCount + 1
			endif
			nsnps(chromCount) = curChromSnpCount
			curChromSnpCount = 0
			chromCount = chromCount + 1

			prevChrom = chrom
			if (chromCount > maxChroms) then
				maxChroms = chromCount
			endif
		endif

		curChromSnpCount = curChromSnpCount + 1

		lengths(chromCount,curChromSnpCount) = chromPos
		basepairs(chromCount,curChromSnpCount) = pos

		call dict%addKey(id, i)

		bimInfo(i)%chrom = chrom
		!  TODO clean this up

		bimInfo(i)%id = id


		bimInfo(i)%chrompos = chrompos
		bimInfo(i)%pos = pos
		bimInfo(i)%ref = ref
		bimInfo(i)%alt = alt

		call chroms(chromCount)%snps%list_add(i)
	end do

	maxChroms = maxChroms -1
	chromCount = chromCount - 1

	if (chromCount /= LARGECHROMNUMBER) then
		allocate(temparray(chromCount))
		temparray(1:chromCount) = nsnps(1:chromCount)
		call move_alloc(temparray,nsnps)

		allocate(tmpLengths(chromCount,totalSnps))
		tmpLengths(1:chromCount,:) = lengths(1:chromCount,:)
		call move_alloc(tmpLengths,lengths)


		allocate(tmpbasePairs(chromCount,totalSnps))
		tmpbasePairs(1:chromCount,:) = basepairs(1:chromCount,:)
		call move_alloc(tmpbasePairs,basepairs)
	endif
	close (unit)


end subroutine readBim

!---------------------------------------------------------------------------
!< @brief Reads BED Files
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2017
!---------------------------------------------------------------------------
subroutine readBED(bed, totalSnps,ped, minor,genotypes, phase)

	use PedigreeModule
	use genotypeModule
	implicit none

	! Arguments
	character(*), intent(in) :: bed !< bed file name
	type(PedigreeHolder), intent(in) :: ped !< pedigree read in from .fam
	integer, intent(in) :: totalSnps, minor !< totalsnps, and if the first allele is minor or major
	integer(kind=1), dimension(:,:), allocatable,intent(out) ::  genotypes !< total genotypes (nanimals, nsnps)
	integer(kind=1), dimension(:,:,:), allocatable,intent(out) ::  phase

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
	allocate(genotypes(ped%pedigreeSize-ped%nDummys,totalSnps))
	allocate(phase(ped%pedigreeSize-ped%nDummys,totalSnps,2))
	genotypes(:,:) = 9


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
				genotypes(j,k) = codes(1)
				phase(j,k,:) = phasecodes(1)
				case (1) ! missing
				genotypes(j,k) = codes(4)
				phase(j,k,:) = phaseCodes(4)
				snpcount = snpcount - 1
				case (2) ! heterozygote
				genotypes(j,k) = codes(2)
				! Set to missing - as we don't know which snp is which
				phase(j,k,1) = MISSINGPHASECODE
				phase(j,k,2) = MISSINGPHASECODE
				majorcount = majorcount + 1
				case (3) ! homozygote, minor
				genotypes(j,k) = codes(3)
				phase(j,k,:) = phasecodes(3)
				majorcount = majorcount + 2
			endselect
			if (j == ped%pedigreeSize-ped%nDummys) then
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
subroutine readPlinkNoneBinary(filePre,ped,outputPaths ,nsnps,sexChrom)
	use HashModule
	use PedigreeModule
	use AlphaHouseMod
	use ifport

	character(len=*),intent(in) :: filePre
	character(len=128), dimension(:), allocatable,intent(out) :: outputPaths
	integer, dimension(:) ,allocatable :: nsnps
	real(kind=real64), dimension(:,:) ,allocatable :: lengths
	integer, dimension(:,:) ,allocatable :: basepairs
	integer :: totalSnps
	type(Chromosome), dimension(:), allocatable :: chroms
	logical,intent(out) :: sexChrom
	character(len=2),dimension(:), allocatable :: referenceAllelePerSnps !<array saying for which snp the reference allele is


	integer(kind=1), dimension(:,:), allocatable ::  genotypes
	integer(kind=1), dimension(:,:,:), allocatable ::  phase

	integer :: maxChroms
	type(DictStructure) :: dict
	type(pedigreeHolder) :: ped

	call readMap(trim(filePre)//".map", dict,chroms,maxChroms, nsnps, totalSnps,sexChrom, lengths, basepairs)
	call readRef(trim(filePre)//".ref", dict, referenceAllelePerSnps)
	call readPedFile(trim(filePre)//".ped",ped, totalSnps, genotypes,phase, referenceAllelePerSnps)
	call createOutputFiles(genotypes,chroms, phase,lengths,basepairs,maxChroms,nsnps,ped, outputPaths,referenceAllelePerSnps)

end subroutine readPlinkNoneBinary


!---------------------------------------------------------------------------
!< @brief Reads plink map file into datastructures
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readMap(filename,dict,chroms,maxChroms, snpCounts, totalSnps,hasSexChrom,lengths, basepairs)
	use HashModule
	use AlphaHouseMod

	character(len=*),intent(in) :: filename
	integer, dimension(:) ,allocatable, intent(out) :: snpCounts !< number of snps for each chromosome
	integer, intent(out) :: maxChroms
	type(DictStructure), intent(out) :: dict !< map dictionary to accompany ref file
	type(Chromosome),dimension(:), allocatable, intent(out) :: chroms !< individual info for each chromosome
	real(real64), dimension(:,:) ,allocatable,intent(out) :: lengths !< array of snp lengths in morgans, in format (chrom, snps)
	integer, dimension(:,:) ,allocatable,intent(out) :: basepairs !< array of basepair positions, in format (chrom, snps)


	integer, dimension(:) ,allocatable :: tempArray
	real(real64), dimension(:,:) ,allocatable ::tmpLengths
	integer, intent(out) :: totalSnps
	integer, dimension(:,:) ,allocatable :: tmpbasePairs
	integer :: unit,i,chromCount,basepair

	real(kind=real64) :: length
	logical :: hasSexChrom
	character(len=2) :: chrom,prevChrom
	character(len=128) :: id

	totalSnps = countLines(fileName)
	print *,"TOTAL SNPS IN MAP:", totalSnps
	open(newunit=unit, file=filename, status='OLD')

	allocate(chroms(LARGECHROMNUMBER))
	allocate(lengths(LARGECHROMNUMBER, totalSnps))
	allocate(basepairs(LARGECHROMNUMBER, totalSnps))

	allocate(snpCounts(LARGECHROMNUMBER))

	call dict%DictStructure()
	hasSexChrom = .false.
	maxChroms = 0
	snpCounts = 0
	chromCount = 0
	prevChrom = 'MT'

	write(*,*) "Start reading map file"
	do i=1,totalSnps

		read(unit, *) chrom, id,length, basepair

		if (chrom == "X" .or. chrom == 'Y') then
			hasSexChrom = .true.
		endif
		if (chrom == 'XY' .or. chrom == 'MT') then
			write(error_unit,*) "WARNING - No support currently for XY or MT chromosomes"
		endif
		if (chrom /=prevChrom) then
			chromCount = chromCount + 1
			prevChrom = chrom
		endif

		snpCounts(chromCount) = snpCounts(chromCount) + 1
		lengths(chromCount,snpCounts(chromCount)) = length
		basepairs(chromCount,snpCounts(chromCount)) = basepair
		call dict%addKey(id, i)
		call chroms(chromCount)%snps%list_add(i)
	enddo

	if (chromCount /= LARGECHROMNUMBER) then
		allocate(temparray(chromCount))
		temparray(1:chromCount) = snpCounts(1:chromCount)
		call move_alloc(temparray,snpCounts)

		allocate(tmpLengths(chromCount,totalSnps))
		tmpLengths(1:chromCount,:) = lengths(1:chromCount,:)
		call move_alloc(tmpLengths,lengths)


		allocate(tmpbasePairs(chromCount,totalSnps))
		tmpbasePairs(1:chromCount,:) = basepairs(1:chromCount,:)
		call move_alloc(tmpbasePairs,basepairs)
	endif

	maxChroms = chromCount
	write(*,*) "Finished reading map file"
end subroutine readMap


!---------------------------------------------------------------------------
!< @brief Sets the individual to be genotyped at high density.
!< @author  David Wilson david.wilson@roslin.ed.ac.uk
!< @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine readPedFile(filename,ped, totalSnps,genotypes,phase, referenceAllelePerSnps)
	use PedigreeModule

	use AlphaHouseMod

	character(len=*), intent(in) :: filename
	type(pedigreeHolder), intent(out) :: ped
	integer, intent(in) :: totalSnps

	character(len=IDLENGTH), dimension(:,:), allocatable :: pedArray

	character(len=2),dimension(:), allocatable, intent(inout) :: referenceAllelePerSnps !<array saying for which snp the reference allele is
	character(len=2),dimension(:,:), allocatable :: alleles !<(nanimals. nsnp*2) size nsnp x2 (for each allele,)
	integer(kind=1), dimension(:,:), allocatable,intent(out) ::  genotypes
	integer(kind=1), dimension(:,:,:), allocatable,intent(out) ::  phase
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



	allocate(alleles(size, totalSnps*2))
	allocate(phase(size,totalSnps,2))
	allocate(genotypes(size,totalSnps))

	open(newunit=fileUnit, file=fileName, status="unknown")

	! minor is always false atm - this is because we always want the same ref allele to be used here - potentially used for future
	minor = 4
	if (minor == 1) then
		codes = (/ 0, 1, 2, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 0, 1, 1, MISSINGGENOTYPECODE /)
	else
		codes = (/ 2, 1, 0, MISSINGGENOTYPECODE /)
		phaseCodes = (/ 1, 1, 0, MISSINGGENOTYPECODE /)
	endif
	write(*,*) "Start Reading Ped File"
	do i=1,size
		read(fileUnit,*) familyID(i),pedArray(1,i),pedArray(2,i),pedArray(3,i),gender,phenotype, alleles(i,:)

		read(gender,*,iostat=stat)  genderArray(i)
		read(phenotype,*,iostat=stat)  phenotypeArray(i)
	enddo

	close(fileUnit)

		! check if reference alleles have been passed in
	if (.not. allocated(referenceAllelePerSnps)) then
		allocate(referenceAllelePerSnps(totalSnps))


		do j=1,totalSnps*2,2
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
					referenceAllelePerSnps(cursnp) = one
				else
					referenceAllelePerSnps(cursnp) = two
				endif
			end block
		enddo
	endif

	do j=1,totalSnps*2,2
		cursnp = (j/2) + 1
		! ref is always the first - allele in snp -
		! referenceAllelePerSnps(cursnp) = alleles(i,j+1)
		if (.not. allocated(referenceAllelePerSnps)) then
			referenceAllelePerSnps(cursnp) = getFirstHetPosition(alleles, cursnp)
		endif
		do i=1,size !< loop through animals
			all1 = alleles(i,j)
			all2 = alleles(i,j+1)
			if (all1 == '0' .or. all2 == '0') then
				genotypes(i,curSnp) = MISSINGGENOTYPECODE

			else if (all1 == all2) then

				if (all1 == referenceAllelePerSnps(curSnp) ) then
					genotypes(i,curSnp) = codes(1)
					phase(i,cursnp,:) = phasecodes(1)
				else
					genotypes(i,curSnp) = codes(3)
					phase(i,cursnp,:) = phasecodes(3)

				endif
			else !< means they are different
				genotypes(i,curSnp) = codes(2)

				if (all1 == referenceAllelePerSnps(cursnp)) then
					phase(i,cursnp,1) = phasecodes(3)
					phase(i,cursnp,2) =  phasecodes(1)
				else
					phase(i,cursnp,1) = phasecodes(1)
					phase(i,cursnp,2) = phasecodes(3)
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
	integer :: i, h

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
subroutine readRef(file, dict, referenceAllelePerSnps)
	use AlphaHouseMod, only :countLines
	character(len=*),intent(in) :: file
	type(DictStructure), intent(in) :: dict
	character(len=2), dimension(:), allocatable, intent(out) :: referenceAllelePerSnps

	character(len=2) :: refAllele, minor,id
	logical :: fileExists
	integer :: lines,i, unit,snpId

	inquire( file=file, exist=fileExists )


	if (fileExists) then

		lines = countLines(file)
		open(newunit=unit, file=file, status="old")
		allocate(referenceAllelePerSnps(lines))
		do i=1,lines
			read(unit, *) id, refAllele, minor
			snpId = dict%getValue(id)
			if (snpId /= DICT_NULL) then
				referenceAllelePerSnps(snpId) = refAllele
			else
				write(error_unit, *) "WARNING: snp " , id, " not found in map file"
			endif
		enddo
		close(unit)
	else
		write(error_unit,*) "Warning: no .ref file found for this dataset"
	endif

end subroutine readRef

end module CompatibilityModule


















