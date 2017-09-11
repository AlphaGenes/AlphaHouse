module CompatibilityModule
	use integerLinkedListModule

	type bimHolder
	character(len=1) :: ref,alt
	character(len=IDLENGTH) :: id
	integer :: chrom
	integer(kind=int64) :: pos, chrompos
end type bimHolder

type Chromosome

type(integerLinkedList) :: snps
contains
	final :: destroyChrom
end type Chromosome
contains
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


	subroutine destroyChrom(chrom)


		type(Chromosome) :: chrom

		call chrom%snps%destroyLinkedList()

	end subroutine





	function readToPedigreeFormat(pedFile) result(ped)
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


			! write(*,'(3a20)') pedArray(1,i),pedArray(2,i),pedArray(3,i)
			read(gender,*,iostat=stat)  genderArray(i)
			read(phenotype,*,iostat=stat)  phenotypeArray(i)
		enddo


		ped = PedigreeHolder(pedArray, genderArray)

		call ped%printPedigreeOriginalFormat("pedigreeOutput.txt")

	end function readToPedigreeFormat


	subroutine readPlink(binaryFilePre, ped, outputPaths)
		use HashModule
		use AlphaHouseMod, only : CountLines
		use PedigreeModule
		use ifport

		character(len=*),intent(in) :: binaryFilePre
		type(pedigreeholder), intent(out) :: ped
		type(DictStructure) :: dict
		integer:: maxChroms
		type(bimHolder) , allocatable, dimension(:) :: bimInfo
		type(Chromosome), dimension(:), allocatable :: chroms
		logical, dimension(:), allocatable :: maskedLogi
		integer, dimension(:), allocatable :: masked
		character(len=100) :: fmt
		integer :: nsnps
		integer(kind=1), dimension(:,:), allocatable ::  allsnps
		character(len=128) :: path, outChrFile
		character(len=128), dimension(:), allocatable :: outputPaths
		integer result,i,p,h,outChrF
		! integer, allocatable


		ped = readToPedigreeFormat(binaryFilePre//".fam")

		call readBim(binaryFilePre//".bim",dict,bimInfo,nsnps,chroms,maxChroms)
		print *,"READ BIM"
		call readplinkSnps(binaryFilePre//".bed",nsnps,ped,1, allsnps)
		print *,"READ BED"
		allocate(maskedLogi(size(allSnps(1,:))))

		path = "chromosomeGenotypes/"
		result=makedirqq(path)

		print *, "MAX CHROMS",maxChroms

		allocate(outputPaths(maxChroms))
		do i =1, maxChroms
			write(outChrFile, '(a,a,i2,a)') trim(path),trim("chr"),i,DASH
			outputPaths(i) = outChrFile
			open(newunit=outChrF, file=outChrFile//"genotypes.txt", status="unknown")
			masked = chroms(i)%snps%convertToArray()
			maskedLogi = .false.
			do h =1, size(masked)
				maskedLogi(masked(h)) = .true.

			enddo
			write(fmt, '(a,i10,a)')  '(a20,', size(allSnps(1,:)), 'i2)'
			do p=1,ped%pedigreeSize-ped%nDummys
				! print *, "in write"
				write(outChrF,fmt) ped%pedigree(p)%originalId,pack(allSnps(p,:), maskedLogi)
			end do
			close(outChrF)
		enddo

		call dict%destroy()
	end subroutine readPlink







	subroutine readBim(bimFile, dict, bimInfo,nsnps,chroms, maxChroms)
		use HashModule
		use AlphaHouseMod

		character(len=*), intent(in) :: bimFile
		type(DictStructure) :: dict
		character :: ref,alt
		character(len=IDLENGTH) :: id

		integer(kind=int64) :: pos, chrompos

		integer,intent(out) :: nsnps
		integer,intent(out) :: maxChroms

		integer :: i, unit,chrom
		type(bimHolder) , allocatable, dimension(:), intent(out) :: bimInfo
		type(Chromosome),dimension(:), allocatable, intent(out) :: chroms
		maxChroms = 0

		dict = DictStructure()
		nsnps = countLines(bimFile)

		open(newUnit=unit, file=bimFile, status='old')
		allocate(chroms(LARGECHROMNUMBER))
		allocate(bimInfo(nsnps))
		do i =1, nsnps

			read(unit, *) chrom, id,chrompos, pos ,ref, alt

			if (chrom > maxChroms) then
				maxChroms = chrom
			endif
			call dict%addKey(id, i)
			bimInfo(i)%chrom = chrom
			!  TODO what does this do!???

			bimInfo(i)%id = id


			bimInfo(i)%chrompos = chrompos
			bimInfo(i)%pos = pos
			bimInfo(i)%ref = ref
			bimInfo(i)%alt = alt

			call chroms(chrom)%snps%list_add(i)
		end do
		close (unit)


	end subroutine readBim

	! Stores entire genotype file in memory
	subroutine readplinkSnps(bed, ncol,ped, minor,snps)

		use PedigreeModule
		use genotypeModule
		implicit none

		! Arguments
		character(*), intent(in) :: bed
		integer, intent(in) :: ncol, minor
		integer :: status
		type(PedigreeHolder), intent(in) :: ped

		!! Types
		INTEGER, PARAMETER :: Byte = SELECTED_INT_KIND(1) ! Byte

		!! Local arguments
		integer(Byte) :: readplinkmode, element, plinkmode
		integer(Byte), dimension(2) :: readmagicnumber, magicnumber
		!logical :: checkmaf
		integer :: stat, i, j, k, snpcount, majorcount
		integer, dimension(4) :: codes
		!integer, dimension(:), allocatable :: domasksnps
		integer(kind=1), dimension(:,:), allocatable,intent(out) ::  snps
		real :: allelefreq
		character(100) :: nChar, fmt
		integer :: na
		integer :: bedInUnit
		! Supported formats as per plink 1.9.
		!data magicnumber/X"6C",X'0000001B' /,  plinkmode/X'01'/
		data magicnumber /108,27/, plinkmode /1/

		
		print *,"start BED read"
		allocate(snps(ped%pedigreeSize-ped%nDummys,ncol))

		na = 9
		snps(:,:) = 9


		if (minor == 1) then
			codes = (/ 0, 1, 2, na /)
		else
			codes = (/ 2, 1, 0, na /)
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
					snps(j,k) = codes(1)
					case (1) ! missing
					snps(j,k) = codes(4)
					snpcount = snpcount - 1
					case (2) ! heterozygote
					snps(j,k) = codes(2)
					majorcount = majorcount + 1
					case (3) ! homozygote, minor
					snps(j,k) = codes(3)
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
		

		! Write output
		! fmt='(i20,'//trim(adjustl(nChar))//'I2)'
		!print *, fmt



		! do i=1,nlines
		! 	ped%pedigree(i)%individualGenotype = NewGenotypeInt(pack(snps(i,:), masksnps))
		! enddo

		! deallocate(snps)

		!print *, 'readplinksimple is done.'

	end subroutine readplinkSnps


end module CompatibilityModule




