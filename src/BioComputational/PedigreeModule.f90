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
!######
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     PedigreeModule.f90
!
! DESCRIPTION:
!> @brief    Module` containing definition of pedigree.
!
!> @details  Fully doubly linked list with useful procedures for operations on the linked list
!
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 Dwilson - Initial Version

!-------------------------------------------------------------------------------




module PedigreeModule
	use IndividualModule
	use IndividualLinkedListModule
	use HashModule
	use constantModule
	use AlphaHouseMod, only : Int2Char
	implicit none


	private addOffspringsAfterReadIn

	type PedigreeHolder

	type(Individual), pointer, dimension(:) :: Pedigree

	type(IndividualLinkedList),allocatable :: Founders !linked List holding all founders
	type(IndividualLinkedList),allocatable, dimension(:) :: generations !linked List holding each generation
	type(DictStructure),allocatable :: dictionary ! hashmap of animal ids to index in pedigree
	integer(kind=int32) :: pedigreeSize, nDummys !pedigree size cannot be bigger than 2 billion animals
	integer(kind=int32) :: addedRealAnimals !< animals that have been added, that aren't dummys
	integer(kind=int32) :: unknownDummys !< dummys that have been set by having one unknown parent
	integer(kind=int32) :: maxPedigreeSize ! maximum size pedigree can be

	integer, dimension(:) , allocatable :: inputMap ! map going from inputMap(1:pedsize) = recID -- this is to maintain the order of the original pedigree

	integer, dimension(:) , allocatable :: genotypeMap ! map going from genotypeMap(1:nAnisG) = recID
	type(DictStructure),allocatable :: genotypeDictionary ! maps id to location in genotype map
	integer(kind=int32) :: nGenotyped ! number of animals that are genotyped

	integer, dimension(:) , allocatable :: hdMap ! map going from genotypeMap(1:nHd) = recID
	type(DictStructure),allocatable :: hdDictionary ! maps id to location in genotype map
	integer(kind=int32) :: nHd ! number of animals that are genotyped hd
	integer :: maxGeneration ! largest generation

	integer :: nsnpsPopulation ! number of snps for
	integer(kind=1) :: isSorted !integer saying how pedigree is sorted. 0 is not sorted. 1 is sorted with all dummys at end, 2 is sorted with unknown dummys at end, 3 is sorted with dummys at top

	type(IndividualLinkedList),allocatable :: sireList, damList !< lists containing all sires and dams
	type(IndividualLinkedList), allocatable :: uniqueParentList

	real(kind=real64), dimension(:), allocatable :: snpLengths !< length in morgans (size nsnp)
	integer, dimension(:), allocatable :: snpBasePairs !< base pairs size nsnp
	type(DictStructure),allocatable :: familyIdDict !< dictionary mapping family id to integer value

	contains
		procedure :: calculatePedigreeCorrelationNoInbreeding
		procedure :: calculatePedigreeCorrelationWithInbreeding
		procedure :: copyPedigree
		final :: destroyPedigree
		procedure :: setPedigreeGenerationsAndBuildArrays
		procedure :: outputSortedPedigree
		procedure :: setOffspringGeneration
		procedure :: addGenotypeInformationFromFile
		procedure :: addPhaseInformationFromFile
		procedure :: addGenotypeInformationFromArray
		generic :: addGenotypeInformation => addGenotypeInformationFromArray, addGenotypeInformationFromFile
		procedure :: outputSortedPedigreeInAlphaImputeFormat
		procedure :: isDummy
		procedure :: sortPedigreeAndOverwrite
		procedure :: sortPedigreeAndOverwriteWithDummyAtTheTop
		procedure :: makeRecodedPedigreeArray
		procedure :: printPedigree
		procedure :: printPedigreeOriginalFormat
		procedure :: getMatePairsAndOffspring
		procedure :: addSequenceFromVCFFile
		procedure :: getAllGenotypesAtPosition
		procedure :: getAllGenotypesAtPositionWithUngenotypedAnimals
		procedure :: getPhaseAtPosition
		procedure :: getPhaseAtPositionUngenotypedAnimals
		procedure :: getPhaseAsArrayWithMissing
		procedure :: setAnimalAsGenotyped
		procedure :: getGenotypesAsArray
		procedure :: getPhaseAsArray
		procedure :: getGenotypesAsArrayWitHMissing
		procedure :: countMissingGenotypesNoDummys
		procedure :: countMissingPhaseNoDummys
		procedure :: setGenotypeFromArray
		procedure :: setPhaseFromArray
		procedure :: getNumGenotypesMissing
		procedure :: getGenotypedFounders
		procedure :: cleanGenotypesBasedOnHaplotypes
		procedure :: getSireDamGenotypeIDByIndex
		procedure :: setAnimalAsHD
		procedure :: getSireDamHDIDByIndex
		procedure :: getGenotypePercentage
		procedure :: writeOutGenotypes
		procedure :: WriteoutPhase
		procedure :: WriteoutPhaseAll
		procedure :: WriteoutPhaseNoDummies
		procedure :: createDummyAnimalAtEndOfPedigree
		procedure :: addAnimalAtEndOfPedigree
		procedure :: addSequenceFromFile
		procedure :: setAnimalAsGenotypedSequence
		procedure :: convertSequenceDataToArray
		procedure :: getSequenceAsArrayWithMissing
		procedure :: writeOutGenders
		procedure :: readInGenders
		procedure :: MakeGenotype
		procedure :: PhaseComplement
		procedure :: homozygoticFillIn
		procedure :: wipeGenotypeAndPhaseInfo
		procedure :: getUniqueParents
		procedure :: findMendelianInconsistencies
		procedure :: writeOutGenotypesAll
		procedure :: writeOutGenotypesNoDummies
		procedure :: getHDPedigree
		procedure :: writeOutPedigree
		procedure :: writeOutGenotypeProbabilities
		procedure :: writeOutPhaseProbabilities
		procedure :: writeOutGenotypesForImputed
		procedure :: setSnpBasePairs
		procedure :: setSnpLengths
		procedure :: deepCheckPedigree
		procedure :: initPedigreeFromOutputFileFolder
		procedure :: setAnimalAsGenotypedFromPhase


#ifdef MPIACTIVE
		procedure:: calculatePedigreeCorrelationWithInBreedingMPI
#endif

	end type PedigreeHolder


	type RecodedPedigreeArray
	integer(kind=int32) :: nInd
	character(len=IDLENGTH), allocatable, dimension(:) :: originalId
	integer(kind=int32), allocatable, dimension(:) :: generation
	integer(kind=int32), allocatable, dimension(:,:) :: id
	contains
		procedure :: init    => initRecodedPedigreeArray
		procedure :: destroy => destroyRecodedPedigreeArray
		procedure :: write   => writeRecodedPedigreeArray
	end type

	! interface PedigreeHolder
	! 	module procedure initPedigree
	! 	module procedure initPedigreeArrays
	! 	module procedure initEmptyPedigree
	! 	module procedure initPedigreeGenotypeFiles
	! 	module procedure initPedigreeIntArrays
	! 	! initPedigreeFromOutputFileFolder initialiser from file
	! end interface PedigreeHolder

	interface assignment (=)
		module procedure copyPedigree
	end interface

	interface operator (==)
		module procedure equality
	end interface
	interface Sort !Sorts into generation list
		module procedure :: setPedigreeGenerationsAndBuildArrays
	end interface Sort
	contains

		!---------------------------------------------------------------------------
		!< @brief Output  of animals that are genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypesForImputed(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit

			open(newUnit=fileUnit,file=filename,status="unknown")
			do i= 1, this%pedigreeSize
				if (this%pedigree(i)%isimputed) then
					write(fileUnit,'(a20,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2,20000i2)')  this%pedigree(this%genotypeMap(i))%originalId, this%pedigree(i)%individualGenotype%toIntegerArray()
				endif
			enddo
			close(fileUnit)
		end subroutine writeOutGenotypesForImputed

		!---------------------------------------------------------------------------
		!< @brief Sets length array and overwrites if already set
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!---------------------------------------------------------------------------
		subroutine setSnpLengths(this, filepath, nsnpIn)
			use AlphaHouseMod, only : countColumns
			class(PedigreeHolder) :: this
			character(len=*),intent(in) :: filePath
			integer, optional :: nsnpIn
			integer :: unit,nsnp

			if (allocated(this%snpLengths)) then
				deallocate(this%snpLengths)
			endif
			if (present(nsnpIn)) then
				allocate(this%snpLengths(nsnpIn))
			else if (this%nsnpsPopulation /= 0) then
				allocate(this%snpLengths(this%nsnpsPopulation))
			else
				nsnp = countColumns(filePath,' ')
			endif

			open(newunit=unit, file=filePath, status="old")

			read(unit,*) this%snpLengths

			close(unit)
		end subroutine setSnpLengths


		!---------------------------------------------------------------------------
		!< @brief Sets base pair array and overwrites if already set
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!---------------------------------------------------------------------------
		subroutine setSnpBasePairs(this, filePath, nsnpIn)
			use AlphaHouseMod, only : countColumns
			class(PedigreeHolder) :: this
			character(len=*),intent(in) :: filePath
			integer, optional :: nsnpIn
			integer :: unit,nsnp
			if (allocated(this%snpBasePairs)) then
				deallocate(this%snpBasePairs)
			endif

			if (present(nsnpIn)) then
				allocate(this%snpBasePairs(nsnpIn))
			else if (this%nsnpsPopulation /= 0) then
				allocate(this%snpBasePairs(this%nsnpsPopulation))
			else
				nsnp = countColumns(filePath,' ')
			endif

			open(newunit=unit, file=filePath, status="old")

			read(unit,*) this%snpBasePairs

			close(unit)
		end subroutine setSnpBasePairs

		!---------------------------------------------------------------------------
		!< @brief Sorts pedigree, and overwrites all fields to new values
		!< @details effectively, does a deep copy to sort pedigree based on generation, but puts dummys at bottom
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine copyPedigree(res,this)
			use iso_fortran_env, only : output_unit, int64
			type(PedigreeHolder),intent(in) :: this
			class(pedigreeHolder), intent(inout) :: res
			integer :: i, tmpId, tmpSire, tmpDam, tmpAnimalArrayCount
			integer(kind=int64) :: sizeDict
			integer, allocatable, dimension(:) :: tmpAnimalArray


			call destroyPedigree(res)
			res%pedigreeSize = this%pedigreeSize
			res%nDummys = this%nDummys
			res%unknownDummys = this%unknownDummys
			res%maxPedigreeSize = this%maxPedigreeSize
			res%inputMap = this%inputMap
			res%genotypeMap = this%genotypeMap
			res%addedRealAnimals = this%addedRealAnimals

			allocate(tmpAnimalArray(this%pedigreeSize))

			res%nGenotyped =this%nGenotyped

			if (allocated(this%hdMap)) then
				if (.not. allocated(res%hdMap)) then
					allocate(res%hdMap(res%pedigreeSize))
				endif
				res%hdMap = this%hdMap
			endif

			! can copy genotype and hd dictionaries as order will be the same,
			! can't do main dictionary as we have to know if parents exist
			if (allocated(this%hdDictionary)) then
				if (allocated(res%hdDictionary)) then
					deallocate(res%hdDictionary)
				endif
				allocate(res%hdDictionary)
				res%hdDictionary = this%hdDictionary
			endif

			if (allocated(this%genotypeDictionary)) then
				if (.not. allocated(res%genotypeDictionary)) then
					allocate(res%genotypeDictionary)
				endif
				call deepCopyHash(res%genotypeDictionary,this%genotypeDictionary)
			endif

			res%nHd = this%nHd
			res%maxGeneration = this%maxGeneration
			res%nsnpsPopulation = this%nsnpsPopulation
			res%isSorted =this%isSorted

			allocate(res%sireList)
			allocate(res%damList)
			sizeDict  =this%pedigreeSize*2

			allocate(res%pedigree(this%maxPedigreeSize))
			allocate(res%dictionary)
			call res%dictionary%DictStructure(sizeDict)

			if (allocated(this%generations)) then
				allocate(res%Generations(0:this%maxGeneration))
			endif

			allocate(res%founders)
			tmpAnimalArrayCount = 0

			do i=1, this%pedigreeSize

				call res%dictionary%addKey(trim(this%pedigree(i)%originalID),i)
				call copyIndividual(res%pedigree(i),this%pedigree(i))
				call res%pedigree(i)%resetOffspringInformation
				if (.not. res%pedigree(i)%Founder) then
					! we know that sire and dam id should be set

					tmpSire= res%dictionary%getValue(res%pedigree(i)%sireId)
					tmpDam = res%dictionary%getValue(res%pedigree(i)%damId)

					! if either sire or dam is null - we need to wait for all the animals are in the new ped to set them
					if (tmpSire == DICT_NULL .or. tmpDam == DICT_NULL) then
						tmpAnimalArrayCount = tmpAnimalArrayCount + 1
						tmpAnimalArray(tmpAnimalArrayCount) = i
					else

						call res%pedigree(tmpSire)%addOffspring(res%pedigree(i))
						res%pedigree(i)%sirePointer =>  res%pedigree(tmpSire)
						if (res%pedigree(tmpSire)%nOffs== 1) then
							call res%sireList%list_add(res%pedigree(tmpSire))
						endif
						call res%pedigree(tmpDam)%addOffspring(res%pedigree(i))
						res%pedigree(i)%damPointer =>  res%pedigree(tmpDam)
						if (res%pedigree(tmpDam)%nOffs== 1) then
							call res%damList%list_add(res%pedigree(tmpDam))
						endif

					endif
				else
					call res%founders%list_add(res%pedigree(i))
				endif

				if (res%isSorted /= 0) then
					call res%generations(res%pedigree(i)%generation)%list_add(res%pedigree(i))
				endif
			enddo

			do i=1, tmpAnimalArrayCount
				tmpId = tmpAnimalArray(i)

				tmpSire= res%dictionary%getValue(res%pedigree(tmpId)%sirePointer%originalId)
				tmpDam = res%dictionary%getValue(res%pedigree(tmpId)%damPointer%originalId)
				if (tmpSire == DICT_NULL .or. tmpDam == DICT_NULL) then
					print *, "WE SHOULD NOT GET HERE IN COPY! PLEASE CONTACT DEVELOPERS"
					print *,res%pedigree(tmpId)%sireId,res%dictionary%getValue(res%pedigree(tmpId)%sireId)
				else
					call res%pedigree(tmpSire)%addOffspring(res%pedigree(tmpId))
					res%pedigree(tmpId)%sirePointer =>  res%pedigree(tmpSire)
					if (res%pedigree(tmpSire)%nOffs== 1) then
						call res%sireList%list_add(res%pedigree(tmpSire))
					endif
					call res%pedigree(tmpDam)%addOffspring(res%pedigree(tmpId))
					res%pedigree(tmpId)%damPointer =>  res%pedigree(tmpDam)
					if (res%pedigree(tmpDam)%nOffs== 1) then
						call res%damList%list_add(res%pedigree(tmpDam))
					endif
				endif
			enddo

		end subroutine copyPedigree

		!---------------------------------------------------------------------------
		!< @brief Checks if all of the pointers line up in the pedigree
		!< details this is done with the LOC intrinsic function
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		logical function deepCheckPedigree(this)

			class(pedigreeHolder) , intent(in) :: this
			integer :: i, h
			type(IndividualLinkedListNode), pointer :: p1
			deepCheckPedigree = .true.
			do i=1, this%pedigreeSize

				if (associated(this%pedigree(i)%sirePointer)) then
					if ( .not. this%pedigree(i)%sirePointer%isUnknownDummy .and. .not. this%pedigree(i)%mendelianError(1)) then
					if (loc(this%pedigree(i)%sirePointer) /= loc(this%pedigree(this%dictionary%getvalue(this%pedigree(i)%sireId))) ) then
						deepCheckPedigree = .false.
						write(error_unit, *) "WARNING: Sire pointer is out of alignment on ind:", this%pedigree(i)%originalId,"  sire: ",this%pedigree(i)%sireId
					endif
					endif
				endif

				if (associated(this%pedigree(i)%damPointer)) then

					if (.not. this%pedigree(i)%damPointer%isUnknownDummy .and. .not. this%pedigree(i)%mendelianError(2)) then
					if (loc(this%pedigree(i)%damPointer) /= loc(this%pedigree(this%dictionary%getvalue(this%pedigree(i)%damId)))) then
						deepCheckPedigree = .false.
						write(error_unit, *) "WARNING: dam pointer is out of alignment on ind:", this%pedigree(i)%originalId,"  dam: ",this%pedigree(i)%damId
					endif
					endif
				endif

				do h=1,this%pedigree(i)%nOffs

					if (loc(this%pedigree(i)%offsprings(h)%p) /= loc(this%pedigree(this%dictionary%getvalue(this%pedigree(i)%offsprings(h)%p%originalId)))) then
						deepCheckPedigree = .false.
						write(error_unit, *) "WARNING: offspring pointer is out of alignment on ind:", this%pedigree(i)%originalId,"  offspring: ",this%pedigree(i)%offsprings(h)%p%originalId
					endif

				enddo

				if (this%pedigree(i)%isDummy .and. this%pedigree(i)%nOffs == 0) then
					write(error_unit, *) "WARNING: Dummy animal does not have any kids attached: "
					deepCheckPedigree = .false.
				endif
			enddo

			p1 => this%sireList%first
			do i=1, this%sireList%length
				if (loc(p1%item) /= loc(this%pedigree(this%dictionary%getvalue(p1%item%originalId)))) then
					deepCheckPedigree = .false.
					write(error_unit, *) "WARNING: damlist is wrong is out of alignment on ind:", p1%item%originalId
					return
				endif
				p1 => p1%next
			enddo

			p1 => this%damList%first
			do i=1, this%damList%length
				if (loc(p1%item) /= loc(this%pedigree(this%dictionary%getvalue(p1%item%originalId)))) then
					deepCheckPedigree = .false.
					write(error_unit, *) "WARNING: damlist is wrong is out of alignment on ind:", p1%item%originalId
				endif
				p1 => p1%next
			enddo
		end function deepCheckPedigree

		!---------------------------------------------------------------------------
		!< @brief Checks for array equality
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		logical function equality(pedOne, pedTwo)

			type(pedigreeHolder) , intent(in) :: pedOne, pedTwo
			integer :: i,h

			equality = .true.

			if (pedOne%pedigreeSize /=  pedTwo%pedigreeSize) then
				equality = .false.
				return
			endif

			if (pedOne%maxPedigreeSize /=  pedTwo%maxPedigreeSize) then
				equality = .false.
				return
			endif

			if (allocated(pedOne%genotypeDictionary)) then
				if (.not. allocated(pedTwo%genotypeDictionary)) then
					equality = .false.
					return
				endif
				if (pedOne%genotypeDictionary == pedTwo%genotypeDictionary) then
					equality = .false.
					return
				endif
			endif

			if (allocated(pedOne%hdDictionary)) then
				if (.not. allocated(pedTwo%hdDictionary)) then
					equality = .false.
					return
				endif
				if (pedOne%hdDictionary == pedTwo%hdDictionary) then
					equality = .false.
					return
				endif
			endif

			if (pedOne%nGenotyped /=  pedTwo%nGenotyped) then
				equality = .false.
				return
			endif

			if (pedOne%unknownDummys /=  pedTwo%unknownDummys) then
				equality = .false.
				return
			endif

			if (pedOne%addedRealAnimals /=  pedTwo%addedRealAnimals) then
				equality = .false.
				return
			endif

			if (pedOne%sireList%length /=  pedTwo%sireList%length) then
				equality = .false.
				return
			endif

			if (pedOne%damList%length /=  pedTwo%damList%length) then
				equality = .false.
				return
			endif

			if (pedOne%founders%length /=  pedTwo%founders%length) then
				equality = .false.
				return
			endif
			! Compare lists
			block
				type(individualLinkedListnode), pointer :: p1,p2

				p1 => pedOne%sireList%first
				p2 => pedTwo%sireList%first
				do i =1, pedOne%sireList%length
					if (.not. compareIndividual(p1%item, p2%item)) then
						equality = .false.
						return
					endif
					p1 => p1%next
					p2 => p2%next
				enddo

				p1 => pedOne%damList%first
				p2 => pedTwo%damList%first
				do i =1, pedOne%damList%length
					if (.not. compareIndividual(p1%item, p2%item)) then
						equality = .false.
						return
					endif
					p1 => p1%next
					p2 => p2%next
				enddo
				p1 => pedOne%founders%first
				p2 => pedTwo%founders%first
				do i =1, pedOne%founders%length
					if (.not. compareIndividual(p1%item, p2%item)) then
						equality = .false.
						return
					endif
					p1 => p1%next
					p2 => p2%next
				enddo

				if (pedOne%maxGeneration /= pedTwo%maxGeneration) then
					equality = .false.
					return
				endif

				do i=1, pedOne%maxGeneration
					if (pedOne%generations(i)%length /=pedTwo%generations(i)%length) then
						equality = .false.
						return
					endif

					p1 => pedOne%generations(i)%first
					p2 => pedTwo%generations(i)%first

					do h=1,pedOne%generations(i)%length

						if (.not. compareIndividual(p1%item, p2%item)) then
							equality = .false.
							return
						endif

						p1 => p1%next
						p2 => p2%next
					enddo
				enddo
			end block

			if (pedOne%nHd /=  pedTwo%nHd) then
				equality = .false.
				return
			endif

			! compare individuals
			do i=1, pedOne%pedigreeSize
				if (.not. (pedOne%pedigree(i) == pedTwo%pedigree(i))) then
					equality = .false.
					return
				endif
			end do

			if (pedOne%founders%length /= pedTwo%founders%length) then
				equality = .false.
				return
			endif

		end function equality



		!---------------------------------------------------------------------------
		!< @brief  writes out phase probabilites
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutPhaseProbabilities(this,filename, indexesToPrint)
			character(len=*),intent(in)  :: filename
			class(PedigreeHolder) :: this
			integer, dimension(:),optional, intent(in) :: indexesToPrint
			character(len=12) :: StrSnp,OutFmt
			integer :: fileUnit,i
			open(newunit=FileUnit, file=filename, status="replace")

			write(StrSnp,*) size(this%pedigree(1)%phaseProbabilities(1,:))
			! OutFmt='(i20,'//trim(adjustl(StrSnp))//'f7.1)'

			OutFmt='(a,'//trim(adjustl(StrSnp))//'f7.1)'
			if (present(indexesToPrint)) then

				do i=1,size(indexesToPrint)
					write(FileUnit, OutFmt) this%pedigree(indexesToPrint(i))%originalId, this%pedigree(indexesToPrint(i))%phaseProbabilities(1,:)
					write(FileUnit, OutFmt) this%pedigree(indexesToPrint(i))%originalId, this%pedigree(indexesToPrint(i))%phaseProbabilities(2,:)
				enddo
			else
				do i=1,this%pedigreeSize
					write(FileUnit, OutFmt) this%pedigree(i)%originalId, this%pedigree(i)%phaseProbabilities(1,:)
					write(FileUnit, OutFmt) this%pedigree(i)%originalId, this%pedigree(i)%phaseProbabilities(2,:)
				enddo
			endif

			close(FileUnit)

		end subroutine writeOutPhaseProbabilities


		!---------------------------------------------------------------------------
		!< @brief  writes out geno probabilites
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypeProbabilities(this,filename, indexesToPrint)
			character(len=*),intent(in)  :: filename
			class(PedigreeHolder) :: this
			integer, dimension(:),optional, intent(in) :: indexesToPrint
			integer :: i, fileUnit
			character(len=30) :: StrSnp,OutFmt

			open(newunit=FileUnit, file=filename, status="replace")

			write(StrSnp,*) size(this%pedigree(1)%genotypeProbabilities)
			! OutFmt='(i20,'//trim(adjustl(StrSnp))//'f7.1)'
			OutFmt='(a,'//trim(adjustl(StrSnp))//'f7.1)'

			if (present(indexesToPrint)) then

				do i=1,size(indexesToPrint)
					write(FileUnit, OutFmt) this%pedigree(indexesToPrint(i))%originalId, this%pedigree(indexesToPrint(i))%genotypeProbabilities
				enddo
			else
				do i=1,this%pedigreeSize
					write(FileUnit, OutFmt) this%pedigree(i)%originalId, this%pedigree(i)%genotypeProbabilities
				enddo
			endif

			close(FileUnit)


		end subroutine writeOutGenotypeProbabilities

		!---------------------------------------------------------------------------
		!< @brief constructor for creating an empty pedgiree
		!< @details Constructs an empty pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initEmptyPedigree(pedStructure,nsnps, minSize)
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure
			integer, optional :: nsnps
			integer, optional :: minSize !< pedigree will be this size * 4

			allocate(pedStructure%dictionary)
			call pedStructure%dictionary%DictStructure()

			pedStructure%pedigreeSize = 0
			pedStructure%nDummys = 0
			pedStructure%nGenotyped = 0
			pedStructure%nHd = 0
			pedStructure%addedRealAnimals = 0

			if (present(minSize)) then
				pedStructure%maxPedigreeSize = minSize
			else
				pedStructure%maxPedigreeSize = DEFAULTDICTSIZE
			endif
			allocate(pedStructure%pedigree(pedStructure%maxPedigreeSize))
			allocate(pedStructure%founders)
			allocate(pedStructure%sireList)
			allocate(pedStructure%damList)
			pedStructure%nsnpsPopulation = 0
			if (present(nsnps)) then
				pedStructure%nsnpsPopulation = nsnps
			endif
		end subroutine initEmptyPedigree




		!---------------------------------------------------------------------------
		!< @brief Returns an HD pedigree
		!< details Creates a new pedigree Object with only the HD animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2018
		!---------------------------------------------------------------------------
		subroutine getHDPedigree(this,new)
			class(pedigreeHolder) :: this
			type(pedigreeHolder) :: new
			integer :: i, sire,dam,tmpAnimalCount, indiv


			integer, allocatable, dimension(:) :: tmpAnimalArray

			! new%pedigreeSize = this%nHd
			call initEmptyPedigree(new,this%nsnpsPopulation, this%pedigreeSize)
			tmpAnimalCount = 0
			new%nhd = 0
			new%nDummys = 0

			allocate(new%hdMap(this%nHd))
			allocate(new%genotypeMap(this%nGenotyped))
			allocate(tmpAnimalArray(this%nHd))
			allocate(new%hdDictionary)
			allocate(new%genotypeDictionary)
			call new%hdDictionary%DictStructure()
			call new%genotypeDictionary%DictStructure()

			if (this%nHd == 0) then
				write(error_unit, * ) "ERROR: No animals defined as HD."
				return
			endif

			do i =1,this%nHd
				sire =new%dictionary%getValue(this%pedigree(this%hdMap(i))%sireId)
				dam = new%dictionary%getValue(this%pedigree(this%hdMap(i))%damId)
				if (this%pedigree(this%hdMap(i))%founder .or.(dam /=DICT_NULL .and. sire /=DICT_NULL) ) then
					new%pedigreeSize = new%pedigreeSize+1
					new%nhd = new%nhd + 1
					new%pedigree(new%pedigreeSize) = this%pedigree(this%hdMap(i))
					call copyIndividual(new%pedigree(new%pedigreeSize),this%pedigree(this%hdMap(i)))
					new%pedigree(i)%id = new%pedigreeSize
					new%hdMap(new%nhd) = new%pedigreeSize
					new%genotypeMap(new%nhd) = new%pedigreeSize
					call new%hdDictionary%addKey(this%pedigree(this%hdMap(i))%originalId, new%nhd)
					call new%genotypeDictionary%addKey(this%pedigree(this%hdMap(i))%originalId, new%nhd)
					call new%dictionary%addKey( this%pedigree(this%hdMap(i))%originalId,new%pedigreeSize)
					call new%pedigree(i)%resetOffspringInformation
					if (sire /= DICT_NULL) then
						call new%pedigree(sire)%AddOffspring(new%pedigree(i))
						! from above condition we know dam is also not null
						call new%pedigree(dam)%AddOffspring(new%pedigree(i))
					endif

				else
					tmpAnimalCount = tmpAnimalCount + 1
					tmpAnimalArray(tmpAnimalCount) = this%hdmap(i)
				endif
			enddo

			do i=1,tmpAnimalCount

			! check if sire and dam are in the new HD pedigree
				sire =new%dictionary%getValue(this%pedigree(tmpAnimalArray(i))%sireId)
				dam = new%dictionary%getValue(this%pedigree(tmpAnimalArray(i))%damId)
				
				new%pedigreeSize = new%pedigreeSize+1
				call new%dictionary%addKey( this%pedigree(tmpAnimalArray(i))%originalId,new%pedigreeSize)
				new%nhd = new%nhd + 1
				new%hdMap(new%nhd) = new%pedigreeSize
				call new%hdDictionary%addKey(this%pedigree(tmpAnimalArray(i))%originalId, new%nhd)
				indiv = new%pedigreeSize
				call copyIndividual(new%pedigree(indiv),this%pedigree(tmpAnimalArray(i)))

				call new%pedigree(i)%resetOffspringInformation

				! check if we need to create new sire/dam dummies
				if (dam ==DICT_NULL) then
					call new%createDummyAnimalAtEndOfPedigree(dam)
					call new%pedigree(dam)%addOffspring(new%pedigree(indiv))
					new%pedigree(indiv)%damId = new%pedigree(dam)%originalId

				else
					call new%pedigree(dam)%addOffspring(new%pedigree(indiv))
				endif

				if (sire ==DICT_NULL) then

					call new%createDummyAnimalAtEndOfPedigree(sire)
					call new%pedigree(sire)%addOffspring(new%pedigree(indiv))
					new%pedigree(indiv)%sireId = new%pedigree(sire)%originalId
				else
					call new%pedigree(sire)%addOffspring(new%pedigree(indiv))
				endif

			enddo
			new%nGenotyped = this%nHd
			deallocate(tmpAnimalArray)
		end subroutine getHDPedigree


		!---------------------------------------------------------------------------
		!< @brief Function to find and correct mendellian inconsistencies in the pedigree
		!< details this is done with the LOC intrinsic function
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		function findMendelianInconsistencies(ped, thresholdIn,file,snpfilePath) result(CountChanges)
			use genotypeModule
			use BitUtilities
			implicit none
			class(pedigreeHolder) :: ped
			type(Mendelian), dimension(:), allocatable :: mend
			real, intent(in), optional :: thresholdIn !< threshold for removing animals
			real :: threshold
			integer :: sireInconsistencies, damInconsistencies
			integer :: outfile,snpFile
			type(individual),pointer :: sire, dam
			integer :: i,j
			character(len=*),intent(in), optional :: file !< path to
			character(len=*),intent(in), optional :: snpfilePath !< path for file to output changes that were made to individual snps.
			integer :: CountChanges !< returns the changes to the pedigree that the function has done
			integer :: dumId
			integer :: snpChanges
			logical :: sireRemoved, damRemoved

			sireRemoved = .false.
			damRemoved = .false.


			if (present(thresholdIn)) then
				threshold = thresholdIn
			else
				threshold = 0.05
			endif
			snpChanges = 0
			if (present(snpfilePath)) then
				open(newunit=snpFile, file=snpfilePath, status='unknown')
				write (snpFile,'(3a30)') "AnimalId","ParentId", "SNP Position of Mismatch"
			endif
			if (present(file)) then
				open(newunit=outfile, file=file, status='unknown')
				write (outfile,'(4a30)') "AnimalId","Sire Or Dam", "ParentId","Number Of Inconsistencies"
			endif

			CountChanges = 0

			allocate(mend(ped%pedigreeSize))

			if (ped%isSorted == 0) then
				write(error_unit, *) "WARNING - mendelian Inconsistency checks are being run without the pedigree being sorted. This could have weird effects"
			endif

			do i=1,ped%pedigreeSize
				! if sire is associated, then dam must be too
				if (associated(ped%pedigree(i)%sirePointer)) then

					sire => ped%pedigree(i)%sirePointer
					dam => ped%pedigree(i)%damPointer

					mend(i) = ped%pedigree(i)%individualGenotype%mendelianInconsistencies(sire%individualGenotype,dam%individualGenotype)

					sireInconsistencies = bitcount(mend(i)%paternalInconsistent)

					! update inconsistencies, maybe do this first
					ped%pedigree(i)%sirePointer%inconsistencyCount = ped%pedigree(i)%sirePointer%inconsistencyCount + sireInconsistencies
					damInconsistencies = bitcount(mend(i)%maternalInconsistent)
					ped%pedigree(i)%damPointer%inconsistencyCount = ped%pedigree(i)%damPointer%inconsistencyCount + damInconsistencies

					if ((float(sireInconsistencies) / ped%pedigree(i)%individualGenotype%length) > threshold) then
						! remove sire link
						if (present(file)) then
							write (outfile,'(3a30,1I)') Ped%pedigree(i)%originalID, "SIRE", Ped%pedigree(i)%sirePointer%originalID, sireInconsistencies
						endif
						CountChanges=CountChanges+1
						! remove offspring link
						ped%pedigree(i)%mendelianError(1) = .true.
						if (ped%pedigree(i)%sirePointer%nOffs ==1) then
							call ped%sireList%list_remove(ped%pedigree(i)%sirePointer)
						endif
						call ped%pedigree(i)%sirePointer%removeOffspring(ped%pedigree(i))

						
						call ped%createDummyAnimalAtEndOfPedigree(dumId, i)

						sireRemoved = .true.
					endif

					if ((float(damInconsistencies) / ped%pedigree(i)%individualGenotype%length) > threshold) then

						! remove dam link
						if (present(file)) then
							write (outfile,'(3a30,1I)') &
							Ped%pedigree(i)%originalID, "DAM", Ped%pedigree(i)%damPointer%originalID, damInconsistencies
						endif
						CountChanges=CountChanges+1
						! remove offspring link
						ped%pedigree(i)%mendelianError(2) = .true.
						if (ped%pedigree(i)%damPointer%nOffs ==1) then
							call ped%damlist%list_remove(ped%pedigree(i)%damPointer)
						endif
						call ped%pedigree(i)%damPointer%removeOffspring(ped%pedigree(i))
						
						call ped%createDummyAnimalAtEndOfPedigree(dumId, i)
						damRemoved =.true.
					endif

					! means we only Calculate inconsistencies for animals that aren't unlinked
					if (.not. sireRemoved) then
						do j=1,ped%pedigree(i)%individualGenotype%length

							if (testBit(mend(i)%paternalInconsistent,j)) then
								ped%pedigree(i)%sirePointer%inconsistencies(j) = ped%pedigree(i)%sirePointer%inconsistencies(j) + 1
								ped%pedigree(i)%inconsistencies(j) = ped%pedigree(i)%inconsistencies(j) + 1
							else if (testBit(mend(i)%paternalconsistent,j))  then
								ped%pedigree(i)%sirePointer%inconsistencies(j) = ped%pedigree(i)%sirePointer%inconsistencies(j) - 1
								ped%pedigree(i)%inconsistencies(j) = ped%pedigree(i)%inconsistencies(j) - 1

							endif
						enddo
					endif

					if (.not. damRemoved) then
						do j=1,ped%pedigree(i)%individualGenotype%length

							if (testBit(mend(i)%maternalInconsistent,j)) then
								ped%pedigree(i)%damPointer%inconsistencies(j) = ped%pedigree(i)%damPointer%inconsistencies(j) + 1
								ped%pedigree(i)%inconsistencies(j) = ped%pedigree(i)%inconsistencies(j) + 1
							else if (testBit(mend(i)%maternalconsistent,j))  then
								ped%pedigree(i)%damPointer%inconsistencies(j) = ped%pedigree(i)%damPointer%inconsistencies(j) - 1
								ped%pedigree(i)%inconsistencies(j) = ped%pedigree(i)%inconsistencies(j) - 1

							endif

						enddo
					endif
				endif

			enddo

			! at this point, we have calculated inconsistenceis

			do i=1, ped%pedigreeSize
				if (ped%pedigree(i)%Founder) cycle

				! if both parents haven't been removed, check most likely one
				! call ped%pedigree(i)%individualGenotype%setMissingBits(mend%individualInconsistencies)
				do j=1,ped%pedigree(i)%individualGenotype%length

					!< if either is a dummy, likely that individualInconsistent is inccorrect
					if (.not. ped%pedigree(i)%sirePointer%isDummy .and. .not. ped%pedigree(i)%damPointer%isDummy) then
						! if both were inconsistent, set to which one is more likely
						if (testBit(mend(i)%individualInconsistent,j) ) then
							snpChanges = snpChanges + 1


							! if dam is more likely correct, set to dam value, and set sire to missing
							if (ped%pedigree(i)%sirePointer%inconsistencies(j) > ped%pedigree(i)%damPointer%inconsistencies(j)) then
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%sirePointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%sirePointer%originalID, j
								endif
							else if(ped%pedigree(i)%sirePointer%inconsistencies(j) < ped%pedigree(i)%damPointer%inconsistencies(j)) then
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%damPointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%damPointer%originalID, j
								endif
							else !< if they are both as unlikely, set all the animals as missing here
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%sirePointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%damPointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%sirePointer%originalID, j
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%damPointer%originalID, j
								endif
							endif
						endif

					else !if animal has a dummy parent

						if (.not. ped%pedigree(i)%sirePointer%isDummy  .and.  testBit(mend(i)%paternalInconsistent,j) ) then
							snpChanges = snpChanges +1
							if (ped%pedigree(i)%inconsistencies(j) > ped%pedigree(i)%sirePointer%inconsistencies(j)) then
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(1a30,I)') & Ped%pedigree(i)%originalID, j
								endif
							else if (ped%pedigree(i)%inconsistencies(j) > ped%pedigree(i)%sirePointer%inconsistencies(j)) then
								call ped%pedigree(i)%sirePointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(1a30,I)') & Ped%pedigree(i)%sirePointer%originalID, j
								endif
							else
								call ped%pedigree(i)%sirePointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%sirePointer%originalID, j
								endif
							endif
						endif

						if (.not. ped%pedigree(i)%damPointer%isDummy  .and.  testBit(mend(i)%maternalInconsistent,j) )  then
							snpChanges = snpChanges +1
							if (ped%pedigree(i)%inconsistencies(j) > ped%pedigree(i)%damPointer%inconsistencies(j)) then
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(1a30,I)') & Ped%pedigree(i)%originalID, j
								endif
							else if (ped%pedigree(i)%inconsistencies(j) < ped%pedigree(i)%damPointer%inconsistencies(j)) then
								call ped%pedigree(i)%dampointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(1a30,I)') & Ped%pedigree(i)%damPointer%originalID, j
								endif
							else
								call ped%pedigree(i)%damPointer%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								call ped%pedigree(i)%individualGenotype%setGenotype(j,MISSINGGENOTYPECODE)
								if (present(snpfilePath)) then
									write (snpFile,'(2a30,I)') & Ped%pedigree(i)%originalID, Ped%pedigree(i)%damPointer%originalID, j
								endif
							endif
						endif
					endif

				enddo
			enddo

			deallocate(mend)

			if (present(snpfilePath)) then
				close(snpFile)
			endif

			if (present(file)) then
				write (outfile,*) CountChanges," changes were made to the pedigree"
				close (outfile)
			endif
			print*, " ",CountChanges," errors in the pedigree due to Mendelian inconsistencies"
			print*, " ",snpChanges," snps changed across individuals"

		end function findMendelianInconsistencies


		!---------------------------------------------------------------------------
		!< @brief Sets phase information from an array
		!< @details sets phase objects for animal when given an array
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setPhaseFromArray(this, array)

			use HaplotypeModule

			class(PedigreeHolder), intent(inout) :: this
			integer(kind=1),dimension(:,:,:),intent(in) :: array !< array should be of format (recodedindId, snp, allele )
			integer :: i

			do i=1, size(array,1)
				call this%pedigree(i)%individualPhase(1)%newHaplotypeInt(array(i,:,1))
				call this%pedigree(i)%individualPhase(2)%newHaplotypeInt(array(i,:,2))
			enddo

		end subroutine setPhaseFromArray



		!---------------------------------------------------------------------------
		!< @brief Returns IndividualLinkedList of animals that are parents (with no overlaps)
		!< @details caches, so avoids code duplication
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getUniqueParents(this) result(res)

			class(PedigreeHolder) :: this

			type(IndividualLinkedList) :: res
			type(DictStructure) ::tmpDict
			integer :: i
			type(IndividualLinkedListNode), pointer :: node


			if (allocated(this%uniqueParentList)) then
				res= this%uniqueParentList
			else

				call tmpDict%DictStructure()
				node => this%sireList%first

				do i=1,this%sireList%length
					call tmpDict%addKey(node%item%originalID,1)
					call res%list_add(node%item)
					node => node%next
				enddo

				node => this%damList%first
				do i=1, this%damList%length

					if (tmpDict%getValue(node%item%originalID) == DICT_NULL) then
						call res%list_add(node%item)
					endif
					node => node%next
				enddo
				allocate(this%uniqueParentList)
				this%uniqueParentList = res
			endif
		end function getUniqueParents


		!---------------------------------------------------------------------------
		!< @brief Sets the phase for homozygotic snps
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		subroutine homozygoticFillIn(this)

			class(PedigreeHolder), intent(inout) :: this
			! Impute phase information for homozygous cases

			integer :: i,j
			!$OMP PARALLEL DO &
			!$OMP PRIVATE(i,j)
			do i=1, this%pedigreeSize
				do j=1, this%pedigree(this%genotypeMap(1))%individualGenotype%length
					if (this%pedigree(i)%individualGenotype%getGenotype(j) == 2) then
						call this%pedigree(i)%individualPhase(1)%setPhase(j,1)
						call this%pedigree(i)%individualPhase(2)%setPhase(j,1)
					else if (this%pedigree(i)%individualGenotype%getGenotype(j) == 0) then
						call this%pedigree(i)%individualPhase(1)%setPhase(j,0)
						call this%pedigree(i)%individualPhase(2)%setPhase(j,0)
					endif
				enddo
			enddo
			!$OMP END PARALLEL DO

		end subroutine homozygoticFillIn

		!---------------------------------------------------------------------------
		!< @brief Calculate correlation between pedigree
		!<Calculates the correlation between individuals in the pedigree, assuming no inbreeding.   Will also handle groups - just give
		!<the group parents a unique ID.   It assumes that the pedigree has been sorted such that the unknown dummies are at the end
		!<(sortPedigreeAndOverwrite(1)).
		!< @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
		!---------------------------------------------------------------------------
		function calculatePedigreeCorrelationNoInbreeding(pedigreeIn, additiveVarianceIn) result(valuesOut)
			class(PedigreeHolder), intent(in):: pedigreeIn
			real(real64), dimension(:,:), allocatable:: valuesOut
			real(real64), intent(in), optional:: additiveVarianceIn
			integer:: damKnown, sireKnown
			integer::numLevels, ID, sireID, damID
			integer:: i
			real(real64):: D

			numLevels = pedigreeIn%pedigreeSize-pedigreeIn%UnknownDummys

			allocate(valuesOut(numLevels, numLevels))
			valuesOut = 0_real64

			do i = 1,numLevels
				damKnown = 0
				sireKnown = 0
				damID=0
				sireID=0
				ID = pedigreeIn%pedigree(i)%ID
				if (associated(pedigreeIn%pedigree(i)%sirePointer)) then
					if (.not. pedigreeIn%pedigree(i)%sirePointer%isUnknownDummy) then
						sireKnown = 1
						sireID = pedigreeIn%pedigree(i)%sirePointer%ID
					end if
				end if

				if (associated(pedigreeIn%pedigree(i)%damPointer)) then
					if (.not. pedigreeIn%pedigree(i)%damPointer%isUnknownDummy) then
						damKnown = 1
						damID = pedigreeIn%pedigree(i)%damPointer%ID
					end if
				end if

				D = 4.0_real64/(2.0_real64+damKnown+sireKnown) !If both are known, then this is set to be 2,if only 1 is known this is 4.0/3.0

				valuesOut(ID, ID) = valuesOut(ID,ID)+D

				if (sireKnown .ne. 0) then
					valuesOut(sireID, ID) = valuesOut(sireID, ID) -0.5*D
					valuesOut(ID, sireID) = valuesOut(ID, sireID) - 0.5*D
					valuesOut(sireID, sireID) = valuesOut(sireId, sireID) +0.25*D
				end if

				if (damKnown .ne. 0) then
					valuesOut(damID, ID) = valuesOut(damID, ID) -0.5*D
					valuesOut(ID, damID) = valuesOut(ID, damID) -0.5*D
					valuesOut(damID, damID) = valuesOut(damId, damID) +0.25*D
				end if

				if ((sireKnown .ne. 0) .and. (damKnown.ne.0)) then
					valuesOut(sireID, damID) = valuesOut(sireID, damID) + 0.25*D
					valuesOut(damID, sireID) = valuesOut(damID, sireID) +0.25*D
				end if
			end do

			if (present(additiveVarianceIn)) then
				valuesOut = valuesOut*additiveVarianceIn
			end if

		end function calculatePedigreeCorrelationNoInbreeding



		!---------------------------------------------------------------------------
		!< @brief Sets pedigree genotype from integer array
		!< @details Sets animals genotype from integer array
		!< NOTE: does not set animals as genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setGenotypeFromArray(this, array)
			use GenotypeModule

			class(PedigreeHolder), intent(inout)  :: this
			integer(kind=1),dimension(:,:) :: array !< array should be of format (recodedindId, snp)
			integer :: i

			do i=1, size(array,1)
				call this%pedigree(i)%individualGenotype%newGenotypeInt(array(i,:))
			enddo

		end subroutine setGenotypeFromArray


		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class
		!< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigree(pedStructure, fileIn, numberInFile, genderFile, nsnps,dontInitAll)
			use AlphaHouseMod, only : countLines
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure
			character(len=*),intent(in) :: fileIn !< path of pedigree file
			character(len=*), intent(in),optional :: genderFile !< path to gender file
			integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
			integer, optional, intent(in) :: nsnps !< number of snps for the population
			integer, optional :: dontInitAll !< don't initialise all animals

			character(len=IDLENGTH) :: tmpId,tmpSire,tmpDam
			integer(kind=int32) :: stat, fileUnit,tmpSireNum, tmpDamNum, tmpGender,tmpIdNum
			integer(kind=int64) :: nIndividuals
			integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
			integer :: tmpAnimalArrayCount,i
			integer(kind=int64) :: sizeDict
			logical :: sireFound, damFound

			pedStructure%isSorted = 0
			pedStructure%nDummys = 0
			pedStructure%nHd = 0
			tmpAnimalArrayCount = 0
			pedStructure%nGenotyped = 0
			pedStructure%nsnpsPopulation = 0

			call destroyPedigree(pedStructure)

			if (present(nsnps)) then
				pedStructure%nsnpsPopulation = nsnps
			endif
			if (present(numberInFile)) then
				nIndividuals = numberInFile
			else
				nIndividuals = countLines(fileIn)
			endif

			pedStructure%addedRealAnimals = nIndividuals

			sizeDict = nIndividuals
			pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
			allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
			pedStructure%pedigreeSize = nIndividuals

			allocate(pedStructure%dictionary)
			call pedStructure%dictionary%DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder

			allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
			allocate(pedStructure%inputMap(pedStructure%maxPedigreeSize))

			pedStructure%inputMap = 0

			pedStructure%maxGeneration = 0
			open(newUnit=fileUnit, file=fileIn, status="old")


			allocate(pedStructure%sireList)
			allocate(pedStructure%damList)
			allocate(pedStructure%Founders)
			do i=1,nIndividuals
				sireFound = .false.
				damFound = .false.
				read(fileUnit,*) tmpId,tmpSire,tmpDam

				if (trim(tmpId) == trim(tmpsire) .or. trim(tmpDam) == trim(tmpId)) then

					write(error_unit,*) "Error: Animal ", trim(tmpId), " has been given itself as a parent. please fix this. Error on line:",i
					stop
				endif
				call pedStructure%dictionary%addKey(tmpId, i)

				if (present(dontInitAll)) then
					call pedStructure%Pedigree(i)%initIndividual(trim(tmpId),trim(tmpSire),trim(tmpDam), i) !Make a new individual based on info from ped
				else
					call pedStructure%Pedigree(i)%initIndividual(trim(tmpId),trim(tmpSire),trim(tmpDam), i, nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped

				endif
				pedStructure%Pedigree(i)%originalPosition = i
				pedStructure%inputMap(i) = i
				if (tmpSire /= EMPTY_PARENT) then !check sire is defined in pedigree
					tmpSireNum = pedStructure%dictionary%getValue(tmpSire)
					if (tmpSireNum /= DICT_NULL) then
						sireFound = .true.
					endif
				endif

				if (tmpDam /= EMPTY_PARENT) then
					tmpDamNum = pedStructure%dictionary%getValue(tmpDam)
					if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
						damFound = .true.
					endif
				endif
				if (tmpSire == EMPTY_PARENT .and. tmpDam == EMPTY_PARENT) then !if animal is a founder
					pedStructure%Pedigree(i)%founder = .true.
					call pedStructure%Founders%list_add(pedStructure%Pedigree(i))
				else if (sireFound == .false. .or. damFound == .false. ) then
					tmpAnimalArrayCount = tmpAnimalArrayCount +1
					tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
				else ! if sire and dam are both found
					pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
					call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))

					if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
						call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum))
					endif
					pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
					call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))
					if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
						call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum))
					endif
				endif
			enddo

			close(fileUnit)
			! if we want gender info read in rather than calculated on the fly, lets do it here
			if (present(genderFile)) then !read in gender here
				open(newUnit=fileUnit, file=genderFile, status="old")
				do
					read (fileUnit,*, IOSTAT=stat) tmpId,tmpGender
					if (stat /=0) exit
					tmpIdNum = pedStructure%dictionary%getValue(tmpId)
					if (tmpIdNum /= DICT_NULL) then
						pedStructure%Pedigree(i)%gender = int(tmpGender)
					else
						write(error_unit, *) "WARNING: Gender  defined for an animal that does not exist in Pedigree!"
						write(error_unit, *) "Amimal:",tmpId
					endif
				end do
			endif
			call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)
			deallocate(tmpAnimalArray)
			!  write(output_unit, *) "Number of animals in Pedigree:",pedStructure%pedigreeSize-pedStructure%nDummys
			!  write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys
		end  subroutine


		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class taking in files that have been written out
		!< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreeFromOutputFile(pedStructure,pedigreeFile, genderFile, pedGenotypeFile, phaseFile, nsnps)
			use AlphaHouseMod, only : countColumns

			character(len=*), intent(in) :: pedigreeFile
			character(len=*), intent(in), optional :: genderFile, pedGenotypeFile, phaseFile
			integer, intent(inout) :: nsnps
			type(PedigreeHolder) :: pedStructure

			call destroyPedigree(pedStructure)
			allocate(pedStructure%sireList)
			allocate(pedStructure%damList)
			allocate(pedStructure%Founders)
			if (nsnps == 0 .and. present(pedGenotypeFile)) then
				nsnps = countColumns(pedGenotypeFile,' ') - 1
			endif
			if (present(genderFile)) then
				call initPedigree(pedStructure,pedigreeFile, nsnps=nsnps, genderFile=genderFile)
			else
				call initPedigree(pedStructure,pedigreeFile, nsnps=nsnps)

			endif

			if (present(pedGenotypeFile)) then
				call pedStructure%addGenotypeInformationFromFile(pedGenotypeFile, nsnps)
			endif

			if (present(phaseFile)) then
				call pedStructure%addPhaseInformationFromFile(phaseFile, nsnps)
			endif


		end subroutine initPedigreeFromOutputFile


		!---------------------------------------------------------------------------
		!< @brief creates a pedigree object from input folder
		!< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreeFromOutputFileFolder(pedStructure,folder, nsnps)
			use AlphaHouseMod, only : countColumns
			character(len=*), intent(in) :: folder !< input folder
			integer, intent(inout) :: nsnps
			class(PedigreeHolder) :: pedStructure
			character(len=:), allocatable :: pedigreeFile, genotypeFile, phaseFile,genderFile

			call destroyPedigree(pedStructure)

			pedigreeFile = "pedigree.txt"
			genotypeFile = "genotypes.txt"
			phaseFile = "phase.txt"
			genderFile = "gender.txt"

			if (nsnps == 0) then
				nsnps = countColumns(genotypeFile,' ') -1
			endif
			call initPedigree(pedStructure, folder//pedigreeFile, nsnps=nsnps, genderFile=folder//genderFile)
			call pedStructure%addGenotypeInformationFromFile(folder//genotypeFile, nsnps)
			call pedStructure%addPhaseInformationFromFile(folder//phaseFile, nsnps)


		end subroutine initPedigreeFromOutputFileFolder


		!---------------------------------------------------------------------------
		!< @brief writes out the pedigree to predefined filess
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!---------------------------------------------------------------------------
		subroutine writeOutPedigree(this,pedigreeFolder)
			class(PedigreeHolder) :: this
			character(len=*), intent(in), optional :: pedigreeFolder
			character(len=:), allocatable :: pedigreeFile, genotypeFile, phaseFile,genderFile


			pedigreeFile = "pedigree.txt"
			genotypeFile = "genotypes.txt"
			phaseFile = "phase.txt"
			genderFile = "gender.txt"


			call this%printPedigreeOriginalFormat(pedigreeFolder//pedigreeFile)
			call this%WriteoutPhase(pedigreeFolder//phaseFile)
			call this%writeOutGenotypes(pedigreeFolder//genotypeFile) !< just output genotype information for animals that are set to genotpyyed
			call this%writeOutGenders(pedigreeFolder//genderFile)

		end subroutine writeOutPedigree




		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class using Genotype File format
		!< @details Constructor builds pedigree, without any sorting being done.
		!< If no pedigree file is supplied, all animals are founders
		!< If an animal is in the pedigree, but not in the genotypeFile, this animal is still created as a dummy!
		!< If the animal is in the genotype file, but not in the pedigree, it is added!
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreeGenotypeFiles(pedStructure,fileIn, numberInFile, nSnp,GenotypeFileFormatIn, pedFile, genderfile,dontInitAll)
			use AlphaHouseMod, only : countLines, countColumns
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure
			integer, intent(inout) :: nSnp !< number of snps to read, if 0, will count columns and return
			integer , intent(in), optional :: GenotypeFileFormatIn
			character(len=*),intent(in) :: fileIn !< path of Genotype file
			integer(kind=int32),optional,intent(in) :: numberInFile !< Number of animals in file
			character(len=*),intent(in), optional :: pedFile !< path of pedigree file
			character(len=*),intent(in), optional :: genderfile !< path of gender file
			integer, intent(in), optional :: dontInitAll !< if genotype and phase of all animals should be initialised

			character(len=IDLENGTH) :: tmpId
			integer(kind=int32) :: fileUnit
			integer(kind=int64) :: nIndividuals
			integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
			integer :: tmpAnimalArrayCount
			integer :: i,j,GenotypeFileFormat
			integer(kind=1), dimension(:), allocatable :: tmpGeno
			integer(kind=int64) :: sizeDict
			integer(kind=1), dimension(nSnp * 2) :: WorkVec


			if (nsnp == 0) then
				nsnp = countColumns(fileIn,' ') - 1
			end if

			pedStructure%isSorted = 0
			pedStructure%nHd = 0
			pedStructure%nGenotyped = 0
			pedStructure%nsnpsPopulation = nsnp
			allocate(tmpGeno(nsnp))
			if (present(numberInFile)) then
				nIndividuals = numberInFile
			else
				nIndividuals = countLines(fileIn)
			endif
			pedStructure%addedRealAnimals = nIndividuals
			if  (present(pedFile)) then
				if (present(genderFile)) then
					if (present(dontInitAll)) then
						call initPedigree(pedStructure,pedFile, genderFile=genderFile, nsnps=nSnp, dontInitAll=dontInitAll)
					else
						call initPedigree(pedStructure,pedFile, genderFile=genderFile, nsnps=nSnp)
					end if
				else
					if (present(dontInitAll)) then
						call initPedigree(pedStructure,pedFile, nsnps=nSnp, dontInitAll=dontInitAll)
					else
						call initPedigree(pedStructure,pedFile, nsnps=nSnp)
					end if
				endif
			else

				allocate(pedStructure%sireList)
				allocate(pedStructure%damList)
				allocate(pedStructure%Founders)
				allocate(pedStructure%dictionary)
				pedStructure%nDummys = 0
				tmpAnimalArrayCount = 0
				sizeDict = nIndividuals
				pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
				allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
				pedStructure%pedigreeSize = nIndividuals
				call pedStructure%dictionary%DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
				allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
				allocate(pedStructure%inputMap(nIndividuals))
				pedStructure%maxGeneration = 0
			endif

			if (present(GenotypeFileFormatIn)) then
				GenotypeFileFormat = GenotypeFileFormatIn
			else
				GenotypeFileFormat = 1
			endif

			open(newUnit=fileUnit, file=fileIn, status="old")

			do i=1,nIndividuals
				select case(GenotypeFileFormat)
				case(1)
					read(fileUnit,*) tmpId,tmpGeno(:)
				case(2)
					read (fileUnit, *) tmpId, tmpGeno(:)
					read (fileUnit, *) tmpId, tmpGeno(:)
				case(3)
					read(fileUnit,*) tmpId,WorkVec(:)
				end select

				if (GenotypeFileFormat == 3) then
					do j=1, nsnp
						tmpGeno(j) = MissingGenotypeCode
						if ((WorkVec(j*2 - 1) == 1).and.(WorkVec(j*2) == 1)) tmpGeno(j) = 0
						if ((WorkVec(j*2 - 1) == 1).and.(WorkVec(j*2) == 2)) tmpGeno(j) = 1
						if ((WorkVec(j*2 - 1) == 2).and.(WorkVec(j*2) == 1)) tmpGeno(j) = 1
						if ((WorkVec(j*2 - 1) == 2).and.(WorkVec(j*2) == 2)) tmpGeno(j) = 2
					enddo
				endif
				if (present(pedFile)) then
					j = pedStructure%dictionary%getValue(tmpID)

					if ( j == DICT_NULL) then
						call pedStructure%addAnimalAtEndOfPedigree(tmpID,tmpGeno)
					else
						call pedStructure%setAnimalAsGenotyped(j,tmpGeno)
					endif
				else
					call pedStructure%dictionary%addKey(tmpId, i)
					call pedStructure%Pedigree(i)%initIndividual(trim(tmpId),"0","0", i, nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
					pedStructure%Pedigree(i)%founder = .true.
					call pedStructure%founders%list_add(pedStructure%pedigree(i))
					pedStructure%Pedigree(i)%originalPosition = i
					pedStructure%inputMap(i) = i
					call pedStructure%setAnimalAsGenotyped(i,tmpGeno)
				endif
			enddo
			if (.not. present(dontInitAll)) then
				do i=1, pedStructure%pedigreeSize
					if (.not. pedStructure%pedigree(i)%genotyped) then
						call pedStructure%pedigree(i)%initPhaseAndGenotypes(pedStructure%nsnpsPopulation)

						if (allocated(pedStructure%pedigree(i)%inconsistencies)) then
							deallocate(pedStructure%pedigree(i)%inconsistencies)
						endif
						allocate(pedStructure%pedigree(i)%inconsistencies(pedStructure%nsnpsPopulation))
						pedStructure%pedigree(i)%inconsistencies = 0
						call pedStructure%pedigree(i)%initPhaseArrays(pedStructure%nsnpsPopulation)
					endif


				enddo
			endif

			close(fileUnit)
		end subroutine initPedigreeGenotypeFiles



		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class using two line phase format
		!< @details Constructor builds pedigree, without any sorting being done.
		!< If no pedigree file is supplied, all animals are founders
		!< If an animal is in the pedigree, but not in the genotypeFile, this animal is still created as a dummy!
		!< If the animal is in the genotype file, but not in the pedigree, it is added!
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreePhaseFiles(pedStructure,fileIn, numberInFile, nSnp,pedFile, genderfile,dontInitAll)
			use AlphaHouseMod, only : countLines, countColumns
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure
			integer, intent(inout) :: nSnp !< number of snps to read, if 0, will count columns and return
			character(len=*),intent(in) :: fileIn !< path of Genotype file
			integer(kind=int32),optional,intent(in) :: numberInFile !< Number of lines in file
			character(len=*),intent(in), optional :: pedFile !< path of pedigree file
			character(len=*),intent(in), optional :: genderfile !< path of gender file
			integer, intent(in), optional :: dontInitAll !< if genotype and phase of all animals should be initialised

			character(len=IDLENGTH) :: tmpId
			integer(kind=int32) :: fileUnit
			integer(kind=int64) :: nIndividuals
			integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
			integer :: tmpAnimalArrayCount
			integer(kind=1), dimension(:), allocatable :: tmpPhase1,tmpphase2
			integer(kind=int64) :: sizeDict
			integer :: i,j


			if (nsnp == 0) then
				nsnp = countColumns(fileIn,' ') - 1
			end if

			pedStructure%isSorted = 0
			pedStructure%nHd = 0
			pedStructure%nGenotyped = 0
			pedStructure%nsnpsPopulation = nsnp
			allocate(tmpphase1(nsnp))
			allocate(tmpphase2(nsnp))
			if (present(numberInFile)) then
				nIndividuals = numberInFile/2
			else
				nIndividuals = countLines(fileIn)/2
			endif
			pedStructure%addedRealAnimals = nIndividuals
			if  (present(pedFile)) then
				if (present(genderFile)) then
					if (present(dontInitAll)) then
						call initPedigree(pedStructure,pedFile, genderFile=genderFile, nsnps=nSnp, dontInitAll=dontInitAll)
					else
						call initPedigree(pedStructure,pedFile, genderFile=genderFile, nsnps=nSnp)
					end if
				else
					if (present(dontInitAll)) then
						call initPedigree(pedStructure,pedFile, nsnps=nSnp, dontInitAll=dontInitAll)
					else
						call initPedigree(pedStructure,pedFile, nsnps=nSnp)
					end if
				endif
			else

				allocate(pedStructure%sireList)
				allocate(pedStructure%damList)
				allocate(pedStructure%Founders)
				allocate(pedStructure%dictionary)
				pedStructure%nDummys = 0
				tmpAnimalArrayCount = 0
				sizeDict = nIndividuals
				pedStructure%maxPedigreeSize = nIndividuals + (nIndividuals * 4)
				allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
				pedStructure%pedigreeSize = nIndividuals
				call pedStructure%dictionary%DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
				allocate(tmpAnimalArray(nIndividuals)) !allocate to nIndividuals in case all animals are in incorrect order of generations
				allocate(pedStructure%inputMap(nIndividuals))
				pedStructure%maxGeneration = 0
			endif


			open(newUnit=fileUnit, file=fileIn, status="old")

			do i=1,nIndividuals
				read (fileUnit, *) tmpId, tmpPhase1(:)
				read (fileUnit, *) tmpId, tmpPhase2(:)


				if (present(pedFile)) then
					j = pedStructure%dictionary%getValue(tmpID)

					if ( j == DICT_NULL) then
						call pedStructure%addAnimalAtEndOfPedigree(tmpID)
						call pedStructure%setAnimalAsGenotypedFromPhase(pedStructure%pedigreeSize,tmpPhase1,tmpPhase2)
					else
						call pedStructure%setAnimalAsGenotypedFromPhase(j,tmpPhase1,tmpPhase2)
					endif
				else
					call pedStructure%dictionary%addKey(tmpId, i)
					call pedStructure%Pedigree(i)%initIndividual(trim(tmpId),"0","0", i, nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
					pedStructure%Pedigree(i)%founder = .true.
					call pedStructure%founders%list_add(pedStructure%pedigree(i))
					pedStructure%Pedigree(i)%originalPosition = i
					pedStructure%inputMap(i) = i
					call pedStructure%setAnimalAsGenotypedFromPhase(i,tmpPhase1,tmpPhase2)
				endif
			enddo
			if (.not. present(dontInitAll)) then
				do i=1, pedStructure%pedigreeSize
					if (.not. pedStructure%pedigree(i)%genotyped) then
						call pedStructure%pedigree(i)%initPhaseAndGenotypes(pedStructure%nsnpsPopulation)

						if (allocated(pedStructure%pedigree(i)%inconsistencies)) then
							deallocate(pedStructure%pedigree(i)%inconsistencies)
						endif
						allocate(pedStructure%pedigree(i)%inconsistencies(pedStructure%nsnpsPopulation))
						pedStructure%pedigree(i)%inconsistencies = 0
						call pedStructure%pedigree(i)%initPhaseArrays(pedStructure%nsnpsPopulation)
					endif


				enddo
			endif

			close(fileUnit)
		end subroutine initPedigreePhaseFiles
		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class
		!< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals.
		!< Takes Two arrays as inputs rather than files
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreeArrays(pedStructure, pedArray, genderArray, nsnps)
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure

			character(len=IDLENGTH), dimension(:,:), intent(in) :: pedArray !< array detailing pedigree of format ped([id, sireId, damId], index)
			integer, dimension(:) ,optional , intent(in):: genderArray !< gender array corresponding to index in pedArray
			integer ,optional , intent(in):: nsnps !< number of snps to initialse to
			integer(kind=int32) :: tmpSireNum, tmpDamNum
			integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
			integer :: tmpAnimalArrayCount
			integer :: i
			integer(kind=int64) :: sizeDict
			logical :: sireFound, damFound

			call destroyPedigree(pedStructure)

			pedStructure%nHd = 0
			pedStructure%nGenotyped = 0
			pedStructure%nDummys = 0
			tmpAnimalArrayCount = 0
			pedStructure%nsnpsPopulation = 0

			allocate(pedStructure%sireList)
			allocate(pedStructure%damList)
			allocate(pedStructure%Founders)
			allocate(pedStructure%dictionary)
			if (present(nsnps)) then
				pedStructure%nsnpsPopulation = nsnps
			endif

			pedStructure%isSorted = 0
			sizeDict = size(pedArray,2)
			pedStructure%maxPedigreeSize = size(pedArray,2) + (size(pedArray,2) * 4)
			allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
			pedStructure%pedigreeSize = size(pedArray,2)
			pedStructure%addedRealAnimals = size(pedArray,2)
			call pedStructure%dictionary%DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
			allocate(tmpAnimalArray(size(pedArray,2))) !allocate to nIndividuals in case all animals are in incorrect order of generations
			allocate(pedStructure%inputMap(size(pedArray,2)))
			pedStructure%maxGeneration = 0

			do i=1,size(pedArray(1,:))

				sireFound = .false.
				damFound = .false.

				call pedStructure%dictionary%addKey(pedArray(1,i), i)

				call pedStructure%Pedigree(i)%initIndividual(pedArray(1,i),pedArray(2,i),pedArray(3,i), i,nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
				pedStructure%Pedigree(i)%originalPosition = i
				pedStructure%inputMap(i) = i
				if (pedArray(2,i) /= EMPTY_PARENT) then !check sire is defined in pedigree
					tmpSireNum = pedStructure%dictionary%getValue(pedArray(2,i))
					if (tmpSireNum /= DICT_NULL) then
						sireFound = .true.
					endif
				endif

				if (pedArray(3,i) /= EMPTY_PARENT) then
					tmpDamNum = pedStructure%dictionary%getValue(pedArray(3,i))
					if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
						damFound = .true.
					endif
				endif
				if (pedArray(2,i) == EMPTY_PARENT .and. pedArray(3,i) == EMPTY_PARENT) then !if animal is a founder
					pedStructure%Pedigree(i)%founder = .true.
					call pedStructure%Founders%list_add(pedStructure%Pedigree(i))

				else if (sireFound == .false. .or. damFound == .false. ) then
					tmpAnimalArrayCount = tmpAnimalArrayCount +1
					tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
				else ! if sire and dam are both found
					pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
					call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))

					if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
						call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
					endif
					pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
					call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))

					if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
						call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum)) ! add animal to dam list
					endif
				endif
			enddo

			if (present(genderArray)) then
				do i=1, size(genderArray)

					if (genderArray(i) /=MISSINGGENDERCODE) then
						pedStructure%pedigree(i)%gender = genderArray(i)
					endif
				enddo
			endif

			call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)

			deallocate(tmpAnimalArray)
			write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys
		end subroutine initPedigreeArrays


		!---------------------------------------------------------------------------
		!< @brief Constructor for pedigree class
		!< @details Constructor builds pedigree, without any sorting being done, but by simply building the linked lists and storing founders, as well as having dummy animals.
		!< Takes Two arrays as inputs rather than files
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine initPedigreeIntArrays(pedStructure,pedArray, genderArray)
			use iso_fortran_env
			type(PedigreeHolder) :: pedStructure

			integer, dimension(:,:), intent(in) :: pedArray !< array detailing pedigree of format ped([id, sireId, damId], index)
			integer, dimension(:) ,optional , intent(in):: genderArray !< gender array corresponding to index in pedArray
			integer(kind=int32) :: tmpSireNum, tmpDamNum
			integer, allocatable, dimension(:) :: tmpAnimalArray !array used for animals which parents are not found
			integer :: tmpAnimalArrayCount
			integer :: i
			integer(kind=int64) :: sizeDict
			character(len=IDLENGTH) :: tmpID,tmpSireId, tmpDamID

			logical :: sireFound, damFound

			pedStructure%nHd = 0
			pedStructure%nGenotyped = 0
			pedStructure%nDummys = 0
			tmpAnimalArrayCount = 0

			pedStructure%isSorted = 0
			sizeDict = size(pedArray,2)
			pedStructure%maxPedigreeSize = size(pedArray,2) + (size(pedArray,2) * 4)
			allocate(pedStructure%Pedigree(pedStructure%maxPedigreeSize))
			allocate(pedStructure%dictionary)
			pedStructure%pedigreeSize = size(pedArray,2)
			call pedStructure%dictionary%DictStructure(sizeDict) !dictionary used to map alphanumeric id's to location in pedigree holder
			allocate(tmpAnimalArray(size(pedArray,2))) !allocate to nIndividuals in case all animals are in incorrect order of generations
			allocate(pedStructure%inputMap(size(pedArray,2)))
			pedStructure%maxGeneration = 0
			pedStructure%addedRealAnimals  =size(PedArray,2)
			do i=1,size(pedArray,2)

				sireFound = .false.
				damFound = .false.
				write(tmpId,*) pedArray(1,i)
				write(tmpSireId,*) pedArray(2,i)
				write(tmpDamID,*) pedArray(3,i)
				call pedStructure%dictionary%addKey(tmpId, i)
				call pedStructure%Pedigree(i)%initIndividual(tmpId,tmpSireId,tmpDamID, i,nsnps=pedStructure%nsnpsPopulation) !Make a new individual based on info from ped
				pedStructure%Pedigree(i)%originalPosition = i
				pedStructure%inputMap(i) = i
				if (pedArray(i,2) /= EMPTY_PARENT) then !check sire is defined in pedigree
					tmpSireNum = pedStructure%dictionary%getValue(tmpSireId)
					if (tmpSireNum /= DICT_NULL) then
						sireFound = .true.
					endif
				endif

				if (pedArray(i,3) /= EMPTY_PARENT) then
					tmpDamNum = pedStructure%dictionary%getValue(tmpSireId)
					if (tmpDamNum /= DICT_NULL) then !check dam is defined in pedigree
						damFound = .true.
					endif
				endif
				if (tmpSireId == EMPTY_PARENT .and. tmpDamID == EMPTY_PARENT) then !if animal is a founder
					pedStructure%Pedigree(i)%founder = .true.
					call pedStructure%Founders%list_add(pedStructure%Pedigree(i))

				else if (sireFound == .false. .or. damFound == .false. ) then
					tmpAnimalArrayCount = tmpAnimalArrayCount +1
					tmpAnimalArray(tmpAnimalArrayCount) = i !Set this animals index to be checked later once all information has been read in
				else ! if sire and dam are both found
					pedStructure%Pedigree(i)%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
					call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(i))

					if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
						call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
					endif
					pedStructure%Pedigree(i)%damPointer =>  pedStructure%Pedigree(tmpDamNum)
					call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(i))

					if (pedStructure%Pedigree(tmpDamNum)%nOffs == 1) then
						call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
						call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamNum)) ! add animal to sire list
					endif
				endif
			enddo

			if (present(genderArray)) then

				do i=1, size(genderArray)

					if (genderArray(i) /=MISSINGGENDERCODE) then
						pedStructure%pedigree(i)%gender = genderArray(i)
					endif
				enddo
			endif
			call addOffspringsAfterReadIn(pedStructure, tmpAnimalArray,tmpAnimalArrayCount)
			deallocate(tmpAnimalArray)
			write (error_unit,*) "NOTE: Number of Dummy Animals: ",pedStructure%nDummys
		end subroutine initPedigreeIntArrays

		!---------------------------------------------------------------------------
		!< @brief Helper function to avoid code duplication
		!< required by constructors to determine animals that need offspring info added
		!< due to it not being initially available in the pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine addOffspringsAfterReadIn(pedStructure, tmpAnimalArray, tmpAnimalArrayCount)
			use ConstantModule, only : IDLENGTH,EMPTY_PARENT
			use IFCORE
			class(PedigreeHolder), intent(inout) :: pedStructure
			integer, dimension(:), intent(in) :: tmpAnimalArray !< array containing indexes of tmp animals
			integer, intent(in) :: tmpAnimalArrayCount !< number of animals actually in tmpAnimalArray
			logical :: sireFound, damFound
			integer(kind=int32) :: tmpSireNum, tmpDamNum,tmpId
			integer(kind=int32) :: i, tmpCounter
			character(len=IDLENGTH) :: tmpSire,tmpDam

			tmpCounter = 0 !< counter for dummy animals

			!check animals that didn't have parental information initially
			! this is done to avoid duplication when a pedigree is sorted
			do i=1,tmpAnimalArrayCount
				sireFound = .false.
				damFound = .false.
				tmpSire = pedStructure%Pedigree(tmpAnimalArray(i))%getSireId()
				tmpSireNum = pedStructure%dictionary%getValue(tmpSire)
				tmpDam = pedStructure%Pedigree(tmpAnimalArray(i))%getDamId()
				tmpDamNum = pedStructure%dictionary%getValue(tmpDam)

				if (tmpSire /= EMPTY_PARENT) then

					! check that we've not already defined the parent above
					if (tmpSireNum /= DICT_NULL .and. .not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer)) then !if sire has been found in hashtable
						pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer =>  pedStructure%Pedigree(tmpSireNum)
						call pedStructure%Pedigree(tmpSireNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
						if (pedStructure%Pedigree(tmpSireNum)%nOffs == 1) then
							call pedStructure%Pedigree(tmpSireNum)%setGender(1) !if its a sire, it should be male
							call pedStructure%sireList%list_add(pedStructure%Pedigree(tmpSireNum)) ! add animal to sire list
						endif
						! check that we've not already defined the parent above
					else if (.not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%sirePointer)) then!if sire is defined but not in the pedigree, create him
						! check if the tmp animal has already been created
						call pedStructure%addAnimalAtEndOfPedigree(originalID=trim(tmpSire),offspringID=tmpAnimalArray(i))
					endif
					sireFound = .true.
				endif

				if (tmpDam /= EMPTY_PARENT) then
					! check that we've not already defined the parent above
					if (tmpDamNum /= DICT_NULL .and. .not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%damPointer)) then !if dam has been found
						pedStructure%Pedigree(tmpAnimalArray(i))%damPointer =>  pedStructure%Pedigree(tmpDamNum)
						call pedStructure%Pedigree(tmpDamNum)%addOffspring(pedStructure%Pedigree(tmpAnimalArray(i)))
						if (pedStructure%Pedigree(tmpDamnum)%nOffs == 1) then
							call pedStructure%Pedigree(tmpDamNum)%setGender(2) !if its a dam, should be female
							call pedStructure%damList%list_add(pedStructure%Pedigree(tmpDamnum)) ! add animal to dam list
						endif
						! check that we've not already defined the parent above
					else if (.not. associated(pedStructure%Pedigree(tmpAnimalArray(i))%damPointer)) then
						! Check for defined animals that have nit been set in pedigree
						call pedStructure%addAnimalAtEndOfPedigree(originalID=trim(tmpDam),offspringID=tmpAnimalArray(i))
						
					endif
					damFound = .true.
				endif


				if (.not. damFound .OR. .not. sireFound) then

					if (.not. damFound) then
						call pedStructure%createDummyAnimalAtEndOfPedigree(tmpId, tmpAnimalArray(i), .false.)
						if (tmpDam == EMPTY_PARENT) then
							pedStructure%unknownDummys = pedStructure%unknownDummys+1
							pedStructure%Pedigree(tmpId)%isUnknownDummy = .true.
						endif

					endif
					if (.not. sireFound) then

						call pedStructure%createDummyAnimalAtEndOfPedigree(tmpId, tmpAnimalArray(i), .true.)
						if (tmpSire == EMPTY_PARENT) then
							pedStructure%unknownDummys = pedStructure%unknownDummys+1
							pedStructure%Pedigree(pedStructure%pedigreeSize)%isUnknownDummy = .true.
						endif
					endif

				endif
			enddo

		end subroutine addOffspringsAfterReadIn





		!---------------------------------------------------------------------------
		!< @brief returns a list of animals that are genotyped, and are classed as founders
		!< Animals are classed as founders if they have no ancestors that are genotyped in a given number of generations
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getGenotypedFounders(this, numberOfGenerations) result(genotypedFounders)
			class(pedigreeHolder) :: this
			integer, intent(in) :: numberOfGenerations
			type(IndividualLinkedList) :: genotypedFounders
			integer :: i
			do i=1, this%pedigreeSize
				if (this%pedigree(i)%isGenotypedNonMissing()) then

					if (this%pedigree(i)%founder) then
						call genotypedFounders%list_add(this%pedigree(i))
					else if (.not. this%pedigree(i)%hasGenotypedAnsestors(numberOfGenerations)) then !< this checks if ancestors are not genotyped given a number
						call genotypedFounders%list_add(this%pedigree(i))
					endif

				endif
			enddo
		end function getGenotypedFounders



		subroutine wipeGenotypeAndPhaseInfo(this)
			class(pedigreeHolder) :: this

			integer :: i

			this%nGenotyped = 0
			! call this%genotypeDictionary%destroy()

			if (allocated(this%genotypeDictionary)) then
				deallocate(this%genotypeDictionary)
			endif
			if (allocated(this%genotypeMap)) then
				deallocate(this%genotypeMap)
			endif
			if (allocated(this%hdDictionary)) then
				deallocate(this%hdDictionary)
			endif
			if (allocated(this%hdMap)) then
				deallocate(this%hdMap)
			endif
			! this%genotypeMap = 0
			do i=1,this%pedigreeSize
				if (allocated(this%pedigree(i)%individualPhase)) then
					deallocate(this%pedigree(i)%individualPhase)
				endif

				if (allocated(this%pedigree(i)%individualGenotype)) then
					deallocate(this%pedigree(i)%individualGenotype)
				endif
			enddo

		end subroutine wipeGenotypeAndPhaseInfo

		!---------------------------------------------------------------------------
		!< @brief distructor for pedigree class
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine destroyPedigree(this)
			type(PedigreeHolder) :: this

			if (ASSOCIATED(this%pedigree)) then
				deallocate(this%pedigree)
				this%pedigree => null()
			endif
			if (allocated(this%generations)) then

				deallocate(this%generations)
			endif

			if (allocated(this%sireList)) then
				deallocate(this%sireList)
			endif

			if (allocated(this%damList)) then
				deallocate(this%damList)
			endif
			if (allocated(this%founders)) then
				deallocate(this%founders)
			endif

			if (allocated(this%inputMap)) then
				deallocate(this%inputMap)
			endif
			if (allocated(this%dictionary)) then
				deallocate(this%dictionary)
			endif

			if (this%nGenotyped > 0) then
				deallocate(this%genotypeDictionary)
				deallocate(this%genotypeMap)
				this%nGenotyped = 0
			endif

			if (this%nHd > 0) then

				deallocate(this%hdMap)
				this%nHd = 0
			endif

			this%pedigreeSize = 0

			if (allocated(this%uniqueParentList)) then

				deallocate(this%uniqueParentList)
			endif
		end subroutine destroyPedigree


		!---------------------------------------------------------------------------
		! DESCRIPTION:
		!< @brief     Adds genotype information to pedigree from a 2 dimensional array
		!
		!< @author     David Wilson, david.wilson@roslin.ed.ac.uk
		!
		!< @date       October 25, 2016
		!---------------------------------------------------------------------------
		subroutine addGenotypeInformationFromArray(this, array, initall)

			use AlphaHouseMod, only : countLines
			implicit none
			class(PedigreeHolder) :: this
			integer(kind=1),allocatable,dimension (:,:), intent(in) :: array !< array should be dimensions nanimals, nsnp
			integer :: i
			integer, intent(in),optional :: initAll !< optional argument- if present initialise whoe pedigree with size of snps
			do i=1,size(array,1) !< assumes dummys are at end, otherwise this will NOT work
				call this%setAnimalAsGenotyped(i, array(i,:))
			enddo
			this%nsnpsPopulation = size(array,1)

			if (present(initAll)) then
				do i=1, this%pedigreeSize
					if (.not. this%pedigree(i)%genotyped) then
						call this%pedigree(i)%initPhaseAndGenotypes(this%nsnpsPopulation)
						allocate(this%pedigree(i)%inconsistencies(this%nsnpsPopulation))
						this%pedigree(i)%inconsistencies = 0
						call this%pedigree(i)%initPhaseArrays(this%nsnpsPopulation)
					endif
				enddo
			endif
		end subroutine addGenotypeInformationFromArray


		!---------------------------------------------------------------------------
		! DESCRIPTION:
		!< @brief     Adds genotype information to pedigree from a file
		!
		!< @author     David Wilson, david.wilson@roslin.ed.ac.uk
		!
		!< @date       October 25, 2016
		!--------------------------------------------------------------------------
		subroutine addGenotypeInformationFromFile(this, genotypeFile, nsnps, nAnnisG, startSnp, endSnp, lockIn, initAll)

			use AlphaHouseMod, only : countLines, countColumns
			implicit none
			class(PedigreeHolder) :: this
			character(len=*) :: genotypeFile
			character(len=IDLENGTH) :: tmpID
			integer,intent(inout) :: nsnps !< if nsnps ==0
			integer,intent(in),optional :: nAnnisG
			integer, intent(in),optional :: startSnp, endSnp
			logical, intent(in), optional :: lockIn
			integer, intent(in), optional :: initAll
			integer(kind=1), allocatable, dimension(:) :: tmpSnpArray
			integer :: i, j,fileUnit, nAnnis,tmpIdNum
			integer :: count, end
			integer :: addCount
			logical :: lock

			addCount = 0
			if (present(lockIn)) then
				lock = lockIn
			else
				lock = .false.

			endif
			if (present(nAnnisG)) then
				nAnnis = nAnnisG
			else
				nAnnis = countLines(genotypeFile)
			endif

			if (nsnps == 0) then
				nsnps = countColumns(genotypeFile, ' ') - 1
			endif


			if (this%nGenotyped /= 0) then
				deallocate(this%genotypeMap)
				deallocate(this%genotypeDictionary)
				this%nGenotyped = 0
			endif

			allocate(tmpSnpArray(nsnps))
			open(newUnit=fileUnit, file=genotypeFile, status="old")
			do i=1, nAnnis
				read (fileUnit,*) tmpId,tmpSnpArray(:)
				do j=1,nsnps
					if ((tmpSnpArray(j)<0).or.(tmpSnpArray(j)>2)) then
						tmpSnpArray(j)=9
					endif
				enddo

				tmpIdNum = this%dictionary%getValue(tmpId)
				if (tmpIdNum == DICT_NULL) then
					addCount = addCount +1
					write(error_unit, *) "WARNING: Genotype info for non existing animal here:",trim(tmpId), " file:", trim(genotypeFile), " line:",i
					write(error_unit, *) "Animal will be added as a founder to pedigree"

					call this%addAnimalAtEndOfPedigree(trim(tmpID), tmpSnpArray)

				else

					if (present(startSnp)) then
						count = 0

						if (present(endSnp)) then
							end =endSnp
						else
							end = nsnps
						endif
						call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray(startSnp:endSnp),lock)
					else if (present(endSnp)) then
						count = 0
						call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray(1:endsnp),lock)
					else
						call this%setAnimalAsGenotyped(tmpIdNum, tmpSnpArray,lock)
					endif
				endif
			enddo

			this%nsnpsPopulation = nsnps
			if (present(initAll)) then
				do i=1, this%pedigreeSize
					if (.not. this%pedigree(i)%genotyped) then
						call this%pedigree(i)%initPhaseAndGenotypes(this%nsnpsPopulation)

						if (allocated(this%pedigree(i)%inconsistencies)) then
							deallocate(this%pedigree(i)%inconsistencies)
						endif
						allocate(this%pedigree(i)%inconsistencies(this%nsnpsPopulation))
						this%pedigree(i)%inconsistencies = 0
						call this%pedigree(i)%initPhaseArrays(this%nsnpsPopulation)
					endif


				enddo
			endif

			write(output_unit,*) "NOTE: Number of Genotyped animals: ",this%nGenotyped
			write(output_unit,*) "NOTE: Number of Founders added to pedigree: ",addCount

		end subroutine addGenotypeInformationFromFile



		subroutine addFamilyIds(this, familyIds)

			class(PedigreeHolder) :: this
			character(len=IDLENGTH),dimension(:),allocatable :: familyIds !< assumed to be same order as pedigree
			integer :: i, familyIdCount

			if (allocated(this%familyIdDict)) then
				deallocate(this%familyIdDict)
			endif
			allocate(this%familyIdDict)
			call this%familyIdDict%dict_create()

			familyIdCount = 0

			do i=1, size(familyIds)

				if (this%familyIdDict%getValue(familyIds(i))==  DICT_NULL) then
					familyIdCount = familyIdCount + 1
					call this%familyIdDict%addKey(familyIds(i), familyIdCount)

				endif
				this%pedigree(i)%familyId = familyids(i)

			enddo

		end subroutine addFamilyIds

		!---------------------------------------------------------------------------
		! DESCRIPTION:
		!< @brief     Adds phase information to pedigree from a file
		!
		!< @author    Daniel Money, daniel.money@roslin.ed.ac.uk
		!
		!< @date       June 19, 2017
		!--------------------------------------------------------------------------
		subroutine addPhaseInformationFromFile(this, phaseFile, nsnps, nAnnisG)

			use AlphaHouseMod, only : countLines
			implicit none
			class(PedigreeHolder) :: this
			character(len=*) :: phaseFile
			character(len=IDLENGTH) :: tmpID
			integer,intent(in) :: nsnps
			integer,intent(in),optional :: nAnnisG
			integer(kind=1), allocatable, dimension(:) :: tmpSnpArray
			integer :: i, j, h, fileUnit, nAnnis,tmpIdNum

			if (present(nAnnisG)) then
				nAnnis = nAnnisG
			else
				nAnnis = int(float(countLines(phaseFile)) / 2.0)
			endif

			allocate(tmpSnpArray(nsnps))
			open(newUnit=fileUnit, file=phaseFile, status="old")
			do i=1, nAnnis
				do h = 1, 2
					read (fileUnit,*) tmpId,tmpSnpArray(:)
					do j=1,nsnps
						if ((tmpSnpArray(j)<0).or.(tmpSnpArray(j)>1)) tmpSnpArray(j)=9
					enddo
					tmpIdNum = this%dictionary%getValue(tmpId)
					if (tmpIdNum == DICT_NULL) then
						write(error_unit, *) "WARNING: Phase info for non existing animal here:",trim(tmpId), " file:", trim(phaseFile), " line:",i
					else
						call this%pedigree(tmpIdNum)%setPhaseArray(h, tmpSnpArray)
						this%pedigree(tmpIdNum)%isPhased = .true.
					endif
				end do
			enddo

		end subroutine addPhaseInformationFromFile


		!---------------------------------------------------------------------------
		! DESCRIPTION:
		!< @brief     Adds sequence information from VCF file - this was taken out of Roberto's HMM
		!< @date       June 19, 2017
		!--------------------------------------------------------------------------
		subroutine addSequenceFromFile(this, seqFile, nsnps, nAnisGIn,maximumReads, startSnp, endSnp)

			use AlphaHouseMod, only : countLines,countColumns
			use ConstantModule, only : IDLENGTH,DICT_NULL
			use iso_fortran_env
			use IFCORE
			implicit none
			class(PedigreeHolder) :: this
			character(len=*) :: seqFile
			integer,intent(in) :: nsnps
			integer, intent(in), optional :: maximumReads
			integer,intent(in),optional :: nAnisGIn
			integer,intent(in),optional :: startSnp, endSnp
			integer :: nanisG, end
			! type(Pedigreeholder), intent(inout) :: genotype
			integer(KIND=1), allocatable, dimension(:) :: tmp
			integer,allocatable, dimension(:) :: ref, alt
			integer :: unit, tmpID,i,j
			character(len=IDLENGTH) :: seqid !placeholder variables
			character(len=1), dimension(3):: delimiter
			integer :: nCol

			if (.not. Present(nAnisGIn)) then
				NanisG = countLines(seqFile)/2
			else
				nanisG = nAnisGIn
			endif

			delimiter(1) = ","
			delimiter(2) = " "
			delimiter(3) = char(9)

			nCol=countColumns(trim(seqFile), delimiter)-1 ! First column is animal id

			open(newunit=unit,FILE=trim(seqFile),STATUS="old") !INPUT FILE
			allocate(ref(nCol))
			allocate(alt(nCol))

			tmp = 9
			print *, "Number of animals in seq file", NanisG
			do i=1,nAnisG
				read (unit,*) seqid, ref(:)
				read (unit,*) seqid, alt(:)

				tmpID = this%dictionary%getValue(seqid)
				if (present(maximumReads)) then
					do j=1,nsnps
						if (ref(j)>=maximumReads) ref(j)=maximumReads-1
						if (alt(j)>=maximumReads) alt(j)=maximumReads-1
					enddo
				endif

				if (tmpID /= DICT_NULL) then

					if (present(startSnp)) then
						if (present(endSnp)) then
							end = endsnp
							if (end > nCol) then
								call TRACEBACKQQ(string= "ERROR - endSNP is greater than number of columns in sequence file",user_exit_code=1)
							endif
						else
							end = nsnps
						endif
						call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref(startSnp:end),alt(startSnp:end))
					else if (present(endSnp)) then
						call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref(1:endsnp),alt(1:endSnp))
					else
						call this%setAnimalAsGenotypedSequence(tmpID,tmp,ref,alt)
					endif
				endif
			end do

			close(unit)
		end subroutine addSequenceFromFile



		!---------------------------------------------------------------------------
		! DESCRIPTION:
		!< @brief     Adds sequence information from VCF file - this was taken out of Mara's code and
		!< adapted to work with the pedigree class
		!< @date       June 19, 2017
		!--------------------------------------------------------------------------
		subroutine addSequenceFromVCFFile(this,seqFile,nSnpsIn,nAnisIn,position,quality,chr,startPos,EndPos)

			use omp_lib
			use AlphaHouseMod, only : countLines,countColumns
			use ConstantModule, only : IDLENGTH,DICT_NULL

			implicit none
			class(PedigreeHolder)        :: this
			character(len=*), intent(in) :: seqFile
			character(len=300),intent(in) :: chr
			integer,intent(in) :: StartPos,EndPos

			integer,intent(in),optional :: nSnpsIn
			integer,intent(in),optional :: nAnisIn


			integer:: nAnis,nSnp,unit,pos,i,j,tmpID
			character(len=1), dimension(3):: delimiter

			character(len=IDLENGTH), allocatable, dimension(:)              :: Ids
			integer(int64), allocatable, dimension(:),intent(out),optional  :: position
			real(real64), allocatable, dimension(:),intent(out),optional    :: quality

			character(100)  :: tCHROM,tREF,tALT
			integer(int64)  :: tPOS
			real(real64)    :: tQUAL

			character(len=100), dimension(:), allocatable:: dumC

			integer(KIND=1), allocatable, dimension(:) :: tmp

			integer, dimension(:,:), allocatable:: tmpSequenceData
			integer, dimension(:,:,:), allocatable:: SequenceData

			real(kind=8)::tstart,tend

			delimiter(1) = ","
			delimiter(2) = " "
			delimiter(3) = char(9)


			if (.not. Present(nSnpsIn)) then
				nSnp = countLines(seqFile)-1 ! First row is the header
			else
				nSnp = nSnpsIn
			endif

			if (.not. Present(nAnisIn)) then
				nAnis = countColumns(seqFile, delimiter)-5 ! First 5 columns are "CHROM POS REF ALT QUAL"
			else
				nAnis = nAnisIn
			endif

			if (Present(position)) allocate(position(nSnp))
			if (Present(quality)) allocate(quality(nSnp))

			open(newunit=unit,FILE=trim(seqFile),STATUS="old") !INPUT FILE

			allocate(Ids(nAnis))
			allocate(dumC(5))
			allocate(tmp(nSnp))
			tmp = 9
			allocate(SequenceData(nAnis, nSnp, 2))
			allocate(tmpSequenceData(nAnis,2))

			! Read header and store animals id
			read (unit,*) dumC,ids
			deallocate(dumC)

			tstart = omp_get_wtime()

			pos=1
			do j = 1, nSnp
				read(unit, *) tCHROM, tPOS, tREF, tALT, tQUAL, (tmpSequenceData(i,1), tmpSequenceData(i,2), i =1, nAnis)
				if (trim(tCHROM).eq.trim(chr)) then
					if (((tPos.ge.StartPos).or.(StartPos.eq.0)).and.((tPos.le.EndPos).or.(EndPos.eq.0))) then
						if (Present(position)) position(pos)=tPOS
						if (Present(quality)) quality(pos)=tQUAL
						SequenceData(:,pos,1)=tmpSequenceData(:,1)
						SequenceData(:,pos,2)=tmpSequenceData(:,2)
						pos=pos+1
					endif
				end if
			end do
			close(unit)

			do i=1, nAnis
				tmpID = this%dictionary%getValue(ids(i))
				if (tmpID /= DICT_NULL) then
					call this%setAnimalAsGenotypedSequence(tmpID,tmp,SequenceData(i,:,1),SequenceData(i,:,2))
				endif
			enddo

			tend = omp_get_wtime()

			write(*,*) "Total wall time for Importing Reads", tend - tstart

		end subroutine addSequenceFromVCFFile



		!---------------------------------------------------------------------------
		!< @brief builds correct generation information by looking at founders
		!< This is effectively a sort function for the pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setPedigreeGenerationsAndBuildArrays(this)

			implicit none
			class(PedigreeHolder) :: this
			integer :: i
			type(IndividualLinkedListNode), pointer :: tmpIndNode
			tmpIndNode => this%Founders%first
			allocate(this%generations(0:generationThreshold))
			do i=1, this%Founders%length
				call this%setOffspringGeneration(tmpIndNode%item)
				tmpIndNode => tmpIndNode%next
			end do

		end subroutine setPedigreeGenerationsAndBuildArrays


		!---------------------------------------------------------------------------
		!< @brief returns true if individual at given index is a isDummy
		!< if 0 is given, return false
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!< @param[in] file path (string)
		!---------------------------------------------------------------------------
		logical function isDummy(this, id)
			implicit none
			class(PedigreeHolder) :: this
			integer, intent(in) :: id !< ID to check if animal is dummy

			if (id == 0) then
				isDummy = .false.
			else if (this%pedigree(id)%isDummy) then
				isDummy = .true.
			else
				isDummy = .false.
			endif
		end function isDummy


		!---------------------------------------------------------------------------
		!< @brief writes sorted pedigree information to either a file or stdout
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!< @param[in] file path (string)
		!---------------------------------------------------------------------------
		subroutine outputSortedPedigree(this,file)
			use iso_fortran_env, only : output_unit
			class(PedigreeHolder) :: this
			character(len=*), intent(in), optional :: file !< output path for sorted pedigree
			integer :: unit, i,h
			type(IndividualLinkedListNode), pointer :: tmpIndNode
			character(len=:), allocatable :: fmt

			if (.not. allocated(this%generations)) then
				call this%setPedigreeGenerationsAndBuildArrays
			endif

			if (present(file)) then
				open(newUnit=unit, file=file, status="unknown")
			else
				unit = output_unit
			endif
			fmt = "(3a"//Int2Char(IDLENGTH)//", i"//Int2Char(IDINTLENGTH)//")"
			do i=0, this%maxGeneration
				tmpIndNode => this%generations(i)%first
				do h=1, this%generations(i)%length
					write(unit, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId, tmpIndNode%item%generation
					! write(*, fmt) tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
					tmpIndNode => tmpIndNode%next
				end do
			enddo
			if (present(file)) then
				close(unit)
			endif
		end subroutine outputSortedPedigree


		!---------------------------------------------------------------------------
		!< @brief Output pedigree to stdout in the format recodedID, recodedSireId, recodedDamId, originalId
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine outputSortedPedigreeInAlphaImputeFormat(this, file)
			use iso_fortran_env, only : output_unit
			class(PedigreeHolder) :: this
			character(len=*), intent(in), optional :: file !< output path for sorted pedigree
			integer :: unit, i,h, sortCounter
			type(IndividualLinkedListNode), pointer :: tmpIndNode
			sortCounter = 0
			if (.not. allocated(this%generations)) then
				call this%setPedigreeGenerationsAndBuildArrays
			endif
			if (present(file)) then
				open(newUnit=unit, file=file, status="unknown")
			else
				unit = output_unit
			endif

			block
				integer :: sireId, damId
				character(len=:), allocatable :: fmt
				fmt = "(3i"//Int2Char(IDINTLENGTH)//", a"//Int2Char(IDLENGTH)//")"
				do i=0, this%maxGeneration
					tmpIndNode => this%generations(i)%first
					do h=1, this%generations(i)%length
						if (associated(tmpIndNode%item%damPointer)) then
							damId = tmpIndNode%item%damPointer%id
						else
							damId = 0
						endif
						if (associated(tmpIndNode%item%sirePointer)) then
							sireId = tmpIndNode%item%sirePointer%id
						else
							sireId = 0
						endif
						sortCounter = sortCounter +1
						write (unit, fmt) tmpIndNode%item%id,sireId,damId, tmpIndNode%item%originalID
						! write(*,'(a,",",a,",",a,",",i8)') tmpIndNode%item%originalID,tmpIndNode%item%sireId,tmpIndNode%item%damId,tmpIndNode%item%generation
						tmpIndNode => tmpIndNode%next
					end do
				enddo
			endblock

			if (present(file)) then !avoids closing stdout
				close(unit)
			endif
		end subroutine outputSortedPedigreeInAlphaImputeFormat




		!---------------------------------------------------------------------------
		!< @brief Sorts pedigree, and overwrites all fields to new values
		!< @details effectively, does a deep copy to sort pedigree based on generation, but puts dummys at bottom
		!< If value is given for unknownDummysAtEnd, then only unknown dummys will be put at the end
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine sortPedigreeAndOverwrite(this, unknownDummysAtEnd)
			use iso_fortran_env, only : output_unit, int64
			use IFCORE
			class(PedigreeHolder) :: this
			integer :: i,h, pedCounter, tmpId,tmpGenotypeMapIndex
			integer(kind=int64) :: sizeDict
			type(IndividualLinkedListNode), pointer :: tmpIndNode
			type(Individual), pointer, dimension(:) :: newPed
			type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList
			type(IndividualLinkedList) :: dummyList
			integer, intent(in) , optional :: unknownDummysAtEnd !< if this option is specified, then only unknown dummies are put at end

			if (this%isSorted /= 0) return

			if (allocated(this%generations)) deallocate(this%generations)
			if (.not. allocated(this%generations)) then
				call this%setPedigreeGenerationsAndBuildArrays
			endif


			pedCounter = 0

			deallocate(this%sireList)
			deallocate(this%damList)
			deallocate(this%founders)
			deallocate(this%dictionary)
			sizeDict = this%pedigreeSize
			allocate(newPed(this%maxPedigreeSize))
			allocate(this%dictionary)
			call this%dictionary%DictStructure(sizeDict)
			allocate(newGenerationList(0:this%maxGeneration))
			allocate(this%sireList)
			allocate(this%damList)
			allocate(this%founders)


			do i=0, this%maxGeneration
				tmpIndNode => this%generations(i)%first
				do h=1, this%generations(i)%length
					if (present(unknownDummysAtEnd)) then
						if (tmpIndNode%item%isUnknownDummy) then
							call dummyList%list_add(tmpIndNode%item)
							tmpIndNode => tmpIndNode%next
							cycle
						endif
					else
						if (tmpIndNode%item%isDummy) then
							call dummyList%list_add(tmpIndNode%item)
							tmpIndNode => tmpIndNode%next
							cycle
						endif
					endif

					pedCounter = pedCounter +1
					! update dictionary index

					call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)

					!  update genotype map
					if (this%nGenotyped > 0) then
						tmpGenotypeMapIndex = this%genotypeDictionary%getValue(tmpIndNode%item%originalID)
						if (tmpGenotypeMapIndex /= DICT_NULL) then
							this%genotypeMap(tmpGenotypeMapIndex) = pedCounter
						endif
					endif
					! Update hd map
					if (this%nHd > 0) then
						tmpGenotypeMapIndex = this%hdDictionary%getValue(tmpIndNode%item%originalID)
						if (tmpGenotypeMapIndex /= DICT_NULL) then
							this%hdMap(tmpGenotypeMapIndex) = pedCounter
						endif
					endif
					! Set the location to the pedigree to the new value
					! newPed(pedCounter) = tmpIndNode%item
					call copyIndividual(newPed(pedCounter),tmpIndNode%item)
					! take the original id, and update it
					if(.not. newPed(pedCounter)%isDummy .and. newPed(pedCounter)%originalPosition /= 0) then
						this%inputMap(newPed(pedCounter)%originalPosition) = pedCounter
					endif
					newPed(pedCounter)%id = pedCounter
					call newPed(pedCounter)%resetOffspringInformation ! reset offsprings
					if (associated(newPed(pedCounter)%sirePointer)) then
						tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
						if (tmpID /= DICT_NULL) then
							call newPed(tmpId)%addOffspring(newPed(pedCounter))
							newPed(pedCounter)%sirePointer=> newPed(tmpId)
							if (newPed(tmpId)%nOffs == 1) then
								call this%sireList%list_add(newPed(tmpId))
							endif
						endif
					endif

					if (associated(newPed(pedCounter)%damPointer)) then
						tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
						if (tmpID /= DICT_NULL) then
							call newPed(tmpId)%addOffspring(newPed(pedCounter))
							newPed(pedCounter)%damPointer=> newPed(tmpId)
							if (newPed(tmpId)%nOffs == 1) then
								call this%damList%list_add(newPed(tmpId))
							endif
						endif
					endif

					if (i ==0 ) then !if object is afounder add to founder array
						call  this%founders%list_add(newPed(pedCounter))
					endif

					call newGenerationList(i)%list_add(newPed(pedCounter))
					tmpIndNode => tmpIndNode%next
				end do
			enddo

			tmpIndNode => dummyList%first

			! add dummys to end of pedigree
			do i=1, dummyList%length
				pedCounter = pedCounter +1
				call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)
				call copyIndividual(newPed(pedCounter),tmpIndNode%item)
				newPed(pedCounter)%id = pedCounter
				call newPed(pedCounter)%resetOffspringInformation() ! reset offsprings

				if (newPed(pedCounter)%generation == NOGENERATIONVALUE) then

					call TRACEBACKQQ(string= "ERROR: Circular pedigree structure given on animal "//newPed(pedCounter)%originalId,user_exit_code=1)

				endif
				do h=1, tmpIndNode%item%nOffs
					tmpId =  this%dictionary%getValue(tmpIndNode%item%offsprings(h)%p%originalID)

					if (tmpId == DICT_NULL .or. tmpIndNode%item%offsprings(h)%p%generation == NOGENERATIONVALUE) then

						call TRACEBACKQQ(string= "ERROR: Circular pedigree structure given on animal "//tmpIndNode%item%offsprings(h)%p%originalID,user_exit_code=1)

					endif
					call newPed(pedCounter)%addOffspring(newPed(tmpID))
					if(associated(tmpIndNode%item%offsprings(h)%p%sirePointer, tmpIndNode%item)) then
						newPed(tmpId)%sirePointer => newPed(pedCounter)
						if (newPed(pedCounter)%nOffs == 1) then
							call this%sireList%list_add(newPed(pedCounter))
						endif
					else
						newPed(tmpId)%damPointer => newPed(pedCounter)
						if (newPed(pedCounter)%nOffs == 1) then
							call this%damList%list_add(newPed(pedCounter))
						endif
					endif
				enddo
				call  this%founders%list_add(newPed(pedCounter))

				call newGenerationList(0)%list_add(newPed(pedCounter))
				tmpIndNode => tmpIndNode%next
			enddo

			deallocate(this%pedigree)

			this%pedigree => newPed

			newPed => null()
			do i = 0, this%maxGeneration
				call destroyLinkedList(this%generations(i))
			enddo
			! this%generations = newGenerationList

			! deallocate(this%generations)
			! allocate(this%generations(0:this%maxGeneration))
			this%generations(0:this%maxGeneration) = newGenerationList

			!
			if (present(unknownDummysAtEnd)) then
				this%isSorted = 2
			else
				this%isSorted = 1
			endif
		end subroutine sortPedigreeAndOverwrite


		!---------------------------------------------------------------------------
		!< @brief Output pedigree to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine sortPedigreeAndOverwriteWithDummyAtTheTop(this)
			use iso_fortran_env, only : output_unit, int64
			class(PedigreeHolder) :: this
			integer :: i,h, pedCounter, tmpId,tmpGenotypeMapIndex
			integer(kind=int64) :: sizeDict
			type(IndividualLinkedListNode), pointer :: tmpIndNode
			type(Individual), pointer, dimension(:) :: newPed
			type(IndividualLinkedList),allocatable, dimension(:) :: newGenerationList
			if (.not. allocated(this%generations)) then
				call this%setPedigreeGenerationsAndBuildArrays
			endif
			pedCounter = 0

			! deallocate to call destructors
			deallocate(this%dictionary)
			deallocate(this%founders)
			deallocate(this%sireList)
			deallocate(this%damList)

			sizeDict = this%pedigreeSize
			allocate(newPed(this%maxPedigreeSize))
			allocate(this%dictionary)
			call this%dictionary%DictStructure(sizeDict)
			allocate(newGenerationList(0:this%maxGeneration))

			! reallocate to use again
			allocate(this%founders)
			allocate(this%sireList)
			allocate(this%damList)
			do i=0, this%maxGeneration
				tmpIndNode => this%generations(i)%first
				do h=1, this%generations(i)%length
					pedCounter = pedCounter +1
					call this%dictionary%addKey(tmpIndNode%item%originalID,pedCounter)

					!  update genotype map
					if (this%nGenotyped > 0) then
						tmpGenotypeMapIndex = this%genotypeDictionary%getValue(tmpIndNode%item%originalID)
						if (tmpGenotypeMapIndex /= DICT_NULL) then
							this%genotypeMap(tmpGenotypeMapIndex) = pedCounter
						endif
					endif
					! Update hd map
					if (this%nHd > 0) then
						tmpGenotypeMapIndex = this%hdDictionary%getValue(tmpIndNode%item%originalID)
						if (tmpGenotypeMapIndex /= DICT_NULL) then
							this%hdMap(tmpGenotypeMapIndex) = pedCounter
						endif
					endif

					newPed(pedCounter) = tmpIndNode%item
					if (.not. newPed(pedCounter)%isDummy .and.  newPed(pedCounter)%originalPosition /= 0) then
						! take the original id, and update it - we don't want dummies in this list
						this%inputMap(newPed(pedCounter)%originalPosition) = pedCounter
					endif
					newPed(pedCounter)%id = pedCounter
					call newPed(pedCounter)%resetOffspringInformation ! reset offsprings

					if (associated(newPed(pedCounter)%sirePointer)) then
						tmpId =  this%dictionary%getValue(newPed(pedCounter)%sirePointer%originalID)
						if (tmpID /= DICT_NULL) then
							call newPed(tmpId)%addOffspring(newPed(pedCounter))
							newPed(pedCounter)%sirePointer=> newPed(tmpId)
							if (newPed(tmpId)%nOffs == 1) then
								call this%sireList%list_add(newPed(tmpId))
							endif

						endif
					endif

					if (associated(newPed(pedCounter)%damPointer)) then
						tmpId =  this%dictionary%getValue(newPed(pedCounter)%damPointer%originalID)
						if (tmpID /= DICT_NULL) then
							call newPed(tmpId)%addOffspring(newPed(pedCounter))
							newPed(pedCounter)%damPointer=> newPed(tmpId)
							if (newPed(tmpId)%nOffs == 1) then
								call this%damList%list_add(newPed(tmpId))
							endif
						endif
					endif
					if (i ==0 ) then !if object is afounder add to founder array
						call  this%founders%list_add(newPed(pedCounter))
					endif

					call newGenerationList(i)%list_add(newPed(pedCounter))
					tmpIndNode => tmpIndNode%next
				end do
			enddo
			do i = 0, this%maxGeneration
				call destroyLinkedList(this%generations(i))
			enddo
			deallocate(this%generations)

			deallocate(this%pedigree)

			this%pedigree => newPed
			allocate(this%generations(0:size(newGenerationList)))
			this%generations(0:size(newGenerationList)) = newGenerationList


			this%isSorted = 3
		end subroutine sortPedigreeAndOverwriteWithDummyAtTheTop

		!---------------------------------------------------------------------------
		!< @brief Output pedigree to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine printPedigree(this)
			class(PedigreeHolder) :: this
			integer ::i
			do i= 1, this%pedigreeSize
				print *, this%pedigree(i)%originalId, this%pedigree(i)%getIntegerVectorOfRecodedIds()
			enddo
		end subroutine printPedigree


		!---------------------------------------------------------------------------
		!< @brief Output pedigree to stdout in the format originalID,sireId,damId
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine printPedigreeOriginalFormat(this, filePath)
			class(PedigreeHolder) :: this
			character(len=*), optional :: filePath
			integer ::i,unit

			if(present(filepath)) then
				open(newunit=unit, file=filePath, status="unknown")
			else
				unit = output_unit
			endif
			do i= 1, this%pedigreeSize
				write(unit,'(3a32)')  trim(this%pedigree(i)%originalId), trim(this%pedigree(i)%sireId),trim(this%pedigree(i)%damId)
			enddo
		end subroutine printPedigreeOriginalFormat


		!---------------------------------------------------------------------------
		!< @brief Output genders to file
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenders(this, filepath)
			class(PedigreeHolder) :: this
			character(len=*), intent(in) :: filepath
			integer ::i, unit

			open(newunit= unit, file= filepath, status="unknown")
			do i= 1, this%pedigreeSize
				write(unit,*) this%pedigree(i)%originalId, this%pedigree(i)%gender
			enddo

			close(unit)
		end subroutine writeOutGenders


		!---------------------------------------------------------------------------
		!< @brief read in gender information
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine readInGenders(this, filepath)
			class(PedigreeHolder) :: this
			character(len=*), intent(in) :: filepath

			integer ::i, tmpGender,tmp,unit
			character(len=IDLENGTH) :: tmpId
			open(newunit= unit, file= filepath, status="old")
			do i= 1, this%pedigreeSize
				read(unit,*) tmpId, tmpGender
				tmp = this%dictionary%getValue(tmpId)
				if (tmp/=DICT_NULL) then
					this%pedigree(tmp)%gender = tmpGender
				else
					write(error_unit,*) "WARNING: Gender info exists for animal not in the pedigree"
				endif
			enddo
		end subroutine readInGenders



		!---------------------------------------------------------------------------
		!< @brief Output  of animals that are genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypes(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt
			open(newUnit=fileUnit,file=filename,status="unknown")
			write(fmt, '(a,i10,a)') '(a20,',this%nsnpsPopulation, 'i2)'
			do i= 1, this%nGenotyped
				write(fileUnit,fmt)  this%pedigree(this%genotypeMap(i))%originalId, this%pedigree(this%genotypeMap(i))%individualGenotype%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine writeOutGenotypes



		!---------------------------------------------------------------------------
		!< @brief Output  of animals with genotype info passed in from array.
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypesArray(this, filename,  array, toPrint)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer, intent(in) :: toPrint !< 1 is all, 2 is genotyped animals, 3 is hd animals
			integer ::i, fileUnit
			integer, dimension(:,:) ,allocatable, intent(in) :: array
			character(len=100) :: fmt
			open(newUnit=fileUnit,file=filename,status="unknown")
			write(fmt, '(a,i10,a)') '(a20,',this%nsnpsPopulation, 'i2)'

			select case (toPrint)


			case(1)
				do i= 1, this%pedigreesize
					write(fileUnit,fmt)  this%pedigree(i)%originalId, array(i,:)
				enddo

			case(2)
				do i= 1, this%nGenotyped
					write(fileUnit,fmt)  this%pedigree(this%genotypeMap(i))%originalId, array(i,:)
				enddo
			case(3)
				do i= 1, this%nHd
					write(fileUnit,fmt)  this%pedigree(this%hdMap(i))%originalId, array(i,:)
				enddo
			end select


			close(fileUnit)
		end subroutine writeOutGenotypesArray


		!---------------------------------------------------------------------------
		!< @brief Output genotypes to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
		!< for all animals not just the ones that are genotyped\
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypesAll(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt

			open(newUnit=fileUnit,file=filename,status="unknown")
			write(fmt, '(a,i10,a)') '(a20,',this%nsnpsPopulation, 'i2)'
			do i= 1, this%pedigreeSize
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualGenotype%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine writeOutGenotypesAll

		!---------------------------------------------------------------------------
		!< @brief Output genotypes to stdout in the format originalID,recodedID,recodedSireID,recodedDamID
		!< for all animals not just the ones that are genotyped\
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeOutGenotypesNoDummies(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt

			open(newUnit=fileUnit,file=filename,status="unknown")
			write(fmt, '(a,i10,a)') '(a20,',this%nsnpsPopulation, 'i2)'
			do i= 1, this%pedigreeSize
				if (this%pedigree(i)%isdummy) cycle
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualGenotype%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine writeOutGenotypesNoDummies


		!---------------------------------------------------------------------------
		!< @brief Outputs phase to file of only animals that are genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine WriteoutPhase(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt

			write(fmt, '(a,i10,a)') '(a20,', this%nsnpsPopulation, 'i2)'
			open(newUnit=fileUnit,file=filename,status="unknown")
			do i= 1, this%nGenotyped
				write(fileUnit,fmt)  this%pedigree(this%genotypeMap(i))%originalId, this%pedigree(this%genotypeMap(i))%individualPhase(1)%toIntegerArray()
				write(fileUnit,fmt)  this%pedigree(this%genotypeMap(i))%originalId, this%pedigree(this%genotypeMap(i))%individualPhase(2)%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine WriteoutPhase

		!---------------------------------------------------------------------------
		!< @brief Outputs phase to file
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine WriteoutPhaseAll(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt

			write(fmt, '(a,i10,a)') '(a20,', this%nsnpsPopulation, 'i2)'
			open(newUnit=fileUnit,file=filename,status="unknown")
			do i= 1, this%pedigreeSize
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(1)%toIntegerArray()
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(2)%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine WriteoutPhaseAll

		!---------------------------------------------------------------------------
		!< @brief Outputs phase to file, excluding dummy animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine WriteoutPhaseNoDummies(this, filename)
			class(PedigreeHolder) :: this
			character(*), intent(in) :: filename
			integer ::i, fileUnit
			character(len=100) :: fmt

			write(fmt, '(a,i10,a)') '(a20,', this%nsnpsPopulation, 'i2)'
			open(newUnit=fileUnit,file=filename,status="unknown")
			do i= 1, this%pedigreeSize
				if (this%pedigree(i)%isDummy) cycle
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(1)%toIntegerArray()
				write(fileUnit,fmt)  this%pedigree(i)%originalId, this%pedigree(i)%individualPhase(2)%toIntegerArray()
			enddo
			close(fileUnit)
		end subroutine WriteoutPhaseNoDummies




		!---------------------------------------------------------------------------
		!< @brief Sets generation of an individual and his children recursively
		!< @details makes assumption that both parents also exist, and that is how generation is got
		!< both parents generation has to be set for this to work
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    Febuary 17, 2016
		!< @param[in] generation (integer)
		!< @param[in] pointer to an individual
		!---------------------------------------------------------------------------
		recursive subroutine setOffspringGeneration(this, indiv)
			type(Individual),pointer, intent(inout) :: indiv
			class(pedigreeHolder):: this


			integer :: i
			if (indiv%generation /= NOGENERATIONVALUE) then !< animal has already been set so return
				return
			endif
			if (.not. indiv%founder) then
				! if the generation of both parents has been set, add one to the greater one
				if (indiv%sirePointer%generation /= NOGENERATIONVALUE .and. indiv%damPointer%generation /= NOGENERATIONVALUE) then

					indiv%generation = max(indiv%sirePointer%generation, indiv%damPointer%generation)+1

				else !otherwise, both parents have not been set so return, as animal will get checked later
					return
				endif

			else
				indiv%generation = 0
			endif

			call this%generations(indiv%generation)%list_add(indiv)
			if(indiv%generation > this%maxGeneration) then
				this%maxGeneration = indiv%generation
			endif

			if ( indiv%nOffs /= 0) then
				do i=1,indiv%nOffs

					call this%setOffspringGeneration(indiv%OffSprings(i)%p)
				enddo
			endif
		end subroutine setOffspringGeneration

		!---------------------------------------------------------------------------
		!< @brief Constructor for recodedPedigreeArray
		!< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
		!< @date   December 22, 2016
		!---------------------------------------------------------------------------
		pure subroutine initRecodedPedigreeArray(this, n)
			implicit none
			class(recodedPedigreeArray), intent(inout) :: This !< @return initialized recoded pedigree array
			integer(int32), intent(in) :: n                    !< number of individuals in pedigree

			this%nInd = n

			if (allocated(this%originalId)) then
				deallocate(this%originalId)
			end if
			allocate(this%originalId(0:n))
			this%originalId = EMPTYID

			if (allocated(this%generation)) then
				deallocate(this%generation)
			end if
			allocate(this%generation(0:n))
			this%generation = 0

			if (allocated(this%id)) then
				deallocate(this%id)
			end if
			allocate(this%id(3, 0:n))
			this%id = 0
		end subroutine initRecodedPedigreeArray

		!---------------------------------------------------------------------------
		!< @brief Destructor for recodedPedigreeArray
		!< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
		!< @date   December 22, 2016
		!---------------------------------------------------------------------------
		pure subroutine destroyRecodedPedigreeArray(this)
			implicit none
			class(recodedPedigreeArray), intent(inout) :: this !< @return recodedPedigreeArray that will be destructed

			if (allocated(this%originalId)) then
				deallocate(this%originalId)
			end if

			if (allocated(this%generation)) then
				deallocate(this%generation)
			end if

			if (allocated(this%id)) then
				deallocate(this%id)
			end if
		end subroutine

		!---------------------------------------------------------------------------
		!< @brief Write recodedPedigreeArray to a file or stdout
		!< @author Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
		!< @date   December 22, 2016
		!---------------------------------------------------------------------------
		subroutine writeRecodedPedigreeArray(this, file)
			use iso_fortran_env, only : output_unit, int32
			implicit none
			class(recodedPedigreeArray), intent(in) :: this !< recodedPedigreeArray that will be written
			character(len=*), intent(in), optional :: file  !< If present File name, else stdout

			integer(int32) :: unit, ind
			character(len=:), allocatable :: fmt

			fmt = "(4i"//Int2Char(IDINTLENGTH)//", a1, 3a"//Int2Char(IDLENGTH)//")"
			if (present(File)) then
				open(newunit=unit, file=trim(file), status="unknown")
			else
				unit = output_unit
			end if
			do ind = 1, this%nInd
				write(unit, fmt) this%id(1:3, ind), this%generation(ind), "", this%originalId(this%id(1:3, ind))
			end do
			if (present(File)) then
				close(unit)
			end if
		end subroutine

		!---------------------------------------------------------------------------
		!< @brief Sorts and recodes pedigree
		!< @details Sorts pedigree such that parents preceede children and recodes ID to 1:n
		!< @author David Wilson david.wilson@roslin.ed.ac.uk & Gregor Gorjanc gregor.gorjanc@roslin.ed.ac.uk
		!< @date   December 20, 2016
		!---------------------------------------------------------------------------
		subroutine makeRecodedPedigreeArray(this, recPed)
			implicit none
			class(pedigreeHolder), intent(in) :: this      !< object to operate on
			type(recodedPedigreeArray), intent(out) :: recPed !< @return recoded pedigree array
			integer :: counter,i,h
			type(IndividualLinkedListNode), pointer :: tmpIndNode

			call RecPed%init(n=this%pedigreeSize)

			call this%sortPedigreeAndOverwriteWithDummyAtTheTop

			counter = 0
			do i=0, this%maxGeneration
				tmpIndNode => this%generations(i)%first
				do h=1, this%generations(i)%length
					counter = counter + 1
					recPed%originalId(counter) = tmpIndNode%item%originalId
					recPed%generation(counter) = tmpIndNode%item%generation
					recPed%id(1:3,counter) = tmpIndNode%item%getIntegerVectorOfRecodedIds()
					tmpIndNode => tmpIndNode%next
				end do
			end do
		end subroutine makeRecodedPedigreeArray


		!---------------------------------------------------------------------------
		!< @brief returns array of genotypes for all animals
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getAllGenotypesAtPosition(this, position) result(res)
			use constantModule, only : MISSINGPHASECODE
			class(pedigreeHolder) :: this
			integer, intent(in) :: position
			integer(KIND=1), allocatable, dimension(:) :: res
			integer :: counter, i
			allocate(res(this%nGenotyped))
			res = MISSINGPHASECODE
			counter = 0

			do i=1, this%nGenotyped

				counter = counter +1
				res(counter) = this%pedigree(this%genotypeMap(i))%individualGenotype%getGenotype(position)
				if (res(counter) /= 0 .and. res(counter) /= 1 .and. res(counter) /= 2 .and. res(counter) /= MISSINGGENOTYPECODE) then
					res(counter) = MISSINGGENOTYPECODE
				endif

			enddo

		end function getAllGenotypesAtPosition




		!---------------------------------------------------------------------------
		!< @brief returns list of mates and offspring for those mate pairs for given pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getPhaseAtPosition(this, position, allele) result(res)
			use constantModule, only : MISSINGPHASECODE
			class(pedigreeHolder) :: this
			integer, intent(in) :: position
			integer, intent(in) :: allele
			integer(KIND=1), allocatable, dimension(:) :: res
			integer :: counter, i
			allocate(res(this%nGenotyped))
			res = MISSINGPHASECODE
			counter = 0

			do i=1, this%nGenotyped

				counter = counter +1
				res(counter) = this%pedigree(this%genotypeMap(i))%individualPhase(allele)%getPhase(position)
				if (res(counter) /= 0 .and. res(counter) /= 1 .and. res(counter) /= 2 .and. res(counter) /= MISSINGPHASECODE) then
					res(counter) = MISSINGPHASECODE
				endif

			enddo

		end function getPhaseAtPosition


		function getPhaseAtPositionUngenotypedAnimals(this, position, allele) result(res)
			use constantModule, only : MISSINGPHASECODE
			class(pedigreeHolder) :: this
			integer, intent(in) :: position
			integer, intent(in) :: allele
			integer(KIND=1), allocatable, dimension(:) :: res
			integer :: counter, i
			allocate(res(this%pedigreeSize))
			res = MISSINGPHASECODE
			counter = 0

			do i=1, this%pedigreeSize

				counter = counter +1
				res(i) = this%pedigree(i)%individualPhase(allele)%getPhase(position)
				if (res(i) /= 0 .and. res(i) /= 1 .and. res(i) /= 2 .and. res(i) /= MISSINGPHASECODE) then
					res(i) = MISSINGPHASECODE
				endif

			enddo

		end function getPhaseAtPositionUngenotypedAnimals



		!---------------------------------------------------------------------------
		!< @brief returns list of mates and offspring for those mate pairs for given pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getAllGenotypesAtPositionWithUngenotypedAnimals(this, position) result(res)
			use constantModule, only : MISSINGPHASECODE
			class(pedigreeHolder) :: this
			integer, intent(in) :: position
			integer(KIND=1), allocatable, dimension(:) :: res
			integer :: i
			allocate(res(this%pedigreeSize))
			res = MISSINGPHASECODE

			do i=1, this%pedigreeSize

				res(i) = this%pedigree(i)%individualGenotype%getGenotype(position)
				if (res(i) /= 0 .and. res(i) /= 1 .and. res(i) /= 2 .and. res(i) /= MISSINGGENOTYPECODE) then
					res(i) = MISSINGGENOTYPECODE
				endif

			enddo

		end function getAllGenotypesAtPositionWithUngenotypedAnimals

		!---------------------------------------------------------------------------
		!< @brief returns array of what percentages an animal has been genotyped
		!<
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getGenotypePercentage(this) result(res)
			use constantModule, only : MISSINGPHASECODE
			class(pedigreeHolder) :: this
			real(KIND=real64), allocatable, dimension(:) :: res
			integer(kind=1), allocatable, dimension(:) :: indGenotypeArray
			logical, dimension(:), allocatable :: genotypedAtMarker
			integer :: i

			allocate(res(this%pedigreeSize))
			res = 0

			do i=1, this%pedigreeSize

				if (this%pedigree(i)%isGenotyped()) then

					! print *, "genotyped", indGenotypeArray
					indGenotypeArray = this%pedigree(i)%individualGenotype%toIntegerArray()
					genotypedAtMarker = ((indGenotypeArray == 0 .or. indGenotypeArray == 1) .or. indGenotypeArray == 2)
					res(i) = count(genotypedAtMarker)*1d0/size(genotypedAtMarker)
				endif
			enddo

		end function getGenotypePercentage



		!---------------------------------------------------------------------------
		!< @brief returns array of genotype information as is used by alphaimpute in format (0:nGenotyped, nSnp)
		!<
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getGenotypesAsArray(this) result(res)

			class(pedigreeHolder) :: this
			integer(kind=1) ,dimension(:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
			integer :: i


			allocate(res(this%nGenotyped, this%pedigree(this%genotypeMap(1))%individualGenotype%length))
			res = 9
			do i=1, this%nGenotyped
				res(i,:) = this%pedigree(this%genotypeMap(i))%individualGenotype%toIntegerArray()
			enddo

		end function getGenotypesAsArray


		!---------------------------------------------------------------------------
		!< @brief returns array of phase information as is used by alphaimpute in format (0:pedSized, nSnp)
		!< only information is populated where animals have been set as genotyped
		!< This takes the genotype info even if an animal is not genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getPhaseAsArray(this) result(res)

			class(pedigreeHolder) :: this
			integer(kind=1) ,dimension(:,:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
			integer :: i


			allocate(res(this%nGenotyped, this%pedigree(this%genotypeMap(1))%individualGenotype%length,2))
			res = 9
			do i=1, this%nGenotyped
				res(i,:,1) = this%pedigree(this%genotypeMap(i))%individualPhase(1)%toIntegerArray()
				res(i,:,2) = this%pedigree(this%genotypeMap(i))%individualPhase(2)%toIntegerArray()
			enddo

		end function getPhaseAsArray



		!---------------------------------------------------------------------------
		!< @brief returns array of phase information as is used by alphaimpute in format (0:pedSized, nSnp)
		!< This takes the genotype info even if an animal is not genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getPhaseAsArrayWithMissing(this) result(res)

			class(pedigreeHolder) :: this
			integer(kind=1) ,dimension(:,:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
			integer :: i

			
			allocate(res(this%pedigreeSize, this%pedigree(1)%individualPhase(1)%length,2))
			
			res = 9
			do i=1, this%pedigreeSize

				res(i,:,1) = this%pedigree(i)%individualPhase(1)%toIntegerArray()
				res(i,:,2) = this%pedigree(i)%individualPhase(2)%toIntegerArray()
			enddo

		end function getPhaseAsArrayWithMissing


		!---------------------------------------------------------------------------
		!< @brief returns array of genotype information as is used by alphaimpute in format (0:pedSized, nSnp)
		!< This takes the genotype info even if an animal is not genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getGenotypesAsArrayWitHMissing(this) result(res)

			class(pedigreeHolder) :: this
			integer(kind=1) ,dimension(:,:), allocatable :: res !indexed from 0 for COMPATIBILITY
			integer :: i


			allocate(res(this%pedigreeSize, this%pedigree(this%genotypeMap(1))%individualGenotype%length))
			do i=1, this%pedigreeSize
				res(i,:) = this%pedigree(i)%individualGenotype%toIntegerArray()
			enddo

		end function getGenotypesAsArrayWitHMissing

		!---------------------------------------------------------------------------
		!< @brief returns integer value of number of missing genotypes accross all genotypes
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getNumGenotypesMissing(this) result(count)

			class(pedigreeHolder) :: this
			integer :: count,i

			count = 0
			!$omp parallel do reduction(+:count)
			do i=1, this%nGenotyped
				count = count + this%pedigree(this%genotypeMap(i))%individualGenotype%numMissing()
			enddo
			!$omp end parallel do
		end function getNumGenotypesMissing



		!---------------------------------------------------------------------------
		!< @brief returns list of mates and offspring for those mate pairs for given pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine getMatePairsAndOffspring(this, offSpringList, listOfParents, nMatingPairs)

			use AlphaHouseMod, only : generatePairing

			class(pedigreeHolder), intent(inout) :: this      !< Pedigree object
			integer, dimension(:, :), allocatable, intent(out) :: listOfParents !< indexed by (sire/dam, mateID) = recodedId
			integer, intent(out) :: nMatingPairs
			type(IndividualLinkedList),allocatable, dimension(:) :: offspringList !< list off spring based on index of parents mateID

			type(IndividualLinkedListNode), pointer :: tmpIndNode
			type(DictStructure) :: dictionary
			integer(kind=int64) :: tmpPairingKey
			character(len=20) :: tmpPairingKeyStr ! 19 is largest characters for int64 so 20 just to be safe [and that its round]
			integer :: i,h,j

			call dictionary%DictStructure()
			nMatingPairs = 0
			if (.not. allocated(this%generations)) then
				call this%setPedigreeGenerationsAndBuildArrays
			endif
			if (allocated(listOfParents)) then
				deallocate(listOfParents)
			endif
			allocate(listOfParents(2,this%pedigreeSize))

			if (allocated(offspringList)) then
				deallocate(offspringList)
			endif
			allocate(offspringList(this%pedigreeSize))

			do i=0,this%maxGeneration
				tmpIndNode => this%generations(i)%first
				do h=1, this%generations(i)%length

					do j=1, tmpIndNode%item%nOffs
						if(associated(tmpIndNode%item,tmpIndNode%item%offsprings(j)%p%sirePointer)) then


							tmpPairingKey = generatePairing(tmpIndNode%item%offsprings(j)%p%sirePointer%id, tmpIndNode%item%offsprings(j)%p%damPointer%id)
							write(tmpPairingKeyStr, '(i0)') tmpPairingKey
							if (.not. dictionary%hasKey(tmpPairingKeyStr)) then
								nMatingPairs = nMatingPairs + 1
								listOfParents(1,nMatingPairs) = tmpIndNode%item%id
								listOfParents(2,nMatingPairs) = tmpIndNode%item%offsprings(j)%p%damPointer%id
								call offspringList(nMatingPairs)%list_add(tmpIndNode%item%offsprings(j)%p)
								call dictionary%addKey(tmpPairingKeyStr, nMatingPairs)
							else
								call offspringList(dictionary%getValue(tmpPairingKeyStr))%list_add(tmpIndNode%item%offsprings(j)%p)
							endif
						endif
					enddo
					tmpIndNode => tmpIndNode%next

				enddo
			enddo

			! call dictionary%destroy()

		end subroutine getMatePairsAndOffspring


		!---------------------------------------------------------------------------
		!< @brief Sets the individual to be genotyped.
		!<If geno array is not given, animal will still be set to genotyped. It is up to the callee
		!<if the animal has enough snps set to actually genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setAnimalAsGenotyped(this, individualIndex, geno, lockIn)

			class(pedigreeHolder) :: this
			integer, intent(in) :: individualIndex !< index of animal to get genotyped
			integer(KIND=1), dimension(:),optional, intent(in) :: geno !< One dimensional array of genotype information
			logical, intent(in), optional :: lockIn
			logical :: lock

			if (present(lockIn)) then
				lock = lockIn
			else
				lock = .false.
			endif

			if (this%nGenotyped == 0) then
				allocate(this%genotypeDictionary)
				call this%genotypeDictionary%DictStructure()
				allocate(this%genotypeMap(this%pedigreeSize))
				this%genotypeMap = 0

			else if (this%nGenotyped > this%pedigreeSize) then
				! Following error should never appear
				write(error_unit,*) "Error: animals being genotyped that are bigger than ped structure size!"
			else if (this%genotypeDictionary%getValue(this%pedigree(individualIndex)%originalID) /= DICT_NULL) then
				! if animal has already been genotyped, overwrite array, but don't increment
				if (present(geno)) then
					call this%pedigree(individualIndex)%setGenotypeArray(geno,lock)
				endif
				return
			endif

			this%nGenotyped = this%nGenotyped+1
			call this%genotypeDictionary%addKey(this%pedigree(individualIndex)%originalID, this%nGenotyped)
			if (present(geno)) then
				call this%pedigree(individualIndex)%setGenotypeArray(geno,lock)
			endif
			this%genotypeMap(this%nGenotyped) = individualIndex

		end subroutine setAnimalAsGenotyped


		!---------------------------------------------------------------------------
		!< @brief Sets the individual to be genotyped.
		!<If geno array is not given, animal will still be set to genotyped. It is up to the callee
		!<if the animal has enough snps set to actually genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setAnimalAsGenotypedFromPhase(this, individualIndex, phase1,phase2, lockIn)
			use HaplotypeModule
			use GenotypeModule
			class(pedigreeHolder) :: this
			integer, intent(in) :: individualIndex !< index of animal to get genotyped
			integer(KIND=1), dimension(:), intent(in) :: phase1,phase2 !< One dimensional array of genotype information
			logical, intent(in), optional :: lockIn
			logical :: lock
			type(Haplotype) :: h1, h2
			type(Genotype) :: geno

			call h1%Haplotype(phase1)
			call h2%Haplotype(phase2)
			call geno%Genotype(h1,h2)
			if (present(lockIn)) then
				lock = lockIn
			else
				lock = .false.
			endif

			if (this%nGenotyped == 0) then
				allocate(this%genotypeDictionary)
				call this%genotypeDictionary%DictStructure()
				allocate(this%genotypeMap(this%pedigreeSize))
				this%genotypeMap = 0

			else if (this%nGenotyped > this%pedigreeSize) then
				! Following error should never appear
				write(error_unit,*) "Error: animals being genotyped that are bigger than ped structure size!"
			else if (this%genotypeDictionary%getValue(this%pedigree(individualIndex)%originalID) /= DICT_NULL) then
				! if animal has already been genotyped, overwrite array, but don't increment
				call this%pedigree(individualIndex)%setGenotypeObject(geno)
				
				this%pedigree(individualIndex)%individualPhase(1) = h1
				this%pedigree(individualIndex)%individualPhase(2) = h2
				return
			endif

			this%nGenotyped = this%nGenotyped+1
			call this%genotypeDictionary%addKey(this%pedigree(individualIndex)%originalID, this%nGenotyped)
			call this%pedigree(individualIndex)%setGenotypeObject(geno)
			this%pedigree(individualIndex)%individualPhase(1) = h1
			this%pedigree(individualIndex)%individualPhase(2) = h2
			this%genotypeMap(this%nGenotyped) = individualIndex

		end subroutine setAnimalAsGenotypedFromPhase



		!---------------------------------------------------------------------------
		!< @brief Sets the individual to be genotyped.
		!<If geno array is not given, animal will still be set to genotyped. It is up to the callee
		!<if the animal has enough snps set to actually genotyped
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setAnimalAsGenotypedSequence(this, individualIndex, geno, referAllele, alterAllele)

			class(pedigreeHolder) :: this
			integer, intent(in) :: individualIndex !< index of animal to get genotyped
			integer(KIND=1), dimension(:),optional, intent(in) :: geno !< One dimensional array of genotype information
			integer, dimension(:), intent(in) :: referAllele, alterAllele
			if (this%nGenotyped == 0) then
				if (.not. allocated(this%genotypeDictionary)) then
					allocate(this%genotypeDictionary)
				endif
				call this%genotypeDictionary%DictStructure()
				allocate(this%genotypeMap(this%pedigreeSize))

			else if (this%nGenotyped > this%pedigreeSize) then
				! Following error should never appear
				write(error_unit,*) "Error: animals being genotyped that are bigger than ped structure size!"
			else if (this%genotypeDictionary%getValue(this%pedigree(individualIndex)%originalID) /= DICT_NULL) then
				! if animal has already been genotyped, overwrite array, but don't increment
				if (present(geno)) then
					call this%pedigree(individualIndex)%setGenotypeArray(geno)
				endif
				return
			endif

			this%nGenotyped = this%nGenotyped+1
			call this%genotypeDictionary%addKey(this%pedigree(individualIndex)%originalID, this%nGenotyped)
			if (present(geno)) then
				call this%pedigree(individualIndex)%setGenotypeArray(geno)
			endif

			call this%pedigree(individualIndex)%setSequenceArray(referAllele,alterAllele)
			this%genotypeMap(this%nGenotyped) = individualIndex

		end subroutine setAnimalAsGenotypedSequence


		!---------------------------------------------------------------------------
		!< @brief  Converts the sequence data counts to a 3D array (similar to phase) of size of pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function convertSequenceDataToArray(this) result(res)

			class(PedigreeHolder), intent(in) :: this
			integer, dimension(:,:,:),allocatable :: res
			integer :: i

			do i =1, this%pedigreeSize

				if (.not. allocated(this%pedigree(i)%referAllele)) cycle

				if(.not. allocated(res)) then
					allocate(res(this%pedigreesize,size(this%pedigree(i)%referAllele),2))
					res = 0
				endif
				res(i,:,1) = this%pedigree(i)%referAllele
				res(i,:,2) = this%pedigree(i)%alterAllele
			enddo

		end function convertSequenceDataToArray

		function getSequenceAsArrayWithMissing(this, index) result(res)

			class(PedigreeHolder), intent(in) :: this
			integer, dimension(:,:),allocatable :: res
			integer :: i, index

			if(.not. allocated(res)) then
				allocate(res(this%pedigreesize,2))
			endif

			res = 0
			do i =1, this%pedigreeSize
				if (.not. allocated(this%pedigree(i)%referAllele)) cycle
				res(i,1) = this%pedigree(i)%referAllele(index)
				res(i,2) = this%pedigree(i)%alterAllele(index)
			enddo

		end function getSequenceAsArrayWithMissing

		!---------------------------------------------------------------------------
		!< @brief Returns either the individuals id, the sires id or dams id based on
		!<which index is passed.

		!<THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function getSireDamGenotypeIDByIndex(this,ind, index) result(v)
			use iso_fortran_env, only : ERROR_UNIT
			class(PedigreeHolder), intent(in) :: this
			type(Individual), intent(in) :: ind
			character(len=IDLENGTH) :: tmp
			integer, intent(in) :: index !< index of geno index to return (1 for this, 2 for sire, 3 for dam)
			integer:: v

			v = 0
			select case (index)
			case(1)
				tmp = ind%originalId
				v = this%genotypeDictionary%getValue(tmp)
				if (v == DICT_NULL) then
					v = 0
				endif
			case(2)
				if (associated(ind%sirePointer)) then
					tmp = ind%sirePointer%originalId
					v = this%genotypeDictionary%getValue(tmp)
					if (v == DICT_NULL) then
						v = 0
					endif
				endif
			case(3)
				if (associated(ind%damPointer)) then
					tmp = ind%damPointer%originalId
					v = this%genotypeDictionary%getValue(tmp)
					if (v == DICT_NULL) then
						v = 0
					endif
				endif
				case default
				write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
			end select
			return
		end function getSireDamGenotypeIDByIndex


		!---------------------------------------------------------------------------
		!< @brief Sets the individual to be genotyped at high density.
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine setAnimalAsHD(this, indId)
			use iso_fortran_env
			class(PedigreeHolder) :: this
			integer, intent(in) :: indId

			! if index not in pedigree return.
			if (indId > this%pedigreeSize .or. indId < 1) then
				write(error_unit, *) "warning - setAnimalAsHD was given an index that was out of range"
				return
			endif
			if (.not. this%pedigree(indId)%genotyped) then
				write(error_unit, *) "warning - setAnimalAsHD was given an index of animal that was not genotyped"
				write(error_unit, *) "animal has ID:", trim(this%pedigree(indId)%originalID), " and recoded ID:", indid
			endif
			if (this%nHd == 0) then
				if (.not. allocated(this%hdDictionary)) then
					allocate(this%hdDictionary)
				endif
				call this%hdDictionary%DictStructure()
				if (allocated(this%hdMap)) deallocate(this%hdMap)

				allocate(this%hdMap(this%pedigreeSize))
				this%hdMap = 0
			endif

			if (this%hdDictionary%getValue(this%pedigree(indId)%originalId) ==DICT_NULL) then
				this%nHd = this%nHd + 1
				this%pedigree(indId)%hd = .true.

				this%hdMap(this%nHd) = indId
				call this%hdDictionary%addKey(this%pedigree(indId)%originalId, this%nHd)
			else
				this%pedigree(indId)%hd = .true.
			endif

		end subroutine setAnimalAsHD



		!---------------------------------------------------------------------------
		!< @brief Returns either the individuals id in hd index, the sires id or dams id based on
		!<which index is passed.
		!<THIS IS DEPRECATED - ONLY MEANT FOR COMPATIBILITY
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		! PARAMETERS:
		!< @param[in] index - the index
		!< @return hdIndex of animal based on index
		!---------------------------------------------------------------------------
		function getSireDamHDIDByIndex(this,ind, index) result(v)
			use iso_fortran_env, only : ERROR_UNIT
			class(PedigreeHolder), intent(in) :: this
			type(Individual), intent(in) :: ind
			character(len=IDLENGTH) :: tmp
			integer, intent(in) :: index !< index of hd index to return (1 for this, 2 for sire, 3 for dam)
			integer:: v

			v = 0
			select case (index)
			case(1)
				tmp = ind%originalId
				v = this%hdDictionary%getValue(tmp)
				if (v == DICT_NULL) then
					v = 0
				endif
			case(2)
				if (associated(ind%sirePointer)) then
					tmp = ind%sirePointer%originalId
					v = this%hdDictionary%getValue(tmp)
					if (v == DICT_NULL) then
						v = 0
					endif
				endif
			case(3)
				if (associated(ind%damPointer)) then
					tmp = ind%damPointer%originalId
					v = this%hdDictionary%getValue(tmp)
					if (v == DICT_NULL) then
						v = 0
					endif
				endif
				case default
				write(error_unit, *) "error: getSireDamByIndex has been given an out of range value"
			end select
			return
		end function getSireDamHDIDByIndex


		!---------------------------------------------------------------------------
		!< @brief creates a new dummy animal at end of pedigree
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		! PARAMETERS:
		!---------------------------------------------------------------------------
		subroutine createDummyAnimalAtEndOfPedigree(this,dummyId, offspringId, sireIn)
			use IFCORE
			class(PedigreeHolder) :: this
			integer, intent(out) :: dummyId
			integer, optional :: offspringId !< offspring recoded id canbe given here
			logical, optional :: sireIn !< if true, assign to sire location
			character(len=IDLENGTH) :: tmpCounterStr

			
			this%pedigreeSize = this%pedigreeSize+1

			if (this%pedigreeSize > this%maxPedigreeSize) then
				write(error_unit,*) "ERROR: too many undefined animals"
				call TRACEBACKQQ(string= "ERROR: too many undefined animals",user_exit_code=1)

			endif

			this%nDummys = this%nDummys + 1

			this%isSorted = 0
			tmpCounterStr = ""
			write(tmpCounterStr, '(I4.4)') this%nDummys
			call this%Pedigree(this%pedigreeSize)%initIndividual(trim(dummyAnimalPrepre)//trim(tmpCounterStr) ,'0','0', this%pedigreeSize,nsnps=this%nsnpsPopulation)
			call this%dictionary%addKey(trim(dummyAnimalPrepre)//trim(tmpCounterStr), this%pedigreeSize)
			this%Pedigree(this%pedigreeSize)%isDummy = .true.
			call this%Founders%list_add(this%Pedigree(this%pedigreeSize))
			this%Pedigree(this%pedigreeSize)%founder = .true.

			if (present(offspringId)) then
				if (offspringId > this%pedigreeSize) then
					write(error_unit,*) "ERROR - dummy list given index larger than pedigree"
				endif

				call this%Pedigree(this%pedigreeSize)%AddOffspring(this%pedigree(offspringId))

				if (present(sireIn)) then
					if (sireIn) then
						this%pedigree(offspringId)%sirePointer => this%Pedigree(this%pedigreeSize)
						call this%sireList%list_add(this%Pedigree(this%pedigreeSize))
						call this%Pedigree(this%pedigreeSize)%setGender(1)
					else
						this%pedigree(offspringId)%damPointer => this%Pedigree(this%pedigreeSize)
						call this%damList%list_add(this%Pedigree(this%pedigreeSize))
						call this%Pedigree(this%pedigreeSize)%setGender(2)
					endif
				else 
					if (.not. associated(this%pedigree(offspringId)%sirePointer)) then
						this%pedigree(offspringId)%sirePointer => this%Pedigree(this%pedigreeSize)
						call this%sireList%list_add(this%Pedigree(this%pedigreeSize))
						call this%Pedigree(this%pedigreeSize)%setGender(1)
					else if (.not. associated(this%pedigree(offspringId)%damPointer)) then
						this%pedigree(offspringId)%damPointer => this%Pedigree(this%pedigreeSize)
						call this%damList%list_add(this%Pedigree(this%pedigreeSize))
						call this%Pedigree(this%pedigreeSize)%setGender(2)
					else
						write(error_unit,*) "ERROR - dummy animal given offspring that already has both parents!"
					end if
				endif
			endif
			dummyId = this%pedigreeSize
		end subroutine createDummyAnimalAtEndOfPedigree




		!---------------------------------------------------------------------------
		!< @brief creates a new animal at end of pedigree
		!< If genotype is supplied, animal is set to hd
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine addAnimalAtEndOfPedigree(this, originalID, geno, offspringID)
			use IFCORE
			class(PedigreeHolder) :: this
			character(len=*) ,intent(in):: OriginalId
			integer(kind=1), dimension(:), intent(in), optional :: geno
			integer, intent(in), optional :: offspringId
			! change pedigree to no longer be sorted

			this%issorted = 0

			this%pedigreeSize = this%pedigreeSize+1
			this%addedRealAnimals = this%addedRealAnimals + 1

			if (this%pedigreeSize > this%maxPedigreeSize) then
				write(error_unit,*) "ERROR: too many undefined animals"
				call TRACEBACKQQ(string= "ERROR: too many undefined animals",user_exit_code=1)

			endif
			call this%Pedigree(this%pedigreeSize)%initIndividual(OriginalId ,'0','0', this%pedigreeSize,nsnps=this%nsnpsPopulation)
			call this%dictionary%addKey(OriginalId, this%pedigreeSize)
			this%Pedigree(this%pedigreeSize)%isDummy = .false.
			this%Pedigree(this%pedigreeSize)%originalPosition = this%addedRealAnimals
			this%inputMap(this%pedigreeSize) = this%addedRealAnimals

			call this%Founders%list_add(this%Pedigree(this%pedigreeSize))
			this%Pedigree(this%pedigreeSize)%founder = .true.
			
			if (present(offspringId)) then
				if (offspringId > this%pedigreeSize) then
					write(error_unit,*) "ERROR - dummy list given index larger than pedigree"
				endif

				call this%Pedigree(this%pedigreeSize)%AddOffspring(this%pedigree(offspringId))
				if (this%pedigree(offspringId)%sireId == originalID) then
					this%pedigree(offspringId)%sirePointer => this%Pedigree(this%pedigreeSize)
					call this%sireList%list_add(this%Pedigree(this%pedigreeSize))
					call this%Pedigree(this%pedigreeSize)%setGender(1)
				else if (this%pedigree(offspringId)%damId == originalID) then
					this%pedigree(offspringId)%damPointer => this%Pedigree(this%pedigreeSize)
					call this%damList%list_add(this%Pedigree(this%pedigreeSize))
					call this%Pedigree(this%pedigreeSize)%setGender(2)
				else
					write(error_unit,*) "ERROR - new animal given offspring which doesn't share its ID!!"
				end if
			endif

			if (present(geno)) then
				call this%setAnimalAsGenotyped(this%pedigreeSize, geno)
				call this%pedigree(this%pedigreeSize)%initPhaseArrays(size(geno))
			else if (this%nsnpsPopulation /=0) then
				call this%pedigree(this%pedigreeSize)%initPhaseAndGenotypes(this%nsnpsPopulation)
			endif
		end subroutine addAnimalAtEndOfPedigree



		!---------------------------------------------------------------------------
		!> @brief Sets all individuals genotype from the Haplotype
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!---------------------------------------------------------------------------
		SUBROUTINE MakeGenotype(this)
			class(PedigreeHolder) :: this
			! Any individual that has a missing genotype information but has both alleles
			! known, has its genotype filled in as the sum of the two alleles
			integer :: i

			!$!OMP PARALLEL DO &
			!$!OMP PRIVATE(i)
			do i=1,this%pedigreeSize
				call this%pedigree(i)%makeIndividualGenotypeFromPhase()
			enddo
			!$!OMP END PARALLEL DO


		END SUBROUTINE MakeGenotype

		!#############################################################################################################################################################################################################################



		!---------------------------------------------------------------------------
		!> @brief Sets the individual haplotypes from the compilement if animal is genotyped for all animals
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine PhaseComplement(this)
			class(PedigreeHolder) :: this
			! If the genotype at a locus for an individual is known and one of its alleles has been determined
			! then impute the missing allele as the complement of the genotype and the known phased allele
			integer :: i

			!$!OMP PARALLEL DO &
			!$!OMP PRIVATE(i)
			do i=1,this%pedigreeSize
				call this%pedigree(i)%makeIndividualPhaseCompliment()
			enddo
			!$!OMP END PARALLEL DO

		end subroutine PhaseComplement



		!---------------------------------------------------------------------------
		!> @brief Sets the individual genotypes from the haplotypes. Overwrites anything that was already there in the genotype
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine cleanGenotypesBasedOnHaplotypes(this)
			use GenotypeModule
			class(PedigreeHolder) :: this
			integer :: i

			!$OMP PARALLEL DO &
			!$OMP PRIVATE(i)
			do i=1,this%pedigreeSize
				call this%pedigree(i)%individualGenotype%Genotype(this%pedigree(i)%individualPhase(1), this%pedigree(i)%individualPhase(2))
			enddo
			!$OMP END PARALLEL DO

		end subroutine cleanGenotypesBasedOnHaplotypes

		!---------------------------------------------------------------------------
		!< @brief  counts missing snps across all animals at every snp
		!<Returns a count of missing snps across all animals at every snp
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2016
		!---------------------------------------------------------------------------
		function countMissingGenotypesNoDummys(this) result(res)
			integer :: res
			integer :: i
			class(PedigreeHolder) :: this

			res = 0
			do i=1, this%pedigreeSize

				if (this%pedigree(i)%isDummy) cycle

				res = res + this%pedigree(i)%individualGenotype%numMissing()
			end do
		end function countMissingGenotypesNoDummys

		!---------------------------------------------------------------------------
		!< @brief  counts missing alleles (in the 2 gametes) across all animals at every snp
		!<Returns a count of missing alleles across all animals at every snp
		!< @author  Mara Battagin mara.battagin@roslin.ed.ac.uk
		!< @date    August 31, 2017
		!---------------------------------------------------------------------------
		function countMissingPhaseNoDummys(this) result(res)
			integer :: res
			integer :: i
			class(PedigreeHolder) :: this

			res = 0

			do i=1, this%pedigreeSize

				if (this%pedigree(i)%isDummy) cycle

				res = res + this%pedigree(i)%individualPhase(1)%numberMissing()
				res = res + this%pedigree(i)%individualPhase(2)%numberMissing()
			end do
		end function countMissingPhaseNoDummys


		!---------------------------------------------------------------------------
		!< @brief Sets the phase for homozygotic snps
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		function calculatePedigreeCorrelationWithInbreeding(pedIn, additVarianceIn) result (values)
			use SortedIntegerLinkedListModule
			class(PedigreeHolder), intent(in):: pedIn
			real(real64), intent(in), optional:: additVarianceIn
			real(real64), dimension(:,:), allocatable:: values
			type(sortedIntegerLinkedList), allocatable:: sireList, damList
			real(real64), dimension(:,:), allocatable:: lValues
			real(real64), dimension(:), allocatable:: F

			integer:: knownDummies
			integer:: i, youngestSire, youngestDam
			integer:: sireI, damI, temp, ID

			logical:: sireKnown, damKnown
			logical:: sireJKnown, damJKnown
			real(real64):: D
			integer:: numLevels

			numLevels = pedIn%pedigreeSize-pedIn%UnknownDummys

			allocate(values(numLevels, numLevels))
			allocate(F(0:numLevels))

			F(0) = -1
			F(1:) = 0

			knownDummies = pedIn%nDummys - pedIn%unKnownDummys

			values = 0

			!$OMP PARALLEL DO PRIVATE(lValues, i, damKnown, sireKnown, damJKnown, sireJKnown, sireI, damI, ID, damList, sireList, temp, youngestSire, youngestDam, D)
			do i = 1, numLevels
				allocate(damList)
				allocate(sireList)
				allocate(lValues(numLevels, numLevels))
				damKnown = .false.
				sireKnown = .false.
				lValues = 0
				F(i) = 0
				sireI = 0
				damI = 0
				ID = pedIn%pedigree(i)%ID
				if (associated(pedIn%pedigree(i)%sirePointer)) then
					if (.not. pedIn%pedigree(i)%sirePointer%isUnknownDummy) then
						sireI = pedIn%pedigree(i)%sirePointer%ID
						call sireList%list_add(pedIn%pedigree(i)%sirePointer%ID)
						lValues(sireI, sireI) = 1
						sireKnown = .true.
					end if
				end if

				if (associated(pedIn%pedigree(i)%damPointer)) then
					if (.not. pedIn%pedigree(i)%damPointer%isUnknownDummy) then
						damI = pedIn%pedigree(i)%damPointer%ID
						call damList%list_add(damI)
						lValues(damI, damI) = 1
						damKnown=.true.
					end if
				end if

				do while (sireList%length>0 .and. damList%length>0)
					youngestSire = sireList%first%item
					youngestDam = damList%first%item

					if (youngestSire>youngestDam) then !i.e. the sire is younger
						call addSireDamToListAndUpdateValues(sireList, pedIn%pedigree(youngestSire), lValues, sireI)
						temp= sireList%pop_first()
					else if (youngestDam>youngestSire) then !i.e. the dam is younger
						call addSireDamToListAndUpdateValues(damList, pedIn%pedigree(youngestDam), lValues, damI)
						temp = damList%pop_first()
					else !youngestSire == youngestDam
						sireJKnown = .false.
						damJKnown = .false.
						if (associated(pedIn%pedigree(youngestSire)%sirePointer)) then
							if (pedIn%pedigree(youngestSire)%sirePointer%isUnknownDummy) then
								sireJKnown = .true.
							end if
						end if

						if (associated(pedIn%pedigree(youngestSire)%damPointer)) then
							if (pedIn%pedigree(youngestSire)%damPointer%isUnknownDummy) then
								damJKnown = .true.
							end if
						end if

						if (sireJKnown .and. damJKnown) then
							D= 2
						else if (sireJKnown .or. damJKnown) then
							D = 4.0_real64/3.0_real64
						else
							D = 1.0_real64
						end if

						call addSireDamToListAndUpdateValues(damList, pedIn%pedigree(youngestSire), lValues, sireI)
						call addSireDamToListAndUpdateValues(sireList, pedIn%pedigree(youngestSire), lValues, damI)
						!$OMP ATOMIC
						F(i) = F(i) +lValues(sireI, youngestSire)*lValues(damI, youngestSire)*0.5*D
						temp = sireList%pop_first()
						temp = damList%pop_first()
					end if
				end do

				if (sireKnown .and. damKnown) then
					D = 0.5_real64 -0.25_real64*(F(sireI)+F(damI))
				else if (sireKnown) then
					D = 0.75_real64 -0.25_real64*F(sireI)
				else if (damKnown) then
					D = 0.75_real64 -0.25_real64*F(damI)
				else
					D = 1_real64
				end if

				D = 1.0_real64/D

				!$OMP ATOMIC
				values(ID, ID) = values(ID, ID) + D
				if (sireKnown) then
					!$OMP ATOMIC
					values(sireI, ID) = values(sireI, ID) -0.5_real64*D
					!$OMP ATOMIC
					values(ID, sireI) = values(ID, sireI) -0.5_real64*D
					!$OMP ATOMIC
					values(sireI, sireI) = values(sireI, sireI) +0.25_real64*D
				end if

				if (damKnown) then
					!$OMP ATOMIC
					values(damI, ID) = values(damI, ID) - 0.5_real64*D
					!$OMP ATOMIC
					values(ID, damI) = values(ID, damI) - 0.5_real64*D
					!$OMP ATOMIC
					values(damI, damI) = values(damI, damI) + 0.25_real64*D
				end if

				if (damKnown .and. sireKnown) then
					!$OMP ATOMIC
					values(sireI, damI) = values(sireI, damI) + 0.25_real64*D
					!$OMP ATOMIC
					values(damI, sireI) = values(damI, sireI) + 0.25_real64*D
				end if
				deallocate(damList)
				deallocate(sireList)
				deallocate(lValues)
			end do
			!    !$OMP END PARALLEL DO

			if (present(additVarianceIn)) then
				values = values*additVarianceIn
			end if

		end function calculatePedigreeCorrelationWithInbreeding

#ifdef MPIACTIVE
		function calculatePedigreeCorrelationWithInBreedingMPI(pedIn, additVarianceIn, communicatorIn) result(values)
			use SortedIntegerLinkedListModule
			use mpi
			use MPIUtilities, only: checkMPI
			class(PedigreeHolder), intent(in):: pedIn
			real(real64), intent(in), optional:: additVarianceIn
			integer, intent(in), optional:: communicatorIn
			type(sortedIntegerLinkedList), allocatable:: sireList, damList
			real(real64), dimension(:,:), allocatable:: lValues
			real(real64), dimension(:), allocatable:: F, F2


			real(real64), dimension(:, :), allocatable:: values !< symetric matrix that is returned is size of animals (with no UNKNOWN dummys)

			integer:: knownDummies, extras, mpiErr
			integer:: i, j, youngestSire, youngestDam
			integer:: sireI, damI, temp, ID
			integer,dimension(:), allocatable:: gatherSizes, offsetLocation

			logical:: sireKnown, damKnown
			logical:: sireJKnown, damJKnown
			real(real64):: D
			type(IndividualLinkedList):: knownAnimals
			type(Individual), dimension(:), pointer:: knownAnimalArray

			integer:: numAnimalsInThisgeneration, startAnimal, endAnimal, animalsPerCore, firstGenAnimal, lastGenAnimal

			integer:: mpiSize, mpiRank, mpiCommunicator
			integer:: numLevels

			if (present(communicatorIn)) then
				mpiCommunicator = communicatorIn
			else
				mpiCommunicator = MPI_COMM_WORLD
			end if

			call MPI_COMM_SIZE(mpiCommunicator, mpiSize, mpiErr)
			call checkMPI(mpiErr)
			call MPI_COMM_RANK(mpiCommunicator, mpiRank, mpiErr)
			call checkMPI(mpiErr)
			numLevels = pedIn%pedigreeSize-pedIn%UnknownDummys

			allocate(values(numLevels, numLevels))

			allocate(F(0:numLevels))
			allocate(gatherSizes(mpiSize))
			allocate(offsetLocation(mpiSize))

			F(0) = -1
			F(1:) = 0

			knownDummies = pedIn%nDummys - pedIn%unKnownDummys
			values = 0
			lastGenAnimal = 0
			GenerationLoop: do j = 0, size(pedIn%generations)-1
				knownAnimals = pedIn%generations(j)%convertToListOfKnownAnimals()
				knownAnimalArray => knownAnimals%convertToArray()
				numAnimalsInThisgeneration = knownAnimals%length

				animalsPerCore = numAnimalsInThisgeneration/mpiSize

				extras = numAnimalsInThisgeneration - mpiSize*(numAnimalsInThisgeneration/mpiSize) !Integer Arithmatic

				startAnimal = mpiRank*animalsPerCore+1
				endAnimal = (mpiRank+1)*animalsPerCore
				if (mpiRank ==mpiSize-1) then
					endAnimal = endAnimal+extras
				end if


				gatherSizes = animalsPerCore
				do i =1, mpiSize
					offsetLocation(i) = animalsPerCore*(i-1)
				end do
				gatherSizes(mpiSize) = gatherSizes(mpiSize)+extras

				allocate(F2(gatherSizes(mpiRank+1)))
				F2=0

				firstGenAnimal = lastGenAnimal+1
				lastGenAnimal = firstGenAnimal+numAnimalsInThisgeneration-1

				F(FirstGenAnimal: lastGenAnimal) = 0

				!$OMP PARALLEL DO PRIVATE(lValues, i, damKnown, sireKnown, damJKnown, sireJKnown, sireI, damI, ID, damList, sireList, temp, youngestSire, youngestDam, D)
				AnimalLoop: do i = startAnimal, endAnimal

					allocate(damList)
					allocate(sireList)
					allocate(lValues(numLevels, numLevels))
					damKnown = .false.
					sireKnown = .false.
					lValues = 0
					F2(i-startAnimal+1) = 0
					sireI = 0
					damI = 0
					ID = knownAnimalArray(i)%ID
					if (associated(knownAnimalArray(i)%sirePointer)) then
						if (.not. knownAnimalArray(i)%sirePointer%isUnknownDummy) then
							sireI = knownAnimalArray(i)%sirePointer%ID
							call sireList%list_add(sireI)!knownAnimalArray(i)%sirePointer%ID)
							lValues(sireI, sireI) = 1
							sireKnown = .true.
						end if
					end if

					if (associated(knownAnimalArray(i)%damPointer)) then
						if (.not. knownAnimalArray(i)%damPointer%isUnknownDummy) then
							damI = knownAnimalArray(i)%damPointer%ID
							call damList%list_add(damI)
							lValues(damI, damI) = 1
							damKnown=.true.
						end if
					end if

					do while (sireList%length>0 .and. damList%length>0)
						youngestSire = sireList%first%item
						youngestDam = damList%first%item

						if (youngestSire>youngestDam) then !i.e. the sire is younger
							call addSireDamToListAndUpdateValues(sireList, pedIn%pedigree(youngestSire), lValues, sireI)
							temp= sireList%pop_first()
						else if (youngestDam>youngestSire) then !i.e. the dam is younger
							call addSireDamToListAndUpdateValues(damList, pedIn%pedigree(youngestDam), lValues, damI)
							temp = damList%pop_first()
						else !youngestOnSireSide == youngestOnDamSide
							sireJKnown = .false.
							damJKnown = .false.
							if (associated(pedIn%pedigree(youngestSire)%sirePointer)) then
								if (pedIn%pedigree(youngestSire)%sirePointer%isUnknownDummy) then
									sireJKnown = .true.
								end if
							end if

							if (associated(pedIn%pedigree(youngestSire)%damPointer)) then
								if (pedIn%pedigree(youngestSire)%damPointer%isUnknownDummy) then
									damJKnown = .true.
								end if
							end if

							if (sireJKnown .and. damJKnown) then
								D= 2
							else if (sireJKnown .or. damJKnown) then
								D = 4.0_real64/3.0_real64
							else
								D = 1.0_real64
							end if

							call addSireDamToListAndUpdateValues(damList, pedIn%pedigree(youngestSire), lValues, sireI)
							call addSireDamToListAndUpdateValues(sireList, pedIn%pedigree(youngestSire), lValues, damI)
							F2(i-startAnimal+1) = F2(i-startAnimal+1) +lValues(sireI, youngestSire)*lValues(damI, youngestSire)*0.5*D
							temp = sireList%pop_first()
							temp = damList%pop_first()
						end if
					end do

					if (sireKnown .and. damKnown) then
						D = 0.5_real64 -0.25_real64*(F(sireI)+F(damI))
					else if (sireKnown) then
						D = 0.75_real64 -0.25_real64*F(sireI)
					else if (damKnown) then
						D = 0.75_real64 -0.25_real64*F(damI)
					else
						D = 1.0_real64
					end if

					D = 1.0_real64/D

					!$OMP ATOMIC
					values(ID, ID) = values(ID, ID) + D!1.0_real64/(1+F(i))
					if (sireKnown) then
						!$OMP ATOMIC
						values(sireI, ID) =values(sireI, ID) -0.5*D!F(i)
						!$OMP ATOMIC
						values(ID, sireI) = values(ID, sireI) -0.5*D!F(i)
						!$OMP ATOMIC
						values(sireI, sireI) = values(sireI, sireI) +0.25*D!F(i)
					end if

					if (damKnown) then
						!$OMP ATOMIC
						values(damI, ID) = values(damI, ID) -0.5*D!F(i)
						!$OMP ATOMIC
						values(ID, damI) = values(ID, damI) - 0.5*D!F(i)
						!$OMP ATOMIC
						values(damI, damI) = values(damI, damI) + 0.25*D!F(i)
					end if

					if (damKnown .and. sireKnown) then
						!$OMP ATOMIC
						values(sireI, damI) = values(sireI, damI) + 0.25*D
						!$OMP ATOMIC
						values(damI, sireI) = values(damI, sireI) + 0.25*D
					end if

					deallocate(sireList)
					deallocate(damList)
					deallocate(lValues)

				end do AnimalLoop
				!$OMP END PARALLEL DO

				call MPI_ALLGATHERV(F2, gatherSizes(mpiRank+1), MPI_DOUBLE, F(firstGenAnimal:lastGenAnimal), gatherSizes, offsetLocation, MPI_DOUBLE, mpiCommunicator, mpiErr)
				knownAnimalArray => null()
				deallocate(F2)
				!      end if
			end do GenerationLoop

			call MPI_ALLREDUCE(MPI_IN_PLACE, values, size(values), MPI_DOUBLE, MPI_SUM, mpiCommunicator, mpiErr)

			if (present(additVarianceIn)) then
				values = values*additVarianceIn
			end if
		end function calculatePedigreeCorrelationWithInBreedingMPI
#endif


		subroutine addSireDamToListAndUpdateValues(listIn, IndividualIn, values, firstValue)
			use SortedIntegerLinkedListModule
			type(sortedIntegerLinkedList), intent(inout):: listIn
			type(Individual), intent(in):: IndividualIn
			real(real64), dimension(:,:), intent(inout):: values
			integer, intent(in):: firstValue
			integer:: sireID, damID
			sireID = 0
			damID = 0
			if (associated(IndividualIn%sirePointer)) then
				if (.not. IndividualIn%sirePointer%isUnknownDummy) then
					sireID = IndividualIn%sirePointer%ID
					call listIn%list_add(SireID)
					!$OMP ATOMIC
					values(firstValue, sireID) = values(firstValue, sireID) + 0.5* values(firstValue, IndividualIn%ID)
				end if
			end if

			if (associated(IndividualIn%damPointer)) then
				if (.not. IndividualIn%damPointer%isUnknownDummy) then
					damID = IndividualIn%damPointer%ID
					call listIn%list_add(IndividualIn%damPointer%ID)
					!$OMP ATOMIC
					values(firstValue, damID) = values(firstValue, damID) + 0.5* values(firstValue, IndividualIn%ID)
				end if
			end if
		end subroutine addSireDamToListAndUpdateValues


		subroutine memoryClearer(this)
			type(PedigreeHolder),intent(inout)  :: this
			integer :: i

			!$OMP Parallel DO
			do i=1, this%pedigreeSize
				this%pedigree(i)%used  = this%pedigree(i)%used  - 1
				if (this%pedigree(i)%used <= 0) then
					call writeOutPhaseAndGenotypeBinary(this%pedigree(i))
					deallocate(this%pedigree(i)%individualGenotype)
					deallocate(this%pedigree(i)%individualPhase)
				endif
			enddo
			!$omp end parallel do

		end subroutine memoryClearer





		subroutine writeOutPhaseAndGenotypeBinary(ind)
			use constantModule, only : storageFolder
			USE IFPORT
			type(individual) :: ind
			logical :: exists, result
			integer :: unit
			inquire(file=storageFolder,EXIST=exists)

			if (.not. exists) then
				result = MAKEDIRQQ(storageFolder)
			endif

			inquire(file=storageFolder//DASH//trim(ind%originalID),EXIST=exists)
			if (.not. exists) then
				result = MAKEDIRQQ(storageFolder//DASH//trim(ind%originalID))
			endif

			inquire(file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase1",EXIST=exists)
			if (.not. exists) then
				result = MAKEDIRQQ(storageFolder//DASH//trim(ind%originalID)//DASH // "phase1")
			endif

			inquire(file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase2",EXIST=exists)
			if (.not. exists) then
				result = MAKEDIRQQ(storageFolder//DASH//trim(ind%originalID)//DASH // "phase2")
			endif

			inquire(file=storageFolder//DASH//trim(ind%originalID)//DASH // "genotype",EXIST=exists)
			if (.not. exists) then
				result = MAKEDIRQQ(storageFolder//DASH//trim(ind%originalID)//DASH // "genotype")
			endif


			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "genotype"// DASH// "genotypeFile", status="unknown", form = 'unformatted')
			write(unit) ind%individualGenotype%sections
			write(unit) ind%individualGenotype%homo
			write(unit) ind%individualGenotype%additional
			write(unit) ind%individualGenotype%hasLock
			if (ind%individualGenotype%hasLock) then
				write(unit) ind%individualGenotype%locked
			endif
			write(unit) ind%individualGenotype%overhang
			write(unit) ind%individualGenotype%length
			close(unit)


			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase1"// DASH// "phaseFile", status="unknown", form = 'unformatted')
			write(unit) ind%individualPhase(1)%sections
			write(unit) ind%individualPhase(1)%phase
			write(unit) ind%individualPhase(1)%missing
			write(unit) ind%individualPhase(1)%hasLock
			if (ind%individualPhase(1)%hasLock) then
				write(unit) ind%individualPhase(1)%locked
			endif
			write(unit) ind%individualPhase(1)%overhang
			write(unit) ind%individualPhase(1)%length
			write(unit) ind%individualPhase(1)%startPosition
			close(unit)

			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase2"// DASH// "phaseFile", status="unknown", form = 'unformatted')
			write(unit) ind%individualPhase(2)%sections
			write(unit) ind%individualPhase(2)%phase
			write(unit) ind%individualPhase(2)%missing
			write(unit) ind%individualPhase(2)%hasLock
			if (ind%individualPhase(2)%hasLock) then
				write(unit) ind%individualPhase(2)%locked
			endif
			write(unit) ind%individualPhase(2)%overhang
			write(unit) ind%individualPhase(2)%length
			write(unit) ind%individualPhase(2)%startPosition
			close(unit)

		end subroutine writeOutPhaseAndGenotypeBinary


		subroutine readInPhaseAndGenotypeBinary(ind)

			type(individual) :: ind
			integer :: unit

			allocate(ind%individualGenotype)
			allocate(ind%individualPhase(2))

			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "genotype"// DASH// "genotypeFile", status="unknown", form = 'unformatted')


			read(unit) ind%individualGenotype%sections

			allocate(ind%individualGenotype%homo(ind%individualGenotype%sections))
			allocate(ind%individualGenotype%additional(ind%individualGenotype%sections))
			allocate(ind%individualGenotype%locked(ind%individualGenotype%sections))
			read(unit) ind%individualGenotype%homo
			read(unit) ind%individualGenotype%additional
			read(unit) ind%individualGenotype%hasLock

			if (ind%individualGenotype%hasLock) then
				read(unit) ind%individualGenotype%locked
			endif
			read(unit) ind%individualGenotype%overhang
			read(unit) ind%individualGenotype%length
			close(unit)


			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase1"// DASH// "phaseFile", status="unknown", form = 'unformatted')
			read(unit) ind%individualPhase(1)%sections

			allocate(ind%individualPhase(1)%phase(ind%individualPhase(1)%sections))
			allocate(ind%individualPhase(1)%missing(ind%individualPhase(1)%sections))
			allocate(ind%individualPhase(1)%locked(ind%individualPhase(1)%sections))
			read(unit) ind%individualPhase(1)%phase
			read(unit) ind%individualPhase(1)%missing
			read(unit) ind%individualPhase(1)%hasLock

			if (ind%individualPhase(1)%hasLock) then
				read(unit) ind%individualPhase(1)%locked
			endif
			read(unit) ind%individualPhase(1)%overhang
			read(unit) ind%individualPhase(1)%length
			read(unit) ind%individualPhase(1)%startPosition
			close(unit)

			open(newunit=unit,file=storageFolder//DASH//trim(ind%originalID)//DASH // "phase2"// DASH// "phaseFile", status="unknown", form = 'unformatted')
			read(unit) ind%individualPhase(2)%sections

			allocate(ind%individualPhase(2)%phase(ind%individualPhase(1)%sections))
			allocate(ind%individualPhase(2)%missing(ind%individualPhase(1)%sections))
			allocate(ind%individualPhase(2)%locked(ind%individualPhase(1)%sections))
			read(unit) ind%individualPhase(2)%phase
			read(unit) ind%individualPhase(2)%missing
			read(unit) ind%individualPhase(2)%hasLock

			if (ind%individualPhase(2)%hasLock) then
				read(unit) ind%individualPhase(2)%locked
			endif
			read(unit) ind%individualPhase(2)%overhang
			read(unit) ind%individualPhase(2)%length
			read(unit) ind%individualPhase(2)%startPosition
			close(unit)

		end subroutine readInPhaseAndGenotypeBinary


		! this should be called in an openmp task
		! We want to call this function only for animals required
		! We also want to know if its phase or genotype they need
		subroutine memGetter(ind)
			type(individual), intent(inout) :: ind !< individual to check what memory needs got from

			if (allocated(ind%individualGenotype)) return !< if info is already there, don't read in

			! spawn new thread here - so other animal jobs can still be done on reading

			call readInPhaseAndGenotypeBinary(ind)

		end subroutine memGetter


end module PedigreeModule

















