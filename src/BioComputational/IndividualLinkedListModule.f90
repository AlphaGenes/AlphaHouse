
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     IndividualLinkedListModule.f90
!
! DESCRIPTION:
!> @brief    Module containing definition of linked List for objects (currently Individuals.
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

module IndividualLinkedListModule
	use iso_fortran_env
	use individualModule

	implicit none

	abstract interface

	! Abstract function required to allow for listing of types
	logical function logicalAbstractFunction(item)
		import :: Individual
		type(Individual),intent(in) :: item
	end function logicalAbstractFunction
end interface


type :: IndividualLinkedList
type(IndividualLinkedListNode),pointer :: first => null()
type(IndividualLinkedListNode),pointer :: last => null()
! TODO - maybe have a middle?
integer :: length = 0

contains
	procedure :: list_add
	procedure :: list_all
	procedure :: list_pop
	procedure :: list_get_nth
	procedure :: list_remove
	procedure :: contains
	procedure :: writeLinkedList
	final :: destroyLinkedList
	procedure :: convertToArray
	procedure :: convertToArrayIDs
	procedure :: getGenotypesAtPosition
	! procedure :: destroyLinkedListFinal
	procedure :: removeIndividualsBasedOnThreshold
	procedure :: convertToArrayOriginalIDs
	procedure :: convertToListOfKnownAnimals
	generic:: write(formatted)=> writeLinkedList

end type IndividualLinkedList

type :: IndividualLinkedListNode
type(individual), pointer :: item =>null()
type(IndividualLinkedListNode),pointer :: next =>null()
type(IndividualLinkedListNode),pointer :: previous =>null()

contains
	final:: destroyIndividualLinkedListNodeFinal
end type IndividualLinkedListNode


contains

	subroutine destroyIndividualLinkedListNodeFinal(this)
		type(IndividualLinkedListNode) :: this

		this%item => null()
		this%next =>null()
		this%previous => null()

	end subroutine destroyIndividualLinkedListNodeFinal

	!---------------------------------------------------------------------------
	!> @brief Destructor for linked list
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine destroyLinkedList(this)
		type(IndividualLinkedList),intent(inout) :: this

		type(IndividualLinkedListNode),pointer :: node, tmp


		node => this%first

        print *, "start destruction"
		if (associated(node)) then
			do while (associated(node%next))
				
                if (ASSOCIATED(node,this%last)) then
                    this%last%next => null()
                    this%last%item => null()
                    this%last%previous => null()
                endif
                node%previous=> null()

				tmp => node%next
				node%next=>null()
				node%item=> null()
				deallocate(node)

				node => tmp
			end do

			node%previous=> null()

			! tmp => node%next
			node%next=>null()


			node%previous=> null()
			node%next=>null()
			node%item=> null()
			deallocate(node)
			print *, "stop destruction"
		endif

    
		
        ! if (associated(this%first)) then
        !     if (ASSOCIATED(this%first%item)) then
        !         this%first%item => null()
        !     endif
        ! endif
        ! this%first => null()

        ! print *, "end destruction"
        this%length = 0
        this%first => null()
        this%last => null()
		! deallocate(this%first)
		! deallocate(this%last)

	end subroutine destroyLinkedList


	!     !---------------------------------------------------------------------------
	! !> @brief Destructor for linked list, note that data set is deallocated as well
	! !> @author  David Wilson david.wilson@roslin.ed.ac.uk
	! !> @date    October 26, 2016
	! !---------------------------------------------------------------------------
	! subroutine destroyLinkedListFinal(this)
	!     class(IndividualLinkedList),intent(inout) :: this
	!     type(IndividualLinkedListNode),pointer :: node
	!     type(individual),pointer :: tmp
	!     if (associated(this%first)) then
	!         node => this%first

	!         do while(associated(node))
	!             call this%list_pop(tmp)


	!             node => this%first
	!         enddo
	!     endif
	!     deallocate(this%first)
	!     deallocate(this%last)

	! end subroutine destroyLinkedListFinal

	!---------------------------------------------------------------------------
	!> @brief output for Linked List
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine writeLinkedList(dtv, unit, iotype, v_list, iostat, iomsg)
		class(IndividualLinkedList), intent(in) :: dtv         !< Object to write.
		integer, intent(in) :: unit         !< Internal unit to write to.
		character(*), intent(in) :: iotype  !< LISTDIRECTED or DTxxx
		integer, intent(in) :: v_list(:)    !< parameters from fmt spec.
		integer, intent(out) :: iostat      !< non zero on error, etc.
		character(*), intent(inout) :: iomsg  !< define if iostat non zero.

		type(IndividualLinkedListNode),pointer :: node
		node => dtv%first

		do while (associated(node))

			if (associated(node%next)) then
				write(unit, *, iostat = iostat, iomsg = iomsg) node%item , char(10) !char(10) is line feed for new line
			else
				write(unit, *, iostat = iostat, iomsg = iomsg) node%item
			endif
			! flush(unit)
			node => node%next
		enddo
	end subroutine writeLinkedList

	!---------------------------------------------------------------------------
	!> @brief Add item to end of the list
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine list_add(this,item)
		class(IndividualLinkedList),intent(inout) :: this
		type(individual),intent(in), target :: item !< item to add

		if (.not.associated(this%last)) then
			allocate(this%first)
			this%last => this%first
			this%first%previous => null()
		else
			allocate(this%last%next)
			this%last%next%previous => this%last
			this%last => this%last%next
		endif
		this%last%item => item
		this%length = this%length + 1
	end subroutine list_add



	recursive logical function list_all(this,proc) result(res)
	class(IndividualLinkedList),intent(inout) :: this
	procedure(logicalAbstractFunction) :: proc
	type(IndividualLinkedListNode),pointer :: node
	res = .true.
	node => this%first

	do while (associated(node))
		res =  proc(node%item)
		if (.not.res) return
		node => node%next
	end do

end function list_all

!---------------------------------------------------------------------------
!> @brief returns and then removes the item at the end of the list
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine list_pop(this, item)
	class(IndividualLinkedList),intent(inout) :: this
	type(individual),pointer,intent(out) :: item !< item at the end of the list

	if (associated(this%last)) then
		item => this%last%item
		this%last => this%last%previous
		if (associated(this%last)) then

			! deallocate(this%last%next)
			this%last%next => null()
		else
			if (associated(this%first)) then

				! deallocate(this%first)
				this%first => null()
			endif
		end if
		this%length = this%length - 1
	else
		item => null()
	end if
end subroutine list_pop


!---------------------------------------------------------------------------
!> @brief Gets value at position N in list
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function list_get_nth(this,n) result(res)
	class(individual),pointer :: res !< item returned
	class(IndividualLinkedList),intent(in) :: this
	integer, intent(in) :: n !< position of item to return
	integer :: i
	type(IndividualLinkedListNode),pointer :: node

	if (associated(this%first).and.this%length>=n) then
		node => this%first
		i = 1
		do
			if (i==n) then
				res => node%item
				exit
			else if (i<n.and..not.associated(node%next)) then
				res => null()
				exit
			else
				i = i + 1
				node => node%next
			end if
		end do
	else
		res => null()
	end if
end function list_get_nth


!---------------------------------------------------------------------------
!> @brief Returns .true. if item is in list, false otherwise
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
logical function contains(this, ind)

	use IndividualModule
	class(IndividualLinkedList),intent(in) :: this
	type(individual),target, intent(in) :: ind !< item to check
	type(IndividualLinkedListNode),pointer :: node

	logical :: tmp
	
	if (associated(this%first)) then
		node => this%first
		if (ASSOCIATED(node%item)) then
			do  
				tmp = compareIndividual(node%item, ind)

				if (tmp) then
					contains = .true.
					return
				else if (.not.associated(node%next)) then
					contains = .false.
					return
				else
					node => node%next
				end if
			end do
		endif
	end if
	contains = .false.
end function contains



!---------------------------------------------------------------------------
!> @brief Removes item from list
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine list_remove(this,item)
	use IndividualModule
	class(IndividualLinkedList),intent(inout) :: this
	type(individual),pointer, intent(in) :: item !< item to remove
	type(individual), pointer :: tmpItem
	type(IndividualLinkedListNode),pointer :: node

    print *, "start ll destroy"
	if (associated(this%first)) then
		node => this%first
		if (associated(node%item,item)) then
			this%first => this%first%next
			this%length = this%length - 1
			if (.not. associated(this%first)) then
				deallocate(this%last)
				!print *, "LIST EMPTY:", this%length
			else
				this%first%previous => null()
			endif
		else
			do while (associated(node%next))
				tmpItem => node%next%item
				! print *,"loop"
				if(associated(tmpItem,item)) then
					if(associated(node%next%next)) then
						deallocate(node%next%next%previous)
						node%next%next%previous => node
						node%next=> node%next%next
					else
						this%last=>node !set last element to node if node to remove is last
						deallocate(this%last%next)
					endif
					this%length = this%length - 1
					return
				endif
				node => node%next
			end do
		endif
	else
		print *, "warning -trying to remove from an empty list"
	endif
    print *, "stop ll destroy"
end subroutine list_remove

!---------------------------------------------------------------------------
!> @brief Converts linked list to 1 dimensional array (vector) of individual objects
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function convertToArray(this) result(res)

	use individualModule
	class(IndividualLinkedList) :: this !< linked list
	type(individual),pointer, dimension(:) :: res !< one dimensional array of animal pointers to return
	integer :: counter
	type(IndividualLinkedListNode),pointer :: node


	counter = 1
	allocate(res(this%length))
	node => this%first

	do while (associated(node))
		res(counter) = node%item

		counter = counter+1
		node => node%next

	enddo
end function convertToArray




!---------------------------------------------------------------------------
!> @brief Converts linked list to 1 dimensional array (vector) of individual objects
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function convertToListOfKnownAnimals(this) result(res)

	use individualModule
	class(IndividualLinkedList) :: this !< linked list
	type(IndividualLinkedList) :: res !< one dimensional array of animal pointers to return
	type(IndividualLinkedListNode),pointer :: node


	node => this%first

	do while (associated(node))

		if (.not. node%item%isUnknownDummy) then
			call res%list_add(node%item)
		end if
		node=> node%next

	enddo
end function convertToListOfKnownAnimals


!---------------------------------------------------------------------------
!> @brief Converts linked list to 1 dimensional array (vector) of integer recoded id's
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function convertToArrayIDs(this) result(res)
	class(IndividualLinkedList) :: this !< linked list
	integer, dimension(:), allocatable :: res !< one dimensional array of recoded id's to return
	integer :: counter
	type(IndividualLinkedListNode),pointer :: node


	counter = 1
	allocate(res(this%length))
	node => this%first

	do while (associated(node))
		res(counter) = node%item%id

		counter = counter+1
		node => node%next

	enddo
end function convertToArrayIDs



!---------------------------------------------------------------------------
!> @brief Converts linked list to 1 dimensional array (vector) of integer recoded id's
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function convertToArrayOriginalIDs(this) result(res)
	use ConstantModule
	class(IndividualLinkedList) :: this !< linked list
	character(len=IDLENGTH), dimension(:), allocatable :: res !< one dimensional array of recoded id's to return
	integer :: counter
	type(IndividualLinkedListNode),pointer :: node


	counter = 1
	allocate(res(this%length))
	node => this%first

	do while (associated(node))

		if (associated(node%item)) then
            print *,"LENGH",this%length
			res(counter) = node%item%originalID

			counter = counter+1
			node => node%next
		endif

	enddo
end function convertToArrayOriginalIDs

!---------------------------------------------------------------------------
!> @brief Across a list of individuals, return all the genotypes at a given position
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function getGenotypesAtPosition(this, pos) result(res)
	class(IndividualLinkedList), intent(in) :: this
	integer, intent(in) :: pos !< snp position
	integer(kind=1), dimension(:), allocatable :: res !< result of all the genotyps
	type(IndividualLinkedListNode), pointer :: ind
	integer :: i
	ind => this%first
	allocate(res(this%length))



	do i=1, this%length
		res(i) = ind%item%individualGenotype%getGenotype(pos)

		ind => ind%next
	enddo

end function getGenotypesAtPosition




!---------------------------------------------------------------------------
!> @brief Across a list of individuals, return all the genotypes at a given position
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine removeIndividualsBasedOnThreshold(this, nOffsThresh, genotyped, hd, genotypedOffspring)

	class(IndividualLinkedList), intent(inout) :: this
	integer, intent(in),optional :: nOffsThresh !< threshold of number of offspring. If less than this these offsprings will be removed
	integer, intent(in), optional :: genotyped,hd !, if either of these are present, then, these are effectively true
	integer, intent(in), optional :: genotypedOffspring !, if present, only care about offspring that are genotpyed
	type(IndividualLinkedListNode), pointer :: ind
	integer :: offspringCount,i

	ind => this%first
	do while(associated(ind))

		if (present(nOffsThresh)) then

			if (present(genotypedOffspring)) then
				offspringCount = 0
				do i=1,ind%item%nOffs

					if (ind%item%offsprings(i)%p%genotyped) then
						offspringCount = offspringCount + 1
					endif
				enddo
			else
				offspringCount = ind%item%nOffs
			endif


			if(offspringCount < nOffsThresh) then
				call this%list_remove(ind%item)
			endif
		endif

		if (present(genotyped)) then
			if(.not. ind%item%Genotyped) then
				call this%list_remove(ind%item)
			endif
		endif

		if (present(hd)) then
			if(.not. ind%item%hd) then
				call this%list_remove(ind%item)
			endif
		endif
		ind => ind%next

	enddo
end subroutine removeIndividualsBasedOnThreshold

end Module IndividualLinkedListModule

