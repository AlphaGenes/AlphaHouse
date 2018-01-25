
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     haplotypeLinkedListModule.f90
!
! DESCRIPTION:
!> @brief    Module containing definition of linked List for objects (currently haplotypes.
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

module haplotypeLinkedListModule
	use iso_fortran_env
	use HaplotypeModule

	implicit none

! 	abstract interface

! 	! Abstract function required to allow for listing of types
! 	logical function logicalAbstractFunction(item)
! 		import :: Haplotype
! 		type(Haplotype),intent(in) :: item
! 	end function logicalAbstractFunction
! end interface


type :: HaplotypeLinkedList
type(HaplotypeLinkedListNode),pointer :: first => null()
type(HaplotypeLinkedListNode),pointer :: last => null()
! TODO - maybe have a middle?
integer :: length = 0

contains
	procedure :: list_add
	procedure :: list_pop
	procedure :: list_get_nth
	procedure :: list_remove
	procedure :: contains
	procedure :: writeLinkedList
	final :: destroyLinkedList
	procedure :: convertToArray
	! procedure :: destroyLinkedListFinal
	generic:: write(formatted)=> writeLinkedList

end type HaplotypeLinkedList

type :: HaplotypeLinkedListNode
type(Haplotype), pointer :: item =>null()
type(HaplotypeLinkedListNode),pointer :: next =>null()
type(HaplotypeLinkedListNode),pointer :: previous =>null()



contains
	final:: destroyhaplotypeLinkedListNodeFinal
end type HaplotypeLinkedListNode
interface assignment (=)
	module procedure deepCopyLinkedList
end interface

contains

	subroutine deepCopyLinkedList(this, listIn)

		class(HaplotypeLinkedList), intent(inout)::this
    	type(HaplotypeLinkedList), intent(in) ::listIn
		type(HaplotypeLinkedListNode) ,pointer :: node, old
		integer :: i
		
		
		if (associated(this%first)) then
			call destroyLinkedList(this)
		endif
		
		this%length = listIn%length
		allocate(node)


		if (associated(listIn%First)) then
			node%item => listin%First%item
			node%previous => null()
			this%first=>node
		endif

		if (listin%length > 1) then
		old => listin%first%next
			do i =1,listIn%length-1

				allocate(node%next)
				
				node%next%previous => node
				node%next%item => old%item

				if (i == listIn%length-1) then
					this%last => node%next
				endif
				old => old%next
				node => node%next
			end do 
		else if (listIn%length > 0) then
			this%last =>node
		endif

	end subroutine deepCopyLinkedList

	subroutine destroyhaplotypeLinkedListNodeFinal(this)
		type(HaplotypeLinkedListNode) :: this

		this%item => null()
		this%next =>null()
		this%previous => null()

	end subroutine destroyhaplotypeLinkedListNodeFinal

	!---------------------------------------------------------------------------
	!> @brief Destructor for linked list
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine destroyLinkedList(this)
		type(HaplotypeLinkedList),intent(inout) :: this

		type(HaplotypeLinkedListNode),pointer :: node, tmp


		node => this%first
		if (this%length == 0) return
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
			node => null()
			
		endif

    
		
        ! if (associated(this%first)) then
        !     if (ASSOCIATED(this%first%item)) then
        !         this%first%item => null()
        !     endif
        ! endif
        ! this%first => null()

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
	!     class(HaplotypeLinkedList),intent(inout) :: this
	!     type(HaplotypeLinkedListNode),pointer :: node
	!     type(Haplotype),pointer :: tmp
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
		class(HaplotypeLinkedList), intent(in) :: dtv         !< Object to write.
		integer, intent(in) :: unit         !< Internal unit to write to.
		character(*), intent(in) :: iotype  !< LISTDIRECTED or DTxxx
		integer, intent(in) :: v_list(:)    !< parameters from fmt spec.
		integer, intent(out) :: iostat      !< non zero on error, etc.
		character(*), intent(inout) :: iomsg  !< define if iostat non zero.

		type(HaplotypeLinkedListNode),pointer :: node
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
		class(HaplotypeLinkedList),intent(inout) :: this
		type(Haplotype),intent(in), target :: item !< item to add

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


!---------------------------------------------------------------------------
!> @brief returns and then removes the item at the end of the list
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
subroutine list_pop(this, item)
	class(HaplotypeLinkedList),intent(inout) :: this
	type(Haplotype),pointer,intent(out) :: item !< item at the end of the list

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
	class(Haplotype),pointer :: res !< item returned
	class(HaplotypeLinkedList),intent(in) :: this
	integer, intent(in) :: n !< position of item to return
	integer :: i
	type(HaplotypeLinkedListNode),pointer :: node

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
logical function contains(this, hap)

	use HaplotypeModule
	class(HaplotypeLinkedList),intent(in) :: this
	type(Haplotype),target, intent(in) :: hap !< item to check
	type(HaplotypeLinkedListNode),pointer :: node

	logical :: tmp
	
	if (associated(this%first)) then
		node => this%first
		if (ASSOCIATED(node%item)) then
			do  
				tmp = compareHaplotype(node%item, hap)

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
	use HaplotypeModule
	class(HaplotypeLinkedList),intent(inout) :: this
	type(Haplotype),pointer, intent(in) :: item !< item to remove
	type(Haplotype), pointer :: tmpItem
	type(HaplotypeLinkedListNode),pointer :: node

	if (associated(this%first)) then
		node => this%first
		if (associated(node%item,item)) then
			this%first => this%first%next
			this%length = this%length - 1
			if (.not. associated(this%first)) then
				deallocate(this%last)
			else
				this%first%previous => null()
			endif
		else
			do while (associated(node%next))
				tmpItem => node%next%item
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
end subroutine list_remove

!---------------------------------------------------------------------------
!> @brief Converts linked list to 1 dimensional array (vector) of Haplotype objects
!> @author  David Wilson david.wilson@roslin.ed.ac.uk
!> @date    October 26, 2016
!---------------------------------------------------------------------------
function convertToArray(this) result(res)

	use HaplotypeModule
	class(HaplotypeLinkedList) :: this !< linked list
	type(Haplotype),pointer, dimension(:) :: res !< one dimensional array of animal pointers to return
	integer :: counter
	type(HaplotypeLinkedListNode),pointer :: node


	counter = 1
	allocate(res(this%length))
	node => this%first

	do while (associated(node))
		res(counter) = node%item

		counter = counter+1
		node => node%next

	enddo
end function convertToArray







end Module haplotypeLinkedListModule

