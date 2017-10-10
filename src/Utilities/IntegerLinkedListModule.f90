

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     IntegerLinkedListModule.f90
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

module IntegerLinkedListModule
	use iso_fortran_env
	use ConstantModule
	implicit none

	type :: IntegerLinkedList
	type(IntegerLinkedListNode),pointer :: first => null()
	type(IntegerLinkedListNode),pointer :: last => null()
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
		generic:: write(formatted)=> writeLinkedList

	end type IntegerLinkedList

	type :: IntegerLinkedListNode
	integer :: item
	type(IntegerLinkedListNode),pointer :: next =>null()
	type(IntegerLinkedListNode),pointer :: previous =>null()
	contains

		final :: destroyIntegerLinkedListNode
	end type IntegerLinkedListNode

	interface assignment (=)
		module procedure deepCopyLinkedList
	end interface

	contains



		subroutine deepCopyLinkedList(this, listIn)

			class(IntegerLinkedList), intent(inout)::this
			type(IntegerLinkedList), intent(in) ::listIn
			type(IntegerLinkedListNode) ,pointer :: node, old
			integer :: i


			if (associated(this%first)) then
				call destroyLinkedList(this)
			endif

			this%length = listIn%length
			allocate(node)

			node%item = listin%First%item
			node%previous => null()
			this%first=>node

			if (listin%length > 1) then
				old => listin%first%next
				do i =1,listIn%length-1

					allocate(node%next)

					node%next%previous => node
					node%next%item = old%item

					if (i == listIn%length-1) then
						this%last => node%next
					endif
					old => old%next
					node => node%next
				end do
			else
				this%last =>node
			endif

		end subroutine deepCopyLinkedList

		!---------------------------------------------------------------------------
		!> @brief Destructor for linked list
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine destroyLinkedList(this)
			type(IntegerLinkedList),intent(inout) :: this

			type(IntegerLinkedListNode),pointer :: node, tmp


			node => this%first
			if (this%length == 0) return
			if (associated(node)) then
				do while (associated(node%next))

					if (ASSOCIATED(node,this%last)) then
						this%last%next => null()
						this%last%previous => null()
					endif
					node%previous=> null()

					tmp => node%next
					node%next=>null()
					deallocate(node)

					node => tmp
				end do

				node%previous=> null()

				! tmp => node%next
				node%next=>null()


				node%previous=> null()
				node%next=>null()
				deallocate(node)
				node => null()

			endif
			this%length = 0
			this%first => null()
			this%last => null()


		end subroutine destroyLinkedList



		subroutine destroyIntegerLinkedListNode(node)
			type(IntegerLinkedListNode) :: node

			node%next => null()
			node%previous => null()

		end subroutine destroyIntegerLinkedListNode
		!---------------------------------------------------------------------------
		!> @brief output for Linked List
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine writeLinkedList(dtv, unit, iotype, v_list, iostat, iomsg)
			class(IntegerLinkedList), intent(in) :: dtv         !< Object to write.
			integer, intent(in) :: unit         !< Internal unit to write to.
			character(*), intent(in) :: iotype  !< LISTDIRECTED or DTxxx
			integer, intent(in) :: v_list(:)    !< parameters from fmt spec.
			integer, intent(out) :: iostat      !< non zero on error, etc.
			character(*), intent(inout) :: iomsg  !< define if iostat non zero.

			type(IntegerLinkedListNode),pointer :: node
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
			class(IntegerLinkedList),intent(inout) :: this
			integer,intent(in) :: item !< item to add

			if (.not.associated(this%last)) then
				allocate(this%first)
				this%last => this%first
				this%first%previous => null()
			else
				allocate(this%last%next)
				this%last%next%previous => this%last
				this%last => this%last%next
			endif
			this%last%item = item
			this%length = this%length + 1
		end subroutine list_add



		!---------------------------------------------------------------------------
		!> @brief returns and then removes the item at the end of the list
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine list_pop(this, item)
			class(IntegerLinkedList),intent(inout) :: this
			integer,intent(out) :: item !< item at the end of the list
			if (associated(this%last)) then
				item = this%last%item
				this%last => this%last%previous
				if (associated(this%last)) then
					deallocate(this%last%next)
				else
					deallocate(this%first)
				end if
				this%length = this%length - 1
			else
				item = DICT_NULL
			end if
		end subroutine list_pop


		!---------------------------------------------------------------------------
		!> @brief Gets value at position N in list
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		function list_get_nth(this,n) result(res)
			integer :: res !< item returned
			class(IntegerLinkedList),intent(in) :: this
			integer, intent(in) :: n !< position of item to return
			integer :: i
			type(IntegerLinkedListNode),pointer :: node

			if (associated(this%first).and.this%length>=n) then
				node => this%first
				i = 1
				do
					if (i==n) then
						res = node%item
						exit
					else if (i<n.and..not.associated(node%next)) then
						res = DICT_NULL
						exit
					else
						i = i + 1
						node => node%next
					end if
				end do
			else
				res = DICT_NULL
			end if
		end function list_get_nth


		!---------------------------------------------------------------------------
		!> @brief Returns .true. if item is in list, false otherwise
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		logical function contains(this, in)
			class(IntegerLinkedList),intent(in) :: this
			integer,target, intent(in) :: in !< item to check
			type(IntegerLinkedListNode),pointer :: node


			if (associated(this%first)) then
				node => this%first

				do
					if (node%item == in) then
						contains = .true.
						return
					else if (.not.associated(node%next)) then
						contains = .false.
						return
					else
						node => node%next
					end if
				end do
			else
				contains = .false.
			end if
		end function contains



		!---------------------------------------------------------------------------
		!> @brief Removes item from list
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		subroutine list_remove(this,item)
			class(IntegerLinkedList),intent(inout) :: this
			integer, intent(in) :: item !< item to remove
			integer, pointer :: tmpItem
			type(IntegerLinkedListNode),pointer :: node
			if (associated(this%first)) then
				node => this%first
				if (node%item == item) then
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
						! print *,"loop"
						if(tmpItem == item) then
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
		!> @brief Converts linked list to 1 dimensional array (vector) of individual objects
		!> @author  David Wilson david.wilson@roslin.ed.ac.uk
		!> @date    October 26, 2016
		!---------------------------------------------------------------------------
		function convertToArray(this) result(res)
			class(IntegerLinkedList) :: this !< linked list
			integer, dimension(:), allocatable :: res !< one dimensional array of integers to return
			integer :: counter
			type(IntegerLinkedListNode),pointer :: node


			counter = 1
			allocate(res(this%length))
			node => this%first

			do while (associated(node))
				res(counter) = node%item

				counter = counter+1
				node => node%next

			enddo
		end function convertToArray


end Module IntegerLinkedListModule

