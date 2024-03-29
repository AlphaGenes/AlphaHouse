
module LinkedListModule




	type LIST_DATA
	character(len=:), allocatable :: key
	! class(*)                :: value
	integer :: value

	contains

		final :: clearData
	end type LIST_DATA




	type LinkedList
	type(LinkedList), pointer :: next
	type(LIST_DATA), allocatable            :: data

	contains
		! final :: list_destroy
	end type LinkedList

	interface operator ( == )
		module procedure LIST_DATAEquals
		!module procedure equalLists
	end interface operator ( == )


	interface assignment ( = )
		module procedure copyLinkedList
	end interface assignment ( = )

	contains



		!---------------------------------------------------------------------------
		!< @brief subroutine to deep copy linked list
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		subroutine copyLinkedList(list, this)
			type(LinkedList), pointer, intent(inout) :: list
			type(LinkedList), pointer, intent(in) :: this
			type(LinkedList), pointer :: cur, t
			type(LinkedList), pointer :: next

			t => this
			if (associated(t)) then

				if (.not. ASSOCIATED(list)) then
					allocate(list)
				endif
				allocate(list%data)
				list%next => null()
				list%data%key = t%data%key
				list%data%value = t%data%value
				cur =>list
				do while (associated(t))
					allocate( next )
					cur%next => next
					allocate( next%data)
					next%data%key = t%data%key
					next%data%value = t%data%value
					cur => cur%next
					t=> t%next
				end do
			endif

		end subroutine copyLinkedList


				!---------------------------------------------------------------------------
		!< @brief Function to check linked list equality 
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		logical function equalLists(left,right)
			type(LinkedList), pointer,intent(in) :: left
			type(LinkedList), pointer,intent(in) :: right
			type(LinkedList), pointer :: l,r

			l => left
			r => right

			do while (associated(l))

				! means lists are different lengths
				if (.not. associated(r)) then
					equalLists = .false.
					Return
				endif

				if (.not. (l%data == r%data)) then
					equalLists = .false.
					Return
				endif
				l => l%next
				r => r%next

			enddo
			! means lists are different lengths
			if (associated(r)) then
				equalLists = .false.
				Return
			endif

			equalLists = .true.

		end function equalLists


		! function copyLinkedList(this) result(list)
		! 	type(LinkedList), pointer :: list
		! 	type(LinkedList), pointer :: cur
		! 	type(LinkedList), pointer :: next, this

		! 	if (associated(this)) then
		! 		allocate(list)
		! 		! allocate(list%data)
		! 		list%next => null()
		! 		list%data%key = this%data%key
		! 		list%data%value = this%data%value
		! 		cur =>list
		! 		do while (associated(this%next))
		! 			allocate( next )
		! 			cur%next => next
		! 			next%data%key = this%data%key
		! 			next%data%value = this%data%value
		! 			cur => cur%next
		! 			this=> this%next
		! 		end do
		! 	else
		! 		list => null()
		! 	endif

		! end function copyLinkedList


				!---------------------------------------------------------------------------
		!< @brief Helper Function to check list equality
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		logical function LIST_DATAEquals(l, r)
			class(LIST_DATA), intent(in) :: l,r
			if (l%value == r%value) then
				LIST_DATAEquals = .true.
			else
				LIST_DATAEquals = .false.
			endif

		end function LIST_DATAEquals


				!---------------------------------------------------------------------------
		!< @brief deallocates list data object
		!< @author  David Wilson david.wilson@roslin.ed.ac.uk
		!< @date    October 26, 2017
		!---------------------------------------------------------------------------
		subroutine clearData(d)
			type(LIST_DATA):: d

			if (allocated(d%key)) then
				deallocate(d%key)
			endif
		end subroutine clearData
		! list_create --
		!     Create and initialise a list
		! Arguments:
		!     list       Pointer to new linked listp
		!     data       The data for the first element
		! Note:
		!     This version assumes a shallow copy is enough
		!     (that is, there are no pointers within the data
		!     to be stored)
		!     It also assumes the argument list does not already
		!     refer to a list. Use list_destroy first to
		!     destroy up an old list.
		!
		pure subroutine list_create( list, data )
			type(LinkedList), pointer  :: list
			type(LIST_DATA), intent(in) :: data

			allocate( list )
			list%next => null()
			list%data =  data
		end subroutine list_create

		! list_destroy --
		!     Destroy an entire list
		! Arguments:
		!     list       Pointer to the list to be destroyed
		! Note:
		!     This version assumes that there are no
		!     pointers within the data that need deallocation
		!
		subroutine list_destroy( list )
			type(LinkedList),target  :: list
			type(LinkedList), pointer  :: current
			type(LinkedList), pointer  :: tmp

			current => list
			if (associated(current)) then
				do while (associated(current%next))

					tmp => current
					current => current%next
					deallocate(tmp%data)
					deallocate(tmp)

				end do


				deallocate(current%data)
				deallocate(current)

			endif

		end subroutine list_destroy

		! list_count --
		!     Count the number of items in the list
		! Arguments:
		!     list       Pointer to the list
		!
		integer function list_count( list )
			type(LinkedList), pointer  :: list

			type(LinkedList), pointer  :: current

			if ( associated(list) ) then
				list_count = 1
				current => list
				do while ( associated(current%next) )
					current => current%next
					list_count = list_count + 1
				enddo
			else
				list_count = 0
			endif
		end function list_count

		! list_next
		!     Return the next element (if any)
		! Arguments:
		!     elem       Element in the linked list
		! Result:
		!
		function list_next( elem ) result(next)
			type(LinkedList), pointer :: elem
			type(LinkedList), pointer :: next

			next => elem%next

		end function list_next

		! list_insert
		!     Insert a new element
		! Arguments:
		!     elem       Element in the linked list after
		!                which to insert the new element
		!     data       The data for the new element
		!
		subroutine list_insert( elem, data )
			type(LinkedList), pointer  :: elem
			type(LIST_DATA), intent(in) :: data

			type(LinkedList), pointer :: next

			allocate(next)

			next%next => elem%next
			elem%next => next
			next%data =  data
			! deallocate(data%key)
		end subroutine list_insert

		! list_insert_head
		!     Insert a new element before the first element
		! Arguments:
		!     list       Start of the list
		!     data       The data for the new element
		!
		subroutine list_insert_head( list, data )
			type(LinkedList), pointer  :: list
			type(LIST_DATA), intent(in) :: data

			type(LinkedList), pointer :: elem

			allocate(elem)
			elem%data =  data

			elem%next => list
			list      => elem
		end subroutine list_insert_head

		! list_delete_element
		!     Delete an element from the list
		! Arguments:
		!     list       Header of the list
		!     elem       Element in the linked list to be
		!                removed
		!
		subroutine list_delete_element( list, elem )
			type(LinkedList), pointer  :: list
			type(LinkedList), pointer  :: elem

			type(LinkedList), pointer  :: current
			type(LinkedList), pointer  :: prev

			if ( associated(list,elem) ) then
				list => elem%next
				deallocate( elem )
			else
				current => list
				prev    => list
				do while ( associated(current) )
					if ( associated(current,elem) ) then
						prev%next => current%next
						deallocate( current ) ! Is also "elem"
						exit
					endif
					prev    => current
					current => current%next
				enddo
			endif
			!    allocate(next)
			!
			!    next%next => elem%next
			!    elem%next => next
			!    next%data =  data
		end subroutine list_delete_element

		! list_get_data
		!     Get the data stored with a list element
		! Arguments:
		!     elem       Element in the linked list
		!
		function list_get_data(this) result(data)
			type(LinkedList), pointer :: this
			type(LIST_DATA) :: data
			data = this%data
		end function list_get_data

		! list_put_data
		!     Store new data with a list element
		! Arguments:
		!     elem       Element in the linked list
		!     data       The data to be stored
		!
		subroutine list_put_data( elem, data )
			type(LinkedList), pointer  :: elem
			type(LIST_DATA), intent(in) :: data

			elem%data = data
		end subroutine list_put_data





end module LinkedListModule

