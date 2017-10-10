
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     HashModule.f90
!
! DESCRIPTION:
!> @brief    Module containing definition for Dictionary
!
!> @details  Fully doubly linked list with useful procedures for operations on the dictionary
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

module HashModule
	use iso_fortran_env
	use LinkedListModule
	use constantModule
	implicit none


	type HASH_LIST
	    type(LinkedList), pointer :: list => null()
    contains 
        final :: destroyHASH_LIST
end type HASH_LIST

type DictStructure

type(HASH_LIST), pointer, dimension(:) :: table => null()
integer(kind=int64) :: hash_size  = 0

contains

	final :: destroy
	procedure :: addKey
	procedure :: deleteKey
	procedure :: getValue
	procedure :: hasKey
	procedure :: getElement
	procedure :: getSize
	procedure :: getAllKeys
	procedure :: hashKey
	procedure :: dict_create
	procedure :: dict_create_val

	generic :: DictStructure => dict_create, dict_create_val
end type DictStructure

interface DictStructure
	module procedure dict_create
	module procedure dict_create_val
end interface DictStructure

interface assignment ( = )
       module procedure deepCopyHash
end interface assignment ( = )


! We do not want everything to be public
!
private :: LIST_DATA
private :: HASH_LIST
private :: LinkedList

private :: getElement
private :: hashKey


! TODO need to write deep copy 


contains

    subroutine destroyHASH_LIST(hashList)

        type(HASH_LIST) :: hashList
        hashList%list => null()
        deallocate(hashList%list)

    end subroutine
	!
	! Routines and functions specific to dictionaries
	!

	!---------------------------------------------------------------------------
	!> @brief Returns size of underlying table of dictionary
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	integer function getSize(this)
		class(DictStructure) :: this
		getSize = size(this%table)

	end function getSize



	subroutine deepCopyHash(hash,  old)
	type(DictStructure) ,intent(inout):: hash
	type(DictStructure), intent(in) :: old
	integer :: i
	
	hash%hash_size = old%hash_size
	allocate(hash%table(hash%hash_size))
	do i=1, hash%hash_size

		hash%table(i) = old%table(i)

	enddo
		
	end subroutine deepCopyHash


	function getAllKeys(this) result(res)
		use CharacterLinkedListModule

		class(DictStructure) :: this
		character(len=IDLENGTH), dimension(:), allocatable :: res !< output array
		type(CharacterLinkedList) :: list
		type(LinkedList), pointer :: tmp
		integer :: i


		do i=1,this%hash_size

			tmp => this%table(i)%list

			do while (associated(tmp))
                
                if (allocated(tmp%data%key)) then
				    call list%list_add(tmp%data%key)
                endif
				tmp => tmp%next

			end do


		enddo

		res = list%convertToArray()


	end function getAllKeys

	!---------------------------------------------------------------------------
	!> @brief Constructor for dictionary
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine dict_create(dict,size)
		class(DictStructure)  :: dict !< Dictionary Object
		integer(kind=int64),intent(in), optional :: size !< size of underlying array datastructure
		integer(kind=int64) :: i

		if (present(size)) then
			dict%hash_size = size
		else 
			dict%hash_size = DEFAULTDICTSIZE
		endif
		! allocate( dict )
		allocate( dict%table(dict%hash_size) )

		do i = 1,dict%hash_size
			dict%table(i)%list => null()
		enddo

	end subroutine dict_create

	!---------------------------------------------------------------------------
	!> @brief Constructor for dictionary withkey and value parameters
	!> @author  David Wilson david.wilson@roslin.ed.ac.uk 
	!> @date    October 26, 2016
	!---------------------------------------------------------------------------
	subroutine dict_create_val(dict,key, value, size )
		class(DictStructure)  :: dict !< Dictionary Object
		character(len=*), intent(in) ::  key !< key for value in dictionary
		integer(kind=int64),intent(in), optional :: size !< size of underlying array datastructure
		integer, intent(in)  :: value !< value to be stored in dictionary
		type(LIST_DATA), allocatable :: data
		integer(kind=int64):: i
		integer(kind=int64):: hash

		if (present(size)) then
			dict%hash_size = size
		endif
		allocate( dict%table(dict%hash_size) )

		do i = 1,dict%hash_size
			dict%table(i)%list => null()
		enddo


		allocate(data)
		data%key   = key
		data%value = value

		hash = dict%hashKey( trim(key ) )
		call list_create( dict%table(hash)%list, data )

	end subroutine dict_create_val

	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief     deinitialises a hashtable and frees all memory it was required to store.
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	!-------------------------------------------
	subroutine destroy(this)
		type(DictStructure) :: this

		integer(kind=int64) :: i
		
		print *,"start destroy"


		if (this%hash_size == 0) return
		if (associated(this%table)) then
			do i = 1,this%hash_size

				call list_destroy(this%table(i)%list )
                this%table(i)%list => null()
				! print *,"here"
				! if ( associated( this%table(i)%list ) ) then
				! 	! call list_destroy( this%table(i)%list )
				! 	this%table(i)%list => null()
				! endif
			enddo
			deallocate(this%table)
			nullify(this%table)
			this%hash_size = 0
		endif
		print *,"STOP"

	end subroutine destroy


	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief     Adds a new key, value pair to hashtable
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	! PARAMETERS:
	!> @param[in] Value. class(*)
	!> @param[in] key. String.
	!---------------------------------------------------------------------------
	subroutine addKey( this, key, value )
		class(DictStructure)  :: this
		character(len=*), intent(in) :: key
		integer, intent(in)  :: value

		type(LIST_DATA),allocatable              :: data
		type(LinkedList), pointer   :: elem
		integer(kind=int64) :: hash


		elem => getElement( this, key )

		if ( associated(elem) ) then
			elem%data%value = value
		else
			allocate(data)
			data%key   = key
			data%value = value
			hash       = this%hashKey( trim(key) )
			if ( associated( this%table(hash)%list ) ) then
				call list_insert( this%table(hash)%list, data )
			else
				call list_create( this%table(hash)%list, data )
			endif
		endif

	end subroutine addKey

	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief     Deletes key from hashtable
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	! PARAMETERS:
	!> @param[in] Value. class(*)
	!> @param[in] key. String.
	!---------------------------------------------------------------------------
	subroutine deleteKey( this, key )
		class(DictStructure) :: this
		character(len=*), intent(in) :: key

		type(LinkedList), pointer   :: elem
		integer(kind=int64) :: hash

		elem => this%getElement(key )

		if ( associated(elem) ) then
			hash = this%hashKey( trim(key) )
			call list_delete_element( this%table(hash)%list, elem )
		endif
	end subroutine deleteKey

	! getValue
	!     Get the value belonging to a key
	! Arguments:
	!     dict       Pointer to the dictionary
	!     key        Key for which the values are sought
	!
	function getValue( this, key ) result(value)
		class(DictStructure)   :: this
		character(len=*), intent(in) :: key
		integer :: value

		type(LinkedList), pointer   :: elem

		elem => this%getElement( key )

		if ( associated(elem) ) then
			value = elem%data%value
		else
			value = DICT_NULL
		endif
	end function getValue

	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief      checks if hashtable contains key
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	! PARAMETERS:
	!> @param[out] boolean relating to if structure contains key.
	!> @param[in] key. String.
	!---------------------------------------------------------------------------
	logical function hasKey( this, key )
		class(DictStructure) :: this
		character(len=*), intent(in) :: key

		type(LinkedList), pointer :: elem

		elem => this%getElement(key )

		hasKey = associated(elem)
	end function hasKey


	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief      gets linked list element a at given key
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	! PARAMETERS:
	!> @param[out] Element. Linked List Node.
	!> @param[in] key. String.
	!---------------------------------------------------------------------------
	function getElement( this, key ) result(elem)
		class(DictStructure) :: this
		character(len=*), intent(in) :: key
		type(LinkedList), pointer :: elem
		integer(kind=int64) :: hash


		hash = this%hashKey(trim(key)) !< if key is empty string, this will return 0 and cause segfault error
        elem => this%table(hash)%list
        do while ( associated(elem) )
            if ( associated(elem)) then
                ! if (allocated(elem%data)) then
                    if (allocated(elem%data%key)) then
                        if ( elem%data%key .eq. key ) then
                            exit
                        else
                            elem => list_next( elem )
                        endif
                    endif
                ! endif
            endif
		enddo

	end function getElement



	!---------------------------------------------------------------------------
	! DESCRIPTION:
	!> @brief      returns a hash of given string
	!
	!> @author     David Wilson, david.wilson@roslin.ed.ac.uk
	!
	!> @date       October 25, 2016
	!
	! PARAMETERS:
	!> @param[in] sizeOut. String.
	!---------------------------------------------------------------------------
	function hashKey(dict, key)
		character(len=*), intent(in) :: key
		class(DictStructure) :: dict

		integer(kind=int64) :: i
		integer(kind=int64) :: hashKey !<hashkey out
		hashKey = 0
		do i = 1,len(key)
			hashKey = KMOD(DICT_MULTIPLIER * hashKey + ichar(key(i:i)), dict%hash_size)
		enddo

		hashKey = 1 + hashKey
	end function hashKey



end module HashModule

