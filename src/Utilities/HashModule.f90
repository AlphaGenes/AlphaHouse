
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

type HASH_LIST
    type(LinkedList), pointer :: list
end type HASH_LIST

type DictStructure
    
    type(HASH_LIST), pointer, dimension(:),private :: table

    contains

    procedure :: destroy
    procedure :: addKey
    procedure :: deleteKey
    procedure :: getValue
    procedure :: hasKey
    procedure :: getElement
    procedure :: getSize
end type DictStructure

interface DictStructure
    module procedure dict_create
    module procedure dict_create_val
end interface DictStructure
!
! We do not want everything to be public
!
private :: LIST_DATA
private :: HASH_LIST
private :: LinkedList

private :: getElement
private :: hashKey

integer(kind=int64) :: hash_size  = DEFAULTDICTSIZE
integer(kind=int64), parameter, private :: multiplier = 17


contains
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


   !---------------------------------------------------------------------------
  !> @brief Constructor for dictionary
  !> @author  David Wilson david.wilson@roslin.ed.ac.uk
  !> @date    October 26, 2016
  !---------------------------------------------------------------------------
function dict_create(size) result(dict)
    type(DictStructure)  :: dict !< Dictionary Object
    integer(kind=int64),intent(in), optional :: size !< size of underlying array datastructure
    integer(kind=int64) :: i

    if (present(size)) then
        hash_size = size
    endif
    ! allocate( dict )
    allocate( dict%table(hash_size) )

    do i = 1,hash_size
        dict%table(i)%list => null()
    enddo

end function dict_create

  !---------------------------------------------------------------------------
  !> @brief Constructor for dictionary withkey and value parameters 
  !> @author  David Wilson david.wilson@roslin.ed.ac.uk
  !> @date    October 26, 2016
  !---------------------------------------------------------------------------
function dict_create_val(key, value, size ) result(dict)
    type(DictStructure)  :: dict !< Dictionary Object
    character(len=*), intent(in) ::  key !< key for value in dictionary
    integer(kind=int64),intent(in), optional :: size !< size of underlying array datastructure
    integer, intent(in)  :: value !< value to be stored in dictionary
    type(LIST_DATA) :: data
    integer(kind=int64):: i
    integer(kind=int64):: hash

    if (present(size)) then
        hash_size = size
    endif
    allocate( dict%table(hash_size) )

    do i = 1,hash_size
        dict%table(i)%list => null()
    enddo

    data%key   = key
    data%value = value

    hash = hashKey( trim(key ) )
    call list_create( dict%table(hash)%list, data )

end function dict_create_val

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
    class(DictStructure) :: this

    integer(kind=int64) :: i

    do i = 1,size(this%table)
        if ( associated( this%table(i)%list ) ) then
            call list_destroy( this%table(i)%list )
        endif
    enddo
    deallocate(this%table )

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

    type(LIST_DATA)              :: data
    type(LinkedList), pointer   :: elem
    integer(kind=int64) :: hash

    elem => getElement( this, key )

    if ( associated(elem) ) then
        elem%data%value = value
    else
        data%key   = key
        data%value = value
        hash       = hashKey( trim(key) )
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
        hash = hashKey( trim(key) )
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

    hash = hashKey( trim(key) ) !< if key is empty string, this will return 0 and cause segfault error
    elem => this%table(hash)%list
    do while ( associated(elem) )
        if ( elem%data%key .eq. key ) then
            exit
        else
            elem => list_next( elem )
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
  !> @param[out] hashKey. Integer( int32).  
  !> @param[in] sizeOut. String.  
  !---------------------------------------------------------------------------
 function hashKey(key)
    character(len=*), intent(in) :: key

    integer(kind=int64) :: i
    integer(kind=int64) :: hashKey
    hashKey = 0
    do i = 1,len(key)
        hashKey = multiplier * hashKey + ichar(key(i:i))
    enddo
   
    hashKey = 1 + KMOD( hashKey-1, hash_size )
end function hashKey



end module HashModule
