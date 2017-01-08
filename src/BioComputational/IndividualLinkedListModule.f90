
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
    use individualModule, only :individual


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
            procedure :: destroyLinkedList
            generic:: write(formatted) => writeLinkedList

    end type IndividualLinkedList

    type :: IndividualLinkedListNode
            type(individual), pointer :: item
            type(IndividualLinkedListNode),pointer :: next =>null()
            type(IndividualLinkedListNode),pointer :: previous =>null()
    end type IndividualLinkedListNode


contains

    !---------------------------------------------------------------------------
    !> @brief Destructor for linked list
    !> @author  David Wilson david.wilson@roslin.ed.ac.uk
    !> @date    October 26, 2016
    !---------------------------------------------------------------------------
    subroutine destroyLinkedList(this)
        class(IndividualLinkedList),intent(inout) :: this
        type(IndividualLinkedListNode),pointer :: node
        type(individual),pointer :: tmp
        if (associated(this%first)) then
            node => this%first

            do while(associated(node))
                call this%list_pop(tmp)
                node => this%first
            enddo
        endif
        deallocate(this%first)
        deallocate(this%last)
        ! deallocate(tmp)

    end subroutine destroyLinkedList

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
                deallocate(this%last%next)
            else
                deallocate(this%first)
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
    logical function contains(this, in)
        class(IndividualLinkedList),intent(in) :: this
        type(individual),target, intent(in) :: in !< item to check
        type(IndividualLinkedListNode),pointer :: node


        if (associated(this%first)) then
          node => this%first

          do
            if (associated(node%item,in)) then
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
        class(IndividualLinkedList),intent(inout) :: this
        type(individual),pointer, intent(in) :: item !< item to remove
        type(individual), pointer :: tmpItem
        type(IndividualLinkedListNode),pointer :: node
        if (associated(this%first)) then
            node => this%first
            if (associated(node%item,item)) then
                this%first => this%first%next
                this%length = this%length - 1 
                if (.not. associated(this%first)) then
                    deallocate(this%last)
                    print *, "LIST EMPTY:", this%length
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

    end subroutine list_remove



end Module IndividualLinkedListModule