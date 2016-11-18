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
            procedure :: list_add_at_n
            procedure :: destroyLinkedList

    end type IndividualLinkedList

    type :: IndividualLinkedListNode
            type(individual), pointer :: item
            type(IndividualLinkedListNode),pointer :: next =>null()
    end type IndividualLinkedListNode


contains

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

    subroutine list_add(this,item)
        class(IndividualLinkedList),intent(inout) :: this
        type(individual),intent(in), target :: item

        if (.not.associated(this%last)) then
        allocate(this%first)
        this%last => this%first
        else
        allocate(this%last%next)
        this%last => this%last%next
        endif
        this%last%item => item
        this%length = this%length + 1
    end subroutine list_add


    ! subroutine will add element item at position n unless n does not exist, at which point will add to the end of the list
    subroutine list_add_at_n(this, item, n)
        class(IndividualLinkedList), intent(inout) :: this
        type(individual), intent(in),target :: item
        integer, intent(in) :: n
        type(IndividualLinkedListNode),pointer :: node, tmp
        integer :: i
        allocate(tmp)
        tmp%item=>item   
        ! add at the beginning
        if (n== 1) then
            tmp%next=>this%first
            this%first=>tmp
                    
        else
            i = 1
            node => this%first
            do
                if (i+1 == this%length) then
                    node%next=>tmp
                    exit
                else if (i+1 == n ) then
                    tmp%item=>item
                    tmp%next=>node%next
                    node%next=>tmp
                    exit
                endif
                i = i+1
                node => node%next
            end do
        endif
        this%length = this%length+1

    end subroutine list_add_at_n


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


    subroutine list_pop(this, item)
        class(IndividualLinkedList),intent(inout) :: this
        type(individual),pointer,intent(out) :: item
        type(IndividualLinkedListNode),pointer :: node, previous

        if (associated(this%last)) then
          item => this%last%item         
          previous => null()
          node => this%first
          do while (associated(node%next))
            previous => node
            node => node%next
          end do
          if (associated(previous)) then
            deallocate(previous%next)
          else
           deallocate(this%first)
          end if
          this%last => previous
          this%length = this%length - 1
        else
          item => null()
        end if
    end subroutine list_pop

    function list_get_nth(this,n) result(res)
        class(individual),pointer :: res
        class(IndividualLinkedList),intent(in) :: this
        integer, intent(in) :: n
        integer :: i
        type(IndividualLinkedListNode),pointer :: node

        if (associated(this%first).and.this%length>=n) then
          node => this%first
          i = 1
          ! TODO work backwards if closer to the end
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

    subroutine list_remove(this,item)
        class(IndividualLinkedList),intent(inout) :: this
        type(individual),pointer, intent(in) :: item
        type(individual), pointer :: tmpItem
        type(IndividualLinkedListNode),pointer :: node
        if (associated(this%first)) then
            node => this%first
            if (associated(node%item,item)) then
                this%first => null()
                this%first => node%next
                this%length = this%length - 1 
            else
                do while (associated(node%next))
                    tmpItem => node%next%item
                    ! print *,"loop"
                    if(associated(tmpItem,item)) then
                        if(associated(node%next%next)) then
                            node%next=> node%next%next
                        else
                            this%last=>node !set last element to node if node to remove is last 
                        endif
                        this%length = this%length - 1 
                        return
                    endif
                    node => node%next
                end do
            endif
        endif
    end subroutine list_remove



end Module IndividualLinkedListModule