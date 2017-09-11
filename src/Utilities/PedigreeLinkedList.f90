module MyLinkedList
  use IntegerLinkedListModule
  implicit none

  type, extends(IntegerLinkedList):: myLL
    contains
      procedure:: list_add => myListlist_add
      procedure:: pop_first
!      procedure:: destroyLinkedList => destroyLinkedListMyLL
  end type

  contains

    subroutine destroyLinkedListMyLL(this)
      class(myLL), intent(inout):: this
      type(IntegerLinkedListNode), pointer:: node

      if (associated(this%last)) then
!        this%last%next => null()
!        this%last%previous => null()
        this%last=> null()
      end if

      if (associated(this%first)) then
        node => this%first
        do while (associated(node))
          node=> node%next
          node%previous => null()
        end do
      end if

      node=> null()

    end subroutine

    subroutine myListlist_add(this, item)
      class(myLL), intent(inout):: this
      integer, intent(in):: item

      type(IntegerLinkedListNode), pointer:: node, previous

      if (associated(this%first)) then
        node => this%first
        do while (item< node%item) 
          if (associated(node%next)) then
            node=> node%next
          else
            allocate(node%next)
            this%last=> node%next
            node%next%item = item
            node%next%previous => node
            node%next%next =>null()
            this%length = this%length+1
          end if
        end do
        if (item > node%item) then
          if (associated(node%previous)) then
            previous => node%previous
            allocate(node%previous)
            node%previous%item = item
            previous%next => node%previous
            node%previous%previous=> previous
            node%previous%next => node
          else
            allocate(node%previous)
            node%previous%item= item
            node%previous%next=> node
            node%previous%previous => null()
            this%first=>node%previous
          end if
            this%length = this%length+1
        end if
      else
        allocate(this%first)
        this%first%previous => null()
        this%first%next => null()
        this%last => this%first
        this%first%item = item
        this%length = 1
      end if
     node=> null()
!     next=> null()
     previous =>null()
   end subroutine myListlist_add

   function pop_first(self) result(itemOut)
     class(myLL), intent(inout):: self
     integer:: itemOut

     if (associated(self%first)) then
       itemOut = self%first%item
       if (associated(self%first%next)) then
         self%first=> self%first%next
         self%first%previous=>null()
       else
         self%first=> null()
       end if
!       write(*,*) "ASSOCIATED", self%length
       self%length = self%length-1
!       write(*,*) "ASSOCIATED", self%length
    else
!      write(*,*) "NOT ASSOCIATED"
      itemOut = DICT_NULL
    end if

    end function pop_first

end module MyLinkedList
