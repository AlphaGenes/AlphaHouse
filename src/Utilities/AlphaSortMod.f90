
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaSortMod.f90
!
! DESCRIPTION:
!> @brief    Module cotaining simple sort procedures
!
!> @details  currently only contains integer and real heap sort procedures 
!
!> @author   David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 DWilson - Initial Version
!
!-------------------------------------------------------------------------------

module AlphaSortMod
    
    use iso_fortran_env
    interface HpSort
        procedure :: HpSortReal
        procedure :: HpSortI
        procedure :: HpSortDS
    end interface HpSort

    interface mergeSortedLists
      module procedure mergeTwoSortedListsInteger, mergeTwoSortedListsReal, mergeTwoSortedListsReal64, mergeTwoSortedListsInteger64
    end interface

  contains
    !> @brief A function that takes in two (sorted) lists and produces a union of the unique elements
    !> @details Takes in two list of integers, and creates a sorted output list of the unique elements.
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    function mergeTwoSortedListsReal64(listOne, listTwo ) result (listOut)
      real(real64), dimension(:),intent(in):: listOne
      real(real64), dimension(:), intent(in):: listTwo

      real(real64), dimension(:), allocatable:: listOut, tempList

      integer:: firstListSize, secondListSize, firstIndex, secondIndex, outputIndex

      firstListSize = size(listOne)
      secondListSize = size(listTwo)

      allocate(tempList(firstListSize+secondListSize))

      firstIndex = 1
      secondIndex = 1
      outputIndex = 1
      do while(firstIndex <= firstListSize .and. secondIndex<= secondListSize)
        if (listOne(firstIndex) < listTwo(secondIndex)) then
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
        else if (listOne(firstIndex) > listTwo(secondIndex)) then
          tempList(outputIndex) = listTwo(secondIndex)
          secondIndex = secondIndex+1
        else 
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
          secondIndex = secondIndex+1
        end if
        outputIndex = outputIndex+1
      end do

      if (firstIndex == firstListSize+1) then
        listRange = secondListSize-secondIndex
        tempList(outputIndex:outputIndex+listRange) = listTwo(secondIndex:)
      else if (secondIndex == secondListSize+1) then
        listRange = firstListSize - firstIndex
        tempList(outputIndex+1: outputIndex+listRange) = listOne(firstIndex:)
      end if
      outputIndex = outputIndex+listRange

      listOut = tempList(:outputIndex)
    end function mergeTwoSortedListsReal64

    !> @brief A function that takes in two (sorted) lists and produces a union of the unique elements
    !> @details Takes in two list of integers, and creates a sorted output list of the unique elements.
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    function mergeTwoSortedListsReal(listOne, listTwo ) result (listOut)
      real, dimension(:),intent(in):: listOne
      real, dimension(:), intent(in):: listTwo

      real, dimension(:), allocatable:: listOut, tempList

      integer:: firstListSize, secondListSize, firstIndex, secondIndex, outputIndex

      firstListSize = size(listOne)
      secondListSize = size(listTwo)

      allocate(tempList(firstListSize+secondListSize))

      firstIndex = 1
      secondIndex = 1
      outputIndex = 1
      do while(firstIndex <= firstListSize .and. secondIndex<= secondListSize)
        if (listOne(firstIndex) < listTwo(secondIndex)) then
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
        else if (listOne(firstIndex) > listTwo(secondIndex)) then
          tempList(outputIndex) = listTwo(secondIndex)
          secondIndex = secondIndex+1
        else 
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
          secondIndex = secondIndex+1
        end if
        outputIndex = outputIndex+1
      end do

      if (firstIndex == firstListSize+1) then
        listRange = secondListSize-secondIndex
        tempList(outputIndex:outputIndex+listRange) = listTwo(secondIndex:)
      else if (secondIndex == secondListSize+1) then
        listRange = firstListSize - firstIndex
        tempList(outputIndex+1: outputIndex+listRange) = listOne(firstIndex:)
      end if
      outputIndex = outputIndex+listRange

      listOut = tempList(:outputIndex)
    end function mergeTwoSortedListsReal

    !> @brief A function that takes in two (sorted) lists and produces a union of the unique elements
    !> @details Takes in two list of integers, and creates a sorted output list of the unique elements.
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    function mergeTwoSortedListsInteger64(listOne, listTwo ) result (listOut)
      integer(int64), dimension(:),intent(in):: listOne
      integer(int64), dimension(:), intent(in):: listTwo

      integer(int64), dimension(:), allocatable:: listOut, tempList

      integer:: firstListSize, secondListSize, firstIndex, secondIndex, outputIndex

      firstListSize = size(listOne)
      secondListSize = size(listTwo)

      allocate(tempList(firstListSize+secondListSize))

      firstIndex = 1
      secondIndex = 1
      outputIndex = 1
      do while(firstIndex <= firstListSize .and. secondIndex<= secondListSize)
        if (listOne(firstIndex) < listTwo(secondIndex)) then
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
        else if (listOne(firstIndex) > listTwo(secondIndex)) then
          tempList(outputIndex) = listTwo(secondIndex)
          secondIndex = secondIndex+1
        else 
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
          secondIndex = secondIndex+1
        end if
        outputIndex = outputIndex+1
      end do

      if (firstIndex == firstListSize+1) then
        listRange = secondListSize-secondIndex
        tempList(outputIndex:outputIndex+listRange) = listTwo(secondIndex:)
      else if (secondIndex == secondListSize+1) then
        listRange = firstListSize - firstIndex
        tempList(outputIndex+1: outputIndex+listRange) = listOne(firstIndex:)
      end if
      outputIndex = outputIndex+listRange

      listOut = tempList(:outputIndex)
    end function mergeTwoSortedListsInteger64

    !> @brief A function that takes in two (sorted) lists and produces a union of the unique elements
    !> @details Takes in two list of integers, and creates a sorted output list of the unique elements.
    !> @author Diarmaid de Búrca, diarmaid.deburca@ed.ac.uk
    function mergeTwoSortedListsInteger(listOne, listTwo ) result (listOut)
      integer, dimension(:),intent(in):: listOne
      integer, dimension(:), intent(in):: listTwo

      integer, dimension(:), allocatable:: listOut, tempList

      integer:: firstListSize, secondListSize, firstIndex, secondIndex, outputIndex

      firstListSize = size(listOne)
      secondListSize = size(listTwo)

      allocate(tempList(firstListSize+secondListSize))

      firstIndex = 1
      secondIndex = 1
      outputIndex = 1
      do while(firstIndex <= firstListSize .and. secondIndex<= secondListSize)
        if (listOne(firstIndex) < listTwo(secondIndex)) then
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
        else if (listOne(firstIndex) > listTwo(secondIndex)) then
          tempList(outputIndex) = listTwo(secondIndex)
          secondIndex = secondIndex+1
        else 
          tempList(outputIndex) = listOne(firstIndex)
          firstIndex = firstIndex+1
          secondIndex = secondIndex+1
        end if
        outputIndex = outputIndex+1
      end do

      if (firstIndex == firstListSize+1) then
        listRange = secondListSize-secondIndex
        tempList(outputIndex:outputIndex+listRange) = listTwo(secondIndex:)
      else if (secondIndex == secondListSize+1) then
        listRange = firstListSize - firstIndex
        tempList(outputIndex+1: outputIndex+listRange) = listOne(firstIndex:)
      end if
      outputIndex = outputIndex+listRange

      listOut = tempList(:outputIndex)
    end function mergeTwoSortedListsInteger
   !  (C) Copr. 1986-92 Numerical Recipes Software 6

   !*****************************************************
   !*  Sorts an array RA of length N in ascending order *
   !*                by the Heapsort method             *
   !* ------------------------------------------------- *
   !* INPUTS:                                           *
   !*      N   size of table RA                   *
   !*          RA    table to be sorted                 *
   !* OUTPUT:                                           *
   !*      RA    table sorted in ascending order    *
   !*                                                   *
   !* NOTE: The Heapsort method is a N Log2 N routine,  *
   !*       and can be used for very large arrays.      *
   !*****************************************************
   subroutine HpSortReal(N,RA)
       !MODifIED FOR REALS BY JOHN HICKEY
       use iso_fortran_env
       implicit none
       integer(kind=int64) :: N
       integer(kind=int32)::IR,J,L,I
       double precision :: RRA
       double precision,intent(inout) ::RA(N)
       if (n.lt.2) return

       L=INT((dble(N)/2)+1)
       IR=N
   !The index L will be decremented from its initial value during the
   !"hiring" (heap creation) phase. Once it reaches 1, the index IR
   !will be decremented from its initial value down to 1 during the
   !"retirement-and-promotion" (heap selection) phase.
10 continue
   if(L > 1)then
       L=L-1
       RRA=RA(L)
   else
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       if(IR.eq.1)then
           RA(1)=RRA
           return
       end if
   end if
   I=L
   J=L+L
20 if(J.le.IR)then
       if(J < IR)then
           if(RA(J) < RA(J+1))  J=J+1
       end if
       if(RRA < RA(J))then
           RA(I)=RA(J)
           I=J; J=J+J
       else
           J=IR+1
       end if
       goto 20
   end if
   RA(I)=RRA
   goto 10
   end subroutine HpSortReal

   subroutine HpSortDS(N,RA)
     !MODifIED FOR integer of kind DOubleSize BY Diarmaid de Burca, July 2016
    use iso_fortran_env
    implicit none
    integer(kind=int64),intent(in) :: N
    integer(kind=int32)::IR,J,L,I
  integer :: RRA
  integer(kind=int64),intent(inout)::RA(N)
    if (n.lt.2) return

    L=INT((dble(N)/2)+1)
    IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
10 continue
   if(L > 1)then
       L=L-1
       RRA=RA(L)
   else
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       if(IR.eq.1)then
           RA(1)=RRA
           return
       end if
   end if
   I=L
   J=L+L
20 if(J.le.IR)then
       if(J < IR)then
           if(RA(J) < RA(J+1))  J=J+1
       end if
       if(RRA < RA(J))then
           RA(I)=RA(J)
           I=J; J=J+J
       else
           J=IR+1
       end if
       goto 20
   end if
   RA(I)=RRA
   goto 10
   end subroutine 

subroutine HpSortI(N,RA)

    use iso_fortran_env
    !MODifIED FOR integerS BY JOHN HICKEY
    implicit none
    integer(kind=int32),intent(in) :: N
    integer (kind=int32) :: IR,J,L,I
  integer :: RRA
  integer(kind=int32),intent(inout) :: RA(N)
    if (n.lt.2) return

    L=INT((dble(N)/2)+1)
    IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
10 continue
   if(L > 1)then
       L=L-1
       RRA=RA(L)
   else
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       if(IR.eq.1)then
           RA(1)=RRA
           return
       end if
   end if
   I=L
   J=L+L
20 if(J.le.IR)then
       if(J < IR)then
           if(RA(J) < RA(J+1))  J=J+1
       end if
       if(RRA < RA(J))then
           RA(I)=RA(J)
           I=J; J=J+J
       else
           J=IR+1
       end if
       goto 20
   end if
   RA(I)=RRA
   goto 10
   end subroutine HpSortI

end module AlphaSortMod
