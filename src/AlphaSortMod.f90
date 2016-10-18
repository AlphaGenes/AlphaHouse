module AlphaSortMod

    interface HpSort
        procedure :: HpSortReal
        procedure :: HpSortI
    end interface HpSort

  contains
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
  integer(kind=8),intent(inout)::RA(N)
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