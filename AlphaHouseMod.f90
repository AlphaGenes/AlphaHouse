
!###############################################################################

module AlphaHouseMod

  use ISO_Fortran_Env

  implicit none

  private
  ! Methods
  public :: CountLines,Int2Char,Real2Char,RandomOrder,ToLower,FindLoc,SetSeed

  ! List of characters for case conversion in ToLower
  CHARACTER(*),PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*),PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface Real2Char
    module procedure RealS2Char,RealD2Char
  end interface

  interface FindLoc
    module procedure FindLocC, FindLocI, FindLocS, FindLocD
  end interface

  contains

    !###########################################################################

    function CountLines(FileName) result(nLines)

      implicit none

      ! Arguments
      character(len=*),intent(in) :: FileName
      integer(int32)              :: nLines

      ! Other
      integer(int32) :: f,Unit

      character(len=300) :: DumC

      nLines=0
      open(newunit=Unit,file=trim(FileName),status="old")
      do
        read(Unit,*,iostat=f) DumC
        nLines=nLines+1
        if (f /= 0) then
          nLines=nLines-1
          exit
        end if
      end do
      close(Unit)
      return
    end function

    !###########################################################################

    function Char2Int(c) result(Res)
      ! From http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90

      implicit none

      character(*), intent(in) :: c
      integer(int32)           :: Res

      read(c, *) Res
      return
    end function

    !###########################################################################

    function Int2Char(i,fmt) result(Res)
      ! From http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran

      implicit none

      integer(int32),intent(in)        :: i
      character(*),intent(in),optional :: fmt

      character(:),allocatable :: Res
      character(range(i)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) i
      else
        write(Tmp,"(i0)") i
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    function RealS2Char(r,fmt) result(Res)
      ! From http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran

      implicit none

      real(real32),intent(in)          :: r
      character(*),intent(in),optional :: fmt

      character(:),allocatable :: Res
      character(range(r)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) r
      else
        write(Tmp,"(f)") r
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    function RealD2Char(r,fmt) result(Res)
      ! From http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran

      implicit none

      real(real64),intent(in)          :: r
      character(*),intent(in),optional :: fmt

      character(:),allocatable :: Res
      character(range(r)+2) :: Tmp

      if (present(fmt)) then
        write(Tmp,fmt) r
      else
        write(Tmp,"(f)") r
      end if
      Res=trim(Tmp)
      return
    end function

    !###########################################################################

    function RandomOrder(n) result(Order)
      ! Generate a random ordering of the integers 1 ... n.

      implicit none

      ! Arguments
      integer(int32),intent(in)  :: n
      integer(int32)             :: Order(n)

      ! Other
      integer(int32) :: i,j,k

      real(real64) :: wk(n)

      do i=1,n
        Order(i)=i
      end do

      ! Starting at the end, swap the current last indicator with one
      ! randomly chosen from those preceeding it.
      call random_number(wk)
      do i=n,2,-1
        j=1 + i * wk(i)
        if (j < i) then
          k=Order(i)
          Order(i)=Order(j)
          Order(j)=k
        end if
      end do

      return
    end function

    !###########################################################################

    function ToLower(StringIn) result(StringOut)
      ! From https://groups.google.com/forum/#!topic/comp.lang.fortran/CKx1L2Ahkxg
      character(len=*),intent(in) :: StringIn
      character(len(StringIn))    :: StringOut
      integer :: i,n

      ! Copy input string
      StringOut=StringIn

      ! Convert case character by character
      do i=1,len(StringOut)
        n=index(UPPER_CASE,StringOut(i:i))
        if (n /= 0) StringOut(i:i)=LOWER_CASE(n:n)
      end do
      return
    end function

    !###########################################################################

    function FindLocI(Val,Vec) result(i)
      implicit none
      integer(int32) :: Val
      integer(int32) :: Vec(:)

      integer(int32) :: i,j
      i=0
      do j=1,size(Vec)
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    function FindLocC(Val,Vec) result(i)
      implicit none
      character(len=*) :: Val
      character(len=*) :: Vec(:)

      integer(int32) :: i,j
      i=0
      do j=1,size(Vec)
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    function FindLocS(Val,Vec) result(i)
      implicit none
      real(real32) :: Val
      real(real32) :: Vec(:)

      integer(int32) :: i,j
      i=0
      do j=1,size(Vec)
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    function FindLocD(Val,Vec) result(i)
      implicit none
      real(real64) :: Val
      real(real64) :: Vec(:)

      integer(int32) :: i,j
      i=0
      do j=1,size(Vec)
        if (Val == Vec(j)) then
          i=j
          exit
        end if
      end do
      return
    end function

    !###########################################################################

    subroutine SetSeed(Seed,SeedFile,Out)

      implicit none

      ! Arguments
      integer(int32),intent(in),optional  :: Seed     ! A number to initialize RNG
      character(len=*),optional           :: SeedFile ! File to save the seed in
      integer(int32),intent(out),optional :: Out      ! Make the seed value available outside

      ! Other
      integer(int32) :: Size,Unit
      integer(int32),allocatable :: SeedList(:)

      ! Get the size of seed array
      call random_seed(size=Size)
      allocate(SeedList(Size))

      ! Set seed
      if (present(Seed)) then ! using the given value
        SeedList(1)=Seed
        SeedList(2:Size)=1
        call random_seed(put=SeedList)
      else                    ! using system/compiler value
        call random_seed
        call random_seed(get=SeedList)
      end if

      ! Save to a file
      if (present(SeedFile)) then
        open(newunit=Unit,file=trim(SeedFile),status="unknown")
        write(Unit,*) SeedList(1)
        close(Unit)
      end if

      ! Output
      if (present(Out)) then
          Out=SeedList(1)
      end if
      deallocate(SeedList)
    end subroutine

    !###########################################################################
end module

!###############################################################################
