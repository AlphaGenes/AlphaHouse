
!###############################################################################

module AlphaStatMod

  ! Basic statistical stuff

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit

  implicit none

  private
  public :: CalcMean,CalcVar,CalcSD
  public :: CalcDescStat,DescStatS,DescStatD
  public :: CalcDescStatSymMatrix,DescStatMatrixS,DescStatMatrixD
  public :: CalcCorrelation,CorrelationS,CorrelationD

  interface CalcMean
    module procedure CalcMeanS, CalcMeanD
  end interface

  interface CalcVar
    module procedure CalcVarS, CalcVarD
  end interface

  interface CalcSD
    module procedure CalcSDS, CalcSDD
  end interface

  interface CalcDescStat
    module procedure CalcDescStatS, CalcDescStatD
  end interface

  type :: DescStatS
    integer(int32) :: n
    real(real32)   :: Mean
    real(real32)   :: Var
    real(real32)   :: SD
    real(real32)   :: Skew ! TODO: do we ever need skewness?
    real(real32)   :: Curt ! TODO: do we ever need curtosis?
    real(real32)   :: Min
    real(real32)   :: Max
  end type

  type :: DescStatD
    integer(int32) :: n
    real(real64)   :: Mean
    real(real64)   :: Var
    real(real64)   :: SD
    real(real64)   :: Skew
    real(real64)   :: Curt
    real(real64)   :: Min
    real(real64)   :: Max
  end type

  interface CalcDescStatSymMatrix
    module procedure CalcDescStatSymMatrixS, CalcDescStatSymMatrixD
  end interface

  type :: DescStatMatrixS
    type(DescStatS) :: Diag
    type(DescStatS) :: OffDiag
  end type

  type :: DescStatMatrixD
    type(DescStatD) :: Diag
    type(DescStatD) :: OffDiag
  end type

  interface CalcCorrelation
    module procedure CalcCorrelationS, CalcCorrelationD
  end interface

  type :: CorrelationS
    real(real32)   :: Cor
    real(real32)   :: Cov
    real(real32)   :: Var1
    real(real32)   :: Var2
  end type

  type :: CorrelationD
    real(real64)   :: Cor
    real(real64)   :: Cov
    real(real64)   :: Var1
    real(real64)   :: Var2
  end type

  contains

    !###########################################################################

    subroutine CalcMeanS(x,Out)
      implicit none

      ! Arguments
      real(real32),intent(in)  :: x(:)
      real(real32),intent(out) :: Out

      Out=sum(x(:))/real(size(x))
    end subroutine

    !###########################################################################

    subroutine CalcMeanD(x,Out)
      implicit none

      ! Arguments
      real(real64),intent(in)  :: x(:)
      real(real64),intent(out) :: Out

      Out=sum(x(:))/dble(size(x))
    end subroutine

    !###########################################################################

    subroutine CalcVarS(x,Out,Mean)
      implicit none

      ! Arguments
      real(real32),intent(in)  :: x(:)
      real(real32),intent(out) :: Out
      real(real32),intent(in),optional :: Mean

      ! Other
      integer(int32) i,n

      real(real32) :: nR,MeanIn,Dev

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcVar"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=real(n)
      if (.not.present(Mean)) then
        MeanIn=sum(x(:))/nR
      else
        MeanIn=Mean
      end if
      Out=0.0

      do i=1,n
        Dev=x(i)-MeanIn
        Out=Out+Dev*Dev
      end do

      Out=Out/(nR-1.0)
    end subroutine

    !###########################################################################

    subroutine CalcVarD(x,Out,Mean)
      implicit none

      ! Arguments
      real(real64),intent(in)  :: x(:)
      real(real64),intent(out) :: Out
      real(real64),intent(in),optional :: Mean

      ! Other
      integer(int32) i,n

      real(real64) :: nR,MeanIn,Dev

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcVar"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=dble(n)
      if (.not.present(Mean)) then
        MeanIn=sum(x(:))/nR
      else
        MeanIn=Mean
      end if
      Out=0.0d0

      do i=1,n
        Dev=x(i)-MeanIn
        Out=Out+Dev*Dev
      end do

      Out=Out/(nR-1.0d0)
    end subroutine

    !###########################################################################

    subroutine CalcSDS(x,Out,Mean)
      implicit none

      ! Arguments
      real(real32),intent(in)  :: x(:)
      real(real32),intent(out) :: Out
      real(real32),intent(in),optional :: Mean

      if (present(Mean)) then
        call CalcVarS(x,Out,Mean)
      else
        call CalcVarS(x,Out)
      end if
      Out=sqrt(Out)
    end subroutine

    !###########################################################################

    subroutine CalcSDD(x,Out,Mean)
      implicit none

      ! Arguments
      real(real64),intent(in)  :: x(:)
      real(real64),intent(out) :: Out
      real(real64),intent(in),optional :: Mean

      if (present(Mean)) then
        call CalcVarD(x,Out,Mean)
      else
        call CalcVarD(x,Out)
      end if
      Out=sqrt(Out)
    end subroutine

    !###########################################################################

    subroutine CalcDescStatS(x,Out)
      implicit none

      ! Arguments
      real(real32),intent(in)   :: x(:)
      type(DescStatS),intent(out) :: Out

      ! Other
      integer(int32) i,n
      real(real32) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcDescStat"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=real(n)
      Out%Mean=sum(x(:))/nR
      Out%Var=0.0
      Out%Skew=0.0
      Out%Curt=0.0

      do i=1,n
        Dev=x(i)-Out%Mean
        Mul=Dev*Dev
        Out%Var=Out%Var+Mul
        Mul=Mul*Dev
        Out%Skew=Out%Skew+Mul
        Mul=Mul*Dev
        Out%Curt=Out%Curt+Mul
      end do

      Out%Var=Out%Var/(nR-1.0)
      Out%SD=sqrt(Out%Var)
      if (Out%Var > 0.0) then
        Out%Skew=Out%Skew/(nR*Out%SD**3)
        Out%Curt=Out%Curt/(nR*Out%Var**2)-3.0
      end if

      Out%Min=minval(x(:))
      Out%Max=maxval(x(:))
    end subroutine

    !###########################################################################

    subroutine CalcDescStatD(x,Out)
      implicit none

      ! Arguments
      real(real64),intent(in)   :: x(:)
      type(DescStatD),intent(out) :: Out

      ! Other
      integer(int32) i,n

      real(real64) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcDescStat"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=dble(n)
      Out%Mean=sum(x(:))/nR
      Out%Var=0.0d0
      Out%Skew=0.0d0
      Out%Curt=0.0d0

      do i=1,n
        Dev=x(i)-Out%Mean
        Mul=Dev*Dev
        Out%Var=Out%Var+Mul
        Mul=Mul*Dev
        Out%Skew=Out%Skew+Mul
        Mul=Mul*Dev
        Out%Curt=Out%Curt+Mul
      end do

      Out%Var=Out%Var/(nR-1.0d0)
      Out%SD=sqrt(Out%Var)
      if (Out%Var > 0.0d0) then
        Out%Skew=Out%Skew/(nR*Out%SD**3)
        Out%Curt=Out%Curt/(nR*Out%Var**2)-3.0d0
      end if

      Out%Min=minval(x(:))
      Out%Max=maxval(x(:))
    end subroutine

    !###########################################################################

    subroutine CalcDescStatSymMatrixS(x,Out)
      implicit none

      ! Arguments
      real(real32),intent(in)           :: x(:,:)
      type(DescStatMatrixS),intent(out) :: Out

      ! Other
      integer(int32) i,j,k,n,p!,MinNP
      real(real32),allocatable :: Diag(:)
      real(real32),allocatable :: OffDiag(:)

      n=size(x,1)
      p=size(x,2)
      if (n /= p) then
        write(STDERR,"(a)") "ERROR: CalcDescStatMatrix work only with symmetric matrices!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! Start of an attempt to work with non-symmetric matrices
      ! MinNP=min([n,p])
      ! allocate(Diag(MinNP))

      allocate(Diag(n))
      do i=1,n
        Diag(i)=x(i,i)
      end do
      call CalcDescStatS(Diag,Out%Diag)
      deallocate(Diag)

      allocate(OffDiag(nint(real(n*n)/2.0-real(n)/2.0))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k=0
      do j=1,(p-1)
        do i=(j+1),n
          k=k+1
          OffDiag(k)=x(i,j)
        end do
      end do
      call CalcDescStatS(OffDiag,Out%OffDiag)
      deallocate(OffDiag)
    end subroutine

    !###########################################################################

    subroutine CalcDescStatSymMatrixD(x,Out)
      implicit none

      ! Arguments
      real(real64),intent(in)           :: x(:,:)
      type(DescStatMatrixD),intent(out) :: Out

      ! Other
      integer(int32) i,j,k,n,p!,MinNP
      real(real64),allocatable :: Diag(:)
      real(real64),allocatable :: OffDiag(:)

      n=size(x,1)
      p=size(x,2)
      if (n /= p) then
        write(STDERR,"(a)") "ERROR: CalcDescStatMatrix work only with symmetric matrices!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      ! Start of an attempt to work with non-symmetric matrices
      ! MinNP=min([n,p])
      ! allocate(Diag(MinNP))

      allocate(Diag(n))
      do i=1,n
        Diag(i)=x(i,i)
      end do
      call CalcDescStatD(Diag,Out%Diag)
      deallocate(Diag)

      allocate(OffDiag(nint(real(n*n)/2.0-real(n)/2.0))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k=0
      do j=1,(p-1)
        do i=(j+1),n
          k=k+1
          OffDiag(k)=x(i,j)
        end do
      end do
      call CalcDescStatD(OffDiag,Out%OffDiag)
      deallocate(OffDiag)
    end subroutine

    !###########################################################################
    subroutine CalcCorrelationS(x,y,Out)
      implicit none
      ! Arguments
      real(real32),intent(in)        :: x(:)
      real(real32),intent(in)        :: y(:)
      type(CorrelationS),intent(out) :: Out

      ! Other
      integer(int32) :: i,n
      real(real32) :: nR,MeanX,MeanY,DevX,DevY,SumXX,SumXY,SumYY

      n=size(x)
      nR=real(n)
      MeanX=sum(x(:))/nR
      MeanY=sum(y(:))/nR
      SumXX=0.0
      SumXY=0.0
      SumYY=0.0
      do i=1,n
        DevX=x(i)-MeanX
        DevY=y(i)-MeanY
        SumXX=SumXX+DevX*DevX
        SumXY=SumXY+DevX*DevY
        SumYY=SumYY+DevY*DevY
      end do
      Out%Cor=SumXY/(sqrt(SumXX*SumYY)+tiny(x))
      Out%Cov=SumXY/(nR-1.0)
      Out%Var1=SumXX/(nR-1.0)
      Out%Var2=SumYY/(nR-1.0)
    end subroutine

    !###########################################################################

    subroutine CalcCorrelationD(x,y,Out)
      implicit none
      ! Arguments
      real(real64),intent(in)        :: x(:)
      real(real64),intent(in)        :: y(:)
      type(CorrelationD),intent(out) :: Out

      ! Other
      integer(int32) :: i,n
      real(real64) :: nR,MeanX,MeanY,DevX,DevY,SumXX,SumXY,SumYY

      n=size(x)
      nR=dble(n)
      MeanX=sum(x(:))/nR
      MeanY=sum(y(:))/nR
      SumXX=0.0d0
      SumXY=0.0d0
      SumYY=0.0d0
      do i=1,n
        DevX=x(i)-MeanX
        DevY=y(i)-MeanY
        SumXX=SumXX+DevX*DevX
        SumXY=SumXY+DevX*DevY
        SumYY=SumYY+DevY*DevY
      end do
      Out%Cor=SumXY/(sqrt(SumXX*SumYY)+tiny(x))
      Out%Cov=SumXY/(nR-1.0d0)
      Out%Var1=SumXX/(nR-1.0d0)
      Out%Var2=SumYY/(nR-1.0d0)
    end subroutine

    !###########################################################################
end module

!###############################################################################
