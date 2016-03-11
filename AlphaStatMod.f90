
!###############################################################################

module AlphaStatMod

  ! Basic statistical stuff

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit

  implicit none

  private
  public :: CalcMoments,MomentsS,MomentsD
  public :: CalcCorrelation,CorrelationS,CorrelationD

  type :: MomentsS
    integer(int32) :: n
    real(real32)   :: Mean
    real(real32)   :: Var
    real(real32)   :: SD
    real(real32)   :: Skew ! TODO: do we ever need skewness?
    real(real32)   :: Curt ! TODO: do we ever need curtosis?
  end type

  type :: MomentsD
    integer(int32) :: n
    real(real64)   :: Mean
    real(real64)   :: Var
    real(real64)   :: SD
    real(real64)   :: Skew
    real(real64)   :: Curt
  end type

  interface CalcMoments
    module procedure CalcMomentsS, CalcMomentsD
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

  interface CalcCorrelation
    module procedure CalcCorrelationS, CalcCorrelationD
  end interface

  contains

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

    subroutine CalcMomentsS(x,Out)
      implicit none

      ! Arguments
      real(real32),intent(in)   :: x(:)
      type(MomentsS),intent(out) :: Out

      ! Other
      integer(int32) i,n
      real(real32) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcMoments"
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
    end subroutine

    !###########################################################################

    subroutine CalcMomentsD(x,Out)
      implicit none

      ! Arguments
      real(real64),intent(in)   :: x(:)
      type(MomentsD),intent(out) :: Out

      ! Other
      integer(int32) i,n

      real(real64) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcMoments"
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
    end subroutine

    !###########################################################################
end module

!###############################################################################
