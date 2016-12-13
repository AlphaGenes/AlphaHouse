
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaStatMod.f90
!
! DESCRIPTION:
!> @brief    Basic statistical functions
!
!> @details  Basic statistical functions such as mean, variance, covariance,
!!           correlation etc. that work with different types
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     September 22, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-22 GGorjanc - Added weights for mean and variance.
!                       Removed skewness and curtosis from DescStat because I
!                       do not know the weighted formulas and we do not use them.
!                       Simplified function names, e.g., CalcMean is now Mean.
! 2016-03-11 GGorjanc - Initial Version
!
!-------------------------------------------------------------------------------
module AlphaStatMod

  use ISO_Fortran_Env, STDIN=>input_unit, STDOUT=>output_unit, STDERR=>error_unit

  implicit none

  private
  ! Types
  public :: DescStatS, DescStatD
  public :: DescStatMatrixS, DescStatMatrixD
  public :: CorrelationS, CorrelationD
  ! Methods
  public :: Mean, Var, StdDev
  public :: DescStat
  public :: DescStatMatrix, DescStatSymMatrix, DescStatLowTriMatrix
  public :: Cov, Cor
  public :: moment, pearsn

  !> @brief Mean interface
  interface Mean
    module procedure MeanS, MeanD
  end interface

  !> @brief Var interface
  interface Var
    module procedure VarS, VarD
  end interface

  !> @brief StdDev interface
  interface StdDev
    module procedure StdDevS, StdDevD
  end interface

  !> @brief DescStat interface
  interface DescStat
    module procedure DescStatSingle, DescStatDouble
  end interface

  !> @brief Descriptive statistics of a variable type (single precision real)
  type :: DescStatS
    integer(int32) :: n
    real(real32)   :: Mean
    real(real32)   :: Var
    real(real32)   :: SD
    real(real32)   :: Min
    real(real32)   :: Max
  end type

  !> @brief Descriptive statistics a variable type (double precision real)
  type :: DescStatD
    integer(int32) :: n
    real(real64)   :: Mean
    real(real64)   :: Var
    real(real64)   :: SD
    real(real64)   :: Min
    real(real64)   :: Max
  end type

  !> @brief DescStatSymMatrix interface
  interface DescStatSymMatrix
    module procedure DescStatSymMatrixSingle, DescStatSymMatrixDouble
  end interface

  !> @brief DescStatLowTriMatrix interface
  interface DescStatLowTriMatrix
    module procedure DescStatSymMatrixSingle, DescStatSymMatrixDouble
  end interface

  !> @brief DescStatMatrix interface
  interface DescStatMatrix
    module procedure DescStatMatrixSingle, DescStatMatrixDouble
  end interface

  !> @brief Descriptive statistics of a matrix type (single precision real)
  type :: DescStatMatrixS
    type(DescStatS) :: All
    type(DescStatS) :: Diag
    type(DescStatS) :: OffDiag
  end type

  !> @brief Descriptive statistics of a matrix type (double precision real)
  type :: DescStatMatrixD
    type(DescStatD) :: All
    type(DescStatD) :: Diag
    type(DescStatD) :: OffDiag
  end type

  !> @brief Covariance interface
  interface Cov
    module procedure CovS, CovD
  end interface

  !> @brief Correlation interface
  interface Cor
    module procedure CorS, CorD
  end interface

  !> @brief Correlation type (single precision real)
  type :: CorrelationS
    real(real32) :: Cor
    real(real32) :: Cov
    real(real32) :: Var1
    real(real32) :: Var2
  end type

  !> @brief Correlation type (double precision real)
  type :: CorrelationD
    real(real64) :: Cor
    real(real64) :: Cov
    real(real64) :: Var1
    real(real64) :: Var2
  end type

  !> @brief Epsilon (single precision real)
  REAL(REAL32), PARAMETER :: EPSILONS = 1.0E-30
  !> @brief Epsilon (double precision real)
  REAL(REAL64), PARAMETER :: EPSILOND = 1.0D-60

  contains

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Mean (single precision real)
    !> @details Compute mean, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function MeanS(x, w) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)           :: x(:) !< values
      real(real32), intent(in), optional :: w(:) !< weights
      real(real32)                       :: Res  !< @return mean

      ! Other
      real(real32) :: SumW

      if (.not.present(w)) then
        Res = sum(x(:)) / size(x)
      else
        if (any(w < 0.0)) then
          write(STDERR, "(a)") "ERROR: Weights must not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILONS) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        Res = sum(x(:) * w(:)) / SumW
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Mean (double precision real)
    !> @details Compute mean, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function MeanD(x, w) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)           :: x(:) !< values
      real(real64), intent(in), optional :: w(:) !< weights
      real(real64)                       :: Res  !< @return mean

      ! Other
      real(real64) :: SumW

      if (.not.present(w)) then
        Res = sum(x(:)) / size(x)
      else
        if (any(w < 0.0)) then
          write(STDERR, "(a)") "ERROR: Weights must not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        Res = sum(x(:) * w(:)) / SumW
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Variance (single precision real)
    !> @details Compute variance, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function VarS(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return variance

      ! Other
      integer(int32) i, n

      real(real32) :: MuIn, Dev, SumW

      n = size(x)
      if (n < 1) then
        write(STDERR, "(a)") "ERROR: number of records must be at least 2 for Var"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (.not.present(Mu)) then
        ! wType does not matter when computing mean
        MuIn = Mean(x, w)
      else
        MuIn = Mu
      end if

      Res = 0.0
      if (.not.present(w)) then
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + Dev * Dev
        end do
        Res = Res / (n - 1)
      else
        if (any(w < 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILONS) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not.present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) /= "count" .and. trim(wType) /= "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + w(i) * Dev * Dev
        end do
        if (trim(wtype) == "count") then
          Res = Res / (SumW - 1.0)
        end if
        if (trim(wtype) == "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Variance (double precision real)
    !> @details Compute variance, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function VarD(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values
      real(real64), intent(in), optional     :: Mu    !< precomputed mean
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real64)                           :: Res   !< @return variance

      ! Other
      integer(int32) i, n

      real(real64) :: MuIn, Dev, SumW

      n = size(x)
      if (n < 1) then
        write(STDERR, "(a)") "ERROR: number of records must be at least 2 for Var"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (.not.present(Mu)) then
        MuIn = Mean(x, w)
      else
        MuIn = Mu
      end if

      Res = 0.0d0
      if (.not.present(w)) then
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + Dev * Dev
        end do
        Res = Res / (n - 1)
      else
        if (any(w < 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not.present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) /= "count" .and. trim(wType) /= "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + w(i) * Dev * Dev
        end do
        if (trim(wtype) == "count") then
          Res = Res / (SumW - 1.0d0)
        end if
        if (trim(wtype) == "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Standard deviation (single precision real)
    !> @details Compute standard deviation, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevS(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return standard deviation

      Res = sqrt(Var(x, Mu, w, wType))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Standard deviation (double precision real)
    !> @details Compute standard deviation, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevD(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values
      real(real64), intent(in), optional     :: Mu    !< precomputed mean
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real64)                           :: Res   !< @return standard deviation

      Res = sqrt(Var(x, Mu, w, wType))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (single precision real)
    !> @details Compute descriptive statistics of a variable, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatSingle(x, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(DescStatS)                        :: Res   !< @return descriptive statistics

      ! Other
      integer(int32) :: n

      n = size(x)
      if (n < 2) then
        write(STDERR, "(a)") "ERROR: number of records must be at least 2 for DescStat"
        write(STDERR, "(a)") " "
        stop 1
      end if

      Res%n    = n
      Res%Mean = Mean(x, w)
      Res%Var  = Var(x, Res%Mean, w, wType)
      Res%SD   = sqrt(Res%Var)
      Res%Min  = minval(x(:))
      Res%Max  = maxval(x(:))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (double precision real)
    !> @details Compute descriptive statistics of a variable, possibly weighted
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatDouble(x, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(DescStatD)                        :: Res   !< @return descriptive statistics

      ! Other
      integer(int32) :: n

      n = size(x)
      if (n < 2) then
        write(STDERR, "(a)") "ERROR: number of records must be at least 2 for DescStat"
        write(STDERR, "(a)") " "
        stop 1
      end if

      Res%n    = n
      Res%Mean = Mean(x, w)
      Res%Var  = Var(x, Res%Mean, w, wType)
      Res%SD   = sqrt(Res%Var)
      Res%Min  = minval(x(:))
      Res%Max  = maxval(x(:))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (single precision real)
    !> @details Compute descriptive statistics of a symetric matrix
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatSymMatrixSingle(x, Diag) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)      :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixS)         :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) i, j, k, n, p
      real(real32), allocatable :: DiagVal(:)
      real(real32), allocatable :: OffDiagVal(:)
      logical                   :: DiagInternal

      n = size(x, 1)
      p = size(x, 2)
      if (n /= p) then
        write(STDERR, "(a)") "ERROR: DescStatSymMatrix work only with symmetric matrices!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (present(Diag)) then
        DiagInternal = Diag
      else
        DiagInternal = .true.
      end if

      ! Diagonal
      if (DiagInternal) then
        allocate(DiagVal(n))
        do i = 1, n
          DiagVal(i) = x(i, i)
        end do
        Res%Diag = DescStat(DiagVal)
      end if

      ! Off-diagonal (lower-triangle only!!!)
      allocate(OffDiagVal(nint(real(n*n)/2-real(n)/2))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k = 0
      do j = 1, (p - 1)
        do i = (j + 1), n
          k = k + 1
          OffDiagVal(k) = x(i, j)
        end do
      end do
      Res%OffDiag = DescStat(OffDiagVal)

      if (DiagInternal) then
        Res%All = DescStat([DiagVal, OffDiagVal])
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (double precision real)
    !> @details Compute descriptive statistics of a symetric matrix
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatSymMatrixDouble(x, Diag) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)      :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixD)        :: Res      !< @return descriptive statistics

      ! Other
      integer(int32) i, j, k, n, p
      real(real64), allocatable :: DiagVal(:)
      real(real64), allocatable :: OffDiagVal(:)
      logical                   :: DiagInternal

      n = size(x, 1)
      p = size(x, 2)
      if (n /= p) then
        write(STDERR, "(a)") "ERROR: DescStatSymMatrix work only with symmetric matrices!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (present(Diag)) then
        DiagInternal = Diag
      else
        DiagInternal = .true.
      end if

      ! Diagonal
      if (DiagInternal) then
        allocate(DiagVal(n))
        do i = 1, n
          DiagVal(i) = x(i, i)
        end do
        Res%Diag = DescStat(DiagVal)
      end if

      ! Off-diagonal (lower-triangle only!!!)
      allocate(OffDiagVal(nint(real(n*n)/2-real(n)/2))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k = 0
      do j = 1, (p - 1)
        do i = (j + 1), n
          k = k + 1
          OffDiagVal(k) = x(i, j)
        end do
      end do
      Res%OffDiag = DescStat(OffDiagVal)

      if (DiagInternal) then
        Res%All = DescStat([DiagVal, OffDiagVal])
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (single precision real)
    !> @details Compute descriptive statistics of a matrix
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixSingle(x) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in) :: x(:, :) !< matrix
      type(DescStatMatrixS)    :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) i, j, k, l, n, p, MinNP
      real(real32), allocatable :: Diag(:)
      real(real32), allocatable :: OffDiag(:)

      n = size(x, 1)
      p = size(x, 2)

      MinNP = minval([n, p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k = 0
      l = 0
      do j = 1, p
        do i = 1, n
          if (i == j) then
            k = k + 1
            Diag(k) = x(i, j)
          else
            l = l + 1
            OffDiag(l) = x(i, j)
          end if
        end do
      end do

      Res%Diag    = DescStat(Diag)
      Res%OffDiag = DescStat(OffDiag)
      Res%All     = DescStat([Diag, OffDiag])
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics (double precision real)
    !> @details Compute descriptive statistics of a matrix
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixDouble(x) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in) :: x(:, :) !< matrix
      type(DescStatMatrixD)    :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) i, j, k, l, n, p, MinNP
      real(real64), allocatable :: Diag(:)
      real(real64), allocatable :: OffDiag(:)

      n = size(x, 1)
      p = size(x, 2)

      MinNP = minval([n, p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k = 0
      l = 0
      do j = 1, p
        do i = 1, n
          if (i == j) then
            k = k + 1
            Diag(k) = x(i, j)
          else
            l = l + 1
            OffDiag(l) = x(i, j)
          end if
        end do
      end do

      Res%Diag    = DescStat(Diag)
      Res%OffDiag = DescStat(OffDiag)
      Res%All     = DescStat([Diag, OffDiag])
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Covariance (single precision real)
    !> @details Compute covariance between two variables
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function CovS(x, y, MuX, MuY, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values for x
      real(real32), intent(in)               :: y(:)  !< values for y
      real(real32), intent(in), optional     :: MuX   !< precomupted mean for x
      real(real32), intent(in), optional     :: MuY   !< precomupted mean for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return covariance

      ! Other
      integer(int32) :: i, n
      real(real32) :: MuXIn, MuYIn, SumW

      n = size(x)

      if (.not.present(MuX)) then
        MuXIn = Mean(x, w)
      else
        MuXIn = MuX
      end if
      if (.not.present(MuY)) then
        MuYIn = Mean(y, w)
      else
        MuYIn = MuY
      end if

      Res = 0.0
      if (.not.present(w)) then
        do i = 1, n
          Res = Res + (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        Res = Res / (n - 1)
      else
        if (any(w < 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not.present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) /= "count" .and. trim(wType) /= "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Res = Res + w(i) * (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        if (trim(wtype) == "count") then
          Res = Res / (SumW - 1.0)
        end if
        if (trim(wtype) == "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Covariance (double precision real)
    !> @details Compute covariance between two variables
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function CovD(x, y, MuX, MuY, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values for x
      real(real64), intent(in)               :: y(:)  !< values for y
      real(real64), intent(in), optional     :: MuX   !< precomupted mean for x
      real(real64), intent(in), optional     :: MuY   !< precomupted mean for y
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real64)                           :: Res   !< @return covariance

      ! Other
      integer(int32) :: i, n
      real(real64) :: MuXIn, MuYIn, SumW

      n = size(x)

      if (.not.present(MuX)) then
        MuXIn = Mean(x, w)
      else
        MuXIn = MuX
      end if
      if (.not.present(MuY)) then
        MuYIn = Mean(y, w)
      else
        MuYIn = MuY
      end if

      Res = 0.0d0
      if (.not.present(w)) then
        do i = 1, n
          Res = Res + (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        Res = Res / (n - 1)
      else
        if (any(w < 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW < EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not.present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) /= "count" .and. trim(wType) /= "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Res = Res + w(i) * (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        if (trim(wtype) == "count") then
          Res = Res / (SumW - 1.0d0)
        end if
        if (trim(wtype) == "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Correlation (single precision real)
    !> @details Compute correlation between two variables
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function CorS(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values for x
      real(real32), intent(in)               :: y(:)  !< values for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationS)                     :: Res   !< @return correlation

      ! Other
      real(real32) :: MuX, MuY

      MuX = Mean(x, w)
      MuY = Mean(y, w)

      Res%Var1 = Var(x, MuX, w, wType)
      Res%Var2 = Var(y, MuY, w, wType)
      Res%Cov  = Cov(x, y, MuX, MuY, w, wType)
      Res%Cor  = Res%Cov / (sqrt(Res%Var1 * Res%Var2) + tiny(x))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Correlation (double precision real)
    !> @details Compute correlation between two variables
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function CorD(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values for x
      real(real64), intent(in)               :: y(:)  !< values for y
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationD)                     :: Res   !< @return correlation

      ! Other
      real(real64) :: MuX, MuY

      MuX = Mean(x, w)
      MuY = Mean(y, w)

      Res%Var1 = Var(x, MuX, w, wType)
      Res%Var2 = Var(y, MuY, w, wType)
      Res%Cov  = Cov(x, y, MuX, MuY, w, wType)
      Res%Cor  = Res%Cov / (sqrt(Res%Var1 * Res%Var2) + tiny(x))
      return
    end function

    !###########################################################################

      subroutine CholDc(a,p,NotPosDefin)    !(DOUBLE PRECISION MODIFIED)
    !SG modified 3/10/15 to return error where user specified matrix is not positive definite
    !DB modified 4/8/2016 to not use Sum as that is a intrinsic procedure
    !DW modified 20/9/2016 to take int64

    use iso_fortran_env, only : int64
    INTEGER(kind=int64) :: n
    real(kind=real64), intent(inout) :: a(:,:),p(:)
    INTEGER(kind=int64) :: i,j,k
    real(kind=real64) :: sum1
    logical, intent(out) :: NotPosDefin

    NotPosDefin = .False.
    n = size(p)
    ! Compute Cholesky factor
    DO i=1,n
      DO j=i,n
        SUM1=a(i,j)
        DO k=i-1,1,-1
          SUM1=SUM1-a(i,k)*a(j,k)
        enddo
        IF(i.eq.j)then
          IF(sum1.le.0.) then
            NotPosDefin = .True.
            return
          endif
          p(i)=sqrt(sum1)
        ELSE
          a(j,i)=sum1/p(i)
        END IF
      enddo
    enddo

    ! Put sd on diagonals so we have complete Cholesky factor
    do i=1,n
      a(i,i)=p(i)
    enddo

  end subroutine choldc
!###############################################################################


  subroutine Pearsn (x,y,n,r)

    implicit none

    integer n
    double precision prob,r,z,x(n),y(n),TINY
    parameter (tiny=1.e-20)
    integer j
    double precision ax,ay,df,sxx,sxy,syy,t,xt,yt
    ! double precision betai

    ax=0.0
    ay=0.0
    DO j=1,n
      ax=ax+x(j)
      ay=ay+y(j)
    END DO
    ax=ax/n                       ! averages
    ay=ay/n

    sxx=0.
    syy=0.
    sxy=0.
    DO j=1,n
      xt=x(j)-ax
      yt=y(j)-ay
      sxx=sxx+xt**2               ! var(x) and var(y)s
      syy=syy+yt**2
      sxy=sxy+xt*yt               ! cov(x,y)
    END DO

    r=sxy/(SQRT(sxx*syy)+TINY)              ! correl coeff
    z=0.5*LOG(((1.+r)+TINY)/((1.-r)+TINY))
    df=n-2
    t=r*SQRT(df/(((1.-r)+TINY)*((1.+r)+TINY)))
    !prob=betai(0.5*df,0.5,df/(df+t**2))
    !prob=erfcc(ABS(z*SQRT(n-1.))/1.4142136)
    prob=0
    return

  end subroutine Pearsn

  !###############################################################################
  SUBROUTINE moment(DATA,n,ave,adev,sdev,var,skew,curt)
    IMPLICIT NONE
    INTEGER n
    DOUBLE PRECISION adev,ave,curt,sdev,skew,var,DATA(n)
    INTEGER j
    DOUBLE PRECISION p,s,ep
    IF (n.le.1) STOP 110003
    s=0
    DO j= 1,n
      s=s+DATA(j)
    END DO

    ave=s/n
    adev=0
    var=0
    skew=0
    curt=0
    ep=0

    DO j=1,n
      s=DATA(j)-ave
      ep=ep+s
      adev=adev+ABS(s)
      p=s*s
      var=var+p
      p=p*s
      skew=skew+p
      p=p*s
      curt=curt+p
    END DO

    adev=adev/n
    var=(var-ep**2/n)/(n-1)
    sdev=SQRT(var)
    IF(var.ne.0)then
      skew=skew/(n*sdev**3)
      curt=curt/(n*var**2)-3
    ELSE
      !PRINT*, 'no skew or kurtosis when zero variance in moment'
      !PAUSE 'no skew or kurtosis when zero variance in moment'
    END IF
    RETURN
  END SUBROUTINE moment
  
!###############################################################################


end module

!###############################################################################
