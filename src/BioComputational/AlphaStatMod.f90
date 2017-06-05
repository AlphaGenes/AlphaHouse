
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
!> @date     January 18, 2017
!
!> @version  0.0.2 (alpha)
!
! REVISION HISTORY:
! 2016-09-22 GGorjanc - IsMissing, IsAnyMissingMatrix, RemoveMissing, and RemoveAnyMissingMatrix added
!                       Added support for Int8 and Int32 arrays to all functions
!                       Renamed types to Type{Int8,Int32,Real32,Real64}
!                       Renamed methods to Method{I8,I32,R32,R64} - otherwise there is a type/method name clash
! 2016-09-22 GGorjanc - Added weights for mean and variance.
!                       Removed skewness and curtosis from DescStat because I
!                       do not know the weighted formulas and we do not use them.
!                       Simplified function names, e.g., CalcMean is now Mean.
! 2016-03-11 GGorjanc - Initial Version
!
!-------------------------------------------------------------------------------
module AlphaStatMod

  use ISO_Fortran_Env, STDIN=>input_unit, STDOUT=>output_unit, STDERR=>error_unit
  use, intrinsic :: IEEE_Arithmetic

  implicit none

  private
  ! Types
  public :: DescStatReal32, DescStatReal64
  public :: DescStatMatrixReal32, DescStatMatrixReal64
  public :: CorrelationReal32, CorrelationReal64
  ! Methods
  public :: IsMissing, IsAnyMissingMatrix, RemoveMissing, RemoveAnyMissingMatrix
  public :: Mean, Var, StdDev
  public :: DescStat
  public :: DescStatMatrix, DescStatLowTriMatrix
  public :: Cov, Cor
  public :: moment, pearsn

  !> @brief IsMissing interface
  interface IsMissing
    module procedure IsMissingI8, IsMissingI32, IsMissingR32, IsMissingR64
  end interface

  !> @brief IsAnyMissingMatrix interface
  interface IsAnyMissingMatrix
    module procedure IsAnyMissingMatrixI8, IsAnyMissingMatrixI32, IsAnyMissingMatrixR32, IsAnyMissingMatrixR64
  end interface

  !> @brief RemoveMissing interface
  interface RemoveMissing
    module procedure RemoveMissingI8, RemoveMissingI32, RemoveMissingR32, RemoveMissingR64
  end interface

  !> @brief RemoveAnyMissingMatrix interface
  interface RemoveAnyMissingMatrix
    module procedure RemoveAnyMissingMatrixI8, RemoveAnyMissingMatrixI32, RemoveAnyMissingMatrixR32, RemoveAnyMissingMatrixR64
  end interface

  !> @brief Mean interface
  interface Mean
    module procedure MeanI8, MeanI32, MeanR32, MeanR64
  end interface

  !> @brief Var interface
  interface Var
    module procedure VarI8, VarI32, VarR32, VarR64
  end interface

  !> @brief StdDev interface
  interface StdDev
    module procedure StdDevI8, StdDevI32, StdDevR32, StdDevR64
  end interface

  !> @brief DescStat interface
  interface DescStat
    module procedure DescStatI8, DescStatI32, DescStatR32, DescStatR64
  end interface

  !> @brief Descriptive statistics of a variable type - real32
  type :: DescStatReal32
    integer(int32) :: n
    real(real32)   :: Mean
    real(real32)   :: Var
    real(real32)   :: SD
    real(real32)   :: Min
    real(real32)   :: Max
  end type

  !> @brief Descriptive statistics a variable type - real64
  type :: DescStatReal64
    integer(int32) :: n
    real(real64)   :: Mean
    real(real64)   :: Var
    real(real64)   :: SD
    real(real64)   :: Min
    real(real64)   :: Max
  end type

  !> @brief DescStatLowTriMatrix interface
  interface DescStatLowTriMatrix
    module procedure DescStatLowTriMatrixI8, DescStatLowTriMatrixI32, DescStatLowTriMatrixR32, DescStatLowTriMatrixR64
  end interface

  !> @brief DescStatMatrix interface
  interface DescStatMatrix
    module procedure DescStatMatrixI8, DescStatMatrixI32, DescStatMatrixR32, DescStatMatrixR64
  end interface

  !> @brief Descriptive statistics of a matrix type - real32
  type :: DescStatMatrixReal32
    type(DescStatReal32) :: All
    type(DescStatReal32) :: Diag
    type(DescStatReal32) :: OffDiag
  end type

  !> @brief Descriptive statistics of a matrix type - real64
  type :: DescStatMatrixReal64
    type(DescStatReal64) :: All
    type(DescStatReal64) :: Diag
    type(DescStatReal64) :: OffDiag
  end type

  !> @brief Covariance interface
  interface Cov
    module procedure CovI8, CovI32, CovR32, CovR64
  end interface

  !> @brief Correlation interface
  interface Cor
    module procedure CorI8, CorI32, CorR32, CorR64
  end interface

  !> @brief Correlation type - real32
  type :: CorrelationReal32
    real(real32) :: Cor
    real(real32) :: Cov
    real(real32) :: Var1
    real(real32) :: Var2
  end type

  !> @brief Correlation type - real64
  type :: CorrelationReal64
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
    !> @brief  Test for missingness - int8
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure elemental function IsMissingI8(x, MissingValue) result(Res)
      implicit none
      integer(int8), intent(in) :: x         !< values
      integer(int8), intent(in) :: MissingValue !< missing value representation
      logical                   :: Res !< @return test for missingness
      Res = x .eq. MissingValue
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for missingness - int32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure elemental function IsMissingI32(x, MissingValue) result(Res)
      implicit none
      integer(int32), intent(in) :: x         !< values
      integer(int32), intent(in) :: MissingValue !< missing value representation
      logical                    :: Res !< @return test for missingness
      Res = x .eq. MissingValue
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for missingness - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure elemental function IsMissingR32(x, MissingValue) result(Res)
      implicit none
      real(real32), intent(in) :: x        !< values
      real(real32), intent(in) :: MissingValue !< missing value representation
      logical                  :: Res !< @return test for missingness
      Res = x .eq. MissingValue
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for missingness - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure elemental function IsMissingR64(x, MissingValue) result(Res)
      implicit none
      real(real64), intent(in) :: x         !< values
      real(real64), intent(in) :: MissingValue !< missing value representation
      logical                  :: Res !< @return test for missingness
      Res = x .eq. MissingValue
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for across row/column missingness in a matrix - int8,
    !> @author Andrew Whalen, awhalen@roslin.ed.ac.uk & Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 18, 2017
    !---------------------------------------------------------------------------
    pure function IsAnyMissingMatrixI8(x, MissingValue, Dim) result(Res)
      implicit none
      integer(int8), intent(in)          :: x(:, :)      !< values
      integer(int8), intent(in)          :: MissingValue !< missing value representation
      integer, intent(in)                :: Dim          !< Dim = 1 test along dimension 1 (i.e. missingnes within each column); Dim = 2, test along dimension 2. (i.e. missingnes within each row)
      logical, allocatable, dimension(:) :: Res          !< @return test for missingness (a 1D array)
      integer :: i, n
      if (Dim .eq. 1) then ! n is the size of the remaining dimension once Dim has been removed.
        n = size(x, 2)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(:, i) .eq. MissingValue)
        end do
      else if (Dim .eq. 2) then
        n = size(x, 1)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(i, :) .eq. MissingValue)
        end do
      endif
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for across row/column missingness in a matrix - int32,
    !> @author Andrew Whalen, awhalen@roslin.ed.ac.uk & Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 18, 2017
    !---------------------------------------------------------------------------
    pure function IsAnyMissingMatrixI32(x, MissingValue, Dim) result(Res)
      implicit none
      integer(int32), intent(in)         :: x(:, :)      !< values
      integer(int32), intent(in)         :: MissingValue !< missing value representation
      integer, intent(in)                :: Dim          !< Dim = 1 test along dimension 1 (i.e. missingnes within each column); Dim = 2, test along dimension 2. (i.e. missingnes within each row)
      logical, allocatable, dimension(:) :: Res          !< @return test for missingness (a 1D array)
      integer :: i, n
      if (Dim .eq. 1) then ! n is the size of the remaining dimension once Dim has been removed.
        n = size(x, 2)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(:, i) .eq. MissingValue)
        end do
      else if (Dim .eq. 2) then
        n = size(x, 1)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(i, :) .eq. MissingValue)
        end do
      endif
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for across row/column missingness in a matrix - real32,
    !> @author Andrew Whalen, awhalen@roslin.ed.ac.uk & Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 18, 2017
    !---------------------------------------------------------------------------
    pure function IsAnyMissingMatrixR32(x, MissingValue, Dim) result(Res)
      implicit none
      real(real32), intent(in)           :: x(:, :)      !< values
      real(real32), intent(in)           :: MissingValue !< missing value representation
      integer, intent(in)                :: Dim          !< Dim = 1 test along dimension 1 (i.e. missingnes within each column); Dim = 2, test along dimension 2. (i.e. missingnes within each row)
      logical, allocatable, dimension(:) :: Res          !< @return test for missingness (a 1D array)
      integer :: i, n
      if (Dim .eq. 1) then ! n is the size of the remaining dimension once Dim has been removed.
        n = size(x, 2)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(:, i) .eq. MissingValue)
        end do
      else if (Dim .eq. 2) then
        n = size(x, 1)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(i, :) .eq. MissingValue)
        end do
      endif
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Test for across row/column missingness in a matrix - real64,
    !> @author Andrew Whalen, awhalen@roslin.ed.ac.uk & Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 18, 2017
    !---------------------------------------------------------------------------
    pure function IsAnyMissingMatrixR64(x, MissingValue, Dim) result(Res)
      implicit none
      real(real64), intent(in)           :: x(:, :)      !< values
      real(real64), intent(in)           :: MissingValue !< missing value representation
      integer, intent(in)                :: Dim          !< Dim = 1 test along dimension 1 (i.e. missingnes within each column); Dim = 2, test along dimension 2. (i.e. missingnes within each row)
      logical, allocatable, dimension(:) :: Res          !< @return test for missingness (a 1D array)
      integer :: i, n
      if (Dim .eq. 1) then ! n is the size of the remaining dimension once Dim has been removed.
        n = size(x, 2)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(:, i) .eq. MissingValue)
        end do
      else if (Dim .eq. 2) then
        n = size(x, 1)
        allocate(Res(n))
        do i = 1, n
          Res(i) = any(x(i, :) .eq. MissingValue)
        end do
      endif
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain a vector without "missing" values - int8
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveMissingI8(x, MissingValue) result(Res)
      implicit none
      integer(int8), intent(in)                :: x(:)         !< values
      integer(int8), intent(in)                :: MissingValue !< missing value representation
      integer(int8), allocatable, dimension(:) :: Res          !< @return x without missing values
      logical :: NotMissing(size(x))
      integer :: nNotMissing, i, j

      NotMissing = .not. IsMissing(x, MissingValue)
      nNotMissing = count(NotMissing)
      allocate(Res(nNotMissing))
      j = 0
      do i = 1, size(x)
        if (NotMissing(i)) then
          j = j + 1
          Res(j) = x(i)
        end if
      enddo
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain a vector without "missing" values - int32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveMissingI32(x, MissingValue) result(Res)
      implicit none
      integer(int32), intent(in)                :: x(:)         !< values
      integer(int32), intent(in)                :: MissingValue !< missing value representation
      integer(int32), allocatable, dimension(:) :: Res          !< @return x without missing values
      logical :: NotMissing(size(x))
      integer :: nNotMissing, i, j

      NotMissing = .not. IsMissing(x, MissingValue)
      nNotMissing = count(NotMissing)
      allocate(Res(nNotMissing))
      j = 0
      do i = 1, size(x)
        if (NotMissing(i)) then
          j = j + 1
          Res(j) = x(i)
        end if
      enddo
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain a vector without "missing" values - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveMissingR32(x, MissingValue) result(Res)
      implicit none
      real(real32), intent(in)                :: x(:)         !< values
      real(real32), intent(in)                :: MissingValue !< missing value representation
      real(real32), allocatable, dimension(:) :: Res          !< @return x without missing values
      logical :: NotMissing(size(x))
      integer :: nNotMissing, i, j

      NotMissing = .not. IsMissing(x, MissingValue)
      nNotMissing = count(NotMissing)
      allocate(Res(nNotMissing))
      j = 0
      do i = 1, size(x)
        if (NotMissing(i)) then
          j = j + 1
          Res(j) = x(i)
        end if
      enddo
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain a vector without "missing" values - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveMissingR64(x, MissingValue) result(Res)
      implicit none
      real(real64), intent(in)                :: x(:)         !< values
      real(real64), intent(in)                :: MissingValue !< missing value representation
      real(real64), allocatable, dimension(:) :: Res          !< @return x without missing values
      logical :: NotMissing(size(x))
      integer :: nNotMissing, i, j

      NotMissing = .not. IsMissing(x, MissingValue)
      nNotMissing = count(NotMissing)
      allocate(Res(nNotMissing))
      j = 0
      do i = 1, size(x)
        if (NotMissing(i)) then
          j = j + 1
          Res(j) = x(i)
        end if
      enddo
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain an array without partially missing rows or columns - int8
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveAnyMissingMatrixI8(x, MissingValue, Dim) result(Res)
      implicit none
      integer(int8), intent(in)                   :: x(:, :)      !< values
      integer(int8), intent(in)                   :: MissingValue !< missing value representation
      integer, intent(in)                         :: Dim          !< Dim = 1 remove columns that contain missing values, if Dim = 2 remove rows that contain missing values
      integer(int8), allocatable, dimension(:, :) :: Res          !< @return x with columns/rows removed
      logical :: NotMissing(size(x, Dim))
      integer :: nNotMissing, i, j
      NotMissing = .not. IsAnyMissingMatrix(x, MissingValue, Dim)
      nNotMissing = count(NotMissing)
      if (Dim .eq. 1) then
        allocate(Res(size(x, 1), nNotMissing))
        j = 0
        do i = 1, size(x, 2)
          if (NotMissing(i)) then
            j = j + 1
            Res(:, j) = x(:, i)
          end if
        enddo
      else if (Dim .eq. 2) then
        allocate(Res(nNotMissing, size(x, 2)))
        j = 0
        do i = 1, size(x, 1)
          if (NotMissing(i)) then
            j = j + 1
            Res(j, :) = x(i, :)
          end if
        enddo
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain an array without partially missing rows or columns - int32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveAnyMissingMatrixI32(x, MissingValue, Dim) result(Res)
      implicit none
      integer(int32), intent(in)                   :: x(:, :)      !< values
      integer(int32), intent(in)                   :: MissingValue !< missing value representation
      integer, intent(in)                          :: Dim          !< Dim = 1 remove columns that contain missing values, if Dim = 2 remove rows that contain missing values
      integer(int32), allocatable, dimension(:, :) :: Res          !< @return x with columns/rows removed
      logical :: NotMissing(size(x, Dim))
      integer :: nNotMissing, i, j
      NotMissing = .not. IsAnyMissingMatrix(x, MissingValue, Dim)
      nNotMissing = count(NotMissing)
      if (Dim .eq. 1) then
        allocate(Res(size(x, 1), nNotMissing))
        j = 0
        do i = 1, size(x, 2)
          if (NotMissing(i)) then
            j = j + 1
            Res(:, j) = x(:, i)
          end if
        enddo
      else if (Dim .eq. 2) then
        allocate(Res(nNotMissing, size(x, 2)))
        j = 0
        do i = 1, size(x, 1)
          if (NotMissing(i)) then
            j = j + 1
            Res(j, :) = x(i, :)
          end if
        enddo
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain an array without partially missing rows or columns - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveAnyMissingMatrixR32(x, MissingValue, Dim) result(Res)
      implicit none
      real(real32), intent(in)                   :: x(:, :)      !< values
      real(real32), intent(in)                   :: MissingValue !< missing value representation
      integer, intent(in)                        :: Dim          !< Dim = 1 remove columns that contain missing values, if Dim = 2 remove rows that contain missing values
      real(real32), allocatable, dimension(:, :) :: Res          !< @return x with columns/rows removed
      logical :: NotMissing(size(x, Dim))
      integer :: nNotMissing, i, j
      NotMissing = .not. IsAnyMissingMatrix(x, MissingValue, Dim)
      nNotMissing = count(NotMissing)
      if (Dim .eq. 1) then
        allocate(Res(size(x, 1), nNotMissing))
        j = 0
        do i = 1, size(x, 2)
          if (NotMissing(i)) then
            j = j + 1
            Res(:, j) = x(:, i)
          end if
        enddo
      else if (Dim .eq. 2) then
        allocate(Res(nNotMissing, size(x, 2)))
        j = 0
        do i = 1, size(x, 1)
          if (NotMissing(i)) then
            j = j + 1
            Res(j, :) = x(i, :)
          end if
        enddo
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Obtain an array without partially missing rows or columns - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   January 17, 2017
    !---------------------------------------------------------------------------
    pure function RemoveAnyMissingMatrixR64(x, MissingValue, Dim) result(Res)
      implicit none
      real(real64), intent(in)                   :: x(:, :)      !< values
      real(real64), intent(in)                   :: MissingValue !< missing value representation
      integer, intent(in)                        :: Dim          !< Dim = 1 remove columns that contain missing values, if Dim = 2 remove rows that contain missing values
      real(real64), allocatable, dimension(:, :) :: Res          !< @return x with columns/rows removed
      logical :: NotMissing(size(x, Dim))
      integer :: nNotMissing, i, j
      NotMissing = .not. IsAnyMissingMatrix(x, MissingValue, Dim)
      nNotMissing = count(NotMissing)
      if (Dim .eq. 1) then
        allocate(Res(size(x, 1), nNotMissing))
        j = 0
        do i = 1, size(x, 2)
          if (NotMissing(i)) then
            j = j + 1
            Res(:, j) = x(:, i)
          end if
        enddo
      else if (Dim .eq. 2) then
        allocate(Res(nNotMissing, size(x, 2)))
        j = 0
        do i = 1, size(x, 1)
          if (NotMissing(i)) then
            j = j + 1
            Res(j, :) = x(i, :)
          end if
        enddo
      end if
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute mean, possibly weighted - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function MeanI8(x, w) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)          :: x(:) !< values
      real(real32), intent(in), optional :: w(:) !< weights
      real(real32)                       :: Res  !< @return mean

      ! Other
      real(real32) :: xR32(size(x))

      xR32 = real(x, kind=real32)
      Res = MeanR32(x=xR32, w=w)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute mean, possibly weighted - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function MeanI32(x, w) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in)         :: x(:) !< values
      real(real32), intent(in), optional :: w(:) !< weights
      real(real32)                       :: Res  !< @return mean

      ! Other
      real(real32) :: xR32(size(x))

      xR32 = real(x, kind=real32)
      Res = MeanR32(x=xR32, w=w)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute mean, possibly weighted - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function MeanR32(x, w) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)           :: x(:) !< values
      real(real32), intent(in), optional :: w(:) !< weights
      real(real32)                       :: Res  !< @return mean

      ! Other
      real(real32) :: SumW

      if (.not. present(w)) then
        Res = sum(x) / size(x)
      else
        if (any(w .lt. 0.0)) then
          write(STDERR, "(a)") "ERROR: Weights must not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILONS) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        Res = sum(x * w) / SumW
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute mean, possibly weighted - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function MeanR64(x, w) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)           :: x(:) !< values
      real(real64), intent(in), optional :: w(:) !< weights
      real(real64)                       :: Res  !< @return mean

      ! Other
      real(real64) :: SumW

      if (.not. present(w)) then
        Res = sum(x) / size(x)
      else
        if (any(w .lt. 0.0)) then
          write(STDERR, "(a)") "ERROR: Weights must not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        Res = sum(x * w) / SumW
      end if
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute variance, possibly weighted - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function VarI8(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)              :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return variance

      ! Other
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = VarR32(x=xR32, Mu=Mu, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute variance, possibly weighted - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function VarI32(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in)             :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return variance

      ! Other
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = VarR32(x=xR32, Mu=Mu, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute variance, possibly weighted - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function VarR32(x, Mu, w, wType) result(Res)
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
      if (n .lt. 2) then
        Res = IEEE_Value(x=Res, class=IEEE_Quiet_NaN)
        return
      end if

      if (.not. present(Mu)) then
        ! wType does not matter when computing mean
        MuIn = Mean(x, w)
      else
        MuIn = Mu
      end if

      Res = 0.0
      if (.not. present(w)) then
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + Dev * Dev
        end do
        Res = Res / (n - 1)
      else
        if (any(w .lt. 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILONS) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not. present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) .ne. "count" .and. trim(wType) .ne. "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + w(i) * Dev * Dev
        end do
        if (trim(wtype) .eq. "count") then
          Res = Res / (SumW - 1.0)
        end if
        if (trim(wtype) .eq. "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute variance, possibly weighted - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function VarR64(x, Mu, w, wType) result(Res)
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
      if (n .lt. 2) then
        Res = IEEE_Value(x=Res, class=IEEE_Quiet_NaN)
        return
      end if

      if (.not. present(Mu)) then
        MuIn = Mean(x, w)
      else
        MuIn = Mu
      end if

      Res = 0.0d0
      if (.not. present(w)) then
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + Dev * Dev
        end do
        Res = Res / (n - 1)
      else
        if (any(w .lt. 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not. present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) .ne. "count" .and. trim(wType) .ne. "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Dev = x(i) - MuIn
          Res = Res + w(i) * Dev * Dev
        end do
        if (trim(wtype) .eq. "count") then
          Res = Res / (SumW - 1.0d0)
        end if
        if (trim(wtype) .eq. "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute standard deviation, possibly weighted - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevI8(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)              :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return standard deviation
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = sqrt(VarR32(x=xR32, Mu=Mu, w=w, wType=wType))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute standard deviation, possibly weighted - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevI32(x, Mu, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in)             :: x(:)  !< values
      real(real32), intent(in), optional     :: Mu    !< precomputed mean
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return standard deviation
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = sqrt(VarR32(x=xR32, Mu=Mu, w=w, wType=wType))
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute standard deviation, possibly weighted - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevR32(x, Mu, w, wType) result(Res)
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
    !> @brief  Compute standard deviation, possibly weighted - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function StdDevR64(x, Mu, w, wType) result(Res)
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
    !> @brief   Compute descriptive statistics of a variable, possibly weighted - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function DescStatI8(x, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)              :: x(:)    !< values
      real(real32), intent(in), optional     :: w(:)    !< weights
      character(len=*), intent(in), optional :: wType   !< type of weights (count, freq)
      type(DescStatReal32)                   :: Res     !< @return descriptive statistics

      ! Other
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = DescStatR32(x=xR32, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Compute descriptive statistics of a variable, possibly weighted - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function DescStatI32(x, w, wType) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in)             :: x(:)    !< values
      real(real32), intent(in), optional     :: w(:)    !< weights
      character(len=*), intent(in), optional :: wType   !< type of weights (count, freq)
      type(DescStatReal32)                   :: Res     !< @return descriptive statistics

      ! Other
      real(real32) :: xR32(size(x))
      xR32 = real(x, kind=real32)
      Res = DescStatR32(x=xR32, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute descriptive statistics of a variable, possibly weighted - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatR32(x, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)               :: x(:)    !< values
      real(real32), intent(in), optional     :: w(:)    !< weights
      character(len=*), intent(in), optional :: wType   !< type of weights (count, freq)
      type(DescStatReal32)                   :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) :: n

      n = size(x)
      Res%n    = n
      Res%Mean = Mean(x, w)
      Res%Var  = Var(x, Res%Mean, w, wType)
      Res%SD   = sqrt(Res%Var)
      Res%Min  = minval(x)
      Res%Max  = maxval(x)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Compute descriptive statistics of a variable, possibly weighted - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatR64(x, w, wType) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)               :: x(:)    !< values
      real(real64), intent(in), optional     :: w(:)    !< weights
      character(len=*), intent(in), optional :: wType   !< type of weights (count, freq)
      type(DescStatReal64)                   :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) :: n

      n = size(x)
      Res%n    = n
      Res%Mean = Mean(x, w)
      Res%Var  = Var(x, Res%Mean, w, wType)
      Res%SD   = sqrt(Res%Var)
      Res%Min  = minval(x)
      Res%Max  = maxval(x)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics of a lower-triangle matrix - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatLowTriMatrixI8(x, Diag) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)     :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixReal32)    :: Res     !< @return descriptive statistics

      ! Other
      real(real32), allocatable, dimension(:, :) :: xR32

      allocate(xR32(size(x, 1), size(x, 2)))
      xR32 = real(x, kind=real32)
      Res = DescStatLowTriMatrixR32(x=xR32, Diag=Diag)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics of a lower-triangle matrix - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatLowTriMatrixI32(x, Diag) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in)    :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixReal32)    :: Res     !< @return descriptive statistics

      ! Other
      real(real32), allocatable, dimension(:, :) :: xR32
      allocate(xR32(size(x, 1), size(x, 2)))
      xR32 = real(x, kind=real32)
      Res = DescStatLowTriMatrixR32(x=xR32, Diag=Diag)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Descriptive statistics of a lower-triangle matrix - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatLowTriMatrixR32(x, Diag) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)      :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixReal32)    :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) i, j, k, n, p
      real(real32), allocatable, dimension(:) :: DiagVal, OffDiagVal
      logical :: DiagInternal

      n = size(x, 1)
      p = size(x, 2)

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
    !> @brief  Descriptive statistics of a lower-triangle matrix - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatLowTriMatrixR64(x, Diag) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)      :: x(:, :) !< matrix
      logical, intent(in), optional :: Diag    !< should diagonal be summarized (default=.true.)
      type(DescStatMatrixReal64)    :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) :: i, j, k, n, p
      real(real64), allocatable, dimension(:) :: DiagVal, OffDiagVal
      logical :: DiagInternal

      n = size(x, 1)
      p = size(x, 2)

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
    !> @brief   Descriptive statistics of a matrix - int8
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixI8(x) result(Res)
      implicit none

      ! Arguments
      integer(int8), intent(in)  :: x(:, :) !< matrix
      type(DescStatMatrixReal32) :: Res     !< @return descriptive statistics

      ! Other
      real(real32), allocatable, dimension(:, :) :: xR32
      allocate(xR32(size(x, 1), size(x, 2)))
      xR32 = real(x, kind=real32)
      Res = DescStatMatrixR32(x=xR32)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Descriptive statistics of a matrix - int32
    !> @details x is cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixI32(x) result(Res)
      implicit none

      ! Arguments
      integer(int32), intent(in) :: x(:, :) !< matrix
      type(DescStatMatrixReal32) :: Res     !< @return descriptive statistics

      ! Other
      real(real32), allocatable, dimension(:, :) :: xR32
      allocate(xR32(size(x, 1), size(x, 2)))
      xR32 = real(x, kind=real32)
      Res = DescStatMatrixR32(x=xR32)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Descriptive statistics of a matrix - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixR32(x) result(Res)
      implicit none

      ! Arguments
      real(real32), intent(in)   :: x(:, :) !< matrix
      type(DescStatMatrixReal32) :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) :: i, j, k, l, n, p, MinNP
      real(real32), allocatable, dimension(:) :: Diag, OffDiag

      n = size(x, 1)
      p = size(x, 2)

      MinNP = minval([n, p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k = 0
      l = 0
      do j = 1, p
        do i = 1, n
          if (i .eq. j) then
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
    !> @brief  Descriptive statistics of a matrix - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function DescStatMatrixR64(x) result(Res)
      implicit none

      ! Arguments
      real(real64), intent(in)   :: x(:, :) !< matrix
      type(DescStatMatrixReal64) :: Res     !< @return descriptive statistics

      ! Other
      integer(int32) :: i, j, k, l, n, p, MinNP
      real(real64), allocatable, dimension(:) :: Diag, OffDiag

      n = size(x, 1)
      p = size(x, 2)

      MinNP = minval([n, p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k = 0
      l = 0
      do j = 1, p
        do i = 1, n
          if (i .eq. j) then
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
    !> @brief   Covariance between two variables - int8
    !> @details x and y are cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function CovI8(x, y, MuX, MuY, w, wType) result(Res)
      implicit none
      ! Arguments
      integer(int8), intent(in)              :: x(:)  !< values for x
      integer(int8), intent(in)              :: y(:)  !< values for y
      real(real32), intent(in), optional     :: MuX   !< precomupted mean for x
      real(real32), intent(in), optional     :: MuY   !< precomupted mean for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return covariance

      ! Other
      real(real32) :: xR32(size(x)), yR32(size(y))

      xR32 = real(x, kind=real32)
      yR32 = real(y, kind=real32)
      Res = CovR32(x=xR32, y=yR32, MuX=MuX, MuY=MuY, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Covariance between two variables - int32
    !> @details x and y are cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function CovI32(x, y, MuX, MuY, w, wType) result(Res)
      implicit none
      ! Arguments
      integer(int32), intent(in)             :: x(:)  !< values for x
      integer(int32), intent(in)             :: y(:)  !< values for y
      real(real32), intent(in), optional     :: MuX   !< precomupted mean for x
      real(real32), intent(in), optional     :: MuY   !< precomupted mean for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      real(real32)                           :: Res   !< @return covariance

      ! Other
      real(real32) :: xR32(size(x)), yR32(size(y))

      xR32 = real(x, kind=real32)
      yR32 = real(y, kind=real32)
      Res = CovR32(x=xR32, y=yR32, MuX=MuX, MuY=MuY, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Covariance between two variables - real32
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function CovR32(x, y, MuX, MuY, w, wType) result(Res)
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

      if (.not. present(MuX)) then
        MuXIn = Mean(x, w)
      else
        MuXIn = MuX
      end if
      if (.not. present(MuY)) then
        MuYIn = Mean(y, w)
      else
        MuYIn = MuY
      end if

      Res = 0.0
      if (.not. present(w)) then
        do i = 1, n
          Res = Res + (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        Res = Res / (n - 1)
      else
        if (any(w .lt. 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not. present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) .ne. "count" .and. trim(wType) .ne. "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Res = Res + w(i) * (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        if (trim(wtype) .eq. "count") then
          Res = Res / (SumW - 1.0)
        end if
        if (trim(wtype) .eq. "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Covariance between two variables - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function CovR64(x, y, MuX, MuY, w, wType) result(Res)
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

      if (.not. present(MuX)) then
        MuXIn = Mean(x, w)
      else
        MuXIn = MuX
      end if
      if (.not. present(MuY)) then
        MuYIn = Mean(y, w)
      else
        MuYIn = MuY
      end if

      Res = 0.0d0
      if (.not. present(w)) then
        do i = 1, n
          Res = Res + (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        Res = Res / (n - 1)
      else
        if (any(w .lt. 0)) then
          write(STDERR, "(a)") "ERROR: Weights should not be negative"
          write(STDERR, "(a)") " "
          stop 1
        end if
        SumW = sum(w)
        if (SumW .lt. EPSILOND) then
          write(STDERR, "(a)") "ERROR: Sum of weights is smaller than EPSILON"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (.not. present(wType)) then
          write(STDERR, "(a)") "ERROR: When weights are given, you must specify the type of weights (wType)"
          write(STDERR, "(a)") " "
          stop 1
        end if
        if (trim(wType) .ne. "count" .and. trim(wType) .ne. "freq") then
          write(STDERR, "(a)") "ERROR: wType must be either count or freq"
          write(STDERR, "(a)") " "
          stop 1
        end if
        ! https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
        do i = 1, n
          Res = Res + w(i) * (x(i) - MuXIn) * (y(i) - MuYIn)
        end do
        if (trim(wtype) .eq. "count") then
          Res = Res / (SumW - 1.0d0)
        end if
        if (trim(wtype) .eq. "freq") then
          Res = Res / (SumW - (sum(w * w) / SumW))
        end if
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Correlation between two variables - int8
    !> @details x and y are cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function CorI8(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      integer(int8), intent(in)              :: x(:)  !< values for x
      integer(int8), intent(in)              :: y(:)  !< values for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationReal32)                :: Res   !< @return correlation

      ! Other
      real(real32) :: xR32(size(x)), yR32(size(y))

      xR32 = real(x, kind=real32)
      yR32 = real(y, kind=real32)
      Res = CorR32(x=xR32, y=yR32, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Correlation between two variables - int32
    !> @details x and y are cast to real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    January 18, 2017
    !---------------------------------------------------------------------------
    function CorI32(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      integer(int32), intent(in)             :: x(:)  !< values for x
      integer(int32), intent(in)             :: y(:)  !< values for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationReal32)                :: Res   !< @return correlation

      ! Other
      real(real32) :: xR32(size(x)), yR32(size(y))

      xR32 = real(x, kind=real32)
      yR32 = real(y, kind=real32)
      Res = CorR32(x=xR32, y=yR32, w=w, wType=wType)
      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Correlation between two variables - real32
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 22, 2016
    !---------------------------------------------------------------------------
    function CorR32(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real32), intent(in)               :: x(:)  !< values for x
      real(real32), intent(in)               :: y(:)  !< values for y
      real(real32), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationReal32)                :: Res   !< @return correlation

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
    !> @brief  Correlation between two variables - real64
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 22, 2016
    !---------------------------------------------------------------------------
    function CorR64(x, y, w, wType) result(Res)
      implicit none
      ! Arguments
      real(real64), intent(in)               :: x(:)  !< values for x
      real(real64), intent(in)               :: y(:)  !< values for y
      real(real64), intent(in), optional     :: w(:)  !< weights
      character(len=*), intent(in), optional :: wType !< type of weights (count, freq)
      type(CorrelationReal64)                :: Res   !< @return correlation

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
