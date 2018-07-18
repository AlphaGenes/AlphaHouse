!###############################################################################
!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaEvolveModule.f90
!
! DESCRIPTION:
!> @brief    Evolutionary algorithms
!
!> @details  Evolutionary algorithms such as random search, differential evolution,
!!           genetic algorithm (not implemented), etc.
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     September 26, 2016
!
!> @version  0.0.2 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 GGorjanc - Initial Version
!
!-------------------------------------------------------------------------------
module AlphaEvolveModule

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use, intrinsic :: IEEE_Arithmetic
  use IntelRNGMod
  use AlphaHouseMod, only : Int2Char, Real2Char, ToLower

  implicit none

  private
  ! Types
  public :: AlphaEvolveSol, AlphaEvolveSpec, AlphaEvolveData
  ! Methods
  public :: DifferentialEvolution, RandomSearch

  !> @brief An evolutionary solution
  type :: AlphaEvolveSol
    real(real64)              :: Objective
    integer(int32)            :: nParam
    real(real64), allocatable :: Chrom(:)
    contains
      procedure :: Initialise => InitialiseAlphaEvolveSol
      procedure :: Assign     => AssignAlphaEvolveSol
      procedure :: UpdateMean => UpdateMeanAlphaEvolveSol
      procedure :: Evaluate   => EvaluateAlphaEvolveSol
      procedure :: Write      => WriteAlphaEvolveSol
      procedure :: Log        => LogAlphaEvolveSol
      procedure :: LogPop     => LogPopAlphaEvolveSol
  end type

  !> @brief Specifications passed to evolutionary algorithm
  type, abstract :: AlphaEvolveSpec
    contains
      procedure :: LogHead    => LogHeadAlphaEvolveSpec
      procedure :: LogPopHead => LogPopHeadAlphaEvolveSpec
  end type

  !> @brief Data passed to evolutionary algorithm
  type, abstract :: AlphaEvolveData
  end type

  contains

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Differential evolution
    !> @details Differential evolution algorithm (Storn and Price) plus additions
    !!          by Brian Kinghorn (vary parameters and genetic algorithm steps).
    !!          This works with continuous representation of a solution.
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  The best evolved solution (BestSol); log on STDOUT and files
    !---------------------------------------------------------------------------
    subroutine DifferentialEvolution(Spec, Data, nParam, nSol, Init, &
                                      nIter, nIterBurnIn, nIterStop, StopTolerance, nIterPrint, &
                                      LogFile, LogStdout, LogPop, LogPopFile, &
                                      CRBurnIn, CRLate1, CRLate2, FBase, FHigh1, FHigh2, Seed, BestSol) ! not pure due to IO & RNG
      
      use IntelRNGMod
      use mkl_vsl_type
      use mkl_vsl
      USE omp_lib
      implicit none

      ! Arguments
      class(AlphaEvolveSpec), intent(in)     :: Spec          !< AlphaEvolveSpec holder
      class(AlphaEvolveData), intent(in)     :: Data          !< AlphaEvolveData holder
      integer(int32), intent(in)             :: nParam        !< No. of parameters in a solution
      integer(int32), intent(in)             :: nSol          !< No. of solutions to test each generation/iteration
      real(real64), intent(in), optional     :: Init(:, :)    !< Initial solutions to start with
      integer(int32), intent(in)             :: nIter         !< No. of generations/iterations to run
      integer(int32), intent(in)             :: nIterBurnIn   !< No. of generations/iterations with more loose parameters
      integer(int32), intent(in)             :: nIterStop     !< Stop after no progress for nIterStop
      real(real64), intent(in)               :: StopTolerance !< Stopping tolerance
      integer(int32), intent(in)             :: nIterPrint    !< Print changed solution every nIterPrint
      character(len=*), intent(in), optional :: LogFile       !< Which file to log best solution into
      logical, intent(in), optional          :: LogStdout     !< Log to STDOUT? (default .true.)
      logical, intent(in), optional          :: LogPop        !< Save all evaluated solutions to a file
      character(len=*), intent(in), optional :: LogPopFile    !< File for the evaluated solutions
      real(real64), intent(in), optional     :: CRBurnIn      !< Crossover rate for nIterBurnIn
      real(real64), intent(in), optional     :: CRLate1       !< Crossover rate (for common small moves)
      real(real64), intent(in), optional     :: CRLate2       !< Crossover rate (for rare large moves)
      real(real64), intent(in), optional     :: FBase         !< F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh1        !< F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh2        !< F is multiplier of difference used to mutate
      integer(int32), intent(in), optional   :: Seed          !< RNG seed value
      class(AlphaEvolveSol), intent(inout)   :: BestSol       !< The best evolved solution

      ! Other
      integer(int32) :: nInit, Param, ParamLoc, Iter, LastIterPrint, LogUnit, LogPopUnit
      integer(int32) :: Sol, a, b, c, RanNumLoc

      real(real32) :: AcceptPct
      real(real64) :: FInt, FBaseInt, FHigh1Int, FHigh2Int, CRInt, CRBurnInInt, CRLateInt1, CRLateInt2
      real(real64) :: Chrom(nParam), OldBestSolObjective
      real(real64), allocatable, dimension(:) :: RanNum

      logical :: DiffOnly, BestSolChanged, LogPopInternal, LogStdoutInternal

      class(AlphaEvolveSol), allocatable :: OldSol(:), NewSol(:)
      type(vsl_stream_state) :: Stream
      integer :: nthreads, ithread,  start, end
      ! --- Trap errors ---

      if (nSol < 4) then
        write(STDERR, "(a)") " ERROR nSol must be at least 4!"
        write(STDERR, "(a)") " "
        stop 1
      end if

      ! --- Allocate and initialise ---

      allocate(OldSol(nSol), source=BestSol)
      allocate(NewSol(nSol), source=BestSol)

      LastIterPrint = 0
      BestSol%Objective = -huge(BestSol%Objective)
      OldBestSolObjective = BestSol%Objective

      

      ! --- Logging ---

      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .true.
      end if

      if (present(LogPop)) then
        LogPopInternal = LogPop
      else
        LogPopInternal = .false.
      end if

      if (LogPopInternal) then
        if (.not. present(LogPopFile)) then
          write (STDERR, "(a)") "ERROR: When LogPop is .true. the LogPopFile must be given!"
          write (STDERR, "(a)") " "
          stop 1
        end if
      end if

      if (present(LogFile)) then
        open(newunit=LogUnit, file=trim(LogFile), status="unknown")
        call Spec%LogHead(LogUnit=LogUnit)
      end if
      if (LogStdoutInternal) then
        call Spec%LogHead
      end if
      if (LogPopInternal) then
        open(newunit=LogPopUnit, file=trim(LogPopFile), status="unknown")
        call Spec%LogPopHead(LogPopUnit=LogPopUnit)
      end if

      ! --- Set parameters ---

      ! Crossover rate
      ! Between 0 and 1, good values are 0.1 (common slow moves) and 0.9 (rare large moves)
      ! ... for first few generations (burn-in)
      if (present(CRBurnIn)) then
        CRBurnInInt = CRBurnIn
      else
        CRBurnInInt = 0.9d0
      end if
      ! ... for later climbs (common slow moves)
      if (present(CRLate1)) then
        CRLateInt1 = CRLate1
      else
        CRLateInt1 = 0.1d0
      end if
      ! ...                  (rare large moves)
      if (present(CRLate2)) then
        CRLateInt2 = CRLate2
      else
        CRLateInt2 = 0.9d0
      end if

      ! F is multiplier of difference used to mutate
      ! Typically between 0.2 and 1.2 (even up to 2.0)
      ! (if alleles should be integer, keep F as integer)
      ! ... conservative moves
      if (present(FBase)) then
        FBaseInt = FBase
      else
        FBaseInt = 0.1d0
      end if
      ! ... adventurous moves
      if (present(FHigh1)) then
        FHigh1Int = FHigh1
      else
        FHigh1Int = 10.0d0 * FBaseInt
      end if
      if (present(FHigh2)) then
        FHigh2Int = FHigh2
      else
        FHigh2Int = 4.0d0 * FHigh1Int
      end if

      ! --- Initialise founder population of solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do Sol = 1, nInit
          call OldSol(Sol)%Evaluate(Chrom=Init(:, Sol), Spec=Spec, Data=Data)
        end do
        nInit = Sol
      else
        nInit = 1
      end if
      do Sol = nInit, nSol
        call OldSol(Sol)%Evaluate(Chrom=SampleIntelUniformD(n=nParam), Spec=Spec, Data=Data)
      end do

      Sol = maxloc(OldSol(:)%Objective, dim=1)
      call BestSol%Assign(OldSol(Sol))
      if (present(LogFile)) then
        call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=0, AcceptPct=IEEE_Value(x=AcceptPct, class=IEEE_Quiet_NaN))
      end if
      if (LogStdoutInternal) then
        call BestSol%Log(Spec=Spec, Iteration=0, AcceptPct=IEEE_Value(x=AcceptPct, class=IEEE_Quiet_NaN))
      end if
      if (LogPopInternal) then
        do Sol = 1, nSol
          call OldSol(Sol)%LogPop(Spec=Spec, LogPopUnit=LogPopUnit, Iteration=0, i=Sol)
        end do
      end if

      ! --- Evolve ---

      do Iter = 1, nIter

        ! Vary differential and non-differential mutation to escape valleys
        if (mod(Iter, 3) == 0) then
          DiffOnly = .true.
        else
          DiffOnly = .false.
        end if

        ! Burn-in
        if (Iter < nIterBurnIn) then
          CRInt = CRBurnInInt
        else
          ! Vary crossover rate every few generations
          if (mod(Iter, 5) == 0) then
            CRInt = CRLateInt2
          else
            CRInt = CRLateInt1
          end if
        end if

        ! Vary mutation rate every few generations
        if (mod(Iter, 4) == 0) then
          FInt = FHigh1Int
        else
          FInt = FBaseInt
        end if
        if (mod(Iter, 7) == 0) then
          FInt = FHigh2Int
        else
          FInt = FBaseInt
        end if

        ! --- Generate competitors ---

        BestSolChanged = .false.
        AcceptPct = 0.0
        !$OMP PARALLEL PRIVATE(Sol, RanNum, RanNumLoc, a, b, c, Param, ParamLoc, Chrom, stream, ithread, start, end)
        

        nthreads = omp_get_num_threads()
        ithread = omp_get_thread_num()
        start = (ithread*nSol/nthreads)+1
        end = ((ithread+1)*nSol/nthreads)

        call IntitialiseIntelRNGStream(stream, Seed=Seed)
        

        ! $OMP PARALLEL DO PRIVATE(Sol, RanNum, RanNumLoc, a, b, c, Param, ParamLoc, Chrom)
        do Sol = start, end

          ! --- Mutate and crossover ---

          RanNum = SampleIntelUniformD(n=5*nParam, stream=stream)
          RanNumLoc = 0
          ! Get three different solutions

          ! a = Sol
          ! do while (a == Sol)
          !   RanNumLoc = RanNumLoc + 1
          !   a = int(RanNum(RanNumLoc) * nSol) + 1
          ! end do
          ! b = Sol
          ! do while ((b == Sol) .or. (b == a))
          !   RanNumLoc = RanNumLoc + 1
          !   b = int(RanNum(RanNumLoc) * nSol) + 1
          ! end do
          ! c = Sol
          ! do while ((c == Sol) .or. (c == a) .or. (c == b))
          !   RanNumLoc = RanNumLoc + 1
          !   c = int(RanNum(RanNumLoc) * nSol) + 1
          ! end do

          RanNumLoc = RanNumLoc + 1
          a = int(RanNum(RanNumLoc) * nSol) + 1
          RanNumLoc = RanNumLoc + 1
          b = int(RanNum(RanNumLoc) * nSol) + 1
          do while (b == a)
            RanNumLoc = RanNumLoc + 1
            b = int(RanNum(RanNumLoc) * nSol) + 1
          end do
          RanNumLoc = RanNumLoc + 1
          c = int(RanNum(RanNumLoc) * nSol) + 1
          do while ((c == a) .or. (c == b))
            RanNumLoc = RanNumLoc + 1
            c = int(RanNum(RanNumLoc) * nSol) + 1
          end do

          ! Mate the solutions to get a new competitor solution
          RanNumLoc = RanNumLoc + 1
          Param = int(RanNum(RanNumLoc) * nParam) + 1 ! Cycle through parameters starting at a random point
          do ParamLoc = 1, nParam
            RanNumLoc = RanNumLoc + 1
            if ((RanNum(RanNumLoc) < CRInt) .or. (ParamLoc == nParam)) then
              ! Crossover
              RanNumLoc = RanNumLoc + 1
              if ((RanNum(RanNumLoc) < 0.8d0) .or. DiffOnly) then
                ! Differential mutation (with prob 0.8 or 1)
                Chrom(Param) = OldSol(c)%Chrom(Param) + FInt * (OldSol(a)%Chrom(Param) - OldSol(b)%Chrom(Param))
              else
                ! Non-differential mutation (to avoid getting stuck)
                RanNumLoc = RanNumLoc + 1
                if (RanNum(RanNumLoc) < 0.5d0) then
                  RanNumLoc = RanNumLoc +1
                  Chrom(Param) = OldSol(c)%Chrom(Param) * (0.9d0 + 0.2d0 * RanNum(RanNumLoc))
                else
                  RanNumLoc = RanNumLoc + 1
                  Chrom(Param) = OldSol(c)%Chrom(Param) + 0.01d0 * FInt * (OldSol(a)%Chrom(Param) + 0.01d0) * (RanNum(RanNumLoc) - 0.5d0)
                end if
              end if
            else
              ! Do not crossover
              Chrom(Param) = OldSol(Sol)%Chrom(Param)
            end if
            Param = Param + 1
            if (Param > nParam) then
              Param = Param - nParam
            end if
          end do ! ParamLoc

          ! --- Evaluate and Select ---

          ! Merit of competitor
          call NewSol(Sol)%Evaluate(Chrom=Chrom, Spec=Spec, Data=Data, Random=RanNum(RanNumLoc))
          RanNumLoc = RanNumLoc + 1
            ! If competitor is better or equal, keep it ("equal" to force evolution)
          if (NewSol(Sol)%Objective >= OldSol(Sol)%Objective) then
            AcceptPct = AcceptPct + 1.0
          else
            ! Keep the old solution
            call NewSol(Sol)%Assign(OldSol(Sol))
          end if
        end do ! Sol
    
        call UnintitialiseIntelRNGStream(stream)
        !$OMP END PARALLEL

        AcceptPct = AcceptPct / nSol * 100.0

        ! --- Promote the new solutions into parental solutions ---

        do Sol = 1, nSol
          call OldSol(Sol)%Assign(NewSol(Sol))
        end do

        ! --- The current best solution ---

        Sol = maxloc(NewSol(:)%Objective, dim=1)
        if (NewSol(Sol)%Objective > BestSol%Objective) then
          call BestSol%Assign(NewSol(Sol))
          BestSolChanged = .true.
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Iter == 1) .or. ((Iter - LastIterPrint) >= nIterPrint)) then
            LastIterPrint = Iter
            if (present(LogFile)) then
              call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=Iter, AcceptPct=AcceptPct)
            end if
            if (LogStdoutInternal) then
              call BestSol%Log(Spec=Spec, Iteration=Iter, AcceptPct=AcceptPct)
            end if
            if (LogPopInternal) then
              do Sol = 1, nSol
                call NewSol(Sol)%LogPop(Spec=Spec, LogPopUnit=LogPopUnit, Iteration=Iter, i=Sol)
              end do
            end if
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Iter, nIterStop) == 0) then
          if ((BestSol%Objective - OldBestSolObjective) > StopTolerance) then
            OldBestSolObjective = BestSol%Objective
          else
            if (LogStdoutInternal) then
              write(STDOUT, "(5a)") "NOTE: Objective did not improve for ", &
                trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nIterStop)), &
                " iterations. Stopping the optimisation."
              write(STDOUT, "(a)") " "
            end if
            exit
          end if
        end if

      end do ! Iter

      ! --- The winner solution ---

      if (present(LogFile)) then
        call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=Iter, AcceptPct=AcceptPct)
        close(LogUnit)
      end if
      if (LogStdoutInternal) then
        call BestSol%Log(Spec=Spec, Iteration=Iter, AcceptPct=AcceptPct)
      end if

      if (LogPopInternal) then
        do Sol = 1, nSol
          call NewSol(Sol)%LogPop(Spec=Spec, LogPopUnit=LogPopUnit, Iteration=Iter, i=Sol)
        end do
        close(LogPopUnit)
      end if

      
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Random search
    !> @details Can either find the best solution (Mode=max) or return mean of
    !!          the evaluated solutions (Mode=avg)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  The best found solution or mean of all evaluated solutions (BestSol);
    !!          log on STDOUT and files
    !---------------------------------------------------------------------------
    subroutine RandomSearch(Mode, Spec, Data, nParam, Init, &
                            nSamp, nSampStop, StopTolerance, &
                            nSampPrint, LogFile, LogStdout, BestSol) ! not pure due to IO & RNG
      implicit none

      ! Arguments
      character(len=*), intent(in)           :: Mode           !< Mode of search (max or avg)
      class(AlphaEvolveSpec), intent(in)     :: Spec           !< AlphaEvolveSpec holder
      class(AlphaEvolveData), intent(in)     :: Data           !< AlphaEvolveData holder
      integer(int32), intent(in)             :: nParam         !< No. of parameters in a solution
      real(real64), intent(inout), optional  :: Init(:, :)     !< Initial solutions to test
      integer(int32), intent(in)             :: nSamp          !< No. of samples to test
      integer(int32), intent(in), optional   :: nSampStop      !< Stop after no progress for nSampStop
      real(real64), intent(in)               :: StopTolerance  !< Stopping tolerance
      integer(int32), intent(in)             :: nSampPrint     !< Print changed solution every nSampPrint
      character(len=*), intent(in), optional :: LogFile        !< Which file to log best solution into
      logical, intent(in), optional          :: LogStdout      !< Log to STDOUT?
      class(AlphaEvolveSol), intent(inout)   :: BestSol        !< The best found solution

      ! Other
      integer(int32) :: nInit, Samp, LastSampPrint, LogUnit

      real(real32) :: AcceptPct
      real(real64) :: RanNum, OldBestSolObjective, BestSolObjective

      logical :: ModeAvg, ModeMax, BestSolChanged, LogStdoutInternal

      class(AlphaEvolveSol), allocatable :: TestSol

      if (present(LogStdout)) then
        LogStdoutInternal = LogStdout
      else
        LogStdoutInternal = .true.
      end if

      ! --- Mode ---

      ModeAvg = .false.
      ModeMax = .false.

      if      (ToLower(trim(Mode)) == "avg") then
        ModeAvg = .true.
      else if (ToLower(trim(Mode)) == "max") then
        ModeMax = .true.
      else
        write (STDERR, "(a)") "ERROR: Mode must be either avg or max!"
        write (STDERR, "(a)") " "
        stop 1
      end if

      ! --- Allocate and Initialise ---

      allocate(TestSol, source=BestSol)

      LastSampPrint = 0

      if (ModeAvg) then
        OldBestSolObjective = 0.0d0
        BestSolObjective = 0.0d0
        BestSolChanged = .true.
        AcceptPct = 100.0
      end if

      if (ModeMax) then
        OldBestSolObjective = -huge(RanNum)
        BestSolObjective = -huge(RanNum)
        BestSolChanged = .false.
        AcceptPct = 0.0
      end if

      call IntitialiseIntelRNG()

      ! --- Printout log header ---

      if (present(LogFile)) then
        open(newunit=LogUnit, file=trim(LogFile), status="unknown")
        call Spec%LogHead(LogUnit=LogUnit)
      end if
      if (LogStdoutInternal) then
        call Spec%LogHead
      end if

      ! --- Initialise with the provided solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do Samp = 1, nInit
          call TestSol%Evaluate(Chrom=Init(:, Samp), Spec=Spec, Data=Data)
          if      (ModeAvg) then
            if (Samp == 1) then
              call BestSol%Assign(TestSol)
            else
              call BestSol%UpdateMean(TestSol, Samp)
            end if
          else if (ModeMax) then
            if (TestSol%Objective > BestSolObjective) then
              call BestSol%Assign(TestSol)
              BestSolObjective = BestSol%Objective
              BestSolChanged = .true.
              AcceptPct = AcceptPct + 1.0
            end if
          end if
        end do
      else
        nInit = 1
      end if

      ! --- Search ---

      ! @todo must make outer and inner loop as in DifferentialEvolution, because
      !       the exit statement must not be within the OMP look
      ! $OMP PARALLEL DO PRIVATE(Samp)
      do Samp = nInit, nSamp

        BestSolChanged = .false.

        ! --- Evaluate a competitor and Select ---

        ! Merit of a competitor
        call TestSol%Evaluate(Chrom=SampleIntelUniformD(n=nParam), Spec=Spec, Data=Data)

        if      (ModeAvg) then
          ! Update the mean
          if (Samp == 1) then
            call BestSol%Assign(TestSol)
          else
            call BestSol%UpdateMean(TestSol, Samp)
          end if
          BestSolChanged = .true.
        else if (ModeMax) then
          ! If the competitor is better, keep it
          if (TestSol%Objective > BestSolObjective) then
            call BestSol%Assign(TestSol)
            BestSolObjective = BestSol%Objective
            BestSolChanged = .true.
            AcceptPct = AcceptPct + 1.0
          end if
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Samp == 1) .or. ((Samp - LastSampPrint) >= nSampPrint)) then
            if      (ModeAvg) then
              call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=Samp, AcceptPct=AcceptPct)
            else if (ModeMax) then
              AcceptPct = AcceptPct / (Samp - LastSampPrint) * 100.0
              if (present(LogFile)) then
                call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=Samp, AcceptPct=AcceptPct)
              end if
              if (LogStdoutInternal) then
                call BestSol%Log(Spec=Spec, Iteration=Samp, AcceptPct=AcceptPct)
              end if
              AcceptPct = 0.0
            end if
            LastSampPrint = Samp
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Samp, nSampStop) == 0) then
          if ((BestSol%Objective - OldBestSolObjective) > StopTolerance) then
            OldBestSolObjective = BestSol%Objective
          else
            if (LogStdoutInternal) then
              write(STDOUT, "(5a)") "NOTE: Objective did not improve for ", &
                trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nSampStop)), &
                " samples. Stopping the random search."
              write(STDOUT, "(a)") " "
            end if
            exit
          end if
        end if

      end do ! Samp
      ! $OMP END PARALLEL DO

      ! --- The winner solution ---

      if (ModeMax) then
        AcceptPct = AcceptPct / (Samp - LastSampPrint) * 100.0
      end if
      if (present(LogFile)) then
        call BestSol%Log(Spec=Spec, LogUnit=LogUnit, Iteration=Samp, AcceptPct=AcceptPct)
        close(LogUnit)
      end if
      if (LogStdoutInternal) then
        call BestSol%Log(Spec=Spec, Iteration=Samp, AcceptPct=AcceptPct)
      end if

      call UnintitialiseIntelRNG
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Initialise AlphaEvolve solution
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure subroutine InitialiseAlphaEvolveSol(This, Chrom, Spec) ! Spec used in future methods
      implicit none

      ! Argument
      class(AlphaEvolveSol), intent(out)           :: This     !< AlphaEvolveSol holder
      real(real64), intent(in)                     :: Chrom(:) !< A solution
      class(AlphaEvolveSpec), intent(in), optional :: Spec     !< Specifications (used in future methods)

      ! Initialisation
      This%Objective = -huge(This%Objective)
      This%nParam = size(Chrom)
      allocate(This%Chrom(This%nParam))
      This%Chrom = Chrom
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Assign one AlphaEvolve solution to another
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    pure subroutine AssignAlphaEvolveSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(out) :: Out !< @return output solution
      class(AlphaEvolveSol), intent(in)  :: In  !< input solution

      ! Assignments
      Out%Objective = In%Objective
      Out%nParam = In%nParam
      Out%Chrom = In%Chrom
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Update (running) mean of AlphaEvolve solutions for random search
    !!          with Mode=avg
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Solution with average metrics
    !---------------------------------------------------------------------------
    pure subroutine UpdateMeanAlphaEvolveSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(inout) :: This !< solution
      class(AlphaEvolveSol), intent(in)    :: Add  !< addition
      integer(int32), intent(in)           :: n    !< number of previous solutions mean was based on

      ! Other
      real(real64) :: kR

      ! Updates
      kR = (dble(n) - 1.0d0) / n

      This%Objective = This%Objective * kR + Add%Objective / n
      ! This%nParam  = This%nParam    * kR + Add%nParam    / n ! the same all the time
      ! This%Chrom   = This%Chrom     * kR + Add%Chrom     / n ! hmm, do we really want to average over chromosomes?
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Evaluate AlphaEvolve solution
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   September 26, 2016
    !---------------------------------------------------------------------------
    subroutine EvaluateAlphaEvolveSol(This, Chrom, Spec, Data, Random) ! Data used in future methods; Not pure as some future Evaluate methods might use RNG
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(inout)         :: This     !< @return solution
      real(real64), intent(in)                     :: Chrom(:) !< A solution
      class(AlphaEvolveSpec), intent(in)           :: Spec     !< AlphaEvolveSpec holder
      class(AlphaEvolveData), intent(in), optional :: Data     !< AlphaEvolveData holder (optional for future methods)
      real(real64), intent(in), optional           :: Random !< random number to avoid reinitialistation
      call This%Initialise(Chrom=Chrom, Spec=Spec)
      This%Objective = sum(This%Chrom) ! just a sum here for simplicity
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief  Write AlphaEvolveSol to a file or standard output
    !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date   March 25, 2017
    !> @return Output to a file or standard output
    !---------------------------------------------------------------------------
    subroutine WriteAlphaEvolveSol(This, File) ! not pure due to IO
      implicit none
      class(AlphaEvolveSol), intent(in)      :: This !< AlphaEvolveSol holder
      character(len=*), intent(in), optional :: File !< Filename, if missing use standard output

      integer(int32) :: Unit
      if (present(File)) then
        open(newunit=Unit, file=File, action="write", status="unknown")
      else
        Unit = STDOUT
      end if

      write(Unit, *) "Objective: ", This%Objective
      write(Unit, *) "nParam: ", This%nParam
      write(Unit, *) "Chrom: ", This%Chrom

      if (present(File)) then
        close(Unit)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print log head
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print log head to unit
    !---------------------------------------------------------------------------
    subroutine LogHeadAlphaEvolveSpec(This, LogUnit, String, StringNum) ! This used in future methods; not pure due to IO
      implicit none
      class(AlphaEvolveSpec), intent(in)     :: This      !< Spec holder
      integer(int32), intent(in), optional   :: LogUnit   !< log file unit (default STDOUT)
      character(len=*), intent(in), optional :: String    !< additional string that will be written before the head
      integer(int32), optional               :: StringNum !< How much space is needed for the String

      character(len=12) :: ColnameLogStdout(3), StringFmt
      character(len=22) :: ColnameLogUnit(3)
      !                      123456789012
      ColnameLogStdout(1) = "   Iteration"
      ColnameLogStdout(2) = "   AcceptPct"
      ColnameLogStdout(3) = "   Objective"
      !                    1234567890123456789012
      ColnameLogUnit(1) = "             Iteration"
      ColnameLogUnit(2) = "             AcceptPct"
      ColnameLogUnit(3) = "             Objective"
      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "(a"//Int2Char(StringNum)//")"
        else
          StringFmt = "(a)"
        end if
      end if
      if (present(LogUnit)) then
        if (present(String)) then
          write(LogUnit, StringFmt, Advance="No") String
        end if
        write(LogUnit, "(3a22)") ColnameLogUnit(:)
      else
        if (present(String)) then
          write(STDOUT, StringFmt, Advance="No") String
        end if
        write(STDOUT,  "(3a12)") ColnameLogStdout(:)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print log
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print log to unit
    !---------------------------------------------------------------------------
    subroutine LogAlphaEvolveSol(This, Spec, LogUnit, Iteration, AcceptPct, String, StringNum) ! Spec used in future methods; not pure due to IO
      implicit none
      class(AlphaEvolveSol), intent(in)      :: This      !< solution
      class(AlphaEvolveSpec), intent(in)     :: Spec      !< spec holder
      integer(int32), intent(in), optional   :: LogUnit   !< log file unit (default STDOUT)
      integer(int32), intent(in)             :: Iteration !< generation/iteration
      real(real32), intent(in)               :: AcceptPct !< acceptance rate
      character(len=*), intent(in), optional :: String    !< additional string that will be written before the head
      integer(int32), optional               :: StringNum !< How much space is needed for the String

      integer(int32) :: Unit
      character(len=20) :: Fmt, StringFmt

      if (present(LogUnit)) then
        Unit = LogUnit
        Fmt = "(a22, 2(1x, es21.14))"
      else
        Unit = STDOUT
        Fmt = "(a12, 1x, f11.1, 1x, f11.5))"
      end if
      if (present(String)) then
        if (present(StringNum)) then
          StringFmt = "(a"//Int2Char(StringNum)//")"
        else
          StringFmt = "(a)"
        end if
      end if
      if (present(String)) then
        write(Unit, StringFmt, Advance="No") String
      end if
      write(Unit, Fmt) Iteration, AcceptPct, This%Objective
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print population log head
    !> @details This is meant to log all the evaluated solutions (the population)
    !!          and not just the best one as LogHeadAlphaEvolveSpec does
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print population log head to unit
    !---------------------------------------------------------------------------
    subroutine LogPopHeadAlphaEvolveSpec(This, LogPopUnit) ! This used in future methods; not pure due to IO
      implicit none
      class(AlphaEvolveSpec), intent(in)   :: This       !< spec
      integer(int32), intent(in), optional :: LogPopUnit !< log file unit (default STDOUT)
      integer(int32) :: Unit
      character(len=22) :: ColnameLogPopUnit(3)
      if (present(LogPopUnit)) then
        Unit = LogPopUnit
      else
        Unit = STDOUT
      end if
      !                       1234567890123456789012
      ColnameLogPopUnit(1) = "             Iteration"
      ColnameLogPopUnit(2) = "              Solution"
      ColnameLogPopUnit(3) = "             Objective"
      write(Unit, "(3a22)") ColnameLogPopUnit(:)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print population log
    !> @details This is meant to log all the evaluated solutions (the population)
    !!          and not just the best one as LogAlphaEvolveSol does
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print population log to unit
    !---------------------------------------------------------------------------
    subroutine LogPopAlphaEvolveSol(This, Spec, LogPopUnit, Iteration, i) ! Spec used in future methods; not pure due to IO
      implicit none
      class(AlphaEvolveSol), intent(in)    :: This       !< solution
      class(AlphaEvolveSpec), intent(in)   :: Spec       !< spec
      integer(int32), intent(in), optional :: LogPopUnit !< population log file unit (default STDOUT)
      integer(int32), intent(in)           :: Iteration  !< generation/iteration
      integer(int32), intent(in)           :: i          !< solution id
      integer(int32) :: Unit
      if (present(LogPopUnit)) then
        Unit = LogPopUnit
      else
        Unit = STDOUT
      end if
      write(Unit, "(2(i22, 1x), es21.14)") Iteration, i, This%Objective
    end subroutine

    !###########################################################################
end module

!###############################################################################
