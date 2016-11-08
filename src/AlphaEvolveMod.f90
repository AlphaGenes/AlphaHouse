
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaEvolveMod.f90
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
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-09-26 GGorjanc - Initial Version
!
!-------------------------------------------------------------------------------
module AlphaEvolveMod

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use AlphaHouseMod, only : int2Char, Real2Char, ToLower

  implicit none

  private
  ! Types
  public :: AlphaEvolveSol
  ! Methods
  public :: DifferentialEvolution, RandomSearch

  !> @brief An evolutionary solution
  type :: AlphaEvolveSol
    real(real64) :: Criterion
    contains
      procedure         :: Initialise => InitialiseAlphaEvolveSol
      procedure         :: Assign     => AssignAlphaEvolveSol
      procedure         :: UpdateMean => UpdateMeanAlphaEvolveSol
      procedure         :: Evaluate   => EvaluateAlphaEvolveSol
      procedure, nopass :: LogHead    => LogHeadAlphaEvolveSol
      procedure         :: Log        => LogAlphaEvolveSol
      procedure, nopass :: LogPopHead => LogPopHeadAlphaEvolveSol
      procedure         :: LogPop     => LogPopAlphaEvolveSol
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
    subroutine DifferentialEvolution(nParam, nSol, Init, nGen, nGenBurnIn, nGenStop,&
      StopTolerance, nGenPrint, LogFile, LogPop, LogPopFile, CritType, CRBurnIn, CRLate, FBase,&
      FHigh1, FHigh2, BestSol)
      implicit none

      ! Arguments
      integer(int32), intent(in)             :: nParam         !< No. of parameters in a solution
      integer(int32), intent(in)             :: nSol           !< No. of solutions to test each generation/Iteration
      real(real64), intent(in), optional     :: Init(:,:)      !< Initial solutions to start with
      integer(int32), intent(in)             :: nGen           !< No. of generations/iterations to run
      integer(int32), intent(in)             :: nGenBurnIn     !< No. of generations/iterations with more loose parameters
      integer(int32), intent(in)             :: nGenStop       !< Stop after no progress for nGenStop
      real(real64), intent(in)               :: StopTolerance  !< Stopping tolerance
      integer(int32), intent(in)             :: nGenPrint      !< Print changed solution every nGenPrint
      character(len=*), intent(in)           :: LogFile        !< Which file to log best solution into
      logical, intent(in), optional          :: LogPop         !< Log all solutions
      character(len=*), intent(in), optional :: LogPopFile     !< Which file to log all solutions into
      character(len=*), intent(in), optional :: CritType       !< Passed to Evaluate
      real(real64), intent(in), optional     :: CRBurnIn       !< Crossover rate for nGenBurnIn
      real(real64), intent(in), optional     :: CRLate         !< Crossover rate
      real(real64), intent(in), optional     :: FBase          !< F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh1         !< F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh2         !< F is multiplier of difference used to mutate
      class(AlphaEvolveSol), intent(out)     :: BestSol        !< The best evolved solution
      !class(*), intent(inout)                :: BestSol        ! The best evolved solution

      ! Other
      integer(int32) :: nInit, Param, ParamLoc, Gen, LastGenPrint, LogUnit, LogPopUnit
      integer(int32) :: i, a, b, c, j
      ! integer(int32) :: OMP_get_num_threads,OMP_get_thread_num

      real(real64) :: RanNum, FInt, FBaseInt, FHigh1Int, FHigh2Int, CRInt, CRBurnInInt, CRLateInt
      real(real64) :: AcceptRate, OldChrom(nParam, nSol), NewChrom(nParam, nSol), Chrom(nParam)
      real(real64) :: BestSolCriterion

      logical :: DiffOnly, BestSolChanged, LogPopInternal

      class(AlphaEvolveSol), allocatable :: Sol(:), HoldSol
      !class(*), allocatable :: Sol(:), HoldSol

      ! --- Allocate and Initialize ---

      allocate(Sol(nSol), source=BestSol)
      allocate(HoldSol, source=BestSol)

      LastGenPrint = 0
      BestSolCriterion = -huge(RanNum)
      BestSol%Criterion = BestSolCriterion

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

      ! --- Printout log header ---

      open(newunit=LogUnit, file=trim(LogFile), status="unknown")
      call HoldSol%LogHead(LogUnit=LogUnit)
      if (LogPopInternal) then
        open(newunit=LogPopUnit, file=trim(LogPopFile), status="unknown")
        call HoldSol%LogPopHead(LogPopUnit)
      end if

      ! --- Set parameters ---

      ! Crossover rate
      ! ... for later climbs
      if (present(CRLate)) then
        CRLateInt = CRLate
      else
        CRLateInt = 0.1d0
      end if
      ! ... for first few generations (burn-in)
      if (present(CRBurnIn)) then
        CRBurnInInt = CRBurnIn
      else
        CRBurnInInt = 2.0d0 * CRLateInt
      end if

      ! F is multiplier of difference used to mutate
      ! Typically between 0.2 and 2.0
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

      ! --- Initialise foundation population of solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do i = 1, nInit
          OldChrom(:,i) = Init(:,i)
          call Sol(i)%Evaluate(OldChrom(:,i), CritType)
        end do
        nInit = i
      else
        nInit = 1
      end if
      do i = nInit, nSol
        call random_number(OldChrom(:,i))
        call Sol(i)%Evaluate(OldChrom(:,i), CritType)
      end do

      ! --- Evolve ---

      do Gen = 1, nGen

        ! Vary differential and non-differential mutation to escape valleys
        if (mod(Gen, 3) == 0) then
          DiffOnly = .true.
        else
          DiffOnly = .false.
        end if

        ! Burn-in
        if (Gen < nGenBurnIn) then
          CRInt = CRBurnInInt
        else
          CRInt = CRLateInt
        end if

        ! Vary mutation rate every few generations
        if (mod(Gen, 4) == 0) then
          FInt = FHigh1Int
        else
          FInt = FBaseInt
        end if

        if (mod(Gen, 7) == 0) then
          FInt = FHigh2Int
        else
          FInt = FBaseInt
        end if

        ! --- Generate competitors ---

        !> @todo: Paralelize this loop?
        !!        The main reason would be to speed up the code as Evaluate() might take quite some time for larger problems
        !!         - some variables are local, say a, b, ...
        !!         - global variable is NewChrom, but is indexed with i so this should not be a problem
        !!         - AcceptRate needs to be in sync between the threads!!!
        !!         - we relly on random_number a lot here and updating the RNG state for each thread can be slow
        !!           and I (GG) am also not sure if we should not have thread specific RNGs
        BestSolChanged = .false.
        AcceptRate = 0.0d0

        ! call OMP_set_num_threads(1)

        ! $OMP PARALLEL DO DEFAULT(PRIVATE)
        do i = 1, nSol

          ! print *, "#Threads: ", OMP_get_num_threads(), "Thread; ", OMP_get_thread_num()+1, ", Solution: ", Sol

          ! --- Mutate and recombine ---

          ! Get three different solutions
          a = i
          do while (a == i)
            call random_number(RanNum)
            a = int(RanNum * nSol) + 1
          end do
          b = i
          do while ((b == i) .or. (b == a))
            call random_number(RanNum)
            b = int(RanNum * nSol) + 1
          end do
          c = i
          do while ((c == i) .or. (c == a) .or. (c == b))
            call random_number(RanNum)
            c = int(RanNum * nSol) + 1
          end do

          ! Mate the solutions
          call random_number(RanNum)
          Param = int(RanNum * nParam) + 1 ! Cycle through parameters starting at a random point
          do ParamLoc = 1, nParam
            call random_number(RanNum)
            if ((RanNum < CRInt) .or. (ParamLoc == nParam)) then
              ! Recombine
              call random_number(RanNum)
              if ((RanNum < 0.8d0) .or. DiffOnly) then
                ! Differential mutation (with prob 0.8 or 1)
                Chrom(Param) = OldChrom(Param, c) + FInt * (OldChrom(Param, a) - OldChrom(Param, b))
              else
                ! Non-differential mutation (to avoid getting stuck)
                call random_number(RanNum)
                if (RanNum < 0.5d0) then
                  call random_number(RanNum)
                  Chrom(Param) = OldChrom(Param, c) * (0.9d0 + 0.2d0 * RanNum)
                else
                  call random_number(RanNum)
                  Chrom(Param) = OldChrom(Param, c) + 0.01d0 * FInt * (OldChrom(Param, a) + 0.01d0) * (RanNum - 0.5d0)
                end if
              end if
            else
              ! Do not recombine
              Chrom(Param) = OldChrom(Param, i)
            end if
            Param = Param + 1
            if (Param > nParam) then
              Param = Param - nParam
            end if
          end do

          ! --- Evaluate and Select ---

          call HoldSol%Evaluate(Chrom, CritType)          ! Merit of competitor
          if (HoldSol%Criterion >= Sol(i)%Criterion) then ! If competitor is better or equal, keep it
            NewChrom(:,i) = Chrom(:)                      !   ("equal" to force evolution)
            call Sol(i)%Assign(HoldSol)
            ! $OMP ATOMIC
            AcceptRate = AcceptRate + 1.0d0
          else
            NewChrom(:,i) = OldChrom(:,i)                 ! Else keep the old solution
          end if
        end do ! i
        ! $OMP END PARALLEL DO

        AcceptRate = AcceptRate / nSol

        ! --- New parents ---

        OldChrom(:,:) = NewChrom(:,:)

        ! --- The current best solution ---

        j = maxloc(Sol(:)%Criterion, dim=1)
        if (Sol(j)%Criterion > BestSol%Criterion) then
          BestSolChanged = .true.
          call BestSol%Assign(Sol(j))
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Gen == 1) .or. ((Gen - LastGenPrint) >= nGenPrint)) then
            LastGenPrint = Gen
            call BestSol%Log(LogUnit, Gen, AcceptRate)
            if (LogPopInternal) then
              do i = 1, nSol
                call Sol(i)%LogPop(LogPopUnit, Gen, i)
              end do
            end if
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Gen, nGenStop) == 0) then
          if ((BestSol%Criterion - BestSolCriterion) > StopTolerance) then
            BestSolCriterion = BestSol%Criterion
          else
            write(STDOUT, "(5a)") "NOTE: Criterion did not improve for ", &
              trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nGenStop)), &
              " iterations. Stopping the optimisation."
            write(STDOUT, "(a)") " "
            exit
          end if
        end if

      end do ! Gen

      ! --- The winner solution ---

      call BestSol%Log(LogUnit, Gen, AcceptRate)
      write(STDOUT, "(a)") " "
      close(LogUnit)
      if (LogPopInternal) then
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
    subroutine RandomSearch(Mode, nParam, Init, nSamp, nSampStop, StopTolerance, &
      nSampPrint, LogFile, CritType, BestSol)
      implicit none

      ! Arguments
      character(len=*), intent(in)           :: Mode           !< Mode of search (max or avg)
      integer(int32), intent(in)             :: nParam         !< No. of parameters in a solution
      real(real64), intent(inout), optional  :: Init(:,:)      !< Initial solutions to test
      integer(int32), intent(in)             :: nSamp          !< No. of samples to test
      integer(int32), intent(in), optional   :: nSampStop      !< Stop after no progress for nSampStop
      real(real64), intent(in)               :: StopTolerance  !< Stopping tolerance
      integer(int32), intent(in)             :: nSampPrint     !< Print changed solution every nSampPrint
      character(len=*), intent(in)           :: LogFile        !< Which file to log best solution into
      character(len=*), intent(in), optional :: CritType       !< Passed to Evaluate
      class(AlphaEvolveSol), intent(out)     :: BestSol        !< The best found solution
      !class(*), intent(inout)                :: BestSol        ! The best found solution

      ! Other
      integer(int32) :: nInit, Samp, LastSampPrint, LogUnit
      ! integer(int32) :: OMP_get_num_threads,OMP_get_thread_num

      real(real64) :: RanNum, AcceptRate, BestSolCriterion, Chrom(nParam)

      logical :: ModeAvg, ModeMax, BestSolChanged

      class(AlphaEvolveSol), allocatable :: HoldSol
      !class(*), allocatable :: HoldSol

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

      allocate(HoldSol, source=BestSol)

      LastSampPrint = 0

      if (ModeAvg) then
        BestSolCriterion = 0.0d0
        BestSolChanged = .true.
        AcceptRate = 1.0d0
      end if

      if (ModeMax) then
        BestSolCriterion = -huge(RanNum)
        BestSol%Criterion = BestSolCriterion
        BestSolChanged = .false.
        AcceptRate = 0.0d0
      end if

      ! --- Printout log header ---

      open(newunit=LogUnit, file=trim(LogFile), status="unknown")
      call HoldSol%LogHead(LogUnit=LogUnit)

      ! --- Initialise with the provided solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do Samp = 1, nInit
          call HoldSol%Evaluate(Init(:, Samp), CritType)
          if      (ModeAvg) then
            if (Samp == 1) then
              call BestSol%Assign(HoldSol)
            else
              call BestSol%UpdateMean(HoldSol, Samp)
            end if
          else if (ModeMax) then
            if (HoldSol%Criterion > BestSolCriterion) then
              call BestSol%Assign(HoldSol)
              AcceptRate = AcceptRate + 1.0d0
            end if
          end if
        end do
      else
        nInit = 1
      end if

      ! --- Search ---

      !> @todo: parallelise this loop?
      do Samp = nInit, nSamp

        BestSolChanged = .false.

        ! --- Generate a competitor ---

        call random_number(Chrom(:))

        ! --- Evaluate and Select ---

        ! Merit of the competitor
        call HoldSol%Evaluate(Chrom(:), CritType)

        if      (ModeAvg) then
          ! Update the mean
          if (Samp == 1) then
            call BestSol%Assign(HoldSol)
          else
            call BestSol%UpdateMean(HoldSol, Samp)
          end if
          BestSolChanged = .true.
        else if (ModeMax) then
          ! If the competitor is better, keep it
          if (HoldSol%Criterion > BestSol%Criterion) then
            call BestSol%Assign(HoldSol)
            BestSolChanged = .true.
            AcceptRate = AcceptRate + 1.0d0
          end if
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Samp == 1) .or. ((Samp - LastSampPrint) >= nSampPrint)) then
            if      (ModeAvg) then
              call BestSol%Log(LogUnit, Samp, AcceptRate)
            else if (ModeMax) then
              AcceptRate = AcceptRate / (Samp - LastSampPrint)
              call BestSol%Log(LogUnit, Samp, AcceptRate)
              AcceptRate = 0.0d0
            end if
            LastSampPrint = Samp
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Samp, nSampStop) == 0) then
          if ((BestSol%Criterion - BestSolCriterion) > StopTolerance) then
            BestSolCriterion = BestSol%Criterion
          else
            write(STDOUT, "(5a)") "NOTE: Criterion did not improve for ", &
              trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nSampStop)), &
              " samples. Stopping the random search."
            write(STDOUT, "(a)") " "
            exit
          end if
        end if

      end do ! Samp

      ! --- The winner solution ---

      if (ModeMax) then
        AcceptRate = AcceptRate / (Samp - LastSampPrint)
      end if
      call BestSol%Log(LogUnit, Samp, AcceptRate)
      write(STDOUT, "(a)") " "
      close(LogUnit)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Initialise AlphaEvolve solution
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    subroutine InitialiseAlphaEvolveSol(This)
      implicit none

      ! Argument
      class(AlphaEvolveSol), intent(out) :: This !< solution

      ! Initialisation
      This%Criterion = 0.0d0
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Assign one AlphaEvolve solution to another
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    subroutine AssignAlphaEvolveSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(out) :: Out !< @return output solution
      class(AlphaEvolveSol), intent(in)  :: In  !< input solution

      ! Assignments
      Out%Criterion = In%Criterion
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Update (running) mean of AlphaEvolve solutions for random search
    !!          with Mode=avg
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Solution with average metrics
    !---------------------------------------------------------------------------
    subroutine UpdateMeanAlphaEvolveSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(inout) :: This !< solution
      class(AlphaEvolveSol), intent(in)    :: Add  !< addition
      integer(int32), intent(in)           :: n    !< number of previous solutions mean was based on

      ! Other
      real(real64) :: kR

      ! Updates
      kR = (dble(n) - 1.0d0) / n

      This%Criterion = This%Criterion * kR + Add%Criterion / n
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Evaluate AlphaEvolve solution
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    subroutine EvaluateAlphaEvolveSol(This, Chrom, CritType) ! Chrom and CritType not used here
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(inout) :: This      !< @return solution
      real(real64), intent(inout)          :: Chrom(:)  !< internal representation of the solution
      character(len=*), intent(in)         :: CritType  !< type of criterion; not used here

      !> @todo is not Chrom now part of solution
      ! Initialize the solution
      call This%Initialise()

      ! Criterion (just a random number here for simplicity)
      call random_number(This%Criterion)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print log head
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print log head to STDOUT and optionally file unit
    !---------------------------------------------------------------------------
    subroutine LogHeadAlphaEvolveSol(LogUnit)
      implicit none
      integer(int32), intent(in), optional :: LogUnit !< log file unit

      character(len=12) :: ColnameLogStdout(3)
      character(len=22) :: ColnameLogUnit(3)
      !                      123456789012
      ColnameLogStdout(1) = "   Iteration"
      ColnameLogStdout(2) = "  AcceptRate"
      ColnameLogStdout(3) = "   Criterion"
      !                    1234567890123456789012
      ColnameLogUnit(1) = "             Iteration"
      ColnameLogUnit(2) = "            AcceptRate"
      ColnameLogUnit(3) = "             Criterion"
      write(STDOUT, "(3a12)")    ColnameLogStdout(:)
      if (present(LogUnit)) then
        write(LogUnit, "(3a22)") ColnameLogUnit(:)
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print log
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print log to STDOUT and optionally file unit
    !---------------------------------------------------------------------------
    subroutine LogAlphaEvolveSol(This, LogUnit, Gen, AcceptRate)
      implicit none
      class(AlphaEvolveSol), intent(in)    :: This       !< solution
      integer(int32), intent(in), optional :: LogUnit    !< log file unit
      integer(int32), intent(in)           :: Gen        !< generation/iteration
      real(real64), intent(in)             :: AcceptRate !< acceptance rate
      write(STDOUT,  "(i12, 2(1x, f11.5))")     Gen, AcceptRate, This%Criterion
      if (present(LogUnit)) then
        write(LogUnit, "(i22, 2(1x, es21.14))") Gen, AcceptRate, This%Criterion
      end if
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print population log head
    !> @details This is meant to log all the evaluated solutions and not just the
    !!          best one as LogHeadAlphaEvolveSol
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print population log head to STDOUT and optionally file unit
    !---------------------------------------------------------------------------
    subroutine LogPopHeadAlphaEvolveSol(LogPopUnit)
      implicit none
      integer(int32), intent(in) :: LogPopUnit !< log file unit

      character(len=22) :: ColnameLogPopUnit(3)
      !                       1234567890123456789012
      ColnameLogPopUnit(1) = "             Iteration"
      ColnameLogPopUnit(2) = "              Solution"
      ColnameLogPopUnit(3) = "             Criterion"
      write(LogPopUnit, "(3a22)") ColnameLogPopUnit(:)
    end subroutine

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Print population log
    !> @details This is meant to log all the evaluated solutions and not just the
    !!          best one as LogAlphaEvolveSol
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Print population log to STDOUT and optionally file unit
    !---------------------------------------------------------------------------
    subroutine LogPopAlphaEvolveSol(This, LogPopUnit, Gen, i)
      implicit none
      class(AlphaEvolveSol), intent(in) :: This       !< solution
      integer(int32), intent(in)        :: LogPopUnit !< population log file unit
      integer(int32), intent(in)        :: Gen        !< generation/iteration
      integer(int32), intent(in)        :: i          !< solution id
      write(LogPopUnit, "(2(i22, 1x), es21.14)") Gen, i, This%Criterion
    end subroutine

    !###########################################################################
end module

!###############################################################################
