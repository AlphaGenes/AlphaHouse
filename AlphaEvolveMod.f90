
!###############################################################################

module AlphaEvolveMod

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use AlphaHouseMod, only : Int2Char, Real2Char, ToLower

  implicit none

  type :: EvolveCrit
    ! TODO: How do I get this done generically? This is a type specific for AlphaMate
    !       I have created an extended type in AlphaMateModule of base EvolveCrit here,
    !       but got bitten through using polymorphic (class() stuff) objects, where I
    !       am not allowed to do a=b etc.
    real(real64)                :: Value
    real(real64)                :: Penalty
    real(real64)                :: Gain
    real(real64)                :: GainStand
    real(real64)                :: PopInb
    real(real64)                :: RatePopInb
    real(real64)                :: PrgInb
    real(real64), allocatable   :: GenericIndVal(:)
    real(real64), allocatable   :: GenericMatVal(:)
    real(real64)                :: Cost
    integer(int32), allocatable :: nVec(:)
    real(real64), allocatable   :: xVec(:)
    integer(int32), allocatable :: MatingPlan(:,:)
    real(real64), allocatable   :: GenomeEdit(:)
    contains
      procedure :: AssignEvolveCrit
      generic   :: assignment(=) => AssignEvolveCrit
      procedure :: UpdateMeanEvolveCrit
      generic   :: UpdateMean => UpdateMeanEvolveCrit
  end type

  private
  ! Types
  public :: EvolveCrit
  ! Methods
  public :: DifferentialEvolution, RandomSearch

  contains

    !###########################################################################

    subroutine AssignEvolveCrit(Out, In)
      implicit none
      class(EvolveCrit), intent(out) :: Out
      class(EvolveCrit), intent(in)  :: In
      Out%Value           = In%Value
      Out%Penalty         = In%Penalty
      Out%Gain            = In%Gain
      Out%GainStand       = In%GainStand
      Out%PopInb          = In%PopInb
      Out%RatePopInb      = In%RatePopInb
      Out%PrgInb          = In%PrgInb
      if (allocated(In%GenericIndVal)) then
        allocate(Out%GenericIndVal(size(In%GenericIndVal)))
        Out%GenericIndVal = In%GenericIndVal
      end if
      if (allocated(In%GenericMatVal)) then
        allocate(Out%GenericMatVal(size(In%GenericMatVal)))
        Out%GenericMatVal = In%GenericMatVal
      end if
      Out%Cost            = In%Cost
      if (allocated(In%nVec)) then
        allocate(Out%nVec(size(In%nVec)))
        Out%nVec          = In%nVec
      end if
      if (allocated(In%xVec)) then
        allocate(Out%xVec(size(In%xVec)))
        Out%xVec          = In%xVec
      end if
      if (allocated(In%MatingPlan)) then
        allocate(Out%MatingPlan(size(In%MatingPlan, dim=1), size(In%MatingPlan, dim=2)))
        Out%MatingPlan    = In%MatingPlan
      end if
      if (allocated(In%GenomeEdit)) then
        allocate(Out%GenomeEdit(size(In%GenomeEdit)))
        Out%GenomeEdit    = In%GenomeEdit
      end if
    end subroutine

    !###########################################################################

    subroutine UpdateMeanEvolveCrit(This, Add, n)
      implicit none
      ! Arguments
      class(EvolveCrit), intent(inout) :: This
      class(EvolveCrit), intent(in)    :: Add
      integer(int32), intent(in)       :: n
      ! Other
      real(real64) :: nR, kR

      nR = dble(n)
      kR = (nR - 1.0d0) / nR

      This%Value           = This%Value         * kR + Add%Value         / nR
      This%Penalty         = This%Penalty       * kR + Add%Penalty       / nR
      This%Gain            = This%Gain          * kR + Add%Gain          / nR
      This%GainStand       = This%GainStand     * kR + Add%GainStand     / nR
      This%PopInb          = This%PopInb        * kR + Add%PopInb        / nR
      This%RatePopInb      = This%RatePopInb    * kR + Add%RatePopInb    / nR
      This%PrgInb          = This%PrgInb        * kR + Add%PrgInb        / nR
      if (allocated(This%GenericIndVal)) then
        This%GenericIndVal = This%GenericIndVal * kR + Add%GenericIndVal / nR
      end if
      if (allocated(This%GenericMatVal)) then
        This%GenericMatVal = This%GenericMatVal * kR + Add%GenericMatVal / nR
      end if
      This%Cost            = This%Cost          * kR + Add%Cost          / nR
      if (allocated(This%nVec)) then
        This%nVec          = This%nVec          * kR + Add%nVec          / nR
      end if
      if (allocated(This%xVec)) then
        This%xVec          = This%xVec          * kR + Add%xVec          / nR
      end if
      if (allocated(This%MatingPlan)) then
        This%MatingPlan    = This%MatingPlan    * kR + Add%MatingPlan    / nR
      end if
      if (allocated(This%GenomeEdit)) then
        This%GenomeEdit    = This%GenomeEdit    * kR + Add%GenomeEdit    / nR
      end if
      return
    end subroutine

    !###########################################################################

    subroutine DifferentialEvolution(nParam, nSol, Init, nGen, nGenBurnIn, nGenStop,&
      StopTolerance, nGenPrint, File, CritType, CRBurnIn, CRLate, FBase, FHigh1, FHigh2,&
      CalcCriterion, LogHead, Log, BestCriterion)

      ! Evolutionary Algorithm - Differential Evolution (continuous representation of solution)

      implicit none

      ! Arguments
      integer(int32), intent(in)          :: nParam         ! No. of parameters in a solution
      integer(int32), intent(in)          :: nSol           ! No. of solutions to test each generation
      real(real64), intent(in), optional  :: Init(:,:)      ! Initial solutions to start with
      integer(int32), intent(in)          :: nGen           ! No. of generations to run
      integer(int32), intent(in)          :: nGenBurnIn     ! No. of generations with more
      integer(int32), intent(in)          :: nGenStop       ! Stop after no progress for nGenStop
      real(real64), intent(in)            :: StopTolerance  ! Stopping tolerance
      integer(int32), intent(in)          :: nGenPrint      ! Print changed solution every nGenPrint
      character(len=*), intent(in)        :: File           ! Which file to write to
      character(len=*), intent(in)        :: CritType       ! Passed to CalcCriterion
      real(real64), intent(in), optional  :: CRBurnIn       ! Crossover rate for nGenBurnIn
      real(real64), intent(in), optional  :: CRLate         ! Crossover rate
      real(real64), intent(in), optional  :: FBase          ! F is multiplier of difference used to mutate
      real(real64), intent(in), optional  :: FHigh1         ! F is multiplier of difference used to mutate
      real(real64), intent(in), optional  :: FHigh2         ! F is multiplier of difference used to mutate
      type(EvolveCrit), intent(out)       :: BestCriterion  ! Criterion for the best solution

      interface
        function CalcCriterion(Sol, CritType) result(Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          real(real64), intent(inout)  :: Sol(:)
          character(len=*), intent(in) :: CritType
          type(EvolveCrit)             :: Criterion
        end function

        subroutine LogHead(LogUnit)
          use ISO_Fortran_Env
          integer(int32), intent(in) :: LogUnit
        end subroutine

        subroutine Log(LogUnit, Gen, AcceptRate, Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          integer(int32), intent(in)   :: LogUnit
          integer(int32), intent(in)   :: Gen
          real(real64), intent(in)     :: AcceptRate
          type(EvolveCrit), intent(in) :: Criterion
        end subroutine
      end interface

      ! Other
      integer(int32) :: nInit, Param, ParamLoc, Sol, Gen, LastGenPrint, Unit
      integer(int32) :: SolA, SolB, SolC, BestSol
      ! integer(int32) :: OMP_get_num_threads,OMP_get_thread_num

      real(real64) :: RanNum, FInt, FBaseInt, FHigh1Int, FHigh2Int, CRInt, CRBurnInInt, CRLateInt
      real(real64) :: AcceptRate, OldChrom(nParam, nSol), NewChrom(nParam, nSol), Chrom(nParam)
      real(real64) :: BestCriterionStopValue

      logical :: DiffOnly, BestSolChanged

      type(EvolveCrit) :: Criterion(nSol), CriterionHold

      ! --- Initialize ---

      LastGenPrint = 0
      BestCriterion%Value = -huge(RanNum)
      BestCriterionStopValue = -huge(RanNum)

      ! --- Printout ---

      open(newunit=Unit, file=trim(File), status="unknown")
      call LogHead(LogUnit=Unit)

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
        do Sol = 1, nInit
          OldChrom(:,Sol) = Init(:,Sol)
          Criterion(Sol) = CalcCriterion(OldChrom(:,Sol), CritType)
        end do
        nInit = Sol
      else
        nInit = 1
      end if
      do Sol = nInit, nSol
        call random_number(OldChrom(:,Sol))
        Criterion(Sol) = CalcCriterion(OldChrom(:,Sol), CritType)
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

        ! TODO: Paralelize this loop?
        !       The main reason would be to speed up the code as CalcCriterion() might take quite some time for larger problems
        !       - some variables are local, say SolA, SolB, ...
        !       - global variable is NewChrom, but is indexed with Sol so this should not be a problem
        !       - AcceptRate needs to be in sync between the threads!!!
        !       - we relly on random_number a lot here and updating the RNG state for each thread can be slow
        !         and I (GG) am also not sure if we should not have thread specific RNGs
        BestSolChanged = .false.
        AcceptRate = 0.0d0

        ! call OMP_set_num_threads(1)

        ! $OMP PARALLEL DO DEFAULT(PRIVATE)
        do Sol = 1, nSol

          ! print *, "#Threads: ", OMP_get_num_threads(), "Thread; ", OMP_get_thread_num()+1, ", Solution: ", Sol

          ! --- Mutate and recombine ---

          ! Get three different solutions
          SolA = Sol
          do while (SolA == Sol)
            call random_number(RanNum)
            SolA = int(RanNum * nSol) + 1
          end do
          SolB = Sol
          do while ((SolB == Sol) .or. (SolB == SolA))
            call random_number(RanNum)
            SolB = int(RanNum * nSol) + 1
          end do
          SolC = Sol
          do while ((SolC == Sol) .or. (SolC == SolA) .or. (SolC == SolB))
            call random_number(RanNum)
            SolC = int(RanNum * nSol) + 1
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
                Chrom(Param) = OldChrom(Param, SolC) + FInt * (OldChrom(Param, SolA) - OldChrom(Param, SolB))
              else
                ! Non-differential mutation (to avoid getting stuck)
                call random_number(RanNum)
                if (RanNum < 0.5d0) then
                  call random_number(RanNum)
                  Chrom(Param) = OldChrom(Param, SolC) * (0.9d0 + 0.2d0 * RanNum)
                else
                  call random_number(RanNum)
                  Chrom(Param) = OldChrom(Param, SolC) + 0.01d0 * FInt * (OldChrom(Param, SolA) + 0.01d0) * (RanNum - 0.5d0)
                end if
              end if
            else
              ! Do not recombine
              Chrom(Param) = OldChrom(Param, Sol)
            end if
            Param = Param + 1
            if (Param > nParam) then
              Param = Param - nParam
            end if
          end do

          ! --- Evaluate and Select ---

          CriterionHold = CalcCriterion(Chrom, CritType)        ! Merit of competitor
          if (CriterionHold%Value >= Criterion(Sol)%Value) then ! If competitor is better or equal, keep it
            NewChrom(:,Sol) = Chrom(:)                          !   ("equal" to force evolution)
            Criterion(Sol) = CriterionHold
            ! $OMP ATOMIC
            AcceptRate = AcceptRate + 1.0d0
          else
            NewChrom(:,Sol) = OldChrom(:,Sol)                   ! Else keep the old solution
          end if
        end do ! Sol
        ! $OMP END PARALLEL DO

        AcceptRate = AcceptRate / dble(nSol)

        ! --- New parents ---

        OldChrom(:,:) = NewChrom(:,:)

        ! --- The current best solution ---

        BestSol = maxloc(Criterion(:)%Value, dim=1)
        if (Criterion(BestSol)%Value > BestCriterion%Value) then
          BestSolChanged = .true.
          BestCriterion = Criterion(BestSol)
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Gen == 1) .or. ((Gen - LastGenPrint) >= nGenPrint)) then
            LastGenPrint = Gen
            call Log(Unit, Gen, AcceptRate, BestCriterion)
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Gen, nGenStop) == 0) then
          if ((BestCriterion%Value - BestCriterionStopValue) > StopTolerance) then
            BestCriterionStopValue = BestCriterion%Value
          else
            write(STDOUT, "(5a)") "NOTE: Objective did not improve for ", &
              trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nGenStop)), &
              " generations. Stopping the optimisation."
            write(STDOUT, "(a)") " "
            exit
          end if
        end if

      end do ! Gen

      ! --- The winner solution ---

      call Log(Unit, Gen, AcceptRate, BestCriterion)
      write(STDOUT, "(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################

    subroutine RandomSearch(Mode, nParam, Init, nSamp, nSampStop, StopTolerance, &
      nSampPrint, File, CritType, CalcCriterion, LogHead, Log, BestCriterion)

      ! Random search

      implicit none

      ! Arguments
      character(len=*), intent(in)          :: Mode           ! Mode of search (max or avg)
      integer(int32), intent(in)            :: nParam         ! No. of parameters in a solution
      real(real64), intent(inout), optional :: Init(:,:)      ! Initial solutions to test
      integer(int32), intent(in)            :: nSamp          ! No. of samples to test
      integer(int32), intent(in), optional  :: nSampStop      ! Stop after no progress for nSampStop
      real(real64), intent(in)              :: StopTolerance  ! Stopping tolerance
      integer(int32), intent(in)            :: nSampPrint     ! Print changed solution every nSampPrint
      character(len=*), intent(in)          :: File           ! Which file to write to
      character(len=*), intent(in)          :: CritType       ! Passed to CalcCriterion
      type(EvolveCrit), intent(out)         :: BestCriterion  ! Criterion for the best solution

      interface
        function CalcCriterion(Sol, CritType) result(Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          real(real64), intent(inout)  :: Sol(:)
          character(len=*), intent(in) :: CritType
          type(EvolveCrit)             :: Criterion
        end function

        subroutine LogHead(LogUnit)
          use ISO_Fortran_Env
          integer(int32), intent(in) :: LogUnit
        end subroutine

        subroutine Log(LogUnit, Gen, AcceptRate, Criterion)
          use ISO_Fortran_Env
          import :: EvolveCrit
          integer(int32), intent(in)   :: LogUnit
          integer(int32), intent(in)   :: Gen
          real(real64), intent(in)     :: AcceptRate
          type(EvolveCrit), intent(in) :: Criterion
        end subroutine
      end interface

      ! Other
      integer(int32) :: nInit, Samp, LastSampPrint, Unit
      ! integer(int32) :: OMP_get_num_threads,OMP_get_thread_num

      real(real64) :: RanNum, AcceptRate, BestCriterionStopValue, Sol(nParam)

      logical :: ModeAvg, ModeMax, BestSolChanged

      type(EvolveCrit) :: CriterionHold

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

      ! --- Initialize ---

      LastSampPrint = 0

      if (ModeAvg) then
        BestCriterion%Value = 0.0d0
        BestCriterionStopValue = 0.0d0
        BestSolChanged = .true.
        AcceptRate = 1.0d0
      end if

      if (ModeMax) then
        BestCriterion%Value = -huge(RanNum)
        BestCriterionStopValue = -huge(RanNum)
        BestSolChanged = .false.
        AcceptRate = 0.0d0
      end if

      ! --- Printout ---

      open(newunit=Unit, file=trim(File), status="unknown")
      call LogHead(LogUnit=Unit)

      ! --- Initialise with the provided solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do Samp = 1, nInit
          CriterionHold = CalcCriterion(Init(:, Samp), CritType)
          if      (ModeAvg) then
            if (Samp == 1) then
              BestCriterion = CriterionHold
            else
              call BestCriterion%UpdateMean(CriterionHold, Samp)
            end if
          else if (ModeMax) then
            if (CriterionHold%Value > BestCriterion%Value) then
              BestCriterion = CriterionHold
              AcceptRate = AcceptRate + 1.0d0
            end if
          end if
        end do
      else
        nInit = 1
      end if

      ! --- Search ---

      ! TODO: parallelise this loop?
      do Samp = nInit, nSamp

        BestSolChanged = .false.

        ! --- Generate a competitor ---

        call random_number(Sol(:))

        ! --- Evaluate and Select ---

        ! Merit of the competitor
        CriterionHold = CalcCriterion(Sol(:), CritType)

        if      (ModeAvg) then
          ! Update the mean
          if (Samp == 1) then
            BestCriterion = CriterionHold
          else
            call BestCriterion%UpdateMean(CriterionHold, Samp)
          end if
          BestSolChanged = .true.
        else if (ModeMax) then
          ! If the competitor is better, keep it
          if (CriterionHold%Value > BestCriterion%Value) then
            BestCriterion = CriterionHold
            BestSolChanged = .true.
            AcceptRate = AcceptRate + 1.0d0
          end if
        end if

        ! --- Monitor ---

        if (BestSolChanged) then
          if ((Samp == 1) .or. ((Samp - LastSampPrint) >= nSampPrint)) then
            if      (ModeAvg) then
              call Log(Unit, Samp, AcceptRate, BestCriterion)
            else if (ModeMax) then
              AcceptRate = AcceptRate / dble(Samp - LastSampPrint)
              call Log(Unit, Samp, AcceptRate, BestCriterion)
              AcceptRate = 0.0d0
            end if
            LastSampPrint = Samp
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Samp, nSampStop) == 0) then
          if ((BestCriterion%Value - BestCriterionStopValue) > StopTolerance) then
            BestCriterionStopValue = BestCriterion%Value
          else
            write(STDOUT, "(5a)") "NOTE: Objective did not improve for ", &
              trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nSampStop)), &
              " samples. Stopping the random search."
            write(STDOUT, "(a)") " "
            exit
          end if
        end if

      end do ! Samp

      ! --- The winner solution ---

      if (ModeMax) then
        AcceptRate = AcceptRate / dble(Samp - LastSampPrint)
      end if
      call Log(Unit, Samp, AcceptRate, BestCriterion)
      write(STDOUT, "(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################
end module
