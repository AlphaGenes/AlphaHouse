
!###############################################################################

module AlphaEvolveMod

  use ISO_Fortran_Env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use AlphaHouseMod, only : Int2Char, Real2Char, ToLower

  implicit none

  type :: AlphaEvolveSol
    real(real64) :: Criterion
    contains
      procedure         :: Initialise    => InitialiseAlphaEvolveSol
      procedure         :: Assign        => AssignAlphaEvolveSol
      procedure         :: UpdateMean    => UpdateMeanAlphaEvolveSol
      procedure         :: CalcCriterion => CalcCriterionAlphaEvolveSol
      procedure, nopass :: LogHead       => LogHeadAlphaEvolveSol
      procedure         :: Log           => LogAlphaEvolveSol
  end type

  private
  ! Types
  public :: AlphaEvolveSol
  ! Methods
  public :: DifferentialEvolution, RandomSearch

  contains
    !###########################################################################

    subroutine DifferentialEvolution(nParam, nSol, Init, nGen, nGenBurnIn, nGenStop,&
      StopTolerance, nGenPrint, File, CritType, CRBurnIn, CRLate, FBase, FHigh1, FHigh2,&
      BestSol)

      ! Evolutionary Algorithm - Differential Evolution (continuous representation of solution)

      implicit none

      ! Arguments
      integer(int32), intent(in)             :: nParam         ! No. of parameters in a solution
      integer(int32), intent(in)             :: nSol           ! No. of solutions to test each generation
      real(real64), intent(in), optional     :: Init(:,:)      ! Initial solutions to start with
      integer(int32), intent(in)             :: nGen           ! No. of generations to run
      integer(int32), intent(in)             :: nGenBurnIn     ! No. of generations with more
      integer(int32), intent(in)             :: nGenStop       ! Stop after no progress for nGenStop
      real(real64), intent(in)               :: StopTolerance  ! Stopping tolerance
      integer(int32), intent(in)             :: nGenPrint      ! Print changed solution every nGenPrint
      character(len=*), intent(in)           :: File           ! Which file to write to
      character(len=*), intent(in), optional :: CritType       ! Passed to CalcCriterion
      real(real64), intent(in), optional     :: CRBurnIn       ! Crossover rate for nGenBurnIn
      real(real64), intent(in), optional     :: CRLate         ! Crossover rate
      real(real64), intent(in), optional     :: FBase          ! F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh1         ! F is multiplier of difference used to mutate
      real(real64), intent(in), optional     :: FHigh2         ! F is multiplier of difference used to mutate
      class(AlphaEvolveSol), intent(out)     :: BestSol        ! The best evolved solution
      !class(*), intent(inout)                :: BestSol        ! The best evolved solution

      ! Other
      integer(int32) :: nInit, Param, ParamLoc, Gen, LastGenPrint, Unit
      integer(int32) :: i, a, b, c, j
      ! integer(int32) :: OMP_get_num_threads,OMP_get_thread_num

      real(real64) :: RanNum, FInt, FBaseInt, FHigh1Int, FHigh2Int, CRInt, CRBurnInInt, CRLateInt
      real(real64) :: AcceptRate, OldChrom(nParam, nSol), NewChrom(nParam, nSol), Chrom(nParam)
      real(real64) :: BestSolCriterion

      logical :: DiffOnly, BestSolChanged

      class(AlphaEvolveSol), allocatable :: Sol(:), HoldSol
      !class(*), allocatable :: Sol(:), HoldSol

      ! --- Allocate and Initialize ---

      allocate(Sol(nSol), source=BestSol)
      allocate(HoldSol, source=BestSol)

      LastGenPrint = 0
      BestSolCriterion = -huge(RanNum)
      BestSol%Criterion = BestSolCriterion

      ! --- Printout log header ---

      open(newunit=Unit, file=trim(File), status="unknown")
      call HoldSol%LogHead(LogUnit=Unit)

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
          call Sol(i)%CalcCriterion(OldChrom(:,i), CritType)
        end do
        nInit = i
      else
        nInit = 1
      end if
      do i = nInit, nSol
        call random_number(OldChrom(:,i))
        call Sol(i)%CalcCriterion(OldChrom(:,i), CritType)
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
        !       - some variables are local, say a, b, ...
        !       - global variable is NewChrom, but is indexed with i so this should not be a problem
        !       - AcceptRate needs to be in sync between the threads!!!
        !       - we relly on random_number a lot here and updating the RNG state for each thread can be slow
        !         and I (GG) am also not sure if we should not have thread specific RNGs
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

          call HoldSol%CalcCriterion(Chrom, CritType)     ! Merit of competitor
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

        AcceptRate = AcceptRate / dble(nSol)

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
            call BestSol%Log(Unit, Gen, AcceptRate)
          end if
        end if

        ! --- Test if solution is improving to continue or stop early ---

        if (mod(Gen, nGenStop) == 0) then
          if ((BestSol%Criterion - BestSolCriterion) > StopTolerance) then ! TODO
            BestSolCriterion = BestSol%Criterion
          else
            write(STDOUT, "(5a)") "NOTE: Criterion did not improve for ", &
              trim(Real2Char(StopTolerance)), " in the last ", trim(Int2Char(nGenStop)), &
              " generations. Stopping the optimisation."
            write(STDOUT, "(a)") " "
            exit
          end if
        end if

      end do ! Gen

      ! --- The winner solution ---

      call BestSol%Log(Unit, Gen, AcceptRate)
      write(STDOUT, "(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################

    subroutine RandomSearch(Mode, nParam, Init, nSamp, nSampStop, StopTolerance, &
      nSampPrint, File, CritType, BestSol)

      ! Random search

      implicit none

      ! Arguments
      character(len=*), intent(in)           :: Mode           ! Mode of search (max or avg)
      integer(int32), intent(in)             :: nParam         ! No. of parameters in a solution
      real(real64), intent(inout), optional  :: Init(:,:)      ! Initial solutions to test
      integer(int32), intent(in)             :: nSamp          ! No. of samples to test
      integer(int32), intent(in), optional   :: nSampStop      ! Stop after no progress for nSampStop
      real(real64), intent(in)               :: StopTolerance  ! Stopping tolerance
      integer(int32), intent(in)             :: nSampPrint     ! Print changed solution every nSampPrint
      character(len=*), intent(in)           :: File           ! Which file to write to
      character(len=*), intent(in), optional :: CritType       ! Passed to CalcCriterion
      class(AlphaEvolveSol), intent(out)     :: BestSol        ! The best found solution
      !class(*), intent(inout)                :: BestSol        ! The best found solution

      ! Other
      integer(int32) :: nInit, Samp, LastSampPrint, Unit
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

      open(newunit=Unit, file=trim(File), status="unknown")
      call HoldSol%LogHead(LogUnit=Unit)

      ! --- Initialise with the provided solutions ---

      if (present(Init)) then
        nInit = size(Init, dim=2)
        do Samp = 1, nInit
          call HoldSol%CalcCriterion(Init(:, Samp), CritType)
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

      ! TODO: parallelise this loop?
      do Samp = nInit, nSamp

        BestSolChanged = .false.

        ! --- Generate a competitor ---

        call random_number(Chrom(:))

        ! --- Evaluate and Select ---

        ! Merit of the competitor
        call HoldSol%CalcCriterion(Chrom(:), CritType)

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
              call BestSol%Log(Unit, Samp, AcceptRate)
            else if (ModeMax) then
              AcceptRate = AcceptRate / dble(Samp - LastSampPrint)
              call BestSol%Log(Unit, Samp, AcceptRate)
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
        AcceptRate = AcceptRate / dble(Samp - LastSampPrint)
      end if
      call BestSol%Log(Unit, Samp, AcceptRate)
      write(STDOUT, "(a)") " "
      close(Unit)
    end subroutine

    !###########################################################################

    subroutine InitialiseAlphaEvolveSol(This)
      implicit none

      ! Argument
      class(AlphaEvolveSol), intent(out) :: This

      ! Initialisation
      This%Criterion = 0.0d0
    end subroutine

    !###########################################################################

    subroutine AssignAlphaEvolveSol(Out, In)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(out) :: Out
      class(AlphaEvolveSol), intent(in)  :: In

      ! Assignments
      Out%Criterion = In%Criterion
    end subroutine

    !###########################################################################

    subroutine UpdateMeanAlphaEvolveSol(This, Add, n)
      implicit none

      ! Arguments
      class(AlphaEvolveSol), intent(inout) :: This
      class(AlphaEvolveSol), intent(in)    :: Add
      integer(int32), intent(in)           :: n

      ! Other
      real(real64) :: nR, kR

      ! Updates
      nR = dble(n)
      kR = (nR - 1.0d0) / nR

      This%Criterion = This%Criterion * kR + Add%Criterion / nR
    end subroutine

    !###########################################################################

    subroutine CalcCriterionAlphaEvolveSol(This, Chrom, CritType) ! Chrom and CritType not used here
      implicit none

      ! Arguments
      class(AlphaEvolveSol)        :: This      ! Solution
      real(real64), intent(inout)  :: Chrom(:)  ! Internal representation of the solution
      character(len=*), intent(in) :: CritType  ! Type of criterion; not used here

      ! Initialize the solution
      call This%Initialise()

      ! Criterion (just a random number here for simplicity)
      call random_number(This%Criterion)
    end subroutine

    !###########################################################################

    subroutine LogHeadAlphaEvolveSol(LogUnit)
      implicit none
      integer(int32), intent(in), optional :: LogUnit

      character(len=12) :: colnamelogstdout(3)
      character(len=22) :: colnamelogunit(3)
      !                      123456789012
      colnamelogstdout(1) = "        Step"
      colnamelogstdout(2) = "  AcceptRate"
      colnamelogstdout(3) = "   Criterion"
      !                    1234567890123456789012
      colnamelogunit(1) = "                  Step"
      colnamelogunit(2) = "            AcceptRate"
      colnamelogunit(3) = "             Criterion"
      write(STDOUT, "(3a12)")    colnamelogstdout(:)
      if (present(LogUnit)) then
        write(LogUnit, "(3a22)") colnamelogunit(:)
      end if
    end subroutine

    !###########################################################################

    subroutine LogAlphaEvolveSol(This, LogUnit, Gen, AcceptRate)
      implicit none
      class(AlphaEvolveSol), intent(in)    :: This
      integer(int32), intent(in), optional :: LogUnit
      integer(int32), intent(in)           :: Gen
      real(real64), intent(in)             :: AcceptRate
      write(STDOUT,  "(i12, 2(1x, f11.5))")     Gen, AcceptRate, This%Criterion
      if (present(LogUnit)) then
        write(LogUnit, "(i22, 2(1x, es21.14))") Gen, AcceptRate, This%Criterion
      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################
