
!###############################################################################

module OrderPackModule
  ! From http://www.fortran-2000.com/rank (2016-02-15)
  ! Modified to work with ISO_Fortran_Env types and converted to
  ! functions where possible

  use ISO_Fortran_Env

  implicit none

  private
  public :: MrgRnk, RnkPar, UniSta, UniInv

  interface MrgRnk
    module procedure D_MrgRnk, R_MrgRnk, I_MrgRnk
  end interface

  interface RnkPar
    module procedure D_RnkPar, R_RnkPar, I_RnkPar
  end interface

  interface UniSta
    module procedure D_UniSta, R_UniSta, I_UniSta
  end interface

  interface UniInv
    module procedure D_UniInv, R_UniInv, I_UniInv
  end interface

  interface Nearless
    module procedure D_Nearless, R_Nearless, I_Nearless
  end interface

  contains

    !###########################################################################

    ! MrgRnk - array ranks

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks array XVALT into index array IRNGT, using merge-sort
      !> @details For performance reasons, the first 2 passes are taken out of the
      !!          standard loop, and use dedicated coding
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function D_MrgRnk (XDONT) result(IRNGT)
          implicit none
          Real(real64), Dimension(:), Intent (In)   :: XDONT !< Vector to rank
          Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

          Real(real64) :: XVALA, XVALB

          Integer(int32), Dimension (SIZE(XDONT)) :: JWRKT
          Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

          ! NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          NVAL = SIZE(XDONT)
          allocate(IRNGT(NVAL))

          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
        !
        !  Fill-in the index array, creating ordered couples
        !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
        !
        !  We will now have ordered subsets A - B - A - B - ...
        !  and merge A and B couples into     C   -   C   - ...
        !
          LMTNA = 2
          LMTNC = 4
        !
        !  First iteration. The length of the ordered subsets goes from 2 to 4
        !
          Do
             If (NVAL <= 2) Exit
        !
        !   Loop on merges of A and B into C
        !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
        !
        !   1 2 3
        !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
        !
        !   1 3 2
        !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
        !
        !   3 1 2
        !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
        !
        !   1 2 3 4
        !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
        !
        !   1 3 x x
        !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
        !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
        !
        !   3 x x x
        !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
        !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
        !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 4
             Exit
          End Do
        !
        !  Iteration loop. Each time, the length of the ordered subsets
        !  is doubled.
        !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
        !
        !   Loop on merges of A and B into C
        !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B. This line may be activated when the
        !   initial set is already close to sorted.
        !
        !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
        !
        !  One steps in the C subset, that we build in the final rank array
        !
        !  Make a copy of the rank array for the merge iteration
        !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
        !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
        !
                Do
                   IWRK = IWRK + 1
        !
        !  We still have unprocessed values in both A and B
        !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
        !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
        !
                End Do
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 2 * LMTNA
          End Do
        !
          Return
        !
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks array XVALT into index array IRNGT, using merge-sort
      !> @details For performance reasons, the first 2 passes are taken out of the
      !!          standard loop, and use dedicated coding
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function R_MrgRnk (XDONT) result(IRNGT)
          implicit none
          Real(real32), Dimension(:), Intent (In)   :: XDONT !< Vector to rank
          Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

          Real(real32) :: XVALA, XVALB

          Integer(int32), Dimension (SIZE(XDONT)) :: JWRKT
          Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

          ! NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          NVAL = SIZE(XDONT)
          allocate(IRNGT(NVAL))

          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
        !
        !  Fill-in the index array, creating ordered couples
        !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
        !
        !  We will now have ordered subsets A - B - A - B - ...
        !  and merge A and B couples into     C   -   C   - ...
        !
          LMTNA = 2
          LMTNC = 4
        !
        !  First iteration. The length of the ordered subsets goes from 2 to 4
        !
          Do
             If (NVAL <= 2) Exit
        !
        !   Loop on merges of A and B into C
        !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
        !
        !   1 2 3
        !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
        !
        !   1 3 2
        !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
        !
        !   3 1 2
        !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
        !
        !   1 2 3 4
        !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
        !
        !   1 3 x x
        !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
        !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
        !
        !   3 x x x
        !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
        !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
        !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 4
             Exit
          End Do
        !
        !  Iteration loop. Each time, the length of the ordered subsets
        !  is doubled.
        !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
        !
        !   Loop on merges of A and B into C
        !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B. This line may be activated when the
        !   initial set is already close to sorted.
        !
        !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
        !
        !  One steps in the C subset, that we build in the final rank array
        !
        !  Make a copy of the rank array for the merge iteration
        !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
        !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
        !
                Do
                   IWRK = IWRK + 1
        !
        !  We still have unprocessed values in both A and B
        !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
        !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
        !
                End Do
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 2 * LMTNA
          End Do
        !
          Return
        !
      End function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks array XVALT into index array IRNGT, using merge-sort
      !> @details For performance reasons, the first 2 passes are taken out of the
      !!          standard loop, and use dedicated coding
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function I_MrgRnk (XDONT) result(IRNGT)
          implicit none
          Integer(int32), Dimension(:), Intent (In) :: XDONT !< Vector to rank
          Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

          Integer(int32) :: XVALA, XVALB

          Integer(int32), Dimension (SIZE(XDONT)) :: JWRKT
          Integer(int32) :: LMTNA, LMTNC, IRNG1, IRNG2
          Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

          ! NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
          NVAL = SIZE(XDONT)
          allocate(IRNGT(NVAL))

          Select Case (NVAL)
          Case (:0)
             Return
          Case (1)
             IRNGT (1) = 1
             Return
          Case Default
             Continue
          End Select
        !
        !  Fill-in the index array, creating ordered couples
        !
          Do IIND = 2, NVAL, 2
             If (XDONT(IIND-1) <= XDONT(IIND)) Then
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             Else
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             End If
          End Do
          If (Modulo(NVAL, 2) /= 0) Then
             IRNGT (NVAL) = NVAL
          End If
        !
        !  We will now have ordered subsets A - B - A - B - ...
        !  and merge A and B couples into     C   -   C   - ...
        !
          LMTNA = 2
          LMTNC = 4
        !
        !  First iteration. The length of the ordered subsets goes from 2 to 4
        !
          Do
             If (NVAL <= 2) Exit
        !
        !   Loop on merges of A and B into C
        !
             Do IWRKD = 0, NVAL - 1, 4
                If ((IWRKD+4) > NVAL) Then
                   If ((IWRKD+2) >= NVAL) Exit
        !
        !   1 2 3
        !
                   If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
        !
        !   1 3 2
        !
                   If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
        !
        !   3 1 2
        !
                   Else
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   End If
                   Exit
                End If
        !
        !   1 2 3 4
        !
                If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
        !
        !   1 3 x x
        !
                If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   Else
        !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   End If
        !
        !   3 x x x
        !
                Else
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                      IRNGT (IWRKD+2) = IRNG1
                      If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
        !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      Else
        !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      End If
                   Else
        !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   End If
                End If
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 4
             Exit
          End Do
        !
        !  Iteration loop. Each time, the length of the ordered subsets
        !  is doubled.
        !
          Do
             If (LMTNA >= NVAL) Exit
             IWRKF = 0
             LMTNC = 2 * LMTNC
        !
        !   Loop on merges of A and B into C
        !
             Do
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                If (IWRKF >= NVAL) Then
                   If (JINDA >= NVAL) Exit
                   IWRKF = NVAL
                End If
                IINDA = 1
                IINDB = JINDA + 1
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B. This line may be activated when the
        !   initial set is already close to sorted.
        !
        !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
        !
        !  One steps in the C subset, that we build in the final rank array
        !
        !  Make a copy of the rank array for the merge iteration
        !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
        !
                XVALA = XDONT (JWRKT(IINDA))
                XVALB = XDONT (IRNGT(IINDB))
        !
                Do
                   IWRK = IWRK + 1
        !
        !  We still have unprocessed values in both A and B
        !
                   If (XVALA > XVALB) Then
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      If (IINDB > IWRKF) Then
        !  Only A still with unprocessed values
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         Exit
                      End If
                      XVALB = XDONT (IRNGT(IINDB))
                   Else
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                      XVALA = XDONT (JWRKT(IINDA))
                   End If
        !
                End Do
             End Do
        !
        !  The Cs become As and Bs
        !
             LMTNA = 2 * LMTNA
          End Do
        !
          Return
        !
      End function

      !#########################################################################

    !###########################################################################

    ! RnkPar - array ranks (for n smallest values)

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks partially XVALT by IRNGT, up to order NORD
      !> @details This routine uses a pivoting strategy such as the one of finding
      !!          the median based on the quicksort algorithm, but we skew the pivot
      !!          choice to try to bring it to NORD as fast as possible. It uses 2
      !!          temporary arrays, one where it stores the indices of the values
      !!          smaller than the pivot, and the other for the indices of values
      !!          larger than the pivot that we might still need later on. It iterates
      !!          until it can bring the number of values in ILOWT to exactly NORD,
      !!          and then uses an insertion sort to rank this set, since it is
      !           supposedly small.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function D_RnkPar (XDONT, NORD) result(IRNGT)
            implicit none
            Real(real64), Dimension(:), Intent(In)    :: XDONT !< Vector to rank
            Integer(int32), Intent(In)                :: NORD  !< Number of smallest values to rank
            Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

            Real(real64) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX

            Integer(int32), Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
            Integer(int32) :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
            Integer(int32) :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
            Integer(int32) :: JLM2, JLM1, JHM2, JHM1
      !
            NDON = SIZE (XDONT)
            allocate(IRNGT(NORD))
      !
      !    First loop is used to fill-in ILOWT, IHIGT at the same time
      !
            If (NDON < 2) Then
              If (NORD >= 1) IRNGT (1) = 1
              Return
            End If
      !
      !  One chooses a pivot, best estimate possible to put fractile near
      !  mid-point of the set of low values.
      !
            If (XDONT(2) < XDONT(1)) Then
              ILOWT (1) = 2
              IHIGT (1) = 1
            Else
              ILOWT (1) = 1
              IHIGT (1) = 2
            End If
      !
            If (NDON < 3) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              Return
            End If
      !
            If (XDONT(3) <= XDONT(IHIGT(1))) Then
              IHIGT (2) = IHIGT (1)
              If (XDONT(3) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = 3
              Else
                  IHIGT (1) = 3
              End If
            Else
              IHIGT (2) = 3
            End If
      !
            If (NDON < 4) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              Return
            End If
      !
            If (XDONT(NDON) <= XDONT(IHIGT(1))) Then
              IHIGT (3) = IHIGT (2)
              IHIGT (2) = IHIGT (1)
              If (XDONT(NDON) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = NDON
              Else
                  IHIGT (1) = NDON
              End If
            Else
              if (XDONT (NDON) < XDONT (IHIGT(2))) Then
                  IHIGT (3) = IHIGT (2)
                  IHIGT (2) = NDON
              else
                  IHIGT (3) = NDON
              endif
            End If
      !
            If (NDON < 5) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              If (NORD >= 4) IRNGT (4) = IHIGT (3)
              Return
            End If
      !
            JDEB = 0
            IDEB = JDEB + 1
            JLOW = IDEB
            JHIG = 3
            XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                        (XDONT(IHIGT(3))-XDONT(ILOWT(IDEB)))
            If (XPIV >= XDONT(IHIGT(1))) Then
              XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                            (XDONT(IHIGT(2))-XDONT(ILOWT(IDEB)))
              If (XPIV >= XDONT(IHIGT(1))) &
                  XPIV = XDONT (ILOWT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                                (XDONT(IHIGT(1))-XDONT(ILOWT(IDEB)))
            End If
            XPIV0 = XPIV
      !
      !  One puts values > pivot in the end and those <= pivot
      !  at the beginning. This is split in 2 cases, so that
      !  we can skip the loop test a number of times.
      !  As we are also filling in the work arrays at the same time
      !  we stop filling in the IHIGT array as soon as we have more
      !  than enough values in ILOWT.
      !
      !
            If (XDONT(NDON) > XPIV) Then
              ICRS = 3
              Do
                  ICRS = ICRS + 1
                  If (XDONT(ICRS) > XPIV) Then
                    If (ICRS >= NDON) Exit
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
      !  One restricts further processing because it is no use
      !  to store more high values
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    Else If (ICRS >= NDON) Then
                        Exit
                    End If
                  End Do
              End If
      !
      !
            Else
      !
      !  Same as above, but this is not as easy to optimize, so the
      !  DO-loop is kept
      !
              Do ICRS = 4, NDON - 1
                  If (XDONT(ICRS) > XPIV) Then
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        If (ICRS >= NDON) Exit
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    End If
                  End Do
              End If
            End If
      !
            JLM2 = 0
            JLM1 = 0
            JHM2 = 0
            JHM1 = 0
            Do
              if (JLOW == NORD) Exit
              If (JLM2 == JLOW .And. JHM2 == JHIG) Then
      !
      !   We are oscillating. Perturbate by bringing JLOW closer by one
      !   to NORD
      !
                If (NORD > JLOW) Then
                      XMIN = XDONT (IHIGT(1))
                      IHIG = 1
                      Do ICRS = 2, JHIG
                        If (XDONT(IHIGT(ICRS)) < XMIN) Then
                            XMIN = XDONT (IHIGT(ICRS))
                            IHIG = ICRS
                        End If
                      End Do
      !
                      JLOW = JLOW + 1
                      ILOWT (JLOW) = IHIGT (IHIG)
                      IHIGT (IHIG) = IHIGT (JHIG)
                      JHIG = JHIG - 1
                  Else
                      ILOW = ILOWT (JLOW)
                      XMAX = XDONT (ILOW)
                      Do ICRS = 1, JLOW
                        If (XDONT(ILOWT(ICRS)) > XMAX) Then
                            IWRK = ILOWT (ICRS)
                            XMAX = XDONT (IWRK)
                            ILOWT (ICRS) = ILOW
                            ILOW = IWRK
                        End If
                      End Do
                      JLOW = JLOW - 1
                  End If
              End If
              JLM2 = JLM1
              JLM1 = JLOW
              JHM2 = JHM1
              JHM1 = JHIG
      !
      !   We try to bring the number of values in the low values set
      !   closer to NORD.
      !
              Select Case (NORD-JLOW)
              Case (2:)
      !
      !   Not enough values in low part, at least 2 are missing
      !
                  Select Case (JHIG)
      !!!!!           CASE DEFAULT
      !!!!!              write (*,*) "Assertion failed"
      !!!!!              STOP
      !
      !   We make a special case when we have so few values in
      !   the high values set that it is bad performance to choose a pivot
      !   and apply the general algorithm.
      !
                  Case (2)
                    If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                    Else
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                    End If
                    Exit
      !
                  Case (3)
      !
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (3)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (3) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
                    JHIG = 0
                    Do ICRS = JLOW + 1, NORD
                        JHIG = JHIG + 1
                        ILOWT (ICRS) = IHIGT (JHIG)
                    End Do
                    JLOW = NORD
                    Exit
      !
                  Case (4:)
      !
      !
                    XPIV0 = XPIV
                    IFIN = JHIG
      !
      !  One chooses a pivot from the 2 first values and the last one.
      !  This should ensure sufficient renewal between iterations to
      !  avoid worst case behavior effects.
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (IFIN)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (IFIN) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
      !
                    JDEB = JLOW
                    NWRK = NORD - JLOW
                    IWRK1 = IHIGT (1)
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = IWRK1
                    XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                            (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
      !
      !  One takes values <= pivot to ILOWT
      !  Again, 2 parts, one where we take care of the remaining
      !  high values because we might still need them, and the
      !  other when we know that we will have more than enough
      !  low values in the end.
      !
                    JHIG = 0
                    Do ICRS = 2, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                          If (JLOW >= NORD) Exit
                        Else
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = IHIGT (ICRS)
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                        End If
                    End Do
                End Select
      !
      !
              Case (1)
      !
      !  Only 1 value is missing in low part
      !
                  XMIN = XDONT (IHIGT(1))
                  IHIG = 1
                  Do ICRS = 2, JHIG
                    If (XDONT(IHIGT(ICRS)) < XMIN) Then
                        XMIN = XDONT (IHIGT(ICRS))
                        IHIG = ICRS
                    End If
                  End Do
      !
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (IHIG)
                  Exit
      !
      !
              Case (0)
      !
      !  Low part is exactly what we want
      !
                  Exit
      !
      !
              Case (-5:-1)
      !
      !  Only few values too many in low part
      !
                  IRNGT (1) = ILOWT (1)
                  Do ICRS = 2, NORD
                    IWRK = ILOWT (ICRS)
                    XWRK = XDONT (IWRK)
                    Do IDCR = ICRS - 1, 1, - 1
                        If (XWRK < XDONT(IRNGT(IDCR))) Then
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        Else
                          Exit
                        End If
                    End Do
                    IRNGT (IDCR+1) = IWRK
                  End Do
      !
                  XWRK1 = XDONT (IRNGT(NORD))
                  Do ICRS = NORD + 1, JLOW
                    If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                        XWRK = XDONT (ILOWT (ICRS))
                        Do IDCR = NORD - 1, 1, - 1
                          If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        End Do
                        IRNGT (IDCR+1) = ILOWT (ICRS)
                        XWRK1 = XDONT (IRNGT(NORD))
                    End If
                  End Do
      !
                  Return
      !
      !
              Case (:-6)
      !
      ! last case: too many values in low part
      !
                  IDEB = JDEB + 1
                  IMIL = (JLOW+IDEB) / 2
                  IFIN = JLOW
      !
      !  One chooses a pivot from 1st, last, and middle values
      !
                  If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                    IWRK = ILOWT (IDEB)
                    ILOWT (IDEB) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                  End If
                  If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
                    IWRK = ILOWT (IFIN)
                    ILOWT (IFIN) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                    If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                        IWRK = ILOWT (IDEB)
                        ILOWT (IDEB) = ILOWT (IMIL)
                        ILOWT (IMIL) = IWRK
                    End If
                  End If
                  If (IFIN <= 3) Exit
      !
                  XPIV = XDONT (ILOWT(1)) + REAL(NORD)/REAL(JLOW+NORD) * &
                                            (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))
                  If (JDEB > 0) Then
                    If (XPIV <= XPIV0) &
                        XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JLOW+NORD) * &
                                        (XDONT(ILOWT(IFIN))-XPIV0)
                  Else
                    IDEB = 1
                  End If
      !
      !  One takes values > XPIV to IHIGT
      !  However, we do not process the first values if we have been
      !  through the case when we did not have enough low values
      !
                  JHIG = 0
                  JLOW = JDEB
      !
                  If (XDONT(ILOWT(IFIN)) > XPIV) Then
                    ICRS = JDEB
                    Do
                      ICRS = ICRS + 1
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                          If (ICRS >= IFIN) Exit
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    If (ICRS < IFIN) Then
                        Do
                          ICRS = ICRS + 1
                          If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                              JLOW = JLOW + 1
                              ILOWT (JLOW) = ILOWT (ICRS)
                          Else
                              If (ICRS >= IFIN) Exit
                          End If
                        End Do
                    End If
                Else
                    Do ICRS = IDEB, IFIN
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                        End If
                    End Do
                  End If
      !
              End Select
      !
            End Do
      !
      !  Now, we only need to complete ranking of the 1:NORD set
      !  Assuming NORD is small, we use a simple insertion sort
      !
            IRNGT (1) = ILOWT (1)
            Do ICRS = 2, NORD
              IWRK = ILOWT (ICRS)
              XWRK = XDONT (IWRK)
              Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                    IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                    Exit
                  End If
              End Do
              IRNGT (IDCR+1) = IWRK
            End Do
          Return
      !
      !
      End Function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks partially XVALT by IRNGT, up to order NORD
      !> @details This routine uses a pivoting strategy such as the one of finding
      !!          the median based on the quicksort algorithm, but we skew the pivot
      !!          choice to try to bring it to NORD as fast as possible. It uses 2
      !!          temporary arrays, one where it stores the indices of the values
      !!          smaller than the pivot, and the other for the indices of values
      !!          larger than the pivot that we might still need later on. It iterates
      !!          until it can bring the number of values in ILOWT to exactly NORD,
      !!          and then uses an insertion sort to rank this set, since it is
      !           supposedly small.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function R_RnkPar (XDONT, NORD) result(IRNGT)
            implicit none
            Real(real32), Dimension(:), Intent(In)    :: XDONT !< Vector to rank
            Integer(int32), Intent(In)                :: NORD  !< Number of smallest values to rank
            Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

            Real(real32) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX

            Integer(int32), Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
            Integer(int32) :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
            Integer(int32) :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
            Integer(int32) :: JLM2, JLM1, JHM2, JHM1

            NDON = SIZE (XDONT)
            allocate(IRNGT(NORD))
      !
      !    First loop is used to fill-in ILOWT, IHIGT at the same time
      !
            If (NDON < 2) Then
              If (NORD >= 1) IRNGT (1) = 1
              Return
            End If
      !
      !  One chooses a pivot, best estimate possible to put fractile near
      !  mid-point of the set of low values.
      !
            If (XDONT(2) < XDONT(1)) Then
              ILOWT (1) = 2
              IHIGT (1) = 1
            Else
              ILOWT (1) = 1
              IHIGT (1) = 2
            End If
      !
            If (NDON < 3) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              Return
            End If
      !
            If (XDONT(3) <= XDONT(IHIGT(1))) Then
              IHIGT (2) = IHIGT (1)
              If (XDONT(3) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = 3
              Else
                  IHIGT (1) = 3
              End If
            Else
              IHIGT (2) = 3
            End If
      !
            If (NDON < 4) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              Return
            End If
      !
            If (XDONT(NDON) <= XDONT(IHIGT(1))) Then
              IHIGT (3) = IHIGT (2)
              IHIGT (2) = IHIGT (1)
              If (XDONT(NDON) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = NDON
              Else
                  IHIGT (1) = NDON
              End If
            Else
              if (XDONT (NDON) < XDONT (IHIGT(2))) Then
                  IHIGT (3) = IHIGT (2)
                  IHIGT (2) = NDON
              else
                  IHIGT (3) = NDON
              endif
            End If
      !
            If (NDON < 5) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              If (NORD >= 4) IRNGT (4) = IHIGT (3)
              Return
            End If
      !
            JDEB = 0
            IDEB = JDEB + 1
            JLOW = IDEB
            JHIG = 3
            XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                        (XDONT(IHIGT(3))-XDONT(ILOWT(IDEB)))
            If (XPIV >= XDONT(IHIGT(1))) Then
              XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                            (XDONT(IHIGT(2))-XDONT(ILOWT(IDEB)))
              If (XPIV >= XDONT(IHIGT(1))) &
                  XPIV = XDONT (ILOWT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                                (XDONT(IHIGT(1))-XDONT(ILOWT(IDEB)))
            End If
            XPIV0 = XPIV
      !
      !  One puts values > pivot in the end and those <= pivot
      !  at the beginning. This is split in 2 cases, so that
      !  we can skip the loop test a number of times.
      !  As we are also filling in the work arrays at the same time
      !  we stop filling in the IHIGT array as soon as we have more
      !  than enough values in ILOWT.
      !
      !
            If (XDONT(NDON) > XPIV) Then
              ICRS = 3
              Do
                  ICRS = ICRS + 1
                  If (XDONT(ICRS) > XPIV) Then
                    If (ICRS >= NDON) Exit
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
      !  One restricts further processing because it is no use
      !  to store more high values
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    Else If (ICRS >= NDON) Then
                        Exit
                    End If
                  End Do
              End If
      !
      !
            Else
      !
      !  Same as above, but this is not as easy to optimize, so the
      !  DO-loop is kept
      !
              Do ICRS = 4, NDON - 1
                  If (XDONT(ICRS) > XPIV) Then
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        If (ICRS >= NDON) Exit
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    End If
                  End Do
              End If
            End If
      !
            JLM2 = 0
            JLM1 = 0
            JHM2 = 0
            JHM1 = 0
            Do
              if (JLOW == NORD) Exit
              If (JLM2 == JLOW .And. JHM2 == JHIG) Then
      !
      !   We are oscillating. Perturbate by bringing JLOW closer by one
      !   to NORD
      !
                If (NORD > JLOW) Then
                      XMIN = XDONT (IHIGT(1))
                      IHIG = 1
                      Do ICRS = 2, JHIG
                        If (XDONT(IHIGT(ICRS)) < XMIN) Then
                            XMIN = XDONT (IHIGT(ICRS))
                            IHIG = ICRS
                        End If
                      End Do
      !
                      JLOW = JLOW + 1
                      ILOWT (JLOW) = IHIGT (IHIG)
                      IHIGT (IHIG) = IHIGT (JHIG)
                      JHIG = JHIG - 1
                  Else
                      ILOW = ILOWT (JLOW)
                      XMAX = XDONT (ILOW)
                      Do ICRS = 1, JLOW
                        If (XDONT(ILOWT(ICRS)) > XMAX) Then
                            IWRK = ILOWT (ICRS)
                            XMAX = XDONT (IWRK)
                            ILOWT (ICRS) = ILOW
                            ILOW = IWRK
                        End If
                      End Do
                      JLOW = JLOW - 1
                  End If
              End If
              JLM2 = JLM1
              JLM1 = JLOW
              JHM2 = JHM1
              JHM1 = JHIG
      !
      !   We try to bring the number of values in the low values set
      !   closer to NORD.
      !
              Select Case (NORD-JLOW)
              Case (2:)
      !
      !   Not enough values in low part, at least 2 are missing
      !
                  Select Case (JHIG)
      !!!!!           CASE DEFAULT
      !!!!!              write (*,*) "Assertion failed"
      !!!!!              STOP
      !
      !   We make a special case when we have so few values in
      !   the high values set that it is bad performance to choose a pivot
      !   and apply the general algorithm.
      !
                  Case (2)
                    If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                    Else
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                    End If
                    Exit
      !
                  Case (3)
      !
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (3)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (3) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
                    JHIG = 0
                    Do ICRS = JLOW + 1, NORD
                        JHIG = JHIG + 1
                        ILOWT (ICRS) = IHIGT (JHIG)
                    End Do
                    JLOW = NORD
                    Exit
      !
                  Case (4:)
      !
      !
                    XPIV0 = XPIV
                    IFIN = JHIG
      !
      !  One chooses a pivot from the 2 first values and the last one.
      !  This should ensure sufficient renewal between iterations to
      !  avoid worst case behavior effects.
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (IFIN)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (IFIN) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
      !
                    JDEB = JLOW
                    NWRK = NORD - JLOW
                    IWRK1 = IHIGT (1)
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = IWRK1
                    XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                            (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
      !
      !  One takes values <= pivot to ILOWT
      !  Again, 2 parts, one where we take care of the remaining
      !  high values because we might still need them, and the
      !  other when we know that we will have more than enough
      !  low values in the end.
      !
                    JHIG = 0
                    Do ICRS = 2, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                          If (JLOW >= NORD) Exit
                        Else
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = IHIGT (ICRS)
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                        End If
                    End Do
                End Select
      !
      !
              Case (1)
      !
      !  Only 1 value is missing in low part
      !
                  XMIN = XDONT (IHIGT(1))
                  IHIG = 1
                  Do ICRS = 2, JHIG
                    If (XDONT(IHIGT(ICRS)) < XMIN) Then
                        XMIN = XDONT (IHIGT(ICRS))
                        IHIG = ICRS
                    End If
                  End Do
      !
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (IHIG)
                  Exit
      !
      !
              Case (0)
      !
      !  Low part is exactly what we want
      !
                  Exit
      !
      !
              Case (-5:-1)
      !
      !  Only few values too many in low part
      !
                  IRNGT (1) = ILOWT (1)
                  Do ICRS = 2, NORD
                    IWRK = ILOWT (ICRS)
                    XWRK = XDONT (IWRK)
                    Do IDCR = ICRS - 1, 1, - 1
                        If (XWRK < XDONT(IRNGT(IDCR))) Then
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        Else
                          Exit
                        End If
                    End Do
                    IRNGT (IDCR+1) = IWRK
                  End Do
      !
                  XWRK1 = XDONT (IRNGT(NORD))
                  Do ICRS = NORD + 1, JLOW
                    If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                        XWRK = XDONT (ILOWT (ICRS))
                        Do IDCR = NORD - 1, 1, - 1
                          If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        End Do
                        IRNGT (IDCR+1) = ILOWT (ICRS)
                        XWRK1 = XDONT (IRNGT(NORD))
                    End If
                  End Do
      !
                  Return
      !
      !
              Case (:-6)
      !
      ! last case: too many values in low part
      !
                  IDEB = JDEB + 1
                  IMIL = (JLOW+IDEB) / 2
                  IFIN = JLOW
      !
      !  One chooses a pivot from 1st, last, and middle values
      !
                  If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                    IWRK = ILOWT (IDEB)
                    ILOWT (IDEB) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                  End If
                  If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
                    IWRK = ILOWT (IFIN)
                    ILOWT (IFIN) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                    If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                        IWRK = ILOWT (IDEB)
                        ILOWT (IDEB) = ILOWT (IMIL)
                        ILOWT (IMIL) = IWRK
                    End If
                  End If
                  If (IFIN <= 3) Exit
      !
                  XPIV = XDONT (ILOWT(1)) + REAL(NORD)/REAL(JLOW+NORD) * &
                                            (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))
                  If (JDEB > 0) Then
                    If (XPIV <= XPIV0) &
                        XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JLOW+NORD) * &
                                        (XDONT(ILOWT(IFIN))-XPIV0)
                  Else
                    IDEB = 1
                  End If
      !
      !  One takes values > XPIV to IHIGT
      !  However, we do not process the first values if we have been
      !  through the case when we did not have enough low values
      !
                  JHIG = 0
                  JLOW = JDEB
      !
                  If (XDONT(ILOWT(IFIN)) > XPIV) Then
                    ICRS = JDEB
                    Do
                      ICRS = ICRS + 1
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                          If (ICRS >= IFIN) Exit
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    If (ICRS < IFIN) Then
                        Do
                          ICRS = ICRS + 1
                          If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                              JLOW = JLOW + 1
                              ILOWT (JLOW) = ILOWT (ICRS)
                          Else
                              If (ICRS >= IFIN) Exit
                          End If
                        End Do
                    End If
                Else
                    Do ICRS = IDEB, IFIN
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                        End If
                    End Do
                  End If
      !
              End Select
      !
            End Do
      !
      !  Now, we only need to complete ranking of the 1:NORD set
      !  Assuming NORD is small, we use a simple insertion sort
      !
            IRNGT (1) = ILOWT (1)
            Do ICRS = 2, NORD
              IWRK = ILOWT (ICRS)
              XWRK = XDONT (IWRK)
              Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                    IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                    Exit
                  End If
              End Do
              IRNGT (IDCR+1) = IWRK
            End Do
          Return
      !
      !
      End Function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Ranks partially XVALT by IRNGT, up to order NORD
      !> @details This routine uses a pivoting strategy such as the one of finding
      !!          the median based on the quicksort algorithm, but we skew the pivot
      !!          choice to try to bring it to NORD as fast as possible. It uses 2
      !!          temporary arrays, one where it stores the indices of the values
      !!          smaller than the pivot, and the other for the indices of values
      !!          larger than the pivot that we might still need later on. It iterates
      !!          until it can bring the number of values in ILOWT to exactly NORD,
      !!          and then uses an insertion sort to rank this set, since it is
      !           supposedly small.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !-------------------------------------------------------------------------
      pure function I_RnkPar (XDONT, NORD) result(IRNGT)
            implicit none
            Integer(int32), Dimension(:), Intent(In)  :: XDONT !< Vector to rank
            Integer(int32), Intent(In)                :: NORD  !< Number of smallest values to rank
            Integer(int32), allocatable, Dimension(:) :: IRNGT !< @return Result

            Integer(int32) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX

            Integer(int32), Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
            Integer(int32) :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
            Integer(int32) :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
            Integer(int32) :: JLM2, JLM1, JHM2, JHM1

            NDON = SIZE (XDONT)
            allocate(IRNGT(NORD))
      !
      !    First loop is used to fill-in ILOWT, IHIGT at the same time
      !
            If (NDON < 2) Then
              If (NORD >= 1) IRNGT (1) = 1
              Return
            End If
      !
      !  One chooses a pivot, best estimate possible to put fractile near
      !  mid-point of the set of low values.
      !
            If (XDONT(2) < XDONT(1)) Then
              ILOWT (1) = 2
              IHIGT (1) = 1
            Else
              ILOWT (1) = 1
              IHIGT (1) = 2
            End If
      !
            If (NDON < 3) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              Return
            End If
      !
            If (XDONT(3) <= XDONT(IHIGT(1))) Then
              IHIGT (2) = IHIGT (1)
              If (XDONT(3) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = 3
              Else
                  IHIGT (1) = 3
              End If
            Else
              IHIGT (2) = 3
            End If
      !
            If (NDON < 4) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              Return
            End If
      !
            If (XDONT(NDON) <= XDONT(IHIGT(1))) Then
              IHIGT (3) = IHIGT (2)
              IHIGT (2) = IHIGT (1)
              If (XDONT(NDON) < XDONT(ILOWT(1))) Then
                  IHIGT (1) = ILOWT (1)
                  ILOWT (1) = NDON
              Else
                  IHIGT (1) = NDON
              End If
            Else
              if (XDONT (NDON) < XDONT (IHIGT(2))) Then
                  IHIGT (3) = IHIGT (2)
                  IHIGT (2) = NDON
              else
                  IHIGT (3) = NDON
              endif
            End If
      !
            If (NDON < 5) Then
              If (NORD >= 1) IRNGT (1) = ILOWT (1)
              If (NORD >= 2) IRNGT (2) = IHIGT (1)
              If (NORD >= 3) IRNGT (3) = IHIGT (2)
              If (NORD >= 4) IRNGT (4) = IHIGT (3)
              Return
            End If
      !
            JDEB = 0
            IDEB = JDEB + 1
            JLOW = IDEB
            JHIG = 3
            XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                        (XDONT(IHIGT(3))-XDONT(ILOWT(IDEB)))
            If (XPIV >= XDONT(IHIGT(1))) Then
              XPIV = XDONT (ILOWT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                            (XDONT(IHIGT(2))-XDONT(ILOWT(IDEB)))
              If (XPIV >= XDONT(IHIGT(1))) &
                  XPIV = XDONT (ILOWT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                                (XDONT(IHIGT(1))-XDONT(ILOWT(IDEB)))
            End If
            XPIV0 = XPIV
      !
      !  One puts values > pivot in the end and those <= pivot
      !  at the beginning. This is split in 2 cases, so that
      !  we can skip the loop test a number of times.
      !  As we are also filling in the work arrays at the same time
      !  we stop filling in the IHIGT array as soon as we have more
      !  than enough values in ILOWT.
      !
      !
            If (XDONT(NDON) > XPIV) Then
              ICRS = 3
              Do
                  ICRS = ICRS + 1
                  If (XDONT(ICRS) > XPIV) Then
                    If (ICRS >= NDON) Exit
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
      !  One restricts further processing because it is no use
      !  to store more high values
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    Else If (ICRS >= NDON) Then
                        Exit
                    End If
                  End Do
              End If
      !
      !
            Else
      !
      !  Same as above, but this is not as easy to optimize, so the
      !  DO-loop is kept
      !
              Do ICRS = 4, NDON - 1
                  If (XDONT(ICRS) > XPIV) Then
                    JHIG = JHIG + 1
                    IHIGT (JHIG) = ICRS
                  Else
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = ICRS
                    If (JLOW >= NORD) Exit
                  End If
              End Do
      !
              If (ICRS < NDON-1) Then
                  Do
                    ICRS = ICRS + 1
                    If (XDONT(ICRS) <= XPIV) Then
                        If (ICRS >= NDON) Exit
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ICRS
                    End If
                  End Do
              End If
            End If
      !
            JLM2 = 0
            JLM1 = 0
            JHM2 = 0
            JHM1 = 0
            Do
              if (JLOW == NORD) Exit
              If (JLM2 == JLOW .And. JHM2 == JHIG) Then
      !
      !   We are oscillating. Perturbate by bringing JLOW closer by one
      !   to NORD
      !
                If (NORD > JLOW) Then
                      XMIN = XDONT (IHIGT(1))
                      IHIG = 1
                      Do ICRS = 2, JHIG
                        If (XDONT(IHIGT(ICRS)) < XMIN) Then
                            XMIN = XDONT (IHIGT(ICRS))
                            IHIG = ICRS
                        End If
                      End Do
      !
                      JLOW = JLOW + 1
                      ILOWT (JLOW) = IHIGT (IHIG)
                      IHIGT (IHIG) = IHIGT (JHIG)
                      JHIG = JHIG - 1
                  Else
                      ILOW = ILOWT (JLOW)
                      XMAX = XDONT (ILOW)
                      Do ICRS = 1, JLOW
                        If (XDONT(ILOWT(ICRS)) > XMAX) Then
                            IWRK = ILOWT (ICRS)
                            XMAX = XDONT (IWRK)
                            ILOWT (ICRS) = ILOW
                            ILOW = IWRK
                        End If
                      End Do
                      JLOW = JLOW - 1
                  End If
              End If
              JLM2 = JLM1
              JLM1 = JLOW
              JHM2 = JHM1
              JHM1 = JHIG
      !
      !   We try to bring the number of values in the low values set
      !   closer to NORD.
      !
              Select Case (NORD-JLOW)
              Case (2:)
      !
      !   Not enough values in low part, at least 2 are missing
      !
                  Select Case (JHIG)
      !!!!!           CASE DEFAULT
      !!!!!              write (*,*) "Assertion failed"
      !!!!!              STOP
      !
      !   We make a special case when we have so few values in
      !   the high values set that it is bad performance to choose a pivot
      !   and apply the general algorithm.
      !
                  Case (2)
                    If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                    Else
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (2)
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = IHIGT (1)
                    End If
                    Exit
      !
                  Case (3)
      !
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (3)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (3) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
                    JHIG = 0
                    Do ICRS = JLOW + 1, NORD
                        JHIG = JHIG + 1
                        ILOWT (ICRS) = IHIGT (JHIG)
                    End Do
                    JLOW = NORD
                    Exit
      !
                  Case (4:)
      !
      !
                    XPIV0 = XPIV
                    IFIN = JHIG
      !
      !  One chooses a pivot from the 2 first values and the last one.
      !  This should ensure sufficient renewal between iterations to
      !  avoid worst case behavior effects.
      !
                    IWRK1 = IHIGT (1)
                    IWRK2 = IHIGT (2)
                    IWRK3 = IHIGT (IFIN)
                    If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                        IHIGT (1) = IWRK2
                        IHIGT (2) = IWRK1
                        IWRK2 = IWRK1
                    End If
                    If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                        IHIGT (IFIN) = IWRK2
                        IHIGT (2) = IWRK3
                        IWRK2 = IWRK3
                        If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                          IHIGT (2) = IHIGT (1)
                          IHIGT (1) = IWRK2
                        End If
                    End If
      !
                    JDEB = JLOW
                    NWRK = NORD - JLOW
                    IWRK1 = IHIGT (1)
                    JLOW = JLOW + 1
                    ILOWT (JLOW) = IWRK1
                    XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                            (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
      !
      !  One takes values <= pivot to ILOWT
      !  Again, 2 parts, one where we take care of the remaining
      !  high values because we might still need them, and the
      !  other when we know that we will have more than enough
      !  low values in the end.
      !
                    JHIG = 0
                    Do ICRS = 2, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                          If (JLOW >= NORD) Exit
                        Else
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = IHIGT (ICRS)
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = IHIGT (ICRS)
                        End If
                    End Do
                End Select
      !
      !
              Case (1)
      !
      !  Only 1 value is missing in low part
      !
                  XMIN = XDONT (IHIGT(1))
                  IHIG = 1
                  Do ICRS = 2, JHIG
                    If (XDONT(IHIGT(ICRS)) < XMIN) Then
                        XMIN = XDONT (IHIGT(ICRS))
                        IHIG = ICRS
                    End If
                  End Do
      !
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (IHIG)
                  Exit
      !
      !
              Case (0)
      !
      !  Low part is exactly what we want
      !
                  Exit
      !
      !
              Case (-5:-1)
      !
      !  Only few values too many in low part
      !
                  IRNGT (1) = ILOWT (1)
                  Do ICRS = 2, NORD
                    IWRK = ILOWT (ICRS)
                    XWRK = XDONT (IWRK)
                    Do IDCR = ICRS - 1, 1, - 1
                        If (XWRK < XDONT(IRNGT(IDCR))) Then
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        Else
                          Exit
                        End If
                    End Do
                    IRNGT (IDCR+1) = IWRK
                  End Do
      !
                  XWRK1 = XDONT (IRNGT(NORD))
                  Do ICRS = NORD + 1, JLOW
                    If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                        XWRK = XDONT (ILOWT (ICRS))
                        Do IDCR = NORD - 1, 1, - 1
                          If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                          IRNGT (IDCR+1) = IRNGT (IDCR)
                        End Do
                        IRNGT (IDCR+1) = ILOWT (ICRS)
                        XWRK1 = XDONT (IRNGT(NORD))
                    End If
                  End Do
      !
                  Return
      !
      !
              Case (:-6)
      !
      ! last case: too many values in low part
      !
                  IDEB = JDEB + 1
                  IMIL = (JLOW+IDEB) / 2
                  IFIN = JLOW
      !
      !  One chooses a pivot from 1st, last, and middle values
      !
                  If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                    IWRK = ILOWT (IDEB)
                    ILOWT (IDEB) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                  End If
                  If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
                    IWRK = ILOWT (IFIN)
                    ILOWT (IFIN) = ILOWT (IMIL)
                    ILOWT (IMIL) = IWRK
                    If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(IDEB))) Then
                        IWRK = ILOWT (IDEB)
                        ILOWT (IDEB) = ILOWT (IMIL)
                        ILOWT (IMIL) = IWRK
                    End If
                  End If
                  If (IFIN <= 3) Exit
      !
                  XPIV = XDONT (ILOWT(1)) + REAL(NORD)/REAL(JLOW+NORD) * &
                                            (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))
                  If (JDEB > 0) Then
                    If (XPIV <= XPIV0) &
                        XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JLOW+NORD) * &
                                        (XDONT(ILOWT(IFIN))-XPIV0)
                  Else
                    IDEB = 1
                  End If
      !
      !  One takes values > XPIV to IHIGT
      !  However, we do not process the first values if we have been
      !  through the case when we did not have enough low values
      !
                  JHIG = 0
                  JLOW = JDEB
      !
                  If (XDONT(ILOWT(IFIN)) > XPIV) Then
                    ICRS = JDEB
                    Do
                      ICRS = ICRS + 1
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                          If (ICRS >= IFIN) Exit
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    If (ICRS < IFIN) Then
                        Do
                          ICRS = ICRS + 1
                          If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                              JLOW = JLOW + 1
                              ILOWT (JLOW) = ILOWT (ICRS)
                          Else
                              If (ICRS >= IFIN) Exit
                          End If
                        End Do
                    End If
                Else
                    Do ICRS = IDEB, IFIN
                        If (XDONT(ILOWT(ICRS)) > XPIV) Then
                          JHIG = JHIG + 1
                          IHIGT (JHIG) = ILOWT (ICRS)
                        Else
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                          If (JLOW >= NORD) Exit
                        End If
                    End Do
      !
                    Do ICRS = ICRS + 1, IFIN
                        If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                          JLOW = JLOW + 1
                          ILOWT (JLOW) = ILOWT (ICRS)
                        End If
                    End Do
                  End If
      !
              End Select
      !
            End Do
      !
      !  Now, we only need to complete ranking of the 1:NORD set
      !  Assuming NORD is small, we use a simple insertion sort
      !
            IRNGT (1) = ILOWT (1)
            Do ICRS = 2, NORD
              IWRK = ILOWT (ICRS)
              XWRK = XDONT (IWRK)
              Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                    IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                    Exit
                  End If
              End Do
              IRNGT (IDCR+1) = IWRK
            End Do
          Return
      !
      !
      End Function

      !#########################################################################

    !###########################################################################

    ! RapKnr - array ranks (for n largest values)
    !###########################################################################

    ! UniSta - unique

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   UniSta - Stable Unique
      !> @details Removes duplicates from an array, leaving unique entries in the
      !!          order of their first appearance in the initial set.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    February, 2000
      !-------------------------------------------------------------------------
      pure function D_unista (XDONT) result(OUT)
            implicit none
            Real(real64), Dimension (:), Intent (In) :: XDONT !< Input Array
            Real(real64), allocatable, Dimension(:) :: OUT    !< @return Unique values in the array

            Real(real64), Dimension (Size(XDONT)) :: XDONTCOPY
            Integer(int32), Dimension (Size(XDONT)) :: IWRKT
            Logical, Dimension (Size(XDONT)) :: IFMPTYT
            Integer(int32) :: ICRS, NUNI

            IWRKT = UNIINV(XDONT)
            IFMPTYT = .True.
            NUNI = 0
            XDONTCOPY = XDONT
            Do ICRS = 1, Size(XDONT)
               If (IFMPTYT(IWRKT(ICRS))) Then
                  IFMPTYT(IWRKT(ICRS)) = .False.
                  NUNI = NUNI + 1
                  XDONTCOPY (NUNI) = XDONTCOPY (ICRS)
               End If
            End Do
            allocate(OUT(NUNI))
            OUT = XDONTCOPY(1:NUNI)
            Return
      End function D_unista

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   UniSta - Stable Unique
      !> @details Removes duplicates from an array, leaving unique entries in the
      !!          order of their first appearance in the initial set.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    February, 2000
      !-------------------------------------------------------------------------
      pure function R_unista (XDONT) result(OUT)
            implicit none
            Real(real32), Dimension (:), Intent (In) :: XDONT !< Input Array
            Real(real32), allocatable, Dimension (:) :: OUT   !< @return Unique values in the array

            Real(real32), Dimension (Size(XDONT)) :: XDONTCOPY
            Integer(int32), Dimension (Size(XDONT)) :: IWRKT
            Logical, Dimension (Size(XDONT)) :: IFMPTYT
            Integer(int32) :: ICRS, NUNI

            IWRKT = UNIINV(XDONT)
            IFMPTYT = .True.
            NUNI = 0
            XDONTCOPY = XDONT
            Do ICRS = 1, Size(XDONT)
               If (IFMPTYT(IWRKT(ICRS))) Then
                  IFMPTYT(IWRKT(ICRS)) = .False.
                  NUNI = NUNI + 1
                  XDONTCOPY (NUNI) = XDONTCOPY (ICRS)
               End If
            End Do
            allocate(OUT(NUNI))
            OUT = XDONTCOPY(1:NUNI)
            Return
      End function R_unista

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   UniSta - Stable Unique
      !> @details Removes duplicates from an array, leaving unique entries in the
      !!          order of their first appearance in the initial set.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    February, 2000
      !-------------------------------------------------------------------------
      pure function I_unista (XDONT) result(OUT)
            implicit none
            Integer(int32), Dimension (:), Intent (In)  :: XDONT !< Input Array
            Integer(int32), allocatable, Dimension(:) :: OUT     !< @return Unique values in the array

            Integer(int32), Dimension (Size(XDONT)) :: XDONTCOPY
            Integer(int32), Dimension (Size(XDONT)) :: IWRKT
            Logical, Dimension (Size(XDONT)) :: IFMPTYT
            Integer(int32) :: ICRS, NUNI

            IWRKT = UNIINV(XDONT)
            IFMPTYT = .True.
            NUNI = 0
            XDONTCOPY = XDONT
            Do ICRS = 1, Size(XDONT)
               If (IFMPTYT(IWRKT(ICRS))) Then
                  IFMPTYT(IWRKT(ICRS)) = .False.
                  NUNI = NUNI + 1
                  XDONTCOPY (NUNI) = XDONTCOPY (ICRS)
               End If
            End Do
            allocate(OUT(NUNI))
            OUT = XDONTCOPY(1:NUNI)
            Return
      End function I_unista

      !#########################################################################

    !###########################################################################

    ! UniInv - unique (inverse?) rank

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   UniInv - Merge-sort inverse ranking of an array, with removal of
      !!          duplicate entries
      !> @details The routine is similar to pure merge-sort ranking, but on the
      !!          last pass, it sets indices in IGOEST to the rank of the value
      !!          in the ordered set with duplicates removed. For performance reasons,
      !!          the first 2 passes are taken out of the standard loop, and use
      !!          dedicated coding.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function D_uniinv (XDONT) result(IGOEST)
            Real(real64), Dimension (:), Intent (In) :: XDONT    !< Vector to rank
            Integer(int32), allocatable, Dimension (:) :: IGOEST !< @return Result

            Real(real64) :: XTST, XDONA, XDONB

            ! Integer(int32), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
            Integer(int32), allocatable, Dimension (:) :: JWRKT, IRNGT
            Integer(int32) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
            Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

            ! NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
            NVAL = SIZE(XDONT)
            allocate(IGOEST(NVAL))
            allocate(JWRKT(NVAL))
            allocate(IRNGT(NVAL))

            Select Case (NVAL)
            Case (:0)
               Return
            Case (1)
               IGOEST (1) = 1
               Return
            Case Default
               Continue
            End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
            Do IIND = 2, NVAL, 2
               If (XDONT(IIND-1) < XDONT(IIND)) Then
                  IRNGT (IIND-1) = IIND - 1
                  IRNGT (IIND) = IIND
               Else
                  IRNGT (IIND-1) = IIND
                  IRNGT (IIND) = IIND - 1
               End If
            End Do
            If (Modulo (NVAL, 2) /= 0) Then
               IRNGT (NVAL) = NVAL
            End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
            LMTNA = 2
            LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
            Do
               If (NVAL <= 4) Exit
      !
      !   Loop on merges of A and B into C
      !
               Do IWRKD = 0, NVAL - 1, 4
                  If ((IWRKD+4) > NVAL) Then
                     If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                     If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                     If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                        IRNG2 = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                     Else
                        IRNG1 = IRNGT (IWRKD+1)
                        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNG1
                     End If
                     Exit
                  End If
      !
      !   1 2 3 4
      !
                  If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
                  If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                        IRNGT (IWRKD+3) = IRNG2
                     Else
      !   1 3 4 2
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+4) = IRNG2
                     End If
      !
      !   3 x x x
      !
                  Else
                     IRNG1 = IRNGT (IWRKD+1)
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                        IRNGT (IWRKD+2) = IRNG1
                        If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                           IRNGT (IWRKD+3) = IRNG2
                        Else
      !   3 1 4 2
                           IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                           IRNGT (IWRKD+4) = IRNG2
                        End If
                     Else
      !   3 4 1 2
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+3) = IRNG1
                        IRNGT (IWRKD+4) = IRNG2
                     End If
                  End If
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 4
               Exit
            End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
            Do
               If (2*LMTNA >= NVAL) Exit
               IWRKF = 0
               LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
               Do
                  IWRK = IWRKF
                  IWRKD = IWRKF + 1
                  JINDA = IWRKF + LMTNA
                  IWRKF = IWRKF + LMTNC
                  If (IWRKF >= NVAL) Then
                     If (JINDA >= NVAL) Exit
                     IWRKF = NVAL
                  End If
                  IINDA = 1
                  IINDB = JINDA + 1
      !
      !  One steps in the C subset, that we create in the final rank array
      !
      !  Make a copy of the rank array for the iteration
      !
                  JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                  XDONA = XDONT (JWRKT(IINDA))
                  XDONB = XDONT (IRNGT(IINDB))
      !
                  Do
                     IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                     If (XDONA > XDONB) Then
                        IRNGT (IWRK) = IRNGT (IINDB)
                        IINDB = IINDB + 1
                        If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                           IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                           Exit
                        End If
                        XDONB = XDONT (IRNGT(IINDB))
                     Else
                        IRNGT (IWRK) = JWRKT (IINDA)
                        IINDA = IINDA + 1
                        If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                        XDONA = XDONT (JWRKT(IINDA))
                     End If
      !
                  End Do
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 2 * LMTNA
            End Do
      !
      !   Last merge of A and B into C, with removal of duplicates.
      !
            IINDA = 1
            IINDB = LMTNA + 1
            NUNI = 0
      !
      !  One steps in the C subset, that we create in the final rank array
      !
            JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
            If (IINDB <= NVAL) Then
              XTST = NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
            Else
              XTST = NEARLESS (XDONT(JWRKT(1)))
            Endif
            Do IWRK = 1, NVAL
      !
      !  We still have unprocessed values in both A and B
      !
               If (IINDA <= LMTNA) Then
                  If (IINDB <= NVAL) Then
                     If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                        IRNG = IRNGT (IINDB)
                        IINDB = IINDB + 1
                     Else
                        IRNG = JWRKT (IINDA)
                        IINDA = IINDA + 1
                     End If
                  Else
      !
      !  Only A still with unprocessed values
      !
                     IRNG = JWRKT (IINDA)
                     IINDA = IINDA + 1
                  End If
               Else
      !
      !  Only B still with unprocessed values
      !
                  IRNG = IRNGT (IWRK)
               End If
               If (XDONT(IRNG) > XTST) Then
                  XTST = XDONT (IRNG)
                  NUNI = NUNI + 1
               End If
               IGOEST (IRNG) = NUNI
      !
            End Do
      !
            Return
      !
      end function

      !-------------------------------------------------------------------------
      !> @brief   UniInv - Merge-sort inverse ranking of an array, with removal of
      !!          duplicate entries
      !> @details The routine is similar to pure merge-sort ranking, but on the
      !!          last pass, it sets indices in IGOEST to the rank of the value
      !!          in the ordered set with duplicates removed. For performance reasons,
      !!          the first 2 passes are taken out of the standard loop, and use
      !!          dedicated coding.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function R_uniinv (XDONT) result(IGOEST)
            Real(real32), Dimension (:), Intent (In) :: XDONT    !< Vector to rank
            Integer(int32), allocatable, Dimension (:) :: IGOEST !< @return Result

            Real(real32) :: XTST, XDONA, XDONB

            ! Integer(int32), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
            Integer(int32), allocatable, Dimension (:) :: JWRKT, IRNGT
            Integer(int32) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
            Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

            ! NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
            NVAL = SIZE(XDONT)
            allocate(IGOEST(NVAL))
            allocate(JWRKT(NVAL))
            allocate(IRNGT(NVAL))

            Select Case (NVAL)
            Case (:0)
               Return
            Case (1)
               IGOEST (1) = 1
               Return
            Case Default
               Continue
            End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
            Do IIND = 2, NVAL, 2
               If (XDONT(IIND-1) < XDONT(IIND)) Then
                  IRNGT (IIND-1) = IIND - 1
                  IRNGT (IIND) = IIND
               Else
                  IRNGT (IIND-1) = IIND
                  IRNGT (IIND) = IIND - 1
               End If
            End Do
            If (Modulo (NVAL, 2) /= 0) Then
               IRNGT (NVAL) = NVAL
            End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
            LMTNA = 2
            LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
            Do
               If (NVAL <= 4) Exit
      !
      !   Loop on merges of A and B into C
      !
               Do IWRKD = 0, NVAL - 1, 4
                  If ((IWRKD+4) > NVAL) Then
                     If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                     If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                     If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                        IRNG2 = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                     Else
                        IRNG1 = IRNGT (IWRKD+1)
                        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNG1
                     End If
                     Exit
                  End If
      !
      !   1 2 3 4
      !
                  If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
                  If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                        IRNGT (IWRKD+3) = IRNG2
                     Else
      !   1 3 4 2
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+4) = IRNG2
                     End If
      !
      !   3 x x x
      !
                  Else
                     IRNG1 = IRNGT (IWRKD+1)
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                        IRNGT (IWRKD+2) = IRNG1
                        If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                           IRNGT (IWRKD+3) = IRNG2
                        Else
      !   3 1 4 2
                           IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                           IRNGT (IWRKD+4) = IRNG2
                        End If
                     Else
      !   3 4 1 2
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+3) = IRNG1
                        IRNGT (IWRKD+4) = IRNG2
                     End If
                  End If
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 4
               Exit
            End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
            Do
               If (2*LMTNA >= NVAL) Exit
               IWRKF = 0
               LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
               Do
                  IWRK = IWRKF
                  IWRKD = IWRKF + 1
                  JINDA = IWRKF + LMTNA
                  IWRKF = IWRKF + LMTNC
                  If (IWRKF >= NVAL) Then
                     If (JINDA >= NVAL) Exit
                     IWRKF = NVAL
                  End If
                  IINDA = 1
                  IINDB = JINDA + 1
      !
      !  One steps in the C subset, that we create in the final rank array
      !
      !  Make a copy of the rank array for the iteration
      !
                  JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                  XDONA = XDONT (JWRKT(IINDA))
                  XDONB = XDONT (IRNGT(IINDB))
      !
                  Do
                     IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                     If (XDONA > XDONB) Then
                        IRNGT (IWRK) = IRNGT (IINDB)
                        IINDB = IINDB + 1
                        If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                           IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                           Exit
                        End If
                        XDONB = XDONT (IRNGT(IINDB))
                     Else
                        IRNGT (IWRK) = JWRKT (IINDA)
                        IINDA = IINDA + 1
                        If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                        XDONA = XDONT (JWRKT(IINDA))
                     End If
      !
                  End Do
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 2 * LMTNA
            End Do
      !
      !   Last merge of A and B into C, with removal of duplicates.
      !
            IINDA = 1
            IINDB = LMTNA + 1
            NUNI = 0
      !
      !  One steps in the C subset, that we create in the final rank array
      !
            JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
            If (IINDB <= NVAL) Then
              XTST = NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
            Else
              XTST = NEARLESS (XDONT(JWRKT(1)))
            Endif
            Do IWRK = 1, NVAL
      !
      !  We still have unprocessed values in both A and B
      !
               If (IINDA <= LMTNA) Then
                  If (IINDB <= NVAL) Then
                     If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                        IRNG = IRNGT (IINDB)
                        IINDB = IINDB + 1
                     Else
                        IRNG = JWRKT (IINDA)
                        IINDA = IINDA + 1
                     End If
                  Else
      !
      !  Only A still with unprocessed values
      !
                     IRNG = JWRKT (IINDA)
                     IINDA = IINDA + 1
                  End If
               Else
      !
      !  Only B still with unprocessed values
      !
                  IRNG = IRNGT (IWRK)
               End If
               If (XDONT(IRNG) > XTST) Then
                  XTST = XDONT (IRNG)
                  NUNI = NUNI + 1
               End If
               IGOEST (IRNG) = NUNI
      !
            End Do
      !
            Return
      !
      end function

      !-------------------------------------------------------------------------
      !> @brief   UniInv - Merge-sort inverse ranking of an array, with removal of
      !!          duplicate entries
      !> @details The routine is similar to pure merge-sort ranking, but on the
      !!          last pass, it sets indices in IGOEST to the rank of the value
      !!          in the ordered set with duplicates removed. For performance reasons,
      !!          the first 2 passes are taken out of the standard loop, and use
      !!          dedicated coding.
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank, modified by
      !!          Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function I_uniinv (XDONT) result(IGOEST)
            Integer(int32), Dimension (:), Intent (In)  :: XDONT !< Vector to rank
            Integer(int32), allocatable, Dimension (:) :: IGOEST !< @return Result

            Integer(int32) :: XTST, XDONA, XDONB

            ! Integer(int32), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
            Integer(int32), allocatable, Dimension (:) :: JWRKT, IRNGT
            Integer(int32) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
            Integer(int32) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

            ! NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
            NVAL = SIZE(XDONT)
            allocate(IGOEST(NVAL))
            allocate(JWRKT(NVAL))
            allocate(IRNGT(NVAL))

            Select Case (NVAL)
            Case (:0)
               Return
            Case (1)
               IGOEST (1) = 1
               Return
            Case Default
               Continue
            End Select
      !
      !  Fill-in the index array, creating ordered couples
      !
            Do IIND = 2, NVAL, 2
               If (XDONT(IIND-1) < XDONT(IIND)) Then
                  IRNGT (IIND-1) = IIND - 1
                  IRNGT (IIND) = IIND
               Else
                  IRNGT (IIND-1) = IIND
                  IRNGT (IIND) = IIND - 1
               End If
            End Do
            If (Modulo (NVAL, 2) /= 0) Then
               IRNGT (NVAL) = NVAL
            End If
      !
      !  We will now have ordered subsets A - B - A - B - ...
      !  and merge A and B couples into     C   -   C   - ...
      !
            LMTNA = 2
            LMTNC = 4
      !
      !  First iteration. The length of the ordered subsets goes from 2 to 4
      !
            Do
               If (NVAL <= 4) Exit
      !
      !   Loop on merges of A and B into C
      !
               Do IWRKD = 0, NVAL - 1, 4
                  If ((IWRKD+4) > NVAL) Then
                     If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
                     If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
                     If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                        IRNG2 = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
                     Else
                        IRNG1 = IRNGT (IWRKD+1)
                        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                        IRNGT (IWRKD+2) = IRNG1
                     End If
                     Exit
                  End If
      !
      !   1 2 3 4
      !
                  If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !
      !   1 3 x x
      !
                  If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   1 3 2 4
                        IRNGT (IWRKD+3) = IRNG2
                     Else
      !   1 3 4 2
                        IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+4) = IRNG2
                     End If
      !
      !   3 x x x
      !
                  Else
                     IRNG1 = IRNGT (IWRKD+1)
                     IRNG2 = IRNGT (IWRKD+2)
                     IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                     If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                        IRNGT (IWRKD+2) = IRNG1
                        If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
      !   3 1 2 4
                           IRNGT (IWRKD+3) = IRNG2
                        Else
      !   3 1 4 2
                           IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                           IRNGT (IWRKD+4) = IRNG2
                        End If
                     Else
      !   3 4 1 2
                        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                        IRNGT (IWRKD+3) = IRNG1
                        IRNGT (IWRKD+4) = IRNG2
                     End If
                  End If
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 4
               Exit
            End Do
      !
      !  Iteration loop. Each time, the length of the ordered subsets
      !  is doubled.
      !
            Do
               If (2*LMTNA >= NVAL) Exit
               IWRKF = 0
               LMTNC = 2 * LMTNC
      !
      !   Loop on merges of A and B into C
      !
               Do
                  IWRK = IWRKF
                  IWRKD = IWRKF + 1
                  JINDA = IWRKF + LMTNA
                  IWRKF = IWRKF + LMTNC
                  If (IWRKF >= NVAL) Then
                     If (JINDA >= NVAL) Exit
                     IWRKF = NVAL
                  End If
                  IINDA = 1
                  IINDB = JINDA + 1
      !
      !  One steps in the C subset, that we create in the final rank array
      !
      !  Make a copy of the rank array for the iteration
      !
                  JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                  XDONA = XDONT (JWRKT(IINDA))
                  XDONB = XDONT (IRNGT(IINDB))
      !
                  Do
                     IWRK = IWRK + 1
      !
      !  We still have unprocessed values in both A and B
      !
                     If (XDONA > XDONB) Then
                        IRNGT (IWRK) = IRNGT (IINDB)
                        IINDB = IINDB + 1
                        If (IINDB > IWRKF) Then
      !  Only A still with unprocessed values
                           IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                           Exit
                        End If
                        XDONB = XDONT (IRNGT(IINDB))
                     Else
                        IRNGT (IWRK) = JWRKT (IINDA)
                        IINDA = IINDA + 1
                        If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                        XDONA = XDONT (JWRKT(IINDA))
                     End If
      !
                  End Do
               End Do
      !
      !  The Cs become As and Bs
      !
               LMTNA = 2 * LMTNA
            End Do
      !
      !   Last merge of A and B into C, with removal of duplicates.
      !
            IINDA = 1
            IINDB = LMTNA + 1
            NUNI = 0
      !
      !  One steps in the C subset, that we create in the final rank array
      !
            JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
            If (IINDB <= NVAL) Then
              XTST = NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
            Else
              XTST = NEARLESS (XDONT(JWRKT(1)))
            Endif
            Do IWRK = 1, NVAL
      !
      !  We still have unprocessed values in both A and B
      !
               If (IINDA <= LMTNA) Then
                  If (IINDB <= NVAL) Then
                     If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                        IRNG = IRNGT (IINDB)
                        IINDB = IINDB + 1
                     Else
                        IRNG = JWRKT (IINDA)
                        IINDA = IINDA + 1
                     End If
                  Else
      !
      !  Only A still with unprocessed values
      !
                     IRNG = JWRKT (IINDA)
                     IINDA = IINDA + 1
                  End If
               Else
      !
      !  Only B still with unprocessed values
      !
                  IRNG = IRNGT (IWRK)
               End If
               If (XDONT(IRNG) > XTST) Then
                  XTST = XDONT (IRNG)
                  NUNI = NUNI + 1
               End If
               IGOEST (IRNG) = NUNI
      !
            End Do
      !
            Return
      !
      end function

      !-------------------------------------------------------------------------
      !> @brief   Nearest value less than given value
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function D_nearless (XVAL) result (D_nl)
        implicit none
        Real(real64), intent (In) :: XVAL
        Real(real64) :: D_nl

        D_nl = nearest (XVAL, -1.0d0)
        return
      end function

      !-------------------------------------------------------------------------
      !> @brief   Nearest value less than given value
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function R_nearless (XVAL) result (R_nl)
        implicit none
        Real(real32), intent (In) :: XVAL
        Real(real32) :: R_nl

        R_nl = nearest (XVAL, -1.0)
        return
      end function

      !-------------------------------------------------------------------------
      !> @brief   Nearest value less than given value
      !> @author  Michel Olagnon, http://www.fortran-2000.com/rank
      !> @date    2010
      !-------------------------------------------------------------------------
      pure function I_nearless (XVAL) result (I_nl)
        implicit none
        Integer(int32), intent (In) :: XVAL
        Integer(int32) :: I_nl

        I_nl = XVAL - 1
        return
      end function

      !#########################################################################

    !###########################################################################
end module

!###############################################################################
