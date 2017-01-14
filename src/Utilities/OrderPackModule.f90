
!###############################################################################

module OrderPackModule
  ! From http://www.fortran-2000.com/rank (2016-02-15)
  ! Modified to work with ISO_Fortran_Env types and converted to
  ! functions where possible

  use ISO_Fortran_Env

  implicit none

  private
  public :: MrgRnk, UniSta, UniInv

  interface MrgRnk
    module procedure D_MrgRnk,R_MrgRnk,I_MrgRnk
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

    ! MrgRnk - ranks array

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
