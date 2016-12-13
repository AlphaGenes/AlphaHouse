module UtilitySubroutines
    implicit none

! Module contains subroutines and functions that are completely independent and can be completely decoupled from other modules.

  !###########################################################################################################################################################
  !###########################################################################################################################################################

!  ____  ____   __   ____  _  _  ____ 
! (  _ \(  __) / _\ (    \( \/ )(  __)
!  )   / ) _) /    \ ) D (/ \/ \ ) _) 
! (__\_)(____)\_/\_/(____/\_)(_/(____)

! THIS IS LEGACY - Please use Alpha*Mod files instead!


  !###########################################################################################################################################################
  !###########################################################################################################################################################


    contains


  !###########################################################################################################################################################
  
  ! same as DescStatSingle in AlphaStatMod.f90
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
  
!#######################################################################

! same as DescStatDouble in AlphaStatMod.f90
  SUBROUTINE momentR4(DATA,n,ave,adev,sdev,var,skew,curt)
    use iso_fortran_env
    IMPLICIT NONE
    INTEGER(kind=int64) :: n
    real(kind=real32) :: adev,ave,curt,sdev,skew,var,DATA(n)
    INTEGER(kind=int64) ::j
    real(kind=real32) :: p,s,ep
    IF (n.le.1) STOP 140001
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
  END SUBROUTINE momentR4


!##############################################################################################################

! pearson correlation that looks similar to cor in alphastatmod
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

! more modern version exists in AlphaStatMod
  subroutine CholDc(a,n,np,p,NotPosDefin)    !(DOUBLE PRECISION MODIFIED)
    !SG modified 3/10/15 to return error where user specified matrix is not positive definite
    !DB modified 4/8/2016 to not use Sum as that is a intrinsic procedure
    !DW modified 20/9/2016 to take int64

    use iso_fortran_env, only : int64
    INTEGER(kind=int64) :: n,np
    DOUBLE PRECISION  :: a(np,np),p(n)
    INTEGER(kind=int64) :: i,j,k
    DOUBLE PRECISION  :: sum1
    logical, intent(inout) :: NotPosDefin

    NotPosDefin = .False.

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

  !  (C) Copr. 1986-92 Numerical Recipes Software 6
  ! converted to fortran 90 syntax for consistency

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !---------------------------------------------------------------------------
    !> @brief   Removes whitespace from a string
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    October 18, 2016
    !---------------------------------------------------------------------------
    subroutine ParseStringWindows(str,newstr)
        implicit none

        integer :: k
        character(len=*), intent(in) :: str
        character(len=:), allocatable, intent(inout) :: newstr

        character(len=512) :: dummyStr

        dummyStr = TLC(str)
        k=1
        newstr=""
        do 
            if (dummyStr(k:k) /= " ") then
                newstr = dummyStr(1:k)
                k=k+1
            else
                exit
            endif
        enddo

    end subroutine ParseStringWindows


    !---------------------------------------------------------------------------
    !> @brief   splits string into initial and second part, where second part is another array
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    October 18, 2016
    !---------------------------------------------------------------------------
    subroutine parseLine(line,first, second)

     implicit none

    integer :: lenin,i
    integer :: sCount1, sCount2, fCount ! SCount1 starts from 0 to avoid extra allocation
    integer, parameter :: numberAfterComma = 10000 ! Can be max 10000 params after first comma
    character(len=*), intent(in) :: line
    character :: c
    character(len=300) :: first
    character (len=300),dimension(:), allocatable,intent(out) :: second
    character (len=300), dimension(:), allocatable :: tmp
    logical :: useSecond    
    first = " "
    useSecond = .false.
    sCount1 = 0
    sCount2 = 1
    fCount = 1
    if (allocated(second)) deallocate(second)
    allocate(second(numberAfterComma)) ! Allocate the second array
    lenin=len(line)
    lenin=len_trim(line(1:lenin))  ! Trim string as no point dealing with trailing whitespace   
    do i=1,lenin
        c=line(i:i)
        if(.not. (ichar(c) == 9 .or. c == " " .or. c=="")) then  !if c is not tab or whitespace
            if (c == ',') then
                sCount1 = sCount1 + 1
                scount2 = 1 ! reset scount 2
                useSecond = .true.
                second(sCount1) = " "
            else if (useSecond) then
                second(sCount1)(sCount2:sCount2) = c
                scount2 = sCount2 + 1
            else
                first(fCount:fCount) = c
                fCount = fCount + 1
            endif
        else
            !We know that it is either a tab or whitespace, so why add them
            continue
        endif
    enddo
    if (sCount1 > 0) then
        allocate(tmp(sCount1)) !Allocate temp to count after
        tmp(1:sCount1) = second(1:sCount1) ! make tmp contain values of second
        call move_alloc(tmp, second)
    else
        deallocate(second) ! Deallocate second if nothing has come afterwards
    endif
    end subroutine parseLine


    !---------------------------------------------------------------------------
    !> @brief   takes in a string and returns the lowercase version
    !> @author  John Hickey, john.hickey@roslin.ed.ac.uk
    !> @date    October 18, 2016
    !---------------------------------------------------------------------------
    ! same as function toLower in AlphaHouseMod
    function TLC(str)
        character(*), intent(in) :: str
        character(len=512) :: TLC
        integer :: i
        TLC = trim(str)
        do i = 1, len(TLC)
            select case(TLC(i:i))
                case("A":"Z")
                    TLC(i:i) = achar(iachar(TLC(i:i))+32)
            end select
        end do
        return
    end function TLC


    ! function ran1(idum) result(Res)
    !   use iso_fortran_env
    !   use IntelRNGMod
    !     implicit none
    !     integer, intent(in), optional:: idum ! this is the seed
    !     real(kind=real64) :: res
    !     real(kind=real64):: temp(1)
    !         if (present(idum)) then
    !             call IntitialiseIntelRNG(idum)
    !         else
    !             call IntitialiseIntelRNG()
    !         endif
    !         temp = SampleIntelUniformD()
    !         res = temp(1)
    !         call UnintitialiseIntelRNG
    !     return
    ! end function ran1


! function corresponds to SampleIntelUniformD in intelRNGMOD
    function ran1(idum)
    use iso_fortran_env
        implicit none
        integer(kind=int32) idum,ia,im,iq,ir,ntab,ndiv
        double precision ran1,am,eps,rnmx
        parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
        integer j,k,iv(ntab),iy
        save iv,iy
        data iv /ntab*0/, iy /0/
        if (idum.le.0.or.iy.eq.0) then
            idum=max(-idum,1)
            do 11 j=ntab+8,1,-1
                k=idum/iq
                idum=ia*(idum-k*iq)-ir*k
                if (idum.lt.0) idum=idum+im
                if (j.le.ntab) iv(j)=idum

    11      continue
            iy=iv(1)
        end if
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        j=1+iy/ndiv
        iy=iv(j)
        iv(j)=idum
        ran1=min(am*iy,rnmx)
        return
    end function ran1


 !---------------------------------------------------------------------------
    !> @brief   Subroutine takes in number n, and returns random ordering of longs based on idum seed
    !---------------------------------------------------------------------------
   
    ! function corresponds to SampleIntelUniformD in intelRNGMOD
    subroutine RandomOrder(order,n,idum)
        use iso_fortran_env
        implicit none

        !     Generate a random ordering of the integers 1 ... n.

        integer(kind=int64), INTENT(IN)  :: n
        integer(kind=int64), INTENT(OUT) :: order(n)
        integer, INTENT(IN)  :: idum
    !    double precision ran1

        !     Local variables

        integer :: i, j, k
        double precision    :: wk

        do i = 1, n
            order(i) = i
        end do

        !     Starting at the end, swap the current last indicator with one
        !     randomly chosen from those preceeding it.

        do i = n, 2, -1
            wk=ran1(idum)
            j = 1 + i * wk
            if (j < i) then
                k = order(i)
                order(i) = order(j)
                order(j) = k
            end if
        end do

        RETURN

    end subroutine RandomOrder


    ! function corresponds to SampleIntelUniformS in intelRNGMOD
    subroutine RandomOrder4byte(order,n,idum)
        use iso_fortran_env
        !modified for 4byte integers
        implicit none

        !     Generate a random ordering of the integers 1 ... n.

        integer(kind=int32), INTENT(IN)  :: n
        integer(kind=int32), INTENT(OUT) :: order(n)
        integer :: idum
    !    double precision ran1

        !     Local variables

        integer :: i, j, k
        double precision    :: wk

        do i = 1, n
            order(i) = i
        end do

        !     Starting at the end, swap the current last indicator with one
        !     randomly chosen from those preceeding it.

        do i = n, 2, -1
            wk=ran1(idum)
            j = 1 + i * wk
            if (j < i) then
                k = order(i)
                order(i) = order(j)
                order(j) = k
            end if
        end do

        RETURN

    end subroutine RandomOrder4byte



   ! function corresponds to SampleIntelGammaD in intelRNGMOD
    function random_gamma(idum,s, b, first) result(fn_val)
        implicit none

        ! adapted from fortran 77 code from the book:
        !     dagpunar, j. 'principles of random variate generation'
        !     clarendon press, oxford, 1988.   isbn 0-19-852202-9

        !     n.b. this version is in `double precision' and includes scaling

        !     function generates a random gamma variate.
        !     calls either random_gamma1 (s > 1.0)
        !     or random_exponential (s = 1.0)
        !     or random_gamma2 (s < 1.0).

        !     s = shape parameter of distribution (0 < real).
        !     b = scale parameter

        !implicit none
        integer, parameter  :: dp = selected_real_kind(12, 60)
        integer :: idum

        double precision, intent(in)  :: s, b
        logical, intent(in)    :: first
        double precision              :: fn_val

        ! local parameters
        double precision, parameter  :: one = 1.0_dp, zero = 0.0_dp

        if (s <= zero) then
            write(*, *) 'shape parameter value must be positive'
            stop
        end if

        if (s >= one) then
            fn_val = random_gamma1(s, first)
        else if (s < one) then
            fn_val = random_gamma2(s, first)
        end if

        ! now scale the random variable
        fn_val = b * fn_val
        return

    contains


        function random_gamma1(s, first) result(fn_val)
            implicit none

            ! adapted from fortran 77 code from the book:
            !     dagpunar, j. 'principles of random variate generation'
            !     clarendon press, oxford, 1988.   isbn 0-19-852202-9

            ! function generates a random variate in [0,infinity) from
            ! a gamma distribution with density proportional to gamma**(s-1)*exp(-gamma),
            ! based upon best's t distribution method

            !     s = shape parameter of distribution
            !          (1.0 < real)

            double precision, intent(in)  :: s
            logical, intent(in)    :: first
            double precision              :: fn_val

            !     local variables
            double precision             :: d, r, g, f, x
            double precision, save       :: b, h
            double precision, parameter  :: sixty4 = 64.0_dp, three = 3.0_dp, pt75 = 0.75_dp,  &
                two = 2.0_dp, half = 0.5_dp
     !       double precision :: ran1

            if (s <= one) then
                write(*, *) 'impermissible shape parameter value'
                stop
            end if

            if (first) then                        ! initialization, if necessary
                b = s - one
                h = sqrt(three*s - pt75)
            end if

            do
                r=ran1(idum)
                g = r - r*r
                if (g <= zero) cycle
                f = (r - half)*h/sqrt(g)
                x = b + f
                if (x <= zero) cycle
                r=ran1(idum)
                d = sixty4*g*(r*g)**2
                if (d <= zero) exit
                if (d*x < x - two*f*f) exit
                if (log(d) < two*(b*log(x/b) - f)) exit
            end do
            fn_val = x

            return
        end function random_gamma1

        function random_gamma2(s, first) result(fn_val)
            implicit none

            ! adapted from fortran 77 code from the book:
            !     dagpunar, j. 'principles of random variate generation'
            !     clarendon press, oxford, 1988.   isbn 0-19-852202-9

            ! function generates a random variate in [0,infinity) from
            ! a gamma distribution with density proportional to
            ! gamma2**(s-1) * exp(-gamma2),
            ! using a switching method.

            !    s = shape parameter of distribution
            !          (real < 1.0)

            double precision, intent(in)  :: s
            logical, intent(in)    :: first
            double precision              :: fn_val

            !     local variables
            double precision            :: r, x, w
            double precision, save       :: a, p, c, uf, vr, d
            double precision, parameter  :: vsmall = epsilon(one)
    !        double precision :: ran1

            if (s <= zero .or. s >= one) then
                write(*, *) 'shape parameter value outside permitted range'
                stop
            end if

            if (first) then                        ! initialization, if necessary
                a = one - s
                p = a/(a + s*exp(-a))
                if (s < vsmall) then
                    write(*, *) 'shape parameter value too small'
                    print*, s
                    stop
                end if
                c = one/s
                uf = p*(vsmall/a)**s
                vr = one - vsmall
                d = a*log(a)
            end if

            do
                r=ran1(idum)
                if (r >= vr) then
                    cycle
                else if (r > p) then
                    x = a - log((one - r)/(one - p))
                    w = a*log(x)-d
                else if (r > uf) then
                    x = a*(r/p)**c
                    w = x
                else
                    fn_val = zero
                    return
                end if

                r=ran1(idum)
                if (one-r <= w .and. r > zero) then
                    if (r*(w + one) >= one) cycle
                    if (-log(r) <= w) cycle
                end if
                exit
            end do

            fn_val = x
            return

        end function random_gamma2

    end function random_gamma



    function gasdev(idum)
        implicit none
        !c uses ran1
        !returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
        !as the source of uniform deviates.

        integer idum
        double precision :: gasdev
        integer iset
        double precision fac,gset,rsq,v1,v2
        save iset,gset
        data iset/0/
        if (idum.lt.0) iset=0
        if (iset.eq.0) then
    1       v1=2.*ran1(idum)-1.
            v2=2.*ran1(idum)-1.
            rsq=v1**2+v2**2
            if(rsq.ge.1..or.rsq.eq.0.)goto 1
            fac=sqrt(-2.*log(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
        else
            gasdev=gset
            iset=0
        endif
        return
    end



! seems to correspond to gamma in intelRNGMOD
   SUBROUTINE cgamma(mo, z, w)
       !-----------------------------------------------------------------------

       !        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

       !                        ---------------

       !     MO IS AN INTEGER, Z A COMPLEX ARGUMENT, AND W A COMPLEX VARIABLE.

       !                 W = GAMMA(Z)       IF MO = 0
       !                 W = LN(GAMMA(Z))   OTHERWISE

       !-----------------------------------------------------------------------
       !     WRITTEN BY ALFRED H. MORRIS, JR.
       !        NAVAL SURFACE WARFARE CENTER
       !        DAHLGREN, VIRGINIA

       !     This version, in a subset of Fortran 90, prepared by
       !     Alan.Miller @ vic.cmis.csiro.au
       !     http://www.ozemail.com.au/~milleraj

       !     This version is accurate to within 5 in the 14th significant
       !     decimal digit.
       !-----------------------------------------------------------------------

       IMPLICIT NONE
       INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

       INTEGER, INTENT(IN)       :: mo
       COMPLEX (dp), INTENT(IN)  :: z
       COMPLEX (dp), INTENT(OUT) :: w

       ! Local variables
       COMPLEX (dp) :: eta, eta2, sum
       REAL (dp), PARAMETER :: c0(12) = (/ .833333333333333E-01_dp,  &
           -.277777777777778E-02_dp, .793650793650794E-03_dp,  &
           -.595238095238095E-03_dp, .841750841750842E-03_dp,  &
           -.191752691752692E-02_dp, .641025641025641E-02_dp,  &
           -.295506535947712E-01_dp, .179644372368831_dp,      &
           -1.39243221690590_dp,     13.4028640441684_dp,      &
           -156.848284626002_dp /), pi = 3.14159265358979_dp,  &
           pi2  = 6.28318530717959_dp, alpi = 1.14472988584940_dp,  &
           hl2p = .918938533204673_dp, half = 0.5_dp
       REAL (dp)  :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
           s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
       INTEGER    :: j, k, l, m, max, n, nm1
       !---------------------------
       !     ALPI = LOG(PI)
       !     HL2P = 0.5 * LOG(2*PI)
       !---------------------------

       !     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
       !            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
       !            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
       !            SUCH THAT 1.0 + EPS > 1.0.

       !                      MAX = IPMPAR(3)
       max = HUGE(3)
       eps = EPSILON(1.0_dp)

       !---------------------------
       x = REAL(z, KIND=dp)
       y = AIMAG(z)
       IF (x < 0.0_dp) THEN
           !-----------------------------------------------------------------------
           !            CASE WHEN THE REAL PART OF Z IS NEGATIVE
           !-----------------------------------------------------------------------
           y = ABS(y)
           t = -pi * y
           et = EXP(t)
           e2t = et * et

           !     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

           a1 = half * (1.0_dp + e2t)
           t2 = t + t
           IF (t2 >= -0.15_dp) THEN
               a2 = -half * rexp(t2)
           ELSE
               a2 = half * (half + (half - e2t))
           END IF

           !     COMPUTE SIN(PI*X) AND COS(PI*X)

           IF (ABS(x) >= MIN(REAL(MAX), 1.0_dp/eps)) GO TO 70
           k = ABS(x)
           u = x + k
           k = MOD(k,2)
           IF (u <= -half) THEN
               u = half + (half + u)
               k = k + 1
           END IF
           u = pi * u
           sn = SIN(u)
           cn = COS(u)
           IF (k == 1) THEN
               sn = -sn
               cn = -cn
           END IF

           !     SET  H1 + H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))

           a1 = sn * a1
           a2 = cn * a2
           a = a1 * a1 + a2 * a2
           IF (a == 0.0_dp) GO TO 70
           IF (mo == 0) THEN

               h1 = a1 / a
               h2 = -a2 / a
               c = pi * et
               h1 = c * h1
               h2 = c * h2
           ELSE

               h1 = (alpi+t) - half * LOG(a)
               h2 = -ATAN2(a2,a1)
           END IF
           IF (AIMAG(z) >= 0.0_dp) THEN
               x = 1.0_dp - x
               y = -y
           ELSE
               h2 = -h2
               x = 1.0_dp - x
           END IF
       END IF
       !-----------------------------------------------------------------------
       !           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
       !-----------------------------------------------------------------------
       w1 = 0.0_dp
       w2 = 0.0_dp
       n = 0
       t = x
       y2 = y * y
       a = t * t + y2
       cut = 36.0_dp
       IF (eps > 1.e-8_dp) cut = 16.0_dp
       IF (a < cut) THEN
           IF (a == 0.0_dp) GO TO 70
10         n = n + 1
           t = t + 1.0_dp
           a = t * t + y2
           IF (a < cut) GO TO 10

           !     LET S1 + S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)

           u1 = (x*t+y2) / a
           u2 = y / a
           s1 = u1
           s2 = n * u2
           IF (n >= 2) THEN
               u = t / a
               nm1 = n - 1
               DO j = 1, nm1
                   v1 = u1 + j * u
                   v2 = (n-j) * u2
                   c = s1 * v1 - s2 * v2
                   d = s1 * v2 + s2 * v1
                   s1 = c
                   s2 = d
               END DO
           END IF

           !     SET  W1 + W2*I = LOG(S1 + S2*I)  WHEN MO IS NONZERO

           s = s1 * s1 + s2 * s2
           IF (mo /= 0) THEN
               w1 = half * LOG(s)
               w2 = ATAN2(s2,s1)
           END IF
       END IF

       !     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

       t1 = half * LOG(a) - 1.0_dp
       t2 = ATAN2(y,t)
       u = x - half
       v1 = (u*t1-half) - y * t2
       v2 = u * t2 + y * t1

       !     LET A1 + A2*I BE THE ASYMPTOTIC SUM

       eta = CMPLX(t/a, -y/a, KIND=dp)
       eta2 = eta * eta
       m = 12
       IF (a >= 289.0_dp) m = 6
       IF (eps > 1.e-8) m = m / 2
       sum = CMPLX(c0(m), 0.0_dp, KIND=dp)
       l = m
       DO j = 2, m
           l = l - 1
           sum = CMPLX(c0(l), 0.0_dp, KIND=dp) + sum * eta2
       END DO
       sum = sum * eta
       a1 = REAL(sum, KIND=dp)
       a2 = AIMAG(sum)
       !-----------------------------------------------------------------------
       !                 GATHERING TOGETHER THE RESULTS
       !-----------------------------------------------------------------------
       w1 = (((a1 + hl2p) - w1) + v1) - n
       w2 = (a2 - w2) + v2
       IF (REAL(z, KIND=dp) < 0.0_dp) GO TO 50
       IF (mo == 0) THEN

           !     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

           a = EXP(w1)
           w1 = a * COS(w2)
           w2 = a * SIN(w2)
           IF (n == 0) GO TO 60
           c = (s1*w1 + s2*w2) / s
           d = (s1*w2 - s2*w1) / s
           w1 = c
           w2 = d
           GO TO 60
       END IF

       !     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
       !     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

40     IF (w2 <= pi) THEN
           k = half - w2 / pi2
           w2 = w2 + pi2 * k
           GO TO 60
       END IF
       k = w2 / pi2 - half
       w2 = w2 - pi2 * REAL(k+1)
       IF (w2 <= -pi) w2 = pi
       GO TO 60

       !     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

50     IF (mo /= 0) THEN
           w1 = h1 - w1
           w2 = h2 - w2
           GO TO 40
       END IF

       !     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO = 0

       a = EXP(-w1)
       t1 = a * COS(-w2)
       t2 = a * SIN(-w2)
       w1 = h1 * t1 - h2 * t2
       w2 = h1 * t2 + h2 * t1
       IF (n /= 0) THEN
           c = w1 * s1 - w2 * s2
           d = w1 * s2 + w2 * s1
           w1 = c
           w2 = d
       END IF

       !     TERMINATION

60     w = CMPLX(w1, w2, KIND=dp)
       RETURN
       !-----------------------------------------------------------------------
       !             THE REQUESTED VALUE CANNOT BE COMPUTED
       !-----------------------------------------------------------------------
70     w = (0.0_dp, 0.0_dp)
       RETURN

       CONTAINS

       FUNCTION rexp(x) RESULT(fn_val)
           !-----------------------------------------------------------------------
           !            EVALUATION OF THE FUNCTION EXP(X) - 1
           !-----------------------------------------------------------------------
           REAL (dp), INTENT(IN) :: x
           REAL (dp)             :: fn_val

           ! Local variables
           REAL (dp), PARAMETER           :: p1 =  .914041914819518E-09_dp,  &
               p2 = .238082361044469E-01_dp, q1 = -.499999999085958_dp,      &
               q2 = .107141568980644_dp,     q3 = -.119041179760821E-01_dp,  &
               q4 = .595130811860248E-03_dp
           REAL (dp) :: e
           !-----------------------
           IF (ABS(x) <= 0.15_dp) THEN
               fn_val = x * (((p2*x + p1)*x + 1.0_dp) /  &
                   ((((q4*x + q3)*x + q2)*x + q1)*x + 1.0_dp))
               RETURN
           END IF

           IF (x >= 0.0_dp) THEN
               e = EXP(x)
               fn_val = e * (half + (half - 1.0_dp/e))
               RETURN
           END IF
           IF (x >= -37.0_dp) THEN
               fn_val = (EXP(x) - half) - half
               RETURN
           END IF
           fn_val = -1.0_dp
           RETURN
       END FUNCTION rexp
END SUBROUTINE cgamma


 

end module UtilitySubroutines
