
!###############################################################################

include "mkl_vsl.f90"

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     IntelRNGMod.f90
!
! DESCRIPTION:
!> @brief    Conventient interfaces to Intel MKL Vector Statistical Library (VSL) RNG capabilities
!
!> @details  See https://software.intel.com/en-us/node/470592 (2014-11-25)
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
module IntelRNGMod

  use mkl_vsl_type
  use mkl_vsl
  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use AlphaHouseMod, only : GetSeed

  implicit none

  integer(int32) RNGErrCode,RNGMethod

  type(vsl_stream_state) :: RNGStream

  private

  ! RNG stream management
  public :: IntitialiseIntelRNG,IntitialiseIntelRNGStream,UnintitialiseIntelRNG,UnintitialiseIntelRNGStream

  ! Discrete
  public :: getIntelUniformI
  public :: SampleIntelUniformI
  public :: SampleIntelBernoulliI
  public :: SampleIntelMultinomialI
  public :: SampleIntelPoissonI

  public :: RandomSample

  ! Continuous
  !>@todo: should we make an interface and have generic for either single or double precision
  !!       and determine single or double based on inputs or???
  public :: getIntelUniformD
  public :: SampleIntelUniformS,SampleIntelUniformD
  public :: SampleIntelGaussS,SampleIntelGaussD
  public :: SampleIntelGammaS,SampleIntelGammaD
  public :: SampleIntelGumbelS,SampleIntelGumbelD

  ! Return a single sample ratger than an array
  public :: SingleSampleIntelUniformI

  contains

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Subroutine that returns an  increasing set of random numbers up until a maximum
    !> @author  David Wilson, david.wilson@roslin.ed.ac.uk
    !> @date    October 21, 2016
    !> @return  array of increasing random numbers
    !---------------------------------------------------------------------------
    subroutine generateIncreasingRandom(amount,maxIndex,randomNumbers)
      integer, intent(in) :: amount,maxIndex
      integer,dimension(:),allocatable, intent(out) :: randomNumbers
      integer :: i
      real(kind=real64) :: delta, tmp
      real(kind=int64) :: res(1)
      allocate(randomNumbers(amount))
      delta = maxIndex / real(amount);
      call IntitialiseIntelRNG()
      tmp = real(maxIndex)
    do i=1, amount
      res = SampleIntelUniformD(b=tmp)
        randomNumbers(i) = int(i*delta + res(1) * delta);
    enddo

    call UnintitialiseIntelRNG()
    end subroutine generateIncreasingRandom
    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Start an Intel RNG stream (including the seed)
    !> @details See https://software.intel.com/en-us/node/470610 (2014-11-25)
    !!          and https://software.intel.com/en-us/node/470612 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !> @return  Started RNG stream in object RNGStream, potentially created file
    !!          (SeedFile), and potentially returned seed value (Out)
    !---------------------------------------------------------------------------
    subroutine IntitialiseIntelRNG(Seed,SeedFile,Out)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional  :: Seed     !< A number to initialize RNG with
      character(len=*),optional           :: SeedFile !< File to save the seed in
      integer(int32),intent(out),optional :: Out      !< Make the seed value available outside

      ! Other
      integer(int32) :: SeedInt,Unit,BRNG

      if (present(Seed)) then
        SeedInt=Seed
      else
        call GetSeed(Out=SeedInt)
      end if

      ! Save to a file
      if (present(SeedFile)) then
        open(newunit=Unit,file=trim(SeedFile),status="unknown")
        write(Unit,*) SeedInt
        close(Unit)
      end if

      ! Output
      if (present(Out)) then
        Out=SeedInt
      end if

      ! Start a RNG stream
      BRNG=VSL_BRNG_MT19937
      RNGErrCode=vslnewstream(RNGStream,BRNG,SeedInt)
      if (RNGErrCode /= VSL_STATUS_OK) then
        write(STDERR,"(a)") "ERROR: IntitialiseIntelRNG failed"
        write(STDERR,"(a)") " "
        stop 1
      end if
    end subroutine


    subroutine IntitialiseIntelRNGStream(Stream, Seed,SeedFile,Out)
      implicit none

      ! Arguments
      type(vsl_stream_state), intent(out) :: Stream
      integer(int32),intent(in),optional  :: Seed     !< A number to initialize RNG with
      character(len=*),optional           :: SeedFile !< File to save the seed in
      integer(int32),intent(out),optional :: Out      !< Make the seed value available outside

      ! Other
      integer(int32) :: SeedInt,Unit,BRNG

      if (present(Seed)) then
        SeedInt=Seed
      else
        call GetSeed(Out=SeedInt)
      end if

      ! Save to a file
      if (present(SeedFile)) then
        open(newunit=Unit,file=trim(SeedFile),status="unknown")
        write(Unit,*) SeedInt
        close(Unit)
      end if

      ! Output
      if (present(Out)) then
        Out=SeedInt
      end if

      ! Start a RNG stream
      BRNG=VSL_BRNG_MT19937
      RNGErrCode=vslnewstream(Stream,BRNG,SeedInt)
      if (RNGErrCode /= VSL_STATUS_OK) then
        write(STDERR,"(a)") "ERROR: IntitialiseIntelRNG failed"
        write(STDERR,"(a)") " "
        stop 1
      end if
    end subroutine
    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   End an Intel RNG stream
    !> @details See https://software.intel.com/en-us/node/470610 (2014-11-25)
    !!          and https://software.intel.com/en-us/node/470612 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    subroutine UnintitialiseIntelRNG
      implicit none

      RNGErrCode=vsldeletestream(RNGStream)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: UnintitialiseIntelRNG failed"
        write(STDERR,"(a)") " "
        stop 1
      end if
    end subroutine


  
  subroutine UnintitialiseIntelRNGStream(stream)
    type(vsl_stream_state), intent(out) :: Stream
      

      RNGErrCode=vsldeletestream(Stream)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: UnintitialiseIntelRNG failed"
        write(STDERR,"(a)") " "
        stop 1
      end if
    end subroutine


    integer function getIntelUniformI()

    integer,allocatable :: sample(:)


    allocate(sample(1))
    sample = SampleIntelUniformI()

    getIntelUniformI = sample(1)

    end function


    function getIntelUniformD() result(res)

    real(real64),allocatable :: sample(:)
    real(real64) :: res

    allocate(sample(1))
    sample = SampleIntelUniformD()

    res = sample(1)

    end function
    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a discrete Uniform(a,b) distribution
    !> @details See https://software.intel.com/en-us/node/470678 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelUniformI(n,a,b,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      integer(int32),intent(in),optional :: a      !< minimal value (inclusive) (default 0)
      integer(int32),intent(in),optional :: b      !< maximal value (inclusive) (default 1)
      type(vsl_stream_state), intent(in),optional :: Stream
      integer(int32),allocatable         :: Res(:) !< @return samples

      ! Other
      integer(int32) :: nOpt,aOpt,bOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(a)) then
        aOpt=a
      else
        aOpt=0
      end if

      if (present(b)) then
        bOpt=b
      else
        bOpt=1
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_UNIFORM_STD
      if (present(stream)) then
      RNGErrCode=virnguniform(RNGMethod,Stream,nOpt,Res,aOpt,bOpt+1)
      else 
        RNGErrCode=virnguniform(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt+1)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelUniformI failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Bernoulli(p) distribution
    !> @details See https://software.intel.com/en-us/node/470686 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelBernoulliI(n,p,Stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: p      !< probability of of a success (default 0.5)
      integer(int32),allocatable         :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt

      real(real64) :: pOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(p)) then
        pOpt=p
      else
        pOpt=0.5d0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_BERNOULLI_ICDF
      if (present(stream)) then
        RNGErrCode=virngbernoulli(RNGMethod,stream,nOpt,Res,pOpt)
      else 
        RNGErrCode=virngbernoulli(RNGMethod,RNGStream,nOpt,Res,pOpt)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelBernoulliI failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Multinomial(p) distribution
    !> @details Not from Intel
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelMultinomialI(n,p,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in)            :: p(:)   !< probabilities for the different categories
      integer(int32),allocatable         :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt,i,j,k
      integer(int32) :: b(1)

      real(real64) :: pi,psum,psumtmp
      real(real64),allocatable :: pInternal(:)

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      k=size(p)
      allocate(pInternal(k))
      pInternal(:)=p(:)

      psum=sum(pInternal)
      if (abs(psum - 1.0d0) > 1e-7) then
        pInternal=pInternal/psum ! rescale
      end if
      psum=1.0d0

      allocate(Res(nOpt))

      ! Over samples
      do j=1,nOpt

        ! Over categories of a sample
        psumtmp=psum
        do i=1,k

          if (pInternal(i) > 0.0d0) then ! sample only if needed (border case)

            if (pInternal(i) < 1.0d0) then ! likewise
              pi=pInternal(i)/psumtmp

              if (pi < 1.0d0) then ! likewise
                if (present(stream)) then
                  b=SampleIntelBernoulliI(p=pi,stream = stream)
                else 
                b=SampleIntelBernoulliI(p=pi)
                endif

                if (b(1) > 0) then
                  Res(j)=i
                  exit
                end if

                psumtmp=psumtmp-pInternal(i)

              else
                Res(j)=i
                exit
              end if

            else
              Res(j)=i
              exit
            end if

          end if

        end do

      end do

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Poisson(lambda) distribution
    !> @details See https://software.intel.com/en-us/node/470694 (2014-12-01)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelPoissonI(n,lambda) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: lambda !< mean and variance of distribution (default 1.0)
      integer(int32),allocatable         :: Res(:) !< @return samples

      ! Other
      integer(int32) :: nOpt

      real(real64) :: lambdaOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(lambda)) then
        lambdaOpt=lambda
      else
        lambdaOpt=1.0d0
      end if

      allocate(Res(nOpt))

      if (lambdaOpt < 27.0d0) then ! based on Intel documentation https://software.intel.com/en-us/node/470598 (2014-12-01)
        RNGMethod=VSL_RNG_METHOD_POISSON_POISNORM
      else
        RNGMethod=VSL_RNG_METHOD_POISSON_PTPE
      end if
      RNGErrCode=virngpoisson(RNGMethod,RNGStream,nOpt,Res,lambdaOpt)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelPoissonI failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Uniform(a,b) distribution (single precision)
    !> @details See https://software.intel.com/en-us/node/470652 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelUniformS(n,a,b,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real32),intent(in),optional   :: a      !< minimal value (default 0.0)
      real(real32),intent(in),optional   :: b      !< maximal value (default 1.0)
      real(real32),allocatable           :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt

      real(real32) :: aOpt,bOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(a)) then
        aOpt=a
      else
        aOpt=0.0
      end if

      if (present(b)) then
        bOpt=b
      else
        bOpt=1.0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_UNIFORM_STD ! should we use here VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
      if (present(stream)) then
        RNGErrCode=vsrnguniform(RNGMethod,stream,nOpt,Res,aOpt,bOpt)
      else
        RNGErrCode=vsrnguniform(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelUniformS failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Uniform(a,b) distribution (double precision)
    !> @details See https://software.intel.com/en-us/node/470652 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelUniformD(n,a,b,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: a      !< minimal value (default 0.0)
      real(real64),intent(in),optional   :: b      !< maximal value (default 1.0)
      real(real64),allocatable           :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt

      real(real64) :: aOpt,bOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(a)) then
        aOpt=a
      else
        aOpt=0.d0
      end if

      if (present(b)) then
        bOpt=b
      else
        bOpt=1.d0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
      if (present(stream)) then
        RNGErrCode=vdrnguniform(RNGMethod,stream,nOpt,Res,aOpt,bOpt)
      else 
        RNGErrCode=vdrnguniform(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelUniformD failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gauss(mu,sigma2) distribution (single precision)
    !> @details See https://software.intel.com/en-us/node/470654 (2016-03-07)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGaussS(n,mu,sigma2,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real32),intent(in),optional   :: mu     !< mean (default 0.0)
      real(real32),intent(in),optional   :: sigma2 !< variance (default 1.0)
      real(real32),allocatable           :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt

      real(real32) :: muOpt,sigma

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(mu)) then
        muOpt=mu
      else
        muOpt=0.0
      end if

      if (present(sigma2)) then
        sigma=sqrt(sigma2)
      else
        sigma=1.0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER

      if (present(stream)) then
        RNGErrCode=vsrnggaussian(RNGMethod,stream,nOpt,Res,muOpt,sigma)
      else 
        RNGErrCode=vsrnggaussian(RNGMethod,RNGStream,nOpt,Res,muOpt,sigma)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGaussS failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gauss(mu,sigma2) distribution (double precision)
    !> @details See https://software.intel.com/en-us/node/470654 (2016-03-07)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGaussD(n,mu,sigma2,stream) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: mu     !< mean (default 0.0)
      real(real64),intent(in),optional   :: sigma2 !< variance (default 1.0)
      real(real64),allocatable           :: Res(:) !< @return samples
      type(vsl_stream_state), intent(in),optional :: Stream
      ! Other
      integer(int32) :: nOpt

      real(real64) :: muOpt,sigma

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(mu)) then
        muOpt=mu
      else
        muOpt=0.0d0
      end if

      if (present(sigma2)) then
        sigma=sqrt(sigma2)
      else
        sigma=1.0d0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    
      if (present(stream)) then
        RNGErrCode=vdrnggaussian(RNGMethod,stream,nOpt,Res,muOpt,sigma)
      else
        RNGErrCode=vdrnggaussian(RNGMethod,RNGStream,nOpt,Res,muOpt,sigma)
      endif
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGaussD failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gamma(alpha,a,beta) distribution (single precision)
    !> @details See https://software.intel.com/en-us/node/470672 (2016-03-07)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGammaS(n,alpha,beta,shape,scale,rate,shift) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real32),intent(in),optional   :: alpha  !< shape
      real(real32),intent(in),optional   :: beta   !< scale
      real(real32),intent(in),optional   :: shape  !< the same as alpha, but to make life easier
      real(real32),intent(in),optional   :: scale  !< the same as beta, but to make life easier
      real(real32),intent(in),optional   :: rate   !< the same as 1/beta, but to make life easier
      real(real32),intent(in),optional   :: shift  !< shift or displacement (default 0.0)
      real(real32),allocatable           :: Res(:) !< @return samples

      ! Other
      integer(int32) :: nOpt

      real(real32) :: shapeOpt,scaleOpt,shiftOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (.not.present(alpha) .and. .not.present(shape)) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires either alpha (shape) or shape argument"
        write(STDERR,"(a)") " "
      end if

      if (.not.present(beta) .and. .not.present(scale) .and. .not.present(rate)) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires either beta (scale), scale, or rate argument"
        write(STDERR,"(a)") " "
      end if

      if (present(alpha)) then
        if (.not. (alpha > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires alpha (shape) parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        shapeOpt=alpha
      end if

      if (present(shape)) then
        if (.not. (shape > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires shape parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        shapeOpt=shape
      end if

      if (present(beta)) then
        if (.not. (beta > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires beta (scale) parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=beta
      end if

      if (present(scale)) then
        if (.not. (scale > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires scale parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=scale
      end if

      if (present(rate)) then
        if (.not. (rate > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires rate parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=1.0/rate
        if (.not. (scaleOpt > 0.0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires rate parameter that gives non-zero scale (=1/rate) parameter"
          write(STDERR,"(a,e)") "ERROR: rate: ",rate
          write(STDERR,"(a,e)") "ERROR: scale: ",scaleOpt
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      if (present(shift)) then
        shiftOpt=shift
      else
        shiftOpt=0.0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE
      RNGErrCode=vsrnggamma(RNGMethod,RNGStream,nOpt,Res,shapeOpt,shiftOpt,scaleOpt)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaS failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gamma(alpha,a,beta) distribution (double precision)
    !> @details See https://software.intel.com/en-us/node/470672 (2016-03-07)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGammaD(n,alpha,beta,shape,scale,rate,shift) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: alpha  !< shape
      real(real64),intent(in),optional   :: beta   !< scale
      real(real64),intent(in),optional   :: shape  !< the same as alpha, but to make life easier
      real(real64),intent(in),optional   :: scale  !< the same as beta, but to make life easier
      real(real64),intent(in),optional   :: rate   !< the same as 1/beta, but to make life easier
      real(real64),intent(in),optional   :: shift  !< shift or displacement (default 0.0)
      real(real64),allocatable           :: Res(:) !< @return samples

      ! Other
      integer(int32) :: nOpt

      real(real64) :: shapeOpt,scaleOpt,shiftOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (.not.present(alpha) .and. .not.present(shape)) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires either alpha (shape) or shape argument"
        write(STDERR,"(a)") " "
      end if

      if (.not.present(beta) .and. .not.present(scale) .and. .not.present(rate)) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires either beta (scale), scale, or rate argument"
        write(STDERR,"(a)") " "
      end if

      if (present(alpha)) then
        if (.not. (alpha > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires alpha (shape) parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        shapeOpt=alpha
      end if

      if (present(shape)) then
        if (.not. (shape > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires shape parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        shapeOpt=shape
      end if

      if (present(beta)) then
        if (.not. (beta > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires beta (scale) parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=beta
      end if

      if (present(scale)) then
        if (.not. (scale > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires scale parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=scale
      end if

      if (present(rate)) then
        if (.not. (rate > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires rate parameter to be greater than zero"
          write(STDERR,"(a)") " "
          stop 1
        end if
        scaleOpt=1.0d0/rate
        if (.not. (scaleOpt > 0.0d0)) then
          write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires rate parameter that gives non-zero scale (=1/rate) parameter"
          write(STDERR,"(a,e)") "ERROR: rate: ",rate
          write(STDERR,"(a,e)") "ERROR: scale: ",scaleOpt
          write(STDERR,"(a)") " "
          stop 1
        end if
      end if

      if (present(shift)) then
        shiftOpt=shift
      else
        shiftOpt=0.0d0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE
      RNGErrCode=vdrnggamma(RNGMethod,RNGStream,nOpt,Res,shapeOpt,shiftOpt,scaleOpt)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGammaD failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gumbel(a,b) (=type-I extreme value) distribution (single precision)
    !> @details See https://software.intel.com/en-us/node/470670 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGumbelS(n,a,b) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real32),intent(in),optional   :: a      !< location (default 0.0)
      real(real32),intent(in),optional   :: b      !< scale (default 1.0)
      real(real32),allocatable           :: Res(:) !< @return samples

      ! Other
      integer(int32) :: nOpt

      real(real32) :: aOpt,bOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(a)) then
        aOpt=a
      else
        aOpt=0.0
      end if

      if (present(b)) then
        bOpt=b
      else
        bOpt=1.0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GUMBEL_ICDF
      RNGErrCode=vsrnggumbel(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGumbelS failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    !###########################################################################

    !---------------------------------------------------------------------------
    !> @brief   Sample from a Gumbel(a,b) (=type-I extreme value) distribution (double precision)
    !> @details See https://software.intel.com/en-us/node/470670 (2014-11-25)
    !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
    !> @date    September 26, 2016
    !---------------------------------------------------------------------------
    function SampleIntelGumbelD(n,a,b) result(Res)
      implicit none

      ! Arguments
      integer(int32),intent(in),optional :: n      !< number of samples to generate (default 1)
      real(real64),intent(in),optional   :: a      !< location (default 0.0)
      real(real64),intent(in),optional   :: b      !< scale (default 1.0)
      real(real64),allocatable           :: Res(:) !< @return samples
      ! Other
      integer(int32) :: nOpt

      real(real64) :: aOpt,bOpt

      if (present(n)) then
        nOpt=n
      else
        nOpt=1
      end if

      if (present(a)) then
        aOpt=a
      else
        aOpt=0.d0
      end if

      if (present(b)) then
        bOpt=b
      else
        bOpt=1.d0
      end if

      allocate(Res(nOpt))

      RNGMethod=VSL_RNG_METHOD_GUMBEL_ICDF
      RNGErrCode=vdrnggumbel(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
      if (RNGErrCode /= vsl_status_ok) then
        write(STDERR,"(a)") "ERROR: SampleIntelGumbelD failed"
        write(STDERR,"(a)") " "
        stop 1
      end if

      return
    end function

    function RandomSample(list,num) result(res)
        implicit none

        integer, dimension(:), intent(in) :: list
        integer, intent(in) :: num

        integer, dimension(num) :: res

        integer, dimension(size(list)) :: tmpA
        integer, dimension(:), allocatable :: r
        integer :: tmp, pos, i


        tmpA = list

        do i = 1, num
            r = SampleIntelUniformI(1,0,size(list) - i)
            pos = size(list) - r(1)
            tmp = tmpA(i)
            tmpA(i) = tmpA(pos)
            tmpA(pos) = tmp
        end do

        res = tmpA(1:num)
    end function RandomSample

    function SingleSampleIntelUniformI(a,b) result(res)
        integer, optional :: a, b
        integer(int32) :: res

        integer(int32), dimension(:), allocatable :: resA

        if (present(a)) then
            if (present(b)) then
                resA = SampleIntelUniformI(1,a,b)
            else
                resA = SampleIntelUniformI(1,a)
            end if
        else
            if (present(b)) then
                resA = SampleIntelUniformI(1,b=b)
            else
                resA = SampleIntelUniformI(1)
            end if
        end if

        res = resA(1)
    end function SingleSampleIntelUniformI

end module

!###############################################################################
