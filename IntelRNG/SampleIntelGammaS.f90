
!###############################################################################

function SampleIntelGammaS(n,alpha,a,beta) result(Res)

  ! Sample n values from a Gamma(alpha,a,beta) distribution; single precision
  ! n input (integer), number of samples to generate (default 1)
  ! alpha input (real), shape
  ! a input (real), displacement (default 0.0)
  ! beta input (real), scale

  ! https://software.intel.com/en-us/node/470672 (2016-03-07)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real32),intent(in)            :: alpha
  real(real32),intent(in),optional   :: a
  real(real32),intent(in)            :: beta
  real(real32),allocatable           :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real32) :: aOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (.not. (alpha > 0.0)) then
    write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires alpha (shape) parameter to be greater than zero"
    write(STDERR,"(a)") " "
    stop 1
  end if

  if (present(a)) then
    aOpt=a
  else
    aOpt=0.0
  end if

  if (.not. (beta > 0.0)) then
    write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires beta (scale) parameter to be greater than zero"
    write(STDERR,"(a)") " "
    stop 1
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_GAMMA_GNORM
  RNGErrCode=vsrnggamma(RNGMethod,RNGStream,nOpt,Res,alpha,aOpt,beta)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelGammaS failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
