
!###############################################################################

function SampleIntelUniformS(n,a,b)

  ! Sample n values values from a Uniform(a,b) distribution; single precision
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), minimal value (default 0.0)
  ! b input (real), maximal value (default 1.0)

  ! https://software.intel.com/en-us/node/470652 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real32),intent(in),optional   :: a
  real(real32),intent(in),optional   :: b

  ! Other
  integer(int32) :: nOpt

  real(real32) :: aOpt,bOpt
  real(real32),allocatable :: SampleIntelUniformS(:)

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

  allocate(SampleIntelUniformS(nOpt))

  RNGMethod=VSL_RNG_METHOD_UNIFORM_STD ! should we use here VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
  RNGErrCode=vsrnguniform(RNGMethod,RNGStream,nOpt,SampleIntelUniformS,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelUniformS failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
