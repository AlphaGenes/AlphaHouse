
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelUniformRS(n,a,b)

  ! Sample n values values from a Uniform(a,b) distribution; kind=4
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), minimal value (default 0.0)
  ! b input (real), maximal value (default 1.0)

  ! https://software.intel.com/en-us/node/470652 (2014-11-25)

  implicit none

  integer(kind=4),optional :: n
  integer(kind=4) :: nOpt

  real(kind=4),optional :: a,b
  real(kind=4) :: aOpt,bOpt
  real(kind=4),dimension(:),allocatable :: SampleIntelUniformRS

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  endif

  if (present(a)) then
    aOpt=a
  else
    aOpt=0.0
  endif

  if (present(b)) then
    bOpt=b
  else
    bOpt=1.0
  endif

  allocate(SampleIntelUniformRS(nOpt))

  RNGMethod=VSL_RNG_METHOD_UNIFORM_STD ! should we use here VSL_RNG_METHOD_UNIFORM_STD_ACCURATE?
  RNGErrCode=vsrnguniform(RNGMethod,RNGStream,nOpt,SampleIntelUniformRS,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelUniformRS failed"
    stop
  endif

  return

end function SampleIntelUniformRS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
