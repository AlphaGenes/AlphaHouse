
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelUniformRD(n,a,b)

  ! Sample n values from a Uniform(a,b) distribution; kind=8
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), minimal value (default 0.0)
  ! b input (real), maximal value (default 1.0)

  ! https://software.intel.com/en-us/node/470652 (2014-11-25)

  implicit none

  integer(kind=4),optional :: n
  integer(kind=4) :: nOpt

  real(kind=8),optional :: a,b
  real(kind=8) :: aOpt,bOpt
  real(kind=8),dimension(:),allocatable :: SampleIntelUniformRD

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

  allocate(SampleIntelUniformRD(nOpt))

  RNGMethod=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
  RNGErrCode=vdrnguniform(RNGMethod,RNGStream,nOpt,SampleIntelUniformRD,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelUniformRD failed"
    stop
  end if

  return

end function SampleIntelUniformRD

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
