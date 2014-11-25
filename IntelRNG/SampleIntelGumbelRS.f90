
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelGumbelRS(n,a,b)

  ! Sample n values from a Gumbel(a,b) (=type-I extreme value) distribution; kind=4
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), location (default 0.0)
  ! b input (real), scale (default 1.0)

  ! https://software.intel.com/en-us/node/470670 (2014-11-25)

  implicit none

  integer(kind=4),optional :: n
  integer(kind=4) :: nOpt

  real(kind=4),optional :: a,b
  real(kind=4) :: aOpt,bOpt
  real(kind=4),dimension(:),allocatable :: SampleIntelGumbelRS

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

  allocate(SampleIntelGumbelRS(nOpt))

  RNGMethod=VSL_RNG_METHOD_GUMBEL_ICDF
  RNGErrCode=vsrnggumbel(RNGMethod,RNGStream,nOpt,SampleIntelGumbelRS,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelGumbelRS failed"
    stop
  endif

  return

end function SampleIntelGumbelRS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
