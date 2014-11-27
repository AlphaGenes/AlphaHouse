
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelGumbelRD(n,a,b)

  ! Sample n values from a Gumbel(a,b) (=type-I extreme value) distribution; kind=8
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), location (default 0.0)
  ! b input (real), scale (default 1.0)

  ! https://software.intel.com/en-us/node/470670 (2014-11-25)

  implicit none

  integer(kind=4),optional,intent(in) :: n
  integer(kind=4) :: nOpt

  real(kind=8),optional,intent(in) :: a,b
  real(kind=8) :: aOpt,bOpt
  real(kind=8),dimension(:),allocatable :: SampleIntelGumbelRD

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

  allocate(SampleIntelGumbelRD(nOpt))

  RNGMethod=VSL_RNG_METHOD_GUMBEL_ICDF
  RNGErrCode=vdrnggumbel(RNGMethod,RNGStream,nOpt,SampleIntelGumbelRD,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelGumbelRD failed"
    stop
  end if

  return

end function SampleIntelGumbelRD

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
