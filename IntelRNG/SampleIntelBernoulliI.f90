
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelBernoulliI(n,p)

  ! Sample n values from a Bernoulli(p) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! p input (real), probability of of a success (default 0.5)

  ! https://software.intel.com/en-us/node/470686 (2014-11-25)

  implicit none

  integer(kind=4),optional,intent(in) :: n
  integer(kind=4) :: nOpt
  integer(kind=4),dimension(:),allocatable :: SampleIntelBernoulliI

  real(kind=8),optional,intent(in) :: p
  real(kind=8) :: pOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(p)) then
    pOpt=p
  else
    pOpt=0.5
  end if

  allocate(SampleIntelBernoulliI(nOpt))

  RNGMethod=VSL_RNG_METHOD_BERNOULLI_ICDF
  RNGErrCode=virngbernoulli(RNGMethod,RNGStream,nOpt,SampleIntelBernoulliI,pOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelBernoulliI failed"
    stop
  end if

  return

end function SampleIntelBernoulliI

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
