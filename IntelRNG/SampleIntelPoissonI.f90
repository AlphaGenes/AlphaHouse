
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelPoissonI(n,lambda)

  ! Sample n values from a Poisson(lambda) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! lambda input (real), mean and variance of distribution (default 1.0)

  ! https://software.intel.com/en-us/node/470694 (2014-12-01)

  implicit none

  integer(kind=4),optional,intent(in) :: n
  integer(kind=4) :: nOpt
  integer(kind=4),dimension(:),allocatable :: SampleIntelPoissonI

  real(kind=8),optional,intent(in) :: lambda
  real(kind=8) :: lambdaOpt

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

  allocate(SampleIntelPoissonI(nOpt))

  if (lambdaOpt < 27.0) then ! based on Intel documentation https://software.intel.com/en-us/node/470598 (2014-12-01)
    RNGMethod=VSL_RNG_METHOD_POISSON_POISNORM
  else
    RNGMethod=VSL_RNG_METHOD_POISSON_PTPE
  end if
  RNGErrCode=virngpoisson(RNGMethod,RNGStream,nOpt,SampleIntelPoissonI,lambdaOpt)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"SampleIntelPoissonI failed"
    stop
  end if

  return

end function SampleIntelPoissonI

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
