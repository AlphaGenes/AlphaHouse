
!###############################################################################

function SampleIntelUniformI(n,a,b)

  ! Sample n values from a discrete Uniform(a,b) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! a input (integer), minimal value (default 0)
  ! b input (integer), maximal value (default 1)

  ! https://software.intel.com/en-us/node/470678 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  integer(int32),intent(in),optional :: a
  integer(int32),intent(in),optional :: b

  ! Other
  integer(int32) :: nOpt,aOpt,bOpt
  integer(int32),allocatable :: SampleIntelUniformI(:)

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

  allocate(SampleIntelUniformI(nOpt))

  RNGMethod=VSL_RNG_METHOD_UNIFORM_STD
  RNGErrCode=virnguniform(RNGMethod,RNGStream,nOpt,SampleIntelUniformI,aOpt,bOpt+1)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelUniformI failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
