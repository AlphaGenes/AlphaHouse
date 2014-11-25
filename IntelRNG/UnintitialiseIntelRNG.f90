
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine UnintitialiseIntelRNG

  ! Delete an RNG stream

  ! https://software.intel.com/en-us/node/470610 (2014-11-25)
  ! https://software.intel.com/en-us/node/470612 (2014-11-25)

  implicit none

  RNGErrCode=vsldeletestream(RNGStream)
  if (RNGErrCode /= vsl_status_ok) then
    print*,"UnintitialiseIntelRNG failed"
    stop
  endif

end subroutine UnintitialiseIntelRNG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
