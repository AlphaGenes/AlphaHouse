
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine IntitialiseIntelRNG

  ! Start an RNG stream

  ! https://software.intel.com/en-us/node/470610 (2014-11-25)
  ! https://software.intel.com/en-us/node/470612 (2014-11-25)

  implicit none

  integer(kind=4) :: BRNG
  integer(kind=4),dimension(1) :: LIdumI

  logical :: LFileExistsL

  character(len=80) :: name

  ! Get seed from a file (if it exists) or else from the current date and time
  call get_command_argument(0, name)
  inquire(file=trim(name)//"Seed.txt", exist=LFileExistsL)
  if (LFileExistsL) then
    open(unit=1,file=trim(name)//"Seed.txt",status="old")
    read(unit=1,fmt=*) LIdumI(1)
    close(unit=1)
    print*,"Seed value taken from the file ",trim(name)//"Seed.txt: ",LIdumI
  else
    call RANDOM_SEED
    call RANDOM_SEED(get=LIdumI(1:1)) ! need to pass a vector
    print*,"Seed value generated from the current date and time: ",LIdumI
  endif

  ! Save the used seed for reproducibility
  open(unit=1,file=trim(name)//"SeedUsed.txt",status="unknown")
  write(unit=1,fmt=*) LIdumI
  close(unit=1)

  ! Start a stream
  BRNG=VSL_BRNG_MT19937
  RNGErrCode=vslnewstream(RNGStream,brng,LIdumI(1))
  if (RNGErrCode /= vsl_status_ok) then
    print*,"IntitialiseIntelRNG failed"
    stop
  endif

end subroutine IntitialiseIntelRNG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
