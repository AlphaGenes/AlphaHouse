module MPIUtilities
use mpi
use ISO_Fortran_Env
implicit none
integer(int32), save:: mpiSize, mpiRank

contains

subroutine checkMPI(errorIn)
  integer(int32):: errorIn
  integer(int32):: mpierr, length
  character(len=1000)::errorString

  if (errorIn .ne. 0) then
  call MPI_Error_string(errorIn, errorString, length, mpierr)
  write(*,*) "MPI Error!"
  write(*,*) trim(errorString)
  stop
  end if

end subroutine checkMPI


function getMPIChunckSize(totalSize) result (chunkSize)
  integer(int32):: totalSize, chunkSize

  chunkSize =totalSize/mpiSize
end function


subroutine initialiseMPI
  integer(int32)::error
  call MPI_Init(error)
  call checkMPI(error)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize,error)
  call checkMPI(error)

  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, error)
  call checkMPI(error)

end subroutine initialiseMPI

subroutine endMPI
  integer(int32)::error
  integer(int32), PARAMETER:: Success=0
  call MPI_FINALIZE(error)

  if (error .ne. Success) then
    write(*,*) "Problem with finishing MPI"
  end if
end subroutine endMPI


end module MPIUtilities

