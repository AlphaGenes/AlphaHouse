module MPIUtilities
	use mpi
	use ISO_Fortran_Env
	implicit none
	integer(int32), save, protected:: mpiSize, mpiRank
	integer(int32), save, protected:: mpiCommunicator
	INTEGER, PARAMETER:: DEFAULTMPIRANK=0

	interface startMPI
		module procedure initialiseMPI
	end interface
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

		subroutine initialiseMPI
			integer(int32)::error
			call MPI_Init(error)
			call checkMPI(error)

			call setMPICommunicator()

			call setMPISizes()

		end subroutine initialiseMPI

		subroutine endMPI
			integer(int32)::error
			integer(int32), PARAMETER:: Success=0
			call MPI_FINALIZE(error)

			if (error .ne. Success) then
				write(*,*) "Problem with finishing MPI"
			end if
		end subroutine endMPI

		subroutine broadcastAllocatableCharacter(characterIn, mpiRankUsed)
			character(len=:), allocatable, intent(inout):: characterIn
			integer, intent(in), optional:: mpiRankUsed
			integer:: messageLength


			integer:: mpiErr

			messageLength=len(characterIn)
			call MPI_BCAST(messageLength, 1, MPI_INTEGER, mpiRankUsed, mpiCommunicator, mpiErr)
			call checkMPI(mpiErr)

			if (mpiRank .ne. mpiRankUsed) then
				if (allocated(characterIn)) then
					deallocate(characterIn)
				end if
				allocate(character(len=messageLength)::characterIn)
			end if

			call MPI_BCAST(characterIn, messageLength, MPI_CHARACTER, mpiRankUsed, mpiCommunicator, mpiErr)
			call checkMPI(mpiErr)
		end subroutine broadcastAllocatableCharacter


		subroutine getMPIChunckSize(totalSize, chunkSize, startObs, endObs)
			integer(int32), intent(in):: totalSize
			integer, intent(out):: chunkSize, startObs, endObs
			integer:: extras

			chunkSize = totalSize/mpiSize

			if (mpiRank==mpiSize-1) then
				extras = totalSize-chunkSize*mpiSize
			else
				extras = 0
			end if

			startObs = mpiRank*chunkSize+1
			endObs = (mpiRank+1)*chunkSize+extras

		end subroutine getMPIChunckSize


		!> @brief Set which communicator to use
		!> @details Mainly here for the tests.   Allows you to set the communicator from
		!> outside of the subroutine
		subroutine setMPICommunicator(commIn)
			integer, intent(in), optional:: commIn

			if (present(commIn)) then
				mpiCommunicator = commIn
			else
				mpiCommunicator = MPI_COMM_WORLD
			end if
		end subroutine setMPICommunicator

		subroutine setMPISizes()
			integer(int32):: error

			call MPI_COMM_SIZE(mpiCommunicator, mpiSize,error)
			call checkMPI(error)

			call MPI_COMM_RANK(mpiCommunicator, mpiRank, error)
			call checkMPI(error)

		end subroutine setMPISizes

		!> @brief writes out on the first MPI rank
		subroutine writeOutToScreen(characterIn)
			use, intrinsic :: iso_fortran_env, only: stdout => output_unit
			character(len=*), intent(in):: characterIn

#ifdef MPIACTIVE
			if (mpiRank ==0) then
				write(stdout,*) characterIn
			end if
#else
			write(stdout,*) characterIn
#endif
		end subroutine writeOutToScreen

end module MPIUtilities



