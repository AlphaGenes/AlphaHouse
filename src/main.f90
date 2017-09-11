program test

	use iso_fortran_env
	use individualModule
	use pedigreeModule
	use compatibilityModule
	implicit none

	type(PedigreeHolder) :: ped
	integer :: nsnp
  character(len=128), dimension(:), allocatable :: chromPaths

	call readPlink("merge_final", ped,chromPaths,nsnp)

end program test

