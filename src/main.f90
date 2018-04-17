program test

	use iso_fortran_env
	use individualModule
	use pedigreeModule
	use alphahousemod
	use compatibilityModule
	implicit none

	character(len=:), allocatable :: path

	call getExecutablePath(path)


end program test

