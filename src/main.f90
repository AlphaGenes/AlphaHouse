program test

	use iso_fortran_env
	use individualModule
	use pedigreeModule
	use compatibilityModule
	implicit none

	type(PedigreeHolder) :: ped
	integer,allocatable,dimension(:) :: nsnp
	integer :: inconsistencies
    character(len=300) :: inputFile,

    inputFile = "pedigree.txt"
    genotypeFile = "genotypes.txt"
    

	call initPedigree(ped,inputFile)
	
    call ped%addGenotypeInformationFromFile(genotypeFile,0)
    call ped%sortPedigreeAndOverwrite()

	inconsistencies = ped%findMendelianInconsistencies(file="mendInfo.txt", snpFilePath="snpinfo.txt")


    call ped%writeOutgenotypes("genotypesAfterMendellian.txt")


end program test

