program test

	use iso_fortran_env
	use individualModule
	use pedigreeModule
	use compatibilityModule
	implicit none

	type(PedigreeHolder) :: ped
	integer,allocatable,dimension(:) :: nsnp
	integer :: inconsistencies, nsnps
    character(len=300) :: inputFile,genotypeFile

    inputFile = "pedigree.txt"
    genotypeFile = "genotypes.txt"
    

	call initPedigree(ped,inputFile)
	nsnps = 0 
    call ped%addGenotypeInformationFromFile(genotypeFile,nsnps,initAll=1)
    call ped%sortPedigreeAndOverwrite()

	inconsistencies = ped%findMendelianInconsistencies(file="mendInfo.txt", snpFilePath="snpinfo.txt")


    call ped%writeOutgenotypes("genotypesAfterMendellian.txt")


end program test

