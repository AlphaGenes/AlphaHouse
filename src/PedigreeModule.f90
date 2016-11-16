module PedigreeModule

type PedigreeHolder

    type(Individual), allocatable, dimension(:) :: Pedigree
end type PedigreeHolder

contains
    

    function Pedigree(fileIn, numberInFile) result(pedStructure)
        use AlphaHouseMod, only : countLines
        type(PedigreeHolder) :: pedStructure
        character(len=*) :: fileIn
        integer(kind=int32),optional :: numberInFile
        integer(kind=int32) :: nIndividuals, fileUnit
        integer()
        if (present(numberInFile)) then
            nIndividuals = numberInFile
        else
            nIndividuals = countLines(fileIn)
        endif
        allocate(pedStructure%Pedigree(nIndividuals))

        open(newUnit=fileUnit, file=fileIn, status="old")
        do i=1,nIndividuals
            read(fileUnit,*) tmpId,tmpSire,tmpDam
        enddo




    end function Pedigree


end module PedigreeModule