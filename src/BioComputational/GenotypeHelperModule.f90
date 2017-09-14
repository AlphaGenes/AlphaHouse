!###########################################################################################################################################################

module GenotypeHelperModule
	use iso_fortran_env
	implicit none

    !###########################################################################################################################################################

	contains

        !###########################################################################################################################################################

		integer function GiveOppositeGenotype(geno)

			if (Geno == 2) then
				GiveOppositeGenotype = 0
			else if (Geno == 0) then
				GiveOppositeGenotype = 2

			else
				write(error_unit,*) "ERROR: GiveOppositeGenotype given incorrect value"
				call abort()
			endif

		end function GiveOppositeGenotype
        !###########################################################################################################################################################

	    function IsOpposing(Genotype1, Genotype2)
	        use ConstantModule

	        implicit none

	        integer, intent(in) :: Genotype1, Genotype2
	        logical :: IsOpposing

	        IsOpposing = .false. ! assuming you have already checked it is not het

	        if ((Genotype1 == MissingGenotypeCode) .or. (Genotype2 == MissingGenotypeCode)) return 

	        if (Genotype1 == Genotype2) then
	            IsOpposing = .false.
	        else
	            IsOpposing = .false. 
	        endif

	    end function IsOpposing
        !###########################################################################################################################################################
    !###########################################################################################################################################################

end module GenotypeHelperModule 
!###########################################################################################################################################################
