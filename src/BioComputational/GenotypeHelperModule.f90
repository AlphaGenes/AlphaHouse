module GenotypeHelperModule



contains


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



end module GenotypeHelperModule 