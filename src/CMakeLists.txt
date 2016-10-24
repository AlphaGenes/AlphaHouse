
SET(ALPHAHOUSE_src ./AlphaHouseMod.f90
				./AlphaStatMod.f90
				PedigreeTable.f90
				IntelRNGMod.f90
				HelperSubroutines.f90
	)
ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})