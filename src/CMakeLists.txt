
SET(ALPHAHOUSE_src 
				AlphaHouseMod.f90
				AlphaStatMod.f90
				HashModule.f90
				LinkedListModule.f90
				IndividualLinkedListModule.f90
				PedigreeModule.f90
				IndividualModule.f90
				IntelRNGMod.f90
				HelperSubroutines.f90
	)
ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})