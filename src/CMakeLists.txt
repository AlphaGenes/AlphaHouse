
SET(ALPHAHOUSE_src 
				AlphaHouseMod.f90
				AlphaStatMod.f90
				HashModule.f90
				LinkedListModule.f90
				
				IndividualModule.f90
				IndividualLinkedListModule.f90
				PedigreeModule.f90				
				IntelRNGMod.f90
				HelperSubroutines.f90
				# test.f90
	)


if (ALPHAHOUSEEXE)
	ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})
else()
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
