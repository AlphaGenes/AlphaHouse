
SET(ALPHAHOUSE_src
        AlphaHouseMod.f90
        AlphaSystemMod.f90
				AlphaStatMod.f90
				HashModule.f90
				LinkedListModule.f90
        DiarmaidLinkedList.f90
				IndividualModule.f90
				IndividualLinkedListModule.f90
				IndividualHelperModule.f90
				PedigreeModule.f90				
				IntelRNGMod.f90
				HelperSubroutines.f90
        StringModule.f90
        LineModule.f90
        PageModule.f90
        Input.f90
        main.f90
	)


if (ALPHAHOUSEEXE)
	ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})
else()
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
