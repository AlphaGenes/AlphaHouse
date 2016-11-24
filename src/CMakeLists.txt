
SET(ALPHAHOUSE_src
        AlphaHouseMod.f90
				AlphaStatMod.f90
				IndividualModule.f90
				IntelRNGMod.f90
				HelperSubroutines.f90
        StringModule.f90
        LineModule.f90
        PageModule.f90
        Input.f90
	)
ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
