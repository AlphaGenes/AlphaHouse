
SET(ALPHAHOUSE_src ./AlphaHouseMod.f90
				./AlphaStatMod.f90
	)
ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})