
SET(ALPHAHOUSE_src ./AlphaHouseMod.f90
				./AlphaStatMod.f90
	)
ADD_LIBRARY(${AHLIB} SHARED ${ALPHAHOUSE_src})