
SET(ALPHAHOUSE_src
		Utilities/ConstantModule.f90
        Utilities/AlphaHouseMod.f90
        Utilities/AlphaSystemMod.f90
		BioComputational/AlphaStatMod.f90
    Utilities/OrderPackModule.f90
		Utilities/HashModule.f90
		Utilities/LinkedListModule.f90
		BioComputational/IndividualModule.f90
		BioComputational/IndividualLinkedListModule.f90
		BioComputational/IndividualHelperModule.f90
		BioComputational/PedigreeModule.f90
		BioComputational/IntelRNGMod.f90
		Utilities/HelperSubroutines.f90
        Utilities/StringModule.f90
        Utilities/LineModule.f90
        Utilities/PageModule.f90
        Utilities/Input.f90
        Utilities/PlinkHandler.f90
	Utilities/AlphaSortMod.f90
	Utilities/IntegerLinkedListModule.f90
	BioComputational/GenotypeModule.f90
	BioComputational/HaplotypeModule.f90
	BioComputational/BitUtilities.f90
	BioComputational/FilterModule.f90
	# BioComputational/CoreIndexModule.f90
        Utilities/SortingModule.f90
	Utilities/ClusteringModule.f90
	Utilities/par_zig_mod.f90
	Utilities/LegacyRandom.f90
        main.f90
	)


if (ALPHAHOUSEEXE)
	ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})
else()
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
