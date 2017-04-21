
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
    Untested/HelperSubroutines.f90
    Utilities/StringModule.f90
    Utilities/LineModule.f90
    Utilities/PageModule.f90
    Utilities/Input.f90
    Untested/CompatibilityModule.f90
    Untested/AlphaSortMod.f90
    Utilities/IntegerLinkedListModule.f90
    Untested/GenotypeModule.f90
    Untested/HaplotypeModule.f90
    Untested/BitUtilities.f90
    Untested/FilterModule.f90
    # Untested/CoreIndexModule.f90
    Untested/SortingModule.f90
    Untested/ClusteringModule.f90
    Untested/par_zig_mod.f90
    Utilities/LegacyRandom.f90
    Untested/AlphaEvolveModule.f90
    main.f90
)


if (ALPHAHOUSEEXE)
	ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})
else()
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
