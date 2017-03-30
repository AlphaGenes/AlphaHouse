
SET(ALPHAHOUSE_src
    Utilities/ConstantModule.f90
    Utilities/AlphaHouseMod.f90
    Utilities/AlphaSystemMod.f90
    Utilities/AlphaStatMod.f90
    Utilities/OrderPackModule.f90
    Utilities/HashModule.f90
    Utilities/LinkedListModule.f90
    BioComputational/IndividualModule.f90
    BioComputational/IndividualLinkedListModule.f90
    BioComputational/IndividualHelperModule.f90
    BioComputational/PedigreeModule.f90
    BioComputational/IntelRNGMod.f90
    NotUnitTested/Utilities/HelperSubroutines.f90
    Utilities/StringModule.f90
    Utilities/LineModule.f90
    Utilities/PageModule.f90
    Utilities/Input.f90
    NotUnitTested/Utilities/CompatibilityModule.f90
    NotUnitTested/Utilities/AlphaSortMod.f90
    Utilities/IntegerLinkedListModule.f90
    NotUnitTested/BioComputational/GenotypeModule.f90
    NotUnitTested/BioComputational/HaplotypeModule.f90
    NotUnitTested/BioComputational/BitUtilities.f90
    NotUnitTested/BioComputational/FilterModule.f90
    # NotUnitTested/BioComputational/CoreIndexModule.f90
    NotUnitTested/Utilities/SortingModule.f90
    NotUnitTested/Utilities/ClusteringModule.f90
    NotUnitTested/Utilities/par_zig_mod.f90
    Utilities/LegacyRandom.f90
    NotUnitTested/BioComputational/AlphaEvolveModule.f90
    main.f90
)


if (ALPHAHOUSEEXE)
	ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})
else()
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
