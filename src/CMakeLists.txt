# Sets source files that are going to be compiled
SET(ALPHAHOUSE_src
    Utilities/ConstantModule.f90
    Utilities/AlphaHouseMod.f90
    Utilities/AlphaSystemMod.f90
    Utilities/OrderPackModule.f90
    Utilities/HashModule.f90
    Utilities/LinkedListModule.f90
    BioComputational/AlphaStatMod.f90
    BioComputational/IndividualModule.f90
    BioComputational/IndividualLinkedListModule.f90
    BioComputational/IndividualHelperModule.f90
    BioComputational/IntelRNGMod.f90
    Utilities/HelperSubroutines.f90
    Utilities/StringModule.f90
    Utilities/LineModule.f90
    Utilities/PageModule.f90
    Utilities/Input.f90
    Utilities/CompatibilityModule.f90
    Utilities/AlphaSortMod.f90
    Utilities/IntegerLinkedListModule.f90
    Utilities/CharacterLinkedListModule.f90
    Utilities/SortedIntegerLinkedListModule.f90
    BioComputational/PedigreeModule.f90
    BioComputational/GenotypeModule.f90
    BioComputational/HaplotypeModule.f90
    BioComputational/HaplotypeLibraryModule.f90
    BioComputational/HaplotypeLinkedListModule.f90
    BioComputational/BitUtilities.f90
    BioComputational/FilterModule.f90
    BioComputational/informationModule.f90
    BioComputational/GenotypeHelperModule.f90
    # BioComputational/CoreIndexModule.f90
    BioComputational/CoreUtils.f90
    Utilities/SortingModule.f90
    Utilities/ClusteringModule.f90
    Utilities/alphafullchrom.f90
    Utilities/par_zig_mod.f90
    Utilities/LegacyRandom.f90
    BioComputational/AlphaEvolveModule.f90
    Utilities/baseSpecFileModule.f90
    main.f90
)

if (USE_MPI)
    set(ALPHAHOUSE_src ${ALPHAHOUSE_src} Utilities/MPIUtilities.f90)
endif()

if (ALPHAHOUSEEXE)
    ADD_EXECUTABLE(${ALPHAHOUSEEXE} ${ALPHAHOUSE_src})

    # link to the mkl libs
    
    TARGET_LINK_LIBRARIES(${ALPHAHOUSEEXE} ${MKL_MIX})

    # link to netcdf if we need
    if (USE_NETCDF)

    target_link_libraries (${ALPHAHOUSEEXE} ${NETCDF_LIBRARIES})

    endif()
else()
# if we don't actually want to build an executable, just add library statically
	ADD_LIBRARY(${AHLIB} STATIC ${ALPHAHOUSE_src})
endif()
