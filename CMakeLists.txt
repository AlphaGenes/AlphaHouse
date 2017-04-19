CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(AlphaHouse)
enable_language (Fortran)
enable_testing()
find_package(Doxygen)
# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../alphacmakefiles/cmake/modules")
set(CMAKE_MACOSX_RPATH 0)

INCLUDE(${CMAKE_MODULE_PATH}/setVersion.cmake)



OPTION(RUN_TESTS "Run and compile tests"
    ON)
OPTION(USE_HDF5 "Use hdf"
       OFF)


# Define the executable name
SET(ALPHAHOUSEEXE AlphaHouse)
MESSAGE("Prog name ${ALPHAHOUSEEXE}")


# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(OBJ ${CMAKE_SOURCE_DIR}/objs)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TESTS ${CMAKE_SOURCE_DIR}/tests)
SET(SRCALPHAHOUSE ${SRC}/)

INCLUDE(${CMAKE_MODULE_PATH}/findDoxygen.cmake)

# Defined libary source

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

# Files below are fortran files that tests link to
if (${RUN_TESTS})
  # file(GLOB testDeps "${SRC}/*.f90")
  SET(testDeps
      ${SRC}/Utilities/ConstantModule.f90
      ${SRC}/Utilities/AlphaHouseMod.f90
      ${SRC}/Utilities/AlphaSystemMod.f90
      ${SRC}/BioComputational/AlphaStatMod.f90
      ${SRC}/Utilities/OrderPackModule.f90
      ${SRC}/Utilities/HashModule.f90
      ${SRC}/Utilities/LinkedListModule.f90
      ${SRC}/BioComputational/IndividualModule.f90
      ${SRC}/BioComputational/IndividualLinkedListModule.f90
      ${SRC}/BioComputational/IndividualHelperModule.f90
      ${SRC}/BioComputational/PedigreeModule.f90
      ${SRC}/BioComputational/IntelRNGMod.f90
      ${SRC}/Utilities/StringModule.f90
      ${SRC}/Utilities/LineModule.f90
      ${SRC}/Utilities/PageModule.f90
      ${SRC}/Utilities/IntegerLinkedListModule.f90
      ${SRC}/Utilities/Input.f90
      ${SRC}/BioComputational/BitUtilities.f90
      ${SRC}/BioComputational/HaplotypeModule.f90
      ${SRC}/BioComputational/GenotypeModule.f90
)

endif (${RUN_TESTS})

INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

INCLUDE(${CMAKE_MODULE_PATH}/findMKL.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache.
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)
if (${USE_HDF5})
	INCLUDE(${CMAKE_MODULE_PATH}/findHDF5.cmake)
endif(${USE_HDF5})

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${OBJ})


# TODO add tests to correct location
ADD_SUBDIRECTORY(${SRCALPHAHOUSE} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake
)

set(TESTFILES
    ${TESTS}/TestAlphaStatMod.pf
    ${TESTS}/TestOrderPackModule.pf
    ${TESTS}/TestHashTable.pf
    ${TESTS}/TestIndividual.pf
    # ${TESTS}/TestIntelRNGMod.pf
    ${TESTS}/TestLinkedList.pf
    ${TESTS}/testAlphaHouseMod.pf
    ${TESTS}/testIndividualHelperModule.pf
    ${TESTS}/testIndividualLinkedList.pf
    ${TESTS}/testInputFile.pf
    ${TESTS}/testLine.pf
    ${TESTS}/testPage.pf
    ${TESTS}/testPedigreeModule.pf
    ${TESTS}/testString.pf
    ${TESTS}/testSystemMod.pf
)

# file(GLOB TESTFILES ${TESTS}/*.pf)
if (${RUN_TESTS})
	INCLUDE(${CMAKE_MODULE_PATH}/findPFUnit.cmake)
endif(${RUN_TESTS})

