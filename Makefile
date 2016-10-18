#Modules
_MODULE_SOURCES=AlphaHouseMod.f90 AlphaEvolveMod.f90 AlphaStatMod.f90 IntelRNGMod.f90 OrderPackMod.f90 PedigreeTable.f90 UtilitySubroutines.f90

_PROGRAM_SOURCES=main2.f90

DEPENDENCIES=
ifeq ($(OS), WINDOWS_NT)
	OSFLAG="OSWIN"

	PROGRAM="ProgramName"

	EXE=".exe"
	RM=del

else
	OSFLAG:="OS_UNIX"
#This gets the name of the current folder and the current commit and sets the program name to that.
	FOLDERNAME:=$(shell basename $(shell pwd)) #Name of program goes here
	NAME:=$(strip $(FOLDERNAME))
	# VERSION:=$(shell git rev-parse --short HEAD)
	SUBVERSION:=0
	#MASTERVERSION:=$(shell git describe --tag | cut -d "-" -f 1)
	PROGRAM:=$(NAME)#_$(VERSION)#$(MASTERVERSION)
	RM=rm -rf

	MAKEDIR=mkdir -p

  LIBNETCDF := -L$(NETCDF)/lib -lnetcdff
	INCNETCDF := -I$(NETCDF)/include
  PFUNIT=/usr/local/pFUnit_serial

	MKLROOT := /opt/intel/mkl
	# On Eddie2
	# MKLROOT:=/exports/applications/apps/intel/ClusterStudio2013/mkl
	# On Eddie3
	# MKLROOT := /exports/applications/apps/SL7/intel/parallel_studio_xe_2016/mkl
	MKLLIB := -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,/opt/intel/lib
	MKLINC := -I$(MKLROOT)/include
endif
#Filenames
OBJECTDIR=objs/
SRCDIR=src/
TESTDIR=tests/
TARGETDIR=bin/
DOXYGENDIR=DoxygenDoc/

PFUNITFLAGS:=-I./$(OBJECTDIR) -I$(PFUNIT)/mod -I$(PFUNIT)/include -lpfunit -module $(OBJECTDIR) -fpp -L$(PFUNIT)/lib
LDFLAGS:=
DRIVER:=$(PFUNIT)/include/driver.f90

HDF5FLAGS:=
#Compiler and flags
CC:=g++
CFLAGS:= -std=c++0x -Wextra -Wall

FC=ifort

FFLAGS= -fpp -module $(OBJECTDIR) -D $(OSFLAG) -DVERS=""commit-$(VERSION)"" -mkl
DEBUG_FLAGS=-traceback -g -debug all -warn all -check bounds -check format \
		-check output_conversion -check pointers
SUPER_DEBUG_FLAGS=$(DEBUG_FLAGS) -ftrapuv -check all -gen-interfaces -warn interfaces

PRODUCTION_FLAG=-fast

MODULE_SOURCES=$(patsubst %, $(SRCDIR)%, $(_MODULE_SOURCES))

PROGRAM_SOURCES=$(patsubst %, $(SRCDIR)%, $(_PROGRAM_SOURCES))

SOURCES=$(MODULE_SOURCES) $(PROGRAM_SOURCES)

MODULE_OBJECTS =$(patsubst $(SRCDIR)%.f90, $(OBJECTDIR)%.o, $(MODULE_SOURCES))
PROGRAM_OBJECTS=$(patsubst $(SRCDIR)%.f90, $(OBJECTDIR)%.o, $(PROGRAM_SOURCES))

OBJECTS= $(MODULE_OBJECTS) $(PROGRAM_OBJECTS)

PFTESTS:=$(wildcard $(TESTDIR)*.pf)

F90TESTS:=$(PFTESTS:.pf=.F90)

F90FINISHED:= $(F90TESTS:.F90=.o)

BUILDDATE=$(shell date +%Y%m%d-%H:%M:%S)

all: directories $(PROGRAM) $(DEPENDENCIES)

build: all doc tests cleanTest



debug: FFLAGS:=$(FFLAGS) $(DEBUG_FLAGS)
debug: $(PROGRAM)

superdebug: FFLAGS:=$(FFLAGS) $(SUPER_DEBUG_FLAGS)
superdebug: all

production: FFLAGS:=$(FFLAGS) $(PRODUCTION_FLAG)
production: PROGRAM:=$(NAME)
production: all

doc:
	doxygen Doxygen.txt > $(DOXYGENDIR)Doxygen.log

tests: $(F90TESTS) $(F90FINISHED)
	cp $(TESTDIR)testSuites.inc .
	$(FC) -o tests.x $(addprefix $(OBJECTDIR), $(notdir $(F90FINISHED))) $(MODULE_OBJECTS) $(DRIVER) $(PFUNITFLAGS) $(MKLINC) $(MKLLIB) $(FFLAGS)
	$(RM) testSuites.inc
	./tests.x -xml test_output

testsDebug : FFLAGS:=$(FFLAGS) $(DEBUG_FLAGS)
testsDebug : tests

clean: cleanIntermediate
	$(RM) $(TARGETDIR)$(NAME)*

cleanIntermediate: cleanTest
	$(RM) $(OBJECTDIR)*

cleanTest:
	$(RM) $(TESTDIR)*.F90
	$(RM) $(OBJECTDIR)tests.x
	$(RM) $(OBJECTDIR)testSuites.inc
	$(RM) tests.x
	$(RM) testSuites.inc

$(PROGRAM):$(OBJECTS)
	$(FC) -o $(TARGETDIR)$(PROGRAM) $^ $(FFLAGS) $(MKLINC) $(MKLLIB) $(LDFLAGS)

$(OBJECTDIR)%.o: $(SRCDIR)%.f90
	$(FC) -o $@ -c $< $(FFLAGS) $(MKLINC) $(MKLLIB) $(LDFLAGS)

directories:
	$(MAKEDIR) $(SRCDIR)
	$(MAKEDIR) $(TESTDIR)
	$(MAKEDIR) $(OBJECTDIR)
	$(MAKEDIR) $(TARGETDIR)
	$(MAKEDIR) $(DOXYGENDIR)

list: 	#Taken from http://stackoverflow.com/questions/4219255/how-do-you-get-the-list-of-targets-in-a-makefile. Answer by mklement0.
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

#$(OBJECTDIR)%.o: %.F90
#	$(FC) -o $@ -c $< $(FFLAGS) $(PFUNITFLAGS) $(MKLLIB) $(MKLINC) $(LDFLAGS)

%.o:%.F90
	$(FC) -c -o $(addprefix $(OBJECTDIR), $(notdir $@))  $< $(PFUNITFLAGS) -module $(OBJECTDIR) $(LDFLAGS)

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

