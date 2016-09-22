#Modules
_MODULE_SOURCES=GenParam.f90 \
								PedigreeTable.f90 \
                UtilitySubroutines.f90 \
								macs.f90
								

_PROGRAM_SOURCES=main.f90

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
	VERSION:=$(shell git rev-parse --short HEAD)
	SUBVERSION:=0
	#MASTERVERSION:=$(shell git describe --tag | cut -d "-" -f 1)
	PROGRAM:=$(NAME)_$(VERSION)#$(MASTERVERSION)
	RM=rm -rf
  
	MAKEDIR=mkdir -p

  LIBNETCDF := -L$(NETCDF)/lib -lnetcdff
	INCNETCDF := -I$(NETCDF)/include 
  PFUNIT=/usr/local/pFUnit_serial
endif
#Filenames
OBJECTDIR=objs/
SRCDIR=src/
TESTDIR=tests/
TARGETDIR=bin/


PFUNITFLAGS:=-I./$(OBJECTDIR) -I$(PFUNIT)/mod -I$(PFUNIT)/include -L$(PFUNIT)/lib -lpfunit -module $(OBJECTDIR) -fpp 
HDF5FLAGS:=
#Compiler and flags
CC:=g++
CFLAGS:= -std=c++0x -Wextra -Wall 

FC=ifort

FFLAGS= -fpp -module $(OBJECTDIR) -D $(OSFLAG) -DVERS=""commit-$(VERSION)""
LDFLAGS:=-I$(OBJECTDIR)
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

debug: FFLAGS:=$(FFLAGS) $(DEBUG_FLAGS)
debug: all 

superdebug: FFLAGS:=$(FFLAGS) $(SUPER_DEBUG_FLAGS)
superdebug: all 
	
production: FFLAGS:=$(FFLAGS) $(PRODUCTION_FLAG)
production: PROGRAM:=$(NAME)
production: all 

build: all tests cleanTest



tests: cleanTest $(F90TESTS) $(F90FINISHED)
	cp $(TESTDIR)testSuites.inc .
	$(FC) -o tests.x $(OTHEROBJECTFILES) $(addprefix $(OBJECTDIR), $(notdir $(F90FINISHED))) $(MODULE_OBJECTS) $(PFUNIT)/include/driver.f90 $(PFUNITFLAGS)
	$(RM) testSuites.inc
	./tests.x -xml test_output

clean: cleanIntermediate
	$(RM) $(TARGETDIR)$(NAME)_*

cleanIntermediate: cleanTest
	$(RM) $(OBJECTDIR)*

cleanTest:
	$(RM) $(TESTDIR)*.F90
	$(RM) $(OBJECTDIR)tests.x
	$(RM) $(OBJECTDIR)testSuites.inc

lib:$(MODULE_OBJECTS)
	$(FC) -o $(TARGETDIR)$(PROGRAM).a -staticlib $(MODULE_OBJECTS) $^

$(PROGRAM):$(OBJECTS)
	$(FC) -o $(TARGETDIR)$(PROGRAM) $(OTHEROBJECTFILES) $^ 





$(OBJECTDIR)%.o: $(SRCDIR)%.f90
	$(FC) -o $@ -c $< $(FFLAGS) $(LDFLAGS)

directories:
	$(MAKEDIR) $(SRCDIR)
	$(MAKEDIR) $(TESTDIR)
	$(MAKEDIR) $(OBJECTDIR)
	$(MAKEDIR) $(TARGETDIR) 

list: 	#Taken from http://stackoverflow.com/questions/4219255/how-do-you-get-the-list-of-targets-in-a-makefile. Answer by mklement0.	
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs


%.o:%.F90
	$(FC) -c -o $(addprefix $(OBJECTDIR), $(notdir $@))  $< $(PFUNITFLAGS)

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

