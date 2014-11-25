.DEFAULT_GOAL:=build

# General vars
FC:=ifort
MKLROOT:=/opt/intel/mkl
bin:=../AlphaHouseBin
optdir:=-I. -I${bin} -I${MKLROOT}/include
#mkl=-L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
opt:=

# Debug flags
ifeq (${debug}, true)
  opt:=${opt} -g -traceback -debug all -check all -ftrapuv -warn all
endif

# List of module directories (should only need to edit this!!!)
mod:=GeneralPurpose IntelRNG Miscellaneous ParameterFile Pedigree # ThirdPartyRoutines

# Get various stuff
src:=$(addsuffix Mod.f90,${mod})
obj:=$(addprefix ${bin}/,$(addsuffix Mod.o,${mod}))
test:=$(addsuffix /Test,${mod}) # TODO: this does not work for TPR!!! Should we reorganize TPR as other modules?

all: build test doc

# --- AlphaHouse library ---

build: .buildecho ${bin}/AlphaHouse.a # Build the library
	@echo "\n * AlphaHouseBin DONE\n"

.buildecho:
	@echo "\n * Building...\n"

${bin}/AlphaHouse.a: ${obj} # Build library
	ar cr ${bin}/AlphaHouse.a ${obj};

# --- Modules ---

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

define make-module

${1}: ${bin}/${1}Mod.o # Build module ${1}

endef

define make-target

# Required modules for the target
reqMod${1}:=$(strip $(addprefix ${bin}/,$(addsuffix .o,$(shell grep -i "use " ${1}/${1}Mod_Start.f90 | awk '{ print $$2 }'))))
# ...remove false positives from the above
reqMod${1}:=$(filter-out,$(filter-out,${mod},${reqMod${1}}),${reqMod${1}})

# Source files of the target
reqSrc${1}:=$(strip $(call rwildcard,,${1}/*.f90))

# Make
${bin}/${1}Mod.o: $${reqMod${1}} $${reqSrc${1}} # Go into module folder, collate all the code in one file, and compile that file
	@echo "\n * Go into module folder ${1}, collate all the code in one file, and compile that file...\n"
	$${MAKE} -C ${1}/;
	$${FC} $${opt} $${optdir} -c ${1}/${1}Mod.f90 -o $${bin}/${1}Mod.o -module $${bin}/

endef

#$(info $(foreach i,${mod},$(call make-module,${i})))
$(foreach i,${mod},$(eval $(call make-module,${i})))
#$(info $(foreach i,${mod},$(call make-target,${i})))
$(foreach i,${mod},$(eval $(call make-target,${i})))

# --- Test ---

test: .testecho ${obj} ${test} # Unit testing - main target
	@echo "\n * Testing DONE \n"

.testecho:
	@echo "\n * Testing... \n"

${test}: # Unit testing - individual targets
	@echo "\n $@ \n"
	${MAKE} -C $@

# --- Documentation ---

doc: docsrc docbin # Create documentation
	@echo "\n * Documentation DONE \n"

docsrc: # Create documentation with the source in this folder
	@echo "\n * Create documentation with the source in this folder...\n";
	mkdir -p DoxygenDoc;
	cat ../Doxygen.txt | sed -e "s|PROJECT_NAME=\"\"|PROJECT_NAME=\"AlphaHouse\"|"\
														-e "s|PROJECT_BRIEF=\"\"|PROJECT_BRIEF=\"A set of housekeeping routines for the Alpha programs\"|"\
														-e "s|SOURCE_BROWSER=NO|SOURCE_BROWSER=YES|" > Doxygen.tmp;
	doxygen Doxygen.tmp > DoxygenDoc/Doxygen.log;
	rm -f Doxygen.tmp;
	cd DoxygenDoc && ln -sf html/index.html .

docbin: # Create documentation without the source in the binary folder
	@echo "\n * Create documentation without the source in the binary folder...\n";
	mkdir -p ${bin}/DoxygenDoc;
	cat ../Doxygen.txt | sed -e "s|PROJECT_NAME=\"\"|PROJECT_NAME=\"AlphaHouse\"|"\
														-e "s|PROJECT_BRIEF=\"\"|PROJECT_BRIEF=\"A set of housekeeping routines for the Alpha programs\"|"\
														-e "s|OUTPUT_DIRECTORY=DoxygenDoc|OUTPUT_DIRECTORY=${bin}/DoxygenDoc|" > Doxygen.tmp;
	doxygen Doxygen.tmp > ${bin}/DoxygenDoc/Doxygen.log;
	rm -f Doxygen.tmp;
	cd ${bin}/DoxygenDoc && ln -sf html/index.html .

# --- Cleanup ---

clean: cleanbin cleansrc # Remove auto-generated source files and compiled files
	@echo "\n * Cleanup DONE \n"

cleanall: cleanbin cleansrc cleantest cleandoc # Remove all auto-generated and compiled files and other files
	@echo "\n * Cleanup-all DONE \n"

cleanbin: # Remove object, module, and library files in the binary folder
	@echo "\n * Remove object, module, and library files in the binary folder...\n"
	rm -f ${bin}/*.{o,mod,a}

cleansrc: # Remove auto-generated source files
	@echo "\n * Remove auto-generated source files...\n"
	rm -f $(addprefix */,${src})

cleantest:=$(addsuffix .clean,${test})
cleantest: .cleantestecho ${cleantest} # Remove test output files - main target

.cleantestecho: # Remove test output files - echo
	@echo "\n * Remove test output files...\n"

${cleantest}: # Remove test output files - individual targets
	${MAKE} -C $(basename $@) clean

cleandoc: cleandocsrc cleandocbin # Remove auto-generated documentation

cleandocsrc: # Remove auto-generated documentation with the source in this folder
	@echo "\n * Remove auto-generated documentation with the source in this folder...\n";
	rm -Rf DoxygenDoc

cleandocbin: # Remove auto-generated documentation without the source in the binary folder
	@echo "\n * Remove auto-generated documentation without the source in the binary folder...\n";
	rm -Rf ${bin}/DoxygenDoc

# --- Utilities ---

help: # Help
	@echo '\nTarget: Dependency # Description';
	@echo '==================================================';
	@egrep -e '^[[:alnum:]+_()%]*: ' Makefile

.PHONY: all build test doc clean help ${mod} ${test}
