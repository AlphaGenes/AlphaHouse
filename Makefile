.DEFAULT_GOAL:=build

# General vars
prog:=AlphaHouse
progDesc:=A set of housekeeping routines for the Alpha programs

FC:=ifort
bin:=../AlphaHouseBin

# MKL
include ../Makefile.MKLRoot
incdir:=${MKLINC}
libdir:=${MKLLIB}

# Options
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

build: .buildecho ${bin}/AlphaHouse.a # Build library
	@printf "\n * ${prog} build DONE\n"

.buildecho:
	@printf "\n * ${prog} build...\n"

${bin}/AlphaHouse.a: ${obj} # Build library
	@printf "\n * ${prog} create library...\n"
	ar cr ${bin}/AlphaHouse.a ${obj};

# --- Modules ---

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

define make-module

${1}: ${bin}/${1}Mod.o # Build module ${1}

endef

define make-target

# Required modules for the target
reqMod${1}:=$(strip $(addprefix ${bin}/,$(addsuffix .o,$(shell grep -i "use " ${1}/${1}Mod.f90 | awk '{ print $$2 }'))))
# ...remove false positives from the above
reqMod${1}:=$(filter-out,$(filter-out,${mod},${reqMod${1}}),${reqMod${1}})

# Source files of the target
reqSrc${1}:=$(strip $(call rwildcard,,${1}/*.f90))

# Make
${bin}/${1}Mod.o: Makefile $${reqMod${1}} $${reqSrc${1}} # Go into module folder and compile, compile, and put object and module files into ${bin}
	@printf "\n * Go into module folder ${1}, compile, and put object and module files into $${bin}...\n"
	$${FC} -c ${1}/${1}Mod.f90 -o $${bin}/${1}Mod.o -module $${bin}/ $${opt} $${incdir}

endef

#$(info $(foreach i,${mod},$(call make-module,${i})))
$(foreach i,${mod},$(eval $(call make-module,${i})))
#$(info $(foreach i,${mod},$(call make-target,${i})))
$(foreach i,${mod},$(eval $(call make-target,${i})))

# --- Test ---

test: .testecho ${obj} ${test} # Unit testing - main target
	@printf "\n * ${prog} test DONE \n"

.testecho:
	@printf "\n * ${prog} test... \n"

${test}: # Unit testing - individual targets
	@printf "\n $@ \n"; \
	if [ -d $@ ]; then \
	  ${MAKE} -C $@; \
	else \
	  printf "\n * ${@} test (no Test folder)...\n" ;\
	fi

# --- Documentation ---

doc: docsrc docbin # Create developer documentation
	@printf "\n * ${prog} create developer documentation DONE \n"

docsrc: # Create developer documentation in this folder
	@printf "\n * ${prog} create developer documentation in this folder...\n"; \
	rm -Rf DoxygenDoc; \
	mkdir -p DoxygenDoc; \
	cat ../Doxygen.txt | sed -e "s|PROJECT_NAME=\"\"|PROJECT_NAME=\"${prog}\"|" \
														-e "s|PROJECT_BRIEF=\"\"|PROJECT_BRIEF=\"${progDesc}\"|" \
														-e "s|FILE_PATTERNS=|FILE_PATTERNS=*.f90|" \
														-e "s|SOURCE_BROWSER=NO|SOURCE_BROWSER=YES|" > Doxygen.tmp; \
	doxygen Doxygen.tmp > DoxygenDoc/Doxygen.log; \
	rm -f Doxygen.tmp; \
	cd DoxygenDoc && ln -sf html/index.html .

docbin: # Create developer documentation in the binary folder
	@printf "\n * ${prog} create developer documentation in the binary folder...\n"; \
	rm -Rf ${bin}/DoxygenDoc; \
	mkdir -p ${bin}/DoxygenDoc; \
	cat ../Doxygen.txt | sed -e "s|PROJECT_NAME=\"\"|PROJECT_NAME=\"${prog}\"|" \
														-e "s|PROJECT_BRIEF=\"\"|PROJECT_BRIEF=\"${progDesc}\"|" \
														-e "s|FILE_PATTERNS=|FILE_PATTERNS=*.f90|" \
														-e "s|OUTPUT_DIRECTORY=DoxygenDoc|OUTPUT_DIRECTORY=${bin}/DoxygenDoc|" > Doxygen.tmp; \
	doxygen Doxygen.tmp > ${bin}/DoxygenDoc/Doxygen.log; \
	rm -f Doxygen.tmp; \
	cd ${bin}/DoxygenDoc && ln -sf html/index.html .

# --- Cleanup ---

clean: cleanbin # Cleanup compiled files
	rm -f */*__genmod.{f90,mod}
	@printf "\n * ${prog} cleanup DONE \n"

cleanall: cleanbin cleantest cleandoc # Cleanup compiled files and other files
	@printf "\n * ${prog} cleanup-all DONE \n"

cleanbin: # Cleanup object, module, and library files in the binary folder
	@printf "\n * ${prog} cleanup object, module, and library files in the binary folder...\n"
	rm -f ${bin}/*.{o,mod,a}

cleantest:=$(addsuffix .clean,${test})
cleantest: .cleantestecho ${cleantest} # Cleanup test output files - main target

.cleantestecho: # Cleanup test output files - echo
	@printf "\n * ${prog} cleanup test output files...\n"

${cleantest}: # Cleanup test output files - individual targets
	@if [ -d $(basename $@) ]; then \
	  ${MAKE} -C $(basename $@) clean; \
	else \
	  printf "\n * $(basename $@) (no Test folder)...\n" ;\
	fi

cleandoc: cleandocsrc cleandocbin # Remove auto-generated developer documentation

cleandocsrc: # Cleanup developer documentation in this folder
	@printf "\n * ${prog} cleanup developer documentation in this folder...\n";
	rm -Rf DoxygenDoc

cleandocbin: # Cleanup developer documentation in the binary folder
	@printf "\n * ${prog} cleanup developer documentation in the binary folder...\n";
	rm -Rf ${bin}/DoxygenDoc

# --- Utilities ---

help: # Help
	@printf '\nTarget: Dependency # Description'; \
	printf '=================================================='; \
	egrep -e '^[[:alnum:]+_()%]*: ' Makefile

.PHONY: all build test doc clean help ${mod} ${test}
