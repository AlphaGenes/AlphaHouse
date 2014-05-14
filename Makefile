comp := ifort
opt :=
bin := ../AlphaHouseBin

# List of module directories (should only need to edit this!!!)
dirs := GeneralPurpose Miscellaneous ParameterFile Pedigree ThirdPartyRoutines

# Get folder names, object file names, and target names
src := $(addsuffix Mod.f90,${dirs})
obj := $(addsuffix Mod.o,${dirs})
binobj := $(addprefix ${bin}/,${obj})

all: ${bin}/AlphaHouse.a # Build the library
	@echo "AlphaHouseBin: DONE"

# --- AlphaHouse library ---

${bin}/AlphaHouse.a: ${binobj} # Build library
	@echo "\n * Build library...\n"
	ar cr ${bin}/AlphaHouse.a ${binobj};

# --- Modules ---

define make-target
${bin}/$(1)Mod.o: $(1)/*.f90 # Go into module folder, collate all the code in one file, and compile that file
	@echo "\n * Go into module folder $(1), collate all the code in one file, and compile that file...\n"
	$$(MAKE) -C $(1)/;
	$${comp} $${opt} -c $(1)/$(1)Mod.f90 -o $${bin}/$(1)Mod.o -module $${bin}/
endef

# $(info $(foreach i,${dirs},$(call make-target,${i})))
$(foreach i,${dirs},$(eval $(call make-target,${i})))

#${bin}/%.o: %/*.f90 # Go into module folder, collate all the code in one file, and compile that file
#	$(MAKE) -C $*/;
#	${comp} ${opt} -c $*/$*Mod.f90 -o $@ -module ${bin}/

#${bin}/GeneralPurposeMod.o: GeneralPurpose/*.f90
#	$(MAKE) -C GeneralPurpose/;
#	${comp} ${opt} -c GeneralPurpose/GeneralPurposeMod.f90 -o ${bin}/GeneralPurposeMod.o -module ${bin}/
#
#${bin}/MiscellaneousMod.o: Miscellaneous/*.f90
#	$(MAKE) -C Miscellaneous/;
#	${comp} ${opt} -I${bin} -c Miscellaneous/MiscellaneousMod.f90 -o ${bin}/MiscellaneousMod.o -module ${bin}/
#
#${bin}/ParameterFileMod.o: ParameterFile/*.f90
#	$(MAKE) -C ParameterFile;
#	${comp} ${opt} -I${bin} -c ParameterFile/ParameterFileMod.f90 -o ${bin}/ParameterFileMod.o -module ${bin}/
#
#${bin}/PedigreeMod.o: Pedigree/*.f90
#	$(MAKE) -C Pedigree/;
#	${comp} ${opt} -c Pedigree/PedigreeMod.f90 -o ${bin}/PedigreeMod.o -module ${bin}/
#
#${bin}/ThirdPartyRoutinesMod.o: ThirdPartyRoutines/*.f90
#	$(MAKE) -C ThirdPartyRoutines/;
#	${comp} ${opt} -c ThirdPartyRoutines/ThirdPartyRoutinesMod.f90 -o ${bin}/ThirdPartyRoutinesMod.o -module ${bin}/

# --- Documentation ---

doc: docsrc docbin # Create documentation
	@echo "\n * Create documentation...\n"

docsrc: # Create documentation with the source
	@echo "\n * Create documentation with the source...\n";
	mkdir -p DoxygenDoc;
	cat Doxygen.txt | sed -e "s/SOURCE_BROWSER=NO/SOURCE_BROWSER=YES/" > Doxygen.tmp;
	doxygen Doxygen.tmp > DoxygenDoc/Doxygen.log;
	rm -f Doxygen.tmp;
	cd DoxygenDoc && ln -sf html/index.html .

docbin: # Create documentation without the source for the binary folder
	@echo "\n * Create documentation without the source for the binary folder...\n";
	mkdir -p ${bin}/DoxygenDoc;
	cat Doxygen.txt | sed -e "s|OUTPUT_DIRECTORY=DoxygenDoc|OUTPUT_DIRECTORY=${bin}/DoxygenDoc|" > Doxygen.tmp;
	doxygen Doxygen.tmp > DoxygenDoc/Doxygen.log;
	rm -f Doxygen.tmp;
	cd ${bin}/DoxygenDoc && ln -sf html/index.html .

# --- Cleanup ---

cleanall: cleanbin cleansrc cleandoc # Cleanup everything
	@echo "\n * Remove...\n"

cleanbin: # Remove object, module, and library files in the binary folder
	@echo "\n * Remove object, module, and library files in the binary folder...\n"
	rm -f ${bin}/*.{o,mod,a}

cleansrc: # Remove auto-generated source files
	@echo "\n * Remove auto-generated source files...\n"
	rm -f $(addprefix */,${src})

cleandoc: cleandocsrc cleandocbin # Remove auto-generated documentation
	@echo "\n * Remove auto-generated documentation...\n"

cleandocsrc: # Remove auto-generated documentation with the source
	@echo "\n * Remove auto-generated documentation with the source...\n";
	rm -Rf DoxygenDoc

cleandocbin: # Remove auto-generated documentation without the source for the binary folder
	@echo "\n * Remove auto-generated documentation without the source for the binary folder...\n";
	rm -Rf ${bin}/DoxygenDoc

help: # Help
	@echo '\nTarget: Dependency # Description';
	@echo '==================================================';
	@egrep -e '^[[:alnum:].+_()%]*:' -e '^\$$' Makefile
