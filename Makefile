comp:=ifort
opt:=
bin:=../AlphaHouseBin

# List of module directories (should only need to edit this!!!)
dirs:=GeneralPurpose Miscellaneous ParameterFile Pedigree ThirdPartyRoutines

# Get folder names, object file names, and target names
src:=$(addsuffix Mod.f90,${dirs})
obj:=$(addsuffix Mod.o,${dirs})
binobj:=$(addprefix ${bin}/,${obj})

all: ${bin}/AlphaHouse.a # Build the library
	@echo "\n * AlphaHouseBin: DONE\n"

# --- AlphaHouse library ---

${bin}/AlphaHouse.a: ${binobj} # Build library
	@echo "\n * Build library...\n"
	ar cr ${bin}/AlphaHouse.a ${binobj};

# --- Modules ---

rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))

define make-target

# Required modules
reqMod$(1):=$(strip $(addprefix ${bin}/,$(addsuffix .o,$(shell grep -i "use " $(1)/$(1)Mod_Start.f90 | awk '{ print $$2 }'))))
# Source files
reqSrc$(1):=$(strip $(call rwildcard,,$(1)/*.f90))
${bin}/$(1)Mod.o: $${reqMod$(1)} $${reqSrc$(1)} # Go into module folder, collate all the code in one file, and compile that file
	@echo "\n * Go into module folder $(1), collate all the code in one file, and compile that file...\n"
	$$(MAKE) -C $(1)/;
	$${comp} $${opt} -I$${bin} -c $(1)/$(1)Mod.f90 -o $${bin}/$(1)Mod.o -module $${bin}/

endef

#$(info $(foreach i,${dirs},$(call make-target,${i})))
$(foreach i,${dirs},$(eval $(call make-target,${i})))

# --- Documentation ---

doc: docsrc docbin # Create documentation

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

cleanall: cleanbin cleansrc cleandoc # Cleanup everything

cleanbin: # Remove object, module, and library files in the binary folder
	@echo "\n * Remove object, module, and library files in the binary folder...\n"
	rm -f ${bin}/*.{o,mod,a}

cleansrc: # Remove auto-generated source files
	@echo "\n * Remove auto-generated source files...\n"
	rm -f $(addprefix */,${src})

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
	@egrep -e '^[[:alnum:].+_()%]*:' -e '^\$$' Makefile
