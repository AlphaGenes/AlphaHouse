comp=ifort
opt=
bin=../AlphaHouseBin

all: ${bin}/AlphaHouse.a
	@echo "AlphaHouseBin: DONE"

# --- AlphaHouse library ---

${bin}/AlphaHouse.a: ${bin}/MiscellaneousMod.o ${bin}/ParameterFileMod.o ${bin}/GeneralPurposeMod.o ${bin}/PedigreeMod.o
	ar cr ${bin}/AlphaHouse.a ${bin}/*.o;

# --- AlphaModules ---

${bin}/GeneralPurposeMod.o: GeneralPurpose/*.f90
	$(MAKE) -C GeneralPurpose/;
	${comp} ${opt} -c GeneralPurpose/GeneralPurposeMod.f90 -o ${bin}/GeneralPurposeMod.o -module ${bin}/

${bin}/MiscellaneousMod.o: Miscellaneous/*.f90
	$(MAKE) -C Miscellaneous/;
	${comp} ${opt} -I${bin} -c Miscellaneous/MiscellaneousMod.f90 -o ${bin}/MiscellaneousMod.o -module ${bin}/

${bin}/ParameterFileMod.o: ParameterFile/*.f90
	$(MAKE) -C ParameterFile;
	${comp} ${opt} -I${bin} -c ParameterFile/ParameterFileMod.f90 -o ${bin}/ParameterFileMod.o -module ${bin}/

${bin}/PedigreeMod.o: Pedigree/*.f90
	$(MAKE) -C Pedigree/;
	${comp} ${opt} -c Pedigree/PedigreeMod.f90 -o ${bin}/PedigreeMod.o -module ${bin}/

# --- Utilities ---

clean: # Cleanup object and module files in binary folder
	rm -f ${bin}/*.{a,o,mod}

help: # Help
	@echo '\nTarget: Dependency # Description'; \
	echo '=================================================='; \
	egrep '^[[:alnum:].+_()%]*:' Makefile
