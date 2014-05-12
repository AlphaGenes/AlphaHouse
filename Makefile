comp=ifort
opt=
bin=../AlphaHouseBin

all: AlphaModule AlphaRoutine
	ar cr ${bin}/AlphaHouse.a ${bin}/*.o;
	@echo "AlphaHouseBin: DONE"

# --- AlphaModule ---

AlphaModule: ${bin}/GeneralPurposeMod.o ${bin}/PedigreeMod.o

${bin}/GeneralPurposeMod.o: AlphaModule/GeneralPurposeMod.f90
	${comp} ${opt} -c AlphaModule/GeneralPurposeMod.f90 -o ${bin}/GeneralPurposeMod.o -module ${bin}/

${bin}/PedigreeMod.o: AlphaModule/PedigreeMod.f90
	${comp} ${opt} -c AlphaModule/PedigreeMod.f90 -o ${bin}/PedigreeMod.o -module ${bin}/

# --- AlphaRoutine ---

AlphaRoutine: ${bin}/MiscellaneousModL1.o ${bin}/ParameterFileModL1.o ${bin}/ParameterFileModL2.o

${bin}/MiscellaneousModL1.o: AlphaRoutine/Miscellaneous/Level1/*.f90
	$(MAKE) -C AlphaRoutine/Miscellaneous/Level1/;
	${comp} ${opt} -I${bin} -c AlphaRoutine/Miscellaneous/Level1/MiscellaneousModL1.f90 -o ${bin}/MiscellaneousModL1.o -module ${bin}/

${bin}/ParameterFileModL1.o: AlphaRoutine/ParameterFile/Level1/*.f90
	$(MAKE) -C AlphaRoutine/ParameterFile/Level1;
	${comp} ${opt} -I${bin} -c AlphaRoutine/ParameterFile/Level1/ParameterFileModL1.f90 -o ${bin}/ParameterFileModL1.o -module ${bin}/

${bin}/ParameterFileModL2.o: AlphaRoutine/ParameterFile/Level2/*.f90
	$(MAKE) -C AlphaRoutine/ParameterFile/Level2;
	${comp} ${opt} -I${bin} -c AlphaRoutine/ParameterFile/Level2/ParameterFileModL2.f90 -o ${bin}/ParameterFileModL2.o -module ${bin}/

clean: # Cleanup object and module files in binary folder
	rm -f ${bin}/*.{a,o,mod}

help: # Help
	@echo '\nTarget: Dependency # Description'; \
	echo '=================================================='; \
	egrep '^[[:alnum:].+_()%]*:' Makefile
