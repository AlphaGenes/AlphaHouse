comp=ifort
opt=
bin=../AlphaHouseBin

all: AlphaModule AlphaRoutine
	ar cr ${bin}/AlphaHouse.a ${bin}/*.o;
	@echo "AlphaHouseBin success"

# --- AlphaModule ---

AlphaModule: GeneralPurposeMod.o PedigreeMod.o

GeneralPurposeMod.o: AlphaModule/GeneralPurposeMod.f90
	${comp} ${opt} -c AlphaModule/GeneralPurposeMod.f90 -o ${bin}/GeneralPurposeMod.o -module ${bin}

PedigreeMod.o: AlphaModule/PedigreeMod.f90
	${comp} ${opt} -c AlphaModule/PedigreeMod.f90 -o ${bin}/PedigreeMod.o -module ${bin}

# --- AlphaRoutine ---

AlphaRoutine: MiscellaneousModL1.o ParameterFileModL1.o ParameterFileModL2.o

MiscellaneousModL1.o: AlphaRoutine/Miscellaneous/Level1/MiscellaneousModL1.f90
	${comp} ${opt} -I${bin} -c AlphaRoutine/Miscellaneous/Level1/MiscellaneousModL1.f90 -o ${bin}/MiscellaneousModL1.o -module ${bin}

AlphaRoutine/Miscellaneous/Level1/MiscellaneousModL1.f90:
	cd AlphaRoutine/Miscellaneous/Level1
	make

ParameterFileModL1.o: AlphaRoutine/ParameterFile/Level1/ParameterFileModL1.f90
	${comp} ${opt} -I${bin} -c AlphaRoutine/ParameterFile/Level1/ParameterFileModL1.f90 -o ${bin}/ParameterFileModL1.o -module ${bin}

AlphaRoutine/ParameterFile/Level1/ParameterFileModL1.f90:
	cd AlphaRoutine/ParameterFile/Level1
	make

ParameterFileModL2.o: AlphaRoutine/ParameterFile/Level2/ParameterFileModL2.f90
	${comp} ${opt} -I${bin} -c AlphaRoutine/ParameterFile/Level2/ParameterFileModL2.f90 -o ${bin}/ParameterFileModL2.o -module ${bin}

AlphaRoutine/ParameterFile/Level2/ParameterFileModL2.f90:
	cd AlphaRoutine/ParameterFile/Level2
	make

clean: # Cleanup in current folder and in ${bin}
	rm -f ${bin}/*.{a,o,mod,lst}

help: # Help
	@echo '\nTarget: Dependency # Description'; \
	echo '=================================================='; \
	egrep '^[[:alnum:].+_()%]*:' Makefile
