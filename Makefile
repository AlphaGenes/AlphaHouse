comp=ifort
opt=
bin=../AlphaHouseBin

AlphaHousePath=../AlphaHouseBin
AlphaHouseVer=
AlphaHouse=$(AlphaHousePath)/AlphaHouse$(AlphaHouseVer).a

all: GeneralPurpose.o Pedigree.o GetRequiredCommandFileInfo.o CheckKeyWordPresence.o ReadParamPedigree.o # Compile everything
	ar cr $(bin)/AlphaHouse.a $(bin)/GeneralPurpose.o $(bin)/Pedigree.o $(bin)/GetRequiredCommandFileInfo.o $(bin)/CheckKeyWordPresence.o $(bin)/ReadParamPedigree.o 

# AlphaModule

GeneralPurpose.o: AlphaModule/GeneralPurpose.f90
	$(comp) $(opt) -c AlphaModule/GeneralPurpose.f90 -o $(bin)/GeneralPurpose.o -module $(bin)

Pedigree.o: AlphaModule/Pedigree.f90
	$(comp) $(opt) -c AlphaModule/Pedigree.f90 -o $(bin)/Pedigree.o -module $(bin)

# AlphaRoutine

ReadParamPedigree.o: AlphaRoutine/ParameterFile/ReadParamPedigree.f90 Pedigree.o
	$(comp) $(opt) -I$(AlphaHousePath) -c AlphaRoutine/ParameterFile/ReadParamPedigree.f90 -o $(bin)/ReadParamPedigree.o

CheckKeyWordPresence.o: AlphaRoutine/ParameterFile/CheckKeyWordPresence.f90
	$(comp) $(opt) -I$(AlphaHousePath) -c AlphaRoutine/ParameterFile/CheckKeyWordPresence.f90 -o $(bin)/CheckKeyWordPresence.o

GetRequiredCommandFileInfo.o: AlphaRoutine/ParameterFile/GetRequiredCommandFileInfo.f90
	$(comp) $(opt) -I$(AlphaHousePath) -c AlphaRoutine/ParameterFile/GetRequiredCommandFileInfo.f90 -o $(bin)/GetRequiredCommandFileInfo.o

clean: # Cleanup in current folder and in $(bin)
	rm -f *.{o,mod,lst}
	rm -f $(bin)/*.{o,mod,lst}

%: # Generic make compile, i.e., use 'make dir/myprogram' to compile 'ifort -c dir/myprogram.f90 -o $(bin)/myprogram.o'
	tmp=$(basename $@)
	$(comp) $(opt) -c $*.f90 -o $(bin)/$(tmp).o -module $(bin)

help: # Help
	@echo '\nTarget: Dependency # Description'; \
	echo '=================================================='; \
	egrep '^[[:alnum:].+_()%]*:' Makefile
