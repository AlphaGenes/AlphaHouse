#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#DEFINE BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"

#else
#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#endif

Module macs
  use iso_fortran_env
  implicit none

type MacsInput
  integer(kind=int64)::nHaps, nChrom,EffecPopSize,ChrLengthBasesMaCS
  real(kind=real64)::MutationRateMaCS, MorgansMaCs
  character(300)::PopHistoryMaCS, MaCSOutputFile="output.txt"
  logical :: InternalMorgansMaCS, InternalChrLenBasesMaCS, InternalMutationRateMaCS, InternalEffecPopSize,RecombHotspotsOnOff
  character(:), allocatable:: AdditionalMaCSParams
  character(:), allocatable:: MacsParams
end type 

  contains

    subroutine getMacsInput(input, MacsSpecFileIn)
      use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts

      type(MacsInput), intent(out):: input

      INTEGER, PARAMETER:: MacsSpecFileID=100
      integer:: MacsIOStatus=0
      character(len=*), optional, intent(in):: MacsSpecFileIn
      character(1000)::MacsSpecFile

      character(len=300) :: first, line
      character(len=:), allocatable::dumStr
      character(len=300),dimension(:),allocatable :: second



      if (.not. present(MacsSpecFileIn)) then
        write(MacsSpecFile, "(A)") "MacsSpec.txt"
      else
        write(MacsSpecFile, "(A)") MacsSpecFileIn
      end if

      open(MacsSpecFileID, file=MacsSpecFile, action="read", status="old")

      READFILE: do while (MacsIOStatus==0)
        read(MacsSpecFileID,"(A)", IOStat=MacsIOStatus)  line
        if (len_trim(line)==0) then
            CYCLE
        end if

        call splitLineIntoTwoParts(trim(line), first, second)

        dumStr = parseToFirstWhitespace(first)
        if (first(1:1)=="#" .or. len(trim(line))==0) then
          cycle
        else
          select case(trim(dumStr))
          case("outputfile")
            if (.not. allocated(second)) then
              write(*, "(A,A)") "No output fileName found. Using default filename: ", input%MaCSOutputFile
            else
              write(input%MaCSOutputFile, "(A)") second(1)
            end if

          case("chromosomelengthbasesmacs")
            if (.not. allocated(second)) goto 57
            dumStr = parseToFirstWhitespace(second(1))
            if (trim(dumStr) == "internal") then
              input%InternalChrLenBasesMaCS = .true.
              input%ChrLengthBasesMaCS = 0 !initialise to 0 just in case user is being stupid
            else
              input%InternalChrLenBasesMaCS = .false.
              read(second(2),*, err=57) input%ChrLengthBasesMaCS
            endif
            cycle
            57                  print *,"ChromosomeLengthBasesMaCS not set properly"
            stop 57

          case("mutationratemacs")
            if (.not. allocated(second)) goto 58
            dumStr = parseToFirstWhitespace(second(1))
            if (trim(dumStr) == "internal") then
              input%InternalMutationRateMaCS = .true.
              input%MutationRateMaCS = 0 !initialise to 0 just in case user is being stupid
            else
              input%InternalMutationRateMaCS = .false.
              read(second(2),*, err=58) input%MutationRateMaCS
            endif
            cycle
            58                  print *,"MutationRateMaCS not set properly"
            stop 58

          case("effectivepopulationsizebasemacs")
            if (.not. allocated(second)) goto 59
            dumStr = parseToFirstWhitespace(second(1))
            if (trim(dumStr) == "internal") then
              input%InternalEffecPopSize = .true.
              input%EffecPopSize = 0 !initialise to 0 just in case user is being stupid
            else
              input%InternalEffecPopSize = .false.
              read(second(2),*, err=59) input%EffecPopSize
            endif
            cycle
            59                  print *,"EffectivePopulationSizeBaseMaCs not set properly"
            stop 59

          case("numberofhaplotypes")
            if (.not. allocated(second)) goto 54
            read(second(1), *, err=54) input%nHaps
            cycle
            54                  print *,"numberofhaplotypes not set properly"
            stop 54


          case("populationhistorymacs")
            if (.not. allocated(second)) then
              print *,"PopulationHistoryMaCS not set properly"
              stop 55
            endif
            input%PopHistoryMaCS = second(1)
            cycle
            55                  print*, "Modelling population history using ", trim(input%PopHistoryMaCS)

          case("recombinationhotspotsonoff")
            if (.not. allocated(second)) goto 56
            dumStr = parseToFirstWhitespace(second(1))
            if (trim(dumStr)=="off") then
              input%RecombHotspotsOnOff = .False.
            else
              input%RecombHotspotsOnOff = .True.
            endif
            cycle
            56                  print *,"recombinationhotspotsonoff not set properly"
            stop 56

          case("chromosomelengthmorgansmacs")
            if (.not. allocated(second)) goto 52
            dumStr = parseToFirstWhitespace(second(1))
            if (trim(dumStr) == "internal") then
              input%InternalMorgansMaCS = .true.
              input%MorgansMaCs = 0 !initialise to 0 just in case user is being stupid
            else
              input%InternalMorgansMaCS = .false.
              read(second(2),*, err=52) input%MorgansMaCs
            endif
            cycle
            52                  print *,"ChromosomeLengthMorgansMaCS not set properly"
            stop 52

          case default
            write(*,"(A,A)") trim(dumStr), "is not a valid input for macs"
            cycle
          end select
        end if
      end do READFILE

      end subroutine getMacsInput


    subroutine runMacs()
      implicit none
      Type(MacsInput):: input

      call getMacsInput(input)
      call RunMacsHelper(input)

    end subroutine runMacs

subroutine RunMacsHelper(input, ChromIn)

    ! use MacsCases
    use ifport  ! required for systemQQ
    implicit none

    type(MacsInput), intent(inout):: input
    integer,optional, intent(in) :: ChromIn
    integer::Chrom
!    character(*),intent(in) :: output ! Folder is full command to copy in to

    character(5) :: fileDirec,fileExtn
    character(100) :: nHapsStr
    character(len=300) :: fileString
    character(1000) :: xR,xO
    character(len=:),allocatable :: x ! allocatable so it can grow automagically

    logical :: success

    if (.not. present(ChromIn)) then
      Chrom = 1
    else
      Chrom=ChromIn
    end if
#ifdef OS_UNIX
    fileDirec = "./"
    fileExtn = " "
#endif
#ifdef OS_WIN
    fileDirec = " "
    fileExtn = ".exe "
#endif

  write(nHapsStr,'(i0)') input%nHaps

  ! Initialise cases for macs
  call InitMacsCases(input)!PopHistoryMaCS,MacsParams,MutationRateMaCS, MorgansMaCS, InternalMorgansMaCS, InternalChrLenBasesMaCS, InternalMutationRateMaCS, ChrLengthBasesMaCS, EffecPopSize, InternalEffecPopSize)

  if (input%RecombHotspotsOnOff) then
    call ReadRecombSpecForMacs(Chrom)
    xR="-R Hotspots.txt"
  else
    xR=""
  endif
  xO="1>"//trim(input%MaCSOutputFile)//" 2>debug.txt"
  x= trim(fileDirec) // "macs" // trim(fileExtn) // " " // trim(nHapsStr) // " " // &
    " " // trim(input%MacsParams) // " " // trim(xR) // " " // trim(xO)
  print*,x
  success = systemQQ(trim(x))
  if(success == .false.) print *, "failure : AlphaSim : RunMacsRewrite A"
!endif

!#ifdef OS_UNIX
!    print *, " Running AlphaFormatter"
!    success = systemQQ("./AlphaFormatter " //" && cp -f SegSites.txt FinishedMacs.txt" )
!    if(success == .false.) print *, "failure : AlphaSim : RunMacsRewrite B"
!#endif
!
!#ifdef OS_WIN
!    print *, " Running AlphaFormatter"
!    success = systemQQ(" ./AlphaFormatter.exe  "// trim(nHapsStr) //" && copy /Y SegSites.txt FinishedMacs.txt > nul" )
!    if(success == .false.) print *, "failure : AlphaSim : RunMacsRewrite C"
!#endif
    
  write(fileString, '("MaCSParametersUsed.txt")') 
  open(unit=5012, file=trim(fileString), status = "unknown")
  write(5012,*) trim(nHapsStr) // " " // trim(input%MacsParams)
  close(5012)

  write(fileString,'("MaCSParametersHumanReadable.txt")')
  open(5100, file =trim(fileString), status = "replace" , action="write")
  write(5100, "(A)") "THIS FILE IS REPLACED EVERY TIME MACS IS RUN"
  write(5100, "(A)") "N. Haplotypes, Chromosome Length(MaCS), Chromosome Length (AlphaSim),      Genome Size, &
    Effective Pop. Size,     MaCS Mutations (-t)/MaCS Recombinations (-r)"
  write(5100, "(i13,f25.4, i18 , i21, A, A, A)") input%nHaps, input%MorgansMaCS, input%ChrLengthBasesMaCs, &
    input%EffecPopSize, "     '", input%MacsParams, "'"
  close(5100)

end subroutine RunMacsHelper

subroutine InitMacsCases(input) !PopHistoryMaCS,MacsParams,MutationRateMaCS, MorgansMaCS, InternalMorgansMaCS, InternalChrLenBasesMaCS, InternalMutationRateMaCS, ChrLengthBasesMaCS, EffecPopSize, InternalEffecPopSize)
  use AlphaHouseMod, only: toLower

    implicit none
    type(MacsInput), intent(inout):: input

    character(len=700) :: dummy
    character(len=:), allocatable :: TmpC, CasePopHistMacs
    character(300) :: buffer
    integer :: size, status,k

    dummy = toLower(input%PopHistoryMaCS)
    k=1
    CasePopHistMacs=""
    do 
        if (dummy(k:k) /= " ") then
            CasePopHistMacs = dummy(1:k)
            k=k+1
        else
            exit
        endif
    enddo

    select case(CasePopHistMacs)

        case("internaltest")
            if (input%InternalEffecPopSize) input%EffecPopSize = 100
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=1.00E+8 !BP
            if (input%InternalMorgansMaCS) input%MorgansMaCS=1.00
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 2.5 * (1.0 / real(input%ChrLengthBasesMaCS)) ! Mu
            input%AdditionalMaCsParams = ""
            ! Old code : MacsParams="-t 4.00E-06 &
            !                        -r 4.00E-06"

        case("internalcattle")
            if (input%InternalEffecPopSize) input%EffecPopSize = 100
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=1.00E+8
            if (input%InternalMorgansMaCS) input%MorgansMaCS=1.00
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 2.5 * (1.0 / real(input%ChrLengthBasesMaCS)) ! Mu
            input%AdditionalMaCSParams="-eN    0.06    2.0 &
                        -eN    0.13    3.0 &
                        -eN    0.25    5.0 &
                        -eN    0.50    7.0 &
                        -eN    0.75    9.0 &
                        -eN    1.00   11.0 &
                        -eN    1.25   12.5 &
                        -eN    1.50   13.0 &
                        -eN    1.75   13.5 &
                        -eN    2.00   14.0 &
                        -eN    2.25   14.5 &
                        -eN    2.50   15.0 &
                        -eN    5.00   20.0 &
                        -eN    7.50   25.0 &
                        -eN   10.00   30.0 &
                        -eN   12.50   35.0 &
                        -eN   15.00   40.0 &
                        -eN   17.50   45.0 &
                        -eN   20.00   50.0 &
                        -eN   22.50   55.0 &
                        -eN   25.00   60.0 &
                        -eN   50.00   70.0 &
                        -eN  100.00   80.0 &
                        -eN  150.00   90.0 &
                        -eN  200.00  100.0 &
                        -eN  250.00  120.0 &
                        -eN  500.00  200.0 &
                        -eN 1000.00  400.0 &
                        -eN 1500.00  600.0 &
                        -eN 2000.00  800.0 &
                        -eN 2500.00 1000.0"            

        case ("internalpig") !TODO - check with GG
            if (input%InternalEffecPopSize) input%EffecPopSize = 100 
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=6.75+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=1.71
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 1.6 * (1.0 / real(input%ChrLengthBasesMaCS)) 
            input%AdditionalMaCSParams = "-eN  25.00    100.0 &
                                    -eN  50.00    200.0 &
                                    -eN  75.00    300.0 &
                                    -eN 100.00    400.0 &
                                    -eN 125.00    500.0 &
                                    -eN 150.00    600.0 &
                                    -eN 175.00    700.0 &
                                    -eN 200.00    800.0 &
                                    -eN 225.00    900.0 &
                                    -eN 250.00   1000.0 &
                                    -eN 275.00   2000.0 &
                                    -eN 300.00   3000.0 &
                                    -eN 325.00   4000.0 &
                                    -eN 350.00   5000.0 &
                                    -eN 375.00   6000.0 &
                                    -eN 400.00   7000.0 &
                                    -eN 425.00   8000.0 &
                                    -eN 450.00   9000.0 &
                                    -eN 475.00  10000.0"            
            print*,"TODO: need to find Ne values in the past for pig!!!"
            stop

        case ("internalchicken") !TODO - check with GG
            if (input%InternalEffecPopSize) input%EffecPopSize = 70                     
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=3.0E+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=0.84
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 2.5 * (1.0 / real(input%ChrLengthBasesMaCS))
            input%AdditionalMaCSParams = "-eN 0.18    0.71 &
                                    -eN 0.36    1.43 &
                                    -eN 0.54    2.14 &
                                    -eN 0.71    2.86 &
                                    -eN 0.89    3.57 &
                                    -eN 1.07    4.29 &
                                    -eN 1.25    5.00 &
                                    -eN 1.43    5.71"            
            print*,"TODO: need to find Ne values in the past for chicken!!!"
            stop

        case("internalrabbit")
            if (input%InternalEffecPopSize) input%EffecPopSize = 100                     
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=1.59E+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=1.36
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 1.74 * (1.0 / real(input%ChrLengthBasesMaCS))
            input%AdditionalMaCSParams = "-eN 0.05    1.25 &
                        -eN 0.08    1.50 &
                        -eN 0.10    1.75 &
                        -eN 0.13    2.00 &
                        -eN 0.15    2.25 &
                        -eN 0.18    2.50 &
                        -eN 0.20    2.75 &
                        -eN 0.23    3.00 &
                        -eN 0.25    3.25 &
                        -eN 0.50    4.00 &
                        -eN 1.00    5.00 &
                        -eN 1.50    6.00 &
                        -eN 2.00    7.00 &
                        -eN 2.50    8.00 &
                        -eN 3.00   90.00 &
                        -eN 3.50   10.00 &
                        -eN 4.00   11.00 &
                        -eN 4.50   12.00 &
                        -eN 5.00 1000.00"            
            ! MacsParams="-t 3.43E-06      &
            !             -r 3.43E-06      &"

        case ("internalmaize")
            if (input%InternalEffecPopSize) input%EffecPopSize = 100                     
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=2.00E+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=2.00
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 2.5 * (1.0 / real(input%ChrLengthBasesMaCS)) 
            input%AdditionalMaCSParams="-eN 0.03   1 &
                        -eN 0.05   2 &
                        -eN 0.10   4 &
                        -eN 0.15   6 &
                        -eN 0.20   8 &
                        -eN 0.25  10 &
                        -eN 0.30  12 &
                        -eN 0.35  14 &
                        -eN 0.40  16 &
                        -eN 0.45  18 &
                        -eN 0.50  20 &
                        -eN 2.00  40 &
                        -eN 3.00  60 &
                        -eN 4.00  80 &
                        -eN 5.00 100"            
            ! MacsParams="-t 2.00E-06  &
            !             -r 2.00E-06  &

        case ("internalmaizelandrace")
            if (input%InternalEffecPopSize) input%EffecPopSize = 100                     
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=2.00E+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=2.00
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 2.5 * (1.0 / real(input%ChrLengthBasesMaCS)) 
            input%AdditionalMaCSParams = "-eN  0.03    1 &
                        -eN  0.05    2 &
                        -eN  0.10    4 &
                        -eN  0.15    6 &
                        -eN  0.20    8 &
                        -eN  0.25   10 &
                        -eN  0.30   12 &
                        -eN  0.35   14 &
                        -eN  0.40   16 &
                        -eN  0.45   18 &
                        -eN  0.50   20 &
                        -eN  2.00   40 &
                        -eN  3.00   60 &
                        -eN  4.00   80 &
                        -eN  5.00  100 &
                        -eN  6.00  120 &
                        -eN  7.00  140 &
                        -eN  8.00  160 &
                        -eN  9.00  180 &
                        -eN 10.00  200 &
                        -eN 12.50  400 &
                        -eN 15.00  600 &
                        -eN 17.50  800 &
                        -eN 20.00 1000 &
                        -eN 22.50 1200 &
                        -eN 25.00 1400 &
                        -eN 27.50 1600 &
                        -eN 30.00 2000"            
            ! MacsParams="-t 2.00E-06    &
            !             -r 2.00E-06    &

        case ("internalwheat")
            if (input%InternalEffecPopSize) input%EffecPopSize = 50                     
            if (input%InternalChrLenBasesMaCS) input%ChrLengthBasesMaCS=8.0E+08
            if (input%InternalMorgansMaCS) input%MorgansMaCS=1.43
            if (input%InternalMutationRateMaCS) input%MutationRateMaCS = 1.6 * (1.0 / real(input%ChrLengthBasesMaCS))
            input%AdditionalMaCSParams = "-eN   0.03   1 &
                        -eN   0.05   2 &
                        -eN   0.10   4 &
                        -eN   0.15   6 &
                        -eN   0.20   8 &
                        -eN   0.25  10 &
                        -eN   0.30  12 &
                        -eN   0.35  14 &
                        -eN   0.40  16 &
                        -eN   0.45  18 &
                        -eN   0.50  20 &
                        -eN   1.00  40 &
                        -eN   2.00  60 &
                        -eN   3.00  80 &
                        -eN   4.00 100 &
                        -eN   5.00 120 &
                        -eN  10.00 140 &
                        -eN  20.00 160 &
                        -eN  30.00 180 &
                        -eN  40.00 200 &
                        -eN  50.00 240 &
                        -eN 100.00 320 &
                        -eN 200.00 400 &
                        -eN 300.00 480 &
                        -eN 400.00 560 &
                        -eN 500.00 640"            
                    ! MacsParams="-t 3.52E-07    &
                    !             -r 3.52E-07    &"

        case default ! Read in param specified by user
            print *, "WARNING: User specified file with -eN options for MaCs are being read in."
            print *, "Please ensure that the following specfile parameters are set to External:"
            print *, "MutationRateMaCS, ChromosomeLengthBasesMaCS, EffectivePopulationSizeBaseMaCs"
            print *, "ChromosomeLengthMorgansMaCS"
            print*, " "
            ! EffecPopSize = 100
            open(unit=7001,file=trim(input%PopHistoryMaCS),status="old")
            status=0
            TmpC=""
            do while (status==0)
                read (7001, "(A)", advance='NO', size=size, iostat=status) buffer
                TmpC = TmpC // buffer
            enddo
            close(7001)
            input%AdditionalMaCSParams=trim(TmpC) ! A hack to be able to make MacsParams grow automagically

    end select

    call CalculateMaCSParams(input%MacsParams, input%ChrLengthBasesMaCS, input%EffecPopSize, input%MutationRateMaCS, input%MorgansMaCS, input%AdditionalMaCsParams)

end subroutine initMacsCases

subroutine CalculateMaCSParams(MacsParams, ChrLengthBasesMaCS, EffecPopSize, MutationRateMaCS, MorgansMaCS, AdditionalMaCsParams)
    ! SR written by SG to calculate parameters to pass to MaCS based on user defined values in specfile Oct 2015

    implicit none

    integer(kind=int64), intent(in) :: ChrLengthBasesMaCS, EffecPopSize
    double precision, intent(in) :: MutationRateMaCS, MorgansMaCS
    character(len=:),allocatable,intent(in) :: AdditionalMaCsParams
    character(len=:),allocatable,intent(inout) :: MacsParams

    double precision :: TOptionMaCS, ROptionMaCS, RecRateMaCS
    character(len=100) :: StrTOptionMaCS, StrROptionMaCS, StrChrLengthBasesMaCs


    RecRateMaCS = real(MorgansMaCS) / real(ChrLengthBasesMaCS)
    TOptionMaCS = 4 * EffecPopSize * MutationRateMaCS
    ROptionMaCS = 4 * EffecPopSize * RecRateMaCS

    write(StrTOptionMaCS, '(E10.2)') TOptionMaCS
    write(StrROptionMaCS, '(E10.2)') ROptionMaCS
    write(StrChrLengthBasesMaCs,'(i0)') ChrLengthBasesMaCS

    MacsParams = trim(StrChrLengthBasesMaCs) // " -t " // trim(StrTOptionMaCS) // " -r " // trim(StrROptionMaCS) // " " // trim(AdditionalMaCsParams)

end subroutine CalculateMaCSParams

subroutine ReadRecombSpecForMacs(Chrom)
  use AlphaHouseMod, only: countLines

    implicit none

    integer(kind=int32),intent(in) :: Chrom

    integer(kind=int64) :: i,nRecBlocks

    double precision :: RecombSpecPropSum
    double precision,allocatable,dimension(:) :: RecombSpecStart,RecombSpecStop,RecombSpecProp

    character(len=512) :: FileString

    write(FileString,'("./RecombinationSpecifications/Chromosome"i0,".txt")') Chrom
    nRecBlocks =  countLines(Filestring)

    allocate(RecombSpecStart(nRecBlocks))
    allocate(RecombSpecStop(nRecBlocks))
    allocate(RecombSpecProp(nRecBlocks))

    ! TODO unify this code and the same code in ReadRecombSpec
    open(unit=201,file=trim(FileString),status="old")
    RecombSpecPropSum=0.d0
    do i=1,nRecBlocks
        read(201,*) RecombSpecStart(i),RecombSpecStop(i),RecombSpecProp(i)
        RecombSpecPropSum=RecombSpecPropSum+RecombSpecProp(i)
        if (i > 1) then
            if ((RecombSpecStart(i) <= RecombSpecStart(i-1)).or.(RecombSpecStop(i) <= RecombSpecStop(i-1))) then
                print *, "ERROR: Variable recombination file should be sorted from start to end of chromosome"
                stop 110031
            endif
        endif
    enddo
    close(201)
    if (RecombSpecPropSum < 1.d0 .or. RecombSpecPropSum > 1.d0) then
        print*,"WARNING: proportion of recombinations over the defined regions must sum to 1.0 - rescale"
        print*, ""
        RecombSpecProp(:)=RecombSpecProp(:)/RecombSpecPropSum
    endif

    ! MaCS hotspots are defined with triplets of start,stop,ratio, where ratio is
    ! (the "standardized" genetic length of a region in cM) divided by
    ! (the physical length of a region in Mb).
    ! For example, a region is 60Mb long physically and 60cM long genetically and chromosome
    ! is R0=2M=200cM long genetically. Then the ratio is (60cM/2M)/60Mb=.5

    write(FileString,'("./Chromosomes/Chromosome"i0,"/Hotspots.txt")')Chrom
    open(unit=201,file=trim(FileString),status="unknown")
    do i=1,nRecBlocks
        !print*,RecombSpecStart(i),RecombSpecStop(i),(RecombSpecProp(i) / (RecombSpecStop(i) - RecombSpecStart(i)))
        write(201,'(3f12.8)') RecombSpecStart(i),RecombSpecStop(i),(RecombSpecProp(i) / (RecombSpecStop(i) - RecombSpecStart(i)))
    enddo
    flush(201)
    close(201)

    deallocate(RecombSpecStart)
    deallocate(RecombSpecStop)
    deallocate(RecombSpecProp)

end subroutine ReadRecombSpecForMacs
end Module macs
