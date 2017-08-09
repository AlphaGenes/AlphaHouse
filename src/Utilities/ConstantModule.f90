module ConstantModule
    use ISO_Fortran_Env

    integer, parameter :: ONEBYTEINT = 1
    integer, parameter :: MissingPhaseCode = 9
    integer, parameter :: ErrorPhaseCode = -1
    integer, parameter :: MissingGenotypeCode = 9
    integer, parameter :: MissingHaplotypeCode = -99
    character, parameter :: EMPTY_PARENT = '0'
    character, parameter :: EMPTYID = '0'
    integer, parameter :: IDLENGTH = 32
    integer, parameter :: IDINTLENGTH = 8
    integer, parameter :: generationThreshold = 1000
    integer, parameter :: OFFSPRINGTHRESHOLD = 20
    integer, parameter :: NOGENERATIONVALUE = -9999
    integer, parameter :: FILELENGTH = 256
    integer, parameter :: SPECOPTIONLENGTH = 300
    integer, parameter :: MISSINGGENDERCODE = -9
    integer , parameter ::  DEFAULTDICTSIZE = 4993
    integer, parameter :: DICT_NULL = -2147483647
    integer, parameter :: DICT_MULTIPLIER = 17
    integer, parameter :: MAXINT32 = 2147483647
    character(len=3), parameter:: DEFAULTNULLVAL = "nan"
    character(len=5), parameter:: DUMMYANIMALPREPRE = "DUMAL"
    character(len=1), parameter:: DEFAULTDELIM =","
!    character(len=1), parameter:: DEFAULTTYPE ='A'
!    character(len=1), dimension(3), parameter:: types = (/"A", "I", "R"/) !A for alphanumeric, I for integer, R for real64 and L for logical
    character(len=1), parameter:: DEFAULTCOMMENT="#"
    real(real64), parameter :: PI = 4.0d0 * atan(1.0d0)
    real(real64), parameter :: RAD2DEG = 180.0d0 / PI
    real(real64), parameter :: DEG2RAD = PI / 180.0d0
end Module ConstantModule
