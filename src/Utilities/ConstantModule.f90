module ConstantModule

    integer, parameter :: MissingPhaseCode = 9
    integer, parameter :: ErrorPhaseCode = -1
    integer, parameter :: MissingGenotypeCode = 9
    integer, parameter :: MissingHaplotypeCode = -99
    character, parameter :: EMPTY_PARENT = '0'
    character, parameter :: EMPTYID = '0'
    integer, parameter :: IDLENGTH = 32
    integer, parameter :: IDINTLENGTH = 8
    integer, parameter :: generationThreshold = 1000
    integer, parameter :: OFFSPRINGTHRESHOLD = 150
    integer, parameter :: NOGENERATIONVALUE = -9999
    integer, parameter :: FILELENGTH = 256
    integer, parameter :: SPECOPTIONLENGTH = 300
    integer, parameter :: MISSINGGENDERCODE = -9
    integer , parameter ::  DEFAULTDICTSIZE = 4993
    integer, parameter :: DICT_NULL = -2147483647
    integer, parameter :: DICT_MULTIPLIER = 17
    integer, parameter :: MAXINT32 = 2147483647
end Module ConstantModule
