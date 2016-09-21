
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

module PedigreeMod

  use ParameterFileMod
  use MiscellaneousMod

  implicit none

  integer :: GnAniPedI
  integer,allocatable,dimension(:,:) :: GPedI

  character(len=2000) :: GPedFileNameC
  character(len=200),allocatable,dimension(:,:) :: GPedOrigC

  contains

  include "ReadPedigree.f90"
  include "ReadParamPedigree.f90"

end module PedigreeMod

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
