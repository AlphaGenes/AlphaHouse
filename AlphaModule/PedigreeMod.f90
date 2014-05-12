
module PedigreeMod

implicit none

public

integer :: GnAniPedI
integer,allocatable,dimension(:,:) :: GPedI

character(len=2000) :: GPedFileNameC
character(len=200),allocatable,dimension(:,:) :: GPedOrigC

end module PedigreeMod
