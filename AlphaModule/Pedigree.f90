module Pedigree

implicit none

integer :: GnAniPedI
integer,allocatable,dimension(:,:) :: GPedI

character(len=2000) :: GPedFileNameC
character(len=200),allocatable,dimension(:,:) :: GPedOrigC

end module Pedigree