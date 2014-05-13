subroutine ReadPedigree

implicit none


call CountLines(GPedFileNameC,GnAniPedI)
print*, GnAniPedI

allocate(GPedOrigC(GnAniPedI,3))

end subroutine ReadPedigree
