program test

  use iso_fortran_env
  use individualModule
  use pedigreeModule
  use compatibilityModule
  implicit none

  type(PedigreeHolder) :: ped


    call readPlink("merge_final", ped)

  end program test
