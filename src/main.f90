program test

  use iso_fortran_env
  use individualModule
  use pedigreeModule
  use individualLinkedListModule
  implicit none

  type(individual) :: temp
  type(individualLinkedList) :: ll
  type(pedigreeHolder) :: pedigree
 integer, allocatable :: x(:)

	type(individualPointerContainer) :: t2  
  temp = initIndividual("1","2","3", 1, 0,1,500)
  ! temp%sirePointer => null()
  ! temp%damPointer => null()

  x = transfer(temp, x )
  temp = transfer(x, temp)
  ! print *, sizeof(pedigree)
  ! print *, sizeOf(ll)
  ! print *, sizeOf(t2)
  !   print *, sizeOf(temp%sirePointer)
  ! print *, sizeOf(temp%sirePointer)

  !  		print *, sizeof(temp%originalID)
  !  		print *, sizeof(temp%sireID)
  !  		print *, sizeof(temp%damID)
  !       print *, sizeof(temp%generation)
  !       print *, sizeof(temp%id)
  !       print *, sizeof(temp%originalPosition)
  !       print *, sizeof(temp%gender)
  !       print *, sizeof(temp%sirePointer)
  !       print *, sizeof(temp%damPointer)
  !       print *, sizeof(temp%OffSprings) !holds array of given size
  !       print *, sizeof(temp%nOffs) !number of offspring
  !       print*, sizeof(temp%Founder)
  !       print*, sizeof(temp%Genotyped)
  !       print*, sizeof(temp%HD)
  !       print*, sizeof(temp%isDummy)  ! if this animal is not in the pedigree, this will be true
  !       print*, sizeof(temp%isUnknownDummy)
  !       print *, sizeof(temp%individualGenotype)
  !       print *, sizeof(temp%individualPhase)
        ! print *, sizeof(temp%referAllele)
       	! print *, sizeof(temp%AlterAllele)
   print *, x 

  end program test
