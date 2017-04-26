module newGenotypeMod
  
  type:: newGenoType
    
  end type

  type:: PhaseType
    private
    integer(kind=1), dimension(:), pointer, public:: phase
    integer::length
    contains


  end type

  interface PhaseType
    module procedure newPhaseType
  end interface

contains
  function newPhaseType(intIn) result(p)
    integer(kind=1), dimension(:), intent(in), target:: intIn
    type(PhaseType):: p

    p%length = size(intIn)
    allocate(p%phase(p%length))

  end function



end module newGenotypeMod

program testNewGenoType
  use newGenotypeMod
  implicit none

  type(newGenotype):: g
  type(PhaseType):: p
  integer, parameter:: N=4
  integer(kind=1), dimension(:), allocatable, target:: temp
  integer:: i
  integer(kind=1), dimension(:), allocatable, target:: temp2
  integer(kind=1):: missing
  integer(kind=1), dimension(:), pointer:: test

  allocate(temp(N))
  allocate(temp2(N))

  temp = 1
  temp2 = 0


  test(1:2) => temp
  test(3:N)=> temp2
  
  p = PhaseType(temp)

  deallocate(temp)
  write(*,*) p%phase

!  forall(i=1:size(p%phase), p%phase(i)-missing==0) 
!    temp2(i) = 1
!  end forall
  write(*,*) temp2
  write(*,*) test


end program testNewGenoType


