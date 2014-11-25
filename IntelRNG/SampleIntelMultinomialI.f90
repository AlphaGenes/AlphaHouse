
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelMultinomialI(n,p)

  ! Sample n values from a Multinomial(p) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! p input (real), probabilities for the different categories

  ! Inspired by genmul() from RANLIB (John Burkardt).

  ! TODO: This is not from Intel so might consider some other name

  implicit none

  integer(kind=4),optional :: n
  integer(kind=4) :: nOpt,i,j,k
  integer(kind=4),dimension(:),allocatable :: SampleIntelMultinomialI
  integer(kind=4),dimension(1) :: b

  real(kind=8) :: pi,ptot
  real(kind=8),dimension(:) :: p

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  k=size(p)

  if (sum(p) > 0.99999e+00 ) then !TODO: is 0.99999e+00 enough?
    print*,"SampleIntelMultinomialI: sum of given probabilities > 1"
    stop
  end if

  allocate(SampleIntelMultinomialI(nOpt))

  do j=1,n
    ptot=1.0d0
    do i=1,k
      pi=p(i)/ptot
      b=SampleIntelBernoulliI(p=pi)
      if (b(1) > 0 .or. i == k) then
        SampleIntelMultinomialI(j)=i
        exit
      end if
      ptot=ptot-p(i)
    end do
  end do

  return

end function SampleIntelMultinomialI

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
