
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function SampleIntelMultinomialI(n,p)

  ! Sample n values from a Multinomial(p) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! p input (real), probabilities for the different categories

  ! TODO: This is not from Intel so might consider some other name

  implicit none

  integer(kind=4),optional,intent(in) :: n
  integer(kind=4) :: nOpt,i,j,k
  integer(kind=4),dimension(:),allocatable :: SampleIntelMultinomialI
  integer(kind=4),dimension(1) :: b

  real(kind=8) :: pi,psum,psumtmp
  real(kind=8),dimension(:),intent(inout) :: p

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  k=size(p)

  psum=sum(p)
  if (abs(psum - 1.0d0) > 1e-7) then
    p=p/psum ! rescale
  end if
  psum=1.0d0

  allocate(SampleIntelMultinomialI(nOpt))

  ! Over samples
  do j=1,nOpt

    ! Over categories of a sample
    psumtmp=psum
    do i=1,k

      if (p(i) > 0.0d0) then ! sample only if needed (border case)

        if (p(i) < 1.0d0) then ! likewise
          pi=p(i)/psumtmp

          if (pi < 1.0d0) then ! likewise
            b=SampleIntelBernoulliI(p=pi)

            if (b(1) > 0) then
              SampleIntelMultinomialI(j)=i
              exit
            end if

            psumtmp=psumtmp-p(i)

          else
            SampleIntelMultinomialI(j)=i
            exit
          end if

        else
          SampleIntelMultinomialI(j)=i
          exit
        end if

      end if

    end do

  end do

  return

end function SampleIntelMultinomialI

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
