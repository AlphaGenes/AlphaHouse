!@test
subroutine TestSampleIntelBernoulliI

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32),dimension(2) :: x

  call IntitialiseIntelRNG

  x=0
#line 15 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelBernoulliI(    p=0.0d0),x(1:1),"Samples from Bernoulli(p=0) should be 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 15) )
  if (anyExceptions()) return
#line 16 "tests/TestIntelRNG.pf"
#line 16 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelBernoulliI(n=2,p=0.0d0),x,     "Samples from Bernoulli(p=0) should be 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 16) )
  if (anyExceptions()) return
#line 17 "tests/TestIntelRNG.pf"

  x=1
#line 19 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelBernoulliI(    p=1.0d0),x(1:1),"Samples from Bernoulli(p=1) should be 1.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 19) )
  if (anyExceptions()) return
#line 20 "tests/TestIntelRNG.pf"
#line 20 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelBernoulliI(n=2,p=1.0d0),x,     "Samples from Bernoulli(p=1) should be 1.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 20) )
  if (anyExceptions()) return
#line 21 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGammaS

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real32),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGammaS(alpha=1.0,beta=1.0)
#line 45 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < 0.0),0,"Samples from Gamma should be positive", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 45) )
  if (anyExceptions()) return
#line 46 "tests/TestIntelRNG.pf"

  x=SampleIntelGammaS(n=n,alpha=1.0,beta=1.0)
#line 48 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < 0.0),0,"Samples from Gamma should be positive", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 48) )
  if (anyExceptions()) return
#line 49 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGammaD

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real64),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGammaD(alpha=1.0d0,beta=1.0d0)
#line 73 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < 0.0),0,"Samples from Gamma should be positive", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 73) )
  if (anyExceptions()) return
#line 74 "tests/TestIntelRNG.pf"

  x=SampleIntelGammaD(n=n,alpha=1.0d0,beta=1.0d0)
#line 76 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < 0.0),0,"Samples from Gamma should be positive", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 76) )
  if (anyExceptions()) return
#line 77 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGaussS

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real32) :: mean1,mean2,sd1,sd2
  real(real32),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGaussS()
#line 102 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < -10.0),0,"Samples from Gauss(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 102) )
  if (anyExceptions()) return
#line 103 "tests/TestIntelRNG.pf"
#line 103 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) > +10.0),0,"Samples from Gauss(0,1) should be < +10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 103) )
  if (anyExceptions()) return
#line 104 "tests/TestIntelRNG.pf"

  x=SampleIntelGaussS(n=n)
#line 106 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < -10.0),0,"Samples from Gauss(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 106) )
  if (anyExceptions()) return
#line 107 "tests/TestIntelRNG.pf"
#line 107 "tests/TestIntelRNG.pf"
  call assertEqual(count(x > +10.0),0,"Samples from Gauss(0,1) should be < +10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 107) )
  if (anyExceptions()) return
#line 108 "tests/TestIntelRNG.pf"
  mean1=sum(x)/real(n)
  x(:)=x(:)-mean1
  x(:)=x(:)*x(:)
  sd1=sqrt(sum(x)/real(n))

  x=SampleIntelGaussS(n=n,mu=10.0)
  mean2=sum(x)/real(n)
#line 115 "tests/TestIntelRNG.pf"
  call assertTrue(mean1 < mean2,"Mean of samples from Gauss(0,1) should be smaller than the mean of samples from Gauss(10,1)", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 115) )
  if (anyExceptions()) return
#line 116 "tests/TestIntelRNG.pf"

  x=SampleIntelGaussS(n=n,sigma2=10.0)
  mean2=sum(x)/real(n)
  x(:)=x(:)-mean2
  x(:)=x(:)*x(:)
  sd2=sqrt(sum(x)/real(n))
#line 122 "tests/TestIntelRNG.pf"
  call assertTrue(sd1 < sd2,"Standard deviation of samples from Gauss(0,1) should be smaller than the standard deviation of samples from Gauss(0,10)", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 122) )
  if (anyExceptions()) return
#line 123 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGaussD

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real64) :: mean1,mean2,sd1,sd2
  real(real64),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGaussD()
#line 148 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < -10.0),0,"Samples from Gauss(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 148) )
  if (anyExceptions()) return
#line 149 "tests/TestIntelRNG.pf"
#line 149 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) > +10.0),0,"Samples from Gauss(0,1) should be < +10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 149) )
  if (anyExceptions()) return
#line 150 "tests/TestIntelRNG.pf"

  x=SampleIntelGaussD(n=n)
#line 152 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < -10.0),0,"Samples from Gauss(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 152) )
  if (anyExceptions()) return
#line 153 "tests/TestIntelRNG.pf"
#line 153 "tests/TestIntelRNG.pf"
  call assertEqual(count(x > +10.0),0,"Samples from Gauss(0,1) should be < +10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 153) )
  if (anyExceptions()) return
#line 154 "tests/TestIntelRNG.pf"
  mean1=sum(x)/dble(n)
  x(:)=x(:)-mean1
  x(:)=x(:)*x(:)
  sd1=sqrt(sum(x)/dble(n))

  x=SampleIntelGaussD(n=n,mu=10.0d0)
  mean2=sum(x)/dble(n)
#line 161 "tests/TestIntelRNG.pf"
  call assertTrue(mean1 < mean2,"Mean of samples from Gauss(0,1) should be smaller than the mean of samples from Gauss(10,1)", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 161) )
  if (anyExceptions()) return
#line 162 "tests/TestIntelRNG.pf"

  x=SampleIntelGaussD(n=n,sigma2=10.0d0)
  mean2=sum(x)/dble(n)
  x(:)=x(:)-mean2
  x(:)=x(:)*x(:)
  sd2=sqrt(sum(x)/dble(n))
#line 168 "tests/TestIntelRNG.pf"
  call assertTrue(sd1 < sd2,"Standard deviation of samples from Gauss(0,1) should be smaller than the standard deviation of samples from Gauss(0,10)", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 168) )
  if (anyExceptions()) return
#line 169 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGumbelS

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real32),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGumbelS()
#line 193 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < -10.0),0,"Samples from Gumbel(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 193) )
  if (anyExceptions()) return
#line 194 "tests/TestIntelRNG.pf"

  x=SampleIntelGumbelS(n=n)
#line 196 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < -10.0),0,"Samples from Gumbel(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 196) )
  if (anyExceptions()) return
#line 197 "tests/TestIntelRNG.pf"
#line 197 "tests/TestIntelRNG.pf"
  call assertEqual(count(x > 100.0),0,"Samples from Gumbel(0,1) should be < 100.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 197) )
  if (anyExceptions()) return
#line 198 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelGumbelD

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real64),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelGumbelD()
#line 222 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) < -10.0),0,"Samples from Gumbel(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 222) )
  if (anyExceptions()) return
#line 223 "tests/TestIntelRNG.pf"

  x=SampleIntelGumbelD(n=n)
#line 225 "tests/TestIntelRNG.pf"
  call assertEqual(count(x < -10.0),0,"Samples from Gumbel(0,1) should be > -10.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 225) )
  if (anyExceptions()) return
#line 226 "tests/TestIntelRNG.pf"
#line 226 "tests/TestIntelRNG.pf"
  call assertEqual(count(x > 100.0),0,"Samples from Gumbel(0,1) should be < 100.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 226) )
  if (anyExceptions()) return
#line 227 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelMultinomialI

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n
  integer(int32),allocatable,dimension(:) :: x

  real(real64),allocatable,dimension(:) :: p,pObs

  call IntitialiseIntelRNG

  allocate(p(3))
  allocate(pObs(3))

  n=100
  allocate(x(n))

  p(1)=0.0d0
  p(2)=1.0d0
  p(3)=0.0d0
  x(1:2)=2
#line 258 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelMultinomialI(    p=p),x(1:1),"Samples from Multinomial(p=(0,1,0)) should be 2.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 258) )
  if (anyExceptions()) return
#line 259 "tests/TestIntelRNG.pf"
#line 259 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelMultinomialI(n=2,p=p),x(1:2),"Samples from Multinomial(p=(0,1,0)) should be 2.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 259) )
  if (anyExceptions()) return
#line 260 "tests/TestIntelRNG.pf"

  p(1)=0.0d0
  p(2)=0.9d0 ! works even is sum(p) not equal to 1, i.e., p's are taken as relative weights
  p(3)=0.0d0
  x(1:2)=2
#line 265 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelMultinomialI(    p=p),x(1:1),"Samples from Multinomial(p=(0,1,0)) should be 2.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 265) )
  if (anyExceptions()) return
#line 266 "tests/TestIntelRNG.pf"
#line 266 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelMultinomialI(n=2,p=p),x(1:2),"Samples from Multinomial(p=(0,1,0)) should be 2.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 266) )
  if (anyExceptions()) return
#line 267 "tests/TestIntelRNG.pf"

  p(1)=0.6d0
  p(2)=0.3d0
  p(3)=0.1d0
  x=SampleIntelMultinomialI(n=n,p=p)
  pObs(1)=count(x == 1)/float(n)
  pObs(2)=count(x == 2)/float(n)
  pObs(3)=count(x == 3)/float(n)
#line 275 "tests/TestIntelRNG.pf"
  call assertGreaterThan(pObs(1),pObs(2),"In Multinomial(p=(0.6,0.3,0.1)), Pr(class 1) > Pr(class 2).", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 275) )
  if (anyExceptions()) return
#line 276 "tests/TestIntelRNG.pf"
#line 276 "tests/TestIntelRNG.pf"
  call assertGreaterThan(pObs(2),pObs(3),"In Multinomial(p=(0.6,0.3,0.1)), Pr(class 2) > Pr(class 3).", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 276) )
  if (anyExceptions()) return
#line 277 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelPoissonI

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n
  integer(int32),allocatable,dimension(:) :: x

  real(real64) :: r

  call IntitialiseIntelRNG

  allocate(x(2))

  x(1:1)=SampleIntelPoissonI()
#line 301 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) >= 0),1,"Samples from Poisson(lambda=1.0) should be >= 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 301) )
  if (anyExceptions()) return
#line 302 "tests/TestIntelRNG.pf"

  x=SampleIntelPoissonI(n=2)
#line 304 "tests/TestIntelRNG.pf"
  call assertEqual(count(x >= 0),2,"Samples from Poisson(lambda=1.0) should be >= 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 304) )
  if (anyExceptions()) return
#line 305 "tests/TestIntelRNG.pf"

  deallocate(x)

  n=100
  allocate(x(n))

  x=SampleIntelPoissonI(n=n)
  r=sum(x)/float(n)
#line 313 "tests/TestIntelRNG.pf"
  call assertTrue(r < 2,"Mean of 100 samples from Poisson(lambda=1.0) should be < 2.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 313) )
  if (anyExceptions()) return
#line 314 "tests/TestIntelRNG.pf"

  x=SampleIntelPoissonI(n=n, lambda=10.0d0)
  r=sum(x)/float(n)
#line 317 "tests/TestIntelRNG.pf"
  call assertTrue(r >  5,"Mean of 100 samples from Poisson(lambda=10.0) should be >  5.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 317) )
  if (anyExceptions()) return
#line 318 "tests/TestIntelRNG.pf"
#line 318 "tests/TestIntelRNG.pf"
  call assertTrue(r < 15,"Mean of 100 samples from Poisson(lambda=10.0) should be < 15.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 318) )
  if (anyExceptions()) return
#line 319 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelUniformI

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n
  integer(int32),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:2)=0
#line 342 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelUniformI(    a=0,b=0),x(1:1),"Samples from Uniform(0,0) should be 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 342) )
  if (anyExceptions()) return
#line 343 "tests/TestIntelRNG.pf"
#line 343 "tests/TestIntelRNG.pf"
  call assertEqual(SampleIntelUniformI(n=2,a=0,b=0),x(1:2),"Samples from Uniform(0,0) should be 0.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 343) )
  if (anyExceptions()) return
#line 344 "tests/TestIntelRNG.pf"

  x=SampleIntelUniformI(n=n,a=5,b=10)
#line 346 "tests/TestIntelRNG.pf"
  call assertEqual(count(x <  5),0,"Samples from Uniform(5,10) should be all >=  5.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 346) )
  if (anyExceptions()) return
#line 347 "tests/TestIntelRNG.pf"
#line 347 "tests/TestIntelRNG.pf"
  call assertEqual(count(x > 10),0,"Samples from Uniform(5,10) should be all <= 10.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 347) )
  if (anyExceptions()) return
#line 348 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelUniformS

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real32),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelUniformS()
#line 372 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) <= 1.0),1,"Samples from Uniform(0,1) should be <= 1.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 372) )
  if (anyExceptions()) return
#line 373 "tests/TestIntelRNG.pf"

  x=SampleIntelUniformS(n=n)
#line 375 "tests/TestIntelRNG.pf"
  call assertEqual(count(x >= 0.0),n,"Samples from Uniform(0,1) should be >= 0.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 375) )
  if (anyExceptions()) return
#line 376 "tests/TestIntelRNG.pf"
#line 376 "tests/TestIntelRNG.pf"
  call assertEqual(count(x <= 1.0),n,"Samples from Uniform(0,1) should be <= 1.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 376) )
  if (anyExceptions()) return
#line 377 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

!@test
subroutine TestSampleIntelUniformD

  use pFUnit_mod
  use IntelRNGMod
  use ISO_Fortran_env

  implicit none

  integer(int32) :: n

  real(real64),allocatable,dimension(:) :: x

  call IntitialiseIntelRNG

  n=100
  allocate(x(n))

  x(1:1)=SampleIntelUniformD()
#line 401 "tests/TestIntelRNG.pf"
  call assertEqual(count(x(1:1) <= 1.0),1,"Samples from Uniform(0,1) should be <= 1.", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 401) )
  if (anyExceptions()) return
#line 402 "tests/TestIntelRNG.pf"

  x=SampleIntelUniformD(n=n)
#line 404 "tests/TestIntelRNG.pf"
  call assertEqual(count(x >= 0.0),n,"Samples from Uniform(0,1) should be >= 0.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 404) )
  if (anyExceptions()) return
#line 405 "tests/TestIntelRNG.pf"
#line 405 "tests/TestIntelRNG.pf"
  call assertEqual(count(x <= 1.0),n,"Samples from Uniform(0,1) should be <= 1.0", &
 & location=SourceLocation( &
 & 'TestIntelRNG.pf', &
 & 405) )
  if (anyExceptions()) return
#line 406 "tests/TestIntelRNG.pf"

  call UnintitialiseIntelRNG

end subroutine

module WrapTestIntelRNG
   use pFUnit_mod
   implicit none
   private

contains


end module WrapTestIntelRNG

function TestIntelRNG_suite() result(suite)
   use pFUnit_mod
   use WrapTestIntelRNG
   type (TestSuite) :: suite

   external TestSampleIntelBernoulliI
   external TestSampleIntelGammaS
   external TestSampleIntelGammaD
   external TestSampleIntelGaussS
   external TestSampleIntelGaussD
   external TestSampleIntelGumbelS
   external TestSampleIntelGumbelD
   external TestSampleIntelMultinomialI
   external TestSampleIntelPoissonI
   external TestSampleIntelUniformI
   external TestSampleIntelUniformS
   external TestSampleIntelUniformD


   suite = newTestSuite('TestIntelRNG_suite')

   call suite%addTest(newTestMethod('TestSampleIntelBernoulliI', TestSampleIntelBernoulliI))

   call suite%addTest(newTestMethod('TestSampleIntelGammaS', TestSampleIntelGammaS))

   call suite%addTest(newTestMethod('TestSampleIntelGammaD', TestSampleIntelGammaD))

   call suite%addTest(newTestMethod('TestSampleIntelGaussS', TestSampleIntelGaussS))

   call suite%addTest(newTestMethod('TestSampleIntelGaussD', TestSampleIntelGaussD))

   call suite%addTest(newTestMethod('TestSampleIntelGumbelS', TestSampleIntelGumbelS))

   call suite%addTest(newTestMethod('TestSampleIntelGumbelD', TestSampleIntelGumbelD))

   call suite%addTest(newTestMethod('TestSampleIntelMultinomialI', TestSampleIntelMultinomialI))

   call suite%addTest(newTestMethod('TestSampleIntelPoissonI', TestSampleIntelPoissonI))

   call suite%addTest(newTestMethod('TestSampleIntelUniformI', TestSampleIntelUniformI))

   call suite%addTest(newTestMethod('TestSampleIntelUniformS', TestSampleIntelUniformS))

   call suite%addTest(newTestMethod('TestSampleIntelUniformD', TestSampleIntelUniformD))


end function TestIntelRNG_suite

