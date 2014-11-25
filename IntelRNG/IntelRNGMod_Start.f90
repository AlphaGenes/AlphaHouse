
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

include 'mkl_vsl.f90'

module IntelRNGMod

  ! Random Number Generation (RNG) module
  !   "interface" to the Intel MKL Vector Statistical Library (VSL) RNG capabilities

  ! https://software.intel.com/en-us/node/470592 (2014-11-25)

  use mkl_vsl_type
  use mkl_vsl

  implicit none

  integer(kind=4) RNGErrCode,RNGMethod

  type(vsl_stream_state) :: RNGStream

  private

  ! RNG stream management
  public :: IntitialiseIntelRNG,UnintitialiseIntelRNG

  ! Discrete
  public :: SampleIntelUniformI
  public :: SampleIntelBernoulliI
  public :: SampleIntelMultinomialI

  ! Continuous
  public :: SampleIntelUniformRS,SampleIntelUniformRD
  public :: SampleIntelGumbelRS,SampleIntelGumbelRD

  contains

