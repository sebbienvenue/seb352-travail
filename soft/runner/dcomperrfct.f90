      function dcomperrfct(x)
      implicit none
      real*8 :: x
      real*8 sqrtpi
      real*8 dcomperrfct
!!
      sqrtpi=1.772453851d0
!!
      dcomperrfct=-2.d0*dexp(-x*x)/sqrtpi
!!
      end
