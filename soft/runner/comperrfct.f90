      function comperrfct(x)
!!
      implicit none
!!
      real*8 :: x,t,z
      real*8 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      real*8 sqrtpi
      real*8 comperrfct
!!
      sqrtpi=1.772453851d0
!!
      a1=-1.26551223d0
      a2=1.00002368d0
      a3=0.37409196d0
      a4=0.09678418d0
      a5=-0.18628806d0
      a6=0.27886807d0
      a7=-1.13520398d0
      a8=1.48851587d0
      a9=-0.82215223d0
      a10=0.17087277d0
!!
      z=abs(x)
      t=1.d0/(1.d0+0.5d0*z)
      comperrfct=t*exp(-z*z+a1 +t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*&
             (a8+t*(a9+t*a10)))))))))
      if(x.lt.0.d0)comperrfct=2.d0-comperrfct
!!
      end
