!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine getstddev(&
        nelem
        ldebug)
!!
      use fileunits
!!
      implicit none
!!
      integer nelem
      
      logical ldebug



      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE avevar(data,n,ave,var)
      INTEGER n
      REAL*8 ave,var,data(n)
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      stddev=dsqrt(var)
      return
      END
