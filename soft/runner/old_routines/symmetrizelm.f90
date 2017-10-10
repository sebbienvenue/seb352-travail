!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine symmetrizelm(&
          dimlm,nelem,alphalm,&
          ldebug)      
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3,i4
      integer nelem
      integer dimlm

      real*8 alphalm(dimlm,dimlm,nelem)   ! in/out

      logical ldebug

      do i1=1,nelem
        do i2=2,dimlm
          do i3=1,i2-1
            alphalm(i3,i2,i1)=alphalm(i2,i3,i1)
          enddo ! i3
        enddo ! i2
      enddo ! i1

      return
      end
