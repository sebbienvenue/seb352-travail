!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_short.f90 (1x)
!! - optimize_ewald.f90 (3x)
!! - fitforcesshort.f90
!!
      subroutine updatesteepest(num_weights,weights,dedw,error,step)
!!
      use fileunits
!!    DON'T use fittingoptions here, because step corresponds to different variables
!!
      implicit none
!!
      integer num_weights
      integer i1
!!
      real*8 weights(num_weights)
      real*8 dedw(num_weights)
      real*8 error
      real*8 step
!!
!!
      do i1=1,num_weights
        weights(i1)=weights(i1)+step*error*dedw(i1)
!!        weights(i1)=weights(i1)+0.01d0*error*dedw(i1)
      enddo ! i1
!!
      return
      end
