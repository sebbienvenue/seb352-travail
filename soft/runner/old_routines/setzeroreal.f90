!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: set a number of real*8 numbers to 0.0d0
!! can also be used for non-arrays

!! called by:
!! fitting.f90
!! fitting_batch.f90
!! fittingpair.f90
!!
      subroutine setzeroreal(ndim,array)
!!  
      implicit none   
!!
      integer ndim
      integer i1
!!
      real*8 array(ndim)
!!
      do i1=1,ndim
        array(i1)=0.0d0
      enddo
!!
      return
      end
