!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by
!! fitting.f90
!! fitting_batch.f90
!! fittingpair.f90

!! Multipurpose subroutine

      subroutine getwshift(ndim,maxnum_weights_local,num_weights_local,&
        weights_local,weights_local_old,weights_local_veryold,&
        wshift_local,wshift_local2)
!!
      implicit none
!!
      integer ndim
      integer maxnum_weights_local
      integer num_weights_local(ndim)
      integer i1,i2
!!
      real*8 wshift_local(ndim)
      real*8 wshift_local2(ndim)
      real*8 weights_local(maxnum_weights_local,ndim)
      real*8 weights_local_old(maxnum_weights_local,ndim)
      real*8 weights_local_veryold(maxnum_weights_local,ndim)
!!
      wshift_local(:)=0.0d0
      wshift_local2(:)=0.0d0
      do i2=1,ndim
        do i1=1,num_weights_local(i2)
          wshift_local(i2)=wshift_local(i2)+abs(weights_local(i1,i2)-weights_local_old(i1,i2))
          wshift_local2(i2)=wshift_local2(i2)+abs(weights_local(i1,i2)-weights_local_veryold(i1,i2))
        enddo ! i1
        wshift_local(i2)=wshift_local(i2)/dble(num_weights_local(i2))
        wshift_local2(i2)=wshift_local2(i2)/dble(num_weights_local(i2))
      enddo ! i2
!!
      return
      end
