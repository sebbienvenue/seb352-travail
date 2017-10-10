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
      subroutine updatesteepest_para(nelem,n_start,n_end,&
                num_weights,weights,dedw,error,step)
!!
      use mpi_mod
      use fileunits
!!
      implicit none
!!
      integer nelem
      integer num_weights
      integer i1,i2
      integer nweights                 ! internal
      integer n_start,n_end            ! internal
!!
      real*8 weights(num_weights) ! in/out
      real*8 dedw(num_weights)    ! internal
      real*8 error                     ! in
      real*8 step                      ! in
!!
!!
      do i1=n_start,n_end
        weights(i1)=weights(i1)+step*error*dedw(i1)
      enddo ! i1
!!
      return
      end
