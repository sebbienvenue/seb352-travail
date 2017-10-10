!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine 

!! called by:
!! - fittingpair.f90
!!
      subroutine systematicweights(ndim,maxnum_layers_local,&
        num_layers_local,nodes_local,&
        maxnum_weights_local,&
        weights_local,weights_min_local,weights_max_local)
!!
      use fileunits
!! Don't use fittingoptions here, because nguyenwidrowweights is called for multiple purposes
!!
      implicit none
!!
      integer i1,i2,i3,i4                             ! internal
      integer ndim                                   ! in
      integer maxnum_layers_local                      ! in
      integer maxnum_weights_local                     ! in
      integer nodes_local(0:maxnum_layers_local,ndim) ! in
      integer icount                                  ! internal
      integer nintervals                              ! internal
      integer num_layers_local(ndim)                  ! in
!!
      real*8 weights_local(maxnum_weights_local,ndim) ! out 
      real*8 weights_min_local                              ! in
      real*8 weights_max_local                              ! out 
      real*8 wstep                                    ! internal
      real*8 wtemp                                    ! internal
!!
!! loop over all elements
      do i1=1,ndim
        icount=0 
!! assign the connecting weights
        do i2=1,num_layers_local(i1)
          nintervals=nodes_local(i2,i1)-1 ! 1 less than number of nodes 
!! loop over all nodes in previous layer
          do i3=1,nodes_local(i2-1,i1)
            wstep=(weights_max_local-weights_min_local)/dble(nintervals)
            wtemp=weights_min_local
!! loop over all nodes in target layer
            do i4=1,nodes_local(i2,i1)
              icount=icount+1
              weights_local(icount,i1)=wtemp
              wtemp=wtemp+wstep
            enddo ! i4
          enddo ! i3
!! bias weights
          wtemp=weights_min_local
          do i3=1,nodes_local(i2,i1)
            icount=icount+1
            weights_local(icount,i1)=wtemp
            wtemp=wtemp+wstep
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      return
      end
