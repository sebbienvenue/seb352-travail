!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
      subroutine nguyenwidrowweights(ndim,maxnum_layers_local,&
        nodes_local,seed_local,windex_local,&
        maxnodes_local,maxnum_weights_local,num_weights_local,&
        weights_local,weights_min_local,weights_max_local,actfunc_local)
!!
      use fileunits
!! Don't use fittingoptions here, because nguyenwidrowweights is called for multiple purposes
!!
      implicit none
!!
      integer ndim                                   ! in
      integer maxnum_layers_local                      ! in
      integer maxnum_weights_local                     ! in
      integer nodes_local(0:maxnum_layers_local,ndim) ! in
      integer num_weights_local(ndim)                 ! in
      integer seed_local                                   ! in
      integer maxnodes_local                          ! in
      integer i1,i2,i3                                ! internal
      integer windex_local(2*maxnum_layers_local,ndim)      ! in
!!
      real*8 weights_min_local                              ! in
      real*8 weights_max_local                              ! in
      real*8 weights_local(maxnum_weights_local,ndim) ! out
      real*8 w_norm                                   ! internal
      real*8 beta(ndim)                              ! internal
      real*8 ran0                                     ! internal
!!
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local,ndim) ! in
!!
!! compute beta for 1. hidden layer
      do i1=1,ndim
        beta(i1)=0.7d0*nodes_local(1,i1)**(1.0d0/real(nodes_local(0,i1)))
      enddo ! i1
!!
!! scale weights in 1. hidden layer
      do i1=1,ndim
        do i2=1,nodes_local(1,i1)
          w_norm=0.0d0
	  !! compute norm for each node in 1. hidden layer
          do i3=1,nodes_local(0,i1)
            w_norm=w_norm+(weights_local(i2+(i3-1)*nodes_local(1,i1),i1))**2
          enddo ! i3
          w_norm=sqrt(w_norm)
          !! scale weights
          do i3=1,nodes_local(0,i1)
            weights_local(i2+(i3-1)*nodes_local(1,i1),i1) &
            =weights_local(i2+(i3-1)*nodes_local(1,i1),i1)*beta(i1)/w_norm
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
!! scale bias weights in 1. hidden layer
      do i1=1,ndim
        do i2=1,nodes_local(1,i1)
          weights_local(i2+(nodes_local(0,i1)*nodes_local(1,i1)),i1) &
          =weights_local(i2+(nodes_local(0,i1)*nodes_local(1,i1)),i1)*beta(i1)
        enddo ! i2
      enddo ! i1
!!
!! if NN has two or more hl
      if(maxnum_layers_local.gt.2)then
   !! compute beta for 2. hidden layer
        do i1=1,ndim
          beta(i1)=0.7d0*nodes_local(2,i1)**(1.0d0/real(nodes_local(1,i1)))
        enddo ! i1
!!
   !! scale weights in 2. hidden layer
       do i1=1,ndim
         do i2=1,nodes_local(2,i1)
           w_norm=0.0d0
	  !! compute norm for each node in 2. hidden layer
           do i3=1,nodes_local(1,i1)
             w_norm=w_norm+(weights_local(windex_local(3,i1)+(i2-1)+(i3-1)*nodes_local(2,i1),i1))**2
           enddo ! i3
           w_norm=sqrt(w_norm)
           !! scale weights
           do i3=1,nodes_local(1,i1)
             weights_local(windex_local(3,i1)+(i2-1)+(i3-1)*nodes_local(2,i1),i1) &
             =weights_local(windex_local(3,i1)+(i2-1)+(i3-1)*nodes_local(2,i1),i1)*beta(i1)/w_norm
           enddo ! i3
         enddo ! i2
       enddo ! i1
!!
   !! scale bias weights in 2. hidden layer
       do i1=1,ndim
         do i2=1,nodes_local(2,i1)
           weights_local(windex_local(4,i1)+(i2-1),i1) &
           =weights_local(windex_local(4,i1)+(i2-1),i1)*beta(i1)
         enddo ! i2
       enddo ! i1
      endif
!!
!! scale weights to output layer
      weights_min_local=-0.5d0
      weights_max_local=0.5d0
      do i1=1,ndim
        do i2=windex_local(2*maxnum_layers_local-1,i1),num_weights_local(i1)
          weights_local(i2,i1)=weights_min_local+(weights_max_local-weights_min_local)*ran0(seed_local) 
        enddo ! i1
      enddo ! i2
!!
!! linearly rescale all weights and bias weights to allow different ranges of activation functions
!! This rescaling should be done only if the activation function in the first hidden layer is tanh
!! FIXME: At the moment we just check the activation function on the first node of the first HL
      do i1=1,ndim
        if(actfunc_local(1,1,i1).eq.'t')then
          do i2=1,num_weights_local(i1)
            weights_local(i2,i1)=weights_local(i2,i1)*2.0d0
          enddo ! i2
        endif
      enddo ! i1
!!
      return
      end

