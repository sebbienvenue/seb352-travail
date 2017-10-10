!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! this is a multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!!
      subroutine initialweights(ndim,&
      maxnum_weights_local,num_weights_local,num_weights_elec,&
      maxnum_layers_local,num_layers_local,windex_local,nodes_local,&
      jseed,nseed,weights_local,weights_elec)
!!
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer ndim                       ! in
      integer maxnum_weights_local       ! in
      integer num_weights_local(ndim)    ! in
      integer num_weights_elec(nelem)    ! in
      integer maxnum_layers_local        ! in
      integer num_layers_local(ndim)     ! in
      integer jseed                      ! in/out
      integer nseed                      ! in/out
      integer i,j,k                      ! internal
      integer windex_local(2*maxnum_layers_local,ndim) ! in
      integer nodes_local(0:maxnum_layers_local,ndim)  ! in
!!
      real*8 weights_local(maxnum_weights_local,ndim)
      real*8 weights_elec(maxnum_weights_elec,nelem)
      real*8 ran0
      real*8 ran1
      real*8 ran2
      real*8 ran3
!!
!!
      if(lshort)then
!! we keep here the somewhat inconsistent looping order to be backwards compatible
        if(nran.eq.1)then
          do j=1,ndim
            do i=1,num_weights_local(j)
              weights_local(i,j)=weights_min+(weights_max-weights_min)*ran0(jseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.2)then
          do j=1,ndim
            do i=1,num_weights_local(j)
              weights_local(i,j)=weights_min+(weights_max-weights_min)*ran1(jseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.3)then
          do j=1,ndim
            do i=1,num_weights_local(j)
              weights_local(i,j)=weights_min+(weights_max-weights_min)*ran2(jseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.4)then
          do j=1,ndim
            do i=1,num_weights_local(j)
              weights_local(i,j)=weights_min+(weights_max-weights_min)*ran3(jseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        endif ! nran
!! overwrite biasweights if requested
        if(lseparatebiasini)then
          do j=1,ndim
            do i=1,num_layers_local(j)
              do k=windex_local(2*i,j),windex_local(2*i,j)+nodes_local(i,j)-1
                weights_local(k,j)=biasweights_min+(biasweights_max-biasweights_min)*ran0(jseed)  
              enddo ! k
            enddo ! i
          enddo ! j
        endif ! lseparatebiasini
      endif ! lshort
!!
      if(lelec.and.(nn_type_elec.eq.1))then
!! we keep here the somewhat inconsistent looping order to be backwards compatible
        if(nran.eq.1)then
          do j=1,nelem
            do i=1,num_weights_elec(j)
              weights_elec(i,j)=weightse_min+(weightse_max-weightse_min)*ran0(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.2)then
          do j=1,nelem
            do i=1,num_weights_elec(j)
              weights_elec(i,j)=weightse_min+(weightse_max-weightse_min)*ran1(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.3)then
          do j=1,nelem
            do i=1,num_weights_elec(j)
              weights_elec(i,j)=weightse_min+(weightse_max-weightse_min)*ran2(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.4)then
          do j=1,nelem
            do i=1,num_weights_elec(j)
              weights_elec(i,j)=weightse_min+(weightse_max-weightse_min)*ran3(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        endif ! nran
      endif ! lelec
!!
!!
      return
      end
!!
