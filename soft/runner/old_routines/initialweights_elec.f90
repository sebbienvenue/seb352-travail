!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine initialweights_elec(&
        num_weights_elec,nseed,weights_elec)
!!
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer num_weights_elec(nelem)    ! in
      integer nseed                      ! in/out
      integer i,j,k                      ! internal
!!
      real*8 weights_elec(maxnum_weights_elec,nelem)
      real*8 ran0
      real*8 ran1
      real*8 ran2
      real*8 ran3
!!
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
!!
!!
      return
      end
!!
