!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine initialweightspair(&
      num_weightsewald,num_weightspair,iseed,nseed,weights_ewald,weights_pair,&
      weights_min,weights_max,weightse_min,weightse_max,&
      lewald,lshort)
!!
      use globaloptions

      implicit none
!!
      integer num_weightsewald(nelem)    ! in
      integer num_weightspair(npairs)
      integer iseed                      ! in/out
      integer nseed                      ! in/out
      integer i,j                        ! internal
!!
      real*8 weights_ewald(maxnum_weightsewald,nelem)
      real*8 weights_pair(maxnum_weightspair,npairs)
      real*8 ran0
      real*8 ran1
      real*8 ran2
      real*8 ran3
      real*8 weights_min                 ! in
      real*8 weights_max                 ! in
      real*8 weightse_min                 ! in
      real*8 weightse_max                 ! in
!!
      logical lewald
      logical lshort
!!
      if(lshort)then
!! we keep here the somewhat inconsistent looping order to be backwards compatible
        if(nran.eq.1)then
          do j=1,npairs 
            do i=1,num_weightspair(j)
              weights_pair(i,j)=weights_min+(weights_max-weights_min)*ran0(iseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.2)then
          do j=1,npairs
            do i=1,num_weightspair(j)
              weights_pair(i,j)=weights_min+(weights_max-weights_min)*ran1(iseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.3)then
          do j=1,npairs
            do i=1,num_weightspair(j)
              weights_pair(i,j)=weights_min+(weights_max-weights_min)*ran2(iseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        elseif(nran.eq.4)then
          do j=1,npairs
            do i=1,num_weightspair(j)
              weights_pair(i,j)=weights_min+(weights_max-weights_min)*ran3(iseed)  !-0.5d0*(weights_max-weights_min)
            enddo
          enddo
        endif ! nran
      endif ! lshort

!!
      if(lewald)then
!! we keep here the somewhat inconsistent looping order to be backwards compatible
        if(nran.eq.1)then
          do j=1,nelem
            do i=1,num_weightsewald(j)
              weights_ewald(i,j)=weightse_min+(weightse_max-weightse_min)*ran0(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.2)then
          do j=1,nelem
            do i=1,num_weightsewald(j)
              weights_ewald(i,j)=weightse_min+(weightse_max-weightse_min)*ran1(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.3)then
          do j=1,nelem
            do i=1,num_weightsewald(j)
              weights_ewald(i,j)=weightse_min+(weightse_max-weightse_min)*ran2(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        elseif(nran.eq.4)then
          do j=1,nelem
            do i=1,num_weightsewald(j)
              weights_ewald(i,j)=weightse_min+(weightse_max-weightse_min)*ran3(nseed)  !-0.5d0*(weightse_max-weightse_min)
            enddo
          enddo
        endif ! nran
      endif ! lewald
!!
!!
      return
      end
