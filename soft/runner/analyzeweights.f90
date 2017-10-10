!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! Purpose: analyze weight vector and print information to output file

!! called by:
!!
      subroutine analyzeweights(iswitch,countepoch,ndim,&
        maxnum_weights_local,num_weights_local,weights_local)
!!
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer i1,i2
      integer iswitch                                    ! in
      integer countepoch                                 ! in
      integer ndim                                       ! in
      integer maxnum_weights_local                       ! in
      integer num_weights_local(ndim)                    ! in
!!
      real*8 weights_local(maxnum_weights_local,ndim)    ! in
      real*8 absweights(ndim)                            ! internal
      real*8 avweights(ndim)                             ! internal
      real*8 stddevweights(ndim)                         ! internal
!!
!! initializations
      absweights(:)   =0.0d0
      avweights(:)    =0.0d0
      stddevweights(:)=0.0d0
!!
!! determine the length of the weight vectors and averages
      do i1=1,ndim
        do i2=1,num_weights_local(i1)
          absweights(i1)=absweights(i1)+weights_local(i2,i1)**2.0d0
          avweights(i1) =avweights(i1) +weights_local(i2,i1)
        enddo
!! absweights is normalized here to avoid very large numbers for many weights
        absweights(i1)=dsqrt(absweights(i1))/dble(num_weights_local(i1))
        avweights(i1) =avweights(i1)/dble(num_weights_local(i1))
      enddo
!!
!! determine the standard deviation of the weights
      do i1=1,ndim
        do i2=1,num_weights_local(i1)
          stddevweights(i1)=stddevweights(i1)+(weights_local(i2,i1)-avweights(i1))**2.0d0
        enddo
        stddevweights(i1)=dsqrt(stddevweights(i1)/dble(num_weights_local(i1)))
      enddo

!! write output to runner.out
      if(iswitch.eq.1)then ! short range atomic weights
        do i1=1,ndim
          write(ounit,'(a10,x,a2,i5,x,3f14.6,a20)') ' WEIGHTS  ',element(i1),countepoch,&
            absweights(i1),avweights(i1),stddevweights(i1),' (abs,avrg,stddev)'
        enddo
      elseif(iswitch.eq.2)then ! electrostatic weights
        do i1=1,ndim
          write(ounit,'(a10,x,a2,i5,x,3f14.6,a20)') ' WEIGHTSE ',element(i1),countepoch,&
            absweights(i1),avweights(i1),stddevweights(i1),' (abs,avrg,stddev)'
        enddo
      elseif(iswitch.eq.3)then ! short range pair NN
        write(ounit,*)'### WARNING ### analyzeweights is not implemented for pair case' 
      endif ! iswitch
!!
      return
      end
