!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90 
!! - fittingpair.f90 
!!
      subroutine writeweightsunformatted(ndim,&
        maxnum_weights,num_weights,&
        maxnum_layers,num_layers,&
        nodes,weights,filename)
!!
      use fileunits
!!
      implicit none
!!
      integer ndim                    ! in
      integer maxnum_weights          ! in
      integer num_weights(ndim)      ! in
      integer maxnum_layers
      integer num_layers(ndim)
      integer i1,i2,i3,i4
      integer nodes(0:maxnum_layers,ndim)
      integer icount,jcount
!!
      real*8 weights(maxnum_weights,ndim)      
!!
      character*23 filename(ndim)
!!
!!
      do i1=1,ndim
        icount=0
        open(woutunit,file=filename(i1),form='unformatted',status='replace')
          do i2=1,num_layers(i1)    !'
!! connecting weights
            do i3=1,nodes(i2-1,i1)
              do i4=1,nodes(i2,i1)
                icount=icount+1
                write(woutunit)&
                  weights(icount,i1)
              enddo ! i4
            enddo ! i3
!! bias weights
            do i3=1,nodes(i2,i1)
              icount=icount+1
              write(woutunit)&
                weights(icount,i1)
            enddo ! i3
          enddo ! i2
        close(woutunit)
      enddo
!!
      return
      end
