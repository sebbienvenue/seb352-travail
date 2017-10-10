!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90 (4x)
!!
      subroutine writetmpweightsunformattedpair(nelem,nucelem,&
        maxnum_weights,num_weights,weights,&
        lswitch)
!!
      use fileunits
!!
      implicit none
!!
      integer nelem
      integer maxnum_weights
      integer num_weights(nelem)
      integer i1,i2
      integer nucelem(nelem)
!!
      real*8 weights(maxnum_weights,nelem)      
!!
      character*22 filename(nelem)
      character*22 filenametemp
!!
      logical lswitch
!!
      if(lswitch)then
       filename(:)='tmpweights.000.000.out'
       do i1=1,nelem
        do i2=i1,nelem
          filenametemp=filename(i1)
          if(nucelem(i1).gt.99)then
            write(filenametemp(12:14),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(13:14),'(i2)')nucelem(i1)
          else
            write(filenametemp(14:14),'(i1)')nucelem(i1)
          endif

          if(nucelem(i2).gt.99)then
            write(filenametemp(16:18),'(i3)')nucelem(i2)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(17:18),'(i2)')nucelem(i2)
          else
            write(filenametemp(18:18),'(i1)')nucelem(i2)
          endif

          filename(i1)=filenametemp
        enddo
       enddo
      endif
!!
      do i1=1,nelem
        open(woutunit,file=filename(i1),form='unformatted',status='replace')
        do i2=1,num_weights(i1)
          write(woutunit)weights(i2,i1)
        enddo
        close(woutunit)
      enddo
!!
!!
      return
      end
