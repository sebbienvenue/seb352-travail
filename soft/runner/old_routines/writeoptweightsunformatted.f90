!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90 
!! - fittingpair.f90 
!!
      subroutine writeoptweightsunformatted(nelem,nucelem,&
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
      character*20 filename(nelem)
      character*20 filenametemp
!!
      logical lswitch
!!
      if(lswitch)then
        filename(:)='optweights.000.out'
        do i1=1,nelem
          filenametemp=filename(i1)
          if(nucelem(i1).gt.99)then
            write(filenametemp(12:14),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(13:14),'(i2)')nucelem(i1)
          else
            write(filenametemp(14:14),'(i1)')nucelem(i1)
          endif
          filename(i1)=filenametemp
        enddo
      else
        filename(:)='optweightse.000.out'
        do i1=1,nelem
          filenametemp=filename(i1)
          if(nucelem(i1).gt.99)then
            write(filenametemp(13:15),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(14:15),'(i2)')nucelem(i1)
          else
            write(filenametemp(15:15),'(i1)')nucelem(i1)
          endif
          filename(i1)=filenametemp
        enddo
      endif
!!
      do i1=1,nelem
        open(symunit,file=filename(i1),form='unformatted',status='replace')
        do i2=1,num_weights(i1)
          write(symunit)weights(i2,i1)
        enddo
        close(symunit)
      enddo
!!
!!
      return
      end
