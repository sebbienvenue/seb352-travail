!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90 
!! - fitting_batch.f90 
!! - fittingpair.f90 
!!
      subroutine writeoptweights(iswitch,ndim,&
        maxnum_weights_local,num_weights_local,weights_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer iswitch
      integer icount
      integer ndim
      integer maxnum_weights_local
      integer num_weights_local(ndim)
      integer i1,i2
!!
      real*8 weights_local(maxnum_weights_local,ndim)      
!!
      character*22 filename(ndim)
      character*22 filenametemp
!!
!!
!! prepare filenames
      if(iswitch.eq.0)then
        filename(:)='optweights.000.out'
        do i1=1,ndim
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
      elseif(iswitch.eq.1)then
        filename(:)='optweightse.000.out'
        do i1=1,ndim
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
      elseif(iswitch.eq.2)then
        icount = 0
        filename(:)='optweights.000.000.out'
        do i1=1,nelem
          do i2=i1,nelem
           icount = icount +1
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
           elseif(nucelem(i2).gt.9)then
             write(filenametemp(17:18),'(i2)')nucelem(i2)
           else
             write(filenametemp(18:18),'(i1)')nucelem(i2)
           endif
           filename(icount)=filenametemp
          enddo
        enddo

      else
        write(ounit,*)'ERROR: unknown iswitch in writeoptweights ',iswitch
        stop
      endif
!!
      do i1=1,ndim
        if(lwriteunformatted)then
          open(symunit,file=filename(i1),form='unformatted',status='replace')
          do i2=1,num_weights_local(i1) !'
            write(symunit)weights_local(i2,i1)
          enddo
        else
          open(symunit,file=filename(i1),form='formatted',status='replace')
          do i2=1,num_weights_local(i1)
            write(symunit,'(f18.10)')weights_local(i2,i1)
          enddo
        endif !'
        close(symunit)
      enddo
!!
!!
      return
      end
