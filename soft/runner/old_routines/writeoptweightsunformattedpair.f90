!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90 
!!
      subroutine writeoptweightsunformattedpair(nelem,npairs,nucelem,&
        maxnum_weights,num_weights,weights,&
        lswitch)
!!
      use fileunits
!!
      implicit none
!!
      integer npairs
      integer nelem
      integer maxnum_weights
      integer num_weights(npairs)
      integer i1,i2
      integer nucelem(nelem)
      integer paircount
!!
      real*8 weights(maxnum_weights,npairs)      
!!
      character*20 filename(nelem)
      character*22 filenamep(npairs)
      character*20 filenametemp
      character*22 filenametempp
!!
      logical lswitch
!!
      if(lswitch)then
        paircount = 0
        filenamep(:)='optweights.000.000.out'
        do i1=1,nelem
          do i2=i1,nelem
           paircount = paircount +1
           filenametempp=filenamep(i1)
           if(nucelem(i1).gt.99)then
             write(filenametempp(12:14),'(i3)')nucelem(i1)
           elseif(nucelem(i1).gt.9)then
             write(filenametempp(13:14),'(i2)')nucelem(i1)
           else
             write(filenametempp(14:14),'(i1)')nucelem(i1)
           endif
           if(nucelem(i2).gt.99)then
             write(filenametempp(16:18),'(i3)')nucelem(i2)
           elseif(nucelem(i2).gt.9)then
             write(filenametempp(17:18),'(i2)')nucelem(i2)
           else
             write(filenametempp(18:18),'(i1)')nucelem(i2)
           endif
           filenamep(paircount)=filenametempp
          enddo
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

      do i1=1,npairs 
        open(symunit,file=filenamep(i1),form='unformatted',status='replace')
        do i2=1,num_weights(i1)
          write(symunit,'(f18.10)')weights(i2,i1)
        enddo
        close(symunit)
      enddo
!!
!!
      return
      end
