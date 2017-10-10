!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine getweightfilenames_elec()
!!
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer icount
      integer i1,i2
!!
      character*20 filenametemp                              ! internal
!!
      filenamewe(:)          ='000000.ewald.000.out'
      do i1=1,nelem
!! electrostatic weights
        filenametemp=filenamewe(i1)
        if(nucelem(i1).gt.99)then
          write(filenametemp(14:16),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(filenametemp(15:16),'(i2)')nucelem(i1)
        else
          write(filenametemp(16:16),'(i1)')nucelem(i1)
        endif
        filenamewe(i1)=filenametemp
      enddo
!!
      return
      end
