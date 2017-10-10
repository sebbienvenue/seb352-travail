!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
      subroutine getfilenamespair(countepoch)
!!
      use fileunits
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer countepoch     ! in
      integer i1
!!
      character*20 filenametemp        ! internal
      character*23 filenametempp        ! internal
!!
!!

      if(countepoch.gt.999999)then
        write(ounit,*)'Number of epochs too large ',countepoch
        stop
      endif
!!
      do i1=1,npairs 
!! short range weights
        filenametempp=filenamewp(i1)
        if(countepoch.gt.99999)then
          write(filenametempp(1:6),'(i6)')countepoch
        elseif(countepoch.gt.9999)then
          write(filenametempp(2:6),'(i5)')countepoch
        elseif(countepoch.gt.999)then
          write(filenametempp(3:6),'(i4)')countepoch
        elseif(countepoch.gt.99)then
          write(filenametempp(4:6),'(i3)')countepoch
        elseif(countepoch.gt.9)then
          write(filenametempp(5:6),'(i2)')countepoch
        else
          write(filenametempp(6:6),'(i1)')countepoch
        endif
        filenamewp(i1)=filenametempp
      enddo

!! charge weights
      do i1=1,nelem
        filenametemp=filenamewe(i1)
        if(countepoch.gt.99999)then
          write(filenametemp(1:6),'(i6)')countepoch
        elseif(countepoch.gt.9999)then
          write(filenametemp(2:6),'(i5)')countepoch
        elseif(countepoch.gt.999)then
          write(filenametemp(3:6),'(i4)')countepoch
        elseif(countepoch.gt.99)then
          write(filenametemp(4:6),'(i3)')countepoch
        elseif(countepoch.gt.9)then
          write(filenametemp(5:6),'(i2)')countepoch
        else
          write(filenametemp(6:6),'(i1)')countepoch
        endif
        filenamewe(i1)=filenametemp
      enddo
!!
      return
      end
