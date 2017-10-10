!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fitting_batch.f90
!!
      subroutine getfilenames(countepoch)
!!
      use fileunits
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer countepoch     ! in
      integer i1
!!
      character*20 filenametemp        ! internal
!!
!!
      if(countepoch.gt.999999)then
        write(ounit,*)'Number of epochs too large ',countepoch
        stop
      endif
!!
!! short range weights
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          filenametemp=filenamews(i1)
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
          filenamews(i1)=filenametemp
        enddo
      endif
!! charge weights
      if(lelec.and.(nn_type_elec.eq.1))then
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
      endif
!!
      return
      end
