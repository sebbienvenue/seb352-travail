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
      subroutine getweightfilenames()
!!
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      integer icount
      integer i1,i2
!!
      character*20 filenametemp                              ! internal
      character*23 filenametempp
      character*30 filenametemppp
!!
      if(lshort.and.(nn_type_short.eq.1))then
        filenamews(:)          ='000000.short.000.out'
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        filenamewp(:)          ='000000.pair.000.000.out'
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        filenamewe(:)          ='000000.ewald.000.out'
      endif
      if(lnntb)then
        if(nntb_flag(3)) then
          write(filenametemppp,'(A15,I3.3,A,I3.3,A,I3.3,A4)')&
          '000000.hextoff.',hextoff_training_triplet(1),'.',&
          hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.out'
          filenamewhextoff = filenametemppp
        endif
 
      endif
      do i1=1,nelem
!! short range weights atomic NNs
        if(lshort.and.(nn_type_short.eq.1))then
          filenametemp=filenamews(i1)
          if(nucelem(i1).gt.99)then
            write(filenametemp(14:16),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(15:16),'(i2)')nucelem(i1)
          else
            write(filenametemp(16:16),'(i1)')nucelem(i1)
          endif
          filenamews(i1)=filenametemp
        endif
!! electrostatic weights
        if(lelec.and.(nn_type_elec.eq.1))then
          filenametemp=filenamewe(i1)
          if(nucelem(i1).gt.99)then
            write(filenametemp(14:16),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filenametemp(15:16),'(i2)')nucelem(i1)
          else
            write(filenametemp(16:16),'(i1)')nucelem(i1)
          endif
          filenamewe(i1)=filenametemp
        endif
      enddo
!!
!! short range weights pair NN
      if(lshort.and.(nn_type_short.eq.2))then
        icount=0
        do i1=1,nelem
          do i2=i1,nelem
             icount= icount +1
             filenametempp=filenamewp(icount)
             if(nucelem(i1).gt.99)then
               write(filenametempp(13:15),'(i3)')nucelem(i1)
             elseif(nucelem(i1).gt.9)then
               write(filenametempp(14:15),'(i2)')nucelem(i1)
             else
               write(filenametempp(15:15),'(i1)')nucelem(i1)
             endif
             if(nucelem(i2).gt.99)then
               write(filenametempp(17:19),'(i3)')nucelem(i2)
             elseif(nucelem(i2).gt.9)then
               write(filenametempp(18:19),'(i2)')nucelem(i2)
             else
               write(filenametempp(19:19),'(i1)')nucelem(i2)
             endif
             filenamewp(icount)=filenametempp
          enddo
        enddo
      endif
!!
      return
      end
