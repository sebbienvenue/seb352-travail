!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: - fitting.f90 (2x)
!!            - prediction.f90 (2x)
!!            - getsymmetryfunctions.f90
!!            - getpairsymmetryfunctions.f90
!!
      subroutine readweightsunformatted(nelem,nucelem,&
                 maxnum_weights,num_weights,weights,&
                 lewald,lshort)
!!
      use fileunits
!!
      implicit none
!!
      integer nelem                         ! in
      integer maxnum_weights                ! in
      integer num_weights(nelem)            ! in
      integer i1,i2                         ! internal
      integer nucelem(nelem)                ! in
!! 
      real*8 weights(maxnum_weights,nelem)     ! out
!!
      character*20 filename
!!
      logical lewald
      logical lshort
      logical lexist
!!
      if(lshort)then
        do i1=1,nelem
          filename='weights.000.data'
          if(nucelem(i1).gt.99)then
            write(filename(9:11),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filename(10:11),'(i2)')nucelem(i1)
          else
            write(filename(11:11),'(i1)')nucelem(i1)
          endif
          inquire(file=filename,exist=lexist)
          if(.not.lexist) then
            write(*,*)'Error: file not found ',filename
            stop
          endif
          open(wunit,file=filename,form='unformatted',status='old')
          rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
          do i2=1,num_weights(i1)
            read(wunit)weights(i2,i1)
          enddo
          close(wunit)
        enddo ! i2
      endif ! lshort
!!
      if(lewald)then
        do i1=1,nelem
          filename='weightse.000.data'
          if(nucelem(i1).gt.99)then
            write(filename(10:12),'(i3)')nucelem(i1)
          elseif(nucelem(i1).gt.9)then
            write(filename(11:12),'(i2)')nucelem(i1)
          else
            write(filename(12:12),'(i1)')nucelem(i1)
          endif
          inquire(file=filename,exist=lexist)
          if(.not.lexist) then
            write(*,*)'Error: file not found ',filename
            stop
          endif
          open(wunit,file=filename,form='unformatted',status='old')
          rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
          do i2=1,num_weights(i1)
            read(wunit)weights(i2,i1)
          enddo
          close(wunit)
        enddo ! i2
      endif ! lewald
!!
!!
      return
      end
