!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by: 
!!            - fitting.f90 
!!            - fitting_batch.f90 
!!            - fittingpair.f90 
!!            - predict.f90 
!!            - getsymmetryfunctions.f90
!!            - getpairsymmetryfunctions.f90
!!
      subroutine readweights(iswitch,ndim,&
                 maxnum_weights_local,num_weights_local,weights_local)
!!
      use fileunits
      use globaloptions
      use nnflags
!!
      implicit none
!!
      integer ndim                                ! in
      integer iswitch                             ! in
      integer icount                              ! internal
      integer maxnum_weights_local                ! in
      integer num_weights_local(ndim)             ! in
      integer i1,i2,i3,i4                         ! internal
!! 
      real*8 weights_local(maxnum_weights_local,ndim)     ! out
!!
      character*40 filename
!!
      logical lexist
!!
      if(iswitch.eq.0)then
        do i1=1,ndim
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
          if(lreadunformatted)then
            open(wunit,file=filename,form='unformatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
            do i2=1,num_weights_local(i1)
              read(wunit)weights_local(i2,i1)
            enddo
          else
            open(wunit,file=filename,form='formatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights_local to be backwards compatible
            do i2=1,num_weights_local(i1)
              read(wunit,*)weights_local(i2,i1)
            enddo
          endif
          close(wunit)
        enddo ! i2
!!
      elseif(iswitch.eq.1)then
        do i1=1,ndim
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
          if(lreadunformatted)then
            open(wunit,file=filename,form='unformatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
            do i2=1,num_weights_local(i1)
              read(wunit)weights_local(i2,i1)
            enddo
          else
            open(wunit,file=filename,form='formatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights_local to be backwards compatible
            do i2=1,num_weights_local(i1)
              read(wunit,*)weights_local(i2,i1)
            enddo
          endif
          close(wunit)
        enddo ! i1
!!
      elseif(iswitch.eq.2)then
        icount=0
        do i1=1,nelem   ! don't put ndim here!
          do i2=i1,nelem
            icount = icount + 1
            filename='weights.000.000.data'
            if(nucelem(i1).gt.99)then
              write(filename(9:11),'(i3)')nucelem(i1)
            elseif(nucelem(i1).gt.9)then
              write(filename(10:11),'(i2)')nucelem(i1)
            else
              write(filename(11:11),'(i1)')nucelem(i1)
            endif
            if(nucelem(i2).gt.99)then
              write(filename(13:15),'(i3)')nucelem(i2)
            elseif(nucelem(i2).gt.9)then
              write(filename(14:15),'(i2)')nucelem(i2)
            else
              write(filename(15:15),'(i1)')nucelem(i2)
            endif
            inquire(file=filename,exist=lexist)
            if(.not.lexist) then
              write(*,*)'Error: file not found ',filename
              stop
            endif
            if(lreadunformatted)then
              open(wunit,file=filename,form='unformatted',status='old')
              rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
              do i3=1,num_weights_local(icount)
                read(wunit)weights_local(i3,icount)
              enddo
            else
              open(wunit,file=filename,form='formatted',status='old')
              rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
              do i3=1,num_weights_local(icount)
                read(wunit,*)weights_local(i3,icount)
              enddo
            endif
            close(wunit)
          enddo ! i2
        enddo ! i1
!!
      elseif(iswitch.eq.3)then !! Hextoff
        if(mode.eq.2)then
          write(filename,'(A,I3.3,A,I3.3,A,I3.3,A)')&
            'weights_hextoff.',hextoff_training_triplet(1),'.',&
            hextoff_training_triplet(2),'.',hextoff_training_triplet(3),&
            '.data'
          if(lreadunformatted)then
            open(wunit,file=filename,form='unformatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
            do i3=1,num_weights_local(1)
              read(wunit)weights_local(i3,1)
            enddo
          else
            open(wunit,file=filename,form='formatted',status='old')
            rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
            do i3=1,num_weights_local(1)
              read(wunit,*)weights_local(i3,1)
            enddo
          endif
          close(wunit)
        elseif(mode.eq.3)then
          do i1=1,nelem
            do i2=1,nelem
              if(nucelem(i2).ge.nucelem(i1))then
                do i3=1,nelem
                  icount = tripletindex(nucelem(i1),nucelem(i2),nucelem(i3))
                  write(filename,'(A,I3.3,A,I3.3,A,I3.3,A)')&
                  'weights_hextoff.',nucelem(i1),'.',nucelem(i2),'.',&
                  nucelem(i3),'.data'
                  if(lreadunformatted)then
                    open(wunit,file=filename,form='unformatted',status='old')
                    rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
                    do i4=1,num_weights_local(icount)
                      read(wunit)weights_local(i4,icount)
                    enddo
                  else
                    open(wunit,file=filename,form='formatted',status='old')
                    rewind(wunit)
!! we keep here an inconsistent order of the weights to be backwards compatible
                    do i4=1,num_weights_local(icount)
                      read(wunit,*)weights_local(i4,icount)
                    enddo
                  endif
                  close(wunit)
                enddo
              endif
            enddo
          enddo
        endif
!!
      else
        write(ounit,*)'ERROR: unknown iswitch in readweights ',iswitch
        stop
      endif ! iswitch 
!!
!!
      return
      end
