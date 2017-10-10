!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by:
!!
      subroutine readscale_nntb(matrixsize,&
        iatom,jatom,katom,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_local,&
        value_min,value_max)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_funcvalues_local            ! in
      integer num_funcvalues_local               ! in
      integer i1,i2,i3                           ! internal
      integer matrixsize                         ! in
      integer i4                                 ! internal
      integer iatom,jatom,katom                  ! in
!!
      real*8 avvalue_local(maxnum_funcvalues_local)     ! out 
      real*8 maxvalue_local(maxnum_funcvalues_local)    ! out 
      real*8 minvalue_local(maxnum_funcvalues_local)    ! out 
      real*8 thres                               ! internal
      real*8 value1, value2                      ! internal
      real*8 value_min(matrixsize)                           ! out
      real*8 value_max(matrixsize)                           ! out
      character*50 dummy                               ! internal
!!
      character*100 filename                      ! internal
!!
      logical lexist                             ! internal
!!
      thres=0.00001d0
      value_min = 0.0d0
      value_max = 0.0d0
!!
      if(nntb_flag(3)) then
         write(filename,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'scaling_hextoff.',iatom,'.'&
           ,jatom,'.',katom,'.data'
      else
        write(ounit,*)'ERROR in readscale_nntb'
        stop
      endif
!!
      inquire(file=filename,exist=lexist)
      if(.not.lexist) then
        write(ounit,*)'Error: could not find ',filename
        stop
      endif
!!
      open(scaleunit,file=filename,form='formatted',status='old')
      rewind(scaleunit)
      do i2=1,num_funcvalues_local
        read(scaleunit,*)i3,minvalue_local(i2),&
            maxvalue_local(i2),avvalue_local(i2)
      enddo ! i2
      do i2=1,matrixsize
        read(scaleunit,'(i4,x,f18.9,f18.9))') i4, value1, value2
        !!value_min(i2), value_max(i2)
        value_min(i2) = value1
        value_max(i2) = value2
!!        write(*,*) value_min(i2), value_max(i2)
      enddo ! i2
      close(scaleunit)
!!
      write(ounit,*)'============================================================='
        write(ounit,*)'Hextoff symmetry function values for triplet ',&
          element(elementindex(hextoff_training_triplet(1))),&
          element(elementindex(hextoff_training_triplet(2))),&
          element(elementindex(hextoff_training_triplet(3)))
        write(ounit,*)'Training set:  min           max       average         range '
!!        write(ounit,*)'-------------------------------------------------------------'
        do i3=1,num_funcvalues_local
!! check if default min and max values are still there => element or pair is not present, then don't write it
          if(minvalue_local(i3).gt.maxvalue_local(i3))then
            write(ounit,*)'No triplets of this type have been present in training set'
          else
            write(ounit,'(i4,x,4f14.8)')i3,minvalue_local(i3),maxvalue_local(i3),&
              avvalue_local(i3),abs(maxvalue_local(i3)-minvalue_local(i3))
            if(abs(minvalue_local(i3)-maxvalue_local(i3)).lt.thres)then
               write(ounit,*)'### WARNING ###: minvalue_local=maxvalue_local ',i3
               stop
            endif
          endif
        enddo ! i3
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)'Hextoff hamiltonian min and max'
      write(ounit,*)'Training set: min                max'
      do i1=1,matrixsize
        write(ounit,'(i4,x,2f14.8)')i1,value_min(i1),value_max(i1)
      enddo
        write(ounit,*)'-------------------------------------------------------------'
!!
      return
      end
