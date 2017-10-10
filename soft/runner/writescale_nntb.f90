!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
      subroutine writescale_nntb(matrixsize,&
        iatom,jatom,katom,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        minvalue_local,maxvalue_local,avvalue_loal,&
        value_min,value_max)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_funcvalues_local            ! in
      integer num_funcvalues_local               ! in
      integer i1,i2                              ! internal
      integer iatom                              ! in
      integer jatom                              ! in
      integer katom                              ! in
      integer matrixsize                         ! in
!!
      real*8 minvalue_local(maxnum_funcvalues_local)           ! in 
      real*8 maxvalue_local(maxnum_funcvalues_local)           ! in 
      real*8 avvalue_loal(maxnum_funcvalues_local)             ! in 
      real*8 value_min(matrixsize)                             ! in
      real*8 value_max(matrixsize)                             ! in
      character*100 filenametemp                               ! internal

!!
!!
      if(nntb_flag(3))then
        write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
        'scaling_hextoff.',iatom,'.'&
         ,jatom,'.',katom,'.data'
      else
        write(ounit,*)'ERROR in writescale_nntb'
        stop
      endif
      open(scaleunit,file=filenametemp,form='formatted',status='replace')
!!     
      do i1=1,num_funcvalues_local
        write(scaleunit,'(i4,x,3f18.9)')i1,minvalue_local(i1),&
          maxvalue_local(i1),avvalue_loal(i1)
      enddo ! i1
!!
      do i1=1,matrixsize
        write(scaleunit,'(i4,x,2f18.9)') i1,value_min(i1),value_max(i1)
      enddo
!!
      close(scaleunit)
!!
      return
      end
