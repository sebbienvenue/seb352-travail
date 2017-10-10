!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine writescalee(iunit,ounit,nelem,&
           maxnum_funcvalues,num_funcvalues,&
           minvalue,maxvalue,avvalue,&
           chargemin,chargemax,&
           ldebug)
!!
      implicit none
!!
      integer iunit                ! in
      integer ounit                ! in
      integer nelem                ! in
      integer maxnum_funcvalues       ! in
      integer num_funcvalues(nelem)       ! in
      integer i1,i2,i3             ! internal
!!
      real*8 minvalue(nelem,maxnum_funcvalues)                 ! in 
      real*8 maxvalue(nelem,maxnum_funcvalues)                 ! in 
      real*8 avvalue(nelem,maxnum_funcvalues)                  ! in 
      real*8 chargemin(nelem)                                  ! in
      real*8 chargemax(nelem)                                  ! in
!!
      logical ldebug               ! in
!!
      open(iunit,file='scalinge.data',form='formatted',status='replace')
        do i2=1,nelem
          do i1=1,num_funcvalues(i2)
          write(iunit,'(i4,x,i4,x,3f18.9)')i2,i1,minvalue(i2,i1),&
            maxvalue(i2,i1),avvalue(i2,i1)
          enddo ! i1
        enddo ! i2
        do i2=1,nelem
          write(iunit,'(2f20.10)')chargemin(i2),chargemax(i2)
        enddo
      close(iunit)
!!
      return
      end
