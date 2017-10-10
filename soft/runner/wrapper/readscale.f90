!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Multipurpose subroutine

!! called by:
!! - predict.f90
!! - fitting.f90
!! - fittingpair.f90
!! - getpairsymfunctions.f90
!! - getsymmetryfunctions.f90
!!
      subroutine readscale(ndim,iswitch,&
           maxnum_funcvalues_local,num_funcvalues_local,&
           minvalue_local,maxvalue_local,avvalue_local,&
           eshortmin,eshortmax,chargemin,chargemax)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ndim                               ! in
      integer maxnum_funcvalues_local                  ! in
      integer num_funcvalues_local(ndim)               ! in
      integer i1,i2,i3                           ! internal
      integer iswitch                            ! in
!!
      real*8 avvalue_local(ndim,maxnum_funcvalues_local)     ! out 
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)    ! out 
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)    ! out 
      real*8 thres                               ! internal
      real*8 eshortmin                           ! out
      real*8 eshortmax                           ! out
      real*8 chargemin(nelem)                    ! out
      real*8 chargemax(nelem)                    ! out
!!
      character*15 filename                      ! internal
!!
      logical lexist                             ! internal
!!
      thres=0.00001d0
!!
      if(iswitch.eq.1)then
        filename='scaling.data'
      elseif(iswitch.eq.2)then
        filename='scaling.data'
      elseif(iswitch.eq.3)then
        filename='scalinge.data'
      else
        write(ounit,*)'ERROR: readscale called for wrong iswitch'
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
      do i1=1,ndim
        do i2=1,num_funcvalues_local(i1)
          read(scaleunit,*)i3,i3,minvalue_local(i1,i2),&
            maxvalue_local(i1,i3),avvalue_local(i1,i2)
        enddo ! i2
      enddo ! i1
      if(iswitch.eq.1)then
        read(scaleunit,*)eshortmin,eshortmax
      elseif(iswitch.eq.2)then
        read(scaleunit,*)eshortmin,eshortmax
      elseif(iswitch.eq.3)then
        do i2=1,nelem
          read(scaleunit,*)chargemin(i2),chargemax(i2)
        enddo ! i2
      endif
      close(scaleunit)
!!
      end
