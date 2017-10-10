!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!! - fitting_batch.f90
!!
      subroutine writescale(ndim,iswitch,&
           maxnum_funcvalues_local,num_funcvalues_local,&
           minvalue_local,maxvalue_local,avvalue_loal,&
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
      integer i1,i2                              ! internal
      integer iswitch                            ! in
!!
      real*8 minvalue_local(ndim,maxnum_funcvalues_local)                  ! in 
      real*8 maxvalue_local(ndim,maxnum_funcvalues_local)                  ! in 
      real*8 avvalue_loal(ndim,maxnum_funcvalues_local)                   ! in 
      real*8 eshortmin                                         ! in
      real*8 eshortmax                                         ! in
      real*8 chargemin(nelem)                                  ! in
      real*8 chargemax(nelem)                                  ! in
!!
!!
      if(iswitch.eq.1)then
        open(scaleunit,file='scaling.data',form='formatted',status='replace')
      elseif(iswitch.eq.2)then
        open(scaleunit,file='scaling.data',form='formatted',status='replace')
      elseif(iswitch.eq.3)then
        open(scaleunit,file='scalinge.data',form='formatted',status='replace')
      else
        write(*,*)'ERROR: wrong iswitch value in writescale ',iswitch
        stop
      endif !'

      do i2=1,ndim
        do i1=1,num_funcvalues_local(i2)
          write(scaleunit,'(i4,x,i4,x,3f18.9)')i2,i1,minvalue_local(i2,i1),&
            maxvalue_local(i2,i1),avvalue_loal(i2,i1)
          enddo ! i1
        enddo ! i2
      if(iswitch.eq.1)then
        write(scaleunit,'(2f20.10)')eshortmin,eshortmax
      elseif(iswitch.eq.2)then
        write(scaleunit,'(2f20.10)')eshortmin,eshortmax
      elseif(iswitch.eq.3)then
        do i2=1,nelem
          write(scaleunit,'(2f20.10)')chargemin(i2),chargemax(i2)
        enddo
      endif
!!      
      close(scaleunit)
!!
      return
      end
