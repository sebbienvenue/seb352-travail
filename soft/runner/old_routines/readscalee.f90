!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - predict.f90
!! - fitting.f90
!! - fitting_batch.f90
!! - fittingpair.f90
!! - getsymmetryfunctions.f90
!! - getpairsymmetryfunctions.f90
!!
      subroutine readscalee(iunit,ounit,nelem,nucelem,&
           maxnum_funcvaluese,num_funcvaluese,&
           minvaluee,maxvaluee,avvaluee,&
           chargemin,chargemax,&
           element,&
           lscalesym,ldebug)
!!
      implicit none
!!
      integer iunit                                            ! in
      integer ounit                                            ! in
      integer nelem                                            ! in
      integer maxnum_funcvaluese                               ! in
      integer num_funcvaluese(nelem)                           ! in
      integer i1,i2,i3                                         ! internal
      integer nucelem(nelem)                                   ! in
!!
      real*8 minvaluee(nelem,maxnum_funcvaluese)               ! out 
      real*8 maxvaluee(nelem,maxnum_funcvaluese)               ! out 
      real*8 avvaluee(nelem,maxnum_funcvaluese)                ! out 
      real*8 thres                                             ! internal
      real*8 chargemin(nelem)                                  ! out
      real*8 chargemax(nelem)                                  ! out
!!
      character*2 element(nelem)                               ! in
!!
      logical lscalesym                                        ! in
      logical ldebug                                           ! in
      logical lexist                                           ! internal
!!
      thres=0.00001d0
!!
      inquire(file='scalinge.data',exist=lexist)
      if(.not.lexist) then
        write(ounit,*)'Error: scalinge.data not found '
        stop
      endif
!!
      open(iunit,file='scalinge.data',form='formatted',status='old')
      rewind(iunit)
      do i2=1,nelem
        do i1=1,num_funcvaluese(i2)
          read(iunit,*)i3,i3,minvaluee(i2,i1),&
            maxvaluee(i2,i1),avvaluee(i2,i1)
        enddo ! i1
      enddo ! i2
      do i2=1,nelem
        read(iunit,*)chargemin(i2),chargemax(i2)
      enddo ! i2
      close(iunit)
!!
      do i1=1,nelem
        write(ounit,*)'============================================================='
        write(ounit,*)'Electrostatic symmetry function values for element ',element(i1)
        write(ounit,*)'Training set:  min           max       average         range '
!!        write(ounit,*)'-------------------------------------------------------------'
        do i2=1,num_funcvaluese(i1)
          write(ounit,'(i4,x,4f14.8)')i2,minvaluee(i1,i2),maxvaluee(i1,i2),&
            avvaluee(i1,i2),abs(maxvaluee(i1,i2)-minvaluee(i1,i2))
          if(abs(minvaluee(i1,i2)-maxvaluee(i1,i2)).lt.thres)then
            write(ounit,*)'### WARNING ###: minvaluee=maxvaluee ',i1,i2,nucelem(i1)
            if(lscalesym)then
              write(ounit,*)'scaling symmetry functions cannot be used with minvaluee=maxvaluee'
              stop
            endif
          endif
        enddo ! i2
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'
      write(ounit,*)' Minimum and Maximum charges from scalinge.data:'
      do i1=1,nelem
        write(ounit,'(x,a2,x,2f20.10)')element(i1),chargemin(i1),chargemax(i1)
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'

!!
      return
      end
