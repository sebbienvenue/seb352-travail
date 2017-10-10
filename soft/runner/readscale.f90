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
            maxvalue_local(i1,i2),avvalue_local(i1,i2)
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
      do i1=1,ndim
        write(ounit,*)'============================================================='
        if(iswitch.eq.1)then
          write(ounit,*)'Short range symmetry function values for element ',element(i1)
        elseif(iswitch.eq.2)then
          write(ounit,*)'Pair symmetry function values for element pair ',elempair(i1,1),elempair(i1,2)
        elseif(iswitch.eq.3)then
          write(ounit,*)'Electrostatic symmetry function values for element ',element(i1)
        endif
        write(ounit,*)'Training set:  min           max       average         range '
!!        write(ounit,*)'-------------------------------------------------------------'
        do i3=1,num_funcvalues_local(i1)
!! check if default min and max values are still there => element or pair is not present, then don't write it
          if(minvalue_local(i1,i3).gt.maxvalue_local(i1,i3))then
            write(ounit,*)'No pairs of this type have been present in training set'
          else
            write(ounit,'(i4,x,4f14.8)')i3,minvalue_local(i1,i3),maxvalue_local(i1,i3),&
              avvalue_local(i1,i3),abs(maxvalue_local(i1,i3)-minvalue_local(i1,i3))
            if(abs(minvalue_local(i1,i3)-maxvalue_local(i1,i3)).lt.thres)then
              if(iswitch.eq.1)then
                write(ounit,*)'### WARNING ###: minvalue=maxvalue ',i1,i3,nucelem(i1)
              elseif(iswitch.eq.2)then
                write(ounit,*)'### WARNING ###: minvalue_short_pair=maxvalue_short_pair ',i1,i3,elempair(i1,1),elempair(i1,2)
              elseif(iswitch.eq.3)then
                write(ounit,*)'### WARNING ###: minvalue_elec=maxvalue_elec ',i1,i3,nucelem(i1)
              endif
              if(lscalesym)then
                if(iswitch.eq.1)then
                  write(ounit,*)'scaling symmetry functions cannot be used with minvalue=maxvalue'
                elseif(iswitch.eq.2)then
                  write(ounit,*)'scaling symmetry functions cannot be used with minvalue_short_pair=maxvalue_short_pair'
                elseif(iswitch.eq.3)then
                  write(ounit,*)'scaling symmetry functions cannot be used with minvalue_elec=maxvalue_elec'
                endif
                stop
              endif
            endif
          endif
        enddo ! i3
      enddo ! i1
      write(ounit,*)'-------------------------------------------------------------'
      if(iswitch.eq.1)then
        write(ounit,'(a,f20.10)')' eshortmin from scaling.data: ',eshortmin
        write(ounit,'(a,f20.10)')' eshortmax from scaling.data: ',eshortmax
      elseif(iswitch.eq.2)then
        write(ounit,'(a,f20.10)')' eshortmin from scaling.data: ',eshortmin
        write(ounit,'(a,f20.10)')' eshortmax from scaling.data: ',eshortmax
      elseif(iswitch.eq.3)then
        write(ounit,*)' Minimum and Maximum charges from scalinge.data:'
        do i1=1,nelem
          write(ounit,'(x,a2,x,2f20.10)')element(i1),chargemin(i1),chargemax(i1)
        enddo ! i1
        write(ounit,*)'-------------------------------------------------------------'
      endif !'
!!
      return
      end
