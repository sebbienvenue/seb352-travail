!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
      subroutine getscale_hextoff(matrixsize,&
        maxnum_funcvalues_local,num_funcvalues_local,&
        ntrain,minvalue_local,&
        maxvalue_local,avvalue_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer maxnum_funcvalues_local                     ! in
      integer num_funcvalues_local                        ! in
      integer matrixsize                                  ! in
      integer npoints                                     ! internal
      integer ncount                                      ! internal
      integer ntrain                                      ! in
      integer counter                                     ! internal
      integer i1,i2,i3,i4                                 ! internal
      integer correlation_stat(maxnum_funcvalues_local,10)                  ! internal
!!
      real*8 symfunction_list_local(maxnum_funcvalues_local,nblock)  ! internal
!!
      real*8 minvalue_local(maxnum_funcvalues_local)                 ! out
      real*8 maxvalue_local(maxnum_funcvalues_local)                 ! out
      real*8 avvalue_local(maxnum_funcvalues_local)                  ! out
      real*8 stddev_local(maxnum_funcvalues_local)                   ! internal
      real*8 hextoffelem_local(matrixsize,nblock)
      real*8 hextoffmin(matrixsize)
      real*8 hextoffmax(matrixsize)
      real*8 thres
      real*8 zdummy
      real*8 pearson_corr(maxnum_funcvalues_local,maxnum_funcvalues_local) ! internal
!!
!! initialization
      minvalue_local(:)    = 100000.d0
      maxvalue_local(:)    =-100000.d0
      avvalue_local(:)     = 0.d0
      stddev_local(:)      = 0.d0
      pearson_corr(:,:)    = 0.0d0
      hextoffmin           = 10000.d10
      hextoffmax           =-10000.d10
      counter              =0
      correlation_stat(:,:)=0
      thres                =0.0001d0
!!
      ncount=ntrain
 10   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! read a block of points
      do i1=1,npoints
        read(symhextoffunit,*) (symfunction_list_local(i3,i1),&
          i3=1,num_funcvalues_local)
        read(symhextoffunit,*) (hextoffelem_local(i2,i1),&
          i2=1,matrixsize)
      enddo ! i1
!!
!! modify minvalue_local and maxvalue_local
      do i1=1,npoints
!! count for later averaging
        counter=counter+1
        do i3=1,num_funcvalues_local
          minvalue_local(i3)&
            =min(minvalue_local(i3),symfunction_list_local(i3,i1))
          maxvalue_local(i3)&
            =max(maxvalue_local(i3),symfunction_list_local(i3,i1))
          avvalue_local(i3) &
            =avvalue_local(i3)+symfunction_list_local(i3,i1)
        enddo ! i3
!!
        do i3=1,matrixsize
          hextoffmin(i3)=min(hextoffmin(i3),hextoffelem_local(i3,i1))
          hextoffmax(i3)=max(hextoffmax(i3),hextoffelem_local(i3,i1))
        enddo
      enddo ! i1
!!
      if(ncount.gt.0) goto 10
!!
!! final calculation of average values
      do i2=1,num_funcvalues_local
        avvalue_local(i2)=avvalue_local(i2)/dble(counter)
      enddo ! i2

!!     
!! Now calculate standard deviations
      rewind(symhextoffunit)
      ncount=ntrain
 20   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! read a block of points
      do i1=1,npoints
        read(symhextoffunit,*) (symfunction_list_local(i3,i1),&
          i3=1,num_funcvalues_local)
        read(symhextoffunit,*)zdummy
      enddo ! i1
!!
      do i1=1,npoints
        do i3=1,num_funcvalues_local
          stddev_local(i3)&
            =stddev_local(i3)+(symfunction_list_local(i3,i1)-avvalue_local(i3))**2
        enddo ! i3
      enddo ! i1
!!
      if(ncount.gt.0) goto 20
!!
!! calculate final part of stddev
      do i2=1,num_funcvalues_local
        stddev_local(i2)=dsqrt(stddev_local(i2)/dble(counter))
      enddo ! i2
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
      write(ounit,*)'============================================================='
        write(ounit,'(a,3a)')' Hextoff symmetry function values for triplet ',&
          element(elementindex(hextoff_training_triplet(1))),&
          element(elementindex(hextoff_training_triplet(2))),&
          element(elementindex(hextoff_training_triplet(3)))
!!
      if(counter.eq.0)then
        write(ounit,*)'WARNING: No triplet present!'
        avvalue_local(:)=0.0d0
      else
        write(ounit,'(a)')' Training set:  min           max       average         range        stddev      range/stddev'
        do i2=1,num_funcvalues_local
          write(ounit,'(i4,x,6f14.8)')i2,minvalue_local(i2),maxvalue_local(i2),&
            avvalue_local(i2),abs(maxvalue_local(i2)-minvalue_local(i2)),&
            stddev_local(i2),abs(maxvalue_local(i2)-minvalue_local(i2))/stddev_local(i2)
          if(abs(minvalue_local(i2)-maxvalue_local(i2)).lt.thres)then
            write(ounit,*)'### WARNING ###: minvalue=maxvalue',i2
            if(lscalesym)then
              write(ounit,*)'scaling symmetry functions cannot be used with minvalue=maxvalue'
              stop
            endif
          endif
        enddo ! i2
      endif
!!
!! determine Pearson's correlation of all pairs of symmetry functions if requested
      if(lpearson_correlation)then
        rewind(symhextoffunit)
        ncount=ntrain
 30     continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!!
!! read a block of points
        do i1=1,npoints
          read(symhextoffunit,*) (symfunction_list_local(i3,i1),&     
            i3=1,num_funcvalues_local)
          read(symhextoffunit,*)zdummy
        enddo ! i1
!!
        do i1=1,npoints
!! get atom or pair type
!! calculate first part of pearson correlation 
          do i3=1,num_funcvalues_local
            do i4=1,num_funcvalues_local
              pearson_corr(i3,i4)&
              =pearson_corr(i3,i4)&
                +(symfunction_list_local(i3,i1)-avvalue_local(i3))&
                *(symfunction_list_local(i4,i1)-avvalue_local(i4))
            enddo ! i4
          enddo ! i3
        enddo ! i1
!!
        if(ncount.gt.0) goto 30
!!
!! calculate final part of pearson correlation 
        do i2=1,num_funcvalues_local
          do i3=1,num_funcvalues_local
            pearson_corr(i2,i3)=pearson_corr(i2,i3)& 
            /(stddev_local(i2)*stddev_local(i3)*counter)
          enddo ! i3
        enddo ! i2
!!
!! write pearson correlation
        write(ounit,*)'============================================================='
        write(ounit,*)'Pearson correlation of symmetry functions'
            write(ounit,'(a,a2)')' Hextoff symmetry functions for triplets ',&
           element(elementindex(hextoff_training_triplet(1))),&
           element(elementindex(hextoff_training_triplet(2))),&
           element(elementindex(hextoff_training_triplet(3)))

        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(8x,300i8)')(i2,i2=1,num_funcvalues_local)
        do i2=1,num_funcvalues_local
          write(ounit,'(a4,i4,300f8.4)')'CORR',i2,(pearson_corr(i2,i3),i3=1,num_funcvalues_local)
        enddo ! i2
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'============================================================='
!!
!! get correlation statistics
        do i2=1,num_funcvalues_local
          do i3=1,num_funcvalues_local
            if(pearson_corr(i2,i3).gt.0.99d0)then
              correlation_stat(i2,1)=correlation_stat(i2,1)+1
            elseif(pearson_corr(i2,i3).gt.0.98d0)then
              correlation_stat(i2,2)=correlation_stat(i2,2)+1
            elseif(pearson_corr(i2,i3).gt.0.95d0)then
              correlation_stat(i2,3)=correlation_stat(i2,3)+1
            elseif(pearson_corr(i2,i3).gt.0.90d0)then
              correlation_stat(i2,4)=correlation_stat(i2,4)+1
            elseif(pearson_corr(i2,i3).gt.0.80d0)then
              correlation_stat(i2,5)=correlation_stat(i2,5)+1
            elseif(pearson_corr(i2,i3).gt.0.70d0)then
              correlation_stat(i2,6)=correlation_stat(i2,6)+1
            elseif(pearson_corr(i2,i3).gt.0.60d0)then
              correlation_stat(i2,7)=correlation_stat(i2,7)+1
            elseif(pearson_corr(i2,i3).gt.0.50d0)then
              correlation_stat(i2,8)=correlation_stat(i2,8)+1
            elseif(pearson_corr(i2,i3).lt.0.00d0)then
              correlation_stat(i2,10)=correlation_stat(i2,10)+1
            else
              correlation_stat(i2,9)=correlation_stat(i2,9)+1
            endif 
          enddo ! i3
        enddo ! i2
!!
!! remove self-correlation
        do i2=1,num_funcvalues_local
          correlation_stat(i2,1)=correlation_stat(i2,1)-1
        enddo ! i2
!!
!! write correlation_stat
        write(ounit,'(a,a2)')' Hextoff symmetry function correlation statistics for triplet ',&
          element(elementindex(hextoff_training_triplet(1))),&
          element(elementindex(hextoff_training_triplet(2))),&
          element(elementindex(hextoff_training_triplet(3)))

        write(ounit,'(15x,10a8)')'   >0.99','   >0.98','   >0.95','   >0.90',&
          '   >0.80','   >0.70','   >0.60','   >0.50','   >0.00','   <0.00'
        write(ounit,*)'-------------------------------------------------------------'
        do i2=1,num_funcvalues_local
          write(ounit,'(a9,x,i4,x,10i8)')' CORRSTAT',i2,(correlation_stat(i2,i3),i3=1,10)
        enddo ! i2
        write(ounit,*)'-------------------------------------------------------------'
!!
      endif !lpearson_correlation
!!
      return
      end
