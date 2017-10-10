!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - fitting.f90 (2x)
!! - fitting_batch.f90 (2x)
!! - fittingpair.f90 (2x)
!!
      subroutine getscale(ndim1,ndim2,iswitch,maxnum_funcvalues_local,num_funcvalues_local,&
        ntrain,symfunction_list_local,minvalue_local,maxvalue_local,avvalue_local)
!!
      use fileunits
      use globaloptions
      use structures
!!
      implicit none
!!
      integer maxnum_funcvalues_local                     ! in
      integer num_funcvalues_local(ndim1)                 ! in
      integer ndim1                                       ! in (=nelem or npairs)
      integer ndim2                                       ! in (=max_num_atoms or max_num_pairs)
      integer icount                                      ! internal
      integer jcount                                      ! internal
      integer npoints                                     ! internal
      integer ncount                                      ! internal
      integer ntrain                                      ! in
      integer counter(ndim1)                              ! internal
      integer i1,i2,i3,i4                                 ! internal
      integer iswitch                                     ! in 
      integer correlation_stat(ndim1,maxnum_funcvalues_local,10)           ! internal
!!
      real*8 symfunction_list_local(maxnum_funcvalues_local,ndim2,nblock)
!!
      real*8 minvalue_local(ndim1,maxnum_funcvalues_local)                 ! out
      real*8 maxvalue_local(ndim1,maxnum_funcvalues_local)                 ! out
      real*8 avvalue_local(ndim1,maxnum_funcvalues_local)                  ! out
      real*8 stddev_local(ndim1,maxnum_funcvalues_local)                   ! internal
      real*8 emin
      real*8 emax
      real*8 eshortmin
      real*8 eshortmax
      real*8 eewaldmin
      real*8 eewaldmax
      real*8 thres
      real*8 zdummy
      real*8 pearson_corr(ndim1,maxnum_funcvalues_local,maxnum_funcvalues_local) ! internal
!!
!! initialization
      minvalue_local(:,:)  = 100000.d0
      maxvalue_local(:,:)  =-100000.d0
      avvalue_local(:,:)   = 0.d0
      stddev_local(:,:)    = 0.d0
      pearson_corr(:,:,:)  = 0.0d0
      emin                 = 10000.d10
      emax                 =-10000.d10
      eshortmin            = 10000.d10
      eshortmax            =-10000.d10
      eewaldmin            = 10000.d10
      eewaldmax            =-10000.d10
      counter(:)           =0
      correlation_stat(:,:,:)=0
!!      thres=0.00000001d0 ! CAUTION: If this is combined with scaling, noise might become very large!!!
      thres        =0.0001d0
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
        if(iswitch.eq.0)then   ! atomic short range NN
          read(symunit,*)num_atoms_list(i1)
          do i2=1,num_atoms_list(i1)
            read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
              i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
          enddo ! i2           ! atomic electrostatic NN
        elseif(iswitch.eq.1)then
          read(symunit,*)num_atoms_list(i1)
          do i2=1,num_atoms_list(i1)
            read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
              i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
          enddo ! i2
        elseif(iswitch.eq.2)then ! pair short range NN
          read(symunit,*)num_atoms_list(i1),num_pairs_list(i1)
          do i2=1,num_pairs_list(i1)
            read(symunit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),&
              (symfunction_list_local(i3,i2,i1),i3=1,num_funcvalues_local(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
          enddo ! i2
        else
          write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
          stop
        endif
        read(symunit,*)totalcharge_list(i1),totalenergy_list(i1),&
          shortenergy_list(i1),elecenergy_list(i1)
      enddo ! i1
!!
!! modify minvalue_local and maxvalue_local
      do i1=1,npoints
!! get number of symfunction lines icount:
        if(iswitch.eq.0)then
          icount=num_atoms_list(i1)
        elseif(iswitch.eq.1)then
          icount=num_atoms_list(i1)
        elseif(iswitch.eq.2)then
          icount=num_pairs_list(i1)
        else
          write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
          stop
        endif
!! get atom or pair type
        do i2=1,icount
          if(iswitch.eq.0)then
            jcount=elementindex(zelem_list(i1,i2))
          elseif(iswitch.eq.1)then
            jcount=elementindex(zelem_list(i1,i2))
          elseif(iswitch.eq.2)then
            jcount=pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
!! count for later averaging
          counter(jcount)&
            =counter(jcount)+1
          do i3=1,num_funcvalues_local(jcount)
            minvalue_local(jcount,i3)&
              =min(minvalue_local(jcount,i3),symfunction_list_local(i3,i2,i1))
            maxvalue_local(jcount,i3)&
              =max(maxvalue_local(jcount,i3),symfunction_list_local(i3,i2,i1))
            avvalue_local(jcount,i3) &
              =avvalue_local(jcount,i3)+symfunction_list_local(i3,i2,i1)
          enddo ! i3
        enddo ! i2
!!
        emin     =min(emin,totalenergy_list(i1))
        emax     =max(emax,totalenergy_list(i1))
        eshortmin=min(eshortmin,shortenergy_list(i1))
        eshortmax=max(eshortmax,shortenergy_list(i1))
        eewaldmin=min(eewaldmin,elecenergy_list(i1))
        eewaldmax=max(eewaldmax,elecenergy_list(i1))
      enddo ! i1
!!
      if(ncount.gt.0) goto 10
!!
!! final calculation of average values
      do i1=1,ndim1
        do i2=1,num_funcvalues_local(i1)
          avvalue_local(i1,i2)=avvalue_local(i1,i2)/dble(counter(i1))
        enddo ! i2
      enddo ! i1
!!     
!! Now calculate standard deviations
      rewind(symunit)
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
        if(iswitch.eq.0)then   ! atomic short range NN
          read(symunit,*)num_atoms_list(i1)
          do i2=1,num_atoms_list(i1)
            read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
              i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
          enddo ! i2           ! atomic electrostatic NN
        elseif(iswitch.eq.1)then
          read(symunit,*)num_atoms_list(i1)
          do i2=1,num_atoms_list(i1)
            read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
              i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
          enddo ! i2
        elseif(iswitch.eq.2)then ! pair short range NN
          read(symunit,*)num_atoms_list(i1),num_pairs_list(i1)
          do i2=1,num_pairs_list(i1)
            read(symunit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),&
              (symfunction_list_local(i3,i2,i1),i3=1,num_funcvalues_local(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
          enddo ! i2
        else
          write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
          stop
        endif
        read(symunit,*)zdummy
      enddo ! i1
!!
      do i1=1,npoints
!! get number of symfunction lines icount:
        if(iswitch.eq.0)then
          icount=num_atoms_list(i1)
        elseif(iswitch.eq.1)then
          icount=num_atoms_list(i1)
        elseif(iswitch.eq.2)then
          icount=num_pairs_list(i1)
        else
          write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
          stop
        endif
!! get atom or pair type
        do i2=1,icount
          if(iswitch.eq.0)then
            jcount=elementindex(zelem_list(i1,i2))
          elseif(iswitch.eq.1)then
            jcount=elementindex(zelem_list(i1,i2))
          elseif(iswitch.eq.2)then
            jcount=pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
!! calculate first part of stddev
          do i3=1,num_funcvalues_local(jcount)
            stddev_local(jcount,i3)&
              =stddev_local(jcount,i3)+(symfunction_list_local(i3,i2,i1)-avvalue_local(jcount,i3))**2
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
      if(ncount.gt.0) goto 20
!!
!! calculate final part of stddev
      do i1=1,ndim1
        do i2=1,num_funcvalues_local(i1)
          stddev_local(i1,i2)=dsqrt(stddev_local(i1,i2)/dble(counter(i1)))
        enddo ! i2
      enddo ! i1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
      do i1=1,ndim1
        write(ounit,*)'============================================================='
        if(iswitch.eq.0)then
          write(ounit,'(a,a2)')' Short range symmetry function values for element ',element(i1)
        elseif(iswitch.eq.1)then
          write(ounit,'(a,a2)')' Electrostatic symmetry function values for element ',element(i1)
        elseif(iswitch.eq.2)then
          write(ounit,'(a,2a3)')' Pair symmetry function values for element pair ',&
            element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2)))
        else
          write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
          stop
        endif
!!
        if(counter(i1).eq.0)then
          write(ounit,*)'WARNING: No species of this type present!'
          avvalue_local(i1,:)=0.0d0
        else
          write(ounit,'(a)')' Training set:  min           max       average         range        stddev      range/stddev'
          do i2=1,num_funcvalues_local(i1)
            write(ounit,'(i4,x,6f14.8)')i2,minvalue_local(i1,i2),maxvalue_local(i1,i2),&
              avvalue_local(i1,i2),abs(maxvalue_local(i1,i2)-minvalue_local(i1,i2)),&
              stddev_local(i1,i2),abs(maxvalue_local(i1,i2)-minvalue_local(i1,i2))/stddev_local(i1,i2)
            if(abs(minvalue_local(i1,i2)-maxvalue_local(i1,i2)).lt.thres)then
              if(iswitch.eq.0)then
                write(ounit,*)'### WARNING ###: minvalue=maxvalue',i1,i2,nucelem(i1)
              elseif(iswitch.eq.1)then
                write(ounit,*)'### WARNING ###: minvalue=maxvalue ',i1,i2,nucelem(i1)
              elseif(iswitch.eq.2)then
                write(ounit,*)'### WARNING ###: minvalue=maxvalue ',i1,i2,&
                  nucelem(elementindex(elempair(i1,1))),nucelem(elementindex(elempair(i1,2)))
              else
                write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
                stop
              endif
              if(lscalesym)then
                write(ounit,*)'scaling symmetry functions cannot be used with minvalue=maxvalue'
                stop
              endif
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! determine Pearson's correlation of all pairs of symmetry functions if requested
      if(lpearson_correlation)then
        rewind(symunit)
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
          if(iswitch.eq.0)then   ! atomic short range NN
            read(symunit,*)num_atoms_list(i1)
            do i2=1,num_atoms_list(i1)
              read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
                i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
            enddo ! i2           ! atomic electrostatic NN
          elseif(iswitch.eq.1)then
            read(symunit,*)num_atoms_list(i1)
            do i2=1,num_atoms_list(i1)
              read(symunit,*)zelem_list(i1,i2),(symfunction_list_local(i3,i2,i1),&
                i3=1,num_funcvalues_local(elementindex(zelem_list(i1,i2))))
            enddo ! i2
          elseif(iswitch.eq.2)then ! pair short range NN
            read(symunit,*)num_atoms_list(i1),num_pairs_list(i1)
            do i2=1,num_pairs_list(i1)
              read(symunit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),&
                (symfunction_list_local(i3,i2,i1),i3=1,num_funcvalues_local(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
            enddo ! i2
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
          read(symunit,*)zdummy
        enddo ! i1
!!
        do i1=1,npoints
!! get number of symfunction lines icount:
          if(iswitch.eq.0)then
            icount=num_atoms_list(i1)
          elseif(iswitch.eq.1)then
            icount=num_atoms_list(i1)
          elseif(iswitch.eq.2)then
            icount=num_pairs_list(i1)
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
!! get atom or pair type
          do i2=1,icount
            if(iswitch.eq.0)then
              jcount=elementindex(zelem_list(i1,i2))
            elseif(iswitch.eq.1)then
              jcount=elementindex(zelem_list(i1,i2))
            elseif(iswitch.eq.2)then
              jcount=pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))
            else
              write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
              stop
            endif
!! calculate first part of pearson correlation 
            do i3=1,num_funcvalues_local(jcount)
              do i4=1,num_funcvalues_local(jcount)
                pearson_corr(jcount,i3,i4)&
                =pearson_corr(jcount,i3,i4)&
                  +(symfunction_list_local(i3,i2,i1)-avvalue_local(jcount,i3))&
                  *(symfunction_list_local(i4,i2,i1)-avvalue_local(jcount,i4))
              enddo ! i4
            enddo ! i3
          enddo ! i2
        enddo ! i1
!!
        if(ncount.gt.0) goto 30
!!
!! calculate final part of pearson correlation 
        do i1=1,ndim1
          do i2=1,num_funcvalues_local(i1)
            do i3=1,num_funcvalues_local(i1)
              pearson_corr(i1,i2,i3)=pearson_corr(i1,i2,i3)& 
              /(stddev_local(i1,i2)*stddev_local(i1,i3)*counter(i1))
            enddo ! i3
          enddo ! i2
        enddo ! i1
!!
!!
!! write pearson correlation
        write(ounit,*)'============================================================='
        write(ounit,*)'Pearson correlation of symmetry functions'
        do i1=1,ndim1
          if(iswitch.eq.0)then
            write(ounit,'(a,a2)')' Short range symmetry functions for element ',element(i1)
          elseif(iswitch.eq.1)then
            write(ounit,'(a,a2)')' Electrostatic symmetry functions for element ',element(i1)
          elseif(iswitch.eq.2)then
            write(ounit,'(a,2a3)')' Pair symmetry functions for element pair ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2)))
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(8x,300i8)')(i2,i2=1,num_funcvalues_local(i1))
          do i2=1,num_funcvalues_local(i1)
            write(ounit,'(a4,i4,300f8.4)')'CORR',i2,(pearson_corr(i1,i2,i3),i3=1,num_funcvalues_local(i1))
          enddo ! i2
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
        write(ounit,*)'============================================================='
!!
!! get correlation statistics
        do i1=1,ndim1
          do i2=1,num_funcvalues_local(i1)
            do i3=1,num_funcvalues_local(i1)
!! JB: 2015 05 09: negative correlation is now treated like positive correlation
!!              if(pearson_corr(i1,i2,i3).gt.0.99d0)then
              if(abs(pearson_corr(i1,i2,i3)).gt.0.99d0)then
                correlation_stat(i1,i2,1)=correlation_stat(i1,i2,1)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.98d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.98d0)then
                correlation_stat(i1,i2,2)=correlation_stat(i1,i2,2)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.95d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.95d0)then
                correlation_stat(i1,i2,3)=correlation_stat(i1,i2,3)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.90d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.90d0)then
                correlation_stat(i1,i2,4)=correlation_stat(i1,i2,4)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.80d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.80d0)then
                correlation_stat(i1,i2,5)=correlation_stat(i1,i2,5)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.70d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.70d0)then
                correlation_stat(i1,i2,6)=correlation_stat(i1,i2,6)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.60d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.60d0)then
                correlation_stat(i1,i2,7)=correlation_stat(i1,i2,7)+1
!!              elseif(pearson_corr(i1,i2,i3).gt.0.50d0)then
              elseif(abs(pearson_corr(i1,i2,i3)).gt.0.50d0)then
                correlation_stat(i1,i2,8)=correlation_stat(i1,i2,8)+1
!!              elseif(pearson_corr(i1,i2,i3).lt.0.00d0)then
!!                correlation_stat(i1,i2,10)=correlation_stat(i1,i2,10)+1
              else
                correlation_stat(i1,i2,9)=correlation_stat(i1,i2,9)+1
              endif 
            enddo ! i3
          enddo ! i2
        enddo ! i1
!!
!! remove self-correlation
        do i1=1,ndim1
          do i2=1,num_funcvalues_local(i1)
            correlation_stat(i1,i2,1)=correlation_stat(i1,i2,1)-1
          enddo ! i2
        enddo ! i1
!!
!! write correlation_stat
        do i1=1,ndim1
          if(iswitch.eq.0)then
            write(ounit,'(a,a2)')' Short range symmetry function correlation statistics for element ',element(i1)
          elseif(iswitch.eq.1)then
            write(ounit,'(a,a2)')' Electrostatic symmetry function correlation statistics for element ',element(i1)
          elseif(iswitch.eq.2)then
            write(ounit,'(a,2a3)')' Pair symmetry function correlation statistics for element pair ',&
              element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2)))
          else
            write(ounit,*)'ERROR: Unknown iswitch in getscale ',iswitch
            stop
          endif
          write(ounit,'(15x,10a8)')'   >0.99','   >0.98','   >0.95','   >0.90',&
            '   >0.80','   >0.70','   >0.60','   >0.50','   >0.00' !,'   <0.00'
!!            '   >0.80','   >0.70','   >0.60','   >0.50','   >0.00','   <0.00'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_local(i1)
            write(ounit,'(a9,x,i4,x,10i8)')' CORRSTAT',i2,(correlation_stat(i1,i2,i3),i3=1,9)
!!            write(ounit,'(a9,x,i4,x,10i8)')' CORRSTAT',i2,(correlation_stat(i1,i2,i3),i3=1,10) !! entry 10 is for negative correlation
          enddo ! i2
          write(ounit,*)'-------------------------------------------------------------'
        enddo ! i1
!!
      endif !lpearson_correlation
!!
      return
      end
