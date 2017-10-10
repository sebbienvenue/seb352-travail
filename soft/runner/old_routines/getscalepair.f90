!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90 
!!
      subroutine getscalepair(nelem,npairs,elempair,&
        nblock,maxnum_funcvaluesp,num_funcvaluesp,&
        max_num_atoms,max_num_pairs,&
        elementindex,pairindex,&
        ntrain,nucelem,&
        symfunctionp_list,minvaluep,maxvaluep,avvaluep,&
        element,lscalesym,lshort,lewald,ldebug)
!!
      use fileunits
!!
      implicit none
!!
      integer maxnum_funcvalues                          ! in
      integer maxnum_funcvaluesp                         ! in
      integer num_funcvalues(nelem)                      ! in
      integer num_funcvaluesp(npairs)
      integer nelem
      integer npairs
      integer nblock
      integer max_num_atoms_
      integer max_num_pairs
      integer npoints
      integer ncount
      integer ntrain
      integer elementindex(102)
      integer pairindex(102,102)
      integer zelemp_list(2,nblock,max_num_pairs)
      integer num_atoms_list(nblock)
      integer size_pairs_list(nblock)
      integer nucelem(nelem)
      integer elemcount(nelem)
      integer paircount(npairs)
      integer i1,i2,i3
      integer elempair(npairs,2)                          ! in
!!
      real*8 symfunctionp_list(maxnum_funcvaluesp,max_num_pairs,nblock)
      real*8 minvaluep(npairs,maxnum_funcvaluesp)         ! out
      real*8 maxvaluep(npairs,maxnum_funcvaluesp)         ! out
      real*8 avvaluep(npairs,maxnum_funcvaluesp)          ! out
      real*8 totalcharge_list(nblock)
      real*8 totalenergy_list(nblock)
      real*8 shortenergy_list(nblock)
      real*8 ewaldenergy_list(nblock)
      real*8 emin
      real*8 emax
      real*8 eshortmin
      real*8 eshortmax
      real*8 eewaldmin
      real*8 eewaldmax
      real*8 thres
!!
      character*2 element(nelem)                           ! in
!!
      logical lscalesym
      logical lshort ! in
      logical lewald ! in
      logical ldebug
!!
!! initialization
      minvaluep(:,:)=100000.d0
      maxvaluep(:,:)=0.d0
      avvaluep(:,:)=0.d0
      emin= 10000.d10
      emax=-10000.d10
      eshortmin= 10000.d10
      eshortmax=-10000.d10
      eewaldmin= 10000.d10
      eewaldmax=-10000.d10
      elemcount(:)=0
      paircount(:)=0
      thres=0.00000001d0
!!
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
        read(symunit,*)num_atoms_list(i1),size_pairs_list(i1)
        do i2=1,size_pairs_list(i1)  
          read(symunit,*)zelemp_list(1,i1,i2),zelemp_list(2,i1,i2),(symfunctionp_list(i3,i2,i1)&
            ,i3=1,num_funcvaluesp(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))))
        enddo ! i2
        read(symunit,*)totalcharge_list(i1),totalenergy_list(i1),shortenergy_list(i1),ewaldenergy_list(i1)
      enddo ! i1

!!
!! modify minvalue and maxvalue
      do i1=1,npoints
        do i2=1,size_pairs_list(i1)
          paircount(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))&
            =paircount(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)))+1
          do i3=1,num_funcvaluesp(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2))) 
            minvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3)&
                =min(minvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3),symfunctionp_list(i3,i2,i1))
            maxvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3)&
                =max(maxvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3),symfunctionp_list(i3,i2,i1))
            avvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3)&
                =avvaluep(pairindex(zelemp_list(1,i1,i2),zelemp_list(2,i1,i2)),i3)+symfunctionp_list(i3,i2,i1)
          enddo ! i3
        enddo ! i2
        emin=min(emin,totalenergy_list(i1))
        emax=max(emax,totalenergy_list(i1))
        eshortmin=min(eshortmin,shortenergy_list(i1))
        eshortmax=max(eshortmax,shortenergy_list(i1))
        eewaldmin=min(eewaldmin,ewaldenergy_list(i1))
        eewaldmax=max(eewaldmax,ewaldenergy_list(i1))
      enddo ! i1
!!
      if(ncount.gt.0) goto 10
!!
!! final calculation of average values
      do i1=1,npairs  
        do i2=1,num_funcvaluesp(i1)
          avvaluep(i1,i2)=avvaluep(i1,i2)/dble(paircount(i1))
        enddo ! i2
      enddo ! i1
!!      
      do i1=1,npairs 

        write(ounit,*)'============================================================='
        write(ounit,'(a,2a3)')' Pair symmetry function values for element pair ',&
          element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2)))
        if(paircount(i1).eq.0)then
          write(ounit,*)'No pairs present in data set'
          avvaluep(i1,:)=0.0d0
        else
          write(ounit,*)'Training set:  min           max       average         range '
          do i2=1,num_funcvaluesp(i1)
            write(ounit,'(i4,x,4f14.8)')i2,minvaluep(i1,i2),maxvaluep(i1,i2),&
              avvaluep(i1,i2),abs(maxvaluep(i1,i2)-minvaluep(i1,i2))
            if(abs(minvaluep(i1,i2)-maxvaluep(i1,i2)).lt.thres)then
              write(ounit,*)'### WARNING ###: minvalue=maxvalue ',i1,i2,nucelem(elementindex(elempair(i1,1))),nucelem(elementindex(elempair(i1,2)))
              if(lscalesym)then
                write(ounit,*)'scaling symmetry functions cannot be used with minvalue=maxvalue'
                stop
              endif
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! write energy ranges only once:
!!      if((lshort.and.(.not.lewald)).or.(lshort.and.lewald.and.(iswitch.eq.1)))then
!!        write(ounit,*)'============================================================='
!!        write(ounit,'(a21,x,f14.8)')' Emin      (Ha/atom) ',emin
!!        write(ounit,'(a21,x,f14.8)')' Emax      (Ha/atom) ',emax
!!        write(ounit,'(a21,x,f14.8)')' Eshortmin (Ha/atom) ',eshortmin
!!        write(ounit,'(a21,x,f14.8)')' Eshortmax (Ha/atom) ',eshortmax
!!        write(ounit,'(a21,x,f14.8)')' Eewaldmin (Ha/atom) ',eewaldmin
!!        write(ounit,'(a21,x,f14.8)')' Eewaldmax (Ha/atom) ',eewaldmax
!!      endif
!!
      return
      end
