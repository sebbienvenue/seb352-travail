!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
!! This routine is called once for each block of points
!!
      subroutine optimize_short(ounit,nblock,npoints,point,idx,&
         nelem,num_weightsshort,max_num_atoms,num_funcvalues,&
         num_layersshort,windex,zelem_list,maxnodes_short,&
         nenergygroup,lseed,&
         num_atoms_list,nodes_short,elementindex,optmode,&
         kalmanthreshold,kalmanlambda,kalmannue,energyrnd,&
         symfunction_list,weights_short,&
         rmse_short,shortenergy_list,&
         corrmatrix_list,kalgainmat_list,&
         actfunc_short,&
         ldebug)
!!
      use mpi_mod
      use kalmandims
!!
      implicit none
!!
      integer ounit                                 ! in
      integer nblock                                ! in
      integer npoints                               ! in
      integer nelem                                 ! in
      integer num_weightsshort                      ! in
      integer max_num_atoms                         ! in
      integer num_funcvalues                        ! in
      integer num_layersshort                       ! in
      integer windex(2*num_layersshort)             ! in
      integer elementindex(102)                     ! in
      integer nodes_short(0:num_layersshort)        ! in
      integer zelem_list(nblock,max_num_atoms)      ! in 
      integer zelem(max_num_atoms)                  ! internal
      integer num_atoms_list(nblock)                ! in
      integer num_atoms                             ! internal
      integer maxnodes_short                        ! in
      integer num_atoms_element(nelem)              ! internal 
      integer optmode                               ! in
      integer point                                 ! in/out
      integer idx(nblock)                           ! in
      integer lseed                                 ! in/out
      integer i1,i2,i3                              ! internal
      integer day                                   ! internal
      integer nenergygroup                          ! in
      integer nenergy                               ! internal
!!
      real*8 kalmanthreshold                        ! in
      real*8 kalmanlambda                           ! in
      real*8 kalmannue                              ! in
      real*8 rmse_short                             ! in
      real*8 eshort                                 ! internal
      real*8 error                                  ! internal
      real*8 errorsum                                  ! internal
      real*8 abserror                               ! internal
      real*8 nnatomenergy(max_num_atoms)            ! internal
      real*8 weights_short(num_weightsshort,nelem)                          ! in/out
      real*8 weights(num_weightsshort)                                      ! internal 
      real*8 symfunction_list(num_funcvalues,max_num_atoms,nblock)          ! in
      real*8 symfunction(num_funcvalues,max_num_atoms)                      ! internal
      real*8 shortenergy_list(nblock)                                       ! in
      real*8 deshortdw(num_weightsshort,nodes_short(num_layersshort),nelem) ! internal
      real*8 deshortdwsum(num_weightsshort,nodes_short(num_layersshort),nelem) ! internal
      real*8 dedw(num_weightsshort)                                         ! internal
      real*8 corrmatrix(corrdim)                                            ! internal
      real*8 corrmatrix_list(corrdim,nelem)                                 ! in/out
      real*8 kalgainmat(kaldim)                                             ! internal
      real*8 kalgainmat_list(kaldim,nelem)                                  ! in/out
      real*8 energyrnd                                                      ! in
      real*8 z,ran0                                                         ! internal
      real*8 timestart
      real*8 timeend
      real*8 timesumkal
      real*8 timesumderiv
!!
      character*1 actfunc_short(num_layersshort)                            ! in
!!
      logical ldebug                                                        ! in
!!
!! initialization
       deshortdw(:,:,:)    = 0.0d0
       deshortdwsum(:,:,:) = 0.0d0
       eshort              = 0.0d0
       timesumkal          = 0.0d0
       timesumderiv        = 0.0d0
       day                 = 0
       nenergy             = 0
       errorsum            = 0.0d0
       num_atoms_element(:)= 0
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! debugging
!!      write(*,*)kalmanthreshold,rmse_short*kalmanthreshold
!!
!! loop over all points
      do i1=1,npoints
        point=point+1
!!
        z=ran0(lseed)
        if(z.lt.energyrnd)then
          write(ounit,*)'optimize_short: point number ',i1
!!
!! predict the short range energy of point i1
        num_atoms=num_atoms_list(idx(i1))
!!        symfunction(:,:)=symfunction_list(:,:,idx(i1))
        zelem(:)=zelem_list(idx(i1),:)
!!
      call calconeshort(ounit,max_num_atoms,num_atoms,num_funcvalues,&
        nelem,num_weightsshort,maxnodes_short,&
        num_layersshort,nodes_short,elementindex,zelem,&
        windex,symfunction_list(1,1,idx(i1)),weights_short,eshort,nnatomenergy,&
        actfunc_short,ldebug)
!! normalize per atom
        eshort=eshort/dble(num_atoms)
!!
!! Caution: Using eshort-shortenergy_list(i1) instead requires sign change of derivatives dedw
        abserror=abs(shortenergy_list(idx(i1))-eshort)
        error=shortenergy_list(idx(i1))-eshort
!!
        if(abserror.gt.kalmanthreshold*rmse_short)then
          nenergy=nenergy+1
!!
!!          write(*,'(i4,a,3f18.6)')point,' point will be used to update short range weights ',&
!!           shortenergy_list(idx(i1)),eshort,error !'
!!
!! calculate the derivative of the short range energy with respect to the weights
!!         call abstime(timestart,day)
         call getdeshortdw(ounit,nelem,num_weightsshort,&
           zelem,elementindex,max_num_atoms,num_atoms,&
           maxnodes_short,&
           windex,num_layersshort,num_funcvalues,nodes_short,&
           symfunction_list(1,1,idx(i1)),&
           weights_short,deshortdw,&
           actfunc_short,&
           ldebug)
!!
         deshortdwsum(:,:,:)=deshortdwsum(:,:,:)+deshortdw(:,:,:)
         errorsum=errorsum+error
!!          call abstime(timeend,day)
!!          timesumderiv=timesumderiv+timeend-timestart
!!
!! calculate the number of atoms per element, because we don't want to 
!! update weights for elements not being present in current example
!!         num_atoms_element(:) = 0
         do i2=1,num_atoms
           num_atoms_element(elementindex(zelem(i2)))=&
             num_atoms_element(elementindex(zelem(i2)))+1
         enddo ! i2
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!         if(mpirank.eq.0)then
         if((mod(nenergy,nenergygroup).eq.0).or.&
         (i1.eq.npoints))then
!!
         error=errorsum/dble(nenergy)
         deshortdw(:,:,:)=deshortdwsum(:,:,:)/dble(nenergy)
!!
!!         write(ounit,*)'weight update is done'
!!
         if(optmode.eq.1)then ! Kalman filter
         do i2=1,nelem
           if(num_atoms_element(i2).gt.0)then ! atoms of element are present
              corrmatrix(:)=corrmatrix_list(:,i2)
              kalgainmat(:)=kalgainmat_list(:,i2)
              weights(:)   =weights_short(:,i2)
              dedw(:)      =deshortdw(:,1,i2)
!!          call abstime(timestart,day)
              call updatekalman(ounit,kaldim,corrdim,&
                num_weightsshort,kalmanlambda,kalmannue,&
                weights,dedw,kalgainmat,corrmatrix,error,ldebug)
!!          call abstime(timeend,day)
!!          timesumkal=timesumkal+timeend-timestart
              corrmatrix_list(:,i2)=corrmatrix(:)
              kalgainmat_list(:,i2)=kalgainmat(:)
              weights_short(:,i2)=weights(:)
           endif
         enddo
         elseif(optmode.eq.2)then ! conjugate gradient'
           write(ounit,*)'CG optimization not yet implemented in optimize_short'
           stop
         elseif(optmode.eq.3)then ! steepest descent'
         do i2=1,nelem
           if(num_atoms_element(i2).gt.0)then ! atoms of element are present
              weights(:)   =weights_short(:,i2)
              dedw(:)      =deshortdw(:,1,i2)
              call updatesteepest(ounit,&
                num_weightsshort,weights,dedw,error,ldebug)
              weights_short(:,i2)=weights(:)
           endif
         enddo
         endif ! optmode
!!
!! reinitializations for next update
           nenergy             = 0
           deshortdwsum(:,:,:) = 0.0d0
           errorsum            = 0.0d0
          num_atoms_element(:) = 0
!!
         endif !((mod(nenergy,nenergygroup).eq.0).or.(i1.eq.npoints))then
!!
!!         endif ! mpirank.eq.0
!!
        else
!!          write(*,'(i6,a,3f18.8)')point,' point is not used for short update ',&
!!          shortenergy_list(i1),eshort,error
        endif
!!
        endif !(z.lt.energyrnd)then
      enddo ! i1=1,npoints
!!
!!      write(ounit,'(a,f14.6,i8)')'deriv time ',timesumderiv,mpirank
!!      write(ounit,'(a,f14.6,i8)')'kalman time ',timesumkal,mpirank
!!
      return
      end
