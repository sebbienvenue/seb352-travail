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
      subroutine optimize_short_para2(ounit,nblock,npoints,point,idx,&
         nelem,num_weightsshort,max_num_atoms,num_funcvalues,&
         num_layersshort,windex,zelem_list,maxnodes_short,&
         nenergygroup,lseed,paramode,&
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
      integer jdx(nblock)                           ! in
      integer kdx(nblock)                           ! in
      integer mdx(nblock)                           ! in
      integer lseed                                 ! in/out
      integer i1,i2,i3,i4                           ! internal
      integer day                                   ! internal
      integer nenergygroup                          ! in
      integer nenergy                               ! internal
      integer nstruct                               ! internal
      integer n_start                               ! internal
      integer n_end                                 ! internal
      integer, dimension(:), allocatable :: atomindex ! internal
      integer icount                                ! internal
      integer jcount                                ! internal
      integer kcount                                ! internal
      integer ncount                                ! internal
      integer npts                                  ! internal
      integer nweights                              ! internal
      integer paramode                              ! in
      integer nbad                                  ! internal
      integer nbad_sum                              ! internal
      integer num_points                            ! internal
!!
      real*8 kalmanthreshold                        ! in
      real*8 kalmanlambda(nelem)                    ! in
      real*8 kalmannue                              ! in
      real*8 rmse_short                             ! in
      real*8 eshort                                 ! internal
      real*8 error                                  ! internal
      real*8 errorsum                               ! internal
      real*8 abserror                               ! internal
      real*8 nnatomenergy(max_num_atoms)            ! internal
      real*8 weights_short(num_weightsshort,nelem)                          ! in/out
      real*8 weights(num_weightsshort)                                      ! internal 
      real*8 symfunction_list(num_funcvalues,max_num_atoms,nblock)          ! in
      real*8, dimension(:,:)  , allocatable :: symfunction
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
      real*8 time1                                              
      real*8 time2
      real*8 klambda                                                        ! internal
!!
      character*1 actfunc_short(num_layersshort)                            ! in
!!
      logical ldebug                                                        ! in
!!
!! initialization
       deshortdw(:,:,:)    = 0.0d0
       deshortdwsum(:,:,:) = 0.0d0
       eshort              = 0.0d0
       day                 = 0
       timestart           = 0.0d0
       timeend             = 0.0d0
       nenergy             = 0
       errorsum            = 0.0d0
       num_atoms_element(:)= 0
       jdx(:)              = 0
       kdx(:)              = 0
       mdx(:)              = 0
       nbad                = 0
!!
      call abstime(timestart,day)
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!! debug
!!      if(mpirank.eq.0)then
!!      write(ounit,*)' number of processes ',mpisize
!!      write(ounit,*)mpirank,' initial points ',npoints
!!      endif
!!
!! decide which points to use if we have a random subset of points
!!    icount counts the points that will be used
!!    jdx() stores the indices of the random points to be used
!!    the array jdx() is smaller than idx() and contains only the relevant npts points
!!
      if(mpirank.eq.0)then 
        npts=0
        do i1=1,npoints 
          z=ran0(lseed)
          if(z.lt.energyrnd)then
            npts=npts+1
            jdx(npts)=idx(i1)
          else
!!            write(ounit,*)'rejected ',i1
          endif
        enddo ! i1
!!
!! debug
!!        write(ounit,*)npts,' points will be used'
!!        do i1=1,npts
!!          write(ounit,*)'point ',i1,jdx(i1)
!!        enddo
      endif ! mpirank.eq.0
!!
!! distribute data to all processes
      call mpi_bcast(npts,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(jdx,nblock,mpi_integer,0,mpi_comm_world,mpierror)
!!
!! debug
!!      if(mpirank.eq.0)then
!!        write(ounit,*)mpirank,' randomly selected points ',npts
!!      endif
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! loop over groups of points
      ncount                 = npts
      jcount                 = 0
!!
 10   continue
      num_points=int(nenergygroup/mpisize)
      if(ncount.gt.nenergygroup)then
        icount=nenergygroup          ! number of structures in this cycle
        ncount=ncount-nenergygroup
      else
        icount=ncount
        ncount=ncount-icount
      endif
!!
!! debug
!!      write(ounit,*)'cycle starts with structures ',icount
!!
!! determine which structures should be calculated by this process
      call mpifitdistribution(ounit,icount,&
              nstruct,n_start,n_end)
!!
!! debug
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' calculates structures ',nstruct,n_start+jcount,n_end+jcount
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!
      errorsum=0.0d0
!!
      do i1=n_start+jcount,n_end+jcount
!!
        num_atoms=num_atoms_list(jdx(i1))
        zelem(:)=zelem_list(jdx(i1),:)
!!
        allocate(symfunction(num_funcvalues,num_atoms))
!!
!! copy symmetry functions for atoms n_start to n_end to mpi array
        symfunction(:,:) = symfunction_list(:,:,jdx(i1))
!!
        nnatomenergy(:)= 0.0d0
        eshort         = 0.0d0
        call calconeshort(ounit,max_num_atoms,num_atoms,num_funcvalues,&
          nelem,num_weightsshort,maxnodes_short,&
          num_layersshort,nodes_short,elementindex,zelem,&
          windex,symfunction,weights_short,&
          eshort,nnatomenergy,actfunc_short,ldebug)
!!
        eshort=eshort/dble(num_atoms)
!! 
!! debug
!!        write(ounit,*)mpirank,i1,eshort
!!
        deallocate(symfunction)
!!
!! Caution: Using eshort-shortenergy_list(i1) instead requires sign change of derivatives dedw
        abserror=abs(shortenergy_list(jdx(i1))-eshort)
        error=shortenergy_list(jdx(i1))-eshort
!!       
        if(abserror.gt.kalmanthreshold*rmse_short)then
!!          write(ounit,*)mpirank,i1,'bad structure found',abserror
          nbad=nbad+1
          kdx(i1)=jdx(i1)
          errorsum=errorsum+error
        endif
      enddo ! i1
!!
!! debug
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' found bad points ',nbad
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
      nbad_sum=nbad
      call mpi_allreduce(nbad_sum,nbad_sum,&
       1,mpi_integer,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(kdx,kdx,&
       nblock,mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!! debug
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!      if(mpirank.eq.0)then
!!        write(ounit,*)mpirank,' bad points ',nbad_sum
!!      endif
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! decide if we have collected enough bad structures to justify update
      jcount=jcount+icount
      if(nbad_sum.ge.nenergygroup)then
!! enough bad structures found: just go ahead
      elseif(ncount.gt.0)then
!! get more bad structures
!!        call mpi_barrier(mpi_comm_world,mpierror)
!!        if(mpirank.eq.0)then
!!          write(ounit,*)mpirank,' needs more bad structures ',nbad_sum
!!        endif
        goto 10
      else
!! not enough bad structures, but no points left: just go ahead
      endif
!!
!! debug
!!      if(mpirank.eq.0)then
!!        write(ounit,*)mpirank,' decided to update with ',nbad_sum
!!      endif
!!
      if(mpirank.eq.0)then
!!        write(ounit,*)'number of bad structures ',nbad_sum
        kcount=0
        do i1=1,npoints
!!          write(ounit,*)'kdx ',i1,kdx(i1)
          if(kdx(i1).ne.0)then
            kcount=kcount+1
            mdx(kcount)=kdx(i1)
          endif
        enddo
      endif ! mpirank.eq.0
      call mpi_bcast(mdx,nblock,mpi_integer,0,mpi_comm_world,mpierror)
!!
!! prepare the calculation of the gradients
!! determine which structures should be calculated by this process
      call mpifitdistribution(ounit,nbad_sum,&
              nstruct,n_start,n_end)
!!
!! debug
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!      write(ounit,*)mpirank,' calculates bad structures ',nstruct,n_start,n_end
!!      if(mpirank.eq.0)then
!!        do i1=1,nbad_sum
!!          write(ounit,*)i1,'using mdx ',mdx(i1)
!!        enddo
!!      endif ! mpirank.eq.0
!!      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! initialization for this update
      deshortdwsum(:,:,:)=0.0d0
!!
      num_atoms_element(:)=0
!!
      do i1=n_start,n_end

        num_atoms=num_atoms_list(mdx(i1))
        zelem(:)=zelem_list(mdx(i1),:)
!!
        allocate(symfunction(num_funcvalues,num_atoms))
!! copy symmetry functions for atoms n_start to n_end to mpi array
        symfunction(:,:) = symfunction_list(:,:,mdx(i1))
!!
!! calculate the number of atoms per element, because we don't want to
!! update weights for elements not being present in current example
            do i2=1,num_atoms
              num_atoms_element(elementindex(zelem(i2)))=&
                num_atoms_element(elementindex(zelem(i2)))+1
            enddo ! i2
!!
        call getdeshortdw(ounit,nelem,num_weightsshort,&
           zelem,elementindex,max_num_atoms,num_atoms,&
           maxnodes_short,&
           windex,num_layersshort,num_funcvalues,nodes_short,&
           symfunction,&
           weights_short,deshortdw,&
           actfunc_short,ldebug)
!!
!! sum derivatives and errors if we group over several points
        deshortdwsum(:,:,:)=deshortdwsum(:,:,:)+deshortdw(:,:,:)


        deallocate(symfunction)
!!
      enddo
!!
      call mpi_allreduce(errorsum,errorsum,&
       1,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(deshortdwsum,deshortdwsum,&
       nelem*num_weightsshort*nodes_short(num_layersshort),&
       mpi_real8,mpi_sum,mpi_comm_world,mpierror)
      call mpi_allreduce(num_atoms_element,num_atoms_element,&
       nelem,mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!!      if(mpirank.eq.0)then
!!        write(ounit,*)mpirank,' errorsum,nbad_sum ',errorsum,nbad_sum
!!      endif

!! normalize the error and the derivatives in case of grouping
!!      if(mpirank.eq.0)then

       error=errorsum/dble(nbad_sum)
       deshortdw(:,:,:)=deshortdwsum(:,:,:)/dble(nbad_sum)

!!
!! chose optimization algorithm
                if(optmode.eq.1)then ! Kalman filter
!! update weights for each element
                  do i2=1,nelem
!! check if atoms of the element are present, only then do update
                    if(num_atoms_element(i2).gt.0)then
!!                      corrmatrix(:)= corrmatrix_list(:,i2)
!!                      kalgainmat(:)= kalgainmat_list(:,i2)
!!                      weights(:)   = weights_short(:,i2)
!!                      dedw(:)      = deshortdw(:,1,i2)
!!                      klambda      = kalmanlambda(i2)
!!
                      call abstime(time1,day)
                      if(mpisize.eq.1)then
                        call updatekalman(ounit,kaldim,corrdim,&
                          num_weightsshort,kalmanlambda(i2),kalmannue,&
                          weights_short(1,i2),deshortdw(1,1,i2),&
                          kalgainmat_list(1,i2),corrmatrix_list(1,i2),&
                          error,ldebug)
                      else ! real parallel case
                        call updatekalman_para(ounit,kaldim,corrdim,&
                          num_weightsshort,kalmanlambda(i2),kalmannue,&
                          weights_short(1,i2),deshortdw(1,1,i2),&
                          kalgainmat_list(1,i2),corrmatrix_list(1,i2),&
                          error,ldebug)
                      endif ! mpisize.eq.1
                      call abstime(time2,day)
                      write(ounit,'(a,f10.3)')'time updatekalman_para ',time2-time1
!!
!!                      corrmatrix_list(:,i2)= corrmatrix(:)
!!                      kalgainmat_list(:,i2)= kalgainmat(:)
!!                      weights_short(:,i2)  = weights(:)
!!                      kalmanlambda(i2)     = klambda
                    endif ! num_atoms_element(i2).gt.0
                  enddo ! i2=1,nelem
!!
                elseif(optmode.eq.2)then ! conjugate gradient'
                  write(ounit,*)'CG optimization not yet implemented in optimize_short'
                  stop
!! '
                elseif(optmode.eq.3)then ! steepest descent'
                  write(ounit,*)'SD optimization not yet implemented in optimize_short'
                  stop
                endif ! optmode'

!!      endif ! mpirank.eq.0

!! if only process 0 is doing the weight updata, the weights have to be distributed
!!      call mpi_bcast(weights_short,nelem*num_weightsshort,mpi_real8,0,mpi_comm_world,mpierror)



!! continue with the next update
      if(ncount.gt.0) then 
        kdx(:)=0.0d0
        nbad=0
        nbad_sum=0
        mdx(:)=0.0d0
        error=0.0d0
        errorsum=0.0d0
!!
        if(mpirank.eq.0)then
!!        write(ounit,*)mpirank,' jumping to next block of points'
        endif
        goto 10
      endif


      call abstime(timeend,day)
      write(ounit,'(i5,a,f10.3)')mpirank,' time optimize_short ',timeend-timestart
!!
      return
      end
