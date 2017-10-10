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
      subroutine optimize_short_para(ounit,nblock,npoints,point,idx,&
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
      integer lseed                                 ! in/out
      integer i1,i2,i3,i4                           ! internal
      integer day                                   ! internal
      integer nenergygroup                          ! in
      integer nenergy                               ! internal
      integer natoms                                ! internal
      integer n_start                               ! internal
      integer n_end                                 ! internal
      integer, dimension(:), allocatable :: atomindex ! internal
      integer icount                                ! internal
      integer nweights                              ! internal
      integer paramode                              ! in
!!
      real*8 kalmanthreshold                        ! in
      real*8 kalmanlambda(nelem)                    ! in
      real*8 kalmannue                              ! in
      real*8 rmse_short                             ! in
      real*8 eshort                                 ! internal
      real*8 errore                                  ! internal
      real*8 erroresum                               ! internal
      real*8 abserrore                               ! internal
      real*8 nnatomenergy(max_num_atoms)            ! internal
      real*8 weights_short(num_weightsshort,nelem)                          ! in/out
      real*8 weights(num_weightsshort)                                      ! internal 
      real*8 symfunction_list(num_funcvalues,max_num_atoms,nblock)          ! in
      real*8 shortenergy_list(nblock)                                       ! in
      real*8 deshortdw(num_weightsshort,nodes_short(num_layersshort),nelem) ! internal
      real*8 deshortdwsum(num_weightsshort,nodes_short(num_layersshort),nelem) ! internal
      real*8 dedw(num_weightsshort)                                         ! internal
      real*8 corrmatrix_list(corrdim,nelem)                                 ! in/out
      real*8 kalgainmat_list(kaldim,nelem)                                  ! in/out
      real*8 energyrnd                                                      ! in
      real*8 z,ran0                                                         ! internal
      real*8 timestart
      real*8 timeend
      real*8 nntotalenergy                                                  ! internal
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
       nenergy             = 0
       erroresum            = 0.0d0
       num_atoms_element(:)= 0
       nntotalenergy       = 0.0d0
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!!
!! for paramode.eq.1 distribute each single point to multiple processes
!! loop over all points
        do i1=1,npoints
          point=point+1
!!
!! if only a random subset of points is used: decide if we use this point
          z=ran0(lseed)
          if(z.lt.energyrnd)then
!!
!! get the number of atoms of this structure
            num_atoms=num_atoms_list(idx(i1))
!!
!! determine which atoms of this structure should be calculated by this process
            call mpifitdistribution(ounit,num_atoms,&
              natoms,n_start,n_end)
!!
!! determine the atomindex array
            allocate(atomindex(natoms))
            do i2=1,natoms
              atomindex(i2)=n_start+i2-1
            enddo
!!
!!
            zelem(:)=zelem_list(idx(i1),:)
!!
!! predict the short-range atomic energies for atoms n_start to n_end
!!            call abstime(timestart,day)
            nnatomenergy(:)= 0.0d0
            nntotalenergy  = 0.0d0
            call calconeshort_para(ounit,natoms,atomindex,&
              max_num_atoms,num_atoms,num_funcvalues,&
              nelem,num_weightsshort,maxnodes_short,&
              num_layersshort,nodes_short,elementindex,zelem,&
              windex,symfunction_list(1,n_start,idx(i1)),weights_short,nnatomenergy,&
              nntotalenergy,actfunc_short,ldebug)
!!
!! merge all atomic energies 
            call mpi_allreduce(nnatomenergy,nnatomenergy,&
              max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
            call mpi_allreduce(nntotalenergy,nntotalenergy,&
              1,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! calculate eshort and normalize per atom
            eshort=nntotalenergy/dble(num_atoms)
!!            call abstime(timeend,day)
!!            if(mpirank.eq.0)then
!!              write(ounit,'(a,i6,f10.3)')'time calconeshort_para: ',i1,timeend-timestart
!!            endif
!!
!! Caution: Using eshort-shortenergy_list(i1) instead requires sign change of derivatives dedw
            abserrore=abs(shortenergy_list(idx(i1))-eshort)
            errore=shortenergy_list(idx(i1))-eshort
!!
!! decide if point is sufficiently bad to justify update
           if(abserrore.gt.kalmanthreshold*rmse_short)then
             nenergy=nenergy+1
!!
!! calculate the derivative of the short range energy with respect to the weights
!!            call abstime(timestart,day)
            call getdeshortdw_para(ounit,natoms,atomindex,&
              nelem,num_weightsshort,&
              zelem,elementindex,max_num_atoms,num_atoms,&
              maxnodes_short,&
              windex,num_layersshort,num_funcvalues,nodes_short,&
              symfunction_list(1,n_start,idx(i1)),weights_short,deshortdw,&
              actfunc_short,ldebug)
!!
!! sum derivatives and errors if we group over several points
            deshortdwsum(:,:,:)=deshortdwsum(:,:,:)+deshortdw(:,:,:)
            erroresum=erroresum+errore
!!            call abstime(timeend,day)
!!            if(mpirank.eq.0)then
!!              write(ounit,'(a,i6,f10.3)')'time getdeshortdw: ',i1,timeend-timestart
!!            endif
!!
!! calculate the number of atoms per element, because we don't want to 
!! update weights for elements not being present in current example
            do i2=1,num_atoms
              num_atoms_element(elementindex(zelem(i2)))=&
                num_atoms_element(elementindex(zelem(i2)))+1
            enddo ! i2
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if we have no parallel updatekalman, we should do the matrix operations only with one process
!!            if(mpirank.eq.0)then
!!
!! check if an update should be done
              if((mod(nenergy,nenergygroup).eq.0).or.&
                (i1.eq.npoints))then
!!
!! normalize the error and the derivatives in case of grouping
                errore=erroresum/dble(nenergy)
                deshortdw(:,:,:)=deshortdwsum(:,:,:)/dble(nenergy)
!!
!! chose optimization algorithm
                if(optmode.eq.1)then ! Kalman filter
!! update weights for each element
                  do i2=1,nelem
!! check if atoms of the element are present, only then do update
                    if(num_atoms_element(i2).gt.0)then 
!!
                      if(mpisize.eq.1)then
                        call updatekalman(ounit,kaldim,corrdim,&
!! updatekalman_para cannot be used for 1 process because the corrmatrix array has a different dimension then
!!                        call updatekalman_para(ounit,kaldim,corrdim,&
                          num_weightsshort,kalmanlambda(i2),kalmannue,&
                          weights_short(1,i2),deshortdw(1,1,i2),&
                          kalgainmat_list(1,i2),corrmatrix_list(1,i2),&
                          errore,ldebug)
                      else ! real parallel case
                        call updatekalman_para(ounit,kaldim,corrdim,&
!! updatekalman cannot be used for more than 1 process because the corrmatrix array has a different dimension then
!!                        call updatekalman(ounit,kaldim,corrdim,&
                          num_weightsshort,kalmanlambda(i2),kalmannue,&
                          weights_short(1,i2),deshortdw(1,1,i2),&
                          kalgainmat_list(1,i2),corrmatrix_list(1,i2),&
                          errore,ldebug)
                      endif ! mpisize.eq.1
!!
                    endif ! num_atoms_element(i2).gt.0
                  enddo ! i2=1,nelem
!!
                elseif(optmode.eq.2)then ! conjugate gradient'
                  write(ounit,*)'CG optimization not yet implemented in optimize_short'
                  stop
!!
!!                   call lbfgs(num_weightsshort,mbfgs,&
!!                     weights_short(1,i2),&
!!                     eshort,-deshortdw(1,1,i2),.false.,diag,&
!!                     -1,1.d-4,1.d-16,work,iflag)
!! '
                elseif(optmode.eq.3)then ! steepest descent'
!! update weights for each element
                  do i2=1,nelem
!! check if atoms of the element are present, only then do update
                    if(num_atoms_element(i2).gt.0)then ! atoms of element are present
                      if(mpisize.gt.1)then
!!
!! determine which atoms of this structure should be calculated by this process
                        call mpifitdistribution(ounit,num_weightsshort,&
                          nweights,n_start,n_end)
!!
                        write(ounit,*)'Warning: steepest descent cannot be parallel if '
                        write(ounit,*)'only mpirank 0 does the update'
                        stop
!! copy arrays'
                        weights(:)   =0.0d0
                        do i3=n_start,n_end
                          weights(i3)   =weights_short(i3,i2)
                          dedw(i3)      =deshortdw(i3,1,i2)
                        enddo
!!
                        call updatesteepest_para(ounit,nelem,n_start,n_end,&
                          num_weightsshort,weights,dedw,errore,ldebug)
!!
!! combine all weights arrays
                        call mpi_allreduce(weights,weights,&
                          num_weightsshort,&
                          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! copy weight array back
                        weights_short(:,i2)=weights(:)
!!
                      else ! mpisize .eq. 1
!! serial original
                        weights(:)   =weights_short(:,i2)
                        dedw(:)      =deshortdw(:,1,i2)
!!
                        call updatesteepest(ounit,&
                          num_weightsshort,weights,dedw,errore,ldebug)
!!
                        weights_short(:,i2)=weights(:)
                      endif ! mpisize
                    endif ! num_atoms_element(i2).gt.0
                  enddo ! i2=1,nelem
                endif ! optmode
!!
!! reinitializations for next update
                nenergy              = 0
                deshortdwsum(:,:,:)  = 0.0d0
                erroresum             = 0.0d0
                num_atoms_element(:) = 0
!!
              endif !((mod(nenergy,nenergygroup).eq.0).or.(i1.eq.npoints))then
!!
!!            endif ! mpirank.eq.0
!!
!! if only process 0 is doing the weight updata, the weights have to be distributed
!!            call mpi_bcast(weights_short,nelem*num_weightsshort,mpi_real8,0,mpi_comm_world,mpierror)
!!
          else
!! point is not used for short update
          endif !(abserror.gt.kalmanthreshold*rmse_short)
!!
        deallocate(atomindex)      

        endif !(z.lt.energyrnd)then
!!
      enddo ! i1=1,npoints
!!
      return
      end
