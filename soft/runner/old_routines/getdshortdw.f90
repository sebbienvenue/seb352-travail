!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting_batch.f90
!!
!! This routine is called once for each block of points
!!
      subroutine getdshortdw(npoints,point,&
         ndone,num_weights_short_atomic_free,&
         wconstraintidx,&
         numbere,numberf,nshort,&
         lseed,ntrain,fitstat,fitstatf,&
         kseed,errorshort,&
         rmse_short,rmse_force_s_ref,&
         minvalue_short_atomic,maxvalue_short_atomic,dshortdw)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_atomic
      use structures
!!
      implicit none
!!
      integer npoints                                  ! in
      integer ndone                                    ! in
      integer num_weights_short_atomic_free(nelem)              ! in
      integer wconstraintidx(maxnum_weights_short_atomic,nelem) ! in
      integer num_atoms_element_e(nelem)               ! internal
      integer zelem(max_num_atoms)                     ! internal
      integer num_atoms                                ! internal
      integer point                                    ! in/out
      integer kseed                                    ! in/out
      integer lseed                                    ! in/out
      integer i1,i2,i3,i4                              ! internal
      integer day                                      ! internal
      integer natoms                                   ! internal
      integer n_start                                  ! internal
      integer n_end                                    ! internal
      integer, dimension(:), allocatable :: atomindex  ! internal
      integer icount                                   ! internal
      integer numbere                                  ! in/out
      integer numberf                                  ! in/out
      integer ntrain                                   ! in
      integer fitstat(ntrain)                          ! in/out
      integer fitstatf(3,max_num_atoms,ntrain)         ! in/out
      integer nshort(nelem)                            ! in/out
      integer, allocatable :: lsta(:,:)                          ! numbers of neighbors
      integer, allocatable :: lstc(:)                            ! identification of atom
      integer, allocatable :: lste(:)                            ! nuclear charge of atom
      integer, allocatable :: num_neighbors_short_atomic(:)              
      integer max_num_neighbors_short_atomic
      integer, allocatable :: neighboridx_short_atomic(:,:)          
      integer, allocatable :: invneighboridx_short_atomic(:,:)  
!!
      real*8 errorshort(nelem)                         ! in/out
      real*8 rmse_short                                ! in
      real*8 eshort                                    ! internal
      real*8 errore                                    ! internal
      real*8 abserrore                                 ! internal
      real*8 errorf                                    ! internal
      real*8 abserrorf                                 ! internal
      real*8 nnatomenergy(max_num_atoms)               ! internal
      real*8 dshortdw(maxnum_weights_short_atomic,nelem)                    ! in/out
      real*8, dimension(:,:)  , allocatable :: symfunctiondummy             ! internal
      real*8 deshortdw(maxnum_weights_short_atomic,1,nelem)                 ! internal
      real*8 z,ran0                                                         ! internal
      real*8 nntotalenergy                                                  ! internal
      real*8 rmse_force_s_ref                                               ! in
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)                 ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)                 ! in
      real*8 shortforce(3,max_num_atoms)                                    ! internal
      real*8 nnshortforce(3,max_num_atoms)                                  ! internal
      real*8 dfshortdw(maxnum_weights_short_atomic)                         ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyz               ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyz_mpi           ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: strs                     ! internal dummy
      real*8 errorthreshold                                                 ! internal 
      real*8 errorthresholdf                                                ! internal 
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
!!
      logical lperiodic                                                     ! internal
      logical lrmin                                                         ! in
!!
!! initializations
      day                    = 0
      eshort                 = 0.0d0
      nntotalenergy          = 0.0d0
      shortforce(:,:)        = 0.0d0
      nnshortforce(:,:)      = 0.0d0
      deshortdw(:,:,:)       = 0.0d0
      dfshortdw(:)           = 0.0d0
      errore                 = 0.0d0
      errorf                 = 0.0d0
!!
!! loop over all points
      do i1=1,npoints
        point=point+1
!!
!! copy structure-specific information on local arrays 
        num_atoms        = num_atoms_list(i1)
        lperiodic        = lperiodic_list(i1)
        shortforce(:,:)  = shortforce_list(:,:,i1)
        zelem(:)         = zelem_list(i1,:)
!!
!! calculate the number of atoms per element
        num_atoms_element_e(:) = 0
        do i2=1,num_atoms
          num_atoms_element_e(elementindex(zelem(i2)))=&
            num_atoms_element_e(elementindex(zelem(i2)))+1
        enddo ! i2
!!
!! determine which atoms of this structure should be calculated by this process
        call mpifitdistribution(num_atoms,natoms,n_start,n_end)
!!
!! determine the atomindex array
        allocate(atomindex(natoms))
        do i2=1,natoms
          atomindex(i2)=n_start+i2-1
        enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! energy part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! decide if we use this point
        z=ran0(lseed)
        if(z.gt.energyrnd) goto 98 ! skip this point
        if(shortenergy_list(i1).gt.maxenergy) goto 98 ! skip this point
!!
!! set errorthreshold criterion for short range energy
        if(lfixederrore)then
          errorthreshold=fixederrore
        else
          errorthreshold=kalmanthreshold*rmse_short
        endif
!!
!! predict the short-range atomic energies for atoms n_start to n_end
        nnatomenergy(:)= 0.0d0
        nntotalenergy  = 0.0d0
        call calconeshort_para(point,natoms,atomindex,&
          zelem,&
          symfunction_short_atomic_list(1,n_start,i1),nnatomenergy,&
          nntotalenergy)
!! merge all atomic energies 
        call mpi_allreduce(mpi_in_place,nnatomenergy,&
          max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
        nntotalenergy=0.0d0
        do i2=1,num_atoms
          nntotalenergy=nntotalenergy+nnatomenergy(i2)
        enddo
!! calculate short range energy eshort and normalize per atom
        eshort=nntotalenergy/dble(num_atoms)
!!
!! calculate the error (per atom) of this training point
        abserrore=abs(shortenergy_list(i1)-eshort)
        errore=shortenergy_list(i1)-eshort
!!
!! decide if point is sufficiently bad to be included in weight update
        if(abserrore.gt.errorthreshold)then
!!
!! calculate the derivative of the short range energy with respect to the weights
          if(paramode.eq.1)then
!! serial version
            call getdeshortdw(&
              num_weights_short_atomic_free,&
              zelem,num_atoms,&
              wconstraintidx,&
              symfunction_short_atomic_list(1,1,i1),deshortdw)
!!
          elseif(paramode.eq.2)then
!! parallel version
!! CAUTION: deshortdw depends slightly on the number of processes in parallel runs
            call getdeshortdw_para(&
              num_atoms,natoms,atomindex,&
              num_weights_short_atomic_free,&
              zelem,wconstraintidx,&
              symfunction_short_atomic_list(1,n_start,i1),deshortdw)
          else 
            write(ounit,*)'Error: paramode not allowed in optimize_short_combined'
            stop !'
          endif ! paramode
!!
!! sum derivatives and errors if we group over several points
          do i2=1,nelem
            dshortdw(:,i2)    = dshortdw(:,i2)+deshortdw(:,1,i2)
            nshort(i2)        = nshort(i2)+1
          enddo ! i2
          numbere          = numbere+1
          fitstat(ndone+i1)= fitstat(ndone+i1)+1
          do i2=1,nelem
            errorshort(i2) = errorshort(i2)&
              +abserrore*dble(num_atoms_element_e(i2))/dble(num_atoms)
          enddo
!!
        endif 
!!
!! jump here if energy update of this point is skipped
 98     continue
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! force part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(luseforces.and.(forcernd.gt.0.0d0))then
!!
!! allocations for parallel use
          allocate(strs(3,3,maxnum_funcvalues_short_atomic,natoms))
          allocate(dsfuncdxyz_mpi(maxnum_funcvalues_short_atomic,natoms,0:max_num_neighbors_short_atomic,3))
          dsfuncdxyz_mpi(:,:,:,:)=0.0d0
          allocate(symfunctiondummy(maxnum_funcvalues_short_atomic,natoms))
!!
          allocate(lsta(2,max_num_atoms))
          allocate(lstc(listdim))
          allocate(lste(listdim))
          allocate(lstb(listdim,4))
          allocate(num_neighbors_short_atomic(natoms))
!!
          call getneighborsatomic(n_start,n_end,natoms,&
            num_atoms,num_neighbors_short_atomic,zelem,max_num_neighbors_short_atomic,&
            lsta,lstc,lste,&
            maxcutoff_short_atomic,lattice_list(1,1,i1),xyzstruct_list(1,1,i1),&
            lstb,lperiodic)
!!
          allocate(dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3))
          allocate(neighboridx_short_atomic(natoms,0:max_num_neighbors_short_atomic))  
          allocate(invneighboridx_short_atomic(natoms,max_num_atoms))  
          call getneighboridxatomic(n_start,n_end,natoms,listdim,&
            max_num_atoms,max_num_neighbors_short_atomic,&
            lsta,lstc,neighboridx_short_atomic,invneighboridx_short_atomic)
!!
!! get dsfuncdxyz:
!! this is independent of the weights and needs to be done only once for each point i1
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!! Caution: symmetry functions from this subroutine cannot be used because they are not scaled
!! parallel version
!! calculate the short range symmetry functions for atoms n_start to n_end
          call calconefunction_atomic(cutoff_type,max_num_neighbors_short_atomic,&
            max_num_atoms,n_start,n_end,natoms,natoms,elementindex,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic, &
            nelem,zelem,listdim,&
            lsta,lstc,lste,invneighboridx_short_atomic,&
            function_type_short_atomic,symelement_short_atomic,&
            xyzstruct_list(1,1,i1),symfunctiondummy,0.0d0,&
            funccutoff_short_atomic,eta_short_atomic,rshift_short_atomic,&
            lambda_short_atomic,zeta_short_atomic,dsfuncdxyz_mpi,strs,lstb,&
            lperiodic,.true.,.false.,lrmin)
          if(.not.lrmin)then
            write(ounit,*)'Error in getdshortdw: lrmin=.false.'
            stop
          endif
!!
!! scale dsfuncdxyz
          if(lscalesym)then
            call scaledsfunc_para(natoms,atomindex,max_num_neighbors_short_atomic,&
              maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
              nelem,minvalue_short_atomic,maxvalue_short_atomic,&
              scmin_short_atomic,scmax_short_atomic,&
              zelem,dsfuncdxyz_mpi,strs)
          endif ! lscalesym
!!
!! combine dsfuncdxyz of all processes
          dsfuncdxyz(:,:,:,:)=0.0d0
          icount=0
          do i2=n_start,n_end
            icount=icount+1
            dsfuncdxyz(:,i2,:,:)=dsfuncdxyz_mpi(:,icount,:,:)
          enddo
          call mpi_allreduce(mpi_in_place,dsfuncdxyz,&
            maxnum_funcvalues_short_atomic*max_num_atoms*max_num_neighbors_short_atomic*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!! 
!! calculate the error of the forces individually for each atom and component
          do i2=1,num_atoms
            do i3=1,3 ! loop over x,y and z
!!
!! if only a random subset of points is used: decide if we use this point
              z=ran0(kseed)
              if(z.gt.forcernd) goto 99 ! don't use force for weight update
              if(abs(shortforce(i3,i2)).gt.maxforce) goto 99 ! don't use force for weight update
!!
!! calculate the NN forces nnshortforce for current weight values
              nnshortforce(:,:)=0.0d0
!!
!! this subroutine calculates only the one force that is needed (specified by i2,i3), not the full array nnshortforce
!! parallel version
!! This scales nearly perfectly! 
              call getoneshortforce_para(i2,i3,max_num_neighbors_short_atomic,&
                num_neighbors_short_atomic,neighboridx_short_atomic,&
                invneighboridx_short_atomic,natoms,atomindex,&
                num_atoms,zelem,&
                symfunction_short_atomic_list(1,1,i1),&
                dsfuncdxyz,nnshortforce)
!!
!! calculate the error for the current force component
              errorf=(shortforce(i3,i2)-nnshortforce(i3,i2))/dble(num_atoms)  ! LABEL1
              abserrorf=abs(errorf)
              errorf=-1.d0*errorf          ! this must be used!
!!
!! set errorthreshold criterion for short range energy
              if(lfixederrorf)then
                errorthresholdf=fixederrorf
              else
                errorthresholdf=kalmanthresholdf*rmse_force_s_ref
              endif
!!
!! decide if force should be used for update
              if(abserrorf.gt.(errorthresholdf/dble(num_atoms)))then ! LABEL1
!!
!! calculate the derivative of the short range forces with respect to the weights
!! dfshortdw is for one atom only, but depends on all atoms
!! calculate dfshortdw
!! serial version
                if(mpirank.eq.0)then
                  call getdfshortdw(i2,i3,natoms,max_num_neighbors_short_atomic,num_neighbors_short_atomic,&
                    neighboridx_short_atomic,invneighboridx_short_atomic,&
                    num_weights_short_atomic_free,zelem,num_atoms,&
                    wconstraintidx,symfunction_short_atomic_list(1,1,i1),dsfuncdxyz,&
                    dfshortdw)
                endif ! mpirank.eq.0
!! CHECK: Do we have to normalize dfshortdw per atom here to be consistent with atomic energies?
                dfshortdw(:)=dfshortdw(:)/dble(num_atoms)
!!
                call mpi_bcast(dfshortdw,maxnum_weights_short_atomic,&
                  mpi_real8,0,mpi_comm_world,mpierror)
!!
                do i4=1,num_weights_short_atomic(elementindex(zelem(i2)))
                  dshortdw(i4,elementindex(zelem(i2)))=&
                    dshortdw(i4,elementindex(zelem(i2)))+dfshortdw(i4)
                enddo ! i4
                errorshort(elementindex(zelem(i2)))=&
                  errorshort(elementindex(zelem(i2)))+abserrorf
                fitstatf(i3,i2,ndone+i1)=fitstatf(i3,i2,ndone+i1)+1
                numberf = numberf+1
                nshort(elementindex(zelem(i2))) = nshort(elementindex(zelem(i2)))+1
!!
!! reinitializations for next update
                dfshortdw(:)          = 0.0d0
                errorf                = 0.0d0
              endif !(abserrorf.gt.kalmanthresholdf*rmse_force_s_ref)
!!
!! jump here if force should be skipped for weight update
 99           continue
!!
            enddo ! i3 loop over x,y and z
          enddo ! i2 loop over atoms
!!
          deallocate(strs)
          deallocate(dsfuncdxyz_mpi)
          deallocate(dsfuncdxyz)
          deallocate(symfunctiondummy)
          deallocate(lsta)
          deallocate(lstc)
          deallocate(lste)
          deallocate(lstb)
          deallocate(neighboridx_short_atomic)  
          deallocate(invneighboridx_short_atomic)  
          deallocate(num_neighbors_short_atomic)
!!
        endif ! luseforces
!!
 101    continue
        deallocate(atomindex)      
!!
      enddo ! i1=1,npoints
!!
!!
      return
      end
