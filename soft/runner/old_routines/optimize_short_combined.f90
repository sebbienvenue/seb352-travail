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
!! Note:
!! maybe the force fits get better by using the 2 lines labeled by LABEL1 ?
!! => this seems to be true!
!!
      subroutine optimize_short_combined(npoints,point,idx,&
         kaldim,maxcorrdim,maxcorrfdim,corrdim,corrfdim,countepoch,&
         ndone,num_weights_short_atomic_free,&
         wconstraintidx,numbere,numberf,&
         lseed,ntrain,fitstat,fitstatf,kseed,&
         kalmanthreshold_temp,kalmanthresholdf_temp,&
         rmse_short,rmse_force_s_ref,&
         corrmatrix_list,corrmatrixf_list,&
         minvalue_short_atomic,maxvalue_short_atomic)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_atomic
      use structures
      use timings
!!
      implicit none
!!
      integer npoints                                  ! in
      integer ndone                                    ! in
      integer countepoch                               ! in
      integer num_weights_short_atomic_free(nelem)     ! in
      integer totnum_weights_short_atomic_free         ! internal
      integer totnum_weights_short_atomic              ! internal
      integer wconstraintidx(maxnum_weights_short_atomic,nelem) ! in
      integer zelem(max_num_atoms)                     ! internal
      integer num_atoms                                ! internal
      integer num_atoms_element(nelem)                 ! internal 
      integer num_atoms_element_struct(nelem)          ! internal 
      integer num_atoms_element_block(nelem)           ! internal 
      integer count_atoms_element_f(nelem)             ! internal 
      integer count_atoms_element_block(nelem)         ! internal 
      integer point                                    ! in/out
      integer idx(nblock)                              ! in
      integer kseed                                    ! in/out
      integer lseed                                    ! in/out
      integer i1,i2,i3,i4,i5                           ! internal
      integer mforcegroup(nelem)                       ! internal
      integer nenergy                                  ! internal
      integer natoms                                   ! internal
      integer n_start                                  ! internal
      integer n_end                                    ! internal
      integer, dimension(:), allocatable :: atomindex  ! internal
      integer icount                                   ! internal
      integer numbere                                  ! in/out
      integer numberf                                  ! in/out
      integer ntrain                                   ! in
      integer, allocatable :: lsta(:,:)                          ! numbers of neighbors
      integer, allocatable :: lstc(:)                            ! identification of atom
      integer, allocatable :: lste(:)                            ! nuclear charge of atom
      integer, allocatable :: num_neighbors(:)                
      integer max_num_neighbors_atomic
      integer, allocatable :: neighboridx(:,:)   
      integer, allocatable :: invneighboridx(:,:)  
!!

!! Kalman matrix dimensions:
      integer maxcorrdim                               ! in
      integer maxcorrfdim                              ! in
      integer corrdim(nelem)                           ! in
      integer corrfdim(nelem)                          ! in
      integer kaldim(nelem)                            ! in
      integer fitstat(ntrain)                          ! in/out
      integer fitstatf(3,max_num_atoms,ntrain)         ! in/out
!!
      real*8 kalmanthreshold_temp                      ! in
      real*8 kalmanthresholdf_temp                     ! in
      real*8 rmse_short                                ! in
      real*8 nneshort                                  ! internal
      real*8 errore                                    ! internal
      real*8 erroresum                                 ! internal
      real*8 abserrore                                 ! internal
      real*8 errorf                                    ! internal
      real*8 errorfsum(nelem)                          ! internal
      real*8 abserrorf                                 ! internal
      real*8 errorjoint(nelem)                         ! internal
      real*8 nnatomenergy(max_num_atoms)               ! internal
      real*8, dimension(:)  , allocatable :: weights                        ! internal
      real*8, dimension(:,:)  , allocatable :: symfunctiondummy             ! internal
      real*8 nneshort_list(nblock)                                          ! internal
      real*8 deshortdw(maxnum_weights_short_atomic,1,nelem)                 ! internal
      real*8 deshortdwsum(maxnum_weights_short_atomic,1,nelem)              ! internal
      real*8, dimension(:) , allocatable :: deshortdw_temp                  ! internal
      real*8 corrmatrix_list(maxcorrdim,nelem)                              ! in/out
      real*8 corrmatrixf_list(maxcorrfdim,nelem)                            ! in/out
      real*8 z,ran0                                                         ! internal
      real*8 nntotalenergy                                                  ! internal
      real*8 nnshortforce_list(3,max_num_atoms,nblock)                      ! internal
      real*8 forceerror_list(3,max_num_atoms,nblock)                        ! internal
      real*8 rmse_force_s_ref                                               ! in
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)                 ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)                 ! in
      real*8 shortforce(3,max_num_atoms)                                    ! internal
      real*8 nnshortforce(3,max_num_atoms)                                  ! internal
      real*8 dfshortdw(maxnum_weights_short_atomic,nelem)                   ! internal
      real*8 dfshortdwsum(maxnum_weights_short_atomic,nelem)                ! internal
      real*8, dimension(:) , allocatable :: dfshortdw_temp                  ! internal
      real*8, dimension(:,:,:,:) , allocatable :: dsfuncdxyz                ! internal 
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyz_mpi           ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: strs                     ! internal dummy
      real*8 debugsum                                                       ! internal debugging
      real*8 debugsum2                                                      ! internal debugging
      real*8 scalefactorftemp(nelem)                                        ! internal
      real*8 energyerror_list(nblock)                                       ! internal
      real*8 sumwsquared                                                    ! internal
      real*8 errorthreshold                                                 ! internal 
      real*8 errorthresholdf                                                ! internal 
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
!!
      logical lperiodic                                                     ! internal
      logical lusee(nblock)                                                 ! internal 
      logical lusef(3,max_num_atoms,nblock)                                 ! internal 
      logical lrmin                                                         ! internal
      logical lenforceupdatee                                               ! internal
      logical lenforceupdatef(nelem)                                        ! internal
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!! initializations
      num_atoms_element(:)       = 0
      num_atoms_element_struct(:)= 0
      num_atoms_element_block(:) = 0
      count_atoms_element_f(:)   = 0
      count_atoms_element_block(:) = 0
      nneshort                   = 0.0d0
      nntotalenergy              = 0.0d0
      shortforce(:,:)            = 0.0d0
      nnshortforce(:,:)          = 0.0d0
      deshortdw(:,:,:)           = 0.0d0
      deshortdwsum(:,:,:)        = 0.0d0
      dfshortdw(:,:)             = 0.0d0
      dfshortdwsum(:,:)          = 0.0d0
      errore                     = 0.0d0
      erroresum                  = 0.0d0
      errorf                     = 0.0d0
      errorfsum(:)               = 0.0d0
      nenergy                    = 0
      lenforceupdatee            = .false.
      lenforceupdatef(:)         = .false.
!! timing variables
      dayefitting                = 0
      dayefittingrepeat          = 0
      dayffitting                = 0
      daydeshortdw               = 0
      daydeshortdwrepeat         = 0
      daydfshortdw               = 0
      dayeupdate                 = 0
      dayeupdaterepeat           = 0
      dayfupdate                 = 0
      dayferror                  = 0
      dayeerror                  = 0
      dayeerrorrepeat            = 0
      dayfsym                    = 0
!!
!! count total number of atoms of each element for this block of points
      do i1=1,npoints
        do i2=1,num_atoms_list(i1)
           num_atoms_element_block(elementindex(zelem_list(i1,i2)))=&
             num_atoms_element_block(elementindex(zelem_list(i1,i2)))+1
        enddo ! i2
      enddo ! i1
!!
!! determine mforcegroup for each element and make sure it is not larger than the number of forces for this element in this block of structures
      mforcegroup(:) = nforcegroup
      do i1=1,nelem
        if(mforcegroup(i1).gt.3*num_atoms_element_block(i1))then
          mforcegroup(i1)=3*num_atoms_element_block(i1)
          if(mforcegroup(i1).eq.0)then    ! in case element is not present
            mforcegroup(i1)=1
          endif
        endif
      enddo
!!
!! initialize total number of free short range weights
      totnum_weights_short_atomic_free = 0
      do i1=1,nelem
        totnum_weights_short_atomic_free=totnum_weights_short_atomic_free+num_weights_short_atomic_free(i1)
      enddo
!!
!! count total number of short range weights
      totnum_weights_short_atomic = 0
      do i1=1,nelem
        totnum_weights_short_atomic=totnum_weights_short_atomic+num_weights_short_atomic(i1)
      enddo
!!
!! determine the weighting factor between energy and force updates
!! CHECK/FIXME: should we put this later when a particular structure is read?
!! This is now done below in force part
!!      do i1=1,nelem
!!        if(scalefactorf.lt.0.0d0)then
!!          scalefactorftemp(i1)=(dble(ntrain*mforcegroup(i1))*energyrnd)/&
!!                    (dble(3*ntrainatoms*nenergygroup)*forcernd)
!!        else ! keep scalefactorf from input.nn
!!          scalefactorftemp(i1)=scalefactorf
!!        endif
!!      enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If requested, determine the points with the worst energies of this block of data
      if(luseworste)then
        if(worste.ge.1.0d0)then
          lusee(:)=.true. ! use all points
        else
!! determine the error of the short range energies
          call geteshort(nblock,npoints,&
            zelem_list,num_atoms_list,&
            symfunction_short_atomic_list,nneshort_list)
!! calculate the energy error of all points
          energyerror_list(:)=0.0d0
          do i1=1,npoints
            energyerror_list(i1)&
              =abs(nneshort_list(i1)-shortenergy_list(i1))
          enddo ! i1
!! sort points by error and determine error array lusee
          call sorteshorterror(npoints,energyerror_list,lusee)
        endif ! worste
      endif ! luseworste
!!
!! If requested, determine the points with the worst forces of this block of data
      if(luseworstf)then
        if(luseforces)then
          if(worstf.ge.1.0d0)then
            lusef(:,:,:)=.true. ! use all points
          else
!! determine the error of all short range forces
!! predict the short range NN forces for the test points here
            call getallshortforces(nblock,npoints,&
              num_atoms_list,zelem_list,&
              symfunction_short_atomic_list,nnshortforce_list,&
              lattice_list,xyzstruct_list,minvalue_short_atomic,maxvalue_short_atomic,&
              lperiodic_list)
!! calculate the error of all forces
!! Caution: in case of lupdatebyelement only the nnshortforces for that element are ok
            forceerror_list(:,:,:)=0.0d0
            do i1=1,npoints
              do i2=1,num_atoms_list(i1)
                do i3=1,3
                  forceerror_list(i3,i2,i1)&
                    =abs(nnshortforce_list(i3,i2,i1)-shortforce_list(i3,i2,i1))
                enddo ! i3
              enddo ! i2
            enddo ! i1
!! sort forces by error and determine array lusef
            call sortforceerror(npoints,forceerror_list,lusef)
          endif ! worstf
        endif ! luseforces
      endif ! luseworstf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! loop over all points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        do i1=1,npoints
          point=point+1
!!
          if((mpirank.eq.0).and.ldebug) then
            write(debugunit,'(a20,i8,a11,i8)')&
              'Short update point ',point,' structure ',idx(i1)+point-i1
          endif
!!
!! copy structure-specific information on local arrays 
          num_atoms        = num_atoms_list(idx(i1))
          lperiodic        = lperiodic_list(idx(i1))
          shortforce(:,:)  = shortforce_list(:,:,idx(i1))
          zelem(:)         = zelem_list(idx(i1),:)
!!
!! count the number of atoms per element, because we don't want to 
!! update weights for elements not being present 
          num_atoms_element_struct(:)=0
          do i2=1,num_atoms
            num_atoms_element(elementindex(zelem(i2)))=&
              num_atoms_element(elementindex(zelem(i2)))+1
            num_atoms_element_struct(elementindex(zelem(i2)))=&
              num_atoms_element_struct(elementindex(zelem(i2)))+1
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! energy part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          call abstime(timeefittingstart,dayefitting)
!! skip this energy if it is not among the worst energies and only these should be used
          if(luseworste.and.(.not.lusee(idx(i1))))then 
            goto 98 
          endif
!!
!! if only a random subset of points is used: decide if we use this point
          z=ran0(lseed)
          if(z.gt.energyrnd)then
            goto 98 
          endif
!!
!! skip this energy if it is too high 
          if(shortenergy_list(idx(i1)).gt.maxenergy)then
            goto 98 ! 
          endif
!!
!! predict the short-range atomic energies for atoms n_start to n_end
          call abstime(timeeerrorstart,dayeerror)
          nnatomenergy(:)= 0.0d0
          nntotalenergy  = 0.0d0
          call calconeshort_para(point,natoms,atomindex,&
            zelem,symfunction_short_atomic_list(1,n_start,idx(i1)),nnatomenergy,&
            nntotalenergy)
!!
!! combine nnatomenergy arrays from all processes 
          call mpi_allreduce(mpi_in_place,nnatomenergy,&
            max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! CAUTION: In parallel runs the result for nntotalenergy very slightly
!! depends on the number of processes because the order of the
!! summation of the nnatomenergies matters. Even invisible changes after the 15th digit
!! still sum up to significant errors in long fits
!! => don't use the nntotalenergy calculated in calconeshort_para but recalculate it in 
!! a well-defined way here 
          nntotalenergy=0.0d0
          do i2=1,num_atoms
            nntotalenergy=nntotalenergy+nnatomenergy(i2)
          enddo
!!
!! calculate short range energy eshort and normalize per atom
          nneshort=nntotalenergy/dble(num_atoms)
!!
!! calculate the error (per atom) of this training point
          errore=shortenergy_list(idx(i1))-nneshort
          abserrore=abs(errore)
          call abstime(timeeerrorend,dayeerror)
          timeeerror=timeeerror+(timeeerrorend-timeeerrorstart)
!!
!! skip energy update if energy error is within noise
          if(abserrore.lt.noisee)then
            goto 98 
          endif
!!
!! if weight damping is requested include normalized sum of squared weights in error 
          if(ldampw)then
            sumwsquared=0.0d0
            do i4=1,nelem
              do i3=1,num_weights_short_atomic(i4)
                sumwsquared=sumwsquared &
                  +(weights_short_atomic(i3,i4))**2.0d0
              enddo ! i3
            enddo ! i4
            errore=(1.d0-dampw)*errore + (dampw*sumwsquared)/dble(totnum_weights_short_atomic)
          endif ! ldampw
!!
!! set errorthreshold criterion for short range energy
          if(lfixederrore)then
            errorthreshold=fixederrore
          else
            errorthreshold=kalmanthreshold_temp*rmse_short
          endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! decide if energy error is sufficiently large to do weight update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          if((abserrore.gt.errorthreshold)&
            .or.&
            (i1.eq.npoints))then ! last point of this block reached, maybe update needs to be done
!!
            if(abserrore.le.errorthreshold) goto 77 ! don't use this gradient because point is good enough
!!
!! do some counting
            fitstat(ndone+idx(i1))=fitstat(ndone+idx(i1))+1
            nenergy=nenergy+1
            numbere=numbere+1
!!
!! calculate the derivative deshortdw of the short range energy with respect to the weights
            daydeshortdw      =0
            call abstime(timedeshortdwstart,daydeshortdw)
            if(paramode.eq.1)then
!! serial version
              call getdeshortdw(num_weights_short_atomic_free,&
                zelem,num_atoms,wconstraintidx,&
                symfunction_short_atomic_list(1,1,idx(i1)),deshortdw)
!!
            elseif(paramode.eq.2)then
!! parallel version
!! CAUTION: deshortdw depends slightly on the number of processes in parallel runs
              call getdeshortdw_para(&
                num_atoms,natoms,atomindex,&
                num_weights_short_atomic_free,zelem,&
                wconstraintidx,&
                symfunction_short_atomic_list(1,n_start,idx(i1)),deshortdw)
            else 
              write(ounit,*)'Error: paramode not allowed in optimize_short_combined'
              stop !'
            endif ! paramode
            call abstime(timedeshortdwend,daydeshortdw)
            timedeshortdw=timedeshortdw+(timedeshortdwend-timedeshortdwstart)
!!
!! debug: write deshortdw to debug.data               
            if((pstring(3:3).eq.'1').and.(mpirank.eq.0))then 
              do i3=1,nelem
                do i4=1,maxnum_weights_short_atomic
                  write(debugunit,'(2i6,a,2i6,f20.10)')&
                    countepoch,point,&
                    ' deshortdw ',i3,i4,deshortdw(i4,1,i3)
                enddo
              enddo
            endif
!!
!! debug: write length of deshortdw in two ways to debug.data 
            if(ldebug)then
              write(debugunit,'(i6,a,i6,x,f14.8)')countepoch,&
                ' errore of point ',point,errore
              do i3=1,nelem !'
                debugsum=0.0d0
                debugsum2=0.0d0
                do i4=1,num_weights_short_atomic(i3)
                  debugsum=debugsum+deshortdw(i4,1,i3)**2
                  debugsum2=debugsum2+deshortdw(i4,1,i3)
                enddo ! i4
                debugsum=dsqrt(debugsum)
                write(debugunit,'(i6,a,i6,x,a2,x,2f14.8)')&
                  countepoch,' deshortdw length point ',point,element(i3),debugsum,debugsum2
              enddo ! i3
            endif
!!
!! add the damping term to the derivatives if weight damping is used
!! CHECK this
            if(ldampw) then
              do i2=1,nelem
                do i3=1,num_weights_short_atomic_free(i2)
                  deshortdw(i3,1,i2)=(1.d0-dampw)*deshortdw(i3,1,i2) &
                    -dampw*2.0d0*weights_short_atomic(wconstraintidx(i3,i2),i2)/dble(totnum_weights_short_atomic_free)
                enddo ! i3
              enddo ! i2
            endif
!!
!! adapt sign, this is needed in case of nenergygroup.gt.1 to avoid cancellation effects upon averaging 
            if(errore.lt.0.0d0)then
              errore=-1.d0*errore
              deshortdw(:,:,:)=-1.d0*deshortdw(:,:,:)
            endif
!!
!! sum derivatives and errors if we group over several points
            deshortdwsum(:,:,:)= deshortdwsum(:,:,:)+deshortdw(:,:,:)
            erroresum          = erroresum+errore
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! energy update of the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! skip update in case of joint_energy_force_update 
            if(ljointefupdate) goto 98 
!!
 77         continue
!! check if there are energies left not used for update
            if((i1.eq.npoints).and.(nenergy.gt.0))then
              lenforceupdatee=.true.
            endif
!!
!! update if we have collected enough structures or there are no more structures left
            call abstime(timeeupdatestart,dayeupdate)
            if(((mod(nenergy,nenergygroup).eq.0)&   ! a sufficient number of energies has been accumulated for update
              .and.(nenergy.gt.0)).or.&             ! CAUTION: this is needed, otherwise mod is also 0 for nenergy=0
              (lenforceupdatee))then                ! enforce update for the remaining energies not yet used
!!
              lenforceupdatee = .false.
!!
!! normalize the error and the derivatives in case grouping has been used
              errore          = erroresum/dble(nenergy)
              deshortdw(:,:,:)= deshortdwsum(:,:,:)/dble(nenergy)
!!
!! Kalman damping
              errore=errore*kalman_dampe
!!
!! chose optimization algorithm
              if(optmodee.eq.1)then ! Kalman filter
!! update weights for each element
                do i2=1,nelem
!!
!! do not update for this element if it is not requested
                  if(lupdatebyelement.and.(nucelem(i2).ne.elemupdate))then
                    goto 198 
                  endif
!!
!! check if atoms of the element are present, only then do update
                  if(num_atoms_element(i2).gt.0)then 
!!
!! prepare temporary arrays
                    allocate(weights(num_weights_short_atomic_free(i2)))
                    allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!!
!! reduce full weight array and deshortdw to arrays containing only the free weights
                    do i3=1,num_weights_short_atomic_free(i2)
                      weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
                      deshortdw_temp(i3)=deshortdw(i3,1,i2)
                    enddo ! i3
!!
                    if((mpisize.eq.1).or.lompmkl)then
!! updatekalman_para cannot be used for 1 process because the corrmatrix array has a different dimension then
                      if(mpirank.eq.0)then
                        call updatekalman(kaldim(i2),corrdim(i2),&
                          num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
                          weights,deshortdw_temp,&
                          corrmatrix_list(1,i2),errore)
                      endif ! mpirank.eq.0
                      call mpi_bcast(weights,num_weights_short_atomic_free(i2),&
                        mpi_real8,0,mpi_comm_world,mpierror)
!!
                    else ! parallel case
!! updatekalman cannot be used for more than 1 process because the corrmatrix array has a different dimension then
                      call updatekalman_para(paramode,kaldim(i2),corrdim(i2),&
                        num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
                        weights,deshortdw_temp,&
                        corrmatrix_list(1,i2),errore)
                    endif ! mpisize.eq.1
!!
!! apply range restrictions to weights if requested
                    if(restrictw.gt.0.0d0)then
                      do i3=1,num_weights_short_atomic_free(i2)
                        if((weights(i3).gt.(-1.d0*restrictw+1.0d0))&
                          .and.(weights(i3).lt.(restrictw-1.0d0)))then 
                        elseif(weights(i3).ge.(restrictw-1.0d0))then
                          weights(i3)=restrictw-1.0d0+tanh(weights(i3)-restrictw+1.0d0)
                        elseif(weights(i3).le.(-1.d0*restrictw+1.0d0))then
                          weights(i3)=-1.d0*restrictw+1.0d0+tanh(weights(i3)+restrictw-1.0d0)
                        endif
                      enddo ! i3
                    endif
!!
!! expand weights array back to original full array
                    do i3=1,num_weights_short_atomic_free(i2)
                      weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
                    enddo ! i3
!!
!! debugging output of weights_short to debug.data
                    if((mpirank.eq.0).and.(pstring(1:1).eq.'1'))then
                      call debugweights(nelem,0,&
                        countepoch,point,i2,&
                        maxnum_weights_short_atomic,&
                        maxnum_layers_short_atomic,num_layers_short_atomic,&
                        nodes_short_atomic,weights_short_atomic)
                    endif
!!
!! deallocate temporary arrays
                    deallocate(weights)
                    deallocate(deshortdw_temp)
!!
                  endif ! num_atoms_element(i2).gt.0
!!
!! jump here if element is not updated due to lupdatebyelement 
 198              continue
!!
                enddo ! i2=1,nelem
!!
              elseif(optmodee.eq.2)then ! conjugate gradient'
                write(ounit,*)'CG optimization not yet implemented in optimize_short'
                stop !'
!!
              elseif(optmodee.eq.3)then ! steepest descent'
!! update weights for each element
                do i2=1,nelem
!!
!! skip this energy if update is not requested for this element
                  if(lupdatebyelement.and.(nucelem(i2).ne.elemupdate))then
                    goto 199 
                  endif
!!
!! check if atoms of the element are present, only then do update
                  if(num_atoms_element(i2).gt.0)then ! atoms of element are present
!!
!! allocate temporary arrays for free weights
                    allocate(weights(num_weights_short_atomic_free(i2)))
                    allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!!
!! serial update only: 
                    if(mpirank.eq.0)then
!!
!! reduce weights array to free weights
                      do i3=1,num_weights_short_atomic_free(i2)
                        weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
                        deshortdw_temp(i3)=deshortdw(i3,1,i2)
                      enddo ! i3
!!
                      call updatesteepest(&
                        num_weights_short_atomic_free(i2),weights,&
                        deshortdw_temp,errore,steepeststepe)
!!
!! apply restrictions to weights is requested
                      if(restrictw.gt.0.0d0)then
                        do i3=1,num_weights_short_atomic_free(i2)
                          if((weights(i3).gt.(-1.d0*restrictw+1.0d0))&
                            .and.(weights(i3).lt.(restrictw-1.0d0)))then
                          elseif(weights(i3).ge.(restrictw-1.0d0))then
                            weights(i3)=restrictw-1.0d0+tanh(weights(i3)-restrictw+1.0d0)
                          elseif(weights(i3).le.(-1.d0*restrictw+1.0d0))then
                            weights(i3)=-1.d0*restrictw+1.0d0+tanh(weights(i3)+restrictw-1.0d0)
                          endif
                        enddo ! i3
                      endif
!!
!! expand weights array back to original full array
                      do i3=1,num_weights_short_atomic_free(i2)
                        weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
                      enddo ! i3
!!
!! debugging output of weights_short
                      if((mpirank.eq.0).and.(pstring(1:1).eq.'1'))then
                        call debugweights(nelem,0,&
                          countepoch,point,i2,&
                          maxnum_weights_short_atomic,&
                          maxnum_layers_short_atomic,num_layers_short_atomic,&
                          nodes_short_atomic,weights_short_atomic)
                      endif
!!
                    endif ! mpirank.eq.0
!!
!! distribute new weights to all processes
                    call mpi_bcast(weights_short_atomic,maxnum_weights_short_atomic*nelem,&
                      mpi_real8,0,mpi_comm_world,mpierror)
!!
                    deallocate(weights)
                    deallocate(deshortdw_temp)
!!
                  endif ! num_atoms_element(i2).gt.0
!!
!! jump here if update for this element was not requested due to lupdatebyelement 
 199              continue
!!
                enddo ! i2=1,nelem
              endif ! optmodee
!!
!! debug: write length of short range weights to debug.data in two ways
              if(ldebug)then
                do i2=1,nelem
                  debugsum=0.0d0
                  debugsum2=0.0d0
                  do i3=1,num_weights_short_atomic(i2)
                    debugsum=debugsum+weights_short_atomic(i3,i2)**2
                    debugsum2=debugsum2+weights_short_atomic(i3,i2)
                  enddo ! i3
                  debugsum=dsqrt(debugsum)
                  write(debugunit,'(i6,a,i6,x,a2,x,2f14.8)')&
                    countepoch,' weight_short length point ',point,element(i2),debugsum,debugsum2
                enddo ! i2
              endif
!!
!! reinitializations for next update
              nenergy                = 0
              num_atoms_element(:)   = 0
              deshortdw(:,:,:)       = 0.0d0
              deshortdwsum(:,:,:)    = 0.0d0
              errore                 = 0.0d0
              erroresum              = 0.0d0
!!
            endif !((mod(nenergy,nenergygroup).eq.0).or.(i1.eq.npoints))then
            call abstime(timeeupdateend,dayeupdate)
            timeeupdate=timeeupdate+(timeeupdateend-timeeupdatestart)
!!
          else
!!
!! if joint_energy_force_update for E and F is used and E is not used, skip point completely (don't use forces)
            if(ljointefupdate)goto 101
!!
!! point is not used for short update
          endif !(abserror.gt.kalmanthreshold_temp*rmse_short)
!!
!!
!! jump here if energy update of this point is skipped
 98       continue
!!
          call abstime(timeefittingend,dayefitting)
          timeefitting=timeefitting+(timeefittingend-timeefittingstart)
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! force part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
        call abstime(timeffittingstart,dayffitting)
        if(luseforces.and.(forcernd.gt.0.0d0))then                    ! for forcernd=0.0d0 only calculate RMSEs but don't do updates
!!
!! automatically determine the number of forces to group for each element in this structure
          if(lfgroupbystruct) then
            mforcegroup(:)=num_atoms_element_struct(:)*3
            do i2=1,nelem
              if(mforcegroup(i2).eq.0)then    ! in case element is not present
                mforcegroup(i2)=1
              endif
            enddo
          endif
!!
!! in case of joint_energy_force_update group all forces of this structure
          if(ljointefupdate) then
            mforcegroup(:)=num_atoms_element_struct(:)*3
            do i2=1,nelem
              if(mforcegroup(i2).eq.0)then    ! in case element is not present
                mforcegroup(i2)=1
              endif
            enddo
          endif
!!
          do i2=1,nelem
            if(scalefactorf.lt.0.0d0)then
              scalefactorftemp(i2)=(dble(mforcegroup(i2))*energyrnd)/&
                        (dble(3*num_atoms_element_struct(i2)*dble(nenergygroup))*forcernd)
            else ! keep scalefactorf from input.nn
              scalefactorftemp(i2)=scalefactorf
            endif
          enddo
!!
!! allocations for parallel use
          call abstime(timefsymstart,dayfsym)
          allocate(strs(3,3,maxnum_funcvalues_short_atomic,natoms))
          allocate(symfunctiondummy(maxnum_funcvalues_short_atomic,natoms))
!!
          allocate(lsta(2,max_num_atoms))
          allocate(lstc(listdim))
          allocate(lste(listdim))
          allocate(lstb(listdim,4))
          allocate(num_neighbors(natoms))
!!
          call getneighborsatomic(n_start,n_end,natoms,&
            num_atoms,num_neighbors,zelem,max_num_neighbors_atomic,&
            lsta,lstc,lste,&
            maxcutoff_short_atomic,lattice_list(1,1,idx(i1)),xyzstruct_list(1,1,idx(i1)),&
            lstb,lperiodic)
!!
          allocate(dsfuncdxyz_mpi(maxnum_funcvalues_short_atomic,natoms,0:max_num_neighbors_atomic,3))
          dsfuncdxyz_mpi(:,:,:,:)=0.0d0
          allocate(dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_atomic,3))
          allocate(neighboridx(natoms,0:max_num_neighbors_atomic))  
          allocate(invneighboridx(natoms,max_num_atoms))  
          call getneighboridx(n_start,n_end,natoms,listdim,&
            max_num_atoms,max_num_neighbors_atomic,&
            lsta,lstc,neighboridx,invneighboridx)
!!
!! get dsfuncdxyz:
!! this is independent of the weights and needs to be done only once for each point i1
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!! Caution: symmetry functions from this subroutine cannot be used because they are not scaled
!! parallel version
!! calculate the short range symmetry functions for atoms n_start to n_end
          call calconefunction_atomic(cutoff_type,max_num_neighbors_atomic,&
            max_num_atoms,n_start,n_end,natoms,elementindex,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic, &
            nelem,zelem,listdim,&
            lsta,lstc,lste,invneighboridx,&
            function_type_short_atomic,symelement_short_atomic,&
            xyzstruct_list(1,1,idx(i1)),symfunctiondummy,0.0d0,&
            funccutoff_short_atomic,eta_short_atomic,rshift_short_atomic,&
            lambda_short_atomic,zeta_short_atomic,dsfuncdxyz_mpi,strs,lstb,&
            lperiodic,.true.,.false.,lrmin)
          if(.not.lrmin)then
            write(ounit,*)'Error: too close atoms in optimize_short_combined'
            stop !'
          endif
!!
          deallocate(lsta)
          deallocate(lstc)
          deallocate(lste)
          deallocate(lstb)
          deallocate(invneighboridx)  
!!
!! scale dsfuncdxyz
          if(lscalesym)then
!!
!! scale dsfuncdxyz and strs for natoms
!! parallel version
            call scaledsfunc_para(natoms,atomindex,max_num_neighbors_atomic,&
              maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
              nelem,minvalue_short_atomic,maxvalue_short_atomic,&
              scmin_short_atomic,scmax_short_atomic,&
              zelem,dsfuncdxyz_mpi,strs)
          endif ! lscalesym
!!
!! combine dsfuncdxyz of all processes
!! TODO: This should be better done with mpi_gatherv
          dsfuncdxyz(:,:,:,:)=0.0d0
          icount=0
          do i2=n_start,n_end
            icount=icount+1
            dsfuncdxyz(:,i2,:,:)=dsfuncdxyz_mpi(:,icount,:,:)
          enddo
          call mpi_allreduce(mpi_in_place,dsfuncdxyz,&
            maxnum_funcvalues_short_atomic*max_num_atoms*max_num_atoms*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!! 
          deallocate(strs)
          deallocate(dsfuncdxyz_mpi)
          deallocate(symfunctiondummy)
          call abstime(timefsymend,dayfsym)
          timefsym=timefsym+(timefsymend-timefsymstart)
!!
!! calculate the error of the forces individually for each atom and component
          do i2=1,num_atoms
            do i3=1,3 ! loop over x,y and z
!!
!! counter for forces (counts all forces no matter if used for updating or not)
              count_atoms_element_block(elementindex(zelem(i2)))=&
                count_atoms_element_block(elementindex(zelem(i2)))+1
!!
!! if this force is not among the worst and only these should be used, skip this force
              if(luseworstf.and.(.not.lusef(i3,i2,idx(i1))))then 
                goto 99 ! skip this point
              endif
!!
!! if only one element should be updated, skip this force if it belongs to the wrong element
              if(lupdatebyelement.and.(zelem_list(idx(i1),i2).ne.elemupdate))then
                goto 99 ! skip this point
              endif
!!
!! if only a random subset of points is used: decide if we use this point
              z=ran0(kseed)
              if(z.gt.forcernd)then 
                goto 99
              endif
!!
!! if the force is larger than maxforce, don't use it
              if(abs(shortforce(i3,i2)).gt.maxforce)then 
                goto 99
              endif
!!
!! calculate the NN forces nnshortforce for current weight values
              nnshortforce(:,:)=0.0d0
!!
!! this subroutine calculates only the one force that is needed (specified by i2,i3), not the full array nnshortforce
!! parallel version
!! This scales nearly perfectly! 
              call abstime(timeferrorstart,dayferror)
              call getoneshortforce_para(i2,i3,max_num_neighbors_atomic,&
                neighboridx,natoms,atomindex,num_atoms,zelem,&
                symfunction_short_atomic_list(1,1,idx(i1)),&
                dsfuncdxyz,nnshortforce)
!!
!! calculate the error for the current force component
!! We have to normalize the force error here per atom to be consistent with atomic energies
              errorf=-1.d0*(shortforce(i3,i2)-nnshortforce(i3,i2))/dble(num_atoms)  ! LABEL1
              abserrorf=abs(errorf)
              call abstime(timeferrorend,dayferror)
              timeferror=timeferror+(timeferrorend-timeferrorstart)
!!
!! skip force update if error is within noise
              if(abserrorf.lt.noisef)then
                goto 99 ! skip this point
              endif
!!
!! if damping is requested, incorporate squared weights into errorf
!! FIXME: This should be done for all elements present in this structure also for one force update 
!! because one force depends on weights of all elements
              if(ldampw)then
                sumwsquared=0.0d0
                do i4=1,num_weights_short_atomic(elementindex(zelem(i2)))
                  sumwsquared=sumwsquared &
                    +(weights_short_atomic(i4,elementindex(zelem(i2))))**2.0d0
                enddo
                errorf=(1.d0-dampw)*errorf+(dampw*sumwsquared)/dble(num_weights_short_atomic(elementindex(zelem(i2))))
              endif
!!
              errorf=errorf*scalefactorftemp(elementindex(zelem(i2)))
!!
!! set errorthreshold criterion for short range energy
              if(lfixederrorf)then
                errorthresholdf=fixederrorf
              else
                errorthresholdf=kalmanthresholdf_temp*rmse_force_s_ref/dble(num_atoms) ! LABEL1
              endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! decide if force should be used for update
!! In case of ljointefupdate we should use all forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
              if((abserrorf.gt.errorthresholdf)&
                .or.&
                ljointefupdate&
                .or.&
                ((count_atoms_element_block(elementindex(zelem(i2))).eq.num_atoms_element_block(elementindex(zelem(i2))))&
                  .and.(i3.eq.3))&                  ! last force of this element in this block, update may be necessary
                .or.&
                (lfgroupbystruct.and.(count_atoms_element_f(elementindex(zelem(i2))).eq.num_atoms_element_struct(elementindex(zelem(i2))))&
                  .and.(i3.eq.3)))then                                  ! last force of this element in this structure, update may be necessary 
!!
                if((abserrorf.le.errorthresholdf).and.(.not.ljointefupdate)) goto 88
!!
!! include this force in fitting statistics report 
                fitstatf(i3,i2,ndone+idx(i1))=fitstatf(i3,i2,ndone+idx(i1))+1
!! identify which elements should be updated
                do i4=1,nelem
                  if(num_atoms_element_struct(i4).gt.0)then
                    count_atoms_element_f(i4)=&
                      count_atoms_element_f(i4)+1
                  endif
                enddo
                numberf = numberf+1
!!
!! calculate the derivative of the short range forces with respect to the weights
!! dfshortdw is for one atom only, but depends on all atoms
!! calculate dfshortdw
!! serial version
                daydfshortdw      =0
                call abstime(timedfshortdwstart,daydfshortdw)
                if(paramode.eq.1)then
                  if(mpirank.eq.0)then
                    call getdfshortdw(i2,i3,max_num_neighbors_atomic,&
                      num_weights_short_atomic_free,zelem,num_atoms,&
                      wconstraintidx,symfunction_short_atomic_list(1,1,idx(i1)),dsfuncdxyz,&
                      dfshortdw)
                  endif ! mpirank.eq.0
!!
!! CHECK: Do we have to normalize dfshortdw per atom here to be consistent with atomic energies?
                  dfshortdw(:,:)=dfshortdw(:,:)/dble(num_atoms)
!!
!! distribute dfshortdw to all processes
                  call mpi_bcast(dfshortdw,maxnum_weights_short_atomic*nelem,&
                    mpi_real8,0,mpi_comm_world,mpierror)
!!
!! debug: write dfshortdw to debug.data
                  if((pstring(4:4).eq.'1').and.(mpirank.eq.0))then 
                    do i4=1,maxnum_weights_short_atomic
                      do i5=1,nelem
                        write(debugunit,'(3i6,a,3i6,f20.10)')&
                          countepoch,point,i2,&
                          ' dfshortdw ',elementindex(zelem(i2)),i4,i5,dfshortdw(i4,i5)
                      enddo ! i5
                    enddo ! i4
                  endif
!!
                elseif(paramode.eq.2)then
!! parallel version
!! CAUTION: in the parallel version dfshortdw depends slightly on the number of processes
!! because the order of summation matters
!! Maybe this is fixed now by introducing a temporary storage array
!! WARNING: this version still crashes after a while
                  call getdfshortdw_para(i2,i3,max_num_neighbors_atomic,&
                    num_weights_short_atomic_free,natoms,atomindex,&
                    zelem,num_atoms,wconstraintidx,&
                    symfunction_short_atomic_list(1,n_start,idx(i1)),dsfuncdxyz,&
                    dfshortdw)
!! CHECK: Do we have to normalize dfshortdw per atom here to be consistent with atomic energies?
                  dfshortdw(:,:)=dfshortdw(:,:)/dble(num_atoms)
!!
                else
                  write(ounit,*)'Error: unsupported paramode in optimize_short_combined'
                  stop !'
                endif ! paramode
                call abstime(timedfshortdwend,daydfshortdw)
                timedfshortdw=timedfshortdw+(timedfshortdwend-timedfshortdwstart)
!!
!! CHECK: should this be done:? this is new in V61
                dfshortdw(:,:)=dfshortdw(:,:)*scalefactorftemp(elementindex(zelem(i2)))
!!
!! adapt sign, this is needed in case of mforcegroup .gt. 1 to avoid cancellations upon averaging
                if(errorf.lt.0.0d0)then
                  errorf=-1.d0*errorf
                  dfshortdw(:,:)=-1.d0*dfshortdw(:,:)
                endif
!!
!! debug:
!! FIXME: debugsum is calculated only for one element here, by dfshortdw is nonzero for all elements for every force
                if(ldebug)then
                  debugsum=0.0d0
                  debugsum2=0.0d0
                  if(i3.eq.1)then
                    write(debugunit,'(i6,a,i6,a6,i5,a3,f14.8)')countepoch,' errorf of point ',point,' atom ',i2,' x ',errorf
                  elseif(i3.eq.2)then
                    write(debugunit,'(i6,a,i6,a6,i5,a3,f14.8)')countepoch,' errorf of point ',point,' atom ',i2,' y ',errorf
                  elseif(i3.eq.3)then
                    write(debugunit,'(i6,a,i6,a6,i5,a3,f14.8)')countepoch,' errorf of point ',point,' atom ',i2,' z ',errorf
                  endif
                  do i4=1,num_weights_short_atomic(elementindex(zelem(i2)))
                    debugsum=debugsum+dfshortdw(i4,elementindex(zelem(i2)))**2
                    debugsum2=debugsum2+dfshortdw(i4,elementindex(zelem(i2)))
                  enddo ! i4
                  debugsum=dsqrt(debugsum)
                  if(i3.eq.1)then
                    write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                    countepoch,' dfshortdw length point ',point,' atom ',i2,'x',element(elementindex(zelem(i2))),debugsum,debugsum2
                  elseif(i3.eq.2)then
                    write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                    countepoch,' dfshortdw length point ',point,' atom ',i2,'y',element(elementindex(zelem(i2))),debugsum,debugsum2
                  elseif(i3.eq.3)then
                    write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                    countepoch,' dfshortdw length point ',point,' atom ',i2,'z',element(elementindex(zelem(i2))),debugsum,debugsum2
                  endif
                endif !'
!!
!! if damping is requested:
                if(ldampw) then
                  do i4=1,num_weights_short_atomic_free(elementindex(zelem(i2)))
                    do i5=1,nelem
                      dfshortdw(i4,i5)=(1.d0-dampw)*dfshortdw(i4,i5) &
                        -dampw*2.0d0*weights_short_atomic(wconstraintidx(i4,i5),i5)/dble(num_weights_short_atomic_free(i5))
                    enddo ! i5'
                  enddo ! i4'
                endif
!!
!! accumulate dfshortdw and errorf for grouping
                dfshortdwsum(:,:)=dfshortdwsum(:,:)+dfshortdw(:,:)
                errorfsum(elementindex(zelem(i2)))=&
                  errorfsum(elementindex(zelem(i2)))+errorf
!!
!! in case of joint_energy_force_update skip individual force updates
                if(ljointefupdate)goto 99 
!!
 88             continue
!! check if there are forces left not used for update (standard grouping case)
                do i4=1,nelem
                  if((count_atoms_element_block(i4).eq.num_atoms_element_block(i4))&   ! is this the last atom of this element in this block?
                     .and.(i3.eq.3)&                                                   ! is this also the last force component for this atom? 
                     .and.(count_atoms_element_f(i4).gt.0))then                        ! there are some so far unused forces 
                       lenforceupdatef(i4)=.true.
!!                       write(ounit,*)'setting lenforceupdatef'
                  endif
                enddo ! i4 
!!
!! check if there are forces left not used for update (group_forces_by_structure case)
                do i4=1,nelem
                  if((lfgroupbystruct)&                                               ! do we group over all atoms within one structure?
                    .and.(count_atoms_element_f(i4).eq.num_atoms_element_struct(i4))& ! is this the last atom of this element in this structure
                    .and.(i3.eq.3))then                                               ! is this also the last force component for this atom?
                    if(count_atoms_element_f(i4).gt.0)then
                      lenforceupdatef(i4)=.true.
!!                      write(ounit,*)'setting lenforceupdatef'
                    endif
                  endif
                enddo ! i4 
!!
!! do force weight update for each element in the present structure (not just for the element of the force component!)
                do i4=1,nelem
!!
                  call abstime(timefupdatestart,dayfupdate)
                  if(((mod(count_atoms_element_f(i4),mforcegroup(i4)).eq.0).or.&     ! number of forces for grouping has been reached for this element
                    (lenforceupdatef(i4)))&                                      ! last force for this element in this block of structures
                    .and.(count_atoms_element_f(i4).gt.0))then                   ! CAUTION: mod is also 0 for count_atoms_element_f=0
!!
                    lenforceupdatef(i4)=.false.
!!                    write(ounit,*)'updating forces ',i1,i2,i3,element(i4)
!!
!! do update for this element i4 only if atoms are present 
!!
!! allocate temporary arrays for free weights
                    allocate(weights(num_weights_short_atomic_free(i4)))
                    allocate(dfshortdw_temp(num_weights_short_atomic_free(i4)))
                    allocate(deshortdw_temp(num_weights_short_atomic_free(i4)))
!! 
!! debug:
                    if((mpirank.eq.0).and.(pstring(1:1).eq.'1')) &
                      write(debugunit,'(a,4i8,f14.6)')&
                      'force update done ',i2,i3,i4,count_atoms_element_f(i4),errorf
!!'
!! normalize by the number of accumulated forces 
                    dfshortdw(:,i4)=dfshortdwsum(:,i4)/&
                      (dble(count_atoms_element_f(i4))) 
                    errorf=errorfsum(i4)/&
                      (dble(count_atoms_element_f(i4))) 
!!
!! Kalman damping
                    errorf=errorf*kalman_dampf
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! reduce weights array to free weights
                    do i5=1,num_weights_short_atomic_free(i4)
                      weights(i5)=weights_short_atomic(wconstraintidx(i5,i4),i4)
                      dfshortdw_temp(i5)=dfshortdw(i5,i4)
                    enddo ! i5
!!
                    if(optmodef.eq.1)then ! Kalman filter
!!
                      if(.not.lsepkalman)then
                        corrmatrixf_list(:,i4)=corrmatrix_list(:,i4)
                      endif
!!
                      if((mpisize.eq.1).or.lompmkl)then
                        if(mpirank.eq.0)then
                          call updatekalman(kaldim(i4),corrfdim(i4),&
                            num_weights_short_atomic_free(i4),kalmanlambda(i4),kalmannue,&
                            weights,dfshortdw_temp,&
                            corrmatrixf_list(1,i4),&
                            errorf)
                        endif ! mpirank.eq.0
                        call mpi_bcast(weights,num_weights_short_atomic_free(i4),&
                          mpi_real8,0,mpi_comm_world,mpierror)
                      else
                        call updatekalman_para(paramode,kaldim(i4),corrfdim(i4),&
                          num_weights_short_atomic_free(i4),kalmanlambda(i4),kalmannue,&
                          weights,dfshortdw_temp,&
                          corrmatrixf_list(1,i4),&
                          errorf)
                      endif ! mpisize.eq.1
!!
                      if(.not.lsepkalman)then
                        corrmatrix_list(:,i4)=corrmatrixf_list(:,i4)
                      endif
!!
!! Update again the weights using the energy after each individual force update if requested:
                      call abstime(timeefittingrepeatstart,dayefittingrepeat)
                      if(lrepeate.and.(shortenergy_list(idx(i1)).lt.maxenergy))then
                        if(lupdatebyelement.and.(nucelem(i4).ne.elemupdate))then
!! Caution: a force on atom of element A depends also on energies of elements B and C
!! => here we update only weights of A, therefore a small force error cannot be reached
!! still we do this here to change only the weights of element A
                          goto 200 ! skip the update for this element
                        endif
                        call abstime(timeeerrorrepeatstart,dayeerrorrepeat)
                        nnatomenergy(:)= 0.0d0
                        nntotalenergy  = 0.0d0
                        do i5=1,num_weights_short_atomic_free(i4)
                          weights_short_atomic(wconstraintidx(i5,i4),i4)=weights(i5)
                        enddo ! i5
                        call calconeshort_para(point,natoms,atomindex,&
                          zelem,&
                          symfunction_short_atomic_list(1,n_start,idx(i1)),nnatomenergy,&
                          nntotalenergy)
                        call mpi_allreduce(mpi_in_place,nnatomenergy,&
                          max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
                        nntotalenergy=0.0d0
                        do i5=1,num_atoms
                          nntotalenergy=nntotalenergy+nnatomenergy(i5)
                        enddo
                        nneshort=nntotalenergy/dble(num_atoms)
                        errore=shortenergy_list(idx(i1))-nneshort
                        call abstime(timeeerrorrepeatend,dayeerrorrepeat)
                        timeeerrorrepeat=timeeerrorrepeat+(timeeerrorrepeatend-timeeerrorrepeatstart)
                        errore=errore*kalman_dampe
                        daydeshortdw      =0
                        call abstime(timedeshortdwrepeatstart,daydeshortdwrepeat)
                        call getdeshortdw(&
                          num_weights_short_atomic_free,&
                          zelem,num_atoms,wconstraintidx,&
                          symfunction_short_atomic_list(1,1,idx(i1)),deshortdw)
                        call abstime(timedeshortdwrepeatend,daydeshortdwrepeat)
                        timedeshortdwrepeat=timedeshortdwrepeat&
                          +(timedeshortdwrepeatend-timedeshortdwrepeatstart)
!! adapt sign 
                        if(errore.lt.0.0d0)then
                          errore=-1.d0*errore
                          deshortdw(:,:,:)=-1.d0*deshortdw(:,:,:)
                        endif
                        do i5=1,num_weights_short_atomic_free(i4)
                          deshortdw_temp(i5)=deshortdw(i5,1,i4)
                        enddo ! i5
                        call abstime(timeeupdaterepeatstart,dayeupdaterepeat)
                        if((mpisize.eq.1).or.lompmkl)then
                          call updatekalman(kaldim(i4),corrdim(i4),&
                            num_weights_short_atomic_free(i4),kalmanlambda(i4),kalmannue,&
                            weights,deshortdw_temp,&
                            corrmatrix_list(1,i4),errore)
                        else ! real parallel case
                          call updatekalman_para(paramode,kaldim(i4),corrdim(i4),&
                            num_weights_short_atomic_free(i4),kalmanlambda(i4),kalmannue,&
                            weights,deshortdw_temp,&
                            corrmatrix_list(1,i4),errore)
                        endif ! mpisize.eq.1
                        call abstime(timeeupdaterepeatend,dayeupdaterepeat)
                        timeeupdaterepeat=timeeupdaterepeat+(timeeupdaterepeatend-timeeupdaterepeatstart)
                        numbere=numbere+1
 200                    continue
                      endif ! lrepeate
                      call abstime(timeefittingrepeatend,dayefittingrepeat)
                      timeefittingrepeat=timeefittingrepeat+(timeefittingrepeatend-timeefittingrepeatstart)
!!
                    elseif(optmodef.eq.2)then ! conjugate gradient'
                      write(ounit,*)'CG optimization not yet implemented in fitforcesshort'
                      stop !'
!!
                    elseif(optmodef.eq.3)then ! steepest descent'
!!
                      if(mpirank.eq.0)then
                        call updatesteepest(&
                          num_weights_short_atomic_free(i4),weights,&
                          dfshortdw_temp,errorf,steepeststepf)
                      endif ! mpirank.eq.0
!!
!! Update again the weights using the energy after each individual force update if requested:
                      call abstime(timeefittingrepeatstart,dayefittingrepeat)
                      if(lrepeate.and.(shortenergy_list(idx(i1)).lt.maxenergy))then
                        if(lupdatebyelement.and.(nucelem(i4).ne.elemupdate))then
!! Caution: a force on atom of element A depends also on energies of elements B and C
!! => here we update only weights of A, therefore a small force error cannot be reached
!! still we do this here to change only the weights of element A
                          goto 201 ! skip the update for this element
                        endif
                        call abstime(timeeerrorrepeatstart,dayeerrorrepeat)
                        nnatomenergy(:)= 0.0d0
                        nntotalenergy  = 0.0d0
                        do i5=1,num_weights_short_atomic_free(i4)
                          weights_short_atomic(wconstraintidx(i5,i4),i4)=weights(i5)
                        enddo ! i5
                        call calconeshort_para(point,natoms,atomindex,&
                          zelem,&
                          symfunction_short_atomic_list(1,n_start,idx(i1)),nnatomenergy,&
                          nntotalenergy)
                        call mpi_allreduce(mpi_in_place,nnatomenergy,&
                          max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
                        nntotalenergy=0.0d0
                        do i5=1,num_atoms
                          nntotalenergy=nntotalenergy+nnatomenergy(i5)
                        enddo
                        nneshort=nntotalenergy/dble(num_atoms)
                        errore=shortenergy_list(idx(i1))-nneshort
                        errore=errore*kalman_dampe
                        call abstime(timeeerrorrepeatend,dayeerrorrepeat)
                        timeeerrorrepeat=timeeerrorrepeat+(timeeerrorrepeatend-timeeerrorrepeatstart)
                        daydeshortdw      =0
                        call abstime(timedeshortdwrepeatstart,daydeshortdwrepeat)
                        call getdeshortdw(&
                          num_weights_short_atomic_free,&
                          zelem,num_atoms,wconstraintidx,&
                          symfunction_short_atomic_list(1,1,idx(i1)),deshortdw)
                        call abstime(timedeshortdwrepeatend,daydeshortdwrepeat)
                        timedeshortdwrepeat=timedeshortdwrepeat&
                          +(timedeshortdwrepeatend-timedeshortdwrepeatstart)
!! adapt sign 
                        if(errore.lt.0.0d0)then
                          errore=-1.d0*errore
                          deshortdw(:,:,:)=-1.d0*deshortdw(:,:,:)
                        endif
                        do i5=1,num_weights_short_atomic_free(i4)
                          deshortdw_temp(i5)=deshortdw(i5,1,i4)
                        enddo ! i5
                        call abstime(timeeupdaterepeatstart,dayeupdaterepeat)
                        call updatesteepest(&
                          num_weights_short_atomic_free(i4),&
                          weights,deshortdw_temp,errore,steepeststepe)
                        call abstime(timeeupdaterepeatend,dayeupdaterepeat)
                        timeeupdaterepeat=timeeupdaterepeat+(timeeupdaterepeatend-timeeupdaterepeatstart)
                        numbere=numbere+1
 201                    continue
                      endif ! lrepeate 
                      call abstime(timeefittingrepeatend,dayefittingrepeat)
                      timeefittingrepeat=timeefittingrepeat+(timeefittingrepeatend-timeefittingrepeatstart)
!!
                    endif ! optmodef
!!
!!
!! apply restrictions to weights is requested
                    if(restrictw.gt.0.0d0)then
                      do i5=1,num_weights_short_atomic_free(i4)
                        if((weights(i5).gt.(-1.d0*restrictw+1.0d0))&
                          .and.(weights(i5).lt.(restrictw-1.0d0)))then
                        elseif(weights(i5).ge.(restrictw-1.0d0))then
                          weights(i5)=restrictw-1.0d0+tanh(weights(i5)-restrictw+1.0d0)
                        elseif(weights(i5).le.(-1.d0*restrictw+1.0d0))then
                          weights(i5)=-1.d0*restrictw+1.0d0+tanh(weights(i5)+restrictw-1.0d0)
                        endif
                      enddo ! i5
                    endif
!!
!! expand weights array back to full arrays 
                    do i5=1,num_weights_short_atomic_free(i4)
                      weights_short_atomic(wconstraintidx(i5,i4),i4)=weights(i5)
                    enddo ! i5
                    call mpi_bcast(weights_short_atomic,maxnum_weights_short_atomic*nelem,&
                      mpi_real8,0,mpi_comm_world,mpierror)
!!
!! debug:
                    if(ldebug)then
                      debugsum=0.0d0
                      debugsum2=0.0d0
                      do i5=1,num_weights_short_atomic(i4)
                        debugsum=debugsum+weights_short_atomic(i5,i4)**2
                        debugsum2=debugsum2+weights_short_atomic(i5,i4)
                      enddo ! i5
                      debugsum=dsqrt(debugsum)
                      if(i3.eq.1)then
                        write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                          countepoch,' weights_short_atomic length point ',point,' atom ',i2,'x',element(i4),debugsum,debugsum2
                      elseif(i3.eq.2)then
                        write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                          countepoch,' weights_short_atomic length point ',point,' atom ',i2,'y',element(i4),debugsum,debugsum2
                      elseif(i3.eq.3)then
                        write(debugunit,'(i6,a,i6,x,a6,x,i5,x,a1,x,a2,x,2f14.8)')&
                          countepoch,' weights_short_atomic length point ',point,' atom ',i2,'z',element(i4),debugsum,debugsum2
                      endif !'
                    endif
!!
!! debugging output of weights_short
                    if((mpirank.eq.0).and.(pstring(1:1).eq.'1'))then
                      call debugweights(nelem,0,&
                        countepoch,point,i4,maxnum_weights_short_atomic,&
                        maxnum_layers_short_atomic,num_layers_short_atomic,&
                        nodes_short_atomic,weights_short_atomic)
                    endif
!!
                    deallocate(weights)
                    deallocate(dfshortdw_temp)
                    deallocate(deshortdw_temp)
!!
!!
!! reinitializations for next update
                    dfshortdw(:,i4)          = 0.0d0
                    dfshortdwsum(:,i4)       = 0.0d0
                    errorf                   = 0.0d0
                    errorfsum(i4)            = 0.0d0
                    count_atoms_element_f(i4)= 0
!!
                  endif !(mod(count_stoms_element_f,mforcegroup).eq.0)then
                  call abstime(timefupdateend,dayfupdate)
                  timefupdate=timefupdate+(timefupdateend-timefupdatestart)
                enddo ! i4=1,nelem
!!
              else ! point is not used for short update
!!
              endif !(abserrorf.gt.kalmanthresholdf_temp*rmse_force_s_ref)
!!
!! jump here if force should be skipped for weight update
 99           continue
!!
            enddo ! i3 loop over x,y and z
!!
          enddo ! i2 loop over atoms
          deallocate(neighboridx)  
          deallocate(dsfuncdxyz)
          deallocate(num_neighbors)
!!
        endif ! luseforces
!!
        call abstime(timeffittingend,dayffitting)
        timeffitting=timeffitting+(timeffittingend-timeffittingstart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end of regular update 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! do joint update of energy and forces here if requested
        if(ljointefupdate)then
!!
!! normalization of force error and derivative
          do i2=1,nelem
            if(num_atoms_element(i2).gt.0)then
              errorfsum(i2)=errorfsum(i2)/dble(3*num_atoms_element(i2))
              dfshortdwsum(:,i2)=dfshortdwsum(:,i2)/dble(3*num_atoms_element(i2))
            endif
          enddo
!!
!! get joint error
          erroresum    =erroresum*kalman_dampe
          errorfsum(:) =errorfsum(:)*kalman_dampf
          errorjoint(:)=0.0d0
          do i2=1,nelem
            errorjoint(i2)=(erroresum+errorfsum(i2))/2.d0
          enddo
!!
!! get joint derivative 
          deshortdw(:,:,:)=deshortdwsum(:,:,:)
          do i2=1,nelem
            deshortdw(:,1,i2)=(deshortdw(:,1,i2)+dfshortdwsum(:,i2))/2.d0
          enddo ! i2
!!
          do i2=1,nelem
            allocate(weights(num_weights_short_atomic_free(i2)))
            allocate(deshortdw_temp(num_weights_short_atomic_free(i2)))
!! 
!! skip update if only weights for another element should be optimized
            if(lupdatebyelement.and.(nucelem(i2).ne.elemupdate))then
              goto 100 ! skip the update for this element
            endif
!!
            errore=errorjoint(i2)
!! chose optimization algorithm
            if(optmodee.eq.1)then ! Kalman filter
!! update weights for each element
!! check if atoms of the element are present, only then do update
              if(num_atoms_element(i2).gt.0)then
!! reduce weights array to free weights
                do i3=1,num_weights_short_atomic_free(i2)
                  weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
                  deshortdw_temp(i3)=deshortdw(i3,1,i2)
                enddo ! i3
!! Kalman damping
                errore=errore*kalman_dampe
!!
                if((mpisize.eq.1).or.lompmkl)then
                  call updatekalman(kaldim(i2),corrdim(i2),&
                    num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
                    weights,deshortdw_temp,&
                    corrmatrix_list(1,i2),errore)
                else ! real parallel case
                  call updatekalman_para(paramode,kaldim(i2),corrdim(i2),&
                    num_weights_short_atomic_free(i2),kalmanlambda(i2),kalmannue,&
                    weights,deshortdw_temp,&
                    corrmatrix_list(1,i2),errore)
                endif ! mpisize.eq.1
!! expand weights array back to original array
                do i3=1,num_weights_short_atomic_free(i2)
                  weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
                enddo ! i3
!!
              endif ! num_atoms_element(i2).gt.0
!!
            elseif(optmodee.eq.2)then ! conjugate gradient'
              write(ounit,*)'CG optimization not yet implemented in optimize_short'
              stop
            elseif(optmodee.eq.3)then ! steepest descent'
!! check if atoms of the element are present, only then do update
              if(num_atoms_element(i2).gt.0)then ! atoms of element are present
                if(mpirank.eq.0)then
!! reduce weights array to free weights
                  do i3=1,num_weights_short_atomic_free(i2)
                    weights(i3)=weights_short_atomic(wconstraintidx(i3,i2),i2)
                    deshortdw_temp(i3)=deshortdw(i3,1,i2)
                  enddo ! i3
                  call updatesteepest(&
                    num_weights_short_atomic_free(i2),weights,&
                    deshortdw_temp,errore,steepeststepe)
!! expand weights array back to original array
                  do i3=1,num_weights_short_atomic_free(i2)
                    weights_short_atomic(wconstraintidx(i3,i2),i2)=weights(i3)
                  enddo ! i3
                endif ! mpirank.eq.0
                call mpi_bcast(weights_short_atomic,maxnum_weights_short_atomic*nelem,&
                  mpi_real8,0,mpi_comm_world,mpierror)
              endif ! num_atoms_element(i2).gt.0
            endif ! optmodee
!!
!! debugging output of weights_short
            if((mpirank.eq.0).and.(pstring(1:1).eq.'1'))then
              call debugweights(nelem,0,&
                countepoch,point,i2,maxnum_weights_short_atomic,&
                maxnum_layers_short_atomic,num_layers_short_atomic,&
                nodes_short_atomic,weights_short_atomic)
            endif
!!
 100        continue
            deallocate(weights)
            deallocate(deshortdw_temp)
          enddo ! i2=1,nelem
!!
!! reinitializations for next energy update
          nenergy                = 0
          deshortdw(:,:,:)       = 0.0d0
          deshortdwsum(:,:,:)    = 0.0d0
          errore                 = 0.0d0
          erroresum              = 0.0d0
!! reinitializations for next force update
          dfshortdw(:,:)        = 0.0d0
          dfshortdwsum(:,:)     = 0.0d0
          errorf                = 0.0d0
          errorfsum(:)          = 0.0d0
          count_atoms_element_f(:)= 0
!!
        endif ! ljointefupdate
!!
!! jump here if joint_energy_force_update is used and due to the energy neither energy nor forces are used for update
 101    continue
!!
!!
        deallocate(atomindex)      
      enddo ! i1=1,npoints
!!
      return
      end
