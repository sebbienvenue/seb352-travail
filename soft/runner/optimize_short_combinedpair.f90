!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90
!!
!! This routine is called once for each block of points
!!
      subroutine optimize_short_combinedpair(npoints,point,idx,&
         kaldim,maxcorrdim,maxcorrfdim,corrdim,corrfdim,countepoch,&
         num_weightspairfree,ndone,&
         wconstraintidxp,numbere,numberf,&
         lseed,ntrain,fitstat,fitstatf,kseed,&
         kalmanthreshold_temp,kalmanthresholdf_temp,&
         rmse_short,rmse_force_s_ref,&
         corrmatrix_list,corrmatrixf_list,&
         minvalue_short_pair,maxvalue_short_pair)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use symfunctions
      use nnshort_pair
      use structures
!!
      implicit none
!!
      integer npoints                                           ! in
      integer countepoch                                        ! in
      integer num_weightspairfree(npairs)                       ! in
      integer totnum_weightspairfree                            ! internal
      integer totnum_weightspair                                ! internal
      integer wconstraintidxp(maxnum_weights_short_pair,npairs)        ! in
      integer zelem(max_num_atoms)                              ! internal
      integer zelemp(2,max_num_pairs)                           ! internal
      integer num_atoms                                         ! internal
      integer num_pairs_para                                    ! internal
      integer num_pairs                                         ! internal
      integer point                                             ! in/out
      integer idx(nblock)                                       ! in
      integer kseed                                             ! in/out
      integer lseed                                             ! in/out
      integer i1,i2,i3,i4,i5                                    ! internal
      integer day                                               ! internal
      integer mforcegroup                                       ! internal
      integer nenergy                                           ! internal
      integer n_start                                           ! internal
      integer n_end                                             ! internal
      integer, dimension(:), allocatable :: pindex              ! internal
      integer icount                                            ! internal
      integer nforce                                            ! internal
      integer numbere                                           ! in/out
      integer numberf                                           ! in/out
      integer ntrain                                            ! in
      integer fitstat(ntrain)                                   ! in/out
      integer fitstatf(3,max_num_atoms,ntrain)                  ! in/out
      integer maxcorrdim                                        ! in
      integer maxcorrfdim                                       ! in
      integer corrdim(npairs)                                   ! in
      integer corrfdim(npairs)                                  ! in
      integer kaldim(npairs)                                    ! in
      integer pairs_charge(2,listdim)                           ! internal 
      integer num_atoms_element_block(nelem)                    ! internal
      integer num_atoms_element_struct(nelem)                   ! internal
      integer num_pairs_element(npairs)                         ! internal 
      integer num_pairs_element_f(npairs)                       ! internal 
      integer num_pairs_element_struct(npairs)                  ! internal 
      integer count_atoms_element_block(nelem)                  ! internal
      integer count_forces                                      ! internal
      integer ndone                                             ! in
      integer totnumatoms                                       ! internal
!!
      real*8 kalmanthreshold_temp                               ! in
      real*8 kalmanthresholdf_temp                              ! in
      real*8 rmse_short                                         ! in
      real*8 nneshort                                           ! internal
      real*8 errore                                             ! internal
      real*8 erroresum                                          ! internal
      real*8 abserrore                                          ! internal
      real*8 errorf                                             ! internal
      real*8 errorfsum                                          ! internal
      real*8 abserrorf                                          ! internal
      real*8 nnpairenergy(max_num_pairs)                        ! internal
      real*8 errorthreshold                                     ! internal
      real*8 errorthresholdf                                    ! internal
      real*8, dimension(:)  , allocatable :: weightsp                       ! internal
      real*8, dimension(:,:)  , allocatable :: symfunctiondummy             ! internal
      real*8 nneshort_list(nblock)                                          ! internal
      real*8 depairdw(maxnum_weights_short_pair,1,npairs)                   ! internal       CHECK
      real*8 depairdwsum(maxnum_weights_short_pair,1,npairs)                ! internal       CHECK
      real*8, dimension(:) , allocatable :: depairdw_temp                   ! internal 
      real*8 corrmatrix_list(maxcorrdim,npairs)                             ! in/out
      real*8 corrmatrixf_list(maxcorrfdim,npairs)                           ! in/out
      real*8 z,ran0                                                         ! internal
      real*8 nntotalenergy                                                  ! internal
      real*8 nnshortforce_list(3,max_num_atoms,nblock)                      ! internal
      real*8 energyerror_list(nblock)                                       ! internal
      real*8 forceerror_list(3,max_num_atoms,nblock)                        ! internal
      real*8 rmse_force_s_ref                                               ! in
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)       ! in
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)       ! in
      real*8 shortforce(3,max_num_atoms)                                    ! internal
      real*8 nnshortforce(3,max_num_atoms)                                  ! internal
      real*8 dfpairdw(maxnum_weights_short_pair,npairs)                            ! internal       CHECK
      real*8 dfpairdwsum(maxnum_weights_short_pair,npairs)                         ! internal       CHECK
      real*8, dimension(:) , allocatable :: dfpairdw_temp                   ! internal
      real*8 dsfuncdxyz_pair(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3)   ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: dsfuncdxyz_mpi           ! internal
      real*8, dimension(:,:,:,:)  , allocatable :: strs                     ! internal dummy
      real*8 debugsum                                                       ! internal debugging
      real*8 debugsum2                                                      ! internal debugging
      real*8 scalefactorftemp(nelem)                                        ! internal
      real*8 sumwsquared                                                    ! internal
!!
      logical lperiodic                                                     ! internal
      logical lusee(nblock)                                                 ! internal 
      logical lusef(3,max_num_atoms,nblock)                                 ! internal 
      logical lrmin                                                         ! internal
      logical lenforceupdatee                                               ! internal
      logical lenforceupdatef                                               ! internal
!!
!! initializations
      num_pairs_element(:)       = 0
      num_pairs_element_f(:)     = 0
      num_pairs_element_struct(:)= 0
      num_atoms_element_struct(:)= 0
      num_atoms_element_block(:) = 0
      count_atoms_element_block(:) = 0
      count_forces               = 0
      nneshort                   = 0.0d0
      nntotalenergy              = 0.0d0
      shortforce(:,:)            = 0.0d0
      nnshortforce(:,:)          = 0.0d0
      depairdw(:,:,:)            = 0.0d0
      depairdwsum(:,:,:)         = 0.0d0
      dfpairdw(:,:)              = 0.0d0
      dfpairdwsum(:,:)           = 0.0d0
      errore                     = 0.0d0
      erroresum                  = 0.0d0
      errorf                     = 0.0d0
      errorfsum                  = 0.0d0
      nenergy                    = 0
      nforce                     = 0
      lenforceupdatee            = .false.
      lenforceupdatef            = .false.
!! timing variables
      day                        = 0
!!
!! count total number of atoms of each element for this block of points
      do i1=1,npoints
        do i2=1,num_atoms_list(i1)
           num_atoms_element_block(elementindex(zelem_list(i1,i2)))=&
             num_atoms_element_block(elementindex(zelem_list(i1,i2)))+1
        enddo ! i2
      enddo ! i1
!!
!! determine mforcegroup and make sure it is not larger than the number of forces in this block of structures
      mforcegroup = nforcegroup
      totnumatoms=0
      do i1=1,nelem
        totnumatoms=totnumatoms+num_atoms_element_block(i1)
      enddo
      mforcegroup=min(3*totnumatoms,nforcegroup)
      if(mforcegroup.eq.0)then    ! in case element is not present
        write(ounit,*)'Error in optimize_short_combinedpair, mforcegroup =0' 
        stop !'
      endif
!!
!! initialize total number of free short range weights
      totnum_weightspairfree = 0
      do i1=1,npairs 
        totnum_weightspairfree = totnum_weightspairfree + num_weightspairfree(i1)
      enddo
!!
!! count total number of short range weights
      totnum_weightspair = 0
      do i1=1,npairs 
        totnum_weightspair = totnum_weightspair + num_weights_short_pair(i1)
      enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If requested, determine the points with the worst energies of this block of data
      if(luseworste)then
        if(worste.ge.1.0d0)then
          lusee(:)=.true. ! use all points
        else
!! determine the error of the short range energies
         call geteshortpair(npoints,&
           zelemp_list,num_atoms_list,num_pairs_list,&
           symfunction_short_pair_list,nneshort_list)
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
          call getallshortforcespair(nblock,npoints,&
            num_atoms_list,zelem_list,zelemp_list,&
            symfunction_short_pair_list,nnshortforce_list,lattice_list,&
            xyzstruct_list,minvalue_short_pair,&
            maxvalue_short_pair,lperiodic_list)
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
        num_pairs        = num_pairs_list(idx(i1))
        lperiodic        = lperiodic_list(idx(i1))
        shortforce(:,:)  = shortforce_list(:,:,idx(i1))
        zelem(:)         = zelem_list(idx(i1),:)
        zelemp(:,:)      = zelemp_list(:,idx(i1),:)
!!
!! calculate the number of pairs, because we don't want to
!! update weights for pairs not being present
        num_pairs_element_struct(:)=0
        do i2=1,num_pairs 
          num_pairs_element_struct(pairindex(zelemp(1,i2),zelemp(2,i2)))=&
            num_pairs_element_struct(pairindex(zelemp(1,i2),zelemp(2,i2)))+1
        enddo
        num_atoms_element_struct(:)=0
        do i2=1,num_atoms 
          num_atoms_element_struct(elementindex(zelem(i2)))=&
            num_atoms_element_struct(elementindex(zelem(i2)))+1
        enddo
!!
!! determine which atoms of this structure should be calculated by this process
        call mpifitdistribution(num_pairs,num_pairs_para,n_start,n_end)
!!
!! determine the pindex array
        allocate(pindex(num_pairs_para))
        do i2=1,num_pairs_para
          pindex(i2)=n_start+i2-1
        enddo
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! energy part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! skip this energy if it is not among the worst energies and only these should be used
        if(luseworste.and.(.not.lusee(idx(i1))))then ! continue
          goto 98 ! skip this point
        endif
!!
!! if only a random subset of points is used: decide if we use this point
        z=ran0(lseed)
        if(z.gt.energyrnd)then
          goto 98 ! skip this point
        endif
!!
!! skip this energy if it is too high
        if(shortenergy_list(idx(i1)).gt.maxenergy)then
          goto 98 !
        endif
!!
!! predict the short-range pair energies for pairs n_start to n_end
        nnpairenergy(:)= 0.0d0
        nntotalenergy  = 0.0d0
        call calconeshort_parapair(point,num_pairs_para,pindex,&
          zelemp,symfunction_short_pair_list(1,1,idx(i1)),nnpairenergy,&
          nntotalenergy)
!!
!! merge all pair energies
        call mpi_allreduce(mpi_in_place,nnpairenergy,&
          max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
!! CAUTION: In parallel runs the result for nntotalenergy very slightly
!! depends on the number of processes because the order of the
!! summation of the nnpairenergies matters. Even invisible changes after the 15th digit
!! still sum up to significant errors in long fits
!! => don't use the nntotalenergy calculated above but recalculate it in
!! a well-defined way here
        nntotalenergy=0.0d0
        do i2=1,num_pairs 
          nntotalenergy=nntotalenergy+nnpairenergy(i2)
        enddo
!!
!! calculate short range energy eshort and normalize per atom
        nneshort=nntotalenergy/dble(num_atoms)
!!
!! calculate the error (per pair) of this training point
        errore=(shortenergy_list(idx(i1))-nneshort)*dble(num_atoms)/dble(num_pairs)
        abserrore=abs(errore)
!!
!! skip energy update if energy error is within noise
        if(abserrore.lt.noisee)then
          goto 98
        endif
!!
!! if weight damping is used modify the energy error (add normalized sum of squared weights)
        if(ldampw)then
          sumwsquared=0.0d0
          do i4=1,npairs
            do i3=1,num_weights_short_pair(i4)
              sumwsquared=sumwsquared +(weights_short_pair(i3,i4))**2.0d0
            enddo ! i3
          enddo ! i4
          errore=(1.d0-dampw)*errore + (dampw*sumwsquared)/dble(totnum_weightspair)
        endif ! ldampw
!!
!! set errorthreshold criterion for short range energy
        if(lfixederrore)then
          errorthreshold=fixederrore*dble(num_atoms)/dble(num_pairs) ! errorthreshold per pair
        else
          errorthreshold=kalmanthreshold_temp*rmse_short *dble(num_atoms)/dble(num_pairs) ! errorthreshold per pair
        endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! decide if point is sufficiently bad to do weight update
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
!! calculate the derivative of the short range energy with respect to the weights
          if(paramode.eq.1)then
!! serial version
            call getdepairdw(num_weightspairfree,&
              zelemp,num_pairs,wconstraintidxp,&
              symfunction_short_pair_list(1,1,idx(i1)),depairdw)
          else
            write(ounit,*)'Error: paramode not allowed in optimize_short_combined'
            stop !'
          endif ! paramode
!!
!! debug: write deshortdw to debug.data
          if((pstring(3:3).eq.'1').and.(mpirank.eq.0))then
            do i3=1,npairs
              do i4=1,maxnum_weights_short_pair
                write(debugunit,'(2i6,a,2i6,f20.10)')&
                  countepoch,point,&
                  ' depairdw ',i3,i4,depairdw(i4,1,i3)
              enddo
            enddo
          endif
!!
!! debug: write length of depairdw in two ways to debug.data
          if(ldebug)then
            write(debugunit,'(i6,a,i6,x,f14.8)')countepoch,&
              ' errore of point ',point,errore
            do i3=1,npairs !'
              debugsum =0.0d0
              debugsum2=0.0d0
              do i4=1,num_weights_short_pair(i3)
                debugsum=debugsum+depairdw(i4,1,i3)**2
                debugsum2=debugsum2+depairdw(i4,1,i3)
              enddo ! i4
              debugsum=dsqrt(debugsum)
              write(debugunit,'(i6,a,i6,x,3f14.8)')&
                countepoch,' depairdw length point ',point,errore,debugsum,debugsum2
            enddo ! i3
          endif ! ledbug
!!
!! add the damping term to the derivatives if weight damping is used
          if(ldampw) then
            do i2=1,npairs 
              do i3=1,num_weightspairfree(i2)
                depairdw(i3,1,i2)=(1.d0-dampw)*depairdw(i3,1,i2)&
                -dampw*2.0d0*weights_short_pair(wconstraintidxp(i3,i2),i2)/dble(totnum_weightspairfree)
              enddo ! i3
            enddo ! i2
          endif !ldampw
!!
!! adapt sign, this is needed in case of nenergygroup.gt.1 to avoid cancellation effects upon averaging
          if(errore.lt.0.0d0)then
            errore=-1.d0*errore
            depairdw(:,:,:)=-1.d0*depairdw(:,:,:)
          endif
!!
          do i2=1,num_pairs 
            num_pairs_element(pairindex(zelemp(1,i2),zelemp(2,i2)))=&
              num_pairs_element(pairindex(zelemp(1,i2),zelemp(2,i2)))+1
          enddo
!!
!! sum derivatives and errors if we group over several points
          depairdwsum(:,:,:)= depairdwsum(:,:,:)+depairdw(:,:,:)
          erroresum         = erroresum+errore
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! energy update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check if an update should be done
          if(ljointefupdate)goto 98
!!
 77       continue
!! check if there are energies left not used for update
          if((i1.eq.npoints).and.(nenergy.gt.0))then
            lenforceupdatee=.true.
          endif
!!
!! update if we have collected enough structures or there are no more structures left
          if(((mod(nenergy,nenergygroup).eq.0)&   ! a sufficient number of energies has been accumulated for update
            .and.(nenergy.gt.0)).or.&             ! CAUTION: this is needed, otherwise mod is also 0 for nenergy=0
            (lenforceupdatee))then                ! enforce update for the remaining energies not yet used
!!
            lenforceupdatee = .false.
!!
!! normalize the error and the derivatives in case grouping has been used
            errore          =erroresum/dble(nenergy)
            depairdw(:,:,:) =depairdwsum(:,:,:)/dble(nenergy)
!!
!! Kalman damping
            errore=errore*kalman_dampe
!!
!! chose optimization algorithm
            if(optmodee.eq.1)then ! KALMANN FILTER
!! update weights for each elemental pair
              do i2=1,npairs
!!
!! check if pair of this type are present, only then do update
                if(num_pairs_element(i2).gt.0)then
!!
!! prepare temporary arrays
                  allocate(weightsp(num_weightspairfree(i2)))
                  allocate(depairdw_temp(num_weightspairfree(i2)))
!!
!! reduce weights array to free weights
                  do i3=1,num_weightspairfree(i2)
                    weightsp(i3)     =weights_short_pair(wconstraintidxp(i3,i2),i2)
                    depairdw_temp(i3)=depairdw(i3,1,i2)
                  enddo ! i3
!!
                  if(mpisize.eq.1)then
!!
!! updatekalman_para cannot be used for 1 process because the corrmatrix array has a different dimension then
                    if(mpirank.eq.0)then
                      call updatekalman(kaldim(i2),corrdim(i2),&
                        num_weightspairfree(i2),kalmanlambdap(i2),kalmannue,&
                        weightsp,depairdw_temp,corrmatrix_list(1,i2),errore)
                    endif ! mpirank.eq.0
!!
                    call mpi_bcast(weightsp,num_weightspairfree(i2),&
                      mpi_real8,0,mpi_comm_world,mpierror)
                  else ! real parallel case
!!
!! updatekalman cannot be used for more than 1 process because the corrmatrix array has a different dimension then
                    call updatekalman_para(paramode,kaldim(i2),corrdim(i2),&
                      num_weightspairfree(i2),kalmanlambdap(i2),kalmannue,&
                      weightsp,depairdw_temp,corrmatrix_list(1,i2),errore)
                  endif ! mpisize.eq.1
!!
!! apply range restrictions to weights if requested
                  if(restrictw.gt.0.0d0)then
                    do i3=1,num_weightspairfree(i2)
                      if((weightsp(i3).gt.(-1.d0*restrictw+1.0d0))&
                        .and.(weightsp(i3).lt.(restrictw-1.0d0)))then
                      elseif(weightsp(i3).ge.(restrictw-1.0d0))then
                        weightsp(i3)=restrictw-1.0d0+tanh(weightsp(i3)-restrictw+1.0d0)
                      elseif(weightsp(i3).le.(-1.d0*restrictw+1.0d0))then
                        weightsp(i3)=-1.d0*restrictw+1.0d0+tanh(weightsp(i3)+restrictw-1.0d0)
                      endif
                    enddo ! i3
                  endif
!!
!! expand weights array back to original array
                  do i3=1,num_weightspairfree(i2)
                    weights_short_pair(wconstraintidxp(i3,i2),i2)=weightsp(i3)
                  enddo ! i3
!!
!! debugging output of weights_pair to debug.data
                  if((mpirank.eq.0).and.(pstring(1:1).eq.'1'))then
                    call debugweights(npairs,2,&
                      countepoch,point,i2,&
                      maxnum_weights_short_pair,&
                      maxnum_layers_short_pair,num_layers_short_pair,&
                      nodes_short_pair,weights_short_pair)
                  endif
!!
!! deallocate temporary arrays
                  deallocate(weightsp)
                  deallocate(depairdw_temp)
!!
                endif ! num_pairs_element(i2).gt.0
!!
              enddo ! i2=1,npairs
!!
            elseif(optmodee.eq.2)then ! conjugate gradient'
              write(ounit,*)'CG optimization not yet implemented in optimize_short'
              stop !'

            elseif(optmodee.eq.3)then ! steepest descent'
!! update weights for each element
              do i2=1,npairs
!!
!! check if pairs of this type are present, only then do update
                if(num_pairs_element(i2).gt.0)then ! pairs of element are present
!!
!! allocate temporary arrays for free weights
                  allocate(weightsp(num_weightspairfree(i2)))
                  allocate(depairdw_temp(num_weightspairfree(i2)))
!!
!! serial original
                  if(mpirank.eq.0)then
!!
!! reduce weights array to free weights
                    do i3=1,num_weightspairfree(i2)
                      weightsp(i3)=weights_short_pair(wconstraintidxp(i3,i2),i2)
                      depairdw_temp(i3)=depairdw(i3,1,i2)
                    enddo ! i3
!!
                    call updatesteepest(&
                      num_weightspairfree(i2),weightsp,&
                      depairdw_temp,errore,steepeststepe)
!!
!! apply restrictions to weights is requested
                    if(restrictw.gt.0.0d0)then
                      do i3=1,num_weightspairfree(i2)
                        if((weightsp(i3).gt.(-1.d0*restrictw+1.0d0))&
                          .and.(weightsp(i3).lt.(restrictw-1.0d0)))then
                        elseif(weightsp(i3).ge.(restrictw-1.0d0))then
                          weightsp(i3)=restrictw-1.0d0+tanh(weightsp(i3)-restrictw+1.0d0)
                        elseif(weightsp(i3).le.(-1.d0*restrictw+1.0d0))then
                          weightsp(i3)=-1.d0*restrictw+1.0d0+tanh(weightsp(i3)+restrictw-1.0d0)
                        endif
                      enddo ! i3
                    endif
!!
!! expand weights array back to original array
                    do i3=1,num_weightspairfree(i2)
                      weights_short_pair(wconstraintidxp(i3,i2),i2)=weightsp(i3)
                    enddo ! i3
!!
                  endif ! mpirank.eq.0
!!
                  call mpi_bcast(weights_short_pair,maxnum_weights_short_pair*npairs,&
                    mpi_real8,0,mpi_comm_world,mpierror)
!!
                  deallocate(weightsp)
                  deallocate(depairdw_temp)
!!
                endif ! num_pairs_element(i2).gt.0
!!
              enddo ! i2=1,npairs

            endif ! optmodee
!!
!!
!! reinitializations for next update
            nenergy                = 0
            num_pairs_element(:)   = 0
            depairdw(:,:,:)        = 0.0d0
            depairdwsum(:,:,:)     = 0.0d0
            errore                 = 0.0d0
            erroresum              = 0.0d0
!!
          endif !((mod(nenergy,nenergygroup).eq.0).or.(i1.eq.npoints))then
!!
        else
!! if joint_energy_force_update for E and F is used and E is not used, skip point completely(don't use forces)
          if(ljointefupdate)goto 101
!!
!! point is not used for short update
        endif !(abserror.gt.kalmanthreshold_temp*rmse_short)
!!
!! jump here if energy update of this point is skipped
 98     continue
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! force part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
        if(luseforces.and.(forcernd.gt.0.0d0))then
!!
!! automatically determine the number of forces to group for each element in this structure
          if(lfgroupbystruct) then
            totnumatoms=0
            do i2=1,nelem
              totnumatoms=totnumatoms+num_atoms_element_struct(i2)
            enddo
            mforcegroup=3*totnumatoms
          endif
!!
!! in case of joint_energy_force_update group all forces of this structure
          if(ljointefupdate) then
            totnumatoms=0
            do i2=1,nelem
              totnumatoms=totnumatoms+num_atoms_element_struct(i2)
            enddo
            mforcegroup=3*totnumatoms
          endif
!!
          do i2=1,nelem
            if(scalefactorf.lt.0.0d0)then
              scalefactorftemp(i2)=(dble(mforcegroup)*energyrnd)/&
                        (dble(3*num_atoms_element_struct(i2)*dble(nenergygroup))*forcernd)
            else ! keep scalefactorf from input.nn
              scalefactorftemp(i2)=scalefactorf
            endif
          enddo
!!
!! allocations for parallel use
          allocate(strs(3,3,maxnum_funcvalues_short_pair,max_num_pairs))
          allocate(dsfuncdxyz_mpi(maxnum_funcvalues_short_pair,max_num_pairs,max_num_atoms,3))
          dsfuncdxyz_mpi(:,:,:,:)=0.0d0
          allocate(symfunctiondummy(maxnum_funcvalues_short_pair,max_num_pairs))
!!
!! get dsfuncdxyz:
!! this is independent of the weights and needs to be done only once for each point i1
!! in principle dsfuncdxyz could be precalculated, but needs too much storage on disk
!! Caution: symmetry functions from this subroutine cannot be used because they are not scaled
!! parallel version
!! calculate the short range symmetry functions for pairs n_start to n_end
          if(mpisize.gt.1)then
            write(ounit,*)'ERROR: optimize_short_combined_pair is not yet parallel'
            stop
          endif
          call calconefunction_pair(n_start,&
            idx(i1),num_atoms,zelem,&
            num_pairs,pairs_charge,&
            lattice_list(1,1,idx(i1)),xyzstruct_list(1,1,idx(i1)),symfunctiondummy,&
            dsfuncdxyz_mpi,strs,&
            lperiodic,.true.,lrmin)
          if(.not.lrmin)then
            write(ounit,*)'Error: too close atoms in optimize_short_combined'
            stop !'
          endif
!'
!! scale dsfuncdxyz
          if(lscalesym)then
!! scale dsfuncdxyz and strs 
!! parallel version
            call scaledsfuncpair_para(num_pairs_para,&
              minvalue_short_pair,maxvalue_short_pair,zelemp,dsfuncdxyz_mpi,strs)
          endif ! lscalesym
!!
!! combine dsfuncdxyz of all processes
!! TODO: This should be better done with mpi_gatherv
          dsfuncdxyz_pair(:,:,:,:)=0.0d0
          icount=0
          do i2=n_start,n_end
            icount=icount+1
            dsfuncdxyz_pair(:,i2,:,:)=dsfuncdxyz_mpi(:,icount,:,:)
          enddo
          call mpi_allreduce(mpi_in_place,dsfuncdxyz_pair,&
            maxnum_funcvalues_short_pair*max_num_atoms*max_num_atoms*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!!
          deallocate(strs)
          deallocate(dsfuncdxyz_mpi)
          deallocate(symfunctiondummy)
!!
!! calculate the error of the forces individually for each atom and component
          do i2=1,num_atoms ! Lopp over all the atom in a structure
            do i3=1,3 ! loop over x,y and z
!!
!! counter for forces (counts all forces no matter if used for updating or not)
              count_atoms_element_block(elementindex(zelem(i2)))=&
                count_atoms_element_block(elementindex(zelem(i2)))+1
!!
              if(luseworstf.and.(.not.lusef(i3,i2,idx(i1))))then
                goto 99 ! skip this point
              endif
!!
!! if only a random subset of points is used: decide if we use this point
              z=ran0(kseed)
              if(z.gt.forcernd)then ! don't use force for weight update
                goto 99 ! skip this point
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
              call getoneshortforcepair_para(i2,i3,&
                pindex,num_pairs,zelemp,&
                symfunction_short_pair_list(1,1,idx(i1)),&
                dsfuncdxyz_pair,nnshortforce)
!!
!! calculate the error for the current force component
               errorf=-1.d0*(shortforce(i3,i2)-nnshortforce(i3,i2)) 
!!             
!! transform units to a per pair basis
              errorf=errorf/dble(num_pairs)   ! error per pairs
              abserrorf=abs(errorf)
!!
!! skip force update if error is within noise
              if(abserrorf.lt.noisef)then
                goto 99 ! skip this point
              endif
!!
              if(ldampw)then
                write(ounit,*)'ERROR: weight damping not implemented for forces'
                stop !'
              endif
!!
              errorf=errorf*scalefactorftemp(elementindex(zelem(i2)))
!!
!! set errorthreshold criterion for short range energy
              if(lfixederrorf)then
                errorthresholdf=fixederrorf/dble(num_pairs)                         !! errorthreshold per pair
              else
                errorthresholdf=kalmanthresholdf_temp*(rmse_force_s_ref)/dble(num_pairs) !! errorthreshold per pair

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
                (lfgroupbystruct.and.(count_forces.eq.3*num_atoms))&
                  )then                                  ! last force of this element in this structure, update may be necessary
!!
!! do some counting
                fitstatf(i3,i2,ndone+idx(i1))=fitstatf(i3,i2,ndone+idx(i1))+1
                count_forces=count_forces+1
                numberf = numberf+1

                do i4=1,num_pairs 
                  num_pairs_element_f(pairindex(zelemp(1,i4),zelemp(2,i4)))=&
                    num_pairs_element_f(pairindex(zelemp(1,i4),zelemp(2,i4)))+1
                enddo
!!
!! calculate the derivative of the short range forces with respect to the weights
!! dfpairdw is for one atom only, but depends on all atoms
!! calculate dfpairdw
!!
!! serial version
                if(paramode.eq.1)then
                  if(mpirank.eq.0)then
                    call getdfpairdw(i2,i3,&
                      num_weightspairfree,zelemp,num_pairs,&
                      wconstraintidxp,&
                      symfunction_short_pair_list(1,1,idx(i1)),dsfuncdxyz_pair,&
                      dfpairdw)
                  endif ! mpirank.eq.0
!!
!! CHECK: Do we have to normalize dfshortdw per atom here to be consistent with atomic energies?
                   dfpairdw(:,:)=dfpairdw(:,:)/dble(num_pairs) ! Jovan 
!!                   dfpairdw(:,:)=dfpairdw(:,:)/dble(num_atoms_element_struct(elementindex(i2))) ! JB test => NaN 
!!                   dfpairdw(:,:)=dfpairdw(:,:)/(dble(num_atoms_element_struct(elementindex(i2)))*dble(num_pairs)) ! JB test => NaN 
!!
                  call mpi_bcast(dfpairdw,maxnum_weights_short_pair*npairs,&
                    mpi_real8,0,mpi_comm_world,mpierror)
!!
               elseif(paramode.eq.2)then
!! parallel version
                 call getdfpairdw_para(i2,i3,&
                   num_weightspairfree,num_pairs_para,pindex,zelemp,num_pairs,&
                   wconstraintidxp,symfunction_short_pair_list(1,1,idx(i1)),dsfuncdxyz_pair,&
                   dfpairdw)
!!
!! CHECK: Do we have to normalize dfshortdw per atom here to be consistent with atomic energies?
                 dfpairdw(:,:)=dfpairdw(:,:)/dble(num_pairs)    !! NEW NEW
!!
               else
                 write(ounit,*)'Error: unsupported paramode'
                 stop
               endif ! paramode
!!
!! scaling for each pair is done here by a factor determined by the element! This should be correct!  ??????
               dfpairdw(:,:)=dfpairdw(:,:)*scalefactorftemp(elementindex(zelem(i2)))
!!
!! adapt sign, this is needed in case of mforcegroup .gt. 1 to avoid cancellations upon averaging  NEWLYADDED
               if(errorf.lt.0.0d0)then
                 errorf=-1.d0*errorf
                 dfpairdw(:,:)=-1.d0*dfpairdw(:,:)         !!CHECK
               endif
!!
!! if damping is requested:
!! CAUTION in pair case!!!
!!                  if(ldampw) then
!!                    do i4=1,num_weightspairfree(pairindex(zelemp(1,i2),zelemp(2,i2)))
!!                      dfpairdw(i4)=(1.d0-dampw)*dfpairdw(i4)-dampw*2.0d0*weights_pair(wconstraintidxp(i4,&
!!                            pairindex(zelemp(1,i2),zelemp(2,i2))),pairindex(zelemp(1,i2),zelemp(2,i2)))/&
!!                            dble(num_weightspairfree(pairindex(zelemp(1,i2),zelemp(2,i2))))
!!                    enddo ! i4'
!!                  endif
!!
!! accumulate dfpairdw and errorf for grouping
               dfpairdwsum(:,:)=dfpairdwsum(:,:)+dfpairdw(:,:)  ! for each pair NN!
               errorfsum=errorfsum+errorf
               nforce=nforce+1
!!
!! In case of joint_energy_force_update skip individual force updates
               if(ljointefupdate)goto 99
!!
 88            continue
!!
!! check if there are forces left not used for update (standard grouping case)
               if((i1.eq.npoints).and.(i2.eq.num_atoms).and.(i3.eq.3)&  ! is this the last force in block of points?
                 .and.(count_forces.gt.0))then                          ! there are some so far unused forces
                 lenforceupdatef=.true.
               endif
!!
!! check if there are forces left not used for update (group_forces_by_structure case)
               if((lfgroupbystruct).and.(i2.eq.num_atoms).and.(i3.eq.3))then ! is this the last force in this structure
                 if(count_forces.gt.0)then
                   lenforceupdatef=.true.
                 endif
               endif
!!
!! do force weight update for each pair in the present structure 
!!
               if(((mod(count_forces,mforcegroup).eq.0).or.&     ! number of forces for grouping has been reached for this element
                 (lenforceupdatef))&                             ! last force for this element in thisblock of structures
                 .and.(count_forces.gt.0))then                   ! CAUTION: mod is also 0 for count_forces=0
!!
                 do i4=1,npairs
!!
                   allocate(weightsp(num_weightspairfree(i4)))
                   allocate(dfpairdw_temp(num_weightspairfree(i4)))
                   allocate(depairdw_temp(num_weightspairfree(i4)))
!!
!! Do update for this element i4 only if pairs are present
                   if(num_pairs_element_f(i4).gt.0)then
!!
!! debug
                     if((mpirank.eq.0).and.(pstring(1:1).eq.'1')) &
                       write(debugunit,'(a,4i8,f14.6)')'force update done ',i2,i3,i4,num_pairs_element_f(i4),errorf
!!'
!! normalize by the number of atoms in this element to reduce the importance of the forces with respect to the energies, should we do this?
                     dfpairdw(:,i4)=dfpairdwsum(:,i4)/(dble(nforce)) 
                     errorf=errorfsum/(dble(nforce)) 
!!
!! Kalman damping
                     errorf=errorf*kalman_dampf
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the weights:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
                     do i5=1,num_weightspairfree(i4)
                       weightsp(i5)=weights_short_pair(wconstraintidxp(i5,i4),i4)
                       dfpairdw_temp(i5)=dfpairdw(i5,i4)
!                      if(i2.eq.2.and.i3.eq.2)
!                       write(*,*)'DFPDW',i1,i2,i3,i5,dfpairdw_temp(i5)
                     enddo ! i5
!!
                     if(optmodef.eq.1)then !                                    KALMAN FILTER
!!
                       if(.not.lsepkalman)then
                         corrmatrixf_list(:,i4)=corrmatrix_list(:,i4)
                       endif
!!
                       if(mpisize.eq.1)then
                         call updatekalman(kaldim(i4),corrfdim(i4),&
                           num_weightspairfree(i4),kalmanlambdap(i4),kalmannue,&
                           weightsp,dfpairdw_temp,corrmatrixf_list(1,i4),errorf)

                       else
                         call updatekalman_para(paramode,kaldim(i4),corrfdim(i4),&
                           num_weightspairfree(i4),kalmanlambdap(i4),kalmannue,&
                           weightsp,dfpairdw_temp,corrmatrixf_list(1,i4),errorf)
                       endif ! mpisize.eq.1
!!
                       if(.not.lsepkalman)then
                         corrmatrix_list(:,i4)=corrmatrixf_list(:,i4)
                       endif
!!
!! Update again the weights using the energy after each individual force update if requested:
                       if(lrepeate.and.(shortenergy_list(idx(i1)).lt.maxenergy))then
                         nnpairenergy(:)= 0.0d0
                         nntotalenergy  = 0.0d0
                         do i5=1,num_weightspairfree(i4)
                           weights_short_pair(wconstraintidxp(i5,i4),i4)=weightsp(i5)
                         enddo ! i5
!!
                         call calconeshort_parapair(point,num_pairs_para,pindex,&
                           zelemp,symfunction_short_pair_list(1,n_start,idx(i1)),nnpairenergy,&
                           nntotalenergy)
!!
                         call mpi_allreduce(mpi_in_place,nnpairenergy,&
                           max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)
                           nntotalenergy=0.0d0
                         do i5=1,num_pairs 
                           nntotalenergy=nntotalenergy+nnpairenergy(i5)
                         enddo
                         nneshort=nntotalenergy/dble(num_atoms)
                         errore=(shortenergy_list(idx(i1))-nneshort)*dble(num_atoms)/dble(num_pairs)
                         abserrore=abs(errore)  
                         errore   =errore*kalman_dampe          
!!
                         call getdepairdw(&
                           num_weightspairfree,&
                           zelemp,num_pairs,wconstraintidxp,&
                           symfunction_short_pair_list(1,1,idx(i1)),depairdw)
!!
                         if(errore.lt.0.0d0)then
                           errore=-1.d0*errore
                           depairdw(:,:,:)=-1.d0*depairdw(:,:,:)
                         endif
                         do i5=1,num_weightspairfree(i4)
                           depairdw_temp(i5)=depairdw(i5,1,i4)
                         enddo ! i5
!!
                         if(mpisize.eq.1)then
                           call updatekalman(kaldim(i4),corrdim(i4),&
                             num_weightspairfree(i4),kalmanlambdap(i4),kalmannue,&
                             weightsp,depairdw_temp,corrmatrix_list(1,i4),&
                             errorf)
                         else ! real parallel case
                           call updatekalman_para(paramode,kaldim(i4),corrdim(i4),&
                             num_weightspairfree(i4),kalmanlambdap(i4),kalmannue,&
                             weightsp,depairdw_temp,corrmatrix_list(1,i4),&
                             errorf)
                         endif ! mpisize.eq.1
                         numbere=numbere+1  
                       endif ! lrepeate
!!
                     elseif(optmodef.eq.2)then !                                             conjugate gradient'
!!
                     elseif(optmodef.eq.3)then !                                              steepest descent'
!!                                                                                                                       
                       if(mpirank.eq.0)then                                                                              
!!
                         weightsp(:)   =weights_short_pair(:,pairindex(zelemp(1,i2),zelemp(2,i2)))
!!                                                                                                                       
                         call updatesteepest(num_weightspairfree(i4),&                                             
                           weightsp,dfpairdw_temp,errorf,steepeststepf)                                           
!!                                                                                                                       
                       endif ! mpirank.eq.0                                                                              

                       call mpi_bcast(weightsp,maxnum_weights_short_pair,&
                         mpi_real8,0,mpi_comm_world,mpierror)      
!!                                                                 
!! Update again the weights using the energy after each individual force update if requested:
                       if(lrepeate.and.(shortenergy_list(idx(i1)).lt.maxenergy))then         
                         nnpairenergy(:)= 0.0d0                                              
                         nntotalenergy  = 0.0d0                                              
                         do i5=1,num_weightspairfree(i4)                                     
                           weights_short_pair(wconstraintidxp(i5,i4),i4)=weightsp(i5)              
                         enddo ! i5                                                          
!!                                                                                           
                         call calconeshort_parapair(point,num_pairs_para,pindex,&          
                           zelemp,symfunction_short_pair_list(1,n_start,idx(i1)),nnpairenergy,& 
                           nntotalenergy)                    
!!                                                                                           
                         call mpi_allreduce(mpi_in_place,nnpairenergy,&                      
                           max_num_atoms,mpi_real8,mpi_sum,mpi_comm_world,mpierror)          
!!                                                                                           
                           nntotalenergy=0.0d0                                               
                         do i5=1,num_pairs ! total energy = sum of pair-energies             
                           nntotalenergy=nntotalenergy+nnpairenergy(i5)                      
                         enddo                                                               
!!                                                                                           
                         nneshort=nntotalenergy/dble(num_atoms)                                
                         errore=(shortenergy_list(idx(i1))-nneshort)*dble(num_atoms)/dble(num_pairs)
                         abserrore=abs(errore)                                                    
                         errore   =errore*kalman_dampe                                            
!!                                                                                                
                         call getdepairdw(&                                                 
                           num_weightspairfree,&                        
                           zelemp,num_pairs,wconstraintidxp,&                                        
                           symfunction_short_pair_list(1,1,idx(i1)),depairdw)                                        
!!                                                                                                
                         if(errore.lt.0.0d0)then                                                  
                           errore=-1.d0*errore                                                    
                           depairdw(:,:,:)=-1.d0*depairdw(:,:,:)                                  
                         endif                                                                    
!!                                                                                                
                         do i5=1,num_weightspairfree(i4)                                          
                           depairdw_temp(i5)=depairdw(i5,1,i4)                                   
                         enddo ! i5                                                               
!!                                                                                                
                         call updatesteepest(&                                              
                           num_weightspairfree(i4),weightsp,&                                     
                           depairdw_temp,errore,steepeststepe)                             
                           numbere=numbere+1  !!1 NEWLY ADDED                                     
                       endif ! lrepeate                   
!!
                     endif ! optmodef                                                      end optmodef
!!
!!
!! apply restrictions to weights is requested
                     if(restrictw.gt.0.0d0)then
                       do i5=1,num_weightspairfree(i4)
                         if((weightsp(i5).gt.(-1.d0*restrictw+1.0d0))&
                           .and.(weightsp(i5).lt.(restrictw-1.0d0)))then
                         elseif(weightsp(i5).ge.(restrictw-1.0d0))then
                           weightsp(i5)=restrictw-1.0d0+tanh(weightsp(i5)-restrictw+1.0d0)
                         elseif(weightsp(i5).le.(-1.d0*restrictw+1.0d0))then
                           weightsp(i5)=-1.d0*restrictw+1.0d0+tanh(weightsp(i5)+restrictw-1.0d0)
                         endif
                       enddo ! i5
                     endif
!!
!! expand weights array back
                     do i5=1,num_weightspairfree(i4)
                       weights_short_pair(wconstraintidxp(i5,i4),i4)=weightsp(i5)
                     enddo ! i5
!!
                     call mpi_bcast(weights_short_pair,maxnum_weights_short_pair*npairs,&
                       mpi_real8,0,mpi_comm_world,mpierror)
!!
                   endif !(num_pairs_element_f(i4).gt.0)then
!!
                   deallocate(weightsp)
                   deallocate(dfpairdw_temp)
                   deallocate(depairdw_temp) 
!!
                 enddo ! i4 loop over all pairs
!!
                 nforce                = 0
                 dfpairdw(:,:)         = 0.0d0
                 dfpairdwsum(:,:)      = 0.0d0
                 errorf                = 0.0d0
                 errorfsum             = 0.0d0
                 lenforceupdatef       =.false.
                 count_forces          = 0
                 num_pairs_element_f(:)= 0
!!
               endif !(((mod(count_forces,mforcegroup).eq.0)
!!
!!
!! reinitializations for next update
!!
           else ! point is not used for short update
!!
           endif !(abserrorf.gt.kalmanthresholdf_temp*rmse_force_s_ref)
!!
!! jump here if force should be skipped for weight update
 99        continue
!!
          enddo ! i3 loop over x,y and z
!!
         enddo ! i2 loop over all atoms
!!
       endif ! luseforces
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end of regular update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!
!! do joint update of energy and forces here if requested          
       if(ljointefupdate)then
         write(ounit,*)'ERROR: jointefupdate not implemented in optimize_short_combinedpair'
         stop !'
       endif ! ljointefupdate
!!
 101   continue
       deallocate(pindex)
!!
!!
      enddo ! i1=1,npoints                           
!!
      return
      end
