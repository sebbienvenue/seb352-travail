!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!! - fittingpair.f90
!!
!! This routine is called once for each block of points
!!
      subroutine optimize_ewald(npoints,pointe,idx,&
         countepoch,ntrain,ndone,&
         maxcorredim,kaledim,kalcdim,corredim,corrcdim,&
         num_weightsewaldfree,fitstatq,&
         wconstraintidxe,mseed,numberq,&
         kalmanthresholde_temp,&
         rmse_charge,rmse_totalcharge,&
         corrmatrixe_list,corrmatrixc)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use globaloptions
      use nnewald
      use structures
      use timings
!!
      implicit none
!!
      integer countepoch                                 ! in
      integer wconstraintidxe(maxnum_weights_elec,nelem) ! in
      integer i1,i2,i3,i4                                ! internal
      integer icount,jcount                              ! internal
      integer idx(nblock)                                ! in
      integer itemp                                      ! internal
      integer mseed                                      ! in/out
      integer ncharge(nelem)                             ! internal
      integer mchargegroup(nelem)                        ! internal
      integer npoints                                    ! in
      integer ndone                                      ! in
      integer num_atoms                                  ! internal
      integer num_atoms_element(nelem)                   ! internal 
      integer num_atoms_element_block(nelem)             ! internal 
      integer num_weightsewaldfree(nelem)                ! in
      integer numberq                                    ! in/out
      integer pointe                                     ! in/out
      integer fitstatq(max_num_atoms,ntrain)             ! in/out
      integer ntrain                                     ! in
      integer ncount(nelem)                              ! internal
!! Kalman matrix dimensions:
      integer maxcorredim                                ! in
      integer corredim(nelem)                            ! in
      integer corrcdim                                   ! in
      integer kaledim(nelem)                             ! in
      integer kalcdim                                    ! in
      integer constraintdim                              ! internal
!!
      real*8 kalmanthresholde_temp                                               ! in
      real*8 rmse_charge                                                    ! in
      real*8 rmse_totalcharge                                               ! in
      real*8 error                                                          ! internal
      real*8 errorsum(nelem)                                                ! internal
      real*8 abserror                                                       ! internal
      real*8 abschargeerror                                                 ! internal
      real*8 chargeerror                                                    ! internal
      real*8 nnchargesum                                                    ! internal
      real*8, dimension(:), allocatable :: weights                          ! internal
!! CAUTION: just one output node is assumed here
      real*8 dqdw(maxnum_weights_elec,1)                                    ! internal
!! CAUTION: just one output node is assumed here
      real*8 dqdwc(maxnum_weights_elec,1)                                   ! internal
      real*8, dimension(:), allocatable :: dqdwc_temp                       ! internal
      real*8 corrmatrixe_list(maxcorredim,nelem)                            ! in/out
      real*8 nnatomcharge                                                   ! internal
      real*8 nnatomcharge_list(nblock,max_num_atoms)                        ! internal 
      real*8 chargeerror_list(nblock,max_num_atoms)                         ! internal 
      real*8 nodes_values(maxnum_layers_elec,maxnodes_elec)                ! internal, just dummy 
      real*8 nodes_sum(maxnum_layers_elec,maxnodes_elec)                   ! internal, just dummy
      real*8 corrmatrixc(corrcdim)                                          ! in/out
      real*8, dimension(:), allocatable :: dqdwconstraint                   ! internal
      real*8, dimension(:), allocatable :: weightsconstraint                ! internal
!! CAUTION: just one output node is assumed here
      real*8 dqdwsum(maxnum_weights_elec,1,nelem)                           ! internal
!! CAUTION: just one output node is assumed here
      real*8 dqdwsumc(maxnum_weights_elec,1,nelem)                          ! internal
      real*8 z,ran0                                                         ! internal
      real*8 sumwsquared                                                    ! internal
!!
      logical luseq(nblock,max_num_atoms)                                   ! internal 
      logical ldebug2                                                       ! internal
      logical lenforceupdate                                                ! internal
!!
!! initialization
      errorsum(:)         = 0.0d0
      dqdwsum(:,:,:)      = 0.0d0
      dqdwsumc(:,:,:)     = 0.0d0
      ncharge(:)          = 0
      luseq(:,:)          = .false.
      ldebug2             = .false. 
      lenforceupdate      = .false.
      ncount(:)           = 0
!! timing variables
      dayqerror           = 0
      daydqdw             = 0
      dayqupdate          = 0
!!
!! count the total number of atoms of each element for the structures in this block of points, no need to use idx here
      num_atoms_element_block(:)=0
      do i1=1,npoints
        do i2=1,num_atoms_list(i1)
          num_atoms_element_block(elementindex(zelem_list(i1,i2)))=&
            num_atoms_element_block(elementindex(zelem_list(i1,i2)))+1
        enddo
      enddo
!!
!! determine mchargegroup individually for each element, make sure mchargegroup is not larger than the number of atoms of this element
      mchargegroup(:)=nchargegroup
      do i1=1,nelem
        if(mchargegroup(i1).gt.num_atoms_element_block(i1))then
          mchargegroup(i1)=num_atoms_element_block(i1)
          if(mchargegroup(i1).eq.0)then ! in case element is not present
            mchargegroup(i1)=1
          endif
        endif
      enddo
!! 
      if(ldebug2) write(ounit,*)mpirank,' optimize_ewald starts '
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if only the worst points shall be used for the update, do the preparations for the block of points here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(luseworstq)then
        if(worstq.ge.1.0d0)then
          luseq(:,:)=.true.
        else
          call getcharges(nblock,npoints,&
            zelem_list,num_atoms_list,&
            symfunction_elec_list,nnatomcharge_list)
!! calculate the error of all charges
!! Caution: in case of lupdatebyelement only the nnatomcharges for that element are ok
          chargeerror_list(:,:)=0.0d0
          do i1=1,npoints
            do i2=1,num_atoms_list(i1)
              chargeerror_list(i1,i2)&
                =abs(nnatomcharge_list(i1,i2)-atomcharge_list(i1,i2))
            enddo ! i2
          enddo ! i1
!! sort points by error and determine error array luseq 
          call sortchargeerror(npoints,&
            chargeerror_list,luseq)
        endif ! worstq
      endif ! luseworstq
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i1=1,npoints
!!
        pointe              =pointe+1                 ! counter of points
        num_atoms           =num_atoms_list(idx(i1))
!!
!! count the total number of atoms of each element in this structure
        num_atoms_element(:)=0 
        do i2=1,num_atoms_list(idx(i1))
          num_atoms_element(elementindex(zelem_list(idx(i1),i2)))&
            =num_atoms_element(elementindex(zelem_list(idx(i1),i2)))+1
        enddo ! i2
!!
!! automatically determine number of atoms to group for this element if requested
        if(lqgroupbystruct)then
          mchargegroup(:)=num_atoms_element(:)
          do i2=1,nelem
            if(mchargegroup(i2).eq.0)then ! in case element is not present
              mchargegroup(i2)=1
            endif
          enddo
        endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all atoms i2 of point idx(i1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i2=1,num_atoms
!!
          itemp=elementindex(zelem_list(idx(i1),i2))               ! abbreviation
!!
!! count atoms of each element done so far in this block (over all points, not just this structure!)
          ncount(itemp)=ncount(itemp)+1
!!
          if(ldebug2) write(ounit,*)& 
            mpirank,' optimize_ewald point ',i1,idx(i1),zelem_list(idx(i1),i2)
!!
!! skip this charge, if it is not one of the worst charges
          if(luseworstq.and.(.not.luseq(idx(i1),i2)))then 
            goto 99 
          endif     
!!
!! skip this charge, if only one element is to be fitted and this is not the correct one 
          if(lupdatebyelement.and.(zelem_list(idx(i1),i2).ne.elemupdate))then
            goto 99 
          endif
!!
!! if only a random subset of charges is used: decide if we use this charge 
          z=ran0(mseed)
          if(z.gt.chargernd)then
            goto 99 
          endif 
!!
!! predict the NN charge of this atom i2, this is not parallel
          call abstime(timeqerrorstart,dayqerror)
          if(mpirank.eq.0)then
            call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
              maxnum_layers_elec,num_layers_elec(itemp),&
              maxnum_weights_elec,nodes_elec(0,itemp),&
              symfunction_elec_list(1,i2,idx(i1)),&
              weights_elec(1,itemp),&
              nodes_values,nodes_sum,&
              nnatomcharge,actfunc_elec(1,1,itemp))
!!         
            error   =    atomcharge_list(idx(i1),i2)-nnatomcharge
            abserror=abs(error)
            call abstime(timeqerrorend,dayqerror)
            timeqerror=timeqerror+(timeqerrorend-timeqerrorstart)
!!
!! if weight damping is requested, add squared sum of weights to error
            if(ldampw)then
              sumwsquared=0.0d0
              do i3=1,num_weights_elec(itemp)
                sumwsquared=sumwsquared &
                  +(weights_elec(i3,itemp))**2.0d0         
              enddo
              error=(1.d0-dampw)*error &
                + (dampw*sumwsquared)/dble(num_weights_elec(itemp))
            endif ! ldampw
!!
            if(ldebug2) write(ounit,*)'abserror ',i1,idx(i1),zelem_list(idx(i1),i2),abserror
          endif ! mpirank.eq.0
!!
!! distribute abserror and error to all processes
          call mpi_bcast(abserror,1,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(error,1,mpi_real8,0,mpi_comm_world,mpierror)
!!
!! skip charge if error is within noise (should not be done within mpirank.eq.0 if statement)
          if(abserror.lt.noiseq)then
            goto 99
          endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! check if error is sufficiently large to justify weight update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          if((abserror.gt.kalmanthresholde_temp*rmse_charge)&
             .or.&
             (ncount(itemp).eq.num_atoms_element_block(itemp))&
             .or.&
             (lqgroupbystruct.and.(i2.eq.num_atoms)))then ! enable update also for last charge of this element
!!
            if(abserror.le.kalmanthresholde_temp*rmse_charge) goto 77 ! jump over gradient directly to update (this charge is too good)
!!
!!            write(ounit,*)'Using charge ',i1,i2,itemp
!! count this charge
            fitstatq(i2,ndone+idx(i1))=fitstatq(i2,ndone+idx(i1))+1
            ncharge(itemp)=ncharge(itemp)+1
            numberq = numberq+1
!!
!! calculate the derivate dqdw with respect to the weights
            dqdw(:,:) =0.0d0 ! initialization of derivative array for this atom
            dqdwc(:,:)=0.0d0 ! initialization for this atom
!!
            if(mpirank.eq.0)then
!! get dqdw, not parallel
              call abstime(timedqdwstart,daydqdw)
              call getonededw(1,&
                maxnum_funcvalues_elec,maxnum_weights_elec,&
                maxnodes_elec,maxnum_layers_elec,&
                num_layers_elec(itemp),&
                windex_elec(1,itemp),&
                nodes_elec(0,itemp),&
                symfunction_elec_list(1,i2,idx(i1)),&
                weights_elec(1,itemp),&
                dqdw,actfunc_elec(1,1,itemp))
              call abstime(timedqdwend,daydqdw)
              timedqdw=timedqdw+(timedqdwend-timedqdwstart)
!!
!! adapt sign, this is needed in case of nchargegroup.gt.1 to avoid cancellation effects upon averaging
              if(error.lt.0.0d0)then
                error=-1.d0*error
                dqdw(:,:)=-1.d0*dqdw(:,:)
              endif
!!
!! compress array dqdw to array dqdwc by skipping all fixed weights
              do i3=1,num_weightsewaldfree(itemp)
                dqdwc(i3,1)=dqdw(wconstraintidxe(i3,itemp),1)
              enddo
            endif ! mpirank.eq.0
!!
            call mpi_bcast(dqdwc,maxnum_weights_elec*1,&
              mpi_real8,0,mpi_comm_world,mpierror)
!!
!! damp the charge weights in dqdwc
            if(ldampw)then
              do i3=1,num_weightsewaldfree(itemp)
                dqdwc(i3,1)=(1.d0-dampw)*dqdwc(i3,1) &
                  -dampw*2.0d0*weights_elec(wconstraintidxe(i3,i2),itemp)&
                  /dble(num_weightsewaldfree(itemp))
              enddo ! i3
            endif
!!
!! sum up dqdwc and error in case we want to group charges for updating
            dqdwsumc(:,:,itemp)=dqdwsumc(:,:,itemp)+dqdwc(:,:)
            errorsum(itemp)=errorsum(itemp)+error
!!
  77        continue
!! check if we have arrived at the final atom of an element in this block of points and if some
!! derivatives have not been used for updating 
!! => enforce a final update for this element
            if((ncharge(itemp).gt.0).and.&                           ! => some derivatives are left
              (ncount(itemp).eq.num_atoms_element_block(itemp)))then ! => present atom is last atom of this element in this block
               lenforceupdate=.true.         ! do a final update with the remaining derivatives for this element
!!               write(ounit,*)'setting lenforceupdate'
            endif
!!
!! check that in case of group_charges_by_structure if some
!! derivatives have not been used for updating (e.g. due to charge_fraction or adaptive Kalman filter some charges may be missing) 
            if(lqgroupbystruct.and.&           ! we MUST do one update for each element for each structure
              (ncharge(itemp).gt.0).and.&     ! => some derivatives are left
              (i2.eq.num_atoms))then ! => present atom is last atom of this structure 
               lenforceupdate=.true.         ! do a final update with the remaining derivatives for this element
!!               write(ounit,*)'setting lenforceupdate'
            endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! decide if weight update should be done here now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call abstime(timequpdatestart,dayqupdate)
            if(((mod(ncharge(itemp),mchargegroup(itemp)).eq.0).or.&    ! enough charges have been accumulated for this element
              (lenforceupdate)).and.&                                  ! update using remaining charges for this element
               (ncharge(itemp).gt.0))then                                ! CAUTION: mod is also 0 for ncharge=0
!!
!!              write(ounit,*)'Charge update done ',i1,i2,itemp
!!
              lenforceupdate=.false.
!!
!! renormalize the error and the derivatives here (division by the number of accumulated charges of this element)
              dqdwc(:,:)=dqdwsumc(:,:,itemp)/dble(ncharge(itemp))
              error     =errorsum(itemp)/dble(ncharge(itemp))
!!
!! reduce arrays to temporary arrays including just free weights
              allocate(weights(num_weightsewaldfree(itemp)))
              allocate(dqdwc_temp(num_weightsewaldfree(itemp)))
              do i3=1,num_weightsewaldfree(itemp)
                weights(i3)    = weights_elec(wconstraintidxe(i3,itemp),itemp)
                dqdwc_temp(i3) = dqdwc(i3,1)
              enddo ! i3
!!
              error=error*kalman_dampq
!!
              if(optmodeq.eq.1)then ! Kalman filter
!!
                if((mpisize.eq.1).or.lompmkl)then
                  call updatekalman(&
                    kaledim(itemp),corredim(itemp),&
                    num_weightsewaldfree(itemp),&
                    kalmanlambdae(itemp),kalmannuee,&
                    weights,dqdwc_temp,&
                    corrmatrixe_list(1,itemp),error)
                else
                  call updatekalman_para(paramode,&
                    kaledim(itemp),corredim(itemp),&
                    num_weightsewaldfree(itemp),&
                    kalmanlambdae(itemp),kalmannuee,&
                    weights,dqdwc_temp,&
                    corrmatrixe_list(1,itemp),error)
                endif
!!
              elseif(optmodeq.eq.2)then ! conjugate gradient'
                write(ounit,*)'Error: CG optimization not yet implemented in optimize_ewald'
                stop !'
!!
              elseif(optmodeq.eq.3)then ! steepest descent'
!!
                call updatesteepest(&
                  num_weightsewaldfree(itemp),&
                  weights,&
                  dqdwc_temp,error,steepeststepq)
!!
              endif ! optmodeq
!!
!! expand weights array back to original array
              do i3=1,num_weightsewaldfree(itemp)
                weights_elec(wconstraintidxe(i3,itemp),itemp)=weights(i3)
              enddo ! i3
!!
!! debugging output of weights_short
              if((mpirank.eq.0).and.(pstring(2:2).eq.'1'))then
                call debugweights(nelem,0,&
                  countepoch,pointe,itemp,&
                  maxnum_weights_elec,&
                  maxnum_layers_elec,num_layers_elec,&
                  nodes_elec,weights_elec)
              endif
!!
!! deallocate temporary arrays
              deallocate(weights)
              deallocate(dqdwc_temp)
!!
              call mpi_barrier(mpi_comm_world,mpierror)
!!
!! reinitialize counters and grouping arrays after update
              ncharge(itemp)      = 0
              dqdwsumc(:,:,itemp) = 0.0d0
              errorsum(itemp)     = 0.0d0
!!
            endif !((mod(ncharge(itemp),mchargegroup(itemp)).eq.0)
            call abstime(timequpdateend,dayqupdate)
            timequpdate=timequpdate+(timequpdateend-timequpdatestart)
!!
          else
!!              write(ounit,*)mpirank,' no update done ',i1,i2,i3
          endif ! (abserror.gt.kalmanthresholde_temp*rmse_charge)
!!
          call mpi_barrier(mpi_comm_world,mpierror)
!!
!! jump here, if charge is not used for updating
 99       continue
!!
        enddo ! i2, loop over all atoms
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! from here on nothing is parallelized 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! CHECK if the constraint here still works!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! apply constraint on total charges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(lchargeconstraint)then
!!
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### using the charge constraint is not'
          write(ounit,*)'implemented efficiently'
          write(ounit,*)'implementation is very old and not checked'
          stop
        endif
!!'
        constraintdim=0
        do i2=1,nelem
          constraintdim=constraintdim+num_weightsewaldfree(i2)
        enddo ! i2
        allocate(dqdwconstraint(constraintdim))
        allocate(weightsconstraint(constraintdim))
!!
!! calculate all atomic charges/derivatives dqdw for the final set of weights
        nnchargesum=0.0d0
        dqdwconstraint(:)=0.0d0 ! initialization
        do i2=1,num_atoms
!!          weights(:)     =weights_ewald(:,elementindex(zelem_list(idx(i1),i2)))
!!
!! predict the NN charge of that atom
          call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
            maxnum_layers_elec,num_layers_elec(elementindex(zelem_list(idx(i1),i2))),&
            maxnum_weights_elec,nodes_elec(0,elementindex(zelem_list(idx(i1),i2))),&
            symfunction_elec_list(1,i2,idx(i1)),&
            weights_elec(1,elementindex(zelem_list(idx(i1),i2))),&
            nodes_values,nodes_sum,&
            nnatomcharge,actfunc_elec(1,1,elementindex(zelem_list(idx(i1),i2))))
!!
          nnchargesum=nnchargesum+nnatomcharge ! sum the NN charges for constraint
!!
          dqdw(:,:)=0.0d0
          call getonededw(1,&
            maxnum_funcvalues_elec,maxnum_weights_elec,&
            maxnodes_elec,maxnum_layers_elec,num_layers_elec(elementindex(zelem_list(idx(i1),i2))),&
            windex_elec(1,elementindex(zelem_list(idx(i1),i2))),&
            nodes_elec(0,elementindex(zelem_list(idx(i1),i2))),&
            symfunction_elec_list(1,i2,idx(i1)),&
            weights_elec(1,elementindex(zelem_list(idx(i1),i2))),&
            dqdw,&
            actfunc_elec(1,1,elementindex(zelem_list(idx(i1),i2))))
!!
!! sum derivatives for charge constraint
          icount=(elementindex(zelem_list(idx(i1),i2))-1)*num_weights_elec(elementindex(zelem_list(idx(i1),i2)))+1
          jcount=1
          do i3=1,num_weights_elec(elementindex(zelem_list(idx(i1),i2)))
            dqdwconstraint(icount)=dqdwconstraint(icount)+dqdw(jcount,1)
            icount=icount+1
            jcount=jcount+1
          enddo ! i3
        enddo ! i2, loop over all atoms
!!
!! normalize dqdwconstraint per atom
!! if we use this, we often get oscillations:
!!         dqdwconstraint(:)=dqdwconstraint(:)/dble(num_atoms)
!!
        chargeerror=totalcharge_list(idx(i1))-nnchargesum
!!         write(*,'(a,3f14.6)')'chargeerror ',chargeerror,totalcharge_list(idx(i1)),nnchargesum
        abschargeerror=abs(totalcharge_list(idx(i1))-nnchargesum)
!!
        if(abschargeerror.gt.kalmanthresholdc*rmse_totalcharge)then
!!
!! update the electrostatic weights for total charge constraint 
         if(optmodeq.eq.1)then ! Kalman filter
!!
           icount=1
           do i3=1,num_weights_elec(elementindex(zelem_list(idx(i1),i2)))
             do i4=1,nelem
               weightsconstraint(icount)=weights_elec(i3,i4)
               icount=icount+1
             enddo
           enddo
!!
           if((mpisize.eq.1).or.lompmkl)then
             call updatekalman(kalcdim,corrcdim,&
               constraintdim,kalmanlambdac,kalmannuec,&
               weightsconstraint,dqdwconstraint,&
               corrmatrixc,chargeerror)
           else
             call updatekalman_para(paramode,kalcdim,corrcdim,&
               constraintdim,kalmanlambdac,kalmannuec,&
               weightsconstraint,dqdwconstraint,&
               corrmatrixc,chargeerror)
           endif
!!
           icount=1
           do i3=1,num_weights_elec(elementindex(zelem_list(idx(i1),i2)))
             do i4=1,nelem
               weights_elec(i3,i4)=weightsconstraint(icount)
               icount=icount+1
             enddo
           enddo
!!
         elseif(optmodeq.eq.2)then ! conjugate gradient'
           write(ounit,*)'CG optimization not yet implemented for constraint'
           stop !'
         elseif(optmodeq.eq.3)then ! steepest descent'
!!
           icount=1
           do i3=1,num_weights_elec(elementindex(zelem_list(idx(i1),i2)))
             do i4=1,nelem
               weightsconstraint(icount)=weights_elec(i3,i4)
               icount=icount+1
             enddo
           enddo
!!
           call updatesteepest(&
                constraintdim,weightsconstraint,dqdwconstraint,&
                chargeerror,steepeststepq)
!!
           icount=1
           do i3=1,num_weights_elec(elementindex(zelem_list(idx(i1),i2)))
             do i4=1,nelem
               weights_elec(i3,i4)=weightsconstraint(icount)
               icount=icount+1
             enddo
           enddo
!!
         endif ! optmodeq
        endif ! adaptive filter for charge constraint
!!
        deallocate(dqdwconstraint)
        deallocate(weightsconstraint)
      endif !! lchargeconstraint
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
      enddo ! i1, loop over all points
!!
      return
      end
