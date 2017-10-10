!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine fitting_batch(iseed)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use symfunctions 
      use nnewald
      use nnshort_atomic
      use structures
!!
      implicit none
!!
!! for dimensions
      integer num_weightsshortfree(nelem)                       ! internal
      integer num_weightsshortfixed(nelem)                      ! internal
      integer num_weightsewaldfree(nelem)                       ! internal
      integer num_weightsewaldfixed(nelem)                      ! internal
!! counters
      integer icount                                            ! internal
      integer npoints                                           ! internal
      integer point                                             ! internal
      integer pointe                                            ! internal
      integer ncount                                            ! internal
      integer ndone                                             ! internal
      integer countepoch                                        ! internal
      integer i1,i2                                             ! internal 
      integer optepoch                                          ! internal
      integer optepoche                                         ! internal
      integer belowmaxenergy                                    ! internal
      integer belowfmax(nelem)                                  ! internal
!!
      integer iseed                                             ! in
      integer jseed                                             ! internal
      integer kseed                                             ! internal
      integer lseed                                             ! internal
      integer mseed                                             ! internal
      integer nseed                                             ! internal
      integer ntrain                                            ! internal
      integer ntest                                             ! internal
      integer day                                               ! internal
      integer idx(nblock)                                       ! internal
      integer nstruct                                           ! internal
      integer n_start,n_end                                     ! internal
      integer numbere                                           ! internal
      integer numberf                                           ! internal
      integer numberq                                           ! internal
      integer ntrainatoms                                       ! internal
      integer ntestatoms                                        ! internal
      integer trainelem(nelem)                                  ! internal
      integer testelem(nelem)                                   ! internal
      integer wconstraint(maxnum_weights_short_atomic,nelem)            ! internal
      integer wconstrainte(maxnum_weights_elec,nelem)           ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidx  ! internal 
      integer, dimension(:,:)  , allocatable :: wconstraintidxe ! internal 
      integer, dimension(:)  , allocatable :: fitstat           ! internal 
      integer, dimension(:,:,:)  , allocatable :: fitstatf      ! internal 
      integer, dimension(:,:)  , allocatable :: fitstatq        ! internal 
!! Kalman matrix dimensions:
      integer corrdim(nelem)                                    ! internal
      integer corrfdim(nelem)                                   ! internal
      integer corredim(nelem)                                   ! internal
      integer corrcdim
      integer kaldim(nelem)
      integer kaledim(nelem)
      integer kalcdim
      integer maxkaldim
      integer maxkaledim
      integer maxcorrdim                                        ! internal
      integer maxcorrfdim                                       ! internal
      integer maxcorredim                                       ! internal

      integer nshort(nelem)                                     ! internal
      integer newald                                            ! internal
      real*8 dshortdw(maxnum_weights_short_atomic,nelem)                ! internal
      real*8 dewalddw(maxnum_weights_elec,nelem)                ! internal
      real*8 errorshort(nelem)                                  ! internal
      real*8 errorewald(nelem)                                  ! internal
      real*8 dummy                                              ! internal
      real*8 rdummy(nelem)                                      ! internal

!! symmetry function related arrays
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)
!! weights
      real*8 optweights_short(maxnum_weights_short_atomic,nelem)
      real*8 optweights_ewald(maxnum_weights_elec,nelem)
      real*8 weights_short_old(maxnum_weights_short_atomic,nelem)
      real*8 weights_ewald_old(maxnum_weights_elec,nelem)
      real*8 weights_short_veryold(maxnum_weights_short_atomic,nelem)
      real*8 weights_ewald_veryold(maxnum_weights_elec,nelem)
!! RMSEs
      real*8 rmse_short
      real*8 rmse_charge
      real*8 rmse_totalcharge
      real*8 rmse_etot
      real*8 rmse_ewald
      real*8 rmse_short_test
      real*8 rmse_short_test_old
      real*8 rmse_charge_test
      real*8 rmse_charge_test_old
      real*8 rmse_totalcharge_test
      real*8 rmse_etot_test
      real*8 rmse_ewald_test
      real*8 rmse_force_s
      real*8 rmse_force_s_test
      real*8 rmse_force_t
      real*8 rmse_force_t_test
      real*8 rmse_force_e
      real*8 rmse_force_e_test
      real*8 rmse_short_ref
      real*8 rmse_charge_ref
      real*8 rmse_force_s_ref
      real*8 rmse_totalcharge_ref
!! MADs
      real*8 mad_short
      real*8 mad_short_test
      real*8 mad_ewald
      real*8 mad_ewald_test
      real*8 mad_etot
      real*8 mad_etot_test
      real*8 mad_charge
      real*8 mad_charge_test
      real*8 mad_totalcharge
      real*8 mad_totalcharge_test
      real*8 mad_force_s
      real*8 mad_force_s_test
      real*8 mad_force_t
      real*8 mad_force_t_test
      real*8 mad_force_e
      real*8 mad_force_e_test
!! for final output
      real*8 optrmse_short
      real*8 optrmse_charge
      real*8 optrmse_totalcharge
      real*8 optrmse_etot
      real*8 optrmse_ewald
      real*8 optrmse_short_test
      real*8 optrmse_charge_test
      real*8 optrmse_totalcharge_test
      real*8 optrmse_etot_test
      real*8 optrmse_ewald_test
      real*8 optrmse_force_s
      real*8 optrmse_force_s_test
      real*8 optrmse_force_e
      real*8 optrmse_force_e_test
      real*8 fmin(nelem)
      real*8 fmax(nelem)

      real*8 eshortmin
      real*8 eshortmax
      real*8 eshortav                                         ! internal

!! Kalman filter 
      real*8 kalmanthreshold_temp                             ! internal
      real*8 kalmanthresholdf_temp                            ! internal
      real*8 kalmanthresholde_temp                            ! internal
      real*8, dimension(:,:), allocatable :: corrmatrix_list
      real*8, dimension(:,:), allocatable :: corrmatrixf_list
      real*8, dimension(:,:), allocatable :: corrmatrixe_list
      real*8, dimension(:)  , allocatable :: corrmatrixc 
!! arrays for convergence vector
      real*8 wshift(nelem)                                    ! internal
      real*8 wshift2(nelem)                                   ! internal
      real*8 wshifte(nelem)                                   ! internal
      real*8 wshifte2(nelem)                                  ! internal

!! miscellaneous
      real*8 epochtime                                        ! internal
      real*8 timestart                                        ! internal
      real*8 timeend                                          ! internal
      real*8 tounit                                           ! internal
      real*8 ztemp                                            ! internal
      real*8 avcharge(nelem)                                  ! internal
      real*8 stddevcharge(nelem)                              ! internal
      real*8 eshortstddev                                     ! internal
      real*8 chargemin(nelem)                                        ! internal
      real*8 chargemax(nelem)                                        ! internal
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! preparations and initializations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! do not fit charges if they are fixed anyway
      if(nn_type_elec.eq.3)then
        lelec=.false.
      endif
!!
!!    allocate arrays of structures module
      call allocatestructures()
!!
      chargemin(:)           = 100000.d0
      chargemax(:)           = -100000.d0
!!
!! initialization of local Kalman filter parameters
      kalmanthreshold_temp   = kalmanthreshold
      kalmanthresholdf_temp  = kalmanthresholdf
      kalmanthresholde_temp  = kalmanthresholde
!! initialization of scaling parameters
      minvalue_short_atomic(:,:)          = 0.0d0
      maxvalue_short_atomic(:,:)          = 0.0d0
      avvalue_short_atomic(:,:)           = 0.0d0
      minvalue_elec(:,:)         = 0.0d0
      maxvalue_elec(:,:)         = 0.0d0
      avvalue_elec(:,:)          = 0.0d0
!! set filenames for weights files
      allocate(filenamewe(nelem))
      allocate(filenamews(nelem))
      allocate(filenamewp(npairs))
      call getweightfilenames()
!!
!! initialization of reference data arrays
      totalcharge_list(:)    = 0.0d0
      totalenergy_list(:)    = 0.0d0
      shortenergy_list(:)    = 0.0d0
      elecenergy_list(:)     = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
      shortforce_list(:,:,:) = 0.0d0
      elecforce_list(:,:,:)  = 0.0d0
      lattice_list(:,:,:)    = 0.0d0
!! setting the energy and force unit converter
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV 
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0   !  stay with Ha 
      endif
!! initialization of reference RMSEs
      rmse_short_ref         =0.0d0
      rmse_charge_ref        =0.0d0
      rmse_totalcharge_ref   =0.0d0
      rmse_force_s_ref       =0.0d0
!!
!! generating local copies of seed for various purposes
!! iseed is just used for random order training below
      jseed                  = iseed ! for initial short range weights 
      kseed                  = iseed ! for forcernd 
      lseed                  = iseed ! for energyrnd
      mseed                  = iseed ! for chargernd 
      nseed                  = iseed ! for initial electrostatic weights
!! initializations of counters
      countepoch             = 0
      point                  = 0
      pointe                 = 0
      npoints                = 0
      numbere                = 0
      numberf                = 0
      numberq                = 0
      optepoch               = 0
      optepoche              = 0
!! initialization of timing variables
      timeend                = 0.0d0
      timestart              = 0.0d0
      day                    = 0
!!
      write(ounit,*)'optmode ',optmodee,optmodef,optmodeq
!!
!! determine starting time
      call abstime(timestart,day)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end of preparations and initializations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! check if we have the same number of structures in function.data and functione.data
!! This might not be the case if we merge separate fits for short range and charge NNs
!! In this case the results would be nonsense, because both fits would use the same trainstruct.data file
      if(mpirank.eq.0)then
        if(lshort.and.lelec.and.(nn_type_elec.eq.1))then
          call checkfunction()
        endif   
      endif
!!
!! get information on fixed weights from input.nn file
      if(mpirank.eq.0)then
        wconstraint(:,:)=0
        wconstrainte(:,:)=0
        num_weightsshortfree(:)=num_weights_short_atomic(:)
        num_weightsshortfixed(:)=0
        num_weightsewaldfree(:)=num_weights_elec(:)
        num_weightsewaldfixed(:)=0
        if(lshort.and.lfixweights)then
          call getfixedweights(0,nelem,&
            maxnum_layers_short_atomic,num_layers_short_atomic,&
            nodes_short_atomic,maxnum_weights_short_atomic,num_weights_short_atomic,windex_short_atomic,&
            num_weightsshortfree,num_weightsshortfixed,wconstraint)
        endif
        if(lelec.and.(nn_type_elec.eq.1).and.lfixweights)then
          call getfixedweights(1,nelem,&
            maxnum_layers_elec,num_layers_elec,&
            nodes_elec,maxnum_weights_elec,num_weights_elec,windex_elec,&
            num_weightsewaldfree,num_weightsewaldfixed,wconstrainte)
        endif
!! write summary
        do i1=1,nelem
          write(ounit,*)' Element ',element(i1)
          write(ounit,*)'  Total number of short range weights       ',num_weights_short_atomic(i1)
          write(ounit,*)'  Number of optimized short range weights   ',num_weightsshortfree(i1)
          write(ounit,*)'  Number of constrained short range weights ',num_weightsshortfixed(i1)
          write(ounit,*)'  Total number of charge weights            ',num_weights_elec(i1)
          write(ounit,*)'  Number of optimized charge weights        ',num_weightsewaldfree(i1)
          write(ounit,*)'  Number of constrained charge weights      ',num_weightsewaldfixed(i1)
        enddo ! i1
      endif ! mpirank.eq.0
!! distribute weight constraints to all nodes here
      call mpi_bcast(wconstraint,maxnum_weights_short_atomic*nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(wconstrainte,maxnum_weights_elec*nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsshortfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsshortfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsewaldfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsewaldfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!! prepare constraint arrays
      allocate(wconstraintidx(maxnum_weights_short_atomic,nelem))  
      allocate(wconstraintidxe(maxnum_weights_elec,nelem))  
!! get wconstraintidx
      do i1=1,nelem
        icount=0
        do i2=1,num_weights_short_atomic(i1)
          if(wconstraint(i2,i1).eq.0)then
            icount=icount+1
            wconstraintidx(icount,i1)=i2
          endif
        enddo
      enddo
!! get wconstraintidxe
      do i1=1,nelem
        icount=0
        do i2=1,num_weights_elec(i1)
          if(wconstrainte(i2,i1).eq.0) then
            icount=icount+1
            wconstraintidxe(icount,i1)=i2
          endif
        enddo
      enddo
!!
!! get dimensions for the Kalman filter
      if((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))then
        call getkaldims(nelem,&
          kaldim,kaledim,kalcdim,corrdim,corrfdim,corredim,corrcdim,&
          maxkaldim,maxkaledim,maxcorrdim,maxcorrfdim,maxcorredim,&
          num_weightsshortfree,num_weightsewaldfree)
      endif
!!
!! allocate Kalman matrices
!! we have to allocate all matrices all the time because of subroutine getkalmanmatrices
      allocate(corrmatrix_list(maxcorrdim,nelem))
      allocate(corrmatrixf_list(maxcorrfdim,nelem))
      allocate(corrmatrixe_list(maxcorredim,nelem))
      allocate(corrmatrixc(corrcdim)) 
!!
!! get the values for the Kalman matrices and parameters (new or from file)
      call getkalmanmatrices(nelem,&
        iseed,kseed,lseed,mseed,&
        kaldim,kaledim,corrdim,corrfdim,corredim,corrcdim,&
        maxcorrdim,maxcorrfdim,maxcorredim,&
        num_weightsshortfree,num_weightsewaldfree,&
        corrmatrix_list,corrmatrixf_list,corrmatrixe_list,corrmatrixc,&
        kalmanlambda)
      call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(kseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mseed,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
!! count the number of training and test points
      if(mpirank.eq.0)then
        call countpoints(1,nelem,&
          ntrain,ntest,ntrainatoms,ntestatoms,&
          trainelem,testelem)
      endif ! mpirank.eq.0
!!
      call mpi_bcast(ntrain,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntest,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntrainatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntestatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(trainelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(testelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
!! allocate arrays for fitting statistics output
      allocate(fitstat(ntrain))
      fitstat(:)=0
      allocate(fitstatf(3,max_num_atoms,ntrain))
      fitstatf(:,:,:)=0
      allocate(fitstatq(max_num_atoms,ntrain))
      fitstatq(:,:)=0
!!
!! if growth mode is not used make sure the full training set is used
      if(.not.lgrowth)then
        ngrowth=ntrain
      endif
!!
!! get maxcutoff of short range symmetry functions
      maxcutoff_short_atomic = 0.0d0
      do i2=1,nelem
        do i1=1,num_funcvalues_short_atomic(i2)
          maxcutoff_short_atomic=max(maxcutoff_short_atomic,funccutoff_short_atomic(i1,i2))
        enddo ! i1
      enddo
!! get maxcutoffe of electrostatic symmetry functions
      maxcutoff_elec = 0.0d0
      do i2=1,nelem
        do i1=1,num_funcvalues_elec(i2)
          maxcutoff_elec=max(maxcutoff_elec,funccutoff_elec(i1,i2))
        enddo ! i1
      enddo ! i2
!!
      if(mpirank.eq.0)then
!! get scaling data for the short range symmetry functions
        if(lshort)then
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit) !'
          call getscale(nelem,max_num_atoms,0,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            ntrain,symfunction_short_atomic_list,minvalue_short_atomic,&
            maxvalue_short_atomic,avvalue_short_atomic)
          close(symunit)
        endif ! lshort
!!
!! get scaling data for the charge symmetry functions 
        if(lelec.and.(nn_type_elec.eq.1))then
          open(symunit,file='functione.data',form='formatted',status='old')
          rewind(symunit) !'
          call getscale(nelem,max_num_atoms,1,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            ntrain,symfunction_elec_list,minvalue_elec,maxvalue_elec,avvalue_elec)
          close(symunit)
        endif ! lelec
!!
!! analyze structures in trainstruct.data and write chemical composition to runner.out
        if(lshort.and.(nelem.le.4))then
          call analyzeinput(0) 
          call analyzeinput(1) 
        elseif(lshort.and.(nelem.gt.3))then
          write(ounit,*)'WARNING: detailed analysis of structures not possible for more than 3 elements'
        endif !'
!!
!! get statistics of energies in training set (just for process 0)
        call getenergystatistics(1,belowmaxenergy,&
          eshortav,eshortstddev,eshortmin,eshortmax)
!!
!! get statistics of forces in training set (just for process 0)
        if((luseforces.and.lshort).or.lfinalforce)then
          call getforcestatistics(belowfmax,fmin,fmax)
        endif
!!
!! get the average charges
        if(lelec)then
          avcharge(:)=0.0d0
          if(mpirank.eq.0)then
            call getavcharge(ntrain,&
              avcharge,stddevcharge,&
              chargemin,chargemax)
          endif ! mpirank.eq.0
          call mpi_bcast(avcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_bcast(stddevcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        endif ! lelec
!!
!! write short range scaling factors to file scaling.data
        if(lshort)then
          rdummy(:)=0.0d0
          call writescale(nelem,1,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
            eshortmin,eshortmax,rdummy,rdummy)
        endif ! lshort
!!
!! write electrostatic scaling factors to file scalinge.data
        if(lelec.and.(nn_type_elec.eq.1))then
          dummy=0.0d0
          call writescale(nelem,3,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            dummy,dummy,chargemin,chargemax)
        endif ! lelec
      endif ! mpirank.eq.0
!!
!! distribute scaling.data to all processes
      if(lshort)then
        call mpi_bcast(minvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortstddev,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(minvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'number of training points ',ntrain !'
        write(ounit,*)'number of training atoms  ',ntrainatoms
        if(luseforces.or.lfinalforce)then
          write(ounit,*)'number of training forces ',3*ntrainatoms
        endif
        write(ounit,*)'number of testing points  ',ntest
        write(ounit,*)'number of testing atoms   ',ntestatoms
        if(luseforces.or.lfinalforce)then
          write(ounit,*)'number of testing forces  ',3*ntestatoms
        endif
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Number of atoms for each element:   '
        write(ounit,*)'           training:    testing:   '
        do i1=1,nelem
          write(ounit,'(i3,x,a3,x,2i12)')i1,element(i1),trainelem(i1),testelem(i1)
        enddo
      endif ! mpirank.eq.0 '
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! preparation of initial weight parameters'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! first initialize all weights randomly no matter if we restart a fit or not
      if(mpirank.eq.0)then
        call initialweights(nelem,&
             maxnum_weights_short_atomic,num_weights_short_atomic,num_weights_elec,&
             maxnum_layers_short_atomic,num_layers_short_atomic,windex_short_atomic,nodes_short_atomic,&
             jseed,nseed,weights_short_atomic,weights_elec)
!!
!! If systematic weights are requested, overwrite random weights for short range NN
        if(lsysweights)then
          call systematicweights(nelem,maxnum_layers_short_atomic,&
            num_layers_short_atomic,nodes_short_atomic,&
            maxnum_weights_short_atomic,&
            weights_short_atomic,weights_min,weights_max)
        endif
!! If systematic weights are requested, overwrite random weights for electrostatic NN
        if(lsysweightse)then
          call systematicweights(nelem,maxnum_layers_elec,&
            num_layers_elec,nodes_elec,&
            maxnum_weights_elec,&
            weights_elec,weightse_min,weightse_max)
        endif
!!
!! if requested overwrite short range weights with weights from weights.XXX.data
        if(luseoldweightsshort.and.lshort)then
          call readweights(0,nelem,&
            maxnum_weights_short_atomic,num_weights_short_atomic,&
            weights_short_atomic)
        endif
!!
!! if requested overwrite charge weights with weights from weightse.XXX.data
        if(luseoldweightscharge.and.lelec.and.(nn_type_elec.eq.1))then
          call readweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            weights_elec)
        endif
!!
!! write initial weights to file
        if(lshort)then
          call writeweights(0,nelem,maxnum_weights_short_atomic,&
            maxnum_layers_short_atomic,num_layers_short_atomic,&
            nodes_short_atomic,weights_short_atomic)
        endif ! lshort
        if(lelec.and.(nn_type_elec.eq.1))then
          call writeweights(1,nelem,maxnum_weights_elec,&
            maxnum_layers_elec,num_layers_elec,&
            nodes_elec,weights_elec)
        endif ! lelec
      endif ! mpirank.eq.0
!!
!! keep copy of weights for the calculation of the weight change wshift/wshifte
      weights_short_old(:,:)=weights_short_atomic(:,:)
      weights_ewald_old(:,:)=weights_elec(:,:)
!!
!! distribute weights_short and weights_ewald to all processes
      if(lshort)then
        call mpi_bcast(weights_short_atomic,nelem*maxnum_weights_short_atomic,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(weights_elec,nelem*maxnum_weights_elec,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end of preparation of initial weight parameters'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!! determine and print initialization time
      call abstime(timeend,day)
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.2)')&
          ' initialization time (min):',(timeend-timestart)/60.d0
      endif ! mpirank.eq.0
!!'
!! set time to zero for epoch 0 time measurement
      timestart              =0.0d0
      day                    =0
      call abstime(timestart,day)
!!
      rmse_short             = 0.0d0
      rmse_charge            = 0.0d0
      rmse_totalcharge       = 0.0d0
      rmse_etot              = 0.0d0
      rmse_ewald             = 0.0d0
      rmse_short_test        = 0.0d0
      rmse_charge_test       = 0.0d0
      rmse_totalcharge_test  = 0.0d0
      rmse_etot_test         = 0.0d0
      rmse_ewald_test        = 0.0d0
      rmse_force_s           = 0.0d0
      rmse_force_s_test      = 0.0d0
      rmse_force_t           = 0.0d0
      rmse_force_t_test      = 0.0d0
      rmse_force_e           = 0.0d0
      rmse_force_e_test      = 0.0d0
      mad_short              = 0.0d0
      mad_short_test         = 0.0d0
      mad_charge             = 0.0d0
      mad_charge_test        = 0.0d0
      mad_etot               = 0.0d0
      mad_etot_test          = 0.0d0
      mad_ewald              = 0.0d0
      mad_ewald_test         = 0.0d0
      mad_totalcharge        = 0.0d0
      mad_totalcharge_test   = 0.0d0
      mad_force_s            = 0.0d0
      mad_force_s_test       = 0.0d0
      mad_force_t            = 0.0d0
      mad_force_t_test       = 0.0d0
      mad_force_e            = 0.0d0
      mad_force_e_test       = 0.0d0
      maxerror_eshort_train    = 0.0d0
      maxerror_eshort_test     = 0.0d0
      imaxerror_eshort_train   = 0
      imaxerror_eshort_test    = 0
      maxerror_elec_train     = 0.0d0
      maxerror_elec_test      = 0.0d0
      imaxerror_elec_train    = 0
      imaxerror_elec_test     = 0
      maxerror_etot_train      = 0.0d0
      maxerror_etot_test       = 0.0d0
      imaxerror_etot_train     = 0
      imaxerror_etot_test      = 0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!! (rmse_short,rmse_charge,rmse_totalcharge,rmse_ewald,rmse_etot,rmse_force_s,rmse_force_e)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror(0,countepoch,ntrain,&
        imaxerror_eshort_train,imaxerror_elec_train,imaxerror_etot_train,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_ewald,rmse_etot,rmse_force_s,rmse_force_t,&
        rmse_force_e,&
        mad_short,mad_charge,mad_totalcharge,&
        mad_ewald,mad_etot,&
        mad_force_s,mad_force_t,mad_force_e,&
        maxerror_eshort_train,maxerror_elec_train,maxerror_etot_train)
!
! CMH
      write(*,*) 'CMH TEST FORCES TEST'
! CMH
!

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of the first training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! get new references for adaptive Kalman filter
      if(lshort) rmse_short_ref      =rmse_short
      if(lshort) rmse_force_s_ref    =rmse_force_s
      if(lelec)  rmse_charge_ref     =rmse_charge
      if(lelec)  rmse_totalcharge_ref=rmse_totalcharge
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the first testing error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror(1,countepoch,ntest,&
        imaxerror_eshort_test,imaxerror_elec_test,imaxerror_etot_test,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_short_test,rmse_charge_test,rmse_totalcharge_test,&
        rmse_ewald_test,rmse_etot_test,rmse_force_s_test,rmse_force_t_test,&
        rmse_force_e_test,&
        mad_short_test,mad_charge_test,mad_totalcharge_test,&
        mad_ewald_test,mad_etot_test,&
        mad_force_s_test,mad_force_t_test,mad_force_e_test,&
        maxerror_eshort_test,maxerror_elec_test,maxerror_etot_test)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of first test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! set first references for determination of optweights files
      rmse_short_test_old  = rmse_short_test
      rmse_charge_test_old = rmse_charge_test
!!
      call abstime(timeend,day)
      epochtime=timeend-timestart
!!
!! write RMSE headers
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------'
        if(lshort)then
          write(ounit,'(a,f14.3,a)')' ### WARNING ### just short range energies below ',&
            maxenergy,' Ha/atom are used for fitting and Eshort RMSE!'
          write(ounit,'(a,f14.3,a)')' => Fitted energy range has width of ',&
            maxenergy-eshortmin,' Ha/atom = ',(maxenergy-eshortmin)*27.211d0,' eV/atom'
          write(ounit,'(a,i10)')' => Number of short range training energies below max_energy: ',belowmaxenergy 
        endif
        if((luseforces.and.lshort).or.lfinalforce)then
          write(ounit,'(a,f14.3,a)')' ### WARNING ### just force components below     ',&
            maxforce,' Ha/Bohr are used for fitting and Fshort RMSE!'
          do i1=1,nelem
            write(ounit,'(x,a2,x,a,f14.3,a,f14.3,a)')element(i1),&
              ' => Fitted force range has width of ',maxforce-fmin(i1),' Ha/Bohr = ',(maxforce-fmin(i1))*27.211d0,' eV/Bohr'
            write(ounit,'(x,a2,x,a,i10)')element(i1),' => Number of short range training forces below max_force: ',belowfmax(i1) 
          enddo
        endif
        write(ounit,*)'-------------------------------------------------------------------------------'
        if(lelec.and.lupdatebyelement)then
          write(ounit,*)'### WARNING ### lupdatebyelement.eq.T => atom charge error refers only to this element'
          write(ounit,*)'                                      => electrostatic energies are not printed'
          write(ounit,*)'                                      => total charges are not printed'
          write(ounit,*)'                                      => printed total energies are meaningless'
          write(ounit,*)'-------------------------------------------------------------------------------'
          rmse_totalcharge      =0.0d0
          rmse_totalcharge_test =0.0d0
          rmse_totalcharge_ref  =0.0d0
          rmse_ewald            =0.0d0
          rmse_ewald_test       =0.0d0
          rmse_force_e          =0.0d0
          rmse_force_e_test     =0.0d0
          rmse_force_t          =0.0d0
          rmse_force_t_test     =0.0d0
          mad_totalcharge       =0.0d0
          mad_totalcharge_test  =0.0d0
          mad_ewald             =0.0d0
          mad_ewald_test        =0.0d0
          mad_force_e           =0.0d0
          mad_force_e_test      =0.0d0
          mad_force_t           =0.0d0
          mad_force_t_test      =0.0d0
        endif
        if(lshort.and.lupdatebyelement.and.luseforces)then
          write(ounit,*)'### WARNING ### lupdatebyelement.eq.T => atom forces error refers only to this element'
          write(ounit,*)'                                      => short range energies are half meaningless'
          write(ounit,*)'                                      => printed total energies are meaningless'
          write(ounit,*)'-------------------------------------------------------------------------------'
        endif
        if(fitting_unit.eq.1)then
          write(ounit,*)'RMSEs (energies: eV/atom, charges: e, forces: eV/Bohr):'
        elseif(fitting_unit.eq.2)then
          write(ounit,*)'RMSEs (energies: Ha/atom, charges: e, forces: Ha/Bohr):'
        else
          write(ounit,*)'Error: unknown energy unit in fitting.f90'
          stop
        endif
        write(ounit,'(6a)')'                      --- E_short: --- ',&
                         '     --- el. energies: ---',&
                         '    --- total energies: ---',&
                         '    --- atom charges: --- ',&
                         '   --- total charges: ---',&
                         ' - time -'
        write(ounit,'(6a)')'                          /atom      ',&
                         '             /atom        ',&
                         '             /atom         ',&
                         '             /e            ',&
                         '             /e           ',&
                         '   min'
        write(ounit,'(5a)')'       epoch         train         test',&
                    '        train         test',&
                    '        train         test',&
                    '        train         test',&
                    '        train         test'
!! CAUTION: If the training data is constructed with short and ewald and
!! then lshort is set false to fit the charges, then rmse_etot and rmse_ewald
!! cannot be the same because they have different DFT reference energies
!! (for etot it is the total DFT energy and for the electrostatic energy
!! it is only the DFT ewald energy, from the NN we in both cases have just the
!! ewald energy)
!!
!! write RMSEs'
        write(ounit,'(a8,i5,x,10f13.6,f8.2)') ' ENERGY ',countepoch,&
            rmse_short*tounit,rmse_short_test*tounit,&
            rmse_ewald*tounit,rmse_ewald_test*tounit,&
            rmse_etot*tounit,rmse_etot_test*tounit,&
            rmse_charge,rmse_charge_test,&
            rmse_totalcharge,rmse_totalcharge_test,&
            epochtime/60.d0
        if(luseforces)then
          write(ounit,'(a8,i5,x,6f13.6)')' FORCES ',countepoch,&
            rmse_force_s*tounit,rmse_force_s_test*tounit,&
            rmse_force_e*tounit,rmse_force_e_test*tounit,& 
            rmse_force_t*tounit,rmse_force_t_test*tounit 
        endif
        if(lprintmad)then
          write(ounit,'(a9,i5,2f14.6,3x,2f20.10)')' MADE   ',countepoch,&
            mad_short*tounit,mad_short_test*tounit,&
            mad_ewald*tounit,mad_ewald_test*tounit,&
            mad_etot*tounit,mad_etot_test*tounit,&
            mad_charge,mad_charge_test,&
            mad_totalcharge,mad_totalcharge_test
          if(luseforces)then
            write(ounit,'(a8,i5,x,6f13.6)')' MADF   ',countepoch,&
              mad_force_s*tounit,mad_force_s_test*tounit,&
              mad_force_e*tounit,mad_force_e_test*tounit,&
              mad_force_t*tounit,mad_force_t_test*tounit
          endif
        endif
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all training epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do countepoch=1,nepochs
        if(mpirank.eq.0) write(debugunit,'(a6,i6)')'epoch ',countepoch
        call abstime(timestart,day)
!!
!! epoch-specific counter initializations
!!
!!        ncount                = ntrain
        ndone                 = 0
        ncount = (int((countepoch-1)/growthstep)+1)*ngrowth    
        ncount = min(ncount,ntrain)
        point                 = 0
        pointe                = 0
        numbere               = 0
        numberf               = 0
        numberq               = 0
!!
!! epoch-specific averaging initializations
        nshort(:)     = 0
        newald        = 0
        dshortdw(:,:) = 0.0d0
        dewalddw(:,:) = 0.0d0
        errorshort(:) = 0.0d0
        errorewald(:) = 0.0d0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mpirank.eq.0)then
        open(symunit,file='function.data',form='formatted',status='old')
        rewind(symunit)
        open(trainstructunit,file='trainstruct.data',form='formatted',status='old')
        rewind(trainstructunit)
        if(luseforces.or.lfinalforce)then
          open(trainfunit,file='trainforces.data',form='formatted',status='old')
          rewind(trainfunit) 
          if(lelec)then
            open(trainfeunit,file='trainforcese.data',form='formatted',status='old')  
            rewind(trainfeunit) !' 
          endif
        endif ! luseforces
        if(lelec.and.(nn_type_elec.eq.1))then
          open(symeunit,file='functione.data',form='formatted',status='old')
          rewind(symeunit) !'
        endif ! lelec
      endif ! mpirank.eq.0
!!
!!
!! loop block-wise over training points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 11   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!! get index for random order of training points => idx
      call getidx(npoints,iseed,idx)
!!
!! read npoint data sets      
!! do all file reading for training data here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      if(mpirank.eq.0)then
        if(lshort)then
          call readfunctions(1,symunit,npoints,nelem,&
            max_num_atoms,maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            symfunction_short_atomic_list)
        endif ! lshort
!!
        if(lelec.and.(nn_type_elec.eq.1))then
          call readfunctions(1,symeunit,npoints,nelem,&
            max_num_atoms,maxnum_funcvalues_elec,num_funcvalues_elec,&
            symfunction_elec_list)
        endif ! lelec
!!
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! is needed for lshort (force fitting) and lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
        call getstructures(trainstructunit,npoints)
!!
!! read short range forces from trainforces.data
        if((luseforces.and.lshort).or.(lfinalforce.and.(nepochs.eq.countepoch)))then
          call readforces(trainfunit,npoints,shortforce_list)
        endif
!!
!! read electrostatic forces from trainforcese.data
        if((luseforces.or.(lfinalforce.and.(nepochs.eq.countepoch))).and.lelec)then
          call readforces(trainfeunit,npoints,elecforce_list)
        endif
      endif ! mpirank.eq.0
!!
!! distribute the data to all processes
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalcharge_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elecenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_short_atomic_list,&
        nblock*max_num_atoms*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_elec_list,nblock*max_num_atoms*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortforce_list,nblock*max_num_atoms*3,&
        mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elecforce_list,nblock*max_num_atoms*3,&
        mpi_real8,0,mpi_comm_world,mpierror)
!!
!! end of file reading for training data 
!!--------------------------------------------------------------------
!!
!! determine which structures of this block should be scaled by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! short range part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lshort)then
!! scale symmetry functions for the short-range interaction
        call scalesymfit_para(1,nelem,&
          nstruct,n_start,n_end,&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,symfunction_short_atomic_list,&
          scmin_short_atomic,scmax_short_atomic)
!!
!! get part of dshortdw for this block of points
        call getdshortdw(npoints,point,&
          ndone,num_weightsshortfree,&
          wconstraintidx,&
          numbere,numberf,nshort,&
          lseed,ntrain,fitstat,fitstatf,&
          kseed,errorshort,&
          rmse_short,rmse_force_s_ref,&
          minvalue_short_atomic,maxvalue_short_atomic,dshortdw)
!!
      endif ! lshort
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end short range part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! electrostatic part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(lelec.and.(nn_type_elec.eq.1))then
!!
!! scale symmetry functions for the charge prediction 
        call scalesymfit_para(3,nelem,&
          nstruct,n_start,n_end,&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          minvalue_elec,maxvalue_elec,avvalue_elec,symfunction_elec_list,&
          scmin_elec,scmax_elec)
!!
!! get part of dewalddw for this block of points
!!        call getdewalddw()
          write(ounit,*)'Error: batch learning not implemented for ewald'
          stop !'
!!
      endif ! lelec
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end electrostatic part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! count finished structures
      ndone=ndone+npoints
!!
      if(ncount.gt.0) goto 11 ! do next block of points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end block-wise loop over training points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! Do weight update here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!      write(ounit,*)'dshortdw'
!!      do i1=1,maxnum_weightsshort
!!        write(ounit,*)dshortdw(i1,1)
!!      enddo

!!
      if(lshort)then
         call updateshort(nshort,&
         maxkaldim,kaldim,corrdim,maxcorrdim,&
         num_weightsshortfree,&
         wconstraintidx,&
         dshortdw,errorshort,&
         corrmatrix_list)
      endif ! lshort
!!
      if(lelec)then
!!        call updateewald()
      endif
!!
!!
!! write Kalman data for restart if requested
      if(lsavekalman)then
        if((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))then
          call writekalman(nelem,&
            iseed,kseed,lseed,mseed,&
            maxcorrdim,maxcorredim,&
            kaldim,kaledim,corrcdim,&
            num_weightsshortfree,num_weightsewaldfree,&
            kalmanlambda,&
            corrmatrix_list,corrmatrixe_list,corrmatrixc)
        endif !lsavekalman
      endif ! optmode
!!
      if(mpirank.eq.0)then
        close(symunit)
        close(trainstructunit)
        if(luseforces.or.lfinalforce)close(trainfunit)
        if(lelec.and.(nn_type_elec.eq.1)) close(symeunit)
        if(lelec.and.(luseforces.or.lfinalforce)) close(trainfeunit)
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! write final weights of this epoch to files
      if(mpirank.eq.0)then
        call getfilenames(countepoch)
        if(lshort)then
          call writeweights(0,nelem,&
            maxnum_weights_short_atomic,&
            maxnum_layers_short_atomic,num_layers_short_atomic,&
            nodes_short_atomic,weights_short_atomic)
        endif
        if(lelec.and.(nn_type_elec.eq.1))then
          call writeweights(1,nelem,&
            maxnum_weights_elec,&
            maxnum_layers_elec,num_layers_elec,&
            nodes_elec,weights_elec)
        endif
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! initialize rmse's for next error calculation 
!!
      rmse_short             = 0.0d0
      rmse_charge            = 0.0d0
      rmse_totalcharge       = 0.0d0
      rmse_etot              = 0.0d0
      rmse_ewald             = 0.0d0
      rmse_short_test        = 0.0d0
      rmse_charge_test       = 0.0d0
      rmse_totalcharge_test  = 0.0d0
      rmse_etot_test         = 0.0d0
      rmse_ewald_test        = 0.0d0
      rmse_force_s           = 0.0d0
      rmse_force_s_test      = 0.0d0
      rmse_force_t           = 0.0d0
      rmse_force_t_test      = 0.0d0
      rmse_force_e           = 0.0d0
      rmse_force_e_test      = 0.0d0
      mad_short              = 0.0d0
      mad_short_test         = 0.0d0
      mad_ewald              = 0.0d0
      mad_ewald_test         = 0.0d0
      mad_etot               = 0.0d0
      mad_etot_test          = 0.0d0
      mad_charge             = 0.0d0
      mad_charge_test        = 0.0d0
      mad_totalcharge        = 0.0d0
      mad_totalcharge_test   = 0.0d0
      mad_force_s            = 0.0d0
      mad_force_s_test       = 0.0d0
      mad_force_t            = 0.0d0
      mad_force_t_test       = 0.0d0
      mad_force_e            = 0.0d0
      mad_force_e_test       = 0.0d0
      maxerror_eshort_train  = 0.0d0
      maxerror_eshort_test   = 0.0d0
      imaxerror_eshort_train = 0
      imaxerror_eshort_test  = 0
      maxerror_elec_train    = 0.0d0
      maxerror_elec_test     = 0.0d0
      imaxerror_elec_train   = 0
      imaxerror_elec_test    = 0
      maxerror_etot_train    = 0.0d0
      maxerror_etot_test     = 0.0d0
      imaxerror_etot_train   = 0
      imaxerror_etot_test    = 0
!!
!! If we are in the final epoch calculate the force error if requested
      if(lfinalforce.and.(nepochs.eq.countepoch))then
        luseforces=.true.
        forcernd=0.0d0
      endif
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate this epoch's training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
      call geterror(0,countepoch,ntrain,&
        imaxerror_eshort_train,imaxerror_elec_train,imaxerror_etot_train,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_ewald,rmse_etot,rmse_force_s,rmse_force_t,&
        rmse_force_e,&
        mad_short,mad_charge,mad_totalcharge,&
        mad_ewald,mad_etot,&
        mad_force_s,mad_force_t,mad_force_e,&
        maxerror_eshort_train,maxerror_elec_train,maxerror_etot_train)
!!
!! get new references for adaptive Kalman filter (this has nothing to do with optweights files)
      if(lshort) rmse_short_ref      =rmse_short
      if(lshort) rmse_force_s_ref    =rmse_force_s
      if(lelec)  rmse_charge_ref     =rmse_charge
      if(lelec)  rmse_totalcharge_ref=rmse_totalcharge
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of the epoch's training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the epoch's testing error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror(1,countepoch,ntest,&
        imaxerror_eshort_test,imaxerror_elec_test,imaxerror_etot_test,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_short_test,rmse_charge_test,rmse_totalcharge_test,&
        rmse_ewald_test,rmse_etot_test,rmse_force_s_test,rmse_force_t_test,&
        rmse_force_e_test,&
        mad_short_test,mad_charge_test,mad_totalcharge_test,&
        mad_ewald_test,mad_etot_test,&
        mad_force_s_test,mad_force_t_test,mad_force_e_test,&
        maxerror_eshort_test,maxerror_elec_test,maxerror_etot_test)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of epoch's test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call abstime(timeend,day)
      epochtime=timeend-timestart
!!
!! calculate change of the short range weights for INFO output
      if(lshort)then
        call getwshift(nelem,maxnum_weights_short_atomic,num_weights_short_atomic,&
          weights_short_atomic,weights_short_old,weights_short_veryold,&
          wshift,wshift2)
      endif ! lshort
!!
!! calculate change of the charge weights for INFO output
      if(lelec.and.(nn_type_elec.eq.1))then
        call getwshift(nelem,maxnum_weights_elec,num_weights_elec,&
          weights_elec,weights_ewald_old,weights_ewald_veryold,&
          wshifte,wshifte2)
      endif
!!
!! keep weights
      if(lshort) weights_short_old(:,:)=weights_short_atomic(:,:)
!! keep weights
      if(lelec.and.(nn_type_elec.eq.1)) weights_ewald_old(:,:)=weights_elec(:,:)
!!
!! write RMSEs
!! If update by element is used: Set some meaningless RMSEs to zero for output
      if(lelec.and.lupdatebyelement)then
        rmse_totalcharge      =0.0d0
        rmse_totalcharge_test =0.0d0
        rmse_totalcharge_ref  =0.0d0
        rmse_ewald            =0.0d0
        rmse_ewald_test       =0.0d0
        rmse_force_e          =0.0d0
        rmse_force_e_test     =0.0d0
        rmse_force_t          =0.0d0
        rmse_force_t_test     =0.0d0
        mad_totalcharge       =0.0d0
        mad_totalcharge_test  =0.0d0
        mad_ewald             =0.0d0
        mad_ewald_test        =0.0d0
        mad_force_e           =0.0d0
        mad_force_e_test      =0.0d0
        mad_force_t           =0.0d0
        mad_force_t_test      =0.0d0
      endif
      if(mpirank.eq.0)then
        write(ounit,'(a8,i5,x,10f13.6,f8.2)')' ENERGY ', countepoch,&
          rmse_short*tounit,rmse_short_test*tounit,&
          rmse_ewald*tounit,rmse_ewald_test*tounit,&
          rmse_etot*tounit,rmse_etot_test*tounit,&
          rmse_charge,rmse_charge_test,&
          rmse_totalcharge,rmse_totalcharge_test,&
          epochtime/60.d0
        if(luseforces)then
          write(ounit,'(a8,i5,x,6f13.6)')' FORCES ',countepoch,&
          rmse_force_s*tounit,rmse_force_s_test*tounit,&
          rmse_force_e*tounit,rmse_force_e_test*tounit,& 
          rmse_force_t*tounit,rmse_force_t_test*tounit 
        endif
        if(lprintmad)then
          write(ounit,'(a9,i5,2f14.6,3x,2f20.10)')' MADE   ',countepoch,&
            mad_short*tounit,mad_short_test*tounit,&
            mad_ewald*tounit,mad_ewald_test*tounit,&
            mad_etot*tounit,mad_etot_test*tounit,&
            mad_charge,mad_charge_test,&
            mad_totalcharge,mad_totalcharge_test
          if(luseforces)then
            write(ounit,'(a8,i5,x,6f13.6)')' MADF   ',countepoch,&
              mad_force_s*tounit,mad_force_s_test*tounit,&
              mad_force_e*tounit,mad_force_e_test*tounit,&
              mad_force_t*tounit,mad_force_t_test*tounit
          endif
        endif
        write(ounit,'(a,i5,3i10,2f20.15)')' INFORMATION USED FOR UPDATE (E,F,Q) ',&
          countepoch,numbere,numberf,numberq !'
      endif ! mpirank.eq.0
!!
!! check if a non-zero number of structures has been used for the energy update
      if((numbere.eq.0).and.(lshort).and.(kalmanthreshold.gt.0.0d0).and.(energyrnd.gt.0.0d0))then
        kalmanthreshold_temp=kalmanthreshold_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthreshold has been adjusted to ',kalmanthreshold_temp
        endif
      elseif((numbere.gt.0).and.(lshort).and.(kalmanthreshold.gt.0.0d0))then
        ztemp=kalmanthreshold_temp
        kalmanthreshold_temp=kalmanthreshold_temp/0.9d0
        kalmanthreshold_temp=min(kalmanthreshold_temp,kalmanthreshold)
        if(abs(ztemp-kalmanthreshold_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthreshold has been adjusted to ',kalmanthreshold_temp
          endif
        endif
      endif
!! check if a non-zero number of structures has been used for the force update
      if((numberf.eq.0).and.(lshort).and.(luseforces).and.(kalmanthresholdf.gt.0.0d0).and.(forcernd.gt.0.0d0))then
        kalmanthresholdf_temp=kalmanthresholdf_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthresholdf has been adjusted to ',kalmanthresholdf_temp
        endif
      elseif((numberf.gt.0).and.(lshort).and.(luseforces).and.(kalmanthresholdf.gt.0.0d0))then
        ztemp=kalmanthresholdf_temp
        kalmanthresholdf_temp=kalmanthresholdf_temp/0.9d0
        kalmanthresholdf_temp=min(kalmanthresholdf_temp,kalmanthresholdf)
        if(abs(ztemp-kalmanthresholdf_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthresholdf has been adjusted to ',kalmanthresholdf_temp
          endif
        endif
      endif
!! check if a non-zero number of charges has been used for the update
      if((numberq.eq.0).and.lelec.and.(nn_type_elec.eq.1).and.(kalmanthresholde.gt.0.0d0).and.(chargernd.gt.0.0d0))then
        kalmanthresholde_temp=kalmanthresholde_temp*0.9d0
        if(mpirank.eq.0)then
          write(ounit,*)'### WARNING ### kalmanthresholde has been adjusted to ',kalmanthresholde_temp
        endif
      elseif((numberq.gt.0).and.lelec.and.(nn_type_elec.eq.1).and.(kalmanthresholde.gt.0.0d0))then
        ztemp=kalmanthresholde_temp
        kalmanthresholde_temp=kalmanthresholde_temp/0.9d0
        kalmanthresholde_temp=min(kalmanthresholde_temp,kalmanthresholde)
        if(abs(ztemp-kalmanthresholde_temp).gt.0.000000001d0)then
          if(mpirank.eq.0)then
            write(ounit,*)'### WARNING ### kalmanthresholde has been adjusted to ',kalmanthresholde_temp
          endif
        endif
      endif
!!'
!! write optimum set of weights to files optweights.out and opweightse.out
      if((rmse_short_test.le.rmse_short_test_old).and.lshort)then
        if(mpirank.eq.0)then
          optepoch=countepoch
          optweights_short(:,:)=weights_short_atomic(:,:)
          call writeoptweights(0,nelem,&
            maxnum_weights_short_atomic,num_weights_short_atomic,&
            optweights_short)
        endif ! mpirank.eq.0
!! define new reference only if we found a new minimum rmse
        rmse_short_test_old=rmse_short_test
!! store best values obtained so far for final output
        optrmse_etot        =rmse_etot
        optrmse_etot_test   =rmse_etot_test
        optrmse_short       =rmse_short
        optrmse_short_test  =rmse_short_test
        optrmse_force_s     =rmse_force_s
        optrmse_force_s_test=rmse_force_s_test
      endif
!! electrostatic NN:
      if((rmse_charge_test.le.rmse_charge_test_old).and.lelec)then
        if(mpirank.eq.0)then
          optepoche=countepoch
          optweights_ewald(:,:)=weights_elec(:,:)
          call writeoptweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            optweights_ewald)
        endif ! mpirank.eq.0
!! define new reference only if we found a new minimum rmse
        rmse_charge_test_old    =rmse_charge_test
!! store best values obtained so far for final output
        optrmse_ewald           =rmse_ewald
        optrmse_ewald_test      =rmse_ewald_test
        optrmse_charge          =rmse_charge
        optrmse_charge_test     =rmse_charge_test
        optrmse_totalcharge     =rmse_totalcharge
        optrmse_totalcharge_test=rmse_totalcharge_test
        optrmse_force_e         =rmse_force_e
        optrmse_force_e_test    =rmse_force_e_test 
      endif
!!
      enddo ! countepoch, loop over all epochs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end loop over all epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! summarize optimum fit
      if((mpirank.eq.0).and.(nepochs.gt.0))then
        write(ounit,*)'============================================================='
        if(lshort)then
          if(optepoch.eq.0)then
            write(ounit,'(a,i5)')' No improvement of the short range fit has been obtained'
          else
            write(ounit,'(a,i5)')' Best short range fit has been obtained in epoch ',optepoch
            if(luseforces)then
              write(ounit,'(2a)')'                    --- E_short: --- ',&
                             '         --- F_short: ---'
              write(ounit,'(2a)')'                   train         test ',&
                             '       train         test'
              write(ounit,'(a11,x,4f13.6)')' OPTSHORT  ',&
                optrmse_short*tounit,optrmse_short_test*tounit,&
                optrmse_force_s*tounit,optrmse_force_s_test*tounit
            else ! no forces:
              write(ounit,'(a)')'                    --- E_short: --- '
              write(ounit,'(a)')'                   train         test '
              write(ounit,'(a11,x,2f13.6)')' OPTSHORT  ',&
                optrmse_short*tounit,optrmse_short_test*tounit
            endif ! luseforces
          endif ! optepoch.eq.0
        endif ! lshort
        if(lshort.and.lelec)then
          write(ounit,*)'-------------------------------------------------------------'
        endif
        if(lelec)then
          if(optepoche.eq.0)then
            write(ounit,'(a,i5)')' No improvement of the charge fit has been obtained'
          else
            write(ounit,'(a,i5)')' Best charge fit has been obtained in epoch ',optepoche
            if(luseforces)then
              write(ounit,'(4a)')'                    --- E_elec: --- ',&
                           '          --- charges: ---',&
                           '      --- total charges: ---',&
                           '        --- F_elec: ---'
              write(ounit,'(4a)')'                   train         test ',&
                           '       train         test',&
                           '        train         test',&
                           '        train         test'
              write(ounit,'(a11,x,8f13.6)')' OPTCHARGE ',&
              optrmse_ewald*tounit,optrmse_ewald_test*tounit,&
              optrmse_charge,optrmse_charge_test,&
              optrmse_totalcharge,optrmse_totalcharge_test,&
              optrmse_force_e*tounit,optrmse_force_e_test*tounit 
            else ! no forces:
              write(ounit,'(3a)')'                    --- E_elec: --- ',&
                           '          --- charges: ---',&
                           '      --- total charges: ---'
              write(ounit,'(3a)')'                   train         test ',&
                           '       train         test',&
                           '        train         test'
              write(ounit,'(a11,x,6f13.6)')' OPTCHARGE ',&
              optrmse_ewald*tounit,optrmse_ewald_test*tounit,&
              optrmse_charge,optrmse_charge_test,&
              optrmse_totalcharge,optrmse_totalcharge_test
            endif ! luseforces
          endif ! optepoche.eq.0
        endif ! lelec
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a)')'                      --- E_tot: --- '
        write(ounit,'(a)')'                   train         test '
        write(ounit,'(a11,x,2f13.6)')' OPTETOT   ',&
          optrmse_etot*tounit,optrmse_etot_test*tounit
      endif ! mpirank.eq.0
!!
!! report final statistics if requested:
      if(lfitstats.and.(mpirank.eq.0))then
        call writefitstat(ntrain,fitstat,fitstatf,fitstatq)
      endif
!!
      if(allocated(wconstraintidx))  deallocate(wconstraintidx)  
      if(allocated(wconstraintidxe)) deallocate(wconstraintidxe)  
!! deallocate Kalman matrices
      if(allocated(corrmatrix_list)) deallocate(corrmatrix_list)
      if(allocated(corrmatrixf_list))deallocate(corrmatrixf_list)
      if(allocated(corrmatrixe_list))deallocate(corrmatrixe_list)
      if(allocated(corrmatrixc))     deallocate(corrmatrixc) 
      deallocate(fitstat)
      deallocate(fitstatf)
      deallocate(fitstatq)
      deallocate(filenamews)
      deallocate(filenamewp)
      deallocate(filenamewe)
!!
!!    deallocate arrays of structures module
      call deallocatestructures()
!!
      return
      end
