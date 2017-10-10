!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: fit the atomic short range part and if requested atomic charges as second output nodes

!! called by:
!!
      subroutine fit_hextoff(iseed)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use symfunctions
      use nnham
      use structures
      use timings
      use basismod
!!
      implicit none
!!
!! for dimensions
      integer num_weights(1)                                    ! internal
      integer num_weightsfixed(1)                               ! internal
!! counters
      integer npoints                                           ! internal
      integer point                                             ! internal
      integer ncount                                            ! internal
      integer ndone                                             ! internal
      integer countepoch                                        ! internal
      integer i0,i1,i2                                          ! internal 
      integer optepoch                                          ! internal
      integer belowmaxenergy                                    ! internal
!!
      integer iseed                                             ! in
      integer jseed                                             ! internal
      integer kseed                                             ! internal
      integer lseed                                             ! internal
      integer mseed                                             ! internal
      integer nseed                                             ! internal
      integer oseed                                             ! internal
      integer ntrain                                            ! internal
      integer ntest                                             ! internal
      integer idx(nblock)                                       ! internal
      integer nstruct                                           ! internal
      integer n_start,n_end                                     ! internal
      integer numbere                                           ! internal
      integer numberf                                           ! internal
      integer numberq                                           ! internal
      integer ntrainatoms                                       ! internal
      integer ntestatoms                                        ! internal
      integer trainelem(1)                                      ! internal
      integer testelem(1)                                       ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidx  ! internal 
      integer, dimension(:)    , allocatable :: fitstate        ! internal 
      integer, dimension(:)  , allocatable :: pointindex        ! internal 
      integer, dimension(:)  , allocatable :: num_atoms_all     ! internal 
!! Kalman matrix dimensions:
      integer corrdim                                           ! internal
      integer corrfdim(nelem)                                   ! internal
      integer kaldim   
      integer maxkaldim
      integer maxcorrdim                                        ! internal
      integer maxcorrfdim                                       ! internal
      integer block_counter                                     ! internal
!!
!! symmetry function related arrays
      real*8 minvalue_hextoff(maxnum_funcvalues_hextoff)
      real*8 maxvalue_hextoff(maxnum_funcvalues_hextoff)
      real*8 avvalue_hextoff(maxnum_funcvalues_hextoff)
!! weights
      real*8 optweights_hextoff(maxnum_weights_hextoff,1)
      real*8 weights_hextoff_veryold(maxnum_weights_hextoff,1)        ! internal for lprintconv
      real*8 weights_hextoff_old(maxnum_weights_hextoff,1)
!! Reference data:
      real*8 chargemin(nelem)                                        ! internal
      real*8 chargemax(nelem)                                        ! internal
      real*8 dummy                                                   ! internal
      real*8 rdummy(nelem)                                           ! internal
!! RMSEs
      real*8 rmse_hextoff
      real*8 rmse_charge
      real*8 rmse_totalcharge
      real*8 rmse_etot
      real*8 rmse_elec
      real*8 rmse_hextoff_test
      real*8 rmse_hextoff_test_old
      real*8 rmse_charge_test
      real*8 rmse_charge_test_old
      real*8 rmse_totalcharge_test
      real*8 rmse_etot_test
      real*8 rmse_elec_test
!!
      real*8 rmse_hextoff_ref
      real*8 rmse_charge_ref
      real*8 rmse_totalcharge_ref
!! MADs
      real*8 mad_hextoff
      real*8 mad_hextoff_test
      real*8 mad_ewald
      real*8 mad_ewald_test
      real*8 mad_etot
      real*8 mad_etot_test
      real*8 mad_charge
      real*8 mad_charge_test
      real*8 mad_totalcharge
      real*8 mad_totalcharge_test
!! for final output
      real*8 optrmse_hextoff
      real*8 optrmse_charge
      real*8 optrmse_totalcharge
      real*8 optrmse_etot
      real*8 optrmse_elec
      real*8 optrmse_hextoff_test
      real*8 optrmse_charge_test
      real*8 optrmse_totalcharge_test
      real*8 optrmse_etot_test
      real*8 optrmse_elec_test

!! Kalman filter 
      real*8 kalmanthreshold_temp                             ! internal
      real*8 kalmanthresholdf_temp                            ! internal
      real*8, dimension(:,:), allocatable :: corrmatrix_list
      real*8, dimension(:,:), allocatable :: corrmatrixf_list
!! arrays for convergence vector
      real*8 wshift(nelem)                                    ! internal
      real*8 wshift2(nelem)                                   ! internal
      real*8 wshiftold(nelem)                                 ! internal
      real*8 convvec(nelem,2,3)                               ! internal
      real*8 convtemp1                                        ! internal
      real*8 convtemp2                                        ! internal
      real*8 convtemp3                                        ! internal
!! miscellaneous
      real*8, dimension(:), allocatable :: hextoff_min             ! internal
      real*8, dimension(:), allocatable :: hextoff_max             ! internal
      real*8 tounit                                           ! internal
      real*8 avcharge(nelem)                                  ! internal
      real*8 stddevcharge(nelem)                              ! internal
      real*8, dimension(:), allocatable :: hextoff_av              ! internal
      real*8, dimension(:), allocatable :: hextoff_stddev          ! internal
      character*200 :: filenametemp                           ! internal

      integer*8 matrixsize                                       ! internal
!!
!!
!!=========================================================
!! some checks for electrostatics  !! But of course we don't use this for NNTB training
!!=========================================================
      if(lelec)then
        write(ounit,*) 'ERROR: fitting electrostatics at the same time as NNTB is not supported'
        stop
      endif ! lelec
      if(lshort)then
        write(ounit,*) 'ERROR: fitting short range at the same time as NNTB is not supported'
        stop
      endif ! lelec
!!
      matrixsize=num_basis(elementindex(hextoff_training_triplet(1)))*num_basis(elementindex(hextoff_training_triplet(2)))
      allocate(hextoff_min(matrixsize))
      allocate(hextoff_max(matrixsize))
      allocate(hextoff_av(matrixsize))
      allocate(hextoff_stddev(matrixsize))
!!
!!=========================================================
!! allocate arrays of structures module
!!=========================================================
      call allocatestructures()
!!
!!=========================================================
!! initialization of local Kalman filter parameters
!!=========================================================
      kalmanthreshold_temp   = kalmanthreshold
!!=========================================================
!! initialization of scaling parameters
!!=========================================================
      minvalue_hextoff(:)          = 0.0d0
      maxvalue_hextoff(:)          = 0.0d0
      avvalue_hextoff(:)           = 0.0d0
!!
!!=========================================================
!! set filenames for weights files
!!=========================================================
      call getweightfilenames()
!!
!!=========================================================
!! initialization of reference data arrays
!!=========================================================
      lattice_list(:,:,:)    = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
!!!!!!!!! need to sort this one      hextoff_list(:,:)      = 0.0d0
!!
!!
!!=========================================================
!! initialization of reference RMSEs
!!=========================================================
      rmse_hextoff_ref      =0.0d0
!!
!!=========================================================
!! initialization of optrmses
!!=========================================================
      optrmse_hextoff            =0.0d0
      optrmse_hextoff_test       =0.0d0
!!
!!=========================================================
!! generating local copies of seed for various purposes
!! iseed is just used for random order training below
!!=========================================================
      jseed                  = iseed ! for initial short range weights 
      kseed                  = iseed ! for forcernd 
      lseed                  = iseed ! for energyrnd
      mseed                  = iseed ! for chargernd 
      nseed                  = iseed ! for initial electrostatic weights
      oseed                  = iseed ! for pointindex 
!!
!!=========================================================
!! initializations of counters
!!=========================================================
      countepoch             = 0
      point                  = 0
      npoints                = 0
      numbere                = 0
      optepoch               = 0
      block_counter          = 1
!!
!!=========================================================
!! initialization of timing variables
!!=========================================================
      timegeterror           = 0.0d0
      timedeshortdw          = 0.0d0
      timeefitting           = 0.0d0
      timeeupdate            = 0.0d0
      timeeerror             = 0.0d0
      timeqerror             = 0.0d0
      timequpdate            = 0.0d0
      timeio                 = 0.0d0
      timemix                = 0.0d0
      timeshortfit           = 0.0d0
      daygeterror            = 0
      dayio                  = 0
      daymix                 = 0
      dayshortfit            = 0
!!
!!=========================================================
!! further initializations
!!=========================================================
      convvec(:,:,:)         = 0.0d0
      wshift(:)              = 0.0d0
      wshift2(:)             = 0.0d0
!!
!!=========================================================
!! determine starting time
!!=========================================================
      call abstime(timeepochinistart,dayepochini)
!!
!!=========================================================
!!=========================================================
!! end of initializations
!!=========================================================
!!=========================================================
!!
!!=========================================================
!! get information on fixed weights from input.nn file
!!=========================================================
      i0 = tripletindex(hextoff_training_triplet(1),hextoff_training_triplet(2),hextoff_training_triplet(3))
      allocate(wconstraintidx(maxnum_weights_hextoff,1)) !! again we only train and so use one set of weights at a time CMH 
      num_weights = num_weights_hextoff  
!!
!!=========================================================
!! get dimensions for the Kalman filter
!!=========================================================
      if(optmodehextoff.eq.1)then
        call getkaldims_hextoff(ntriplets,&
          kaldim,corrdim,&
          maxkaldim,maxcorrdim,&
          num_weights)
      endif
!!
!!=========================================================
!! allocate Kalman matrices
!!=========================================================
      allocate(corrmatrix_list(maxcorrdim,1))
!!
!!=========================================================
!! get the values for the Kalman matrices and parameters (new or from file)
!!=========================================================
      call getkalmanmatrices_hextoff(1,&
        iseed,kseed,lseed,mseed,&
        corrdim,&
        maxcorrdim,num_weights,&
        corrmatrix_list,kalmanlambdahextoff)
!!=========================================================
!! distribute random seeds possibly read from file
!!=========================================================
      call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(kseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mseed,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!=========================================================
!! count the number of training and test points
!!=========================================================
       
      if(mpirank.eq.0)then
        call countpoints(1,1,&
          ntrain,ntest,ntrainatoms,ntestatoms,trainelem,testelem)
      endif ! mpirank.eq.0
      call mpi_bcast(ntrain,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntest,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntrainatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntestatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(trainelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(testelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!=========================================================
!! allocate pointindex array for random order of structures 
!!=========================================================
      allocate(pointindex(ntrain))
!!
!!=========================================================
!! get array num_atoms_all, contains number of atoms in each training point
!!=========================================================
      allocate(num_atoms_all(ntrain))
      num_atoms_all = 3
!!      call getnumatomsall(ntrain,num_atoms_all)
!!
!!=========================================================
!! allocate arrays for fitting statistics output
!!=========================================================
!!
!! CMH No idea which of these or the sizes we need
!!
      allocate(fitstate(ntrain))
      fitstate(:)=0
!!
!!=========================================================
!! if growth mode is not used make sure the full training set is used
!!=========================================================
      if(.not.lgrowth)then
        ngrowth=ntrain
      endif
!!
!!=========================================================
!! get maxcutoff of short range symmetry functions
!!=========================================================
      maxcutoff_hextoff = 0.0d0
!!      do i2=1,ntriplets
        i2 = 1
        do i1=1,num_funcvalues_hextoff(i2)
          maxcutoff_hextoff=max(maxcutoff_hextoff,funccutoff_hextoff(i1,i2))
        enddo ! i1
!!      enddo
!!
      if(mpirank.eq.0)then
!!=========================================================
!! get scaling data for the short range symmetry functions
!!=========================================================
        if(luseoldscaling)then
          call readscale_nntb(matrixsize,&
            hextoff_training_triplet(1),hextoff_training_triplet(2),&
            hextoff_training_triplet(3),&
            maxnum_funcvalues_hextoff,num_funcvalues_hextoff(1),&
            minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
            hextoff_min, hextoff_max)
        else ! luseoldscaling
          write(filenametemp,'(A,I3.3,A,I3.3,A,I3.3,A)') &
          'function_hextoff.',hextoff_training_triplet(1),'.',hextoff_training_triplet(2),'.',hextoff_training_triplet(3),'.data'
          open(symhextoffunit,file=filenametemp,form='formatted',status='old')
          rewind(symhextoffunit) !'
          call getscale_hextoff(matrixsize,&
            maxnum_funcvalues_hextoff,num_funcvalues_hextoff(1),&
            ntrain,minvalue_hextoff,&
            maxvalue_hextoff,avvalue_hextoff)
          close(symhextoffunit)
        endif ! luseoldscaling
!!
!!=========================================================
!! get statistics of energies in training set (just for process 0)
!!=========================================================
        call gethamiltonianstatistics(1,hextoff_av,hextoff_stddev,&
          hextoff_min,hextoff_max,matrixsize)
!!
!!=========================================================
!! write short range scaling factors to file scaling.data
!!=========================================================
        if(.not.luseoldscaling)then
          call writescale_nntb(matrixsize,&
            hextoff_training_triplet(1),hextoff_training_triplet(2),&
            hextoff_training_triplet(3),&
            maxnum_funcvalues_hextoff,num_funcvalues_hextoff(1),&
            minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
            hextoff_min,hextoff_max)
        endif ! luseoldscaling
      endif ! mpirank.eq.0
!!
!!=========================================================
!! distribute scaling.data to all processes
!!=========================================================
      call mpi_bcast(minvalue_hextoff,maxnum_funcvalues_hextoff,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(maxvalue_hextoff,maxnum_funcvalues_hextoff,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(avvalue_hextoff,maxnum_funcvalues_hextoff,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(hextoff_stddev,maxnum_funcvalues_hextoff,mpi_real8,0,mpi_comm_world,mpierror)
!!
!!=========================================================
!! write information about training data to output 
!!=========================================================
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'number of training points ',ntrain !'
        write(ounit,*)'number of testing points  ',ntest
      endif ! mpirank.eq.0 '
!!
!!=========================================================
!! preparation of initial weight parameters'
!!=========================================================
!!
!!=========================================================
!! first initialize all weights randomly no matter if we restart a fit or not
!!=========================================================
      if(mpirank.eq.0)then
        call initialweights(1,&
          maxnum_weights_hextoff,num_weights_hextoff,&
          maxnum_layers_hextoff,num_layers_hextoff,&
          windex_hextoff,nodes_hextoff,&
          jseed,nseed,weights_hextoff)
!!
!!=========================================================
!! if weights according to Nguyen Widrow are requested for short range NN, overwrite random weights:
!!=========================================================
        if(lnwweights)then
          write(*,*) 'Warning: ngugen widrow weights not checked for nntb_flag(3)'
          call nguyenwidrowweights(1,maxnum_layers_hextoff,&
            nodes_hextoff,jseed,windex_hextoff,&
            maxnodes_hextoff,maxnum_weights_hextoff,num_weights_hextoff,&
            weights_hextoff,weights_min,weights_max,actfunc_hextoff)
        endif
!!
!!=========================================================
!! If systematic weights are requested, overwrite random weights for short range NN
!!=========================================================
        if(lsysweights)then
          write(*,*) 'Warning: systematic weights not checked for nntb_flag(3)'
          call systematicweights(1,maxnum_layers_hextoff,&
            num_layers_hextoff,nodes_hextoff,&
            maxnum_weights_hextoff,&
            weights_hextoff,weights_min,weights_max)
        endif
!!
!!=========================================================
!! if requested overwrite hextoff weights with weights from weights.XXX.data
!!=========================================================
        luseoldweights_hextoff = .true.
        if(luseoldweights_hextoff)then
          call readweights(3,1,&
            maxnum_weights_hextoff,num_weights_hextoff,&
            weights_hextoff)
        endif
      endif ! mpirank.eq.0
!!
!!=========================================================
!! distribute weights_hextoff to all processes
!!=========================================================
      if(lnntb)then
        call mpi_bcast(weights_hextoff,1*maxnum_weights_hextoff,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!=========================================================
!! write initial weights to file
!!=========================================================
      if(mpirank.eq.0)then
        call writeweights(3,1,maxnum_weights_hextoff,&
          maxnum_layers_hextoff,num_layers_hextoff,&
          nodes_hextoff,weights_hextoff)
      endif ! mpirank.eq.0
!!
!!=========================================================
!! keep copy of weights for the calculation of the weight change wshift
!!=========================================================
      weights_hextoff_veryold(:,:)=0.0d0
      weights_hextoff_old(:,:)    =weights_hextoff(:,:)
!!
!!=========================================================
!! end of preparation of initial weight parameters'
!!=========================================================
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!=========================================================
!! determine and print initialization time
!!=========================================================
      call abstime(timeepochiniend,dayepochini)
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.2)')&
          ' initialization time (min):',(timeepochiniend-timeepochinistart)/60.d0
      endif ! mpirank.eq.0
!!'
!!=========================================================
!! terminate RuNNer here if just the initialization but no fit and no error calculation were requested
!!=========================================================
      if(linionly)then
        write(ounit,*)'Initialization done - terminating RuNNer as requested'
        write(ounit,*)'-------------------------------------------------------------'
        return !'
      endif
!!
!!=========================================================
!! set time to zero for epoch 0 time measurement
!!=========================================================
      timeepochstart =0.0d0
      dayepoch       =0
      call abstime(timeepochstart,dayepoch)
!!
!!=========================================================
!! initialize RMSEs
!!=========================================================
      rmse_hextoff          =0.0d0
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_etot             =0.0d0
      rmse_elec             =0.0d0
      rmse_hextoff_test     =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_etot_test        =0.0d0
      rmse_elec_test        =0.0d0
!!
!!=========================================================
!! initialize MADs
!!=========================================================
      mad_hextoff           =0.0d0
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_etot              =0.0d0
      mad_ewald             =0.0d0
      mad_hextoff_test      =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_etot_test         =0.0d0
      mad_ewald_test        =0.0d0
!!
      maxerror_hextoff_train  = 0.0d0 
      maxerror_hextoff_test   = 0.0d0 
      imaxerror_hextoff_train = 0 
      imaxerror_hextoff_test  = 0 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!! (rmse_short,rmse_charge,rmse_totalcharge,rmse_elec,rmse_etot,rmse_force_s,rmse_force_e)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call abstime(timegeterrorstart,daygeterror)
      call geterror_hextoff(0,countepoch,ntrain,&
        imaxerror_hextoff_train,minvalue_hextoff,&
        maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff,mad_hextoff,maxerror_hextoff_train)
!!
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of the first training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! get new references for adaptive Kalman filter
      if(lnntb.and.nntb_flag(3))  rmse_hextoff_ref    =rmse_hextoff
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the first testing error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror_hextoff(1,countepoch,ntest,&
        imaxerror_hextoff_test,minvalue_hextoff,&
        maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff_test,mad_hextoff_test,maxerror_hextoff_test)
      call abstime(timegeterrorend,daygeterror)
      timegeterror=timegeterrorend-timegeterrorstart
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of first test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! set first references for determination of optweights files
      rmse_hextoff_test_old  = rmse_hextoff_test
!!
      call abstime(timeepochend,dayepoch)
      timeepoch=timeepochend-timeepochstart

!!
!! write RMSE headers
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'RMSEs (energies: Ha/atom, charges: e, forces: Ha/Bohr):'
        write(ounit,'(6a)')'                      --- E_hextoff: --- ',&
                         ' - time -'
        write(ounit,'(6a)')'                          /atom      ',&
                         '   min'
        write(ounit,'(5a)')'       epoch         train         test'
!!
!! CAUTION: If the training data is constructed with short and ewald and
!! then lshort is set false to fit the charges, then rmse_etot and rmse_elec
!! cannot be the same because they have different DFT reference energies
!! (for etot it is the total DFT energy and for the electrostatic energy
!! it is only the DFT ewald energy, from the NN we in both cases have just the
!! ewald energy)
!!
!! write RMSEs'
        write(ounit,'(a8,i5,x,10f13.6,f8.2)') ' ENERGY ',countepoch,&
            rmse_hextoff,rmse_hextoff_test,&
            timeepoch/60.d0
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)') ' MADE   ',countepoch,&
            mad_hextoff,mad_hextoff_test
        endif
!!
!! Convergence - might need this
!! 
!        if(lprintconv)then
!          do i1=1,nelem
!            write(ounit,'(a9,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',element(i1),&
!              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshift(i1),wshift2(i1)
!            convvec(i1,:,3)=convvec(i1,:,2)
!            convvec(i1,:,2)=convvec(i1,:,1)
!            convvec(i1,:,1)=0.0d0
!          enddo ! i1
!        endif
!!
!!
!!
        if(lweightanalysis)then
          call analyzeweights(1,countepoch,1,&
            maxnum_weights_hextoff,num_weights_hextoff,weights_hextoff)
        endif ! lweightanalysis
      endif ! mpirank.eq.0

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all training epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do countepoch=1,nepochs
        write(ounit,*)'-------------------------------------------------------------------------------'
        if(mpirank.eq.0) write(debugunit,'(a6,i6)')'epoch ',countepoch
        call abstime(timeepochstart,dayepoch) !'
!!
!! epoch-specific initializations
!!
        ndone                 = 0
        ncount = (int((countepoch-1)/growthstep)+1)*ngrowth    
        ncount = min(ncount,ntrain)
        point                 = 0
        numbere               = 0
        numberf               = 0
        numberq               = 0
!!
!! get array pointindex for randomly mixing points in training
        call abstime(timemixstart,daymix)
        call getpointindex(ncount,oseed,pointindex)
        call abstime(timemixend,daymix)
        timemix=timemixend-timemixstart
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
!! reinitialize correlation matrices of Kalman filter if requested
      lresetkalman = .true. 
      if(lresetkalman)then
        do i1=1,1
          if(lnntb.and.nntb_flag(3).and.(optmodehextoff.eq.1))then
            call initialcorrmatrix(num_weights_hextoff(i1),corrmatrix_list(1,i1))
          endif
        enddo
      endif ! lresetkalman
!!
!! get index for random order of training points => idx
      call getidx(npoints,iseed,idx)
!!
!! read npoint data sets      
!! do all file reading for training data here at one place to allow for parallelization
!!--------------------------------------------------------------------------------------
!!
      call abstime(timeiostart,dayio)
      write(*,*) npoints, ntrain
      if(mpirank.eq.0)then
        if(lnntb.and.nntb_flag(3))then
          call readfunctions_hextoff_mixed(0,npoints,1,ntrain,&
            pointindex,block_counter,max_num_atoms,maxnum_funcvalues_hextoff,&
            num_funcvalues_hextoff,&
            symfunction_hextoff_list)
        endif ! lnntb
!!
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! is needed for lshort (force fitting) and lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
        call getstructures_mixed(npoints,&
          ntrain,block_counter,pointindex,num_atoms_all)
!!
      endif ! mpirank.eq.0
      call abstime(timeioend,dayio)
      timeio=timeio+timeioend-timeiostart
!!
!! end of file reading for training data 
!!--------------------------------------------------------------------
!!
!! determine which structures of this block should be scaled by this process
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! start optimization 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call abstime(timeshortfitstart,dayshortfit)
!!
!! scale symmetry functions for the short-range interaction
      call scalesym_hextoff(npoints,&
        maxnum_funcvalues_hextoff,num_funcvalues_hextoff,&
        symfunction_hextoff_list,&
        minvalue_hextoff,maxvalue_hextoff,avvalue_hextoff,&
        scmin_hextoff,scmax_hextoff)
!!
!! update using energies and forces together
!!      call optimize_hextoff(npoints,point,idx,&
!!        kaldim,maxcorrdim,maxcorrfdim,corrdim,corrfdim,countepoch,&
!!        ndone,num_weights_hextoff,&
!!        wconstraintidx,numbere,numberf,numberq,&
!!        lseed,ntrain,&
!!        fitstate,&
!fitstatf,fitstatq,
!!        kseed,&
!!        kalmanthreshold_temp,kalmanthresholdf_temp,&
!!        rmse_hextoff_ref,&
!!        corrmatrix_list,corrmatrixf_list,&
!!        minvalue_hextoff,maxvalue_hextoff)
!!
      STOP
      call abstime(timeshortfitend,dayshortfit)
      timeshortfit=timeshortfit+timeshortfitend-timeshortfitstart
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end optimization of the short range part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! write temporary weights if requested
      if(mpirank.eq.0)then
        if(lwritetmpweights)then
          call writetmpweights(0,nelem,&
            maxnum_weights_hextoff,num_weights_hextoff,&
            weights_hextoff)
        endif
      endif ! mpirank.eq.0
!! 
!! count finished structures
      ndone        =ndone+npoints
      block_counter=block_counter+npoints
!!
      if(ncount.gt.0) goto 11 ! do next block of points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end block-wise loop over training points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! write Kalman data for restart if requested
      if(lsavekalman)then
        if((optmodee.eq.1).or.(optmodef.eq.1))then
          call writekalman_short(nelem,&
            iseed,kseed,lseed,mseed,&
            maxcorrdim,kaldim,&
            num_weights_hextoff,&
            kalmanlambda,corrmatrix_list)
        endif ! optmode
      endif !lsavekalman
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! write final weights of this epoch to files
      if(mpirank.eq.0)then
        if(mod(countepoch,iwriteweight).eq.0) then
          call getfilenames(countepoch)
          call writeweights(0,nelem,&
            maxnum_weights_hextoff,&
            maxnum_layers_hextoff,num_layers_hextoff,&
            nodes_hextoff,weights_hextoff)
        endif
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! initialize RMSEs for next epoch
      rmse_hextoff            =0.0d0
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_etot             =0.0d0
      rmse_elec             =0.0d0
      rmse_hextoff_test       =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_etot_test        =0.0d0
      rmse_elec_test        =0.0d0
!!
!! initialize MADs for next epoch
      mad_hextoff             =0.0d0
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_etot              =0.0d0
      mad_ewald             =0.0d0
      mad_hextoff_test        =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_etot_test         =0.0d0
      mad_ewald_test        =0.0d0
!!
      maxerror_hextoff_train  = 0.0d0 
      maxerror_hextoff_test   = 0.0d0 
      imaxerror_hextoff_train = 0 
      imaxerror_hextoff_test  = 0 
      maxerror_elec_train    = 0.0d0 
      maxerror_elec_test     = 0.0d0 
      imaxerror_elec_train   = 0 
      imaxerror_elec_test    = 0 
      maxerror_etot_train    = 0.0d0 
      maxerror_etot_test     = 0.0d0 
      imaxerror_etot_train   = 0 
      imaxerror_etot_test    = 0 
      block_counter          = 1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate this epoch's training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
      timegeterrorstart=0.0d0
      daygeterror      =0
      call abstime(timegeterrorstart,daygeterror)
!!
      call geterror_hextoff(0,countepoch,ntrain,&
        imaxerror_hextoff_train,minvalue_hextoff,&
        maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff,mad_hextoff,maxerror_hextoff_train)
!!
!! get new references for adaptive Kalman filter (this has nothing to do with optweights files)
      if(lnntb.and.nntb_flag(3)) rmse_hextoff_ref      =rmse_hextoff
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of the epoch's training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the epoch's testing error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror_hextoff(1,countepoch,ntest,&
        imaxerror_hextoff_test,minvalue_hextoff,&
        maxvalue_hextoff,avvalue_hextoff,&
        rmse_hextoff_test,mad_hextoff_test,maxerror_hextoff_test)

      write(*,*) rmse_hextoff, rmse_hextoff_test
      STOP
!!
      call abstime(timegeterrorend,daygeterror)
      timegeterror=timegeterrorend-timegeterrorstart
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of epoch's test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call abstime(timeepochend,dayepoch)
      timeepoch=timeepochend-timeepochstart
!!
!! calculate change of the short range weights for INFO and CONVVEC output
      if(lshort)then
        call getwshift(nelem,maxnum_weights_hextoff,num_weights_hextoff,&
          weights_hextoff,weights_hextoff_old,weights_hextoff_veryold,&
          wshift,wshift2)
      endif ! lshort
!!
!! keep weights
      weights_hextoff_veryold(:,:)=weights_hextoff_old(:,:)
      weights_hextoff_old(:,:)    =weights_hextoff(:,:)
!!
!! write RMSEs
      if(mpirank.eq.0)then
        write(ounit,'(a8,i5,x,10f13.6,f8.2)')' ENERGY ', countepoch,&
          rmse_hextoff*tounit,rmse_hextoff_test*tounit,&
          rmse_elec*tounit,rmse_elec_test*tounit,&
          rmse_etot*tounit,rmse_etot_test*tounit,&
          rmse_charge,rmse_charge_test,&
          rmse_totalcharge,rmse_totalcharge_test,&
          timeepoch/60.d0
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)')' MADE   ', countepoch,&
            mad_hextoff*tounit,mad_hextoff_test*tounit,&
            mad_ewald*tounit,mad_ewald_test*tounit,&
            mad_etot*tounit,mad_etot_test*tounit,&
            mad_charge,mad_charge_test,&
            mad_totalcharge,mad_totalcharge_test
        endif
        write(ounit,'(a,i5,3i10,2f20.15)')' INFORMATION USED FOR UPDATE (E,F,Q) ',&
          countepoch,numbere,numberf,numberq !'
        if(lprintconv)then
          do i1=1,nelem
            if(countepoch.eq.1)then
              convvec(i1,1,1)=wshift(i1)
              convvec(i1,2,1)=0.0d0
            else
              convtemp1=atan((convvec(i1,2,2)-convvec(i1,2,3))/(convvec(i1,1,2)-convvec(i1,1,3)))
              convtemp2=wshiftold(i1) !abs((convvec(i1,1,3)-convvec(i1,1,2))/dcos(convtemp1))
              convtemp3=acos((-wshift(i1)**2+convtemp2**2+wshift2(i1)**2)/(2.d0*convtemp2*wshift2(i1)))
              convtemp1=convtemp1+convtemp3
              convvec(i1,1,1)=convvec(i1,1,3)+wshift2(i1)*dcos(convtemp1)
              convvec(i1,2,1)=convvec(i1,2,3)+wshift2(i1)*dsin(convtemp1)
            endif
            write(ounit,'(a9,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',element(i1),&
              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshift(i1),wshift2(i1)
            convvec(i1,:,3)=convvec(i1,:,2)
            convvec(i1,:,2)=convvec(i1,:,1)
            convvec(i1,:,1)=0.0d0
            wshiftold(i1) =wshift(i1)
          enddo ! i1
        endif ! lprintconv
        if(lweightanalysis)then
          call analyzeweights(1,countepoch,nelem,&
            maxnum_weights_hextoff,num_weights_hextoff,weights_hextoff)
        endif ! lweightanalysis
!! print detailed timing for epoch if requested
        if(lfinetimeepoch)then
          call writeepochtime(countepoch)
        endif ! lfinetimeepoch
      endif ! mpirank.eq.0
!!
      call adjustkalman_short(numbere,numberf,&
        kalmanthreshold_temp,kalmanthresholdf_temp)
!!'
!! write optimum set of weights to files optweights.out and opweightse.out
!! FIXME: WE NEED BETTER CRITERION TO DEFINE OPTIMUM EPOCH NOW
      if((rmse_hextoff_test.le.rmse_hextoff_test_old).and.lshort)then
        if(mpirank.eq.0)then
          optepoch=countepoch
          optweights_hextoff(:,:)=weights_hextoff(:,:)
          call writeoptweights(0,nelem,&
            maxnum_weights_hextoff,num_weights_hextoff,&
            optweights_hextoff)
        endif ! mpirank.eq.0
!! define new reference only if we found a new minimum rmse
        rmse_hextoff_test_old=rmse_hextoff_test
!! store best values obtained so far for final output
        optrmse_etot        =rmse_etot
        optrmse_etot_test   =rmse_etot_test
        optrmse_hextoff     =rmse_hextoff
        optrmse_hextoff_test =rmse_hextoff_test
        optrmse_elec            =rmse_elec
        optrmse_elec_test       =rmse_elec_test
        optrmse_charge          =rmse_charge
        optrmse_charge_test     =rmse_charge_test
        optrmse_totalcharge     =rmse_totalcharge
        optrmse_totalcharge_test=rmse_totalcharge_test
      endif
!!
      enddo ! countepoch, loop over all epochs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end loop over all epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! summarize optimum fit
      if((mpirank.eq.0).and.(nepochs.gt.0))then
!! FIXME: This routine has not been adjusted
!        call writeoptfit_atomic(optepoch,&
!          tounit,optrmse_hextoff,optrmse_hextoff,optrmse_force_hextoff,&
!          optrmse_force_hextoff_test,optrmse_elec,optrmse_elec_test,&
!          optrmse_charge,optrmse_charge_test,optrmse_totalcharge,&
!          optrmse_totalcharge_test,optrmse_force_e,optrmse_force_e_test,&
!          optrmse_etot,optrmse_etot_test)
      endif ! mpirank.eq.0
!!
      if(allocated(wconstraintidx))  deallocate(wconstraintidx)  
!! deallocate Kalman matrices
      if(allocated(corrmatrix_list)) deallocate(corrmatrix_list)
      if(allocated(fitstate))        deallocate(fitstate)
      if(allocated(pointindex))      deallocate(pointindex)
      if(allocated(num_atoms_all))   deallocate(num_atoms_all)
      deallocate(filenamews)
!!
!!    deallocate arrays of structures module
      call deallocatestructures()
!!
!!
      return
      end
