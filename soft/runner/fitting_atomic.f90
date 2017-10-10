!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: fit the atomic short range part and if requested atomic charges as second output nodes

!! called by:
!!
      subroutine fitting_atomic(iseed)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use symfunctions
      use nnshort_atomic
      use structures
      use timings
!!
      implicit none
!!
!! for dimensions
      integer num_weightsshortfree(nelem)                       ! internal
      integer num_weightsshortfixed(nelem)                      ! internal
!! counters
      integer npoints                                           ! internal
      integer point                                             ! internal
      integer ncount                                            ! internal
      integer ndone                                             ! internal
      integer countepoch                                        ! internal
      integer i1,i2                                             ! internal 
      integer optepoch                                          ! internal
      integer belowmaxenergy                                    ! internal
      integer belowfmax(nelem)                                  ! internal
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
      integer trainelem(nelem)                                  ! internal
      integer testelem(nelem)                                   ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidx  ! internal 
      integer, dimension(:)    , allocatable :: fitstate        ! internal 
      integer, dimension(:,:,:), allocatable :: fitstatf        ! internal 
      integer, dimension(:,:)  , allocatable :: fitstatq        ! internal 
      integer, dimension(:)  , allocatable :: pointindex        ! internal 
      integer, dimension(:)  , allocatable :: num_atoms_all     ! internal 
!! Kalman matrix dimensions:
      integer corrdim(nelem)                                    ! internal
      integer corrfdim(nelem)                                   ! internal
      integer kaldim(nelem)
      integer maxkaldim
      integer maxcorrdim                                        ! internal
      integer maxcorrfdim                                       ! internal
      integer block_counter                                     ! internal
!!
!! symmetry function related arrays
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)
!! weights
      real*8 optweights_short(maxnum_weights_short_atomic,nelem)
      real*8 weights_short_veryold(maxnum_weights_short_atomic,nelem)        ! internal for lprintconv
      real*8 weights_short_old(maxnum_weights_short_atomic,nelem)
!! Reference data:
      real*8 chargemin(nelem)                                        ! internal
      real*8 chargemax(nelem)                                        ! internal
      real*8 dummy                                                   ! internal
      real*8 rdummy(nelem)                                           ! internal
      real*8 fmin(nelem)                                             ! internal
      real*8 fmax(nelem)                                             ! internal
      real*8 fvecmin(nelem)                                             ! internal
      real*8 fvecmax(nelem)                                             ! internal
!! RMSEs
      real*8 rmse_short
      real*8 rmse_charge
      real*8 rmse_totalcharge
      real*8 rmse_etot
      real*8 rmse_elec
      real*8 rmse_short_test
      real*8 rmse_short_test_old
      real*8 rmse_charge_test
      real*8 rmse_charge_test_old
      real*8 rmse_totalcharge_test
      real*8 rmse_etot_test
      real*8 rmse_elec_test
      real*8 rmse_force_s
      real*8 rmse_force_s_test
      real*8 rmse_force_t
      real*8 rmse_force_t_test
      real*8 rmse_force_e
      real*8 rmse_force_e_test
!!
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
      real*8 optrmse_elec
      real*8 optrmse_short_test
      real*8 optrmse_charge_test
      real*8 optrmse_totalcharge_test
      real*8 optrmse_etot_test
      real*8 optrmse_elec_test
      real*8 optrmse_force_s
      real*8 optrmse_force_s_test
      real*8 optrmse_force_e
      real*8 optrmse_force_e_test

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
      real*8 eshortmin                                        ! internal
      real*8 eshortmax                                        ! internal
      real*8 tounit                                           ! internal
      real*8 avcharge(nelem)                                  ! internal
      real*8 stddevcharge(nelem)                              ! internal
      real*8 eshortav                                         ! internal
      real*8 eshortstddev                                     ! internal
      real*8, dimension(:,:)  , allocatable :: sens           ! internal
!!
!!
!!=========================================================
!! initializations
!!=========================================================
      if(.not.lshort)then
        write(ounit,*)'ERROR: without lshort fitting_atomic does not make sense'
        stop !'
      endif
!!
!!=========================================================
!! some checks for electrostatics 
!!=========================================================
      if(lelec)then
        if(nn_type_elec.eq.1)then
          write(ounit,*)'ERROR: fitting_atomic is not the right subroutine here'
          stop !'
!! do not fit charges if they are fixed anyway
        elseif(nn_type_elec.eq.3)then
          lelec=.false.
        endif
      endif ! lelec
!!
!!=========================================================
!! allocation of sensitivity array 
!!=========================================================
      if(lsens)then
!! FIXME: sensitivity needs to take into account several output nodes
        allocate(sens(nelem,maxnum_funcvalues_short_atomic))
        sens(:,:)=0.0d0
      endif
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
      kalmanthresholdf_temp  = kalmanthresholdf
!!=========================================================
!! initialization of scaling parameters
!!=========================================================
      minvalue_short_atomic(:,:)          = 0.0d0
      maxvalue_short_atomic(:,:)          = 0.0d0
      avvalue_short_atomic(:,:)           = 0.0d0
!!
!!=========================================================
!! set filenames for weights files
!!=========================================================
      allocate(filenamews(nelem))
      call getweightfilenames()
!!
!!=========================================================
!! initialization of reference data arrays
!!=========================================================
      totalenergy_list(:)    = 0.0d0
      shortenergy_list(:)    = 0.0d0
      elecenergy_list(:)     = 0.0d0
      atomcharge_list(:,:)   = 0.0d0
      totalcharge_list(:)    = 0.0d0
      shortforce_list(:,:,:) = 0.0d0
      elecforce_list(:,:,:)  = 0.0d0
      totalforce_list(:,:,:) = 0.0d0
      lattice_list(:,:,:)    = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
!!
!!=========================================================
!! setting the energy and force unit converter
!!=========================================================
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV 
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0   !  stay with Ha 
      endif
!!
!!=========================================================
!! initialization of reference RMSEs
!!=========================================================
      rmse_short_ref      =0.0d0
      rmse_charge_ref     =0.0d0
      rmse_totalcharge_ref=0.0d0
      rmse_force_s_ref    =0.0d0
!!
!!=========================================================
!! initialization of optrmses
!!=========================================================
      optrmse_short              =0.0d0
      optrmse_charge             =0.0d0
      optrmse_totalcharge        =0.0d0
      optrmse_etot               =0.0d0
      optrmse_elec               =0.0d0
      optrmse_short_test         =0.0d0
      optrmse_charge_test        =0.0d0
      optrmse_totalcharge_test   =0.0d0
      optrmse_etot_test          =0.0d0 
      optrmse_elec_test          =0.0d0
      optrmse_force_s            =0.0d0
      optrmse_force_s_test       =0.0d0
      optrmse_force_e            =0.0d0
      optrmse_force_e_test       =0.0d0
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
      numberf                = 0
      numberq                = 0
      optepoch               = 0
      block_counter          = 1
!!
!!=========================================================
!! initialization of timing variables
!!=========================================================
      timegeterror           = 0.0d0
      timedeshortdw          = 0.0d0
      timedeshortdwrepeat    = 0.0d0
      timedfshortdw          = 0.0d0
      timeefitting           = 0.0d0
      timeefittingrepeat     = 0.0d0
      timeffitting           = 0.0d0
      timeeupdate            = 0.0d0
      timeeupdaterepeat      = 0.0d0
      timefupdate            = 0.0d0
      timeeerror             = 0.0d0
      timeeerrorrepeat       = 0.0d0
      timeferror             = 0.0d0
      timefsym               = 0.0d0
      timeqerror             = 0.0d0
      timedqdw               = 0.0d0
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
      chargemin(:)           = 100000.d0
      chargemax(:)           = -100000.d0
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
      allocate(wconstraintidx(maxnum_weights_short_atomic,nelem))  
      if(lshort)then
        call getwconstraintidx(0,nelem,windex_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic,maxnum_weights_short_atomic,&
          num_weights_short_atomic,num_weightsshortfree,num_weightsshortfixed,&
          nodes_short_atomic,wconstraintidx)
      endif ! lshort
      call mpi_bcast(num_weightsshortfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsshortfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!=========================================================
!! get dimensions for the Kalman filter
!!=========================================================
      if((optmodee.eq.1).or.(optmodef.eq.1))then
        call getkaldims_short(nelem,&
          kaldim,corrdim,corrfdim,&
          maxkaldim,maxcorrdim,maxcorrfdim,&
          num_weightsshortfree)
      endif
!!
!!=========================================================
!! allocate Kalman matrices
!!=========================================================
      allocate(corrmatrix_list(maxcorrdim,nelem))
      allocate(corrmatrixf_list(maxcorrfdim,nelem))
!!
!!=========================================================
!! get the values for the Kalman matrices and parameters (new or from file)
!!=========================================================
      call getkalmanmatrices_short(nelem,&
        iseed,kseed,lseed,mseed,&
        corrdim,corrfdim,&
        maxcorrdim,maxcorrfdim,num_weightsshortfree,&
        corrmatrix_list,corrmatrixf_list,kalmanlambda)
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
        call countpoints(1,nelem,&
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
      call getnumatomsall(ntrain,num_atoms_all)
!!
!!=========================================================
!! allocate arrays for fitting statistics output
!!=========================================================
      allocate(fitstate(ntrain))
      fitstate(:)=0
      allocate(fitstatf(3,max_num_atoms,ntrain))
      fitstatf(:,:,:)=0
      allocate(fitstatq(max_num_atoms,ntrain))
      fitstatq(:,:)=0
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
      maxcutoff_short_atomic = 0.0d0
      do i2=1,nelem
        do i1=1,num_funcvalues_short_atomic(i2)
          maxcutoff_short_atomic=max(maxcutoff_short_atomic,funccutoff_short_atomic(i1,i2))
        enddo ! i1
      enddo
!!
      if(mpirank.eq.0)then
!!=========================================================
!! get scaling data for the short range symmetry functions
!!=========================================================
        if(luseoldscaling)then
          call readscale(nelem,1,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
            dummy,dummy,rdummy,rdummy)
        else ! luseoldscaling
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit) !'
          call getscale(nelem,max_num_atoms,0,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            ntrain,symfunction_short_atomic_list,minvalue_short_atomic,&
            maxvalue_short_atomic,avvalue_short_atomic)
          close(symunit)
        endif ! luseoldscaling
!!
!!=========================================================
!! analyze structures in trainstruct.data and write chemical composition to runner.out
!!=========================================================
        if(lshort.and.(nelem.le.5).and.(mpirank.eq.0).and.lanalyzecomposition)then
          call analyzeinput(0) 
          call analyzeinput(1) 
        elseif(lshort.and.(nelem.gt.5))then
          write(ounit,*)'WARNING: detailed analysis of structures not possible for more than 5 elements'
        endif !'
!!
!!=========================================================
!! analyze atomic environments in trainstruct.data and teststruct.data
!!=========================================================
        if(lenvironmentanalysis.and.(mpirank.eq.0))then
          if(nelem.le.4)then
            call environmentanalysis(0,maxcutoff_short_atomic) 
            call environmentanalysis(1,maxcutoff_short_atomic) 
          elseif(lshort.and.(nelem.gt.4))then
            write(ounit,*)'WARNING: detailed analysis of atomic environments not possible for more than 4 elements'
          endif !'
        endif ! lenvironmentanalysis
!!
!!=========================================================
!! check for contradictory training data from atomic short range symmetry functions 
!!=========================================================
        if(lfindcontradictions.and.(mpirank.eq.0))then
          call findcontradictions(ntrain,maxnum_funcvalues_short_atomic,&
            num_funcvalues_short_atomic,&
            minvalue_short_atomic,&
            maxvalue_short_atomic,avvalue_short_atomic,&
            scmin_short_atomic,scmax_short_atomic)
        endif
!!
!!=========================================================
!! get statistics of energies in training set (just for process 0)
!!=========================================================
        call getenergystatistics(1,belowmaxenergy,eshortav,eshortstddev,eshortmin,eshortmax)
!!
!!=========================================================
!! get statistics of forces in training set (just for process 0)
!!=========================================================
        if((luseforces.and.lshort).or.lfinalforce)then
          call getforcestatistics(belowfmax,fmin,fmax,fvecmin,fvecmax)
        endif
!!
!!=========================================================
!! get the charge statistics, works for all nn_type_elec cases
!!=========================================================
        if(lelec)then
          avcharge(:)=0.0d0
          if(mpirank.eq.0)then
            call getavcharge(ntrain,avcharge,stddevcharge,chargemin,chargemax)
          endif ! mpirank.eq.0
        else ! set dummies for writescale
          chargemin(:)=0.0d0
          chargemax(:)=0.0d0
        endif ! lelec
!!
!!=========================================================
!! write short range scaling factors to file scaling.data
!!=========================================================
        if(.not.luseoldscaling)then
          call writescale(nelem,1,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
            eshortmin,eshortmax,chargemin,chargemax)
        endif ! luseoldscaling
      endif ! mpirank.eq.0
!!
!!=========================================================
!! distribute scaling.data to all processes
!!=========================================================
      call mpi_bcast(minvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(maxvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(avvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(eshortstddev,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(avcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(stddevcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(chargemin,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(chargemax,nelem,mpi_real8,0,mpi_comm_world,mpierror)
!!
!!=========================================================
!! write information about training data to output 
!!=========================================================
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
!!=========================================================
!! preparation of initial weight parameters'
!!=========================================================
!!
!!=========================================================
!! first initialize all weights randomly no matter if we restart a fit or not
!!=========================================================
      if(mpirank.eq.0)then
        call initialweights(nelem,&
          maxnum_weights_short_atomic,num_weights_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic,&
          windex_short_atomic,nodes_short_atomic,&
          jseed,nseed,weights_short_atomic)
!!
!!=========================================================
!! if weights according to Nguyen Widrow are requested for short range NN, overwrite random weights:
!!=========================================================
        if(lnwweights)then
          call nguyenwidrowweights(nelem,maxnum_layers_short_atomic,&
            nodes_short_atomic,jseed,windex_short_atomic,&
            maxnodes_short_atomic,maxnum_weights_short_atomic,num_weights_short_atomic,&
            weights_short_atomic,weights_min,weights_max,actfunc_short_atomic)
        endif
!!
!!=========================================================
!! If systematic weights are requested, overwrite random weights for short range NN
!!=========================================================
        if(lsysweights)then
          call systematicweights(nelem,maxnum_layers_short_atomic,&
            num_layers_short_atomic,nodes_short_atomic,&
            maxnum_weights_short_atomic,&
            weights_short_atomic,weights_min,weights_max)
        endif
!!
!!=========================================================
!! if requested overwrite short range weights with weights from weights.XXX.data
!!=========================================================
        if(luseoldweightsshort)then
          call readweights(0,nelem,&
            maxnum_weights_short_atomic,num_weights_short_atomic,&
            weights_short_atomic)
        endif
      endif ! mpirank.eq.0
!!
!!=========================================================
!! distribute weights_short to all processes
!!=========================================================
      if(lshort)then
        call mpi_bcast(weights_short_atomic,nelem*maxnum_weights_short_atomic,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!=========================================================
!! Caution: precondition is parallel inside
!! Get a better set of initial weights by preconditioning
!!=========================================================
      if(lprecond)then
        if(lelec.and.(nn_type_elec.eq.2))then
!! FIXME: works only without electrostatics, charges are not yet taken into account to find optimum initial weights
          write(ounit,*)'ERROR: precondition_short_atomic does not yet work for nn_type_elec 2'
          stop !'
        endif
        call precondition_atomic(ntrain,trainelem,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          eshortmin,eshortmax,eshortav,eshortstddev,avcharge,stddevcharge)
      endif
!!
!!=========================================================
!! write initial weights to file
!!=========================================================
      if(mpirank.eq.0)then
        call writeweights(0,nelem,maxnum_weights_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic,&
          nodes_short_atomic,weights_short_atomic)
      endif ! mpirank.eq.0
!!
!!=========================================================
!! keep copy of weights for the calculation of the weight change wshift
!!=========================================================
      weights_short_veryold(:,:)=0.0d0
      weights_short_old(:,:)    =weights_short_atomic(:,:)
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
      rmse_short            =0.0d0
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_etot             =0.0d0
      rmse_elec             =0.0d0
      rmse_short_test       =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_etot_test        =0.0d0
      rmse_elec_test        =0.0d0
      rmse_force_s          =0.0d0
      rmse_force_s_test     =0.0d0
      rmse_force_e          =0.0d0
      rmse_force_e_test     =0.0d0
      rmse_force_t          =0.0d0
      rmse_force_t_test     =0.0d0
!!
!!=========================================================
!! initialize MADs
!!=========================================================
      mad_short             =0.0d0
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_etot              =0.0d0
      mad_ewald             =0.0d0
      mad_short_test        =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_etot_test         =0.0d0
      mad_ewald_test        =0.0d0
      mad_force_s           =0.0d0
      mad_force_s_test      =0.0d0
      mad_force_e           =0.0d0
      mad_force_e_test      =0.0d0
      mad_force_t           =0.0d0
      mad_force_t_test      =0.0d0
!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!! (rmse_short,rmse_charge,rmse_totalcharge,rmse_elec,rmse_etot,rmse_force_s,rmse_force_e)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call abstime(timegeterrorstart,daygeterror)
      call geterror_atomic(0,countepoch,ntrain,&
        imaxerror_eshort_train,imaxerror_elec_train,imaxerror_etot_train,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_elec,rmse_etot,&
        rmse_force_s,rmse_force_t,rmse_force_e,&
        mad_short,mad_charge,mad_totalcharge,&
        mad_ewald,mad_etot,&
        mad_force_s,mad_force_t,mad_force_e,&
        maxerror_eshort_train,maxerror_elec_train,maxerror_etot_train)
!!
!! calculate sensitivity for initial weights if requested
      if(lsens)then
!! FIXME: this subroutine has not been adjusted:
        call getsensitivity_atomic(nelem,max_num_atoms,ntrain,trainelem,&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          scmin_short_atomic,scmax_short_atomic,sens)
      endif
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
      call geterror_atomic(1,countepoch,ntest,&
        imaxerror_eshort_test,imaxerror_elec_test,imaxerror_etot_test,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        rmse_short_test,rmse_charge_test,rmse_totalcharge_test,&
        rmse_elec_test,rmse_etot_test,&
        rmse_force_s_test,rmse_force_t_test,rmse_force_e_test,&
        mad_short_test,mad_charge_test,mad_totalcharge_test,&
        mad_ewald_test,mad_etot_test,&
        mad_force_s_test,mad_force_t_test,mad_force_e_test,&
        maxerror_eshort_test,maxerror_elec_test,maxerror_etot_test)
      call abstime(timegeterrorend,daygeterror)
      timegeterror=timegeterrorend-timegeterrorstart
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of first test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! set first references for determination of optweights files
      rmse_short_test_old  = rmse_short_test
      rmse_charge_test_old = rmse_charge_test
!!
      call abstime(timeepochend,dayepoch)
      timeepoch=timeepochend-timeepochstart
!!
!! write RMSE headers
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------'
        if(lshort)then
          write(ounit,'(a,f14.3,a)')' ### WARNING ### just short range energies below ',&
            maxenergy,' Ha/atom are used for fitting and Eshort RMSE!'
          write(ounit,'(a,f14.3,a,f14.3,a)')' => Fitted energy range has width of ',&
            maxenergy-eshortmin,' Ha/atom = ',(maxenergy-eshortmin)*27.211d0,' eV/atom'
          write(ounit,'(a,i10)')' => Number of short range training energies below max_energy: ',belowmaxenergy 
        endif
        if((luseforces.and.lshort).or.lfinalforce)then
          write(ounit,'(a,f14.3,a)')' ### WARNING ### just force components below     ',&
            maxforce,' Ha/Bohr are used for fitting and Fshort RMSE!'
          do i1=1,nelem
            write(ounit,'(x,a2,x,a,f14.3,a,f14.3,a)')element(i1),&
              ' => Fitted force range has width of ',maxforce-fmin(i1),&
              ' Ha/Bohr = ',(maxforce-fmin(i1))*27.211d0,' eV/Bohr'
            write(ounit,'(x,a2,x,a,i10)')element(i1),&
              ' => Number of short range training force components below max_force: ',belowfmax(i1) 
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
          rmse_elec             =0.0d0
          rmse_elec_test        =0.0d0
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
!! then lshort is set false to fit the charges, then rmse_etot and rmse_elec
!! cannot be the same because they have different DFT reference energies
!! (for etot it is the total DFT energy and for the electrostatic energy
!! it is only the DFT ewald energy, from the NN we in both cases have just the
!! ewald energy)
!!
!! write RMSEs'
        if(lprintdateandtime)then
          call printdateandtime(countepoch)
        endif
        write(ounit,'(a8,i5,x,10f13.6,f8.2)') ' ENERGY ',countepoch,&
            rmse_short*tounit,rmse_short_test*tounit,&
            rmse_elec*tounit,rmse_elec_test*tounit,&
            rmse_etot*tounit,rmse_etot_test*tounit,&
            rmse_charge,rmse_charge_test,&
            rmse_totalcharge,rmse_totalcharge_test,&
            timeepoch/60.d0
        if(luseforces)then
          write(ounit,'(a8,i5,x,6f13.6)')' FORCES ',countepoch,&
            rmse_force_s*tounit,rmse_force_s_test*tounit,&
            rmse_force_e*tounit,rmse_force_e_test*tounit,& 
            rmse_force_t*tounit,rmse_force_t_test*tounit 
        endif
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)') ' MADE   ',countepoch,&
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
        if(lprintconv)then
          do i1=1,nelem
            write(ounit,'(a9,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',element(i1),&
              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshift(i1),wshift2(i1)
            convvec(i1,:,3)=convvec(i1,:,2)
            convvec(i1,:,2)=convvec(i1,:,1)
            convvec(i1,:,1)=0.0d0
          enddo ! i1
        endif
        if(lweightanalysis)then
          call analyzeweights(1,countepoch,nelem,&
            maxnum_weights_short_atomic,num_weights_short_atomic,weights_short_atomic)
        endif ! lweightanalysis
!! print sensitivity if requested
!! FIXME: take sensitivity with respect to charges into account
        if(lsens)then
          write(ounit,*)'Short range NN sensitivity: '
          do i1=1,nelem
            do i2=1,num_funcvalues_short_atomic(i1)
              write(ounit,'(i5,a15,x,a2,x,i5,f16.8)')countepoch,' NNsensitivity ',&
                element(i1),i2,sens(i1,i2)  
            enddo ! i2
            write(ounit,*)'-------------------------------------------------------------'
            enddo ! i1
        endif ! lsens
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
      if(lresetkalman)then
        if(mpisize.gt.1)then
          write(ounit,*)'Error: bcast for reinitialized corrmatrix still missing'
          stop !'
        endif
        do i1=1,nelem
          if(lshort.and.(optmodee.eq.1))then
            call initialcorrmatrix(num_weightsshortfree(i1),corrmatrix_list(1,i1))
          endif
          if(lshort.and.(optmodef.eq.1))then
            call initialcorrmatrix(num_weightsshortfree(i1),corrmatrixf_list(1,i1))
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
      if(mpirank.eq.0)then
        if(lshort)then
          call readfunctions_mixed(npoints,nelem,max_num_atoms,&
            ntrain,block_counter,pointindex,0,&
            maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
            symfunction_short_atomic_list)
        endif ! lshort
!!
!! read the structures needed for the calculation of the electrostatic energy
!! and get reference forces from DFT
!! is needed for lshort (force fitting) and lelec (structure for electrostatics)
!! must be called after readfunctions because it needs num_atoms_list
        call getstructures_mixed(npoints,&
          ntrain,block_counter,pointindex,num_atoms_all)
!!
!! read short range forces from trainforces.data
        if((luseforces.and.lshort).or.(lfinalforce.and.(nepochs.eq.countepoch)))then
          call readforces_mixed(npoints,&
            ntrain,block_counter,pointindex,num_atoms_all)
        endif
!!
!! read electrostatic forces from trainforcese.data
!! we don't need this because we don't fit electrostatic forces
!!        if((luseforces.or.(lfinalforce.and.(nepochs.eq.countepoch))).and.lelec)then
!!          call readforces(trainfeunit,npoints,&
!!            elecforce_list)
!!        endif
      endif ! mpirank.eq.0
      call abstime(timeioend,dayio)
      timeio=timeio+timeioend-timeiostart
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
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(shortforce_list,nblock*max_num_atoms*3,&
        mpi_real8,0,mpi_comm_world,mpierror)
!!      call mpi_bcast(elecforce_list,nblock*max_num_atoms*3,&
!!        mpi_real8,0,mpi_comm_world,mpierror)
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
!! start optimization 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      call abstime(timeshortfitstart,dayshortfit)
!!
!! scale symmetry functions for the short-range interaction
      call scalesymfit_para(1,nelem,&
        nstruct,n_start,n_end,&
        maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        symfunction_short_atomic_list,&
        scmin_short_atomic,scmax_short_atomic)
!!
!! update using energies and forces together
      call optimize_atomic(npoints,point,idx,&
        kaldim,maxcorrdim,maxcorrfdim,corrdim,corrfdim,countepoch,&
        ndone,num_weightsshortfree,&
        wconstraintidx,numbere,numberf,numberq,&
        lseed,ntrain,&
        fitstate,fitstatf,fitstatq,kseed,&
        kalmanthreshold_temp,kalmanthresholdf_temp,&
        rmse_short_ref,rmse_force_s_ref,rmse_charge_ref,&
        corrmatrix_list,corrmatrixf_list,&
        minvalue_short_atomic,maxvalue_short_atomic)
!!
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
            maxnum_weights_short_atomic,num_weights_short_atomic,&
            weights_short_atomic)
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
            num_weightsshortfree,&
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
            maxnum_weights_short_atomic,&
            maxnum_layers_short_atomic,num_layers_short_atomic,&
            nodes_short_atomic,weights_short_atomic)
        endif
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! initialize RMSEs for next epoch
      rmse_short            =0.0d0
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_etot             =0.0d0
      rmse_elec             =0.0d0
      rmse_short_test       =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_etot_test        =0.0d0
      rmse_elec_test        =0.0d0
      rmse_force_s          =0.0d0
      rmse_force_s_test     =0.0d0
      rmse_force_e          =0.0d0
      rmse_force_e_test     =0.0d0
      rmse_force_t          =0.0d0
      rmse_force_t_test     =0.0d0
!!
!! initialize MADs for next epoch
      mad_short             =0.0d0
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_etot              =0.0d0
      mad_ewald             =0.0d0
      mad_short_test        =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_etot_test         =0.0d0
      mad_ewald_test        =0.0d0
      mad_force_s           =0.0d0
      mad_force_s_test      =0.0d0
      mad_force_e           =0.0d0
      mad_force_e_test      =0.0d0
      mad_force_t           =0.0d0
      mad_force_t_test      =0.0d0
!!
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
      block_counter          = 1
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
      timegeterrorstart=0.0d0
      daygeterror      =0
      call abstime(timegeterrorstart,daygeterror)
!!
      call geterror_atomic(0,countepoch,ntrain,&
        imaxerror_eshort_train,imaxerror_elec_train,imaxerror_etot_train,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_elec,rmse_etot,&
        rmse_force_s,rmse_force_t,rmse_force_e,&
        mad_short,mad_charge,mad_totalcharge,&
        mad_ewald,mad_etot,&
        mad_force_s,mad_force_t,mad_force_e,&
        maxerror_eshort_train,maxerror_elec_train,maxerror_etot_train)
!!
!! calculate sensitivity for this epoch if requested
      if(lsens)then
!! FIXME: this subroutine has not been adjusted:
        call getsensitivity_atomic(nelem,max_num_atoms,ntrain,trainelem,&
          maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
          minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
          scmin_short_atomic,scmax_short_atomic,sens)
      endif
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
      call geterror_atomic(1,countepoch,ntest,&
        imaxerror_eshort_test,imaxerror_elec_test,imaxerror_etot_test,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        rmse_short_test,rmse_charge_test,rmse_totalcharge_test,&
        rmse_elec_test,rmse_etot_test,&
        rmse_force_s_test,rmse_force_t_test,rmse_force_e_test,&
        mad_short_test,mad_charge_test,mad_totalcharge_test,&
        mad_ewald_test,mad_etot_test,&
        mad_force_s_test,mad_force_t_test,mad_force_e_test,&
        maxerror_eshort_test,maxerror_elec_test,maxerror_etot_test)
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
        call getwshift(nelem,maxnum_weights_short_atomic,num_weights_short_atomic,&
          weights_short_atomic,weights_short_old,weights_short_veryold,&
          wshift,wshift2)
      endif ! lshort
!!
!! keep weights
      weights_short_veryold(:,:)=weights_short_old(:,:)
      weights_short_old(:,:)    =weights_short_atomic(:,:)
!!
!! If update by element is used: Set some meaningless RMSEs to zero for output
      if(lelec.and.lupdatebyelement)then
        rmse_totalcharge       =0.0d0
        rmse_totalcharge_test  =0.0d0
        rmse_totalcharge_ref   =0.0d0
        rmse_elec              =0.0d0
        rmse_elec_test         =0.0d0
        rmse_force_e           =0.0d0
        rmse_force_e_test      =0.0d0
        rmse_force_t           =0.0d0
        rmse_force_t_test      =0.0d0
        mad_totalcharge        =0.0d0
        mad_totalcharge_test   =0.0d0
        mad_ewald              =0.0d0
        mad_ewald_test         =0.0d0
        mad_force_e            =0.0d0
        mad_force_e_test       =0.0d0
        mad_force_t            =0.0d0
        mad_force_t_test       =0.0d0
      endif
!!
!! write RMSEs
      if(mpirank.eq.0)then
        if(lprintdateandtime)then
          call printdateandtime(countepoch)
        endif
        write(ounit,'(a8,i5,x,10f13.6,f8.2)')' ENERGY ', countepoch,&
          rmse_short*tounit,rmse_short_test*tounit,&
          rmse_elec*tounit,rmse_elec_test*tounit,&
          rmse_etot*tounit,rmse_etot_test*tounit,&
          rmse_charge,rmse_charge_test,&
          rmse_totalcharge,rmse_totalcharge_test,&
          timeepoch/60.d0
        if(luseforces)then
          write(ounit,'(a8,i5,x,6f13.6)')' FORCES ',countepoch,&
          rmse_force_s*tounit,rmse_force_s_test*tounit,&
          rmse_force_e*tounit,rmse_force_e_test*tounit,& 
          rmse_force_t*tounit,rmse_force_t_test*tounit 
        endif
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)')' MADE   ', countepoch,&
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
            maxnum_weights_short_atomic,num_weights_short_atomic,weights_short_atomic)
        endif ! lweightanalysis
!! print detailed timing for epoch if requested
        if(lfinetimeepoch)then
          call writeepochtime(countepoch)
        endif ! lfinetimeepoch
!! print sensitivity if requested
        if(lsens)then
          write(ounit,*)'Short range NN sensitivity: '
          do i1=1,nelem
            do i2=1,num_funcvalues_short_atomic(i1)
              write(ounit,'(i5,a15,x,a2,x,i5,f16.8)')countepoch,' NNsensitivity ',&
                element(i1),i2,sens(i1,i2) 
            enddo ! i2
            write(ounit,*)'-------------------------------------------------------------'
          enddo ! i1
        endif ! lsens
      endif ! mpirank.eq.0
!!
      call adjustkalman_short(numbere,numberf,&
        kalmanthreshold_temp,kalmanthresholdf_temp)
!!'
!! write optimum set of weights to files optweights.out and opweightse.out
!! FIXME: WE NEED BETTER CRITERION TO DEFINE OPTIMUM EPOCH NOW
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
        optrmse_elec            =rmse_elec
        optrmse_elec_test       =rmse_elec_test
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
!! summarize optimum fit
      if((mpirank.eq.0).and.(nepochs.gt.0))then
!! FIXME: This routine has not been adjusted
        call writeoptfit_atomic(optepoch,&
          tounit,optrmse_short,optrmse_short_test,optrmse_force_s,&
          optrmse_force_s_test,optrmse_elec,optrmse_elec_test,&
          optrmse_charge,optrmse_charge_test,optrmse_totalcharge,&
          optrmse_totalcharge_test,optrmse_force_e,optrmse_force_e_test,&
          optrmse_etot,optrmse_etot_test)
      endif ! mpirank.eq.0
!!
!! report final statistics if requested:
      if(lfitstats.and.(mpirank.eq.0))then
        call writefitstat(ntrain,fitstate,fitstatf,fitstatq)
      endif
!!
!! analyze error distribution if requested
      if(lanalyzeerror.and.(mpirank.eq.0))then
        call erroranalysis(ntrain,ntest)
      endif
!!
      if(allocated(wconstraintidx))  deallocate(wconstraintidx)  
!! deallocate Kalman matrices
      if(allocated(corrmatrix_list)) deallocate(corrmatrix_list)
      if(allocated(corrmatrixf_list))deallocate(corrmatrixf_list)
      if(allocated(fitstate))        deallocate(fitstate)
      if(allocated(fitstatf))        deallocate(fitstatf)
      if(allocated(fitstatq))        deallocate(fitstatq)
      if(allocated(pointindex))      deallocate(pointindex)
      if(allocated(num_atoms_all))   deallocate(num_atoms_all)
      deallocate(filenamews)
!!
!!    deallocate arrays of structures module
      call deallocatestructures()
!!
      if(lsens)then
        deallocate(sens)
      endif
!!
      return
      end
