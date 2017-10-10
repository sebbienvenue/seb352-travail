!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine fitting_electrostatic(iseed)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use symfunctions
      use nnewald
      use structures
      use timings
!!
      implicit none
!!
!! for dimensions
      integer num_weightsewaldfree(nelem)                       ! internal
      integer num_weightsewaldfixed(nelem)                      ! internal
!! counters
      integer npoints                                           ! internal
      integer point                                             ! internal
      integer pointe                                            ! internal
      integer ncount                                            ! internal
      integer ndone                                             ! internal
      integer countepoch                                        ! internal
      integer i1,i2                                             ! internal 
      integer optepoche                                         ! internal
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
      integer numberq                                           ! internal
      integer ntrainatoms                                       ! internal
      integer ntestatoms                                        ! internal
      integer trainelem(nelem)                                  ! internal
      integer testelem(nelem)                                   ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidxe ! internal 
      integer, dimension(:,:)  , allocatable :: fitstatq        ! internal 
      integer, dimension(:)  , allocatable :: pointindex        ! internal 
      integer, dimension(:)  , allocatable :: num_atoms_all     ! internal 
      integer isum                                              ! internal
!! Kalman matrix dimensions:
      integer corredim(nelem)                                   ! internal
      integer corrcdim
      integer kaledim(nelem)
      integer kalcdim
      integer maxkaledim
      integer maxcorredim                                       ! internal
      integer block_counter                                     ! internal
!!
!! symmetry function related arrays
      real*8 minvalue_elec(nelem,maxnum_funcvalues_elec)
      real*8 maxvalue_elec(nelem,maxnum_funcvalues_elec)
      real*8 avvalue_elec(nelem,maxnum_funcvalues_elec)
!! weights
      real*8 optweights_ewald(maxnum_weights_elec,nelem)
      real*8 weights_ewald_veryold(maxnum_weights_elec,nelem)        ! internal for lprintconv
      real*8 weights_ewald_old(maxnum_weights_elec,nelem)
!! Reference data:
      real*8 chargemin(nelem)                                        ! internal
      real*8 chargemax(nelem)                                        ! internal
      real*8 dummy                                                   ! internal
!! RMSEs
      real*8 rmse_charge
      real*8 rmse_totalcharge
      real*8 rmse_elec
      real*8 rmse_charge_test
      real*8 rmse_charge_test_old
      real*8 rmse_totalcharge_test
      real*8 rmse_elec_test
      real*8 rmse_charge_ref
      real*8 rmse_totalcharge_ref
!! MADs
      real*8 mad_elec
      real*8 mad_elec_test
      real*8 mad_charge
      real*8 mad_charge_test
      real*8 mad_totalcharge
      real*8 mad_totalcharge_test
!! for final output
      real*8 optrmse_charge
      real*8 optrmse_totalcharge
      real*8 optrmse_elec
      real*8 optrmse_charge_test
      real*8 optrmse_totalcharge_test
      real*8 optrmse_elec_test

!! Kalman filter 
      real*8 kalmanthresholde_temp                            ! internal
      real*8, dimension(:,:), allocatable :: corrmatrixe_list
      real*8, dimension(:)  , allocatable :: corrmatrixc 
!! arrays for convergence vector
      real*8 wshifte(nelem)                                   ! internal
      real*8 wshifteold(nelem)                                ! internal
      real*8 wshifte2(nelem)                                  ! internal
      real*8 convvec(nelem,2,3)                               ! internal
      real*8 convtemp1                                        ! internal
      real*8 convtemp2                                        ! internal
      real*8 convtemp3                                        ! internal

!! miscellaneous
      real*8 tounit                                           ! internal
      real*8 avcharge(nelem)                                  ! internal
      real*8 stddevcharge(nelem)                              ! internal
      real*8, dimension(:,:)  , allocatable :: sense          ! internal
!!
!!
!!=============================================================
!! allocate sensitivity array if requested 
!!=============================================================
      if(lsens)then
        allocate(sense(nelem,maxnum_funcvalues_elec))
        sense(:,:)=0.0d0
      endif
!!
!!=============================================================
!!    allocate arrays of structures module
!!=============================================================
      call allocatestructures()
!!
!!=============================================================
!! initialization of local Kalman filter parameters
!!=============================================================
      kalmanthresholde_temp  = kalmanthresholde
!!=============================================================
!! initialization of scaling parameters
!!=============================================================
      minvalue_elec(:,:)         = 0.0d0
      maxvalue_elec(:,:)         = 0.0d0
      avvalue_elec(:,:)          = 0.0d0
!!
!!=============================================================
!! set filenames for weights files
!!=============================================================
      allocate(filenamewe(nelem))
      call getweightfilenames_elec()
!!
!!=============================================================
!! initialization of reference data arrays
!!=============================================================
      totalcharge_list(:)    = 0.0d0
      totalenergy_list(:)    = 0.0d0
      shortenergy_list(:)    = 0.0d0
      elecenergy_list(:)     = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
      shortforce_list(:,:,:) = 0.0d0
      lattice_list(:,:,:)    = 0.0d0
!!
!!=============================================================
!! setting the energy and force unit converter
!!=============================================================
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV 
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0   !  stay with Ha 
      endif
!!
!!=============================================================
!! initialization of reference RMSEs
!!=============================================================
      rmse_charge_ref     =0.0d0
      rmse_totalcharge_ref=0.0d0
!!
!!=============================================================
!! initialization of optrmses
!!=============================================================
      optrmse_charge             =0.0d0
      optrmse_charge_test        =0.0d0
      optrmse_elec               =0.0d0
      optrmse_elec_test          =0.0d0
      optrmse_totalcharge        =0.0d0
      optrmse_totalcharge_test   =0.0d0
!!
!!=============================================================
!! generating local copies of seed for various purposes
!! iseed is just used for random order training below
!!=============================================================
      jseed                  = iseed ! for initial short range weights 
      kseed                  = iseed ! for forcernd 
      lseed                  = iseed ! for energyrnd
      mseed                  = iseed ! for chargernd 
      nseed                  = iseed ! for initial electrostatic weights
      oseed                  = iseed ! for pointindex 
!!=============================================================
!! initializations of counters
!!=============================================================
      countepoch             = 0
      point                  = 0
      pointe                 = 0
      npoints                = 0
      numberq                = 0
      optepoche              = 0
      block_counter          = 1
!!=============================================================
!! initialization of timing variables
!!=============================================================
      timegeterror           = 0.0d0
      timeqerror             = 0.0d0
      timedqdw               = 0.0d0
      timequpdate            = 0.0d0
      timeio                 = 0.0d0
      timemix                = 0.0d0
      timeelecfit            = 0.0d0
      daygeterror            = 0
      dayio                  = 0
      daymix                 = 0
      dayelecfit             = 0
!!=============================================================
!! further initializations
!!=============================================================
      chargemin(:)           = 100000.d0
      chargemax(:)           = -100000.d0
      convvec(:,:,:)         = 0.0d0
      wshifte(:)             = 0.0d0
      wshifte2(:)            = 0.0d0
!!
!! determine starting time
      call abstime(timeepochinistart,dayepochini)
!!
!!=============================================================
!!=============================================================
!! end of initializations
!!=============================================================
!!=============================================================
!!
!!=============================================================
!! get information on fixed weights from input.nn file
!!=============================================================
      allocate(wconstraintidxe(maxnum_weights_elec,nelem))  
      call getwconstraintidx(1,nelem,windex_elec,&
        maxnum_layers_elec,num_layers_elec,maxnum_weights_elec,&
        num_weights_elec,num_weightsewaldfree,num_weightsewaldfixed,&
        nodes_elec,wconstraintidxe)
      call mpi_bcast(num_weightsewaldfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsewaldfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!=============================================================
!! get dimensions for the Kalman filter
!!=============================================================
      if(optmodeq.eq.1)then
        call getkaldims_elec(&
          kaledim,kalcdim,corredim,corrcdim,&
          maxkaledim,maxcorredim,&
          num_weightsewaldfree)
      endif
!!
!!=============================================================
!! allocate Kalman matrices
!!=============================================================
      allocate(corrmatrixe_list(maxcorredim,nelem))
      allocate(corrmatrixc(corrcdim)) 
!!
!!=============================================================
!! get the values for the Kalman matrices and parameters (new or from file)
!!=============================================================
      call getkalmanmatrices_elec(&
        iseed,kseed,lseed,mseed,&
        corredim,corrcdim,&
        maxcorredim,&
        num_weightsewaldfree,&
        corrmatrixe_list,corrmatrixc)
      call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(kseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mseed,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!=============================================================
!! count the number of training and test points
!!=============================================================
      if(mpirank.eq.0)then
        call countpoints(1,nelem,&
          ntrain,ntest,ntrainatoms,ntestatoms,trainelem,testelem)
      endif ! mpirank.eq.0
!!
      call mpi_bcast(ntrain,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntest,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntrainatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntestatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(trainelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(testelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
      allocate(pointindex(ntrain))
      allocate(num_atoms_all(ntrain))
!!
!!=============================================================
!! get array num_atoms_all, contains number of atoms in each training point
!!=============================================================
      call getnumatomsall(ntrain,num_atoms_all)
!!
!!=============================================================
!! allocate arrays for fitting statistics output
!!=============================================================
      allocate(fitstatq(max_num_atoms,ntrain))
      fitstatq(:,:)=0
!!
!!=============================================================
!! if growth mode is not used make sure the full training set is used
!!=============================================================
      if(.not.lgrowth)then
        ngrowth=ntrain
      endif
!!
!!=============================================================
!! get maxcutoff_elec of electrostatic symmetry functions
!!=============================================================
      maxcutoff_elec = 0.0d0
      do i2=1,nelem
        do i1=1,num_funcvalues_elec(i2)
          maxcutoff_elec=max(maxcutoff_elec,funccutoff_elec(i1,i2))
        enddo ! i1
      enddo ! i2
!!
      if(mpirank.eq.0)then
!!=============================================================
!! get scaling data
!!=============================================================
        if(luseoldscaling)then
          call readscale(nelem,3,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            dummy,dummy,chargemin,chargemax) ! FIXME: CHECK: chargemin and chargemax should be a dummy? 
        else ! luseoldscaling
!! get scaling data for the charge symmetry functions 
          open(symunit,file='functione.data',form='formatted',status='old')
          rewind(symunit) !'
          call getscale(nelem,max_num_atoms,1,&
            maxnum_funcvalues_elec,num_funcvalues_elec,&
            ntrain,symfunction_elec_list,minvalue_elec,maxvalue_elec,avvalue_elec)
          close(symunit)
        endif ! luseoldscaling
      endif
!!
!!=============================================================
!! analyze structures in trainstruct.data and write chemical composition to runner.out
!!=============================================================
        if((nelem.le.5).and.(mpirank.eq.0).and.lanalyzecomposition)then
          call analyzeinput(2) 
          call analyzeinput(3) 
        elseif(nelem.gt.5)then
          write(ounit,*)'WARNING: detailed analysis of structures not possible for more than 5 elements'
        endif !'
!!
!!=============================================================
!! analyze atomic environments in trainstruct.data and teststruct.data
!!=============================================================
        if(lenvironmentanalysis.and.(mpirank.eq.0))then
          if(nelem.le.4)then
            call environmentanalysis(0,maxcutoff_elec) 
            call environmentanalysis(1,maxcutoff_elec) 
          elseif(nelem.gt.4)then
            write(ounit,*)'ERROR: detailed analysis of atomic environments not possible for more than 4 elements'
            stop
          endif !'
        endif ! lenvironmentanalysis
!!
        if(lfindcontradictions.and.(mpirank.eq.0))then
          write(ounit,*)'ERROR: find_contradictions works for short range case only'
          stop !' 
        endif
!!
!!=============================================================
!! get the charge statistics
!!=============================================================
        avcharge(:)=0.0d0
        if(mpirank.eq.0)then
          call getavcharge(ntrain,avcharge,stddevcharge,chargemin,chargemax)
        endif ! mpirank.eq.0
        call mpi_bcast(avcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(stddevcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
!!
!!=============================================================
!! write electrostatic scaling factors to file scalinge.data
!!=============================================================
      if(.not.luseoldscaling)then
        call writescale(nelem,3,&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          dummy,dummy,chargemin,chargemax)
      endif ! mpirank.eq.0
!!
!!=============================================================
!! distribute scalinge.data to all processes
!!=============================================================
      call mpi_bcast(minvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(maxvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(avvalue_elec,nelem*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
!!
!!=============================================================
!! print results 
!!=============================================================
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'number of training points ',ntrain !'
        write(ounit,*)'number of training atoms  ',ntrainatoms
        write(ounit,*)'number of testing points  ',ntest
        write(ounit,*)'number of testing atoms   ',ntestatoms
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Number of atoms for each element:   '
        write(ounit,*)'           training:    testing:   '
        do i1=1,nelem
          write(ounit,'(i3,x,a3,x,2i12)')i1,element(i1),trainelem(i1),testelem(i1)
        enddo
      endif ! mpirank.eq.0 '
!!
!!=============================================================
!! preparation of initial weight parameters'
!!=============================================================
!!
!!=============================================================
!! first initialize all weights randomly no matter if we restart a fit or not
!!=============================================================
      if(mpirank.eq.0)then
        call initialweights(nelem,&
          maxnum_weights_elec,num_weights_elec,&
          maxnum_layers_elec,num_layers_elec,windex_elec,nodes_elec,&
          jseed,nseed,weights_elec)
!!
!!=============================================================
!! if weights according to Nguyen Widrow are requested for electrostatic NN, overwrite random weights:
!!=============================================================
        if(lnwweightse)then
          call nguyenwidrowweights(nelem,maxnum_layers_elec,&
            nodes_elec,nseed,windex_elec,&
            maxnodes_elec,maxnum_weights_elec,num_weights_elec,&
            weights_elec,weightse_min,weightse_max,actfunc_elec)
        endif
!!
!!=============================================================
!! If systematic weights are requested, overwrite random weights for electrostatic NN
!!=============================================================
        if(lsysweightse)then
          call systematicweights(nelem,maxnum_layers_elec,&
            num_layers_elec,nodes_elec,&
            maxnum_weights_elec,&
            weights_elec,weightse_min,weightse_max)
        endif
!!
!!=============================================================
!! if requested overwrite charge weights with weights from weightse.XXX.data
!!=============================================================
        if(luseoldweightscharge)then
          call readweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            weights_elec)
        endif
      endif ! mpirank.eq.0
!!
!!=============================================================
!! distribute weights to all processes
!!=============================================================
      call mpi_bcast(weights_elec,nelem*maxnum_weights_elec,&
        mpi_real8,0,mpi_comm_world,mpierror)
!!
      if(lprecond)then
        call precondition_electrostatic(ntrain,trainelem,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          avcharge,stddevcharge)
      endif
!!
!!=============================================================
!! write initial weights to file
!!=============================================================
      if(mpirank.eq.0)then
        call writeweights(1,nelem,maxnum_weights_elec,&
          maxnum_layers_elec,num_layers_elec,&
          nodes_elec,weights_elec)
      endif ! mpirank.eq.0
!!
!!=============================================================
!! keep copy of weights for the calculation of the weight change wshifte
!!=============================================================
      weights_ewald_veryold(:,:)=0.0d0
      weights_ewald_old(:,:)    =weights_elec(:,:)
!!
!!=============================================================
!! end of preparation of initial weight parameters'
!!=============================================================
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!=============================================================
!! determine and print initialization time
!!=============================================================
      call abstime(timeepochiniend,dayepochini)
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.2)')&
          ' initialization time (min):',(timeepochiniend-timeepochinistart)/60.d0
      endif ! mpirank.eq.0
!!'
!!=============================================================
!! terminate RuNNer here if just the initialization but no fit and no error calculation were requested
!!=============================================================
      if(linionly)then
        write(ounit,*)'Initialization done - terminating RuNNer as requested'
        write(ounit,*)'-------------------------------------------------------------'
        return !'
      endif
!!
!!=============================================================
!! set time to zero for epoch 0 time measurement
!!=============================================================
      timeepochstart =0.0d0
      dayepoch       =0
      call abstime(timeepochstart,dayepoch)
!!
!!=============================================================
!! initialize RMSEs
!!=============================================================
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_elec             =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_elec_test        =0.0d0
!!
!!=============================================================
!! initialize MADs
!!=============================================================
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_elec             =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_elec_test        =0.0d0
!!
      maxerror_elec_train     = 0.0d0 
      maxerror_elec_test      = 0.0d0 
      imaxerror_elec_train   = 0 
      imaxerror_elec_test    = 0 
!!
!!=============================================================
!! calculate the initial training error
!!=============================================================
      call abstime(timegeterrorstart,daygeterror)
      call geterror_elec(0,countepoch,ntrain,&
        imaxerror_elec_train,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_charge,rmse_totalcharge,&
        rmse_elec,mad_charge,mad_totalcharge,&
        mad_elec,maxerror_elec_train)
!!
!!=============================================================
!! calculate sensitivity for initial weights if requested
!!=============================================================
      if(lsens)then
        call getsensitivity_elec(ntrain,trainelem,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          sense)
      endif
!!
!!=============================================================
!! get new references for adaptive Kalman filter
!!=============================================================
      rmse_charge_ref      = rmse_charge
      rmse_totalcharge_ref = rmse_totalcharge
!!
!!=============================================================
!! calculate the first test error
!!=============================================================
      call geterror_elec(1,countepoch,ntest,&
        imaxerror_elec_test,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_charge_test,rmse_totalcharge_test,&
        rmse_elec_test,mad_charge_test,mad_totalcharge_test,&
        mad_elec_test,maxerror_elec_test)
      call abstime(timegeterrorend,daygeterror)
      timegeterror=timegeterrorend-timegeterrorstart
!!
!!=============================================================
!! set first references for determination of optweightse files
!!=============================================================
      rmse_charge_test_old = rmse_charge_test
!!
!!=============================================================
!! timing for epoch 0 
!!=============================================================
      call abstime(timeepochend,dayepoch)
      timeepoch=timeepochend-timeepochstart
!!
!!=============================================================
!! write RMSE headers
!!=============================================================
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------'
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
          mad_totalcharge       =0.0d0
          mad_totalcharge_test  =0.0d0
          mad_elec             =0.0d0
          mad_elec_test        =0.0d0
        endif
        if(fitting_unit.eq.1)then
          write(ounit,*)'RMSEs (energies: eV/atom, charges: e, forces: eV/Bohr):'
        elseif(fitting_unit.eq.2)then
          write(ounit,*)'RMSEs (energies: Ha/atom, charges: e, forces: Ha/Bohr):'
        else
          write(ounit,*)'Error: unknown energy unit in fitting.f90'
          stop
        endif
        write(ounit,'(6a)')'                      --- el. energies --- ',&
                         '    --- atom charges: --- ',&
                         '   --- total charges: ---',&
                         ' - time -'
        write(ounit,'(6a)')'                          /atom      ',&
                         '             /e            ',&
                         '             /e           ',&
                         '   min'
        write(ounit,'(5a)')'       epoch         train         test',&
                    '        train         test',&
                    '        train         test'
!!
!!=============================================================
!! write RMSEs'
!!=============================================================
        if(lprintdateandtime)then
          call printdateandtime(countepoch)
        endif
        write(ounit,'(a8,i5,x,6f13.6,f8.2)') ' ENERGY ',countepoch,&
          rmse_elec*tounit,rmse_elec_test*tounit,&
          rmse_charge,rmse_charge_test,&
          rmse_totalcharge,rmse_totalcharge_test,&
          timeepoch/60.d0
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)') ' MADE   ',countepoch,&
            mad_elec*tounit,mad_elec_test*tounit,&
            mad_charge,mad_charge_test,&
            mad_totalcharge,mad_totalcharge_test
        endif
        if(lprintconv)then
          do i1=1,nelem
            write(ounit,'(a9,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',element(i1),&
              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshifte(i1),wshifte2(i1)
            convvec(i1,:,3)=convvec(i1,:,2)
            convvec(i1,:,2)=convvec(i1,:,1)
            convvec(i1,:,1)=0.0d0
          enddo ! i1
        endif
!!
!!=============================================================
!! perform analysis of weight vector if requested 
!!=============================================================
        if(lweightanalysis)then
          call analyzeweights(2,countepoch,nelem,&
            maxnum_weights_elec,num_weights_elec,weights_elec)
        endif ! lweightanalysis
!!
!!=============================================================
!! print sensitivity if requested
!!=============================================================
        if(lsens)then
          write(ounit,*)'Electrostatic NN sensitivity: '
          do i1=1,nelem
            do i2=1,num_funcvalues_elec(i1)
              write(ounit,'(i5,a16,x,a2,x,i5,f16.8)')countepoch,' NNsensitivitye ',&
                element(i1),i2,sense(i1,i2)  
            enddo ! i2
            write(ounit,*)'-------------------------------------------------------------'
          enddo ! i1
        endif ! lsens
      endif ! mpirank.eq.0
!!
!!=============================================================
!!=============================================================
!! loop over all training epochs
!!=============================================================
!!=============================================================
      do countepoch=1,nepochs
        write(ounit,*)'-------------------------------------------------------------------------------'
        if(mpirank.eq.0) write(debugunit,'(a6,i6)')'epoch ',countepoch 
        call abstime(timeepochstart,dayepoch) !'
!!
!!=============================================================
!! epoch-specific initializations
!!=============================================================
        ndone                 = 0
        ncount = (int((countepoch-1)/growthstep)+1)*ngrowth    
        ncount = min(ncount,ntrain)
        point                 = 0
        pointe                = 0
        numberq               = 0
!!
!!=============================================================
!! get array pointindex for randomly mixing points of training set
!!=============================================================
        call abstime(timemixstart,daymix)
        call getpointindex(ncount,oseed,pointindex)
        call abstime(timemixend,daymix)
        timemix=timemixend-timemixstart
!!
!!=============================================================
!! loop block-wise over training points
!!=============================================================
 11   continue
      if(ncount.gt.nblock)then
        npoints=nblock
        ncount=ncount-nblock
      else
        npoints=ncount
        ncount=ncount-npoints
      endif
!!
!!=============================================================
!! reinitialize correlation matrices of Kalman filter if requested
!!=============================================================
      if(lresetkalman)then
        if(mpisize.gt.1)then
          write(ounit,*)'Error: bcast for reinitialized corrmatrix still missing'
          stop !'
        endif
        do i1=1,nelem
          if(optmodeq.eq.1)then
            call initialcorrmatrix(num_weightsewaldfree(i1),corrmatrixe_list(1,i1))
          endif
        enddo
        if(lchargeconstraint.and.(optmodeq.eq.1))then
          isum=0
          do i1=1,nelem
            isum=isum+num_weightsewaldfree(i1)
          enddo
          call initialcorrmatrix(isum,corrmatrixc)
        endif
      endif ! lresetkalman
!!
!!=============================================================
!! get index for random order of training points => idx
!!=============================================================
      call getidx(npoints,iseed,idx)
!!
!!=============================================================
!! read npoint data sets      
!! do all file reading for training data here at one place to allow for parallelization
!!=============================================================
      call abstime(timeiostart,dayio)
      if(mpirank.eq.0)then
!!=============================================================
!! read electrostatic symmetry functions 
!!=============================================================
        call readfunctions_mixed(npoints,nelem,max_num_atoms,&
          ntrain,block_counter,pointindex,1,&
          maxnum_funcvalues_elec,num_funcvalues_elec,&
          symfunction_elec_list)
!!
!!=============================================================
!! read the structures needed for the calculation of the electrostatic energy
!! must be called after readfunctions because it needs num_atoms_list
!!=============================================================
        call getstructures_mixed(npoints,&
          ntrain,block_counter,pointindex,num_atoms_all)
      endif ! mpirank.eq.0
      call abstime(timeioend,dayio)
      timeio=timeio+timeioend-timeiostart
!!
!!=============================================================
!! distribute the training data to all processes
!!=============================================================
      call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalcharge_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(totalenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(elecenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(symfunction_elec_list,nblock*max_num_atoms*maxnum_funcvalues_elec,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(atomcharge_list,nblock*max_num_atoms,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
!!
!!=============================================================
!! end of io part/file reading for training data 
!!=============================================================
!!
!!=============================================================
!! determine which structures of this block should be scaled by this process
!!=============================================================
      call mpifitdistribution(npoints,nstruct,n_start,n_end)
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!=============================================================
!! start optimization of the charges
!!=============================================================
      call abstime(timeelecfitstart,dayelecfit)
!!
!!=============================================================
!! scale electrostatic symmetry functions for a subset of structures (parallel) 
!!=============================================================
      call scalesymfit_para(3,nelem,&
        nstruct,n_start,n_end,&
        maxnum_funcvalues_elec,num_funcvalues_elec,&
        minvalue_elec,maxvalue_elec,avvalue_elec,symfunction_elec_list,&
        scmin_elec,scmax_elec)
!!
!!=============================================================
!! optimize the atomic charges point by point and atom by atom
!!=============================================================
      if(optmodeq.le.3)then
        call optimize_ewald(npoints,pointe,idx,&
          countepoch,ntrain,ndone,&
          maxcorredim,kaledim,kalcdim,corredim,corrcdim,&
          num_weightsewaldfree,fitstatq,&
          wconstraintidxe,mseed,numberq,&
          kalmanthresholde_temp,&
          rmse_charge_ref,rmse_totalcharge_ref,&
          corrmatrixe_list,corrmatrixc)
      else
        write(ounit,*)'ERROR: unknown optmodeq in fitting ',optmodeq
        stop
      endif ! optmodeq
!!
      call abstime(timeelecfitend,dayelecfit)
      timeelecfit=timeelecfit+timeelecfitend-timeelecfitstart
!!
!!=============================================================
!! end optimization of the charges
!!=============================================================
!!
!!=============================================================
!! write temporary weights if requested
!!=============================================================
      if(mpirank.eq.0)then
        if(lwritetmpweights)then
          call writetmpweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            weights_elec)
        endif
      endif ! mpirank.eq.0
!! 
!! count finished structures
      ndone        =ndone+npoints
      block_counter=block_counter+npoints
!!
!!=============================================================
!! if there are structures left go to next group of structures 
!!=============================================================
      if(ncount.gt.0) goto 11 
!!=============================================================
!! end block-wise loop over training points
!!=============================================================
!!
!!=============================================================
!! write Kalman data for restart if requested
!!=============================================================
      if(lsavekalman)then
        if(optmodeq.eq.1)then
          call writekalman_elec(&
            iseed,kseed,lseed,mseed,&
            maxcorredim,kaledim,corrcdim,&
            num_weightsewaldfree,&
            corrmatrixe_list,corrmatrixc)
        endif ! optmode
      endif !lsavekalman
!!
!!=============================================================
!! write final weights of this epoch to files
!!=============================================================
      if(mpirank.eq.0)then
        if(mod(countepoch,iwriteweight).eq.0) then
          call getfilenames(countepoch)
          call writeweights(1,nelem,&
            maxnum_weights_elec,&
            maxnum_layers_elec,num_layers_elec,&
            nodes_elec,weights_elec)
        endif
      endif ! mpirank.eq.0
!!
!!=============================================================
!! initialize RMSEs for next epoch
!!=============================================================
      rmse_charge           =0.0d0
      rmse_totalcharge      =0.0d0
      rmse_elec             =0.0d0
      rmse_charge_test      =0.0d0
      rmse_totalcharge_test =0.0d0
      rmse_elec_test        =0.0d0
!!
!!=============================================================
!! initialize MADs for next epoch
!!=============================================================
      mad_charge            =0.0d0
      mad_totalcharge       =0.0d0
      mad_elec             =0.0d0
      mad_charge_test       =0.0d0
      mad_totalcharge_test  =0.0d0
      mad_elec_test        =0.0d0
!!
      maxerror_elec_train     = 0.0d0 
      maxerror_elec_test      = 0.0d0 
      imaxerror_elec_train   = 0 
      imaxerror_elec_test    = 0 
      block_counter          = 1
!!
!!=============================================================
!! calculate this epoch's training error
!!=============================================================
      timegeterrorstart=0.0d0
      daygeterror      =0
      call abstime(timegeterrorstart,daygeterror)
!!
      call geterror_elec(0,countepoch,ntrain,&
        imaxerror_elec_train,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_charge,rmse_totalcharge,&
        rmse_elec,mad_charge,mad_totalcharge,&
        mad_elec,maxerror_elec_train)
!!
!!=============================================================
!! calculate sensitivity for this epoch if requested
!!=============================================================
      if(lsens)then
        call getsensitivity_elec(ntrain,trainelem,&
          minvalue_elec,maxvalue_elec,avvalue_elec,&
          sense)
      endif
!!
!!=============================================================
!! get new references for adaptive Kalman filter (this has nothing to do with optweights files)
!!=============================================================
      if(lelec)  rmse_charge_ref     =rmse_charge
      if(lelec)  rmse_totalcharge_ref=rmse_totalcharge
!!
!!=============================================================
!! calculate the epoch's testing error
!!=============================================================
      call geterror_elec(1,countepoch,ntest,&
        imaxerror_elec_test,&
        minvalue_elec,maxvalue_elec,avvalue_elec,&
        rmse_charge_test,rmse_totalcharge_test,&
        rmse_elec_test,mad_charge_test,mad_totalcharge_test,&
        mad_elec_test,maxerror_elec_test)
!!
      call abstime(timegeterrorend,daygeterror)
      timegeterror=timegeterrorend-timegeterrorstart
!!
!!=============================================================
!! calculate time of this epoch 
!!=============================================================
      call abstime(timeepochend,dayepoch)
      timeepoch=timeepochend-timeepochstart
!!
!!=============================================================
!! calculate change of the charge weights for INFO and CONVVEC output
!!=============================================================
      if(lelec.and.(nn_type_elec.eq.1))then
        call getwshift(nelem,maxnum_weights_elec,num_weights_elec,&
          weights_elec,weights_ewald_old,weights_ewald_veryold,&
          wshifte,wshifte2)
      endif
!!
!!=============================================================
!! keep copy of weights for the calculation of the weight change wshifte
!!=============================================================
      weights_ewald_veryold(:,:)=weights_ewald_old(:,:)
      weights_ewald_old(:,:)    =weights_elec(:,:)
!!
!!=============================================================
!! If update by element is used: Set some meaningless RMSEs and MADs to zero for output
!!=============================================================
      if(lupdatebyelement)then
        rmse_totalcharge       =0.0d0
        rmse_totalcharge_test  =0.0d0
        rmse_totalcharge_ref   =0.0d0
        rmse_elec              =0.0d0
        rmse_elec_test         =0.0d0
        mad_totalcharge        =0.0d0
        mad_totalcharge_test   =0.0d0
        mad_elec              =0.0d0
        mad_elec_test         =0.0d0
      endif
!!
!!=============================================================
!! write RMSEs
!!=============================================================
      if(mpirank.eq.0)then
        if(lprintdateandtime)then
          call printdateandtime(countepoch)
        endif
        write(ounit,'(a8,i5,x,10f13.6,f8.2)')' ENERGY ', countepoch,&
          rmse_elec*tounit,rmse_elec_test*tounit,&
          rmse_charge,rmse_charge_test,&
          rmse_totalcharge,rmse_totalcharge_test,&
          timeepoch/60.d0
        if(lprintmad)then
          write(ounit,'(a8,i5,x,10f13.6)')' MADE   ', countepoch,&
            mad_elec*tounit,mad_elec_test*tounit,&
            mad_charge,mad_charge_test,&
            mad_totalcharge,mad_totalcharge_test
        endif
        write(ounit,'(a,i5,3i10,2f20.15)')' INFORMATION USED FOR UPDATE (Q) ',&
          countepoch,numberq !'
        if(lprintconv)then
          do i1=1,nelem
            if(countepoch.eq.1)then
              convvec(i1,1,1)=wshifte(i1)
              convvec(i1,2,1)=0.0d0
            else
              convtemp1=atan((convvec(i1,2,2)-convvec(i1,2,3))/(convvec(i1,1,2)-convvec(i1,1,3)))
              convtemp2=wshifteold(i1) !abs((convvec(i1,1,3)-convvec(i1,1,2))/dcos(convtemp1))
              convtemp3=acos((-wshifte(i1)**2+convtemp2**2+wshifte2(i1)**2)/(2.d0*convtemp2*wshifte2(i1)))
              convtemp1=convtemp1+convtemp3
              convvec(i1,1,1)=convvec(i1,1,3)+wshifte2(i1)*dcos(convtemp1)
              convvec(i1,2,1)=convvec(i1,2,3)+wshifte2(i1)*dsin(convtemp1)
            endif
            write(ounit,'(a9,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',element(i1),&
              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshifte(i1),wshifte2(i1)
            convvec(i1,:,3)=convvec(i1,:,2)
            convvec(i1,:,2)=convvec(i1,:,1)
            convvec(i1,:,1)=0.0d0
            wshifteold(i1)=wshifte(i1)
          enddo ! i1
        endif ! lprintconv
!!
!!=============================================================
!! perform analysis of weight vector if requested 
!!=============================================================
        if(lweightanalysis)then
          call analyzeweights(2,countepoch,nelem,&
            maxnum_weights_elec,num_weights_elec,weights_elec)
        endif ! lweightanalysis
!!
!!=============================================================
!! print detailed timing for epoch if requested
!!=============================================================
        if(lfinetimeepoch)then
          call writeepochtime(countepoch)
        endif ! lfinetimeepoch
!!
!!=============================================================
!! print sensitivity if requested
!!=============================================================
        if(lsens)then
          write(ounit,*)'Electrostatic NN sensitivity: '
          do i1=1,nelem
            do i2=1,num_funcvalues_elec(i1)
              write(ounit,'(i5,a16,x,a2,x,i5,f16.8)')countepoch,' NNsensitivitye ',&
                element(i1),i2,sense(i1,i2)
            enddo ! i2
            write(ounit,*)'-------------------------------------------------------------'
          enddo ! i1
        endif ! lsens
      endif ! mpirank.eq.0
!!
!!=============================================================
!! adjust Kalman threshold if needed 
!!=============================================================
      call adjustkalman_elec(numberq,kalmanthresholde_temp)
!!'
!!=============================================================
!! write optimum set of weights opweightse.out
!!=============================================================
      if(rmse_charge_test.le.rmse_charge_test_old)then
        if(mpirank.eq.0)then
          optepoche=countepoch
          optweights_ewald(:,:)=weights_elec(:,:)
          call writeoptweights(1,nelem,&
            maxnum_weights_elec,num_weights_elec,&
            optweights_ewald)
        endif ! mpirank.eq.0
!!
!!=============================================================
!! define new reference only if we found a new minimum rmse
!!=============================================================
        rmse_charge_test_old    =rmse_charge_test
!!=============================================================
!! store best values obtained so far for final output
!!=============================================================
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
!!=============================================================
!!=============================================================
!! end loop over all epochs
!!=============================================================
!!=============================================================
!!
!!=============================================================
!! summarize optimum fit
!!=============================================================
      if((mpirank.eq.0).and.(nepochs.gt.0))then
        call writeoptfit_elec(optepoche,tounit,&
          optrmse_elec,optrmse_elec_test,&
          optrmse_charge,optrmse_charge_test,optrmse_totalcharge,&
          optrmse_totalcharge_test)
      endif ! mpirank.eq.0
!!
!!=============================================================
!! report final statistics if requested:
!!=============================================================
      if(lfitstats.and.(mpirank.eq.0))then
        call writefitstat_elec(ntrain,fitstatq)
      endif
!!
!!=============================================================
!! analyze error distribution if requested
!!=============================================================
      if(lanalyzeerror.and.(mpirank.eq.0))then
        call erroranalysis(ntrain,ntest)
      endif
!!
!!=============================================================
!! final cleanup 
!!=============================================================
      if(allocated(wconstraintidxe)) deallocate(wconstraintidxe)  
      if(allocated(corrmatrixe_list))deallocate(corrmatrixe_list)
      if(allocated(corrmatrixc))     deallocate(corrmatrixc) 
      if(allocated(fitstatq))        deallocate(fitstatq)
      if(allocated(pointindex))      deallocate(pointindex)
      if(allocated(num_atoms_all))   deallocate(num_atoms_all)
      deallocate(filenamewe)
!!
!!    deallocate arrays of structures module
      call deallocatestructures()
!!
      if(lsens)then
        deallocate(sense)
      endif
!!
      return
      end
