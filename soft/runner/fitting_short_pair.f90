!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine fitting_short_pair(iseed)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags 
      use globaloptions
      use symfunctions 
      use nnshort_pair
      use structures
      use timings
!!
      implicit none
!!
!! for dimensions
      integer num_weightspairfree(npairs)                       ! internal
      integer num_weightspairfixed(npairs)                      ! internal
!! seeds 
      integer iseed                                    
      integer jseed                                    
      integer kseed
      integer lseed
      integer mseed
      integer nseed
      integer oseed
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
      integer ntrain                                            ! internal
      integer ntest                                             ! internal
      integer idx(nblock)                                       ! internal
      integer nstruct                                           ! internal
      integer n_start,n_end                                     ! internal
      integer numbere                                           ! internal
      integer numberf                                           ! internal
      integer ntrainatoms                                       ! internal
      integer ntestatoms                                        ! internal
      integer trainpair(npairs)                                 ! internal
      integer testpair(npairs)                                  ! internal
      integer block_counter                                     ! internal
      integer day
!!
      integer npaircount                                        ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidxp ! internal 
      integer, dimension(:)  , allocatable :: fitstat           ! internal
      integer, dimension(:,:,:)  , allocatable :: fitstatf      ! internal
      integer, dimension(:)  , allocatable :: pointindex        ! internal
      integer, dimension(:)  , allocatable :: num_atoms_all     ! internal
!! Kalman matrix dimensions:
      integer corrdim(npairs)
      integer corrfdim(npairs)
      integer kaldim(npairs)
      integer maxkaldim
      integer maxcorrdim
      integer maxcorrfdim
!!
!! weights
      real*8 optweights_pair(maxnum_weights_short_pair,npairs)      ! internal
      real*8 weights_pair_old(maxnum_weights_short_pair,npairs)     ! internal
      real*8 weights_pair_veryold(maxnum_weights_short_pair,npairs) ! internal
!! symmetry function related arrays
      real*8 minvalue_short_pair(npairs,maxnum_funcvalues_short_pair)         ! internal
      real*8 maxvalue_short_pair(npairs,maxnum_funcvalues_short_pair)         ! internal
      real*8 avvalue_short_pair(npairs,maxnum_funcvalues_short_pair)          ! internal
!! Reference data:
      real*8 dummy                                                  ! internal
      real*8 rdummy(nelem)                                          ! internal
      real*8 fmin(nelem)
      real*8 fmax(nelem)
      real*8 fvecmin(nelem)
      real*8 fvecmax(nelem)
!! RMSEs
      real*8 rmse_short
      real*8 rmse_short_test
      real*8 rmse_short_test_old
      real*8 rmse_force_s
      real*8 rmse_force_s_test
      real*8 rmse_short_ref
      real*8 rmse_force_s_ref
!! MADs
      real*8 mad_short
      real*8 mad_short_test
      real*8 mad_force_s
      real*8 mad_force_s_test
!! for final output
      real*8 optrmse_short
      real*8 optrmse_short_test
      real*8 optrmse_force_s
      real*8 optrmse_force_s_test
!! Kalman filter
      real*8 kalmanthreshold_temp                             ! internal
      real*8 kalmanthresholdf_temp                            ! internal
      real*8, dimension(:,:), allocatable :: corrmatrix_list
      real*8, dimension(:,:), allocatable :: corrmatrixf_list
!! arrays for convergence vector
      real*8 wshiftp(npairs)                                  ! internal
      real*8 wshiftp2(npairs)                                 ! internal
      real*8 wshiftpold(npairs)                               ! internal
      real*8 convvec(npairs,2,3)                              ! internal
      real*8 convtemp1                                        ! internal
      real*8 convtemp2                                        ! internal
      real*8 convtemp3                                        ! internal
!! miscellaneous
      real*8 timestart                                        ! internal
      real*8 timeend                                          ! internal
      real*8 tounit                                           ! internal
      real*8 eshortmin                                        ! internal
      real*8 eshortmax                                        ! intenral
      real*8 eshortav                                         ! internal
      real*8 eshortstddev                                     ! internal
      real*8, dimension(:,:)  , allocatable :: sens           ! internal
!!
!!================================================
!! initializations
!!================================================
      if(lsens)then
        allocate(sens(npairs,maxnum_funcvalues_short_pair))
        sens(:,:)=0.0d0
      endif
!!
!!================================================
!! allocate arrays of structures module
!!================================================
      call allocatestructures()
!!
!!================================================
!! initialization of local Kalman filter parameters
!!================================================
      kalmanthreshold_temp   = kalmanthreshold
      kalmanthresholdf_temp  = kalmanthresholdf
!!================================================
!! initialization of scaling parameters
!!================================================
      minvalue_short_pair(:,:)   = 0.0d0
      maxvalue_short_pair(:,:)   = 0.0d0
      avvalue_short_pair(:,:)    = 0.0d0
!!================================================
!! set filenames for weights files
!!================================================
      allocate(filenamewp(npairs))
      call getweightfilenames()
!!
!!================================================
!! initialization of reference data arrays
!!================================================
      totalcharge_list(:)    = 0.0d0
      totalenergy_list(:)    = 0.0d0
      shortenergy_list(:)    = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
      shortforce_list(:,:,:) = 0.0d0
      lattice_list(:,:,:)    = 0.0d0
!!
!!================================================
!! setting the energy and force unit converter
!!================================================
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV 
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0       ! stay with Ha 
      endif
!!
!!================================================
!! initialization of reference RMSEs
!!================================================
      rmse_short_ref         =0.0d0
      rmse_force_s_ref       =0.0d0
!!
!!================================================
!! initialization of optrmses
!!================================================
      optrmse_short           =0.0d0
      optrmse_short_test      =0.0d0
      optrmse_force_s         =0.0d0
      optrmse_force_s_test    =0.0d0
!!
!!================================================
!! generating local copies of seed for various purposes
!! iseed is just used for random order training below
!!================================================
      jseed                  = iseed ! for initial short range weights  
      kseed                  = iseed ! for forcernd 
      lseed                  = iseed ! for energyrnd
      mseed                  = iseed ! for chargernd
      nseed                  = iseed ! for initial electrostatic weights
      oseed                  = iseed ! for pointindex
!!================================================
!! initializations of counters
!!================================================
      countepoch             = 0
      point                  = 0
      npoints                = 0
      numbere                = 0
      numberf                = 0
      optepoch               = 0
      block_counter          = 1
!!================================================
!! initialization of timing variables
!!================================================
      timeend                = 0.0d0
      timestart              = 0.0d0
      day                    = 0
!!
!!================================================
!! initialization of arrays for fit analysis 
!!================================================
      convvec(:,:,:)         = 0.0d0
      wshiftp(:)             = 0.0d0
      wshiftp2(:)            = 0.0d0
!!
!!================================================
!! determine starting time for epoch
!!================================================
      call abstime(timestart,day)
!!
!!================================================
!! get information on fixed weights from input.nn file
!!================================================
      allocate(wconstraintidxp(maxnum_weights_short_pair,npairs))
      call getwconstraintidx(0,npairs,windex_short_pair,&
        maxnum_layers_short_pair,num_layers_short_pair,maxnum_weights_short_pair,&
        num_weights_short_pair,num_weightspairfree,num_weightspairfixed,&
        nodes_short_pair,wconstraintidxp)
      call mpi_bcast(num_weightspairfree,npairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightspairfixed,npairs,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!================================================
!! get dimensions for the Kalman filter
!!================================================
      if((optmodee.eq.1).or.(optmodef.eq.1))then
        call getkaldims_short(npairs,&
          kaldim,corrdim,corrfdim,&
          maxkaldim,maxcorrdim,maxcorrfdim,&
          num_weightspairfree)
      endif
!!
!!================================================
!! allocate Kalman matrices
!!================================================
      allocate(corrmatrix_list(maxcorrdim,npairs))
      allocate(corrmatrixf_list(maxcorrfdim,npairs))
!!
!!================================================
!! get the values for the Kalman matrices and parameters (new or from file)
!!================================================
      call getkalmanmatrices_short(npairs,&
        iseed,kseed,lseed,mseed,&
        corrdim,corrfdim,&
        maxcorrdim,maxcorrfdim,num_weightspairfree,&
        corrmatrix_list,corrmatrixf_list,&
        kalmanlambdap)
!!
!!================================================
!! distribute random numbers (possibly read from file) 
!!================================================
      call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(kseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lseed,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mseed,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
!!================================================
!! count and distribute the number of training and test points
!!================================================
      if(mpirank.eq.0)then
        call countpoints(2,npairs,&
          ntrain,ntest,ntrainatoms,ntestatoms,trainpair,testpair)
      endif ! mpirank.eq.0
      call mpi_bcast(ntrain,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntest,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntrainatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntestatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(trainpair,npairs,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(testpair,npairs,mpi_integer,0,mpi_comm_world,mpierror)
!!
      allocate(pointindex(ntrain))
!!
!!================================================
!! get array num_atoms_all, contains number of atoms in each training point
!!================================================
      allocate(num_atoms_all(ntrain))
      call getnumatomsall(ntrain,num_atoms_all)
!!
!!================================================
!! allocate arrays for fitting statistics output
!!================================================
      allocate(fitstat(ntrain))
      fitstat(:)=0
      allocate(fitstatf(3,max_num_atoms,ntrain))
      fitstatf(:,:,:)=0
!!
!!================================================
!! if growth mode is not used make sure the full training set is used
!!================================================
      if(.not.lgrowth)then
        ngrowth=ntrain
      endif
!!
!!================================================
!! get maxcutoff_short_pair of short range symmetry function
!!================================================
      maxcutoff_short_pair=0.0d0
      do i2=1,npairs
        do i1=1,num_funcvalues_short_pair(i2)
          maxcutoff_short_pair=max(maxcutoff_short_pair,funccutoff_short_pair(i1,i2))
        enddo ! i1
      enddo
!!
!!================================================
!! get scaling.data 
!!================================================
      if(mpirank.eq.0)then
        if(luseoldscaling)then
          call readscale(npairs,2,&
            maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
            dummy,dummy,rdummy,rdummy) 
        else ! not luseoldscaling
!! get scaling data for the short range symmetry functions
          open(symunit,file='function.data',form='formatted',status='old')
          rewind(symunit) !'
          call getscale(npairs,max_num_pairs,2,maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            ntrain,symfunction_short_pair_list,minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
          close(symunit)
        endif ! luseoldscaling
!!
!!================================================
!! analyze structures in trainstruct.data and write chemical composition to runner.out
!!================================================
        if(lshort.and.(nelem.le.5).and.lanalyzecomposition)then
          call analyzeinput(0)
          call analyzeinput(1)
        elseif(lshort.and.(nelem.gt.5))then
          write(ounit,*)'WARNING: detailed analysis of structures not possible for more than 5 elements'
        endif !'
!!
!!================================================
!! get statistics of energies in training set (just for process 0)
!!================================================
        call getenergystatistics(2,belowmaxenergy,&
          eshortav,eshortstddev,eshortmin,eshortmax)
!!
!!================================================
!! get statistics of forces in training set (just for process 0)
!!================================================
        if((luseforces.and.lshort).or.lfinalforce)then
          call getforcestatistics(belowfmax,fmin,fmax,fvecmin,fvecmax)
        endif
!!
!!================================================
!! write scaling factors to file scaling.data
!!================================================
        if(lshort.and.(.not.luseoldscaling))then
          call writescale(npairs,2,&
            maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
            eshortmin,eshortmax,rdummy,rdummy)
        endif ! lshort
      endif ! mpirank.eq.0
!!
!!================================================
!! distribute scaling.data to all processes
!!================================================
      if(lshort)then
        call mpi_bcast(minvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue_short_pair,npairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortstddev,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!!================================================
!! write information on data sets to output file
!!================================================
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'number of training points ',ntrain !'
        write(ounit,*)'number of training atoms  ',ntrainatoms
        if(luseforces)then
          write(ounit,*)'number of training forces ',3*ntrainatoms
        endif
        write(ounit,*)'number of testing points  ',ntest
        write(ounit,*)'number of testing atoms   ',ntestatoms
        if(luseforces)then
          write(ounit,*)'number of testing forces  ',3*ntestatoms
        endif
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Number of pairs for each element combination:   '
        write(ounit,*)'                  training:    testing:   '
        npaircount=0
        do i1=1,nelem
         do i2=i1,nelem
          npaircount=npaircount+1
          write(ounit,'(i3,x,a3,x,a3,x,2i12)')npaircount,element(i1),element(i2),trainpair(npaircount),testpair(npaircount)
         enddo !'
        enddo
      endif ! mpirank.eq.0
!!
!!================================================
!! preparation of initial weight parameters'
!!================================================
!!
!!================================================
!! first initialize all weights randomly no matter if we restart a fit or not
!!================================================
      if(mpirank.eq.0)then
        call initialweights(npairs,&
          maxnum_weights_short_pair,num_weights_short_pair,&
          maxnum_layers_short_pair,num_layers_short_pair,windex_short_pair,nodes_short_pair,&
          jseed,nseed,weights_short_pair)
!!
!!================================================
!! if weights according to Nguyen Widrow are requested for short range NN, overwrite random weights:
!!================================================
        if(lnwweights)then
          call nguyenwidrowweights(npairs,maxnum_layers_short_pair,&
            nodes_short_pair,jseed,windex_short_pair,&
            maxnodes_short_pair,maxnum_weights_short_pair,num_weights_short_pair,&
            weights_short_pair,weights_min,weights_max,actfunc_short_pair)
        endif
!!
!!================================================
!! If systematic weights are requested, overwrite random weights for short range NN
!!================================================
        if(lsysweights)then
          call systematicweights(npairs,maxnum_layers_short_pair,&
            num_layers_short_pair,nodes_short_pair,&
            maxnum_weights_short_pair,&
            weights_short_pair,weights_min,weights_max)
        endif
!!
!!================================================
!! if requested overwrite short range weights with weights from weights.data
!!================================================
        if(luseoldweightsshort.and.lshort)then
          call readweights(2,npairs,&
            maxnum_weights_short_pair,num_weights_short_pair,weights_short_pair)
        endif
      endif ! mpirank.eq.0
!!
!!================================================
!! distribute weights_short_pair to all processes
!!================================================
      call mpi_bcast(weights_short_pair,npairs*maxnum_weights_short_pair,&
        mpi_real8,0,mpi_comm_world,mpierror)
!!
!!================================================
!! Caution: precondition is parallel inside
!!================================================
      if(lprecond)then
        call precondition_short_pair(&
         ntrain,trainpair,&
         minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
         eshortmin,eshortmax,eshortav,eshortstddev)
      endif ! lprecond
!!
!!================================================
!! write initial weights of epoch 0 to file
!!================================================
      if(mpirank.eq.0)then
        call writeweights(2,npairs,&
          maxnum_weights_short_pair,&
          maxnum_layers_short_pair,num_layers_short_pair,&
          nodes_short_pair,weights_short_pair)
      endif ! mpirank.eq.0

!!
!!================================================
!! keep copy of weights for the calculation of the weight change wshiftp
!!================================================
      weights_pair_veryold(:,:) =0.0d0
      weights_pair_old(:,:)     =weights_short_pair(:,:)
!!
!!================================================
!! end of preparation of initial weight parameters'
!!================================================
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!================================================
!! determine and print initialization time
!!================================================
      call abstime(timeend,day)
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.2)')&
          ' initialization time (min):',(timeend-timestart)/60.d0
      endif ! mpirank.eq.0
!!'
!!================================================
!! terminate RuNNer here if only initialization is requested 
!!================================================
      if(linionly)then
        write(ounit,*)'Initialization done - terminating RuNNer as requested'
        write(ounit,*)'-------------------------------------------------------------'
        return !'
      endif
!!
!!================================================
!! set time to zero for epoch 0 time measurement'
!!================================================
      timestart              =0.0d0
      day                    =0
      call abstime(timestart,day)
!!
!!================================================
!! initialize RMSEs
!!================================================
      rmse_short            =0.0d0
      rmse_short_test       =0.0d0
      rmse_force_s          =0.0d0
      rmse_force_s_test     =0.0d0
!!
!!================================================
!! initialize MADs
!!================================================
      mad_short             =0.0d0
      mad_short_test        =0.0d0
      mad_force_s           =0.0d0
      mad_force_s_test      =0.0d0
!!
      maxerror_eshort_train  = 0.0d0
      maxerror_eshort_test   = 0.0d0
      imaxerror_eshort_train = 0
      imaxerror_eshort_test  = 0
!!
!!================================================
!! calculate the initial training error
!!================================================
!!
      call geterror_short_pair(0,countepoch,ntrain,&
        imaxerror_eshort_train,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        rmse_short,rmse_force_s,mad_short,&
        mad_force_s,maxerror_eshort_train)
!!
!!================================================
!! calculate sensitivity for initial weights if requested
!!================================================
      if(lsens)then
        call getsensitivity_short(npairs,max_num_pairs,ntrain,trainpair,&
          maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
          scmin_short_pair,scmax_short_pair,sens)
      endif
!!
!!================================================
!! get new references for adaptive Kalman filter
!!================================================
      rmse_short_ref      =rmse_short
      rmse_force_s_ref    =rmse_force_s
!!
!!================================================
!! calculate the first test error
!!================================================
      call geterror_short_pair(1,countepoch,ntest,&
        imaxerror_eshort_test,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        rmse_short_test,rmse_force_s_test,mad_short_test,&
        mad_force_s_test,maxerror_eshort_test)
!!
!!================================================
!! set first references for determination of optweights files
!!================================================
      rmse_short_test_old = rmse_short_test
!!
      call abstime(timeend,day)
      timeepoch=timeend-timestart
!!
!!================================================
!! write headers
!!================================================
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------------------------'
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
              ' => Number of short range training force componentss below max_force: ',belowfmax(i1) 
          enddo
        endif
        write(ounit,*)'-------------------------------------------------------------------------------'
        if(lshort.and.lupdatebyelement.and.luseforces)then
          write(ounit,*)'### WARNING ### lupdatebyelement.eq.T => atom forces error refers only to this element'
          write(ounit,*)'                                      => short range energies are half meaningless'
          write(ounit,*)'-------------------------------------------------------------------------------'
        endif
        if(fitting_unit.eq.1)then
          write(ounit,*)'RMSEs (energies: eV/atom, forces: eV/Bohr):'
        elseif(fitting_unit.eq.2)then
          write(ounit,*)'RMSEs (energies: Ha/atom, forces: Ha/Bohr):'
        else
          write(ounit,*)'Error: please add new energy unit in fitting.f90'
          stop
        endif
        write(ounit,'(6a)')'                      --- E_short: --- ',&
                         ' - time -'
        write(ounit,'(6a)')'                          /atom      ',&
                         '   min'
        write(ounit,'(5a)')'       epoch         train         test' !'
!!================================================
!! write RMSEs 
!!================================================
        if(lprintdateandtime)then
          call printdateandtime(countepoch)
        endif
        write(ounit,'(a8,i5,x,2f13.6,f8.2)') ' ENERGY ',countepoch,&
          rmse_short*tounit,rmse_short_test*tounit,&
          timeepoch/60.d0
        if(luseforces)then
          write(ounit,'(a8,i5,x,2f13.6)')' FORCES ',countepoch,&
            rmse_force_s*tounit,rmse_force_s_test*tounit
        endif
        if(lprintmad)then
          write(ounit,'(a9,i5,2f14.6)')' MADE   ',countepoch,&
            mad_short*tounit,mad_short_test*tounit
          if(luseforces)then
            write(ounit,'(a8,i5,x,2f13.6)')' MADF   ',countepoch,&
              mad_force_s*tounit,mad_force_s_test*tounit
          endif
        endif
!!================================================
!! write convergence vector 
!!================================================
        if(lprintconv)then
          do i1=1,npairs
            write(ounit,'(a9,a2,x,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',&
              element(elempair(i1,1)),element(elempair(i1,2)),&
              countepoch,convvec(i1,1,1),convvec(i1,2,1),wshiftp(i1),wshiftp2(i1)
            convvec(i1,:,3)=convvec(i1,:,2)
            convvec(i1,:,2)=convvec(i1,:,1)
            convvec(i1,:,1)=0.0d0
          enddo ! i1'
        endif
!!================================================
!! print sensitivity if requested
!!================================================
        if(lsens)then
          write(ounit,*)'Short range NN sensitivity: '
          do i1=1,npairs
            do i2=1,num_funcvalues_short_pair(i1)
              write(ounit,'(i5,a15,x,a2,x,a2,x,i5,f16.8)')countepoch,' NNsensitivity ',&
                element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),i2,sens(i1,i2)  
            enddo ! i2
            write(ounit,*)'-------------------------------------------------------------'
          enddo ! i1
        endif ! lsens
      endif ! mpirank.eq.0
!!
!!================================================
!!================================================
!! loop over all training epochs
!!================================================
!!================================================
      do countepoch=1,nepochs
        if(mpirank.eq.0) write(debugunit,'(a6,i6)')'epoch ',countepoch
        call abstime(timestart,day)
!!
!!================================================
!! epoch-specific initializations
!!================================================
        ndone                 = 0 
        ncount = (int((countepoch-1)/growthstep)+1)*ngrowth    
        ncount = min(ncount,ntrain)
        point                 = 0
        numbere               = 0
        numberf               = 0
!!
!!================================================
!! get pointindex array for mixing structures randomly 
!!================================================
        call getpointindex(ncount,oseed,pointindex)
!!
!!================================================
!! loop block-wise over training points
!!================================================
11      continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!!
!!================================================
!! reset Kalman matrices if requested 
!!================================================
        if(lresetkalman)then
          if(mpisize.gt.1)then
            write(ounit,*)'Error: bcast for reinitialized corrmatrix still missing'
            stop !'
          endif
          do i1=1,npairs
            if(optmodee.eq.1)then
              call initialcorrmatrix(num_weightspairfree(i1),corrmatrix_list(1,i1))
            endif
            if(optmodef.eq.1)then
              call initialcorrmatrix(num_weightspairfree(i1),corrmatrixf_list(1,i1))
            endif
          enddo
        endif ! lresetkalman
!!
!!================================================
!! get index for random order of training points => idx
!!================================================
        call getidx(npoints,iseed,idx)
!!
!!================================================
!! read npoint data sets      
!! do all file reading for training data here at one place to allow for parallelization
!!================================================
        if(mpirank.eq.0)then
          call readfunctions_mixed(npoints,npairs,max_num_pairs,&
            ntrain,block_counter,pointindex,2,&
            maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
            symfunction_short_pair_list)
!!
!!================================================
!! read the structures needed for the calculation of the forces 
!! must be called after readfunctions because it needs num_atoms_list
!!================================================
          call getstructures_mixed(npoints,&
            ntrain,block_counter,pointindex,num_atoms_all)
!!
!!================================================
!! read short range forces from trainforces.data
!!================================================
          if(luseforces.or.(lfinalforce.and.(nepochs.eq.countepoch)))then
            call readforces_mixed(npoints,&
              ntrain,block_counter,pointindex,num_atoms_all)
          endif
!!
        endif ! mpirank.eq.0
!!
!!================================================
!! distribute the data to all processes
!!================================================
        call mpi_bcast(num_atoms_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_pairs_list,nblock,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(zelem_list,nblock*max_num_atoms,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(zelemp_list,2*nblock*max_num_pairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(shortenergy_list,nblock,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(symfunction_short_pair_list,&
          nblock*max_num_pairs*maxnum_funcvalues_short_pair,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(xyzstruct_list,nblock*max_num_atoms*3,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lattice_list,nblock*9,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lperiodic_list,nblock,mpi_logical,0,mpi_comm_world,mpierror)
        call mpi_bcast(shortforce_list,nblock*max_num_atoms*3,&
          mpi_real8,0,mpi_comm_world,mpierror)
!!
!!================================================
!! end of file reading for training data 
!!================================================
!!
!!================================================
!! determine which nstruct structures of this block should be scaled by this process
!!================================================
        call mpifitdistribution(npoints,nstruct,n_start,n_end)
        call mpi_barrier(mpi_comm_world,mpierror)
!!
!!================================================
!! start optimization of the short range part 
!!================================================
!!
!!================================================
!! scale symmetry functions for the short-range interaction
!!================================================
        call scalesymfit_parapair(&
          nstruct,n_start,n_end,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair)
!!
!!================================================
!! update weights using energies and forces
!!================================================
        call optimize_short_combinedpair(npoints,point,idx,&
          kaldim,maxcorrdim,maxcorrfdim,corrdim,corrfdim,countepoch,&
          num_weightspairfree,ndone,&
          wconstraintidxp,numbere,numberf,&
          lseed,ntrain,fitstat,fitstatf,kseed,&
          kalmanthreshold_temp,kalmanthresholdf_temp,&
          rmse_short_ref,rmse_force_s_ref,&
          corrmatrix_list,corrmatrixf_list,&
          minvalue_short_pair,maxvalue_short_pair)
!!
!!================================================
!! end optimization of the short range part 
!!================================================
!!
!!================================================
!! write temporary weights if requested
!!================================================
        if(mpirank.eq.0)then
          if(lwritetmpweights)then
            call writetmpweights(2,npairs,&
              maxnum_weights_short_pair,num_weights_short_pair,&
              weights_short_pair)
          endif
        endif ! mpirank.eq.0
!! 
!!================================================
!! count finished structures
!!================================================
        ndone=ndone+npoints
        block_counter=block_counter+npoints
!!
!!================================================
!! if there are structures left go back to next group of structures
!!================================================
        if(ncount.gt.0) goto 11 
!!================================================
!! end block-wise loop over training points
!!================================================
!!
!!================================================
!! write Kalman data for restart if requested
!!================================================
        if(lsavekalman)then
          if((optmodee.eq.1).or.(optmodef.eq.1))then
            call writekalman_short(npairs,&
              iseed,kseed,lseed,mseed,&
              maxcorrdim,kaldim,&
              num_weightspairfree,&
              kalmanlambdap,corrmatrix_list)
          endif ! optmode
        endif !lsavekalman
!!
!!================================================
!! write final weights of this epoch to files
!!================================================
        if(mpirank.eq.0)then
          if(mod(countepoch,iwriteweight).eq.0) then
            call getfilenamespair(countepoch)
            call writeweights(2,npairs,&
              maxnum_weights_short_pair,&
              maxnum_layers_short_pair,num_layers_short_pair,&
              nodes_short_pair,weights_short_pair)
          endif
        endif ! mpirank.eq.0
!!
!!================================================
!! initialize RMSEs for next epoch
!!================================================
        rmse_short            =0.0d0
        rmse_short_test       =0.0d0
        rmse_force_s          =0.0d0
        rmse_force_s_test     =0.0d0
!!
!!================================================
!! initialize MADs for next epoch
!!================================================
        mad_short             =0.0d0
        mad_short_test        =0.0d0
        mad_force_s           =0.0d0
        mad_force_s_test      =0.0d0
!!
        maxerror_eshort_train  = 0.0d0
        maxerror_eshort_test   = 0.0d0
        imaxerror_eshort_train = 0
        imaxerror_eshort_test  = 0
        block_counter          = 1
!!
!!================================================
!! If we are in the final epoch calculate the force error if requested
!!================================================
        if(lfinalforce.and.(nepochs.eq.countepoch))then
          luseforces=.true.
          forcernd=0.0d0
        endif
!!
!!================================================
!! calculate this epoch's training error
!!================================================
      call geterror_short_pair(0,countepoch,ntrain,&
        imaxerror_eshort_train,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        rmse_short,rmse_force_s,mad_short,&
        mad_force_s,maxerror_eshort_train)
!!
!!================================================
!! calculate sensitivity for this epoch if requested
!!================================================
      if(lsens)then
        call getsensitivity_short(npairs,max_num_pairs,ntrain,trainpair,&
          maxnum_funcvalues_short_pair,num_funcvalues_short_pair,&
          minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
          scmin_short_pair,scmax_short_pair,sens)
      endif
!!
!!================================================
!! get new references for adaptive Kalman filter
!!================================================
      rmse_short_ref   =rmse_short
      rmse_force_s_ref =rmse_force_s
!!
!!================================================
!! calculate the epoch's test error
!!================================================
      call geterror_short_pair(1,countepoch,ntest,&
        imaxerror_eshort_test,&
        minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        rmse_short_test,rmse_force_s_test,mad_short_test,&
        mad_force_s_test,maxerror_eshort_test)
!!
        call abstime(timeend,day)
        timeepoch=timeend-timestart
!!
!!================================================
!! calculate change of the short range weights for INFO output
!!================================================
        call getwshift(npairs,maxnum_weights_short_pair,num_weights_short_pair,&
          weights_short_pair,weights_pair_old,weights_pair_veryold,&
          wshiftp,wshiftp2)
!!
!!================================================
!! keep weights for calculation of wshiftp
!!================================================
        weights_pair_veryold(:,:)=weights_pair_old(:,:)
        weights_pair_old(:,:)=weights_short_pair(:,:)
!!
!!================================================
!! write RMSEs
!!================================================
        if(mpirank.eq.0)then
          if(lprintdateandtime)then
            call printdateandtime(countepoch)
          endif
          write(ounit,'(a8,i5,x,2f13.6,f8.2)')' ENERGY ', countepoch,&
            rmse_short*tounit,rmse_short_test*tounit,&
            timeepoch/60.d0
          if(luseforces)then
            write(ounit,'(a8,i5,x,2f13.6)')' FORCES ',countepoch,&
              rmse_force_s*tounit,rmse_force_s_test*tounit
          endif
!!================================================
!! write MADs
!!================================================
          if(lprintmad)then
            write(ounit,'(a9,i5,2f14.6)')' MADE   ',countepoch,&
              mad_short*tounit,mad_short_test*tounit
            if(luseforces)then
              write(ounit,'(a8,i5,x,2f13.6)')' MADF   ',countepoch,&
                mad_force_s*tounit,mad_force_s_test*tounit
            endif
          endif ! lrpintmad
          write(ounit,'(a,i5,2i10)')' INFORMATION USED FOR UPDATE (E,F) ',&
            countepoch,numbere,numberf  !'
!!================================================
!! write convergence vector 
!!================================================
          if(lprintconv)then
            do i1=1,npairs
              if(countepoch.eq.1)then
                convvec(i1,1,1)=wshiftp(i1)
                convvec(i1,2,1)=0.0d0
              else
                convtemp1=atan((convvec(i1,2,2)-convvec(i1,2,3))/(convvec(i1,1,2)-convvec(i1,1,3)))
                convtemp2=wshiftpold(i1) !abs((convvec(i1,1,3)-convvec(i1,1,2))/dcos(convtemp1))
                convtemp3=acos((-wshiftp(i1)**2+convtemp2**2+wshiftp2(i1)**2)/(2.d0*convtemp2*wshiftp2(i1)))
                convtemp1=convtemp1+convtemp3
                convvec(i1,1,1)=convvec(i1,1,3)+wshiftp2(i1)*dcos(convtemp1)
                convvec(i1,2,1)=convvec(i1,2,3)+wshiftp2(i1)*dsin(convtemp1)
              endif
              write(ounit,'(a9,a2,x,a2,x,i5,2f14.6,3x,2f20.10)')' CONVVEC ',&
                element(elempair(i1,1)),element(elempair(i1,2)),&
                countepoch,convvec(i1,1,1),convvec(i1,2,1),wshiftp(i1),wshiftp2(i1)
              convvec(i1,:,3)=convvec(i1,:,2)
              convvec(i1,:,2)=convvec(i1,:,1)
              convvec(i1,:,1)=0.0d0
              wshiftpold(i1)=wshiftp(i1)
            enddo ! i1'
          endif ! lprintconv
!!================================================
!! print sensitivity if requested
!!================================================
          if(lsens)then
            write(ounit,*)'Short range NN sensitivity: '
            do i1=1,npairs
              do i2=1,num_funcvalues_short_pair(i1)
                write(ounit,'(i5,a15,x,a2,x,a2,x,i5,f16.8)')countepoch,' NNsensitivity ',&
                  element(elementindex(elempair(i1,1))),element(elementindex(elempair(i1,2))),i2,sens(i1,i2)  
              enddo ! i2
              write(ounit,*)'-------------------------------------------------------------'
            enddo ! i1
          endif ! lsens
        endif ! mpirank.eq.0
!!'
!!================================================
!! if needed adjust Kalman filter parameters 
!!================================================
        call adjustkalman_short(numbere,numberf,&
          kalmanthreshold_temp,kalmanthresholdf_temp)
!!
!!================================================
!! write optimum set of weights to files optweights.out 
!!================================================
        if(rmse_short_test.le.rmse_short_test_old)then
          if(mpirank.eq.0)then
            optepoch=countepoch
            optweights_pair(:,:)=weights_short_pair(:,:)
            call writeoptweights(2,npairs,&
              maxnum_weights_short_pair,num_weights_short_pair,&
              optweights_pair)
          endif ! mpirank.eq.0
!!================================================
!! define new reference only if we found a new minimum rmse
!!================================================
          rmse_short_test_old= rmse_short_test
!!================================================
!! store best values obtained so far for final output
!!================================================
          optrmse_short       =rmse_short
          optrmse_short_test  =rmse_short_test
          optrmse_force_s     =rmse_force_s
          optrmse_force_s_test=rmse_force_s_test
        endif ! rmse_short_test
!!
      enddo ! countepoch, loop over all epochs
!!
!!================================================
!!================================================
!! end loop over all epochs
!!================================================
!!================================================
!!
!!================================================
!! summarize optimum fit
!!================================================
      if((mpirank.eq.0).and.(nepochs.gt.0))then
        call writeoptfit_short(optepoch,&
          tounit,optrmse_short,optrmse_short_test,optrmse_force_s,&
          optrmse_force_s_test)
      endif ! mpirank.eq.0
!!
!!================================================
!! report final statistics if requested:
!!================================================
      if(lfitstats.and.(mpirank.eq.0))then
        call writefitstat_short(ntrain,fitstat,fitstatf) 
      endif
!!
!!================================================
!! analyze error distribution if requested
!!================================================
      if(lanalyzeerror.and.(mpirank.eq.0))then
        call erroranalysis(ntrain,ntest)
      endif
!!
!!================================================
!! cleanup 
!!================================================
      if(allocated(wconstraintidxp)) deallocate(wconstraintidxp)  
      if(allocated(corrmatrix_list)) deallocate(corrmatrix_list)
      if(allocated(corrmatrixf_list))deallocate(corrmatrixf_list)
      if(allocated(fitstat))         deallocate(fitstat)
      if(allocated(fitstatf))        deallocate(fitstatf)
      if(allocated(pointindex))      deallocate(pointindex)
      if(allocated(num_atoms_all))   deallocate(num_atoms_all)
      deallocate(filenamewp)
      call deallocatestructures()
      if(lsens)then
        deallocate(sens)
      endif
!!
      return
      end
