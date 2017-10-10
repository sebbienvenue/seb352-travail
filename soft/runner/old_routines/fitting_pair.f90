!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine fitting(iunit,junit,kunit,munit,nunit,ounit,punit,&
         qunit,runit,sunit,tunit,uunit,vunit,dunit,max_num_atoms,nelem,listdim,&
         maxnum_weightsshort,maxnum_weightsewald,&
         num_weightsshort,num_weightsewald,iseed,nblock,paramode,&
         nepochs,nran,&
         maxnum_funcvalues,maxnum_funcvaluese,&
         num_funcvalues,num_funcvaluese,&
         function_type,function_typee,symelement,symeelement,&
         optmodee,optmodef,optmodeq,&
         maxnum_layersshort,maxnum_layersewald,&
         num_layersshort,num_layersewald,nodes_short,nodes_ewald,&
         maxnodes_short,maxnodes_ewald,ewaldkmax,elemupdate,&
         ngrowth,growthstep,fitting_unit,&
         windex,windexe,nenergygroup,nforcegroup,nchargegroup,&
         elementindex,nucelem,symfunction_list,symfunctione_list,&
         weights_short,weights_ewald,&
         kalmanthreshold,kalmanthresholdf,kalmanlambda,kalmannue,&
         kalmanthresholde,kalmanlambdae,kalmannuee,&
         kalmanthresholdc,kalmanlambdac,kalmannuec,&
         kalman_damp,kalman_dampf,kalman_dampq,&
         funccutoff,funccutoffe,eta,etae,zeta,zetae,&
         rshift,rshifte,lambda,lambdae,&
         worste,worstf,worstq,dampw,&
         weights_min,weights_max,weightse_min,weightse_max,&
         ewaldalpha,ewaldcutoff,forcernd,energyrnd,chargernd,&
         steepeststepe,steepeststepf,steepeststepq,scalefactorf,&
         actfunc_short,actfunc_ewald,element,pstring,&
         luseoldweightsshort,luseoldweightscharge,&
         lewald,lshort,lglobalfit,ljointefupdate,&
         lscalesym,lcentersym,lchargeconstraint,&
         lfinetime,lupdatebyelement,&
         lwritetrainpoints,lwritetraincharges,lrandomtrain,&
         luseforces,lwritetmpweights,lsavekalman,lrestkalman,&
         lwritetrainforces,luseworste,luseworstf,luseworstq,&
         lgrowth,ldampw,lfixweights,lfixedcharges,&
         lsysweights,lsysweightse,lreadunformatted,&
         lwriteunformatted,lresetkalman,lsepkalman,&
         lrepeate,ldebug)
!!
      use mpi_mod
      use times_mod
!!
      implicit none
!!
      integer icount
      integer dunit
      integer iunit
      integer junit
      integer kunit
      integer munit
      integer nunit
      integer ounit
      integer punit
      integer qunit
      integer runit
      integer sunit
      integer tunit
      integer uunit
      integer vunit
      integer max_num_atoms
      integer maxnum_weightsshort
      integer num_weightsshort(nelem)
      integer num_weightsshortfree(nelem)
      integer num_weightsshortfixed(nelem)
      integer maxnum_weightsewald
      integer num_weightsewald(nelem)
      integer num_weightsewaldfree(nelem)
      integer num_weightsewaldfixed(nelem)
      integer nelem
      integer iseed                                    !
      integer jseed                                    !
      integer kseed
      integer lseed
      integer mseed
      integer nblock
      integer npoints
      integer ncount
      integer point                                    ! internal
      integer pointe                                   ! internal
      integer point1
      integer nepochs
      integer maxnum_funcvalues                        ! in
      integer maxnum_funcvaluese                       ! in
      integer num_funcvalues(nelem)                    ! in
      integer num_funcvaluese(nelem)                   ! in
      integer elementindex(102)
      integer ntrain
      integer ntest
      integer nucelem(nelem)
      integer num_layersshort(nelem)
      integer num_layersewald(nelem)
      integer maxnum_layersshort
      integer maxnum_layersewald
      integer windex(2*maxnum_layersshort,nelem)
      integer windexe(2*maxnum_layersewald,nelem)
      integer zelem_list(nblock,max_num_atoms) 
      integer num_atoms_list(nblock)
      integer nodes_short(0:maxnum_layersshort,nelem)
      integer nodes_ewald(0:maxnum_layersewald,nelem)
      integer maxnodes_short
      integer maxnodes_ewald
      integer countepoch
      integer optmodee
      integer optmodef
      integer optmodeq
      integer day
      integer ewaldkmax                                  ! in
      integer listdim                                    ! in
      integer idx(nblock)                                ! internal
      integer nenergygroup                               ! in
      integer nforcegroup                                ! in
      integer nchargegroup                               ! in
      integer i1,i2,i3                                   ! internal 
      integer nstruct                                    ! internal
      integer function_type(maxnum_funcvalues,nelem)     ! in
      integer function_typee(maxnum_funcvaluese,nelem)   ! in
      integer symelement(maxnum_funcvalues,2,nelem)      ! in
      integer symeelement(maxnum_funcvaluese,2,nelem)    ! in
      integer n_start,n_end                              ! internal
      integer paramode                                   ! in
      integer numbere                                    ! internal
      integer numberf                                    ! internal
      integer numberq                                    ! internal
      integer ntrainatoms                                ! internal
      integer ntestatoms                                 ! internal
      integer trainelem(nelem)                           ! internal
      integer testelem(nelem)                            ! internal
      integer elemupdate                                 ! in
      integer dimlmq                                     ! internal
      integer ngrowth                                    ! in
      integer growthstep                                 ! in
      integer wconstraint(maxnum_weightsshort,nelem)     ! internal
      integer wconstrainte(maxnum_weightsewald,nelem)    ! internal
      integer, dimension(:,:)  , allocatable :: wconstraintidx ! internal 
      integer, dimension(:,:)  , allocatable :: wconstraintidxe ! internal 
      integer fitting_unit                               ! in
      integer isum                                       ! internal
      integer nran                                       ! in
!! Kalman matrix dimensions:
      integer corrdim(nelem)
      integer corrfdim(nelem)
      integer corredim(nelem)
      integer corrcdim
      integer kaldim(nelem)
      integer kaledim(nelem)
      integer kalcdim
      integer maxkaldim
      integer maxkaledim
      integer maxcorrdim
      integer maxcorrfdim
      integer maxcorredim
!!
      real*8 weights_short(maxnum_weightsshort,nelem)
      real*8 weights_short_old(maxnum_weightsshort,nelem)
      real*8 weights_ewald(maxnum_weightsewald,nelem)
      real*8 weights_ewald_save(maxnum_weightsewald,nelem)
      real*8 weights_ewald_old(maxnum_weightsewald,nelem)
      real*8 optweights_short(maxnum_weightsshort,nelem)
      real*8 optweights_ewald(maxnum_weightsewald,nelem)
      real*8 symfunction_list(maxnum_funcvalues,max_num_atoms,nblock)
      real*8 symfunctione_list(maxnum_funcvaluese,max_num_atoms,nblock)
      real*8 minvalue(nelem,maxnum_funcvalues)
      real*8 minvaluee(nelem,maxnum_funcvaluese)
      real*8 maxvalue(nelem,maxnum_funcvalues)
      real*8 maxvaluee(nelem,maxnum_funcvaluese)
      real*8 avvalue(nelem,maxnum_funcvalues)
      real*8 avvaluee(nelem,maxnum_funcvaluese)
      real*8 xyzstruct_list(3,max_num_atoms,nblock)                  ! internal
      real*8 lattice_list(3,3,nblock)
      real*8 weights_min                              ! in
      real*8 weights_max                              ! in
      real*8 weightse_min                             ! in
      real*8 weightse_max                             ! in
!! symmetry function parameters
      real*8 funccutoff(maxnum_funcvalues,nelem)    ! in
      real*8 funccutoffe(maxnum_funcvaluese,nelem)  ! in
      real*8 eta(maxnum_funcvalues,nelem)           ! in
      real*8 etae(maxnum_funcvaluese,nelem)         ! in
      real*8 lambda(maxnum_funcvalues,nelem)        ! in
      real*8 lambdae(maxnum_funcvaluese,nelem)      ! in
      real*8 zeta(maxnum_funcvalues,nelem)          ! in
      real*8 zetae(maxnum_funcvaluese,nelem)        ! in
      real*8 rshift(maxnum_funcvalues,nelem)        ! in
      real*8 rshifte(maxnum_funcvaluese,nelem)      ! in
      real*8 maxcutoff                                ! internal
      real*8 maxcutoffe                               ! internal
!! DFT data:
      real*8 totalcharge_list(nblock)
      real*8 totalenergy_list(nblock)
      real*8 shortenergy_list(nblock)                                ! internal
      real*8 ewaldenergy_list(nblock)
      real*8 atomcharge_list(nblock,max_num_atoms)
      real*8 xyzforce_list(3,max_num_atoms,nblock)
      real*8 ewaldforce_list(3,max_num_atoms,nblock)
!! NN data:
      real*8 nneshort_list(nblock)
      real*8 nnewald_list(nblock)
      real*8 nnetot_list(nblock)
      real*8 nnewaldforce_list(3,max_num_atoms,nblock)               ! internal
!!
      real*8 rmse_short
      real*8 rmse_charge
      real*8 rmse_charge_old
      real*8 rmse_totalcharge
      real*8 rmse_etot
      real*8 rmse_ewald
      real*8 rmse_short_test
      real*8 rmse_charge_test
      real*8 rmse_charge_test_old
      real*8 rmse_totalcharge_test
      real*8 rmse_etot_test
      real*8 rmse_ewald_test
      real*8 rmse_etot_old
      real*8 rmse_etot_test_old
      real*8 rmse_force_s
      real*8 rmse_force_s_test
      real*8 rmse_force_e
      real*8 rmse_force_e_test
      real*8 rmse_short_ref
      real*8 rmse_charge_ref
      real*8 rmse_force_s_ref
      real*8 rmse_totalcharge_ref
      real*8 epochtime
      real*8 timestart
      real*8 timeend
      real*8 kalmanthreshold
      real*8 kalmanthresholdf
      real*8 kalmanlambda(nelem)
      real*8 kalmanlambdaf(nelem)
      real*8 kalmannue
      real*8 kalmanthresholde
      real*8 kalmanlambdae(nelem)
      real*8 kalmannuee
      real*8 kalmanthresholdc
      real*8 kalmanlambdac
      real*8 kalmannuec
      real*8 kalman_damp
      real*8 kalman_dampf
      real*8 kalman_dampq
      real*8, dimension(:,:), allocatable :: corrmatrix_list
      real*8, dimension(:,:), allocatable :: corrmatrixf_list
      real*8, dimension(:,:), allocatable :: corrmatrixe_list
      real*8, dimension(:)  , allocatable :: corrmatrixc 
      real*8 tounit
      real*8 ewaldalpha                                       ! in
      real*8 ewaldcutoff                                      ! in
      real*8 forcernd                                         ! in
      real*8 energyrnd                                        ! in
      real*8 chargernd                                        ! in
      real*8 wshift                                           ! internal
      real*8 wshifte                                          ! internal
      real*8 steepeststepe                                    ! in
      real*8 steepeststepf                                    ! in
      real*8 steepeststepq                                    ! in
      real*8 scalefactorf                                     ! in(out)
      real*8 worste                                           ! in
      real*8 worstf                                           ! in
      real*8 worstq                                           ! in
      real*8 dampw                                            ! in
!! Levenberg-Marquardt variables
      real*8, dimension(:,:,:)  , allocatable :: alphalmq 
      real*8, dimension(:,:,:)  , allocatable :: covarlmq 
      real*8, dimension(:,:)  , allocatable :: betalmq 
      real*8, dimension(:,:)  , allocatable :: dalmq 
      real*8 chisqq(nelem)
      real*8 ochisqq(nelem)
      real*8 stddevinvq
      real*8 alambdalmq(nelem)
      real*8 avcharge(nelem)
      real*8 stddevcharge(nelem)
      real*8 eshortstddev                                    ! internal
!!
      character*1 actfunc_short(maxnodes_short,maxnum_layersshort,nelem) ! in
      character*1 actfunc_ewald(maxnodes_ewald,maxnum_layersewald,nelem) ! in
      character*20 filenamews(nelem)
      character*20 filenamewe(nelem)
      character*20 filenametemp
      character*20 pstring                                   ! in
      character*2 element(nelem)                             ! in
!!
      logical luseoldweightsshort                            ! in
      logical luseoldweightscharge                           ! in
      logical lewald                                         ! in
      logical lshort                                         ! in
      logical lscalesym                                      ! in
      logical lcentersym                                     ! in
      logical ldebug                                         ! in
      logical lperiodic_list(nblock)
      logical lchargeconstraint                              ! in
      logical lglobalfit                                     ! in
      logical lfinetime                                      ! in
      logical lwritetrainpoints                              ! in
      logical lwritetraincharges                             ! in
      logical lwritetrainforces                              ! in
      logical lrandomtrain                                   ! in
      logical luseforces                                     ! in
      logical lwritetmpweights                               ! in
      logical lsavekalman                                    ! in
      logical lrestkalman                                    ! in
      logical lupdatebyelement                               ! in
      logical luseworste                                     ! in
      logical luseworstf                                     ! in
      logical luseworstq                                     ! in
      logical lgrowth                                        ! in
      logical ldampw                                         ! in
      logical lfixweights                                    ! in
      logical lfixedcharges                                  ! in
      logical ljointefupdate                                 ! in
      logical lsysweights                                    ! in
      logical lsysweightse                                   ! in
      logical lrmin                                          ! in
      logical lreadunformatted                               ! in
      logical lwriteunformatted                              ! in
      logical lresetkalman                                   ! in
      logical lsepkalman                                     ! in
      logical lrepeate                                       ! in
!!
!! initialization
      kalmanlambdaf(:)       = kalmanlambda(:)
      minvalue(:,:)          = 0.0d0
      minvaluee(:,:)         = 0.0d0
      maxvalue(:,:)          = 0.0d0
      maxvaluee(:,:)         = 0.0d0
      avvalue(:,:)           = 0.0d0
      avvaluee(:,:)          = 0.0d0
      filenamews(:)          ='000000.short.000.out'
      filenamewp(:)          ='000000.pair.000.000.out'      ! Jovan: fix here
      filenamewe(:)          ='000000.ewald.000.out'
      do i1=1,nelem
        filenametemp=filenamews(i1)
        if(nucelem(i1).gt.99)then
          write(filenametemp(14:16),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(filenametemp(15:16),'(i2)')nucelem(i1)
        else
          write(filenametemp(16:16),'(i1)')nucelem(i1)
        endif
        filenamews(i1)=filenametemp
        filenametemp=filenamewe(i1)
        if(nucelem(i1).gt.99)then
          write(filenametemp(14:16),'(i3)')nucelem(i1)
        elseif(nucelem(i1).gt.9)then
          write(filenametemp(15:16),'(i2)')nucelem(i1)
        else
          write(filenametemp(16:16),'(i1)')nucelem(i1)
        endif
        filenamewe(i1)=filenametemp
      enddo
      countepoch             =0
      totalcharge_list(:)    =0.0d0
      totalenergy_list(:)    =0.0d0
      shortenergy_list(:)    =0.0d0
      ewaldenergy_list(:)    =0.0d0
      nnewald_list(:)        =0.0d0
      nneshort_list(:)       =0.0d0
      nnetot_list(:)         =0.0d0
      point                  =0
      pointe                 =0
      point1                 =0
      if(fitting_unit.eq.1)then
        tounit               = 27.211d0   ! energy conversion Ha to eV 
      elseif(fitting_unit.eq.2)then
        tounit               = 1.d0   !  stay with Ha 
      endif
      timeend                = 0.0d0
      xyzstruct_list(:,:,:)  = 0.0d0
      xyzforce_list(:,:,:)   = 0.0d0
      ewaldforce_list(:,:,:) = 0.0d0
      nnewaldforce_list(:,:,:)=0.0d0
      lattice_list(:,:,:)    = 0.0d0
      npoints                = 0
      kseed                  = iseed 
      lseed                  = iseed 
      mseed                  = iseed 
      jseed                  = iseed 
      maxcutoff              = 0.0d0 
      maxcutoffe             = 0.0d0 
      rmse_short_ref         = 0.0d0
      rmse_charge_ref        = 0.0d0
      rmse_totalcharge_ref   = 0.0d0
      rmse_force_s_ref       = 0.0d0
      numbere                = 0
      numberf                = 0
      numberq                = 0
      time1                  = 0.0d0
      time2                  = 0.0d0
      time3                  = 0.0d0
      time4                  = 0.0d0
!!
!! get process id mpirank
      call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)
!!
!! check parallel job size
      call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
!!
!!      write(ounit,*)mpirank,' fitting starts'
!!
!! get information on fixed weights
      if(mpirank.eq.0)then
        wconstraint(:,:)=0
        wconstrainte(:,:)=0
        num_weightsshortfree(:)=num_weightsshort(:)
        num_weightsshortfixed(:)=0
        num_weightsewaldfree(:)=num_weightsewald(:)
        num_weightsewaldfixed(:)=0
!!
!! FIXME later: read weight constraints here later
!! write summary
        do i1=1,nelem
          write(ounit,*)' Element ',element(i1)
          write(ounit,*)'  Total number of short range weights       ',num_weightsshort(i1)
          write(ounit,*)'  Number of optimized short range weights   ',num_weightsshortfree(i1)
          write(ounit,*)'  Number of constrained short range weights ',num_weightsshortfixed(i1)
          write(ounit,*)'  Total number of charge weights            ',num_weightsewald(i1)
          write(ounit,*)'  Number of optimized charge weights        ',num_weightsewaldfree(i1)
          write(ounit,*)'  Number of constrained charge weights      ',num_weightsewaldfixed(i1)
        enddo ! i1
      endif ! mpirank.eq.0
!! distribute constraints to all nodes here
      call mpi_bcast(wconstraint,maxnum_weightsshort*nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(wconstrainte,maxnum_weightsewald*nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsshortfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsshortfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsewaldfree,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(num_weightsewaldfixed,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!! prepare constraint arrays
      allocate(wconstraintidx(maxnum_weightsshort,nelem))  
      allocate(wconstraintidxe(maxnum_weightsewald,nelem))  
!! get wconstraintidx
      do i1=1,nelem
        icount=0
        do i2=1,num_weightsshort(i1)
          if(wconstraint(i2,i1).eq.0)then
            icount=icount+1
            wconstraintidx(icount,i1)=i2
          endif
        enddo
      enddo
!! get wconstraintidxe
      do i1=1,nelem
        icount=0
        do i2=1,num_weightsewald(i1)
          if(wconstrainte(i2,i1).eq.0) then
            icount=icount+1
            wconstraintidxe(icount,i1)=i2
          endif
        enddo
      enddo
!!
!! get dimensions for the Kalman filter
      if((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))then
        call getkaldims(ounit,nelem,optmodee,optmodef,optmodeq,&
          kaldim,kaledim,kalcdim,corrdim,corrfdim,corredim,corrcdim,&
          maxkaldim,maxkaledim,maxcorrdim,maxcorrfdim,maxcorredim,&
          num_weightsshortfree,num_weightsewaldfree,&
          lshort,lewald,lchargeconstraint)
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
        call getkalmanmatrices(ounit,iunit,nucelem,&
          iseed,kseed,lseed,mseed,&
          kaldim,kaledim,kalcdim,corrdim,corrfdim,corredim,corrcdim,&
          maxkaldim,maxkaledim,maxcorrdim,maxcorrfdim,maxcorredim,&
          optmodee,optmodef,optmodeq,&
          nelem,num_weightsshortfree,num_weightsewaldfree,&
          corrmatrix_list,corrmatrixf_list,corrmatrixe_list,corrmatrixc,&
          kalmanlambda,kalmanlambdae,kalmanlambdac,&
          lshort,lewald,lchargeconstraint,lrestkalman)
!!
!! count the number of training and test points
      if(mpirank.eq.0)then
        call countpoints(iunit,ounit,nelem,elementindex,&
          ntrain,ntest,ntrainatoms,ntestatoms,&
          trainelem,testelem,ldebug)
      endif ! mpirank.eq.0
!!
      call mpi_bcast(ntrain,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntest,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntrainatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ntestatoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(trainelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(testelem,nelem,mpi_integer,0,mpi_comm_world,mpierror)
!!
!! if growth mode is not used make sure the full training set is used
      if(.not.lgrowth)then
        ngrowth=ntrain
      endif
!!
!! get maxcutoff
      maxcutoff  =0.0d0
      do i2=1,nelem
        do i1=1,num_funcvalues(i2)
          maxcutoff=max(maxcutoff,funccutoff(i1,i2))
        enddo ! i1
      enddo
!! get maxcutoffe
      maxcutoffe =0.0d0
      do i2=1,nelem
        do i1=1,num_funcvaluese(i2)
          maxcutoffe=max(maxcutoffe,funccutoffe(i1,i2))
        enddo ! i1
      enddo ! i2
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! get scaling.data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! get scaling data for the short range symmetry functions
      if(mpirank.eq.0)then
        if(lshort)then
          open(iunit,file='function.data',form='formatted',status='old')
          rewind(iunit)
          call getscale(iunit,ounit,nelem,0,&
            nblock,maxnum_funcvalues,num_funcvalues,&
            max_num_atoms,elementindex,&
            ntrain,nucelem,symfunction_list,minvalue,&
            maxvalue,avvalue,element,lscalesym,lshort,lewald,ldebug)
          close(iunit)
        endif ! lshort
!!
!! get scaling data for the charge symmetry functions 
        if(lewald)then
          open(iunit,file='functione.data',form='formatted',status='old')
          rewind(iunit)
          call getscale(iunit,ounit,nelem,1,&
            nblock,maxnum_funcvaluese,num_funcvaluese,&
            max_num_atoms,elementindex,&
            ntrain,nucelem,symfunctione_list,minvaluee,&
            maxvaluee,avvaluee,element,lscalesym,lshort,lewald,ldebug)
          close(iunit)
        endif ! lewald
!!
!! get statistics of energies in training set (just for process 0)
        call getenergystatistics(iunit,ounit,nelem,nucelem,&
          elementindex,eshortstddev,element,lshort,lewald,ldebug)
!!
!! get statistics of forces in training set (just for process 0)
        if(luseforces.and.lshort)then
          call getforcestatistics(iunit,junit,ounit,nelem,nucelem,&
            elementindex,element,ldebug)
        endif
!!
!! write scaling factors to file scaling.data
        if(lshort)then
          call writescale(iunit,ounit,nelem,&
             maxnum_funcvalues,num_funcvalues,&
             minvalue,maxvalue,avvalue,&
             .true.,ldebug)
        endif ! lshort
        if(lewald)then
          call writescale(iunit,ounit,nelem,&
             maxnum_funcvaluese,num_funcvaluese,&
             minvaluee,maxvaluee,avvaluee,&
             .false.,ldebug)
         endif ! lewald
      endif ! mpirank.eq.0
!!
!! distribute scaling.data to all processes
      if(lshort)then
        call mpi_bcast(minvalue,nelem*maxnum_funcvalues,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvalue,nelem*maxnum_funcvalues,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvalue,nelem*maxnum_funcvalues,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eshortstddev,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lewald)then
        call mpi_bcast(minvaluee,nelem*maxnum_funcvaluese,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxvaluee,nelem*maxnum_funcvaluese,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(avvaluee,nelem*maxnum_funcvaluese,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
!! get the average charges
      if(lewald)then
        avcharge(:)=0.0d0
        if(mpirank.eq.0)then
          call getavcharge(ounit,iunit,junit,nelem,ntrain,&
            elementindex,nucelem,avcharge,stddevcharge,&
            ldebug)
        endif ! mpirank.eq.0
        call mpi_bcast(avcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(stddevcharge,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif ! lewald
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
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
        write(ounit,*)'Number of atoms for each element:   '
        write(ounit,*)'           training:    testing:   '
        do i1=1,nelem
          write(ounit,'(i3,x,a3,x,2i12)')i1,element(i1),trainelem(i1),testelem(i1)
        enddo
      endif ! mpirank.eq.0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! get initial weight parameters'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! first initialize all weights randomly no matter if we restart a fit or not
      if(mpirank.eq.0)then
        call initialweights(nelem,nran,&
             maxnum_weightsshort,maxnum_weightsewald,&
             num_weightsshort,num_weightsewald,&
             jseed,weights_short,weights_ewald,&
             weights_min,weights_max,weightse_min,weightse_max,&
             lewald,lshort)
!!
!! if requested overwrite short range weights with weights from weights.data
        if(luseoldweightsshort.and.lshort)then
            call readweights(iunit,nelem,nucelem,&
              maxnum_weightsshort,num_weightsshort,&
              weights_short,.false.,.true.)
        endif
!!
!! if requested overwrite charge weights with weights from weightse.data
        if(luseoldweightscharge.and.lewald)then
            call readweights(iunit,nelem,nucelem,&
              maxnum_weightsewald,num_weightsewald,&
              weights_ewald,.true.,.false.)
        endif
!!
!! write initial weights to file
        if(lshort)then
            call writeweights(iunit,nelem,&
              maxnum_weightsshort,num_weightsshort,&
              maxnum_layersshort,num_layersshort,&
              nodes_short,weights_short,filenamews)
        endif ! lshort
        if(lewald)then
            call writeweights(iunit,nelem,&
              maxnum_weightsewald,num_weightsewald,&
              maxnum_layersewald,num_layersewald,&
              nodes_ewald,weights_ewald,filenamewe)
        endif ! lewald
      endif ! mpirank.eq.0
!!
!! keep copy of weights
      weights_short_old(:,:)=weights_short(:,:)
      weights_ewald_old(:,:)=weights_ewald(:,:)
!!
!! distribute weights_short and weights_ewald to all processes
      if(lshort)then
        call mpi_bcast(weights_short,nelem*maxnum_weightsshort,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lewald)then
        call mpi_bcast(weights_ewald,nelem*maxnum_weightsewald,&
          mpi_real8,0,mpi_comm_world,mpierror)
      endif
      call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,mpierror)
!!
      call mpi_barrier(mpi_comm_world,mpierror)
!!
!!'
      call abstime(timeend,day)
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.2)')&
          ' initialization time (min):',(timeend-timestart)/60.d0
      endif ! mpirank.eq.0
!!'
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
      rmse_force_e           = 0.0d0
      rmse_force_e_test      = 0.0d0
!!



!! if you are here => very good!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the initial training error
!! (rmse_short,rmse_charge,rmse_totalcharge,rmse_ewald,rmse_etot,rmse_etot_old,rmse_charge_old,rmse_force_s,rmse_force_e)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror(0,ounit,iunit,countepoch,&
        junit,kunit,munit,nunit,qunit,runit,sunit,&
        ntrain,nblock,max_num_atoms,&
        maxnum_funcvalues,maxnum_funcvaluese,&
        num_funcvalues,num_funcvaluese,&
        nelem,elementindex,function_type,&
        function_typee,symelement,symeelement,windex,windexe,&
        maxnum_layersshort,maxnum_layersewald,&
        num_layersshort,num_layersewald,&
        maxnodes_short,maxnodes_ewald,nodes_short,nodes_ewald,&
        maxnum_weightsshort,maxnum_weightsewald,&
        num_weightsshort,num_weightsewald,listdim,nucelem,&
        ewaldkmax,elemupdate,&
        funccutoff,funccutoffe,eta,etae,zeta,zetae,&
        rshift,rshifte,lambda,lambdae,maxcutoff,maxcutoffe,&
        minvalue,maxvalue,avvalue,minvaluee,maxvaluee,avvaluee,&
        weights_short,weights_ewald,ewaldalpha,ewaldcutoff,&
        rmse_short,rmse_charge,rmse_totalcharge,&
        rmse_ewald,rmse_etot,rmse_force_s,&
        rmse_force_e,&
        rmse_etot_old,rmse_charge_old,&
        actfunc_short,actfunc_ewald,pstring,&
        luseforces,ldebug,lscalesym,lcentersym,&
        lshort,lewald,lwritetraincharges,lwritetrainpoints,&
        lwritetrainforces,lfixedcharges,lupdatebyelement)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of the first training error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!! get new references for adaptive Kalman filter
      if(lshort) rmse_short_ref      =rmse_short
      if(lewald) rmse_charge_ref     =rmse_charge
      if(lewald) rmse_totalcharge_ref=rmse_totalcharge
      rmse_force_s_ref               =rmse_force_s
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! calculate the first testing error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call geterror(1,ounit,iunit,countepoch,&
        junit,kunit,munit,nunit,qunit,runit,sunit,&
        ntest,nblock,max_num_atoms,&
        maxnum_funcvalues,maxnum_funcvaluese,&
        num_funcvalues,num_funcvaluese,&
        nelem,elementindex,function_type,&
        function_typee,symelement,symeelement,windex,windexe,&
        maxnum_layersshort,maxnum_layersewald,&
        num_layersshort,num_layersewald,&
        maxnodes_short,maxnodes_ewald,nodes_short,nodes_ewald,&
        maxnum_weightsshort,maxnum_weightsewald,&
        num_weightsshort,num_weightsewald,listdim,nucelem,&
        ewaldkmax,elemupdate,&
        funccutoff,funccutoffe,eta,etae,zeta,zetae,&
        rshift,rshifte,lambda,lambdae,maxcutoff,maxcutoffe,&
        minvalue,maxvalue,avvalue,minvaluee,maxvaluee,avvaluee,&
        weights_short,weights_ewald,ewaldalpha,ewaldcutoff,&
        rmse_short_test,rmse_charge_test,rmse_totalcharge_test,&
        rmse_ewald_test,rmse_etot_test,rmse_force_s_test,&
        rmse_force_e_test,&
        rmse_etot_test_old,rmse_charge_test_old,&
        actfunc_short,actfunc_ewald,pstring,&
        luseforces,ldebug,lscalesym,lcentersym,&
        lshort,lewald,lwritetraincharges,lwritetrainpoints,&
        lwritetrainforces,lfixedcharges,lupdatebyelement)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! end calculation of first test error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!
!!      write(ounit,*)mpirank,' rmse_force_s ',rmse_force_s
!!      write(ounit,*)mpirank,' rmse_force_s_test ',rmse_force_s_test
!!
      call abstime(timeend,day)
      epochtime=timeend-timestart
!!
!! write RMSE headers
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Did you check your output file for warnings? ;-)             '
        write(ounit,*)'-------------------------------------------------------------------------------'
        if(lewald.and.lupdatebyelement)then
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
          write(ounit,*)'Error: please add new energy unit in fitting.f90'
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
          write(ounit,'(a8,i5,x,4f13.6)')' FORCES ',countepoch,&
            rmse_force_s*tounit,rmse_force_s_test*tounit,&
            rmse_force_e*tounit,rmse_force_e_test*tounit 
        endif
      endif ! mpirank.eq.0
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop over all training epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do countepoch=1,nepochs
!!
        write(ounit,*)'Fitting starts now'
        stop

      enddo ! countepoch, loop over all epochs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end loop over all epochs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(allocated(wconstraintidx))  deallocate(wconstraintidx)  
      if(allocated(wconstraintidxe)) deallocate(wconstraintidxe)  
!! deallocate Kalman matrices
      if(allocated(corrmatrix_list)) deallocate(corrmatrix_list)
      if(allocated(corrmatrixf_list))deallocate(corrmatrixf_list)
      if(allocated(corrmatrixe_list))deallocate(corrmatrixe_list)
      if(allocated(corrmatrixc))     deallocate(corrmatrixc) 
!!
      return
      end
