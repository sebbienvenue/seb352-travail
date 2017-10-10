!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine inputnndefaults()
!!
      use nnflags 
      use globaloptions 
      use mode1options
      use predictionoptions
      use fittingoptions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use runneripi
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        nodes_short_atomic(:,:)=0
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        nodes_short_pair(:,:)=0
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        nodes_elec(:,:)=0
      endif
      if(lnntb)then
        if(nntb_flag(0)) then
          nodes_ham(:,:)=0
        endif
        if(nntb_flag(1)) then
          nodes_s(:,:)=0
        endif
        if(nntb_flag(2)) then
          nodes_hexton(:,:)=0
        endif
        if(nntb_flag(3)) then
          if(mode.ge.2)then
            nodes_hextoff(:,:)=0
          endif
        endif
        if(nntb_flag(4)) then
          nodes_dens(:,:)=0
        endif

      endif
      analyze_error_energy_step = 0.01d0
      analyze_error_force_step = 0.01d0
      analyze_error_charge_step = 0.001d0
      ldebug        =.false.
      paramode      =1
      ewaldalpha=0.0d0
      ewaldcutoff=0.0d0
      ewaldkmax=0
      nenergygroup=1
      nforcegroup=1
      nchargegroup=1
      luseforces=.false.
      energyrnd=1.0d0
      forcernd=1.0d0
      chargernd=1.0d0
      luseatomcharges=.false.
      luseatomenergies=.false.
      lremoveatomenergies=.false.
      lchargeconstraint=.false.
      lfitethres=.false.
      fitethres=0.0d0
      lfitfthres=.false.
      fitfthres=0.0d0
      rmin=0.5d0
      optmodee=1
      optmodef=1
      optmodeq=1
      optmodeham=1
      optmodes=1
      optmodedens=1
      optmodehexton=1
      optmodehextoff=1
      nblock=200
      nepochs=0
      iwriteweight=1
      lwritetmpweights=.false.
      kalmanthreshold=0.0d0
      kalmanthresholdf=0.0d0
      kalmanthresholde=0.0d0
      kalmanthresholdc=0.0d0
      kalman_dampe=1.0d0
      kalman_dampf=1.0d0
      kalman_dampq=1.0d0
      steepeststepe=0.01d0
      steepeststepf=0.01d0
      steepeststepq=0.01d0
      scalefactorf=1.d0
      lrandomtrain=.false.
      lscalesym=.false.
      lcentersym=.false.
      luseoldweightsshort=.false.
      luseoldweightscharge=.false.
      luseoldweights_s=.false.
      luseoldweights_hexton=.false.
      luseoldweights_hextoff=.false.
      luseoldweights_dens=.false.
      lglobalfit=.false.
      lsavekalman=.false.
      lrestkalman=.false.
      lupdatebyelement=.false.
      luseworste=.false.
      luseworstf=.false.
      luseworstq=.false.
      lgrowth=.false.
      ngrowth=0
      growthstep=1
      ldampw=.false.
      dampw=0.0d0
      lfixweights=.false.
      ldoforces=.false.
      ldostress=.false.
      lfinetime=.false.
      lfinetimeepoch=.false.
      lwritepdb=.false.
      lwritexyz=.false.
      lwritepov=.false.
      lwritepw=.false.
      lwritetrainpoints=.false.
      lwritetrainforces=.false.
      lwritetraincharges=.false.
      atomrefenergies(:)=0.0d0
      weights_min=-1.d0
      weights_max=1.d0
      biasweights_min=-1.d0
      biasweights_max=1.d0
      weightse_min=-1.d0
      weightse_max=1.d0
      fitting_unit=1
      ljointefupdate=.false.
      pstring='00000000000000000000'
!!      nn_type_short =0
!!      nn_type_elec=0             ! no electrostatics
!!      nn_type_nntb =0
      nran    =1
      enforcetotcharge=0
      fixedcharge(:)=99.0d0
      lsysweights=.false.
      lsysweightse=.false.
      lsens=.false.
      lreadunformatted=.false.
      lwriteunformatted=.false.
      lresetkalman=.false.
      lsepkalman=.false.
      lrepeate=.false.
      maxforce=10000.d0
      maxenergy=10000.d0
      lfinalforce=.false.
      lcheckf=.false.
      lfitstats=.false.
      lfixederrore=.false.
      lfixederrorf=.false.
      lompmkl=.false.
      lnormnodes=.false.
      restrictw=-100000.d0
      fitmode=1            ! default is online learning
      lanalyzeerror=.false.
      lnwweights=.false.
      lnwweightse=.false.
      scmin_short_atomic=0.0d0  
      scmax_short_atomic=1.0d0
      scmin_short_pair=0.0d0  
      scmax_short_pair=1.0d0
      scmin_ham=0.0d0  
      scmax_ham=1.0d0
      scmin_s=0.0d0
      scmax_s=1.0d0
      scmin_hexton=0.0d0
      scmax_hexton=1.0d0
      scmin_hextoff=0.0d0
      scmax_hextoff=1.0d0
      scmin_dens=0.0d0
      scmax_dens=1.0d0

      scmin_elec=0.0d0     
      scmax_elec=1.0d0    
      luseoldscaling=.false.
      lprecond=.false.
      linionly=.false.
      noisee=0.0d0
      noisef=0.0d0
      noiseq=0.0d0
      lprintconv=.false.
      lprintmad=.false.
      lfgroupbystruct=.false.
      lqgroupbystruct=.false.
      cutoff_type=1
      lmixpoints=.false.
      lscreen=.false.
      rscreen_cut=0.0d0
      rscreen_onset=0.0d0
      lsilent=.false.
      lpreparemd=.false.
      lseparatebiasini=.false.
      nn_type_vdw=0
      lpearson_correlation=.false.
      lweightanalysis=.false.
      lenvironmentanalysis=.false.
      lfindcontradictions=.false.
      lmd=.false.
      ldynforcegroup=.false.
      dynforcegroup_start=20
      dynforcegroup_step=2
      lshuffle_weights_short_atomic=.false.
      nshuffle_weights_short_atomic=10
      shuffle_weights_short_atomic=0.1d0
      ldetect_saturation=.false.
      saturation_threshold=0.99d0
      dataclusteringthreshold1=1.0d0
      dataclusteringthreshold2=0.0d0
      ldataclustering=.false.
      lprintdateandtime=.false.
      lenableontheflyinput=.false.
      lcheckinputforces=.false.
      lionforcesonly=.false.
      inputforcethreshold=0.001d0
      luseedkalman=.false.
      ledforcesv2=.false.
      luseipi=.false.
      ipisocket=" "
      ipistring=" "
      desc_str=" "
!!
      return
      end 
