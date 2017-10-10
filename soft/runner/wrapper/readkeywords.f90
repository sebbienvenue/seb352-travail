!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine readkeywords(iseed,&
        nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,nodes_ham_temp,&
        kalmanlambda_local,kalmanlambdae_local)
!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use mode1options
      use predictionoptions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use inputnncounters
!!    
      implicit none
!!
      integer nodes_short_atomic_temp(0:maxnum_layers_short_atomic)  ! out 
      integer nodes_elec_temp(0:maxnum_layers_elec)                  ! out 
      integer nodes_short_pair_temp(0:maxnum_layers_short_pair)      ! out 
      integer nodes_ham_temp(0:maxnum_layers_ham)                    ! out 
      integer iseed                                     ! out, seed for weight initialization
      integer i
      integer i1,i2                                     ! internal
      integer mode_local                                ! internal
      integer nn_type_short_local                       ! internal
      integer nn_type_elec_local                        ! internal
      integer nn_type_nntb_local                        ! internal
!!
      real*8 kalmanlambda_local                         ! internal
      real*8 kalmanlambdae_local                        ! internal
!!
      character*40 dummy                                ! internal
      character*40 dummy2                               ! internal
      character*40 keyword                              ! internal
!!
      logical lshort_local                              ! internal
      logical lelec_local                               ! internal
      logical lnntb_local                               ! internal
!!
!! local initializations
      lshort_local=.false.
      lelec_local =.false.
      lnntb_local =.false.
      nn_type_short_local=0
      nn_type_elec_local=0
      nn_type_nntb_local=0
      mode_local  =0
!!
      open(nnunit,file='input.nn',form='formatted',status='old')
      rewind(nnunit)
 10   continue
!!
      read(nnunit,*,END=20) keyword
!!
!! check for comment lines
      if(keyword(1:1).eq.'#')then
        goto 10
      endif
!! check for blank lines
      if(keyword.eq.'')then
        goto 10
      endif
!!
!! analyze keywords
!!
      if(keyword.eq.'runner_mode')then
        if(count_mode.gt.0)then
          write(ounit,*)'Error: runner_mode specified twice'
          stop
        endif
        count_mode=count_mode+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,mode_local
        if(mode.ne.mode_local)then
          write(ounit,*)'ERROR in determination of mode'
          stop
        endif
!!
      elseif((keyword.eq.'electrostatic_type').or.(keyword.eq.'nn_type_elec'))then
        if(count_nn_type_elec.gt.0)then
          write(ounit,*)'Error: electrostatic_type specified twice'
          stop
        endif
        count_nn_type_elec=count_nn_type_elec+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nn_type_elec_local
        if(nn_type_elec.ne.nn_type_elec_local)then
          write(ounit,*)'ERROR in determination of nn_type_elec'
          stop
        endif
!!
      elseif((keyword.eq.'vdw_type').or.(keyword.eq.'nn_type_vdw'))then
        if(count_vdw_type.gt.0)then
          write(ounit,*)'Error: vdw_type specified twice'
          stop
        endif
        count_vdw_type=count_vdw_type+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nn_type_vdw
!!
      elseif(keyword.eq.'cutoff_vdw')then
        if(count_cutoffvdw.gt.0)then
          write(ounit,*)'Error: cutoff_vdw specified twice'
          stop
        endif
        count_cutoffvdw=count_cutoffvdw+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,cutoffvdw
!!
      elseif(keyword.eq.'vdw_param')then
!!      let this just pass here, values are read in subroutine getvdwparams
!!
      elseif(keyword.eq.'use_electrostatic_nn')then
        write(ounit,*)'keyword use_electrostatic_nn is obsolete, use electrostatic_type and use_electrostatics now'
        stop
!!        
      elseif(keyword.eq.'debug_mode')then
        if(count_ldebug.gt.0)then
          write(ounit,*)'Error: ldebug specified twice'
          stop
        endif
        count_ldebug=count_ldebug+1
        ldebug=.true.
!!
      elseif(keyword.eq.'cutoff_type')then
        if(count_cutoff_type.gt.0)then
          write(ounit,*)'Error: cutoff_type specified twice'
          stop
        endif
        count_cutoff_type=count_cutoff_type+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,cutoff_type
!!
      elseif(keyword.eq.'dynamic_force_grouping')then
        if(count_ldynforcegroup.gt.0)then
          write(ounit,*)'Error: dynamic_force_grouping specified twice'
          stop
        endif
        count_ldynforcegroup=count_ldynforcegroup+1
        ldynforcegroup=.true.
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,dynforcegroup_start,dynforcegroup_step
!!
      elseif(keyword.eq.'analyze_error_energy_step')then
        if(count_analyze_error_energy_step.gt.0)then
          write(ounit,*)'Error: analyze_error_energy_step specified twice'
          stop
        endif
        count_analyze_error_energy_step=count_analyze_error_energy_step+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,analyze_error_energy_step
!!
      elseif(keyword.eq.'analyze_error_force_step')then
        if(count_analyze_error_force_step.gt.0)then
          write(ounit,*)'Error: analyze_error_force_step specified twice'
          stop
        endif
        count_analyze_error_force_step=count_analyze_error_force_step+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,analyze_error_force_step
!!
      elseif(keyword.eq.'analyze_error_charge_step')then
        if(count_analyze_error_charge_step.gt.0)then
          write(ounit,*)'Error: analyze_error_charge_step specified twice'
          stop
        endif
        count_analyze_error_charge_step=count_analyze_error_charge_step+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,analyze_error_charge_step
!!
      elseif(keyword.eq.'parallel_mode')then
        if(count_paramode.gt.0)then
          write(ounit,*)'Error: parallel_mode specified twice'
          stop
        endif
        count_paramode=count_paramode+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,paramode
!!
      elseif(keyword.eq.'use_short_nn')then
        if(count_lshort.gt.0)then
          write(ounit,*)'Error: use_short_nn specified twice'
          stop
        endif
        count_lshort=count_lshort+1
        lshort_local=.true.
        if(lshort.neqv.lshort_local)then
          write(ounit,*)'ERROR in determination of lshort'
          stop
        endif
!!
      elseif(keyword.eq.'symfunction_correlation')then
        if(count_lpearson_correlation.gt.0)then
          write(ounit,*)'Error: symfunction_correlation specified twice'
          stop
        endif
        count_lpearson_correlation=count_lpearson_correlation+1
        lpearson_correlation=.true.
!!
      elseif(keyword.eq.'weight_analysis')then
        if(count_lweightanalysis.gt.0)then
          write(ounit,*)'Error: weight_analysis specified twice'
          stop
        endif
        count_lweightanalysis=count_lweightanalysis+1
        lweightanalysis=.true.
!!
      elseif(keyword.eq.'environment_analysis')then
        if(count_lenvironmentanalysis.gt.0)then
          write(ounit,*)'Error: environment_analysis specified twice'
          stop
        endif
        count_lenvironmentanalysis=count_lenvironmentanalysis+1
        lenvironmentanalysis=.true.
!!
      elseif(keyword.eq.'find_contradictions')then
        if(count_lfindcontradictions.gt.0)then
          write(ounit,*)'Error: find_contradictions specified twice'
          stop
        endif
        count_lfindcontradictions=count_lfindcontradictions+1
        lfindcontradictions=.true.
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,deltagthres,deltafthres
!!
      elseif(keyword.eq.'use_electrostatics')then
        if(count_lelec.gt.0)then
          write(ounit,*)'Error: use_electrostatics specified twice'
          stop
        endif
        count_lelec=count_lelec+1
        lelec_local=.true.
        if(lelec.neqv.lelec_local)then
          write(ounit,*)'ERROR in determination of lelec'
          stop
        endif
!!
      elseif(keyword.eq.'use_hamiltonian')then
        if(count_lnntb.gt.0)then
          write(ounit,*)'Error: use_hamiltonian specified twice'
          stop
        endif
        count_lnntb=count_lnntb+1
        lnntb_local=.true.
        if(lnntb.neqv.lnntb_local)then
          write(ounit,*)'ERROR in determination of lnntb'
          stop
        endif
!!
      elseif(keyword.eq.'use_old_scaling')then
        if(count_luseoldscaling.gt.0)then
          write(ounit,*)'Error: use_old_scaling specified twice'
          stop
        endif
        count_luseoldscaling=count_luseoldscaling+1
        luseoldscaling=.true.
!!
      elseif(keyword.eq.'md_mode')then
        if(count_lmd.gt.0)then
          write(ounit,*)'Error: md_mode specified twice'
          stop
        endif
        count_lmd=count_lmd+1
        lmd=.true.
!!
      elseif(keyword.eq.'global_hidden_layers_short')then
        if(lshort.and.(nn_type_short.eq.1))then
          if(count_num_layers_short_atomic.gt.0)then
            write(ounit,*)'Error: global_hidden_layers_short specified twice'
            stop !'
          endif
          count_num_layers_short_atomic=count_num_layers_short_atomic+1
        endif
!!
      elseif(keyword.eq.'global_hidden_layers_electrostatic')then
        if(count_num_layers_elec.gt.0)then
          write(ounit,*)'Error: global_hidden_layers_electrostatic specified twice'
          stop !'
        endif
        count_num_layers_elec=count_num_layers_elec+1
!!      maxnum_layersewald is read before in getdimensions
!!'
      elseif(keyword.eq.'global_hidden_layers_pair')then
        if(count_num_layers_short_pair.gt.0)then
          write(ounit,*)'Error: global_hidden_layers_pair specified twice'
          stop !'
        endif
        count_num_layers_short_pair=count_num_layers_short_pair+1
!!      maxnum_layerspair is read before in getdimensions
!!
      elseif(keyword.eq.'global_hidden_layers_ham')then
        if(count_num_layers_ham.gt.0)then
          write(ounit,*)'Error: global_hidden_layers_ham specified twice'
          stop !'
        endif
        count_num_layers_ham=count_num_layers_ham+1
!!      maxnum_layersham is read before in getdimensions
!!
      elseif((keyword.eq.'global_nodes_short')&
         .or.(keyword.eq.'global_nodes_short_atomic'))then
        if(lshort.and.(nn_type_short.eq.1))then
          if(count_nodes_short_atomic.gt.0)then
            write(ounit,*)'Error: global_nodes_short specified twice'
            stop
          endif
          count_nodes_short_atomic=count_nodes_short_atomic+1
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,(nodes_short_atomic_temp(i1),i1=1,maxnum_layers_short_atomic-1)
          do i1=1,nelem
            do i2=1,maxnum_layers_short_atomic
              nodes_short_atomic(i2,i1)=nodes_short_atomic_temp(i2)
            enddo
          enddo
        endif
!!
      elseif(keyword.eq.'global_nodes_electrostatic')then
        if(lelec.and.(nn_type_elec.eq.1))then
          if(count_nodes_elec.gt.0)then
            write(ounit,*)'Error: global_nodes_electrostatic specified twice'
            stop !'
          endif
          count_nodes_elec=count_nodes_elec+1
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,(nodes_elec_temp(i1),i1=1,maxnum_layers_elec-1)
          do i1=1,nelem
            do i2=1,maxnum_layers_elec
              nodes_elec(i2,i1)=nodes_elec_temp(i2)
            enddo
          enddo
        endif
!!
      elseif((keyword.eq.'global_nodes_pair')&
         .or.(keyword.eq.'global_nodes_short_pair'))then
        if(lshort.and.(nn_type_short.eq.2))then
          if(count_nodes_short_pair.gt.0)then
            write(ounit,*)'Error: global_nodes_pair specified twice'
            stop
          endif
          count_nodes_short_pair=count_nodes_short_pair+1
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,(nodes_short_pair_temp(i1),i1=1,maxnum_layers_short_pair-1)
          do i1=1,npairs
            do i2=1,maxnum_layers_short_pair
              nodes_short_pair(i2,i1)=nodes_short_pair_temp(i2)
            enddo
          enddo
        endif
!!
      elseif(keyword.eq.'global_nodes_ham')then
        if(lnntb)then
          if(count_nodes_ham.gt.0)then
            write(ounit,*)'Error: global_nodes_ham specified twice'
            stop
          endif
          count_nodes_ham=count_nodes_ham+1
          backspace(nnunit)
          read(nnunit,*,ERR=99)dummy,(nodes_ham_temp(i1),i1=1,maxnum_layers_ham-1)
          do i1=1,npairs
            do i2=1,maxnum_layers_ham
              nodes_ham(i2,i1)=nodes_ham_temp(i2)
            enddo
          enddo
        endif
!!
      elseif(keyword.eq.'global_output_nodes_short')then
        write(ounit,*)'Error: global_output_nodes_short keyword is obsolete, please remove it'
        stop
!!
      elseif(keyword.eq.'global_output_nodes_electrostatic')then
        write(ounit,*)'Error: global_output_nodes_electrostatic keyword is obsolete, please remove'
        stop
!!
      elseif(keyword.eq.'global_output_nodes_pair')then
        write(ounit,*)'Error: global_output_nodes_pair keyword is obsolete, please remove'
        stop !'
!!
      elseif(keyword.eq.'global_output_nodes_ham')then
        write(ounit,*)'Error: global_output_nodes_ham keyword is obsolete, please remove'
        stop !'
!!
      elseif(keyword.eq.'ewald_alpha')then
        if(count_ewaldalpha.gt.0)then
          write(ounit,*)'Error: ewald_alpha specified twice'
          stop
        endif
        count_ewaldalpha=count_ewaldalpha+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,ewaldalpha
!!
      elseif(keyword.eq.'ewald_cutoff')then
        if(count_ewaldcutoff.gt.0)then
          write(ounit,*)'Error: ewald_cutoff specified twice'
          stop
        endif
        count_ewaldcutoff=count_ewaldcutoff+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,ewaldcutoff
!!
      elseif(keyword.eq.'ewald_kmax')then
        if(count_ewaldkmax.gt.0)then
          write(ounit,*)'Error: ewald_kmax specified twice'
          stop
        endif
        count_ewaldkmax=count_ewaldkmax+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,ewaldkmax
!!
      elseif(keyword.eq.'precondition_weights')then
        if(count_lprecond.gt.0)then
          write(ounit,*)'Error: precondition_weights specified twice'
          stop
        endif
        count_lprecond=count_lprecond+1
        lprecond=.true.
!!
      elseif(keyword.eq.'initialization_only')then
        if(count_linionly.gt.0)then
          write(ounit,*)'Error: initialization_only specified twice'
          stop
        endif
        count_linionly=count_linionly+1
        linionly=.true.
!!
      elseif(keyword.eq.'force_grouping_by_structure')then
        if(count_lfgroupbystruct.gt.0)then
          write(ounit,*)'Error: force_grouping_by_structure specified twice'
          stop
        endif
        count_lfgroupbystruct=count_lfgroupbystruct+1
        lfgroupbystruct=.true.
!!
      elseif(keyword.eq.'charge_grouping_by_structure')then
        if(count_lqgroupbystruct.gt.0)then
          write(ounit,*)'Error: charge_grouping_by_structure specified twice'
          stop
        endif
        count_lqgroupbystruct=count_lqgroupbystruct+1
        lqgroupbystruct=.true.
!!
      elseif(keyword.eq.'mix_all_points')then
        if(count_lmixpoints.gt.0)then
          write(ounit,*)'Error: mix_all_points specified twice'
          stop
        endif
        count_lmixpoints=count_lmixpoints+1
        lmixpoints=.true.
!!
      elseif(keyword.eq.'print_convergence_vector')then
        if(count_lprintconv.gt.0)then
          write(ounit,*)'Error: print_convergence_vector specified twice'
          stop
        endif
        count_lprintconv=count_lprintconv+1
        lprintconv=.true.
!!
      elseif(keyword.eq.'print_mad')then
        if(count_lprintmad.gt.0)then
          write(ounit,*)'Error: print_mad specified twice'
          stop
        endif
        count_lprintmad=count_lprintmad+1
        lprintmad=.true.
!!
      elseif(keyword.eq.'noise_energy')then
        if(count_noisee.gt.0)then
          write(ounit,*)'Error: noise_energy specified twice'
          stop
        endif
        count_noisee=count_noisee+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,noisee
!!
      elseif(keyword.eq.'noise_force')then
        if(count_noisef.gt.0)then
          write(ounit,*)'Error: noise_force specified twice'
          stop
        endif
        count_noisef=count_noisef+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,noisef
!!
      elseif(keyword.eq.'noise_charge')then
        if(count_noiseq.gt.0)then
          write(ounit,*)'Error: noise_charge specified twice'
          stop
        endif
        count_noiseq=count_noiseq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,noiseq
!!
      elseif(keyword.eq.'short_energy_group')then
        if(count_short_energy_group.gt.0)then
          write(ounit,*)'Error: short_energy_group specified twice'
          stop
        endif
        count_short_energy_group=count_short_energy_group+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nenergygroup
!!
      elseif(keyword.eq.'short_force_group')then
        if(count_short_force_group.gt.0)then
          write(ounit,*)'Error: short_force_group specified twice'
          stop
        endif
        count_short_force_group=count_short_force_group+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nforcegroup
!!
      elseif(keyword.eq.'charge_group')then
        if(count_charge_group.gt.0)then
          write(ounit,*)'Error: charge_group specified twice'
          stop
        endif
        count_charge_group=count_charge_group+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nchargegroup
!!
      elseif(keyword.eq.'use_short_forces')then
        if(count_luseforces.gt.0)then
          write(ounit,*)'Error: use_short_forces specified twice'
          stop
        endif
        count_luseforces=count_luseforces+1
        luseforces=.true.
!!
      elseif(keyword.eq.'short_energy_fraction')then
        if(count_energyrnd.gt.0)then
          write(ounit,*)'Error: short_energy_fraction specified twice'
          stop
        endif
        count_energyrnd=count_energyrnd+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,energyrnd
!!
      elseif(keyword.eq.'short_force_fraction')then
        if(count_forcernd.gt.0)then
          write(ounit,*)'Error: short_force_fraction specified twice'
          stop
        endif
        count_forcernd=count_forcernd+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,forcernd
!!
      elseif(keyword.eq.'charge_fraction')then
        if(count_chargernd.gt.0)then
          write(ounit,*)'Error: charge_fraction specified twice'
          stop
        endif
        count_chargernd=count_chargernd+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,chargernd
!!
      elseif(keyword.eq.'use_atom_charges')then
        if(count_luseatomcharges.gt.0)then
          write(ounit,*)'Error: use_atom_charges specified twice'
          stop
        endif
        count_luseatomcharges=count_luseatomcharges+1
        luseatomcharges=.true.
!!
      elseif(keyword.eq.'use_atom_energies')then
        if(count_luseatomenergies.gt.0)then
          write(ounit,*)'Error: use_atom_energies specified twice'
          stop
        endif
        count_luseatomenergies=count_luseatomenergies+1
        luseatomenergies=.true.
!!
      elseif(keyword.eq.'remove_atom_energies')then
        if(count_lremoveatomenergies.gt.0)then
          write(ounit,*)'Error: remove_atom_energies specified twice'
          stop
        endif
        count_lremoveatomenergies=count_lremoveatomenergies+1
        lremoveatomenergies=.true.
!!
      elseif(keyword.eq.'analyze_error')then
        if(count_lanalyzeerror.gt.0)then
          write(ounit,*)'Error: analyze_error specified twice'
          stop
        endif
        count_lanalyzeerror=count_lanalyzeerror+1
        lanalyzeerror=.true.
!!
      elseif(keyword.eq.'use_charge_constraint')then
        if(count_lchargeconstraint.gt.0)then
          write(ounit,*)'Error: use_charge_constraint specified twice'
          stop
        endif
        count_lchargeconstraint=count_lchargeconstraint+1
        lchargeconstraint=.true.
!!
      elseif(keyword.eq.'fitmode')then
        if(count_fitmode.gt.0)then
          write(ounit,*)'Error: mode specified twice'
          stop
        endif
        count_fitmode=count_fitmode+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,fitmode
!!
      elseif(keyword.eq.'energy_threshold')then
        if(count_fitethres.gt.0)then
          write(ounit,*)'Error: energy_threshold specified twice'
          stop
        endif
        count_fitethres=count_fitethres+1
        backspace(nnunit)
        lfitethres=.true.
        read(nnunit,*,ERR=99)dummy,fitethres
!!
      elseif(keyword.eq.'bond_threshold')then
        if(count_rmin.gt.0)then
          write(ounit,*)'Error: bond_threshold specified twice'
          stop
        endif
        count_rmin=count_rmin+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,rmin
!!
      elseif(keyword.eq.'optmode_short_energy')then
        if(count_optmodee.gt.0)then
          write(ounit,*)'Error: optmode_short_energy specified twice'
          stop
        endif
        count_optmodee=count_optmodee+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,optmodee
!!
      elseif(keyword.eq.'optmode_short_force')then
        if(count_optmodef.gt.0)then
          write(ounit,*)'Error: optmode_short_force specified twice'
          stop
        endif
        count_optmodef=count_optmodef+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,optmodef
!!
      elseif(keyword.eq.'optmode_charge')then
        if(count_optmodeq.gt.0)then
          write(ounit,*)'Error: optmode_charge specified twice'
          stop
        endif
        count_optmodeq=count_optmodeq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,optmodeq
!!
      elseif(keyword.eq.'random_seed')then
        if(count_iseed.gt.0)then
          write(ounit,*)'Error: random_seed specified twice'
          stop
        endif
        count_iseed=count_iseed+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,iseed
!!
      elseif((keyword.eq.'points_in_memory').or.(keyword.eq.'nblock'))then
        if(count_nblock.gt.0)then
          write(ounit,*)'Error: points_in_memory specified twice'
          stop
        endif
        count_nblock=count_nblock+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nblock
!!
      elseif(keyword.eq.'epochs')then
        if(count_epochs.gt.0)then
          write(ounit,*)'Error: epochs specified twice'
          stop
        endif
        count_epochs=count_epochs+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nepochs
!!
      elseif(keyword.eq.'number_of_elements')then
        if(count_nelem.gt.0)then
          write(ounit,*)'Error: number_of_elements specified twice'
          stop
        endif
        count_nelem=count_nelem+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nelem
!!
      elseif(keyword.eq.'elements')then
        if(count_element.gt.0)then
          write(ounit,*)'Error: elements specified twice'
          stop
        endif
        count_element=count_element+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,(element(i),i=1,nelem)
!!
      elseif(keyword.eq.'write_weights_epoch')then
        if(count_iwriteweight.gt.0)then
          write(ounit,*)'Error: write_weights_epoch specified twice'
          stop
        endif
        count_iwriteweight=count_iwriteweight+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,iwriteweight
!!
      elseif(keyword.eq.'write_temporary_weights')then
        if(count_lwritetmpweights.gt.0)then
          write(ounit,*)'Error: write_temporary_weights specified twice'
          stop
        endif
        count_lwritetmpweights=count_lwritetmpweights+1
        lwritetmpweights=.true.
!!
      elseif(keyword.eq.'test_fraction')then
        if(count_splitthres.gt.0)then
          write(ounit,*)'Error: test_fraction specified twice'
          stop
        endif
        count_splitthres=count_splitthres+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,splitthres
!!
      elseif(keyword.eq.'scale_min_short_atomic')then
        if(count_scmin_short_atomic.gt.0)then
          write(ounit,*)'Error: scale_min_short/scale_min_short_atomic specified twice'
          stop
        endif
        count_scmin_short_atomic=count_scmin_short_atomic+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmin_short_atomic
!!
      elseif(keyword.eq.'scale_max_short_atomic')then
        if(count_scmax_short_atomic.gt.0)then
          write(ounit,*)'Error: scale_max_short/scale_max_short_atomic specified twice'
          stop
        endif
        count_scmax_short_atomic=count_scmax_short_atomic+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmax_short_atomic
!!
      elseif(keyword.eq.'scale_min_short_pair')then
        if(count_scmin_short_pair.gt.0)then
          write(ounit,*)'Error: scale_min_pair/scale_min_short_pair specified twice'
          stop
        endif
        count_scmin_short_pair=count_scmin_short_pair+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmin_short_pair
!!
      elseif(keyword.eq.'scale_max_short_pair')then
        if(count_scmax_short_pair.gt.0)then
          write(ounit,*)'Error: scale_max_pair/scale_max_short_pair specified twice'
          stop
        endif
        count_scmax_short_pair=count_scmax_short_pair+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmax_short_pair
!!
      elseif(keyword.eq.'scale_min_ham')then
        if(count_scmin_ham.gt.0)then
          write(ounit,*)'Error: scale_min_ham specified twice'
          stop
        endif
        count_scmin_ham=count_scmin_ham+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmin_ham
!!
      elseif(keyword.eq.'scale_max_ham')then
        if(count_scmax_ham.gt.0)then
          write(ounit,*)'Error: scale_max_ham specified twice'
          stop
        endif
        count_scmax_ham=count_scmax_ham+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmax_ham
!!
      elseif(keyword.eq.'scale_min_elec')then
        if(count_scmin_elec.gt.0)then
          write(ounit,*)'Error: scale_min_ewald/scale_min_elec specified twice'
          stop
        endif
        count_scmin_elec=count_scmin_elec+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmin_elec
!!
      elseif(keyword.eq.'scale_max_elec')then
        if(count_scmax_elec.gt.0)then
          write(ounit,*)'Error: scale_max_ewald/scale_max_elec specified twice'
          stop
        endif
        count_scmax_elec=count_scmax_elec+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scmax_elec
!!
      elseif(keyword.eq.'short_energy_error_threshold')then
        if(count_kalmanthreshold.gt.0)then
          write(ounit,*)'Error: short_energy_error_threshold specified twice'
          stop
        endif
        count_kalmanthreshold=count_kalmanthreshold+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanthreshold
!!
      elseif(keyword.eq.'short_force_error_threshold')then
        if(count_kalmanthresholdf.gt.0)then
          write(ounit,*)'Error: short_force_error_threshold specified twice'
          stop
        endif
        count_kalmanthresholdf=count_kalmanthresholdf+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanthresholdf
!!
      elseif(keyword.eq.'charge_error_threshold')then
        if(count_kalmanthresholde.gt.0)then
          write(ounit,*)'Error: charge_error_threshold specified twice'
          stop
        endif
        count_kalmanthresholde=count_kalmanthresholde+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanthresholde
!!
      elseif(keyword.eq.'total_charge_error_threshold')then
        if(count_kalmanthresholdc.gt.0)then
          write(ounit,*)'Error: total_charge_error_threshold specified twice'
          stop !'
        endif
        count_kalmanthresholdc=count_kalmanthresholdc+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanthresholdc
!!
      elseif(keyword.eq.'kalman_damp_short')then
        if(count_kalman_dampe.gt.0)then
          write(ounit,*)'Error: kalman_damp_short specified twice'
          stop
        endif
        count_kalman_dampe=count_kalman_dampe+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalman_dampe
!!
      elseif(keyword.eq.'kalman_damp_force')then
        if(count_kalman_dampf.gt.0)then
          write(ounit,*)'Error: kalman_damp_force specified twice'
          stop
        endif
        count_kalman_dampf=count_kalman_dampf+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalman_dampf
!!
      elseif(keyword.eq.'kalman_damp_charge')then
        if(count_kalman_dampq.gt.0)then
          write(ounit,*)'Error: kalman_damp_charge specified twice'
          stop
        endif
        count_kalman_dampq=count_kalman_dampq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalman_dampq
!!
      elseif(keyword.eq.'kalman_lambda_short')then
        if(count_kalmanlambda.gt.0)then
          write(ounit,*)'Error: kalman_lambda_short specified twice'
          stop
        endif
        count_kalmanlambda=count_kalmanlambda+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanlambda_local
!!
      elseif(keyword.eq.'kalman_lambda_charge')then
        if(count_kalmanlambdae.gt.0)then
          write(ounit,*)'Error: kalman_lambda_charge specified twice'
          stop
        endif
        count_kalmanlambdae=count_kalmanlambdae+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanlambdae_local
!!
      elseif(keyword.eq.'kalman_lambda_charge_constraint')then
        if(count_kalmanlambdac.gt.0)then
          write(ounit,*)'Error: kalman_lambda_charge_constraint specified twice'
          stop
        endif
        count_kalmanlambdac=count_kalmanlambdac+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmanlambdac
!!
      elseif(keyword.eq.'kalman_nue_short')then
        if(count_kalmannue.gt.0)then
          write(ounit,*)'Error: kalman_nue_short specified twice'
          stop
        endif
        count_kalmannue=count_kalmannue+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmannue
!!
      elseif(keyword.eq.'kalman_nue_charge')then
        if(count_kalmannuee.gt.0)then
          write(ounit,*)'Error: kalman_nue_charge specified twice'
          stop
        endif
        count_kalmannuee=count_kalmannuee+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmannuee
!!
      elseif(keyword.eq.'kalman_nue_charge_constraint')then
        if(count_kalmannuec.gt.0)then
          write(ounit,*)'Error: kalman_nue_charge_constraint specified twice'
          stop
        endif
        count_kalmannuec=count_kalmannuec+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,kalmannuec
!!
      elseif(keyword.eq.'steepest_descent_step_energy_short')then
        if(count_steepeststepe.gt.0)then
          write(ounit,*)'Error: steepest_descent_step_energy_short specified twice'
          stop
        endif
        count_steepeststepe=count_steepeststepe+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,steepeststepe
!!
      elseif(keyword.eq.'steepest_descent_step_force_short')then
        if(count_steepeststepf.gt.0)then
          write(ounit,*)'Error: steepest_descent_step_force_short specified twice'
          stop
        endif
        count_steepeststepf=count_steepeststepf+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,steepeststepf
!!
      elseif(keyword.eq.'steepest_descent_step_charge')then
        if(count_steepeststepq.gt.0)then
          write(ounit,*)'Error: steepest_descent_step_charge specified twice'
          stop
        endif
        count_steepeststepq=count_steepeststepq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,steepeststepq
!!
      elseif(keyword.eq.'force_update_scaling')then
        if(count_scalefactorf.gt.0)then
          write(ounit,*)'Error: force_update_scaling specified twice'
          stop
        endif
        count_scalefactorf=count_scalefactorf+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scalefactorf
!!
      elseif(keyword.eq.'charge_update_scaling')then
        if(count_scalefactorq.gt.0)then
          write(ounit,*)'Error: charge_update_scaling specified twice'
          stop
        endif
        count_scalefactorq=count_scalefactorq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,scalefactorq
!!
      elseif(keyword.eq.'random_order_training')then
        if(count_lrandomtrain.gt.0)then
          write(ounit,*)'Error: random_order_training specified twice'
          stop
        endif
        count_lrandomtrain=count_lrandomtrain+1
        lrandomtrain=.true.
!!
      elseif(keyword.eq.'scale_symmetry_functions')then
        if(count_lscalesym.gt.0)then
          write(ounit,*)'Error: scale_symmetry_functions specified twice'
          stop
        endif
        count_lscalesym=count_lscalesym+1
        lscalesym=.true.
!!
      elseif(keyword.eq.'center_symmetry_functions')then
        if(count_lcentersym.gt.0)then
          write(ounit,*)'Error: center_symmetry_functions specified twice'
          stop
        endif
        count_lcentersym=count_lcentersym+1
        lcentersym=.true.
!!
      elseif(keyword.eq.'use_old_weights_short')then
        if(count_luseoldweightsshort.gt.0)then
          write(ounit,*)'Error: use_old_weights_short specified twice'
          stop
        endif
        count_luseoldweightsshort=count_luseoldweightsshort+1
        luseoldweightsshort=.true.
!!
      elseif(keyword.eq.'use_old_weights_charge')then
        if(count_luseoldweightscharge.gt.0)then
          write(ounit,*)'Error: use_old_weights_charge specified twice'
          stop
        endif
        count_luseoldweightscharge=count_luseoldweightscharge+1
        luseoldweightscharge=.true.
!!
      elseif(keyword.eq.'save_kalman_matrices')then
        if(count_lsavekalman.gt.0)then
          write(ounit,*)'Error: save_kalman_matrices specified twice'
          stop
        endif
        count_lsavekalman=count_lsavekalman+1
        lsavekalman=.true.
!!
      elseif(keyword.eq.'read_kalman_matrices')then
        if(count_lrestkalman.gt.0)then
          write(ounit,*)'Error: read_kalman_matrices specified twice'
          stop
        endif
        count_lrestkalman=count_lrestkalman+1
        lrestkalman=.true.
!!
      elseif(keyword.eq.'update_single_element')then
        if(count_elemupdate.gt.0)then
          write(ounit,*)'Error: update_single_element specified twice'
          stop
        endif
        count_elemupdate=count_elemupdate+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,elemupdate
        lupdatebyelement=.true.
!!
      elseif(keyword.eq.'update_worst_short_energies')then
        if(count_luseworste.gt.0)then
          write(ounit,*)'Error: update_worst_short_energies specified twice'
          stop
        endif
        count_luseworste=count_luseworste+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,worste
        luseworste=.true.
!!
      elseif(keyword.eq.'update_worst_short_forces')then
        if(count_luseworstf.gt.0)then
          write(ounit,*)'Error: update_worst_short_forces specified twice'
          stop
        endif
        count_luseworstf=count_luseworstf+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,worstf
        luseworstf=.true.
!!
      elseif(keyword.eq.'update_worst_charges')then
        if(count_luseworstq.gt.0)then
          write(ounit,*)'Error: update_worst_charges specified twice'
          stop
        endif
        count_luseworstq=count_luseworstq+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,worstq
        luseworstq=.true.
!!
      elseif(keyword.eq.'growth_mode')then
        if(count_growth.gt.0)then
          write(ounit,*)'Error: growth_mode specified twice'
          stop
        endif
        count_growth=count_growth+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,ngrowth,growthstep
        lgrowth=.true.
!!
      elseif(keyword.eq.'use_damping')then
        if(count_dampw.gt.0)then
          write(ounit,*)'Error: use_damping specified twice'
          stop
        endif
        count_dampw=count_dampw+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,dampw
        ldampw=.true.
!!
      elseif(keyword.eq.'fix_weights')then
        if(count_lfixweights.gt.0)then
          write(ounit,*)'Error: fix_weights specified twice'
          stop
        endif
        count_lfixweights=count_lfixweights+1
        lfixweights=.true.
!!
      elseif(keyword.eq.'calculate_forces')then
        if(count_ldoforces.gt.0)then
          write(ounit,*)'Error: calculate_forces specified twice'
          stop
        endif
        count_ldoforces=count_ldoforces+1
        ldoforces=.true.
!!
      elseif(keyword.eq.'calculate_stress')then
        if(count_ldostress.gt.0)then
          write(ounit,*)'Error: calculate_stress specified twice'
          stop
        endif
        count_ldostress=count_ldostress+1
        ldostress=.true.
!!
      elseif(keyword.eq.'enforce_max_num_neighbors_atomic')then
        if(count_lenforcemaxnumneighborsatomic.gt.0)then
          write(ounit,*)'Error: enforce_max_num_neighbors_atomic specified twice'
          stop
        endif
        count_lenforcemaxnumneighborsatomic=count_lenforcemaxnumneighborsatomic+1
        lenforcemaxnumneighborsatomic=.true.
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,max_num_neighbors_atomic_input
!!
      elseif(keyword.eq.'detailed_timing')then
        if(count_lfinetime.gt.0)then
          write(ounit,*)'Error: detailed_timing specified twice'
          stop
        endif
        count_lfinetime=count_lfinetime+1
        lfinetime=.true.
!!
      elseif(keyword.eq.'detailed_timing_epoch')then
        if(count_lfinetimeepoch.gt.0)then
          write(ounit,*)'Error: detailed_timing_epoch specified twice'
          stop
        endif
        count_lfinetimeepoch=count_lfinetimeepoch+1
        lfinetimeepoch=.true.
!!
      elseif(keyword.eq.'write_pdb')then
        if(count_lwritepdb.gt.0)then
          write(ounit,*)'Error: write_pdb specified twice'
          stop
        endif
        count_lwritepdb=count_lwritepdb+1
        lwritepdb=.true.
!!
      elseif(keyword.eq.'write_xyz')then
        if(count_lwritexyz.gt.0)then
          write(ounit,*)'Error: write_xyz specified twice'
          stop
        endif
        count_lwritexyz=count_lwritexyz+1
        lwritexyz=.true.
!!
      elseif(keyword.eq.'write_pov')then
        if(count_lwritepov.gt.0)then
          write(ounit,*)'Error: write_pov specified twice'
          stop
        endif
        count_lwritepov=count_lwritepov+1
        lwritepov=.true.
!!
      elseif(keyword.eq.'write_pwscf')then
        if(count_lwritepw.gt.0)then
          write(ounit,*)'Error: write_pwscf specified twice'
          stop
        endif
        count_lwritepw=count_lwritepw+1
        lwritepw=.true.
!!
      elseif(keyword.eq.'write_trainpoints')then
        if(count_lwritetrainpoints.gt.0)then
          write(ounit,*)'Error: write_trainpoints specified twice'
          stop
        endif
        count_lwritetrainpoints=count_lwritetrainpoints+1
        lwritetrainpoints=.true.
!!
      elseif(keyword.eq.'write_trainforces')then
        if(count_lwritetrainforces.gt.0)then
          write(ounit,*)'Error: write_trainforces specified twice'
          stop
        endif
        count_lwritetrainforces=count_lwritetrainforces+1
        lwritetrainforces=.true.
!!
      elseif(keyword.eq.'write_traincharges')then
        if(count_lwritetraincharges.gt.0)then
          write(ounit,*)'Error: write_traincharges specified twice'
          stop
        endif
        count_lwritetraincharges=count_lwritetraincharges+1
        lwritetraincharges=.true.
!!
      elseif(keyword.eq.'max_force')then
        if(count_maxforce.gt.0)then
          write(ounit,*)'Error: max_force specified twice'
          stop
        endif
        count_maxforce=count_maxforce+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,maxforce
!!
      elseif(keyword.eq.'max_energy')then
        if(count_maxenergy.gt.0)then
          write(ounit,*)'Error: max_energy specified twice'
          stop
        endif
        count_maxenergy=count_maxenergy+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,maxenergy
!!
      elseif(keyword.eq.'nn_type')then
        write(ounit,*)'Error: nn_type keyword is outdated, use nn_type_short instead'
        stop
!!
      elseif(keyword.eq.'nn_type_short')then
        if(count_nn_type_short.gt.0)then
          write(ounit,*)'Error: nn_type_short specified twice'
          stop
        endif
        count_nn_type_short=count_nn_type_short+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nn_type_short_local
        if(nn_type_short.ne.nn_type_short_local)then
          write(ounit,*)'ERROR in determination of nn_type_short'
          stop
        endif
!!
      elseif(keyword.eq.'nn_type_nntb')then
        if(count_nn_type_nntb.gt.0)then
          write(ounit,*)'Error: nn_type_nntb specified twice'
          stop
        endif
        count_nn_type_nntb=count_nn_type_nntb+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nn_type_nntb_local
        if(nn_type_nntb.ne.nn_type_nntb_local)then
          write(ounit,*)'ERROR in determination of nn_type_nntb'
          stop
        endif
!!
      elseif(keyword.eq.'random_number_type')then
        if(count_nran.gt.0)then
          write(ounit,*)'Error: random_number_type specified twice'
          stop
        endif
        count_nran=count_nran+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,nran
!!
      elseif(keyword.eq.'calculate_final_force')then
        if(count_lfinalforce.gt.0)then
          write(ounit,*)'Error: calculate_final_force specified twice'
          stop
        endif
        count_lfinalforce=count_lfinalforce+1
        lfinalforce=.true.
!!
      elseif(keyword.eq.'normalize_nodes')then
        if(count_lnormnodes.gt.0)then
          write(ounit,*)'Error: normalize_nodes specified twice'
          stop
        endif
        count_lnormnodes=count_lnormnodes+1
        lnormnodes=.true.
!!
      elseif(keyword.eq.'atom_energy')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'weight_constraint')then
        count_wconstraint=count_wconstraint+1
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'weighte_constraint')then
        count_wconstraint=count_wconstraint+1
!!      don't do anything here, just let it pass
!!
      elseif((keyword.eq.'global_symfunction_short')&
         .or.(keyword.eq.'global_symfunction_short_atomic'))then
!!      don't do anything here, just let it pass
!!
      elseif((keyword.eq.'global_symfunction_electrostatic')&
         .or.(keyword.eq.'global_symfunction_elec'))then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_symfunction_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_symfunction_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'symfunction_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'symfunction_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif((keyword.eq.'global_pairsymfunction_short')&
         .or.(keyword.eq.'global_symfunction_short_pair'))then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_symfunction_nntb_Hpair')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_symfunction_nntb_Hatomic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_symfunction_nntb_Hext')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_symfunction_nntb_S')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_pairsymfunction_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_hamsymfunction')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'pairsymfunction_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'symfunction_nntb_Hatomic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'basis')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'weights_min')then
        if(count_weights_min.gt.0)then
          write(ounit,*)'Error: weights_min specified twice'
          stop
        endif
        count_weights_min=count_weights_min+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,weights_min
!!
      elseif(keyword.eq.'weights_max')then
        if(count_weights_max.gt.0)then
          write(ounit,*)'Error: weights_max specified twice'
          stop
        endif
        count_weights_max=count_weights_max+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,weights_max
!!
      elseif(keyword.eq.'separate_bias_ini_short')then
        if(count_lseparatebiasini.gt.0)then
          write(ounit,*)'Error: separate_bias_ini_short specified twice'
          stop
        endif
        count_lseparatebiasini=count_lseparatebiasini+1
        lseparatebiasini=.true.
!!
      elseif(keyword.eq.'biasweights_min')then
        if(count_biasweights_min.gt.0)then
          write(ounit,*)'Error: biasweights_min specified twice'
          stop
        endif
        count_biasweights_min=count_biasweights_min+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,biasweights_min
!!
      elseif(keyword.eq.'biasweights_max')then
        if(count_biasweights_max.gt.0)then
          write(ounit,*)'Error: biasweights_max specified twice'
          stop
        endif
        count_biasweights_max=count_biasweights_max+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,biasweights_max
!!
      elseif(keyword.eq.'weightse_min')then
        if(count_weightse_min.gt.0)then
          write(ounit,*)'Error: weightse_min specified twice'
          stop
        endif
        count_weightse_min=count_weightse_min+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,weightse_min
!!
      elseif(keyword.eq.'weightse_max')then
        if(count_weightse_max.gt.0)then
          write(ounit,*)'Error: weightse_max specified twice'
          stop
        endif
        count_weightse_max=count_weightse_max+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,weightse_max
!!
      elseif(keyword.eq.'use_systematic_weights_short')then
        if(count_lsysweights.gt.0)then
          write(ounit,*)'Error: use_systematic_weights_short specified twice'
          stop
        endif
        count_lsysweights=count_lsysweights+1
        lsysweights=.true.
!!
      elseif(keyword.eq.'use_systematic_weights_electrostatic')then
        if(count_lsysweightse.gt.0)then
          write(ounit,*)'Error: use_systematic_weights_electrostatic specified twice'
          stop
        endif
        count_lsysweightse=count_lsysweightse+1
        lsysweightse=.true.
!!
      elseif(keyword.eq.'print_sensitivity')then
        if(count_lsens.gt.0)then
          write(ounit,*)'Error: print_sensitivity specified twice'
          stop
        endif
        count_lsens=count_lsens+1
        lsens=.true.
!!
      elseif(keyword.eq.'read_unformatted')then
        if(count_lreadunformatted.gt.0)then
          write(ounit,*)'Error: read_unformatted specified twice'
          stop
        endif
        count_lreadunformatted=count_lreadunformatted+1
        lreadunformatted=.true.
!!
      elseif(keyword.eq.'write_unformatted')then
        if(count_lwriteunformatted.gt.0)then
          write(ounit,*)'Error: write_unformatted specified twice'
          stop
        endif
        count_lwriteunformatted=count_lwriteunformatted+1
        lwriteunformatted=.true.
!!
      elseif(keyword.eq.'reset_kalman')then
        if(count_lresetkalman.gt.0)then
          write(ounit,*)'Error: reset_kalman specified twice'
          stop
        endif
        count_lresetkalman=count_lresetkalman+1
        lresetkalman=.true.
!!
      elseif(keyword.eq.'separate_kalman_short')then
        if(count_lsepkalman.gt.0)then
          write(ounit,*)'Error: separate_kalman_short specified twice'
          stop
        endif
        count_lsepkalman=count_lsepkalman+1
        lsepkalman=.true.
!!
      elseif(keyword.eq.'repeated_energy_update')then
        if(count_lrepeate.gt.0)then
          write(ounit,*)'Error: repeated_energy_update specified twice'
          stop
        endif
        count_lrepeate=count_lrepeate+1
        lrepeate=.true.
!!
      elseif(keyword.eq.'enforce_totcharge')then
        if(count_enforcetotcharge.gt.0)then
          write(ounit,*)'Error: enforce_totcharge specified twice'
          stop
        endif
        count_enforcetotcharge=count_enforcetotcharge+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,enforcetotcharge
!!
      elseif(keyword.eq.'check_forces')then
        if(count_lcheckf.gt.0)then
          write(ounit,*)'Error: check_forces specified twice'
          stop
        endif
        count_lcheckf=count_lcheckf+1
        lcheckf=.true.
!!
      elseif(keyword.eq.'write_fit_statistics')then
        if(count_lfitstats.gt.0)then
          write(ounit,*)'Error: write_fit_statistics specified twice'
          stop
        endif
        count_lfitstats=count_lfitstats+1
        lfitstats=.true.
!!
      elseif(keyword.eq.'fixed_short_energy_error_threshold')then
        if(count_lfixederrore.gt.0)then
          write(ounit,*)'Error: fixed_short_energy_error_threshold specified twice'
          stop
        endif
        count_lfixederrore=count_lfixederrore+1
        lfixederrore=.true.
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,fixederrore
!!
      elseif(keyword.eq.'fixed_short_force_error_threshold')then
        if(count_lfixederrorf.gt.0)then
          write(ounit,*)'Error: fixed_short_force_error_threshold specified twice'
          stop
        endif
        count_lfixederrorf=count_lfixederrorf+1
        lfixederrorf=.true.
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,fixederrorf
!!
      elseif(keyword.eq.'restrict_weights')then
        if(count_restrictw.gt.0)then
          write(ounit,*)'Error: restrict_weights specified twice'
          stop
        endif
        count_restrictw=count_restrictw+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,restrictw
!!
      elseif(keyword.eq.'screen_electrostatics')then
        if(count_lscreen.gt.0)then
          write(ounit,*)'Error: screen_electrostatics specified twice'
          stop
        endif
        count_lscreen=count_lscreen+1
        backspace(nnunit)
        lscreen=.true.
        read(nnunit,*,ERR=99)dummy,rscreen_onset,rscreen_cut
!!
      elseif(keyword.eq.'silent_mode')then
        if(count_lsilent.gt.0)then
          write(ounit,*)'Error: silent_mode specified twice'
          stop
        endif
        count_lsilent=count_lsilent+1
        backspace(nnunit)
        lsilent=.true.
!!
      elseif(keyword.eq.'prepare_md')then
        if(count_lpreparemd.gt.0)then
          write(ounit,*)'Error: prepare_md specified twice'
          stop
        endif
        count_lpreparemd=count_lpreparemd+1
        backspace(nnunit)
        lpreparemd=.true.
!!
      elseif(keyword.eq.'fitting_unit')then
        if(count_fittingunit.gt.0)then
          write(ounit,*)'Error: fitting_unit specified twice'
          stop
        endif
        count_fittingunit=count_fittingunit+1
        backspace(nnunit)
        read(nnunit,*,ERR=99)dummy,dummy2
        if(dummy2.eq.'eV')then
          fitting_unit=1
        elseif(dummy2.eq.'Ha')then
          fitting_unit=2
        else
          write(*,*)'Error: unknown energy unit specified for fitting ',dummy2
          stop
        endif
!!
      elseif(keyword.eq.'global_activation_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_activation_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_activation_pair')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'global_activation_ham')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_hidden_layers_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_hidden_layers_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_hidden_layers_pair')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_hidden_layers_ham')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_nodes_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_nodes_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_nodes_pair')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_nodes_ham')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_activation_short')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_activation_electrostatic')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_activation_pair')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'element_activation_ham')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'joint_energy_force_update')then
        if(count_ljointefupdate.gt.0)then
          write(ounit,*)'Error: joint_energy_force_update specified twice'
          stop !'
        endif
        count_ljointefupdate=count_ljointefupdate+1
        ljointefupdate=.true.
!!
      elseif(keyword.eq.'use_fixed_charges')then
        write(ounit,*)'Error: use_fixed_charges is obsolete, use electrostatic_type 3 instead '
        stop
!!
      elseif(keyword.eq.'use_omp_mkl')then
        if(count_lompmkl.gt.0)then
          write(ounit,*)'Error: use_omp_mkl specified twice'
          stop
        endif
        count_lompmkl=count_lompmkl+1
        lompmkl=.true.
!!
      elseif(keyword.eq.'nguyen_widrow_weights_short')then
        if(count_lnwweights.gt.0)then
          write(ounit,*)'Error: nguyen_widrow_weights_short specified twice'
          stop
        endif
        count_lnwweights=count_lnwweights+1
        lnwweights=.true.
!!
      elseif(keyword.eq.'nguyen_widrow_weights_ewald')then
        if(count_lnwweightse.gt.0)then
          write(ounit,*)'Error: nguyen_widrow_weights_ewald specified twice'
          stop
        endif
        count_lnwweightse=count_lnwweightse+1
        lnwweightse=.true.
!!
      elseif(keyword.eq.'fixed_charge')then
!!      don't do anything here, just let it pass
!!
      elseif(keyword.eq.'print_all_short_weights')then
        write(pstring(1:1),'(a1)')'1'
!!
      elseif(keyword.eq.'print_all_electrostatic_weights')then
        write(pstring(2:2),'(a1)')'1'
!!
      elseif(keyword.eq.'print_all_deshortdw')then
        write(pstring(3:3),'(a1)')'1'
!!
      elseif(keyword.eq.'print_all_dfshortdw')then
        write(pstring(4:4),'(a1)')'1'
!!
      else
        write(ounit,*)'Error: unknown keyword in input.nn: ',keyword
        stop
      endif
!!
      goto 10
!!
 20   continue
      close(nnunit)
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
 99   continue
      write(ounit,*)'Error: keyword ',keyword
      write(ounit,*)'is missing arguments '
      stop

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
