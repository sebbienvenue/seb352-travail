!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine checkinputnn()
!!
      use nnflags
      use globaloptions
      use inputnncounters
      use fittingoptions
      use predictionoptions
      use mode1options
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use fileunits
!!
      implicit none
!!
      integer i1
!!
      if(ldynforcegroup.and.(mode.eq.2).and.(nforcegroup.gt.1))then
        write(ounit,*)'ERROR: dynamic_force_grouping cannot be used together with force grouping'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.2).and.lshort.and.(nn_type_short.eq.2))then
        write(ounit,*)'ERROR: charges can be derived from second output node'
        write(ounit,*)'of short range NN in atomic case only!'   
        stop
      endif

      if(lscreen)then
        if(rscreen_onset.gt.rscreen_cut)then
          write(ounit,*)'ERROR: rscreen_onset .gt. rscreen_cut in screen_electrostatics'
          stop
        endif
        if(rscreen_onset.lt.0.0d0)then
          write(ounit,*)'ERROR: rscreen_onset .lt. 0 in screen_electrostatics'
          stop
        endif
        if(rscreen_cut.lt.0.0d0)then
          write(ounit,*)'ERROR: rscreen_cut .lt. 0 in screen_electrostatics'
          stop
        endif
      endif
!!
      if(lrandomtrain.and.lmixpoints)then
        write(ounit,*)'ERROR: random_order_training cannot be combined with mix_all_points'
        write(ounit,*)'Keyword mix_all_points is more general!'
        stop
      endif
!!
      if((cutoff_type.ne.1).and.(nn_type_short.ne.1))then
        write(ounit,*)'ERROR: cutoff_type .ne.1 can be used only for nn_type_short 1'
        stop
      endif
!!
      if((mode.eq.2).and.(lshort).and.(luseforces).and.(lfgroupbystruct))then
        if(nforcegroup.gt.1)then
          write(ounit,*)'ERROR: force_grouping_by_structure cannot be combined with short_force_group>1 '
          stop
        endif
      endif
!!
      if(noisee.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_energy must not be negative ',noisee
        stop
      endif
!!
      if(lwritetrainforces.and.(mode.eq.2).and.lelec.and.(.not.lshort))then
        write(ounit,*)'ERROR: writing force error does not make sense for charge fitting only'
        stop
      endif
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.luseforces)then
        write(ounit,*)'ERROR: electrostatic forces cannot be used for charge fitting'
        stop
      endif
!!
      if(noisef.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_force must not be negative ',noisef
        stop
      endif
!!
      if(noiseq.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_charge must not be negative ',noiseq
        stop
      endif
!!
      if(lsysweights.and.lnwweights)then
        write(ounit,'(a)')'Error: Cannot use systematic_weights_short and nguyen_widrow_weights_short together!' 
        stop
      endif
!!
      if(lsysweightse.and.lnwweightse)then
        write(ounit,'(a)')'Error: Cannot use systematic_weights_ewald and nguyen_widrow_weights_ewald together!' 
        stop
      endif
!!
      if(lnormnodes.and.lnwweights)then
        write(ounit,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_short together!' 
        stop
      endif
!!
      if(lnormnodes.and.lnwweightse)then
        write(ounit,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_ewald together!' 
        stop
      endif
!!
      if((fitmode.ne.1).and.(fitmode.ne.2))then
        write(ounit,'(a)')'Error: Unknown fitmode in readinput ',fitmode
        stop
      endif
!!
      if((count_kalmanthreshold.eq.1).and.(count_lfixederrore.eq.1))then
        write(ounit,'(2a)')'Error: short_energy_error_threshold cannot be used ',&
          'in combination with fixed_short_energy_error_threshold'
        stop
      endif
!!
      if((count_kalmanthresholdf.eq.1).and.(count_lfixederrorf.eq.1))then
        write(ounit,'(a)')'Error: short_force_error_thresholdf cannot be used in combination with fixed_short_force_error_threshold'
        stop
      endif
!!
      if(lfinalforce.and.luseforces.and.(mode.eq.2))then
        write(ounit,'(a)')'Error: calculate_final_forces cannot be used in combination with use_short_forces'
        stop
      endif
!!
      if(count_mode.eq.0)then
        write(ounit,*)'Error: runner_mode is not specified'
        stop
      endif
!!
      if((.not.lshort).and.(.not.lelec))then
        write(ounit,*)'Error: short range and electrostatic NNs are switched off'
        stop
      endif    
!!
      if(lshort.and.(maxnum_layers_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(ounit,*)'Error: global_hidden_layers_short is not specified'
        stop
      endif
!!
      if(lshort.and.(maxnum_layers_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(ounit,*)'Error: global_hidden_layers_pair is not specified'
        stop
      endif
!!
      if((maxnum_layers_ham.eq.0).and.(nn_type_nntb.eq.1))then
        write(ounit,*)'Error: global_hidden_layers_ham is not specified'
        stop
      endif
!!
      if(lelec.and.(nn_type_elec.eq.0))then
        write(ounit,*)'Error: electrostatic_type is not specified'
        stop
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(maxnum_layers_elec.eq.0))then
        write(ounit,*)'Error: global_hidden_layers_electrostatic is not specified'
        stop
      endif
!!
      if(lshort.and.(count_nodes_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(ounit,*)'Error: global_nodes_short is not specified'
        stop
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(count_nodes_elec.eq.0))then
        write(ounit,*)'Error: global_nodes_electrostatic is not specified'
        stop
      endif
!!
      if(lshort.and.(count_nodes_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(ounit,*)'Error: global_nodes_pair is not specified'
        stop
      endif
!!
      if((count_nodes_ham.eq.0).and.(nn_type_nntb.eq.1))then
        write(ounit,*)'Error: global_nodes_ham is not specified'
        stop
      endif
!!
      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).gt.1)then
            write(ounit,*)'Error: More than 1 output node currently does '
            write(ounit,*)'make sense in short range NN'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).eq.0)then
            write(ounit,*)'Error: output_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).gt.1)then
            write(ounit,*)'Error: More than 1 output node currently does '
            write(ounit,*)'make sense in electrostatic NN'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).eq.0)then
            write(ounit,*)'Error: output_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,npairs
        if(lshort.and.(nn_type_short.eq.2))then
          if(nodes_short_pair(maxnum_layers_short_pair,i1).gt.1)then
            write(ounit,*)'Error: More than 1 output node currently does '
            write(ounit,*)'make sense in short range pair NN'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,npairs
        if(lshort.and.(nn_type_short.eq.2))then
          if(nodes_short_pair(maxnum_layers_short_pair,i1).eq.0)then
            write(ounit,*)'Error: output_nodes_pair is 0'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).eq.0)then
            write(ounit,*)'Error: input_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).ne.num_funcvalues_short_atomic(i1))then
            write(ounit,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)',&
              num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            write(ounit,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).eq.0)then
            write(ounit,*)'Error: input_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).ne.num_funcvalues_elec(i1))then
            write(ounit,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)',&
              num_funcvalues_elec(i1),nodes_elec(0,i1)
            write(ounit,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,npairs
        if(lshort.and.(nn_type_short.eq.2))then
          if(nodes_short_pair(0,i1).eq.0)then
            write(ounit,*)'ERROR: input_nodes_pair is 0'
            write(ounit,*)'If this pair is not needed for your structures,'
            write(ounit,*)'please still define one dummy pair symmetry function.'
            stop
          endif
        endif
      enddo ! i1
!!
      do i1=1,npairs
        if(lshort.and.(nn_type_short.eq.2))then
          if(nodes_short_pair(0,i1).ne.num_funcvalues_short_pair(i1))then
            write(ounit,*)'Error: num_funcvalues_short_pair .ne. nodes_short_pair(0)',&
              num_funcvalues_short_pair(i1),nodes_short_pair(0,i1)
            write(ounit,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1
!!
      if(lshort.and.(nn_type_short.eq.1))then
        if(count_global_activation_short_atomic.eq.0)then
          write(ounit,*)'Error: global_activation_short is not specified'
          stop
        endif
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        if(count_global_activation_elec.eq.0)then
          write(ounit,*)'Error: global_activation_ewald is not specified'
          stop
        endif
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        if(count_global_activation_short_pair.eq.0)then
          write(ounit,*)'Error: global_activation_pair is not specified'
          stop
        endif
      endif
!!
      if(lnntb)then
        if((count_global_activation_ham.eq.0).and.(nn_type_nntb.eq.1))then
          write(ounit,*)'Error: global_activation_ham is not specified'
          stop
        endif
      endif
!!
      if(lelec.and.(count_ewaldalpha.eq.0))then
        write(ounit,*)'Error: ewald_alpha must be specified for electrostatic NN'
        stop
      endif
!!
      if(lelec.and.(ewaldalpha.le.0))then
        write(ounit,*)'Error: ewald_alpha must be positive ',ewaldalpha
        stop
      endif
!!
      if(lelec.and.(count_ewaldcutoff.eq.0))then
        write(ounit,*)'Error: ewald_cutoff must be specified for electrostatic NN'
        stop
      endif
!!
      if(lelec.and.(ewaldcutoff.le.0))then
        write(ounit,*)'Error: ewald_cutoff must be positive ',ewaldcutoff
        stop
      endif
!!
      if(lelec.and.(count_ewaldkmax.eq.0))then
        write(ounit,*)'Error: ewald_kmax must be specified for electrostatic NN'
        stop
      endif
!!
      if((.not.lshort).and.(luseforces))then
        write(ounit,*)'### WARNING ### switching off use_short_forces because no short range NN is used'
        luseforces=.false. 
      endif
!!
      if(lelec.and.(.not.luseatomcharges))then
        write(ounit,*)'### WARNING ### use_atom_charges is switched on for electrostatic NN'
        luseatomcharges=.true. 
      endif
!!
      if(lshort.and.(luseatomenergies))then
        write(ounit,*)'### WARNING ### use_atom_energies is switched off (not implemented)'
        luseatomenergies=.false. 
      endif
!!
      if((.not.lshort).and.(lremoveatomenergies))then
        write(ounit,*)'### WARNING ### remove_atom_energies is switched on without short range NN'
      endif
!!
      if(lelec.and.(lchargeconstraint))then
        write(ounit,'(a)')' ### WARNING ### use_charge_constraint is not maintained at the moment and might fail'
      endif
!!
      if(count_iseed.eq.0)then
        write(ounit,*)'### WARNING ### no random_seed specified, using default '
      endif
!!
      if(nenergygroup.gt.nblock)then
        nenergygroup=nblock
        write(ounit,*)'### WARNING ### reducing nenergygroup to nblock'
      endif
!!
      if((mode.eq.2).and.(count_epochs.eq.0))then
        write(ounit,*)'Error: points_in_memory not specified'
        stop
      endif
!!
      if(count_nelem.eq.0)then
        write(ounit,*)'Error: number_of_elements not specified'
        stop
      endif
!!
      if(count_element.eq.0)then
        write(ounit,*)'Error: elements not specified'
        stop
      endif
!!
      if((mode.eq.1).and.(count_splitthres.eq.0))then
        write(ounit,*)'Error: test_fraction not specified'
        stop
      endif
!!
      if((mode.eq.2).and.(lshort).and.(optmodee.eq.1).and.(count_kalmanlambda.eq.0))then
        write(ounit,*)'Error: kalman_lambda_short not specified'
        stop
      endif
!!
      if((mode.eq.2).and.(lelec).and.(nn_type_elec.eq.1).and.(optmodeq.eq.1).and.(count_kalmanlambdae.eq.0))then
        write(ounit,*)'Error: kalman_lambda_charge not specified'
        stop
      endif
!!
      if((mode.eq.2).and.(lelec).and.(optmodeq.eq.1).and.(lchargeconstraint).and.(count_kalmanlambdac.eq.0))then
        write(ounit,*)'Error: kalman_lambda_charge_constraint not specified'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_atomic.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_atomic.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_pair.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_pair.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_ham.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_ham keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_ham.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_ham keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_elec.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_elec keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_elec.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_elec keyword'
        stop
      endif
!!
      if((count_scmin_short_atomic.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_atomic.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_short_pair.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_pair.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_ham.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_ham requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_ham.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_ham requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_elec.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_elec.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if(scmin_short_atomic.ge.scmax_short_atomic)then
        write(ounit,'(a)')'Error: scale_min_short .ge. scale_max_short'
        stop
      endif
!!
      if(scmin_short_pair.ge.scmax_short_pair)then
        write(ounit,'(a)')'Error: scale_min_short_pair .ge. scale_max_short_pair'
        stop
      endif
!!
      if(scmin_ham.ge.scmax_ham)then
        write(ounit,'(a)')'Error: scale_min_ham .ge. scale_max_ham'
        stop
      endif
!!
      if(scmin_elec.ge.scmax_elec)then
        write(ounit,'(a)')'Error: scale_min_elec .ge. scale_max_elec'
        stop
      endif
!!
      if((mode.eq.2).and.(lshort).and.(optmodee.eq.1).and.(count_kalmannue.eq.0))then
        write(ounit,*)'Error: kalman_nue_short not specified'
        stop
      endif
!!
      if((mode.eq.2).and.(lelec).and.(nn_type_elec.eq.1).and.(optmodeq.eq.1).and.(count_kalmannuee.eq.0))then
        write(ounit,*)'Error: kalman_nue_charge not specified'
        stop
      endif
!!
      if((mode.eq.2).and.(lelec).and.(optmodeq.eq.1).and.(lchargeconstraint).and.(count_kalmannuec.eq.0))then
        write(ounit,*)'Error: kalman_nue_charge_constraint not specified'
        stop
      endif
!!
      if((mode.eq.2).and.lsavekalman.and.(.not.((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))))then
        write(ounit,*)'### WARNING ### switching off save_kalman_matrices, because no Kalman filter is used'
        lsavekalman=.false.
      endif
!!
      if((mode.eq.2).and.lrestkalman.and.(.not.((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))))then
        write(ounit,*)'### WARNING ### switching off read_kalman_matrices, because no Kalman filter is used'
        lrestkalman=.false.
      endif
!!
      if((nn_type_vdw.ne.0).and.(count_cutoffvdw.eq.0))then
        write(ounit,*)'ERROR: Please specify cutoff_vdw when using vdW ',cutoffvdw
        stop
      endif
!!
      if(lupdatebyelement.and.lchargeconstraint)then
        lchargeconstraint=.false.
        if(mode.eq.2)then
          write(ounit,*)'### WARNING ### lchargeconstraint is switched off because of lupdatebyelement'
        endif
      endif
!!
      if(lshort.and.lupdatebyelement.and.(mode.eq.2))then
        write(ounit,*)'### WARNING ### lupdatebyelement works only for charges and forces'
      endif
!!
      if(luseworste.and.lshort.and.(energyrnd.lt.1.0d0))then
        energyrnd=1.0d0
        write(ounit,*)'### WARNING ### luseworste overrides energyrnd: ',energyrnd
      endif
!!
      if(luseworstf.and.lshort.and.luseforces.and.(forcernd.lt.1.0d0))then
        forcernd=1.0d0
        write(ounit,*)'### WARNING ### luseworstf overrides forcernd: ',forcernd
      endif
!!
      if(luseworstq.and.lelec.and.(nn_type_elec.eq.1).and.(chargernd.lt.1.0d0))then
        chargernd=1.0d0
        write(ounit,*)'### WARNING ### luseworstq overrides chargernd: ',chargernd
      endif
!!
      if(dampw.gt.1.0d0)then
        write(ounit,*)'Error: dampw must not be larger than 1.0d0 ',dampw
        stop
      endif
!!
      if(ldostress.and.(.not.ldoforces))then
        write(ounit,*)'### WARNING ### Analytic stress is requested without forces'
        write(ounit,*)'Switching on calculation of analytic forces'
        ldoforces=.true.
      endif
!!
      if(ldostress.and.(mode.eq.1))then
        write(ounit,*)'### WARNING ### switching off stress calculation in mode 1 for increased performance'
        ldostress=.false.
      endif
!!
      if((count_wconstraint.gt.0).and.(.not.lfixweights))then
        write(ounit,*)'Error: weight constraints are specified without fix_weights keyword'
        stop
      endif
!!
      if((mode.eq.2).and.(lprecond).and.(luseoldweightsshort).and.(lshort))then
        write(ounit,*)'WARING: switching off precondition_weights because old short range weights are used'
        lprecond=.false.
      endif
!!
      if((mode.eq.2).and.(lprecond).and.(luseoldweightscharge).and.(lelec).and.(nn_type_elec.eq.1))then
        write(ounit,*)'WARING: switching off precondition_weights because old charge weights are used'
        lprecond=.false.
      endif
!!
      if((count_wconstraint.eq.0).and.(lfixweights))then
        write(ounit,*)'Error: no weights constrained but keyword fix_weights has been selected'
        stop
      endif
!!
      if(weights_min.ge.weights_max)then
        write(ounit,*)'Error: weights_min > weights_max'
        stop
      endif
!!
      if(biasweights_min.ge.biasweights_max)then
        write(ounit,*)'Error: biasweights_min > biasweights_max'
        stop
      endif
!!
      if(weightse_min.ge.weightse_max)then
        write(ounit,*)'Error: weightse_min > weightse_max'
        stop
      endif
!!
      if(kalman_dampe.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_short must be non-negative ',kalman_dampe
        stop
      endif
!!
      if(kalman_dampf.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_force must be non-negative ',kalman_dampf
        stop
      endif
!!
      if(kalman_dampq.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_charge must be non-negative ',kalman_dampq
        stop
      endif
!!
      if(ljointefupdate)then
        write(ounit,*)'ERROR: joint_energy_force_update is not implemented for lelec and nn_type_elec 2'
        stop
      endif
!!
      if(ljointefupdate)then
        if(optmodee.ne.optmodef)then
          write(ounit,*)'Error: joint_energy_force_update requires to use the'
          write(ounit,*)'same optimization algorithm for energy and forces'
          stop
        endif
        if(.not.luseforces)then
          write(ounit,*)'Error: switch on use_short_forces for joint_energy_force_update'
          stop
        endif
        if(lrepeate)then
          write(ounit,*)'ERROR: repeated energy update cannot be combined with joint energy and force update'
          stop 
        endif
        if(forcernd.lt.1.0d0)then
          write(ounit,*)'ERROR: joint energy and force update requires force_fraction = 1.0d0'
          stop
        endif
        if(luseworste)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_energies'
          stop
        endif
        if(luseworstf)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_forces'
          stop
        endif
        if(nenergygroup.gt.1)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_energy_group > 1'
          stop
        endif
        if(nforcegroup.gt.1)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_force_group > 1'
          stop
        endif
        if(kalmanthresholdf.gt.0.0d0)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_force_error_threshold > 0.0'
          stop
        endif
!!        write(ounit,*)'### WARNING ### joint_energy_force_update enforces these settings:'
!!        write(ounit,*)'short_force_fraction = 1.0'
!!        forcernd=1.0d0
!!        write(ounit,*)'update_worst_short_energies = .false.'
!!        luseworste=.false.
!!        write(ounit,*)'update_worst_short_forces = .false.'
!!        luseworstf=.false.
!!        write(ounit,*)'short_energy_group = 1'
!!        nenergygroup=1
!!        write(ounit,*)'short_force_group = 1'
!!        nforcegroup=1
!!        write(ounit,*)'kalmanthresholdf = 0.0'
!!        kalmanthresholdf=0.0d0
      endif
!!
      if(maxforce.le.0.0d0)then
        write(ounit,*)'Error: max_force must not be negative ',maxforce
        stop
      endif
!!
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_short_atomic(i1).ne.nodes_short_atomic(0,i1))then
            write(ounit,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)'
            write(ounit,*)i1,num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            stop
          endif
        enddo
      elseif(lshort.and.(nn_type_short.eq.2))then
        do i1=1,npairs
          if(num_funcvalues_short_pair(i1).ne.nodes_short_pair(0,i1))then
            write(ounit,*)'Error: num_funcvalues_short_pair .ne. nodes_short_pair(0)'
            write(ounit,*)i1,num_funcvalues_short_pair(i1),nodes_short_pair(0,i1)
            stop
          endif
        enddo ! i1
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_elec(i1).ne.nodes_elec(0,i1))then
            write(ounit,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)'
            write(ounit,*)i1,num_funcvalues_elec(i1),nodes_elec(0,i1)
            stop
          endif
        enddo ! i1
      endif
!!
      if((mode.eq.3).and.(max_num_atoms.lt.nblock).and.(nn_type_short.eq.1).and.lshort)then
        write(ounit,*) 'WARNING: reducing points_in_memory to max_num_atoms ',max_num_atoms
        nblock=max_num_atoms
      endif
!!
      return
      end
