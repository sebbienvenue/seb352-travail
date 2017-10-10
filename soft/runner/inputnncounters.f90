      module inputnncounters 

      implicit none

      integer count_analyze_error_energy_step
      integer count_analyze_error_force_step
      integer count_analyze_error_charge_step
      integer count_mode
      integer count_nn_type_elec
      integer count_vdw_type
      integer count_cutoffvdw
      integer count_ldebug
      integer count_paramode
      integer count_lshort
      integer count_lelec
      integer count_lnntb
      integer count_lham
      integer count_loverlap
      integer count_lhexton
      integer count_lhextoff
      integer count_ldens
      integer count_num_layers_short_atomic
      integer count_num_layers_elec
      integer count_num_layers_short_pair
      integer count_num_layers_ham
      integer count_num_layers_s
      integer count_num_layers_hextoff
      integer count_num_layers_hexton
      integer count_num_layers_dens
      integer count_nodes_short_atomic
      integer count_nodes_elec
      integer count_nodes_short_pair
      integer count_nodes_ham
      integer count_nodes_s
      integer count_nodes_hexton
      integer count_nodes_hextoff
      integer count_nodes_dens
      integer count_global_activation_short_atomic
      integer count_global_activation_elec
      integer count_global_activation_short_pair
      integer count_global_activation_ham
      integer count_global_activation_s
      integer count_global_activation_hexton
      integer count_global_activation_hextoff
      integer count_global_activation_dens
      integer count_ewaldalpha
      integer count_ewaldkmax
      integer count_ewaldcutoff
      integer count_short_energy_group
      integer count_short_force_group
      integer count_charge_group
      integer count_luseforces
      integer count_energyrnd
      integer count_forcernd
      integer count_chargernd
      integer count_luseatomcharges
      integer count_luseatomenergies
      integer count_lremoveatomenergies
      integer count_lchargeconstraint
      integer count_fitethres
      integer count_fitfthres
      integer count_rmin
      integer count_optmodee
      integer count_optmodef
      integer count_optmodeq
      integer count_optmodedens
      integer count_optmodes
      integer count_optmodehexton
      integer count_optmodehextoff
      integer count_optmodeham
      integer count_iseed
      integer count_nblock
      integer count_epochs
      integer count_nelem
      integer count_element
      integer count_iwriteweight
      integer count_lwritetmpweights
      integer count_splitthres
      integer count_kalmanthreshold
      integer count_kalmanthresholdf
      integer count_kalmanthresholde
      integer count_kalmanthresholdc
      integer count_kalmanlambda
      integer count_kalmanlambdae
      integer count_kalmanlambdac
      integer count_kalmannue
      integer count_kalmannuee
      integer count_kalmannuec
      integer count_kalman_dampe
      integer count_kalman_dampf
      integer count_kalman_dampq
      integer count_steepeststepe
      integer count_steepeststepf
      integer count_steepeststepq
      integer count_scalefactorf
      integer count_scalefactorq
      integer count_lrandomtrain
      integer count_lscalesym
      integer count_lcentersym
      integer count_luseoldweightsshort
      integer count_luseoldweightscharge
      integer count_luseoldweights_s
      integer count_luseoldweights_hexton
      integer count_luseoldweights_hextoff
      integer count_luseoldweights_dens 
      integer count_lglobalfit
      integer count_lsavekalman
      integer count_lrestkalman
      integer count_elemupdate
      integer count_luseworste
      integer count_luseworstf
      integer count_luseworstq
      integer count_growth
      integer count_dampw
      integer count_lfixweights
      integer count_ldoforces
      integer count_ldostress
      integer count_lfinetime
      integer count_lfinetimeepoch
      integer count_lwritepdb
      integer count_lwritexyz
      integer count_lwritepov
      integer count_lwritepw
      integer count_lwritetrainpoints
      integer count_lwritetrainhextoff
      integer count_lwritetrainforces
      integer count_lwritetraincharges
      integer count_wconstraint
      integer count_weights_min
      integer count_weights_max
      integer count_lseparatebiasini
      integer count_biasweights_min
      integer count_biasweights_max
      integer count_weightse_min
      integer count_weightse_max
      integer count_fittingunit
      integer count_ljointefupdate
      integer count_nn_type_short
      integer count_nn_type_nntb
      integer count_nran
      integer count_enforcetotcharge  
      integer count_lsysweights 
      integer count_lsysweightse 
      integer count_lsens 
      integer count_lreadunformatted 
      integer count_lwriteunformatted 
      integer count_lresetkalman 
      integer count_lsepkalman 
      integer count_lrepeate 
      integer count_lfinalforce 
      integer count_maxforce 
      integer count_maxenergy 
      integer count_lcheckf 
      integer count_lfitstats 
      integer count_lfixederrore 
      integer count_lfixederrorf 
      integer count_lompmkl 
      integer count_lnormnodes
      integer count_restrictw
      integer count_fitmode
      integer count_lanalyzeerror
      integer count_lnwweights
      integer count_lnwweightse
      integer count_scmin_short_atomic  
      integer count_scmax_short_atomic
      integer count_scmin_short_pair  
      integer count_scmax_short_pair
      integer count_scmin_ham  
      integer count_scmax_ham
      integer count_scmin_s
      integer count_scmax_s
      integer count_scmin_hexton
      integer count_scmax_hexton
      integer count_scmin_hextoff
      integer count_scmax_hextoff
      integer count_scmin_dens
      integer count_scmax_dens
      integer count_scmin_elec     
      integer count_scmax_elec    
      integer count_luseoldscaling    
      integer count_lprecond    
      integer count_linionly    
      integer count_noisee    
      integer count_noisef    
      integer count_noiseq    
      integer count_lprintconv
      integer count_lprintmad
      integer count_lfgroupbystruct
      integer count_lqgroupbystruct
      integer count_cutoff_type
      integer count_lmixpoints
      integer count_lscreen
      integer count_lsilent
      integer count_lpreparemd
      integer count_lpearson_correlation
      integer count_lweightanalysis
      integer count_lenvironmentanalysis
      integer count_lfindcontradictions
      integer count_lenforcemaxnumneighborsatomic
      integer count_lmd
      integer count_ldynforcegroup
      integer count_shuffle_weights_short_atomic
      integer count_hextoff_training_triplet
      integer count_ldetect_saturation
      integer count_ldataclustering
      integer count_luseedkalman !! MG: ED-Kalman counter
      integer count_ledforcesv2  !! MG: 2nd variant of ED-Kalman force fitting
      integer count_lanalyzecomposition
      integer count_lprintdateandtime
      integer count_lenableontheflyinput
      integer count_lcheckinputforces
      integer count_lionforcesonly
      integer count_useipi

      end module inputnncounters 

