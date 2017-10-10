!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine distribute_fittingoptions()
!!
      use mpi_mod
      use nnflags
      use globaloptions
      use fittingoptions
!!
      implicit none
!!
      call mpi_bcast(ngrowth,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(growthstep,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nenergygroup,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nforcegroup,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nchargegroup,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(optmodee,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(optmodef,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(optmodeq,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nepochs,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(elemupdate,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(fitmode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(fitting_unit,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nshuffle_weights_short_atomic,1,mpi_integer,0,mpi_comm_world,mpierror)

      call mpi_bcast(analyze_error_energy_step,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(analyze_error_force_step,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(analyze_error_charge_step,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(forcernd,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(energyrnd,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(chargernd,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(worste,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(worstf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(worstq,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(fixederrore,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(fixederrorf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(noisee,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(noisef,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(noiseq,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(steepeststepe,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(steepeststepf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(steepeststepq,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalman_dampe,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalman_dampf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalman_dampq,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(scalefactorf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(scalefactorq,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmanthreshold,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmanthresholdf,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmanthresholde,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmanthresholdc,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmannue,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmannuee,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(kalmannuec,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(dampw,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(weights_min,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(weights_max,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(biasweights_min,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(biasweights_max,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(weightse_min,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(weightse_max,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(shuffle_weights_short_atomic,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(dataclusteringthreshold1,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(dataclusteringthreshold2,1,mpi_real8,0,mpi_comm_world,mpierror)

      call mpi_bcast(restrictw,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(maxforce,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(maxenergy,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(deltagthres,1,mpi_real8,0,mpi_comm_world,mpierror)
      call mpi_bcast(deltafthres,1,mpi_real8,0,mpi_comm_world,mpierror)

      call mpi_bcast(lwritetmpweights,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsavekalman,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lrestkalman,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lgrowth,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseworste,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseworstf,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseworstq,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfixederrore,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfixederrorf,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfgroupbystruct,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lqgroupbystruct,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lmixpoints,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lupdatebyelement,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lrandomtrain,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldampw,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ljointefupdate,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsysweights,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsysweightse,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweightsshort,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweightscharge,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweights_ham,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweights_s,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweights_hexton,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweights_hextoff,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(luseoldweights_dens,1,mpi_logical,0,mpi_comm_world,mpierror)

      call mpi_bcast(luseoldscaling,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(linionly,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lprecond,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lprintconv,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfinalforce,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lnwweights,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lnwweightse,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lglobalfit,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lanalyzeerror,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritetrainpoints,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritetraincharges,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lwritetrainforces,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfixweights,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lrepeate,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lsepkalman,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfitstats,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lchargeconstraint,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lresetkalman,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lseparatebiasini,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lweightanalysis,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lfindcontradictions,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lshuffle_weights_short_atomic,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(ldataclustering,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lenableontheflyinput,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lionforcesonly,1,mpi_logical,0,mpi_comm_world,mpierror)

      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(kalmanlambda,nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(kalmanlambdap,npairs,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(kalmanlambdae,nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(kalmanlambdac,1,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      if(lnntb) then
        call mpi_bcast(kalmanlambdahextoff,ntriplets,mpi_real8,0,mpi_comm_world,mpierror)
      endif

      return
      end
