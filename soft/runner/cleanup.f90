!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine cleanup()
!!
      use mpi_mod
      use fileunits
      use timings
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use basismod
      use symfunctions
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      call zerotime(dayfinalize,timefinalizestart,timefinalizeend)
      call abstime(timefinalizestart,dayfinalize)
!!
!! deallocate everything
      if(lshort.and.(nn_type_short.eq.1))then
        deallocate(weights_short_atomic)
        deallocate(symfunction_short_atomic_list)
        deallocate(num_funcvalues_short_atomic)
        deallocate(windex_short_atomic)
        deallocate(num_layers_short_atomic)
        deallocate(actfunc_short_atomic)
        deallocate(nodes_short_atomic)
        deallocate(num_weights_short_atomic)
        deallocate(function_type_short_atomic)
        deallocate(symelement_short_atomic)
        deallocate(funccutoff_short_atomic)
        deallocate(eta_short_atomic)
        deallocate(zeta_short_atomic)
        deallocate(lambda_short_atomic)
        deallocate(rshift_short_atomic)
        if(mode.eq.2)deallocate(kalmanlambda)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        deallocate(weights_short_pair)
        deallocate(symfunction_short_pair_list)
        deallocate(num_funcvalues_short_pair)
        deallocate(windex_short_pair)
        deallocate(num_layers_short_pair)
        deallocate(actfunc_short_pair)
        deallocate(nodes_short_pair)
        deallocate(num_weights_short_pair)
        deallocate(function_type_short_pair)
        deallocate(symelement_short_pair)
        deallocate(funccutoff_short_pair)
        deallocate(eta_short_pair)
        deallocate(zeta_short_pair)
        deallocate(lambda_short_pair)
        deallocate(rshift_short_pair)
        if(mode.eq.2)deallocate(kalmanlambdap)
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        deallocate(weights_elec)
        deallocate(symfunction_elec_list)
        deallocate(num_funcvalues_elec)
        deallocate(windex_elec)
        deallocate(num_layers_elec)
        deallocate(actfunc_elec)
        deallocate(nodes_elec)
        deallocate(num_weights_elec)
        deallocate(function_type_elec)
        deallocate(symelement_elec)
        deallocate(funccutoff_elec)
        deallocate(eta_elec)
        deallocate(zeta_elec)
        deallocate(lambda_elec)
        deallocate(rshift_elec)
        if(mode.eq.2)deallocate(kalmanlambdae)
      endif
      if(lnntb)then
        if(nntb_flag(0))then
          deallocate(weights_ham)
          deallocate(symfunction_ham_list)
          deallocate(num_funcvalues_ham)
          deallocate(windex_ham)
          deallocate(num_layers_ham)
          deallocate(actfunc_ham)
          deallocate(nodes_ham)
          deallocate(num_weights_ham)
          deallocate(function_type_ham)
          deallocate(symelement_ham)
          deallocate(funccutoff_ham)
          deallocate(eta_ham)
          deallocate(zeta_ham)
          deallocate(lambda_ham)
          deallocate(rshift_ham)
        endif
        if(nntb_flag(3))then
          deallocate(weights_hextoff)
          deallocate(symfunction_hextoff_list)
          deallocate(num_funcvalues_hextoff)
          if(mode.ge.2)then
            deallocate(windex_hextoff)
            deallocate(num_layers_hextoff)
            deallocate(actfunc_hextoff)
            deallocate(nodes_hextoff)
            deallocate(num_weights_hextoff)
          endif
          deallocate(function_type_hextoff)
          deallocate(symelement_hextoff)
          deallocate(funccutoff_hextoff)
          deallocate(eta_hextoff)
          deallocate(zeta_hextoff)
          deallocate(lambda_hextoff)
          deallocate(rshift_hextoff)
!!          deallocate(tripletindex)
        endif
!! todo:        if(mode.eq.2)deallocate(kalmanlambdae)
      endif

      deallocate(nucelem)
      deallocate(element)
      if(allocated(atomrefenergies))deallocate(atomrefenergies)
      if(allocated(fixedcharge))deallocate(fixedcharge)
      if(allocated(elempair))deallocate(elempair)

!! deallocate symmetry functions

      if(allocated(basis))deallocate(basis)
      if(allocated(num_basis))deallocate(num_basis)

      if(allocated(vdw_param))deallocate(vdw_param)
!!
!!    for timing let all processes wait here
      call mpi_barrier(mpi_comm_world,mpierror)
!!
      call abstime(timefinalizeend,dayfinalize)
      if(lfinetime)then
        call printtimings()
      endif
!! Print final summary
!!
      if(mpirank.eq.0)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f14.3)')' Total runtime (s)  : ',runtimeend-runtimestart
        write(ounit,'(a,f14.3)')' Total runtime (min): ',(runtimeend-runtimestart)/60.d0
        write(ounit,'(a,f14.3)')' Total runtime (h)  : ',(runtimeend-runtimestart)/3600.d0
        write(ounit,*)'Normal termination of RuNNer'
        write(ounit,*)'-------------------------------------------------------------'
        if(ounit.ne.6)then
          close(ounit)
        endif
        close(debugunit)
      endif ! mpirank.eq.0
!!
      return
      end
