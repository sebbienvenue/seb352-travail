!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initnn.f90
!!
      subroutine distribute_symfunctions()
!!
      use mpi_mod
      use globaloptions
      use symfunctions
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(function_type_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_atomic,maxnum_funcvalues_short_atomic*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_atomic,maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(function_type_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_pair,maxnum_funcvalues_short_pair*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_pair,maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(function_type_elec,maxnum_funcvalues_elec*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_elec,maxnum_funcvalues_elec*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_elec,maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_elec,maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_elec,maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_elec,maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_elec,maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lnntb)then
        call mpi_bcast(function_type_ham,maxnum_funcvalues_ham*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_ham,maxnum_funcvalues_ham*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_ham,maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_ham,maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_ham,maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_ham,maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_ham,maxnum_funcvalues_ham*npairs,mpi_real8,0,mpi_comm_world,mpierror)
      endif
      end
