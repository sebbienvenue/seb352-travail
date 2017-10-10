!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - initnn.f90
!!
      subroutine allocatesymfunctions()
!!
      use globaloptions
      use symfunctions
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        allocate(function_type_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        function_type_short_atomic(:,:)=0
        allocate(symelement_short_atomic(maxnum_funcvalues_short_atomic,2,nelem))
        symelement_short_atomic(:,:,:)=0
        allocate(funccutoff_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        funccutoff_short_atomic(:,:)=0.0d0
        allocate(eta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        eta_short_atomic(:,:)=0.0d0
        allocate(zeta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        zeta_short_atomic(:,:)=0.0d0
        allocate(lambda_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        lambda_short_atomic(:,:)=0.0d0
        allocate(rshift_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        rshift_short_atomic(:,:)=0.0d0
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        allocate(function_type_short_pair(maxnum_funcvalues_short_pair,npairs))
        function_type_short_pair(:,:)=0
        allocate(symelement_short_pair(maxnum_funcvalues_short_pair,2,npairs))
        symelement_short_pair(:,:,:)=0
        allocate(funccutoff_short_pair(maxnum_funcvalues_short_pair,npairs))
        funccutoff_short_pair(:,:)=0.0d0
        allocate(eta_short_pair(maxnum_funcvalues_short_pair,npairs))
        eta_short_pair(:,:)=0.0d0
        allocate(zeta_short_pair(maxnum_funcvalues_short_pair,npairs))
        zeta_short_pair(:,:)=0.0d0
        allocate(lambda_short_pair(maxnum_funcvalues_short_pair,npairs))
        lambda_short_pair(:,:)=0.0d0
        allocate(rshift_short_pair(maxnum_funcvalues_short_pair,npairs))
        rshift_short_pair(:,:)=0.0d0
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        allocate(function_type_elec(maxnum_funcvalues_elec,nelem))
        function_type_elec(:,:)=0
        allocate(symelement_elec(maxnum_funcvalues_elec,2,nelem))
        symelement_elec(:,:,:)=0
        allocate(funccutoff_elec(maxnum_funcvalues_elec,nelem))
        funccutoff_elec(:,:)=0.0d0
        allocate(eta_elec(maxnum_funcvalues_elec,nelem))
        eta_elec(:,:)=0.0d0
        allocate(zeta_elec(maxnum_funcvalues_elec,nelem))
        zeta_elec(:,:)=0.0d0
        allocate(lambda_elec(maxnum_funcvalues_elec,nelem))
        lambda_elec(:,:)=0.0d0
        allocate(rshift_elec(maxnum_funcvalues_elec,nelem))
        rshift_elec(:,:)=0.0d0
      endif
!!
      if(lnntb)then
        allocate(function_type_ham(maxnum_funcvalues_ham,npairs))
        function_type_ham(:,:)=0
        allocate(symelement_ham(maxnum_funcvalues_ham,2,npairs))
        symelement_ham(:,:,:)=0
        allocate(funccutoff_ham(maxnum_funcvalues_ham,npairs))
        funccutoff_ham(:,:)=0.0d0
        allocate(eta_ham(maxnum_funcvalues_ham,npairs))
        eta_ham(:,:)=0.0d0
        allocate(zeta_ham(maxnum_funcvalues_ham,npairs))
        zeta_ham(:,:)=0.0d0
        allocate(lambda_ham(maxnum_funcvalues_ham,npairs))
        lambda_ham(:,:)=0.0d0
        allocate(rshift_ham(maxnum_funcvalues_ham,npairs))
        rshift_ham(:,:)=0.0d0
      endif
      end
