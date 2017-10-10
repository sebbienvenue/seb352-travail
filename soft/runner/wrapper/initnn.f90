!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: main.f90
!!
      subroutine initnn(iseed, ielem)
!!
      use mpi_mod
      use fileunits 
      use fittingoptions
      use mode1options
      use predictionoptions
      use nnflags
      use globaloptions
      use symfunctions
      use timings
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use nnham
      use basismod
!!
      implicit none
!!
      integer ielem ! number of elements in input.data (just used for checking in initialization part)
      integer iseed ! seed for random numbers
!!
      logical lelement(102)                   ! internal for readinput
!!
      call zerotime(dayinitnn,timeinitnnstart,timeinitnnend)
      call abstime(timeinitnnstart,dayinitnn)
!!
      listdim    =100000 ! preliminary value, is overwritten by better estimate below
!!
!! analyze input.nn file for the first time to get dimensions of arrays
      if(mpirank.eq.0)then
        call initialization(ielem,lelement)
      endif ! mpirank.eq.0
!!
!! distribute contents of nnflags module 
      call distribute_nnflags()
!!
!! distribute information from initialization
      call mpi_bcast(totnum_structures,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(max_num_atoms,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(max_num_pairs,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(ielem,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nelem,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(npairs,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lelement,102,mpi_logical,0,mpi_comm_world,mpierror)
!!
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(maxnum_layers_short_atomic,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnum_funcvalues_short_atomic,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnodes_short_atomic,1,mpi_integer,0,mpi_comm_world,mpierror)
        if(mode.eq.2)allocate (kalmanlambda(nelem))
        allocate (num_funcvalues_short_atomic(nelem))
        num_funcvalues_short_atomic(:)=0
        allocate (windex_short_atomic(2*maxnum_layers_short_atomic,nelem))
        allocate (num_layers_short_atomic(nelem))
        num_layers_short_atomic(:)=maxnum_layers_short_atomic
        allocate (actfunc_short_atomic(maxnodes_short_atomic,maxnum_layers_short_atomic,nelem))
        allocate (nodes_short_atomic(0:maxnum_layers_short_atomic,nelem))
        nodes_short_atomic(:,:)=0
        allocate (num_weights_short_atomic(nelem))
        num_weights_short_atomic(:)=0
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(maxnum_layers_short_pair,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnum_funcvalues_short_pair,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnodes_short_pair,1,mpi_integer,0,mpi_comm_world,mpierror)
        if(mode.eq.2)allocate (kalmanlambdap(npairs))
        allocate (num_funcvalues_short_pair(npairs))
        num_funcvalues_short_pair(:)=0
        allocate (windex_short_pair(2*maxnum_layers_short_pair,npairs))
        allocate (num_layers_short_pair(npairs))
        num_layers_short_pair(:)=maxnum_layers_short_pair
        allocate (actfunc_short_pair(maxnodes_short_pair,maxnum_layers_short_pair,npairs))
        allocate (nodes_short_pair(0:maxnum_layers_short_pair,npairs))
        nodes_short_pair(:,:)=0
        allocate (num_weights_short_pair(npairs))
        num_weights_short_pair(:)=0
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(maxnum_layers_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnum_funcvalues_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnodes_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
        if(mode.eq.2)allocate (kalmanlambdae(nelem))
        allocate (num_funcvalues_elec(nelem))
        num_funcvalues_elec(:)=0
        allocate (windex_elec(2*maxnum_layers_elec,nelem))
        allocate (num_layers_elec(nelem))
        num_layers_elec(:)=maxnum_layers_elec
        allocate (actfunc_elec(maxnodes_elec,maxnum_layers_elec,nelem))
        allocate (nodes_elec(0:maxnum_layers_elec,nelem))
        nodes_elec(:,:)=0
        allocate (num_weights_elec(nelem))
        num_weights_elec(:)=0
      endif
      if(lnntb)then
        call mpi_bcast(maxnum_layers_ham,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnum_funcvalues_ham,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(maxnodes_ham,1,mpi_integer,0,mpi_comm_world,mpierror)
!! todo:        if(mode.eq.2)allocate (kalmanlambdae(nelem))
        allocate (num_funcvalues_ham(npairs))
        num_funcvalues_ham(:)=0
        allocate (windex_ham(2*maxnum_layers_ham,npairs))
        allocate (num_layers_ham(npairs))
        num_layers_ham(:)=maxnum_layers_ham
        allocate (actfunc_ham(maxnodes_ham,maxnum_layers_ham,npairs))
        allocate (nodes_ham(0:maxnum_layers_ham,npairs))
        nodes_ham(:,:)=0
        allocate (num_weights_ham(npairs))
        num_weights_ham(:)=0
      endif
!!'
!! allocate globaloptions
      allocate (fixedcharge(nelem))
      fixedcharge(:)=0.0d0
      allocate (nucelem(nelem))
      allocate (element(nelem))
      allocate (atomrefenergies(nelem))
      allocate (elempair(npairs,2))
      elempair(:,:)=0
!!
!! allocate arrays for symmetry functions
      call allocatesymfunctions()
!!
!!======================================================================
!! read input.nn file for the second time, now all keywords are read
!!======================================================================
      call zerotime(dayreadinput,timereadinputstart,timereadinputend) 
      call abstime(timereadinputstart,dayreadinput)
      if(mpirank.eq.0)then
        call readinput(ielem,iseed,lelement)
      endif ! mpirank.eq.0
      call abstime(timereadinputend,dayreadinput)
!!
!! get optimum value for dimension parameter listdim
      call getlistdim()
!!
!! distribute mode 1 options
      if(mode.eq.1)then
        call mpi_bcast(splitthres,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(fitethres,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lfitethres,1,mpi_logical,0,mpi_comm_world,mpierror)
      endif
!!
!! distribute mode 2 options 
      if(mode.eq.2)then
        call distribute_fittingoptions()
      endif
!!
!! distribute mode 3 options
      if(mode.eq.3)then
        call distribute_predictionoptions()
      endif
!!
!! distribute symmetry function parameters
      call distribute_symfunctions()
!!
!! distribute remaining parameters of globaloptions module
      call distribute_globaloptions()
!!
!! distribute NN specific parameters and arrays
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(num_layers_short_atomic,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(nodes_short_atomic,(maxnum_layers_short_atomic+1)*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(windex_short_atomic,2*maxnum_layers_short_atomic*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_weights_short_atomic,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_funcvalues_short_atomic,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(actfunc_short_atomic,&
          maxnodes_short_atomic*maxnum_layers_short_atomic*nelem,mpi_character,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmin_short_atomic,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmax_short_atomic,1,mpi_real8,0,mpi_comm_world,mpierror)
        allocate (weights_short_atomic(maxnum_weights_short_atomic,nelem))
        weights_short_atomic(:,:)=0.0d0
        allocate (symfunction_short_atomic_list(maxnum_funcvalues_short_atomic,max_num_atoms,nblock))
        symfunction_short_atomic_list(:,:,:)=0.0d0
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(num_layers_short_pair,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(nodes_short_pair,(maxnum_layers_short_pair+1)*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(windex_short_pair,2*maxnum_layers_short_pair*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_weights_short_pair,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_funcvalues_short_pair,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(actfunc_short_pair,&
          maxnodes_short_pair*maxnum_layers_short_pair*npairs,mpi_character,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmin_short_pair,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmax_short_pair,1,mpi_real8,0,mpi_comm_world,mpierror)
        allocate (weights_short_pair(maxnum_weights_short_pair,npairs))
        weights_short_pair(:,:)=0.0d0
        allocate (symfunction_short_pair_list(maxnum_funcvalues_short_pair,max_num_pairs,nblock))
        symfunction_short_pair_list(:,:,:)=0.0d0
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(num_layers_elec,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(nodes_elec,(maxnum_layers_elec+1)*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(windex_elec,2*maxnum_layers_elec*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_weights_elec,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_funcvalues_elec,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(actfunc_elec,&
          maxnodes_elec*maxnum_layers_elec*nelem,mpi_character,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmin_elec,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmax_elec,1,mpi_real8,0,mpi_comm_world,mpierror)
        allocate (weights_elec(maxnum_weights_elec,nelem))
        weights_elec(:,:)=0.0d0
        allocate (symfunction_elec_list(maxnum_funcvalues_elec,max_num_atoms,nblock))
        symfunction_elec_list(:,:,:)=0.0d0
      endif
!!
      if(lnntb)then
        call mpi_bcast(num_layers_ham,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(nodes_ham,(maxnum_layers_ham+1)*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(windex_ham,2*maxnum_layers_ham*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_weights_ham,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_funcvalues_ham,npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(actfunc_ham,&
          maxnodes_ham*maxnum_layers_ham*npairs,mpi_character,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmin_ham,1,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(scmax_ham,1,mpi_real8,0,mpi_comm_world,mpierror)
        allocate (weights_ham(maxnum_weights_ham,npairs))
        weights_ham(:,:)=0.0d0
        allocate (symfunction_ham_list(maxnum_funcvalues_ham,max_num_pairs,nblock))
        symfunction_ham_list(:,:,:)=0.0d0
      endif
!!
!! distribute the rest
      call mpi_bcast(ielem,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(iseed,1,mpi_real8,0,mpi_comm_world,mpierror)
      if(nn_type_nntb.eq.1)then
        call mpi_bcast(maxnum_basis,1,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(num_basis,nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(basis,nelem*maxnum_basis*3,mpi_integer,0,mpi_comm_world,mpierror)
      endif
!!
      call abstime(timeinitnnend,dayinitnn)
!!
      return
      end
