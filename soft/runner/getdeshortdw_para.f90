!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - optimize_short_combined.f90
!!
!! called once for each point
!!
      subroutine getdeshortdw_para(&
           num_atoms,natoms,atomindex,&
           num_weights_short_atomic_free,&
           zelem,wconstraintidx,&
           symfunction,deshortdw)
!!
      use mpi_mod
      use fileunits
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer num_atoms                                                     ! in
      integer num_weights_short_atomic_free(nelem)                          ! in
      integer zelem(max_num_atoms)                                          ! in
      integer numelement(nelem)                                             ! internal
      integer i1,i2                                                         ! internal
      integer natoms                                                        ! in
      integer atomindex(natoms)                                             ! in
      integer wconstraintidx(maxnum_weights_short_atomic,nelem)             ! in
      integer itemp                                                         ! internal
!!
      real*8 weights(maxnum_weights_short_atomic)                           ! internal
!! CAUTION: just one output node is assumed here
      real*8 deshortdw(maxnum_weights_short_atomic,1,nelem) ! out 
!! CAUTION: just one output node is assumed here
      real*8 dedw(maxnum_weights_short_atomic,1)            ! internal
      real*8 symfunction(maxnum_funcvalues_short_atomic,natoms)                          ! in
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)                            ! internal
!!
!!
!! initializations
      deshortdw(:,:,:)=0.0d0
      numelement(:)=0 ! counts the number of atoms for each element
!!
!! loop over all atoms of the structure
      do i1=1,natoms
        itemp=elementindex(zelem(atomindex(i1)))
        dedw(:,:)=0.0d0
        weights(:)=weights_short_atomic(:,itemp)
        symfunction_atom(:)=symfunction(:,i1)
        numelement(itemp) = numelement(itemp)+1
!! calculate the derivative dedw for one specific atom/element/output node
!! even for weight constraints we calculate all dedw for simplicity
        call getonededw(1,&
          maxnum_funcvalues_short_atomic,maxnum_weights_short_atomic,&
          maxnodes_short_atomic,maxnum_layers_short_atomic,num_layers_short_atomic(itemp),&
          windex_short_atomic(1,itemp),nodes_short_atomic(0,itemp),&
          symfunction_atom,weights,dedw,&
          actfunc_short_atomic(1,1,itemp))
!! sum up the total derivative array for each element
        do i2=1,num_weights_short_atomic_free(itemp)
          deshortdw(i2,:,itemp)&
            =deshortdw(i2,:,itemp)&
            +dedw(wconstraintidx(i2,itemp),:)
        enddo ! i2
      enddo ! i1
!!
!! combine the deshortdw arrays of all processes
!! CAUTION: just one output node is assumed here
      call mpi_allreduce(mpi_in_place,deshortdw,&
        nelem*maxnum_weights_short_atomic*1,& 
        mpi_real8,mpi_sum,mpi_comm_world,mpierror)
!! 
!! combine the numelement arrays of all processes
      call mpi_allreduce(mpi_in_place,numelement,&
        nelem,mpi_integer,mpi_sum,mpi_comm_world,mpierror)
!!
!! normalization of the derivatives
      do i1=1,nelem
        if(numelement(i1).gt.0)then
!!          deshortdw(:,:,i1)=deshortdw(:,:,i1)/dble(numelement(i1))
          deshortdw(:,:,i1)=deshortdw(:,:,i1)/dble(num_atoms)
        endif
      enddo
!!
      return
      end
