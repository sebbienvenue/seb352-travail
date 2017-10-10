!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!!
      subroutine getatom(ndim,npoints,&
        zelem_local,num_atoms_local,&
        symfunction_local,nneshort_local,&
        nnatomcharge_local,nnchargesum_local)
!!
      use fileunits
      use nnflags
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer ndim                                                    ! in
      integer npoints                                                 ! in
      integer zelem_local(ndim,max_num_atoms)                         ! in
      integer zelem(max_num_atoms)                                    ! internal
      integer num_atoms_local(ndim)                                   ! in
      integer nndim                                                   ! internal
!!
      integer i1,i2
!!
      real*8 symfunction_local(maxnum_funcvalues_short_atomic,max_num_atoms,ndim)  ! in
      real*8 nneshort_local(ndim)                                     ! out
      real*8 nnatomcharge_local(ndim,max_num_atoms)                   ! out
      real*8 nnchargesum_local(ndim)                                  ! out
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)         ! internal
      real*8 weights(maxnum_weights_short_atomic)                     ! internal
      real*8, dimension(:)  , allocatable :: nnoutput                 ! internal
      real*8 nodes_values_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic) ! just dummy in this routine
      real*8 nodes_sum_dummy(maxnum_layers_short_atomic,maxnodes_short_atomic)    ! just dummy in this routine
      real*8 nneshort                                                 ! internal
!!
!!============================================
!! determine number of output nodes
!!============================================
      if(lelec.and.(nn_type_elec.eq.2))then 
        nndim=2
      else
        nndim=1
      endif
      allocate(nnoutput(nndim))
!!
!!============================================
!! loop over all structures 
!!============================================
      do i1=1,npoints
        zelem(:)=zelem_local(i1,:)
!!
!!============================================
!! loop over all atoms of a given structure 
!!============================================
        nnchargesum_local(i1)=0.0d0
        nneshort=0.0d0
        do i2=1,num_atoms_local(i1)
          symfunction_atom(:)=symfunction_local(:,i2,i1)
          weights(:)=weights_short_atomic(:,elementindex(zelem(i2)))
!!
!! calculate nodes_values and nnoutput for atom i2
          call calconenn(nndim,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
            maxnum_layers_short_atomic,num_layers_short_atomic(elementindex(zelem(i2))),&
            maxnum_weights_short_atomic,nodes_short_atomic(0,elementindex(zelem(i2))),&
            symfunction_atom,weights,nodes_values_dummy,nodes_sum_dummy,&
            nnoutput,actfunc_short_atomic(1,1,elementindex(zelem(i2))))
!!
          nneshort=nneshort+nnoutput(1)
          if(lelec.and.(nn_type_elec.eq.2))then 
            nnatomcharge_local(i1,i2)=nnoutput(2)
            nnchargesum_local(i1)=nnchargesum_local(i1)+nnoutput(2) 
          endif
!!
        enddo ! i2
!!
!! normalize to energy per atom
        nneshort_local(i1)=nneshort/dble(num_atoms_local(i1))
!!
      enddo ! i1
!!
      deallocate(nnoutput)
!!
      return
      end
