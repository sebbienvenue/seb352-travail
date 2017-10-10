!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! called once for each point
!!
      subroutine getdnnoutdw(&
           num_weights_short_atomic_free,&
           zelem,num_atoms,wconstraintidx,&
           symfunction,deshortdw,dchargedw)
!!
      use fileunits
      use nnflags
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer num_weights_short_atomic_free(nelem)                          ! in
      integer zelem(max_num_atoms)                                          ! in
      integer num_atoms                                                     ! in
      integer numelement(nelem)                                             ! internal
      integer i1,i2                                                         ! internal
      integer wconstraintidx(maxnum_weights_short_atomic,nelem)             ! in
      integer itemp                                                         ! internal
      integer nndim                                                         ! internal
!!
      real*8 weights(maxnum_weights_short_atomic)                           ! internal
      real*8 deshortdw(maxnum_weights_short_atomic,nelem)                   ! out 
      real*8 dchargedw(maxnum_weights_short_atomic,max_num_atoms)           ! out 
      real*8, dimension(:,:), allocatable :: dedw                           ! internal
      real*8 symfunction(maxnum_funcvalues_short_atomic,max_num_atoms)      ! in
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)               ! internal
!!
      nndim=1
      if(lelec.and.(nn_type_elec.eq.2))then
        nndim=2
      endif
      allocate(dedw(maxnum_weights_short_atomic,nndim))
!!
!! initializations
      deshortdw(:,:) = 0.0d0
      dchargedw(:,:) = 0.0d0
      numelement(:)  = 0 ! counts the number of atoms for each element
!!
!! loop over all atoms of the structure
      do i1=1,num_atoms
        itemp=elementindex(zelem(i1))
        dedw(:,:)=0.0d0
        weights(:)=weights_short_atomic(:,itemp)
        symfunction_atom(:)=symfunction(:,i1)
        numelement(itemp) = numelement(itemp)+1
!! calculate the derivative dedw for one specific atom/element/output node
!! even for weight constraints we calculate all dedw for simplicity
        call getonededw(nndim,&
          maxnum_funcvalues_short_atomic,maxnum_weights_short_atomic,&
          maxnodes_short_atomic,maxnum_layers_short_atomic,&
          num_layers_short_atomic(itemp),&
          windex_short_atomic(1,itemp),nodes_short_atomic(0,itemp),&
          symfunction_atom,weights,dedw,&
          actfunc_short_atomic(1,1,itemp))
!! sum up the total derivative array for each element
        do i2=1,num_weights_short_atomic_free(itemp)
          deshortdw(i2,itemp)&
            =deshortdw(i2,itemp)+dedw(wconstraintidx(i2,itemp),1)
          if(lelec.and.(nn_type_elec.eq.2))then
            dchargedw(i2,i1)&
              =dchargedw(i2,i1)+dedw(wconstraintidx(i2,itemp),2)
          endif
        enddo ! i2
      enddo ! i1
!!
!! normalization of the derivatives
      do i1=1,nelem
        if(numelement(i1).gt.0)then
          deshortdw(:,i1)=deshortdw(:,i1)/dble(num_atoms)
        endif
      enddo
!!
      deallocate(dedw)
!!
      return
      end
