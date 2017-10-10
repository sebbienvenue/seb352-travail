!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - getallshortforces.f90
!!
      subroutine getshortforces(max_num_neighbors_short_atomic,num_atoms,&
        num_neighbors_short_atomic,neighboridx_short_atomic,zelem,symfunction,&
        dsfuncdxyz,deshortdsfunc,nnshortforce)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use nnshort_atomic
!!
      implicit none
!!
      integer max_num_neighbors_short_atomic                                 ! in
      integer num_atoms                                                      ! in
      integer num_neighbors_short_atomic(num_atoms)                          ! in
      integer neighboridx_short_atomic(num_atoms,0:max_num_neighbors_short_atomic) ! in
      integer zelem(max_num_atoms)                                           ! in
      integer i1,i2,i3,i4,i5                                                 ! internal
      integer icount                                                         ! internal
      integer itemp                                                          ! internal
!!
      real*8 dsfuncdxyz(maxnum_funcvalues_short_atomic,max_num_atoms,0:max_num_neighbors_short_atomic,3) ! in
      real*8 nnshortforce(3,max_num_atoms)                                   ! out 
      real*8 deshortdsfunc(max_num_atoms,maxnum_funcvalues_short_atomic)     ! out (for stress) 
      real*8 weights(maxnum_weights_short_atomic)                            ! internal
      real*8 symfunction(maxnum_funcvalues_short_atomic,max_num_atoms)       ! in
      real*8 symfunction_atom(maxnum_funcvalues_short_atomic)                ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                        ! internal 
      real*8 nodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic)  ! internal
      real*8 nodes_sum(maxnum_layers_short_atomic,maxnodes_short_atomic)     ! internal
      real*8 dnodes_values(maxnum_layers_short_atomic,maxnodes_short_atomic) ! internal
      real*8 tempderivative(maxnum_layers_short_atomic,maxnodes_short_atomic,maxnum_funcvalues_short_atomic) ! internal
!!   
!!===============================================
!! initializations
!!===============================================
      nnshortforce(:,:) = 0.0d0
      deshortdsfunc(:,:)= 0.0d0
!!
!!===============================================
!! get deshortdsfunc for each atom
!!===============================================
      do i1=1,num_atoms
!! don't calculate forces that we don't need
        if(lupdatebyelement.and.(zelem(i1).ne.elemupdate)) goto 98
!!
        symfunction_atom(:)=symfunction(:,i1)
        weights(:)=weights_short_atomic(:,elementindex(zelem(i1)))
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
        itemp = elementindex(zelem(i1))
!!
!!===============================================
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
!!===============================================
        call calconenn(1,maxnum_funcvalues_short_atomic,maxnodes_short_atomic,&
          maxnum_layers_short_atomic,num_layers_short_atomic(itemp),&
          maxnum_weights_short_atomic,nodes_short_atomic(0,itemp),&
          symfunction_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_short_atomic(1,1,itemp))
!!
!!===============================================
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
!!===============================================
!!
        call getdnodes_values(maxnum_layers_short_atomic,num_layers_short_atomic,maxnodes_short_atomic,&
          nodes_short_atomic,itemp,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_short_atomic)
!!
!!===============================================
!! calculate the full derivative deshortdsfunc of E_i with respect to G_j^i
!!===============================================
        tempderivative(:,:,:)=0.0d0
!!===============================================
!! for layer 1
!!===============================================
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_short_atomic(1) values for each input node
        do i2=1,nodes_short_atomic(1,itemp) ! over all nodes in layer 1 ("target layer") 
          do i3=1,num_funcvalues_short_atomic(itemp) ! over all nodes in previous layer
            icount=(i3-1)*nodes_short_atomic(1,itemp)+i2 ! set pointer in weights array, don't need windex_short_atomic for first weights
            tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
          enddo ! i3
        enddo ! i2
!!===============================================
!! for layers 2 and beyond (if present) 
!!===============================================
        if(num_layers_short_atomic(itemp).gt.1)then
          do i2=2,num_layers_short_atomic(itemp) ! over all hidden and output layers
            do i3=1,nodes_short_atomic(i2,itemp) ! over all nodes in the target layer
              do i5=1,num_funcvalues_short_atomic(itemp)
!! we have to sum over the nodes in the previous layer (i4)
                do i4=1,nodes_short_atomic(i2-1,itemp) ! sum over all nodes in previous layer
                  icount=windex_short_atomic(2*i2-1,itemp)&
                    +(i4-1)*nodes_short_atomic(i2,itemp)+i3-1 ! set pointer in weight array
                  tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5) + &
                    dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
                enddo ! i4
              enddo ! i5
            enddo ! i3
          enddo ! i2
        endif
!!
!!===============================================
!! get deshortdsfunc from tempderivate of output node
!!===============================================
        deshortdsfunc(i1,:)=tempderivative(num_layers_short_atomic(itemp),1,:)
!!
 98     continue
      enddo ! i1 ! num_atoms
!!
!!===========================================================
!! calculate the final forces 
!!===========================================================
      do i2=1,num_atoms ! over all atoms in structure
        do i1=0,num_neighbors_short_atomic(i2)  
!! don't calculate forces that we don't need
        if(lupdatebyelement.and.(zelem(neighboridx_short_atomic(i2,i1)).ne.elemupdate)) goto 99
          do i3=1,3 ! x,y,z
            do i4=1,num_funcvalues_short_atomic(elementindex(zelem(i2))) ! over all symmetry functions
              nnshortforce(i3,neighboridx_short_atomic(i2,i1))=nnshortforce(i3,neighboridx_short_atomic(i2,i1)) &
                -deshortdsfunc(i2,i4)*dsfuncdxyz(i4,i2,i1,i3)
            enddo ! i4
          enddo ! i3
 99     continue
        enddo ! i1
      enddo ! i2
!!
      return
      end
