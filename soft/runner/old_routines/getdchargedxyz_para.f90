!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate dchargedxyz

!! called by: - prediction.f90 
!!
      subroutine getdchargedxyz_para(max_num_neighbors,&
         num_neighbors,natoms,atomindex,zelem,&
         dsfuncdxyze,&
         dchargedxyz,symfunctione)
!!
      use fileunits
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer ielem                                                     ! in
      integer zelem(max_num_atoms)                                      ! in
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer max_num_neighbors                                         ! in
      integer num_neighbors(natoms)                                     ! in
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                                    ! internal
!!
      real*8 dsfuncdxyze(maxnum_funcvaluese,natoms,0:max_num_neighbors,3) ! in
      real*8 dchargedsfunc(natoms,maxnum_funcvaluese)                   ! internal 
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors,3)           ! out 
      real*8 symfunctione(maxnum_funcvaluese,natoms)                    ! in
      real*8 symfunctione_atom(maxnum_funcvaluese)                      ! internal
      real*8 weights(maxnum_weightsewald)                               ! internal
      real*8 nodes_values(maxnum_layersewald,maxnodes_ewald)            ! internal
      real*8 nodes_sum(maxnum_layersewald,maxnodes_ewald)               ! internal
      real*8 dnodes_values(maxnum_layersewald,maxnodes_ewald)           ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                   ! internal 
      real*8 tempderivative(maxnum_layersewald,maxnodes_ewald,maxnum_funcvaluese) ! internal
      real*8 alphagaussian                                              ! internal
!!
!!
!! initializations
!!----------------
      alphagaussian     =0.5d0
!! 
!! get dchargesfunc for each atom
!!--------------------------------
      do i1=1,natoms
        symfunctione_atom(:)=symfunctione(:,i1)
        weights(:)=weights_ewald(:,elementindex(zelem(atomindex(i1))))
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
!!
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
        call calconenn(1,maxnum_funcvaluese,maxnodes_ewald,&
          maxnum_layersewald,num_layersewald(elementindex(zelem(atomindex(i1)))),&
          maxnum_weightsewald,nodes_ewald(0,elementindex(zelem(atomindex(i1)))),&
          symfunctione_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_ewald(1,1,elementindex(zelem(atomindex(i1)))))
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
        ielem=elementindex(zelem(atomindex(i1)))
!!
        call getdnodes_values(maxnum_layersewald,num_layersewald,maxnodes_ewald,&
          nodes_ewald,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_ewald)
!!
       tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_ewald(1) values for each input node
       do i2=1,nodes_ewald(1,elementindex(zelem(atomindex(i1)))) ! over all nodes in layer 1 ("target layer") 
         do i3=1,num_funcvaluese(elementindex(zelem(atomindex(i1)))) ! over all nodes in previous layer
           icount=(i3-1)*nodes_ewald(1,elementindex(zelem(atomindex(i1))))+i2 ! set pointer in weights array, don't need windex for first weights
           tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
         enddo ! i3
       enddo ! i2
!! for layers 2 and beyond (if present) 
      if(num_layersewald(elementindex(zelem(atomindex(i1)))).gt.1)then
        do i2=2,num_layersewald(elementindex(zelem(atomindex(i1)))) ! over all hidden and output layers
          do i3=1,nodes_ewald(i2,elementindex(zelem(atomindex(i1)))) ! over all nodes in the target layer
            do i5=1,num_funcvaluese(elementindex(zelem(atomindex(i1))))
!! we have to sum over the nodes in the previous layer (i4)
              do i4=1,nodes_ewald(i2-1,elementindex(zelem(atomindex(i1)))) ! sum over all nodes in previous layer
                icount=windexe(2*i2-1,elementindex(zelem(atomindex(i1))))&
                  +(i4-1)*nodes_ewald(i2,elementindex(zelem(atomindex(i1))))+i3-1 ! set pointer in weight array
                tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5) + &
                dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
              enddo ! i4
            enddo ! i5
          enddo ! i3
        enddo ! i2
      endif
!!
        dchargedsfunc(i1,:)=tempderivative(num_layersewald(elementindex(zelem(atomindex(i1)))),1,:)
!!
      enddo ! i1
!!
!! calculation of dchargedxyz 
!!---------------------------
!! summation over all symmetry functions for each atom  
      do i4=1,3
        do i1=1,natoms
          do i2=1,num_funcvaluese(elementindex(zelem(atomindex(i1))))
            do i3=0,num_neighbors(i1)
              dchargedxyz(atomindex(i1),i3,i4)=dchargedxyz(atomindex(i1),i3,i4) &
                + dchargedsfunc(i1,i2)*dsfuncdxyze(i2,i1,i3,i4)        ! CHECK!
            enddo ! i3
          enddo ! i2
        enddo ! i1
      enddo ! i4
!!
      return
      end
