!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! Purpose: calculate dchargedsfunc from natoms atoms 

!! called by: 
!!
      subroutine getdchargedsfunc_para(natoms,atomindex,&
        zelem,symfunctione,dchargedsfunc)
!!
      use fileunits
      use fittingoptions
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer natoms                                                    ! in
      integer atomindex(natoms)                                         ! in
      integer ielem                                                     ! internal
      integer zelem(max_num_atoms)                                      ! in
      integer i1,i2,i3,i4,i5                                            ! internal
      integer icount                                                    ! internal
!!
      real*8 dchargedsfunc(natoms,maxnum_funcvalues_elec)                   ! out  
      real*8 weights(maxnum_weights_elec)                               ! internal
      real*8 symfunctione(maxnum_funcvalues_elec,natoms)                    ! in
      real*8 symfunctione_atom(maxnum_funcvalues_elec)                      ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                   ! internal 
      real*8 nodes_values(maxnum_layers_elec,maxnodes_elec)            ! internal
      real*8 nodes_sum(maxnum_layers_elec,maxnodes_elec)               ! internal
      real*8 dnodes_values(maxnum_layers_elec,maxnodes_elec)           ! internal
      real*8 tempderivative(maxnum_layers_elec,maxnodes_elec,maxnum_funcvalues_elec) ! internal
!!   
!!
!! initialization
      dchargedsfunc(:,:)= 0.0d0
!!
!! get dchargedsfunc for each atom
!!--------------------------------
      do i1=1,natoms
!! don't calculate forces that we don't need
        if(lupdatebyelement.and.(zelem(atomindex(i1)).ne.elemupdate)) goto 98
!!
        ielem=elementindex(zelem(atomindex(i1)))
        symfunctione_atom(:)=symfunctione(:,i1)
        weights(:)=weights_elec(:,ielem)
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
!!
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
        call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
          maxnum_layers_elec,num_layers_elec(ielem),&
          maxnum_weights_elec,nodes_elec(0,ielem),&
          symfunctione_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_elec(1,1,ielem))
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
!!
        call getdnodes_values(maxnum_layers_elec,num_layers_elec,maxnodes_elec,&
          nodes_elec,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_elec)
!!
!! calculate the full derivative of E_i with respect to G_j^i
!!-----------------------------------------------------------
!!
        tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_ewald(1) values for each input node
        do i2=1,nodes_elec(1,ielem) ! over all nodes in layer 1 ("target layer") 
          do i3=1,num_funcvalues_elec(ielem) ! over all nodes in previous layer
            icount=(i3-1)*nodes_elec(1,ielem)+i2 ! set pointer in weights array, don't need windexe for first weights
            tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
          enddo ! i3
        enddo ! i2
!! for layers 2 and beyond (if present) 
        if(num_layers_elec(ielem).gt.1)then
          do i2=2,num_layers_elec(ielem) ! over all hidden and output layers
            do i3=1,nodes_elec(i2,ielem) ! over all nodes in the target layer
              do i5=1,num_funcvalues_elec(ielem)
!! we have to sum over the nodes in the previous layer (i4)
                do i4=1,nodes_elec(i2-1,ielem) ! sum over all nodes in previous layer
                  icount=windex_elec(2*i2-1,ielem)&
                    +(i4-1)*nodes_elec(i2,ielem)+i3-1 ! set pointer in weight array
                  tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5) + &
                    dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
                enddo ! i4
              enddo ! i5
            enddo ! i3
          enddo ! i2
        endif
!!
        do i2=1,maxnum_funcvalues_elec
          dchargedsfunc(i1,i2)=tempderivative(num_layers_elec(ielem),1,i2)
        enddo
!!
 98     continue
      enddo ! i1 ! natoms
!!
!!
      return
      end
