!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by: 
!! - getallelectrostatic.f90
!!
!! this subroutine is for the non-periodic case
!!
      subroutine getcoulombforces(max_num_neighbors_elec,&
         num_neighbors_elec,neighboridx_elec,num_atoms,zelem,&
         dsfuncdxyze,symfunctione,&
         nnatomcharge,xyzstruct,nnelecforce)
!!
      use fileunits
      use globaloptions
      use nnewald
!!
      implicit none
!!
      integer num_atoms                                                    ! in
      integer ielem                                                        ! internal
      integer zelem(max_num_atoms)                                         ! in
      integer max_num_neighbors_elec                                       ! in
      integer num_neighbors_elec(num_atoms)                                ! in
      integer neighboridx_elec(num_atoms,0:max_num_neighbors_elec)         ! in
      integer i1,i2,i3,i4,i5                                               ! internal
      integer icount                                                       ! internal
!!
      real*8 nnatomcharge(max_num_atoms)                                   ! in
      real*8 xyzstruct(3,max_num_atoms)                                    ! in
      real*8 nnelecforce(3,max_num_atoms)                                  ! out
      real*8 dsfuncdxyze(maxnum_funcvalues_elec,max_num_atoms,0:max_num_neighbors_elec,3) ! in
      real*8 dchargedsfunc(max_num_atoms,maxnum_funcvalues_elec)           ! internal
      real*8 dchargedxyz(max_num_atoms,0:max_num_neighbors_elec,3)         ! internal
      real*8 weights(maxnum_weights_elec)                                  ! internal
      real*8 symfunctione(maxnum_funcvalues_elec,max_num_atoms)            ! in
      real*8 symfunctione_atom(maxnum_funcvalues_elec)                     ! internal
      real*8 tempderivative(maxnum_layers_elec,maxnodes_elec,maxnum_funcvalues_elec) ! internal
      real*8 nodes_values(maxnum_layers_elec,maxnodes_elec)               ! internal
      real*8 nodes_sum(maxnum_layers_elec,maxnodes_elec)                  ! internal
      real*8 dnodes_values(maxnum_layers_elec,maxnodes_elec)              ! internal
!! CAUTION: just one output node is assumed here
      real*8 nnoutput                                                      ! internal 
      real*8 distance                                                      ! internal
      real*8 invrij2                                                       ! internal
      real*8 drijdxyz(max_num_atoms,3)                                     ! internal
      real*8 alphagaussian
      real*8 fscreen                                                       ! internal
      real*8 fscreenderiv                                                  ! internal
!!
!!
!! initializations
!!----------------
      dchargedxyz(:,:,:)= 0.0d0
      alphagaussian     = 0.5d0
      nnelecforce(:,:)  = 0.0d0
!! 
!! get dchargesfunc for each atom
!!--------------------------------
      do i1=1,num_atoms
        symfunctione_atom(:)=symfunctione(:,i1)
        weights(:)=weights_elec(:,elementindex(zelem(i1)))
        nodes_sum(:,:)    =0.0d0 ! initialization
        nodes_values(:,:) =0.0d0 ! initialization
        dnodes_values(:,:)=0.0d0 ! initialization
!!
!! calculation of the values on all nodes (nodes_values) in the NN (needed below for the derivatives)
!! This is done on-the-fly here although it is also calculated before
!! in order to save memory
        call calconenn(1,maxnum_funcvalues_elec,maxnodes_elec,&
          maxnum_layers_elec,num_layers_elec(elementindex(zelem(i1))),&
          maxnum_weights_elec,&
          nodes_elec(0,elementindex(zelem(i1))),&
          symfunctione_atom,weights,nodes_values,nodes_sum,&
          nnoutput,actfunc_elec(1,1,elementindex(zelem(i1))))
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
        ielem=elementindex(zelem(i1))
!!
        call getdnodes_values(maxnum_layers_elec,num_layers_elec,maxnodes_elec,&
          nodes_elec,ielem,nelem,nodes_sum,nodes_values,dnodes_values,actfunc_elec)
!!
        tempderivative(:,:,:)=0.0d0
!! for layer 1
!! calculate \frac{\partial f^1(x_j^1)}{\partial G_i}*a_{ij}^{01}
!! This is the derivative of the values of the nodes in the first layer with respect to G_i
!! => nodes_ewald(1) values for each input node
        do i2=1,nodes_elec(1,elementindex(zelem(i1))) ! over all nodes in layer 1 ("target layer") 
          do i3=1,num_funcvalues_elec(elementindex(zelem(i1))) ! over all nodes in previous layer
            icount=(i3-1)*nodes_elec(1,elementindex(zelem(i1)))+i2 ! set pointer in weights array, don't need windex_short_atomic for first weights
            tempderivative(1,i2,i3)=dnodes_values(1,i2)*weights(icount)
          enddo ! i3
        enddo ! i2
!! for layers 2 and beyond (if present) 
       if(num_layers_elec(elementindex(zelem(i1))).gt.1)then
         do i2=2,num_layers_elec(elementindex(zelem(i1))) ! over all hidden and output layers
           do i3=1,nodes_elec(i2,elementindex(zelem(i1))) ! over all nodes in the target layer
             do i5=1,num_funcvalues_elec(elementindex(zelem(i1)))
!! we have to sum over the nodes in the previous layer (i4)
               do i4=1,nodes_elec(i2-1,elementindex(zelem(i1))) ! sum over all nodes in previous layer
                 icount=windex_elec(2*i2-1,elementindex(zelem(i1)))&
                  +(i4-1)*nodes_elec(i2,elementindex(zelem(i1)))+i3-1 ! set pointer in weight array
                 tempderivative(i2,i3,i5)=tempderivative(i2,i3,i5) + &
                 dnodes_values(i2,i3)*weights(icount)*tempderivative(i2-1,i4,i5)
               enddo ! i4
             enddo ! i5
           enddo ! i3
         enddo ! i2
       endif
!!
       dchargedsfunc(i1,:)=tempderivative(num_layers_elec(elementindex(zelem(i1))),1,:)
!!
      enddo ! i1
!!
!! calculation of dchargedxyz 
!!---------------------------
!! summation over all symmetry functions for each atom  
      do i1=1,num_atoms
        do i2=1,num_funcvalues_elec(elementindex(zelem(i1)))
          do i3=0,num_neighbors_elec(i1)
            do i4=1,3
          dchargedxyz(i1,i3,i4)=dchargedxyz(i1,i3,i4) &
            + dchargedsfunc(i1,i2)*dsfuncdxyze(i2,i1,i3,i4)        ! CHECK!
            enddo ! i4
          enddo ! i3
        enddo ! i2
      enddo ! i1
!!
!! calculation of dE/dx
!!----------------------
!!
      if(num_atoms.gt.1) then
        do i3=0,num_neighbors_elec(i1) !! FIXME
          do i4=1,3
!!
!! caution: do we need double counting for the forces???
            do i1=1,num_atoms
              do i2=i1+1,num_atoms
                distance=(xyzstruct(1,i1)-xyzstruct(1,i2))**2 + &
                  (xyzstruct(2,i1)-xyzstruct(2,i2))**2 + &
                  (xyzstruct(3,i1)-xyzstruct(3,i2))**2
                invrij2 =1.d0/distance
                distance=dsqrt(distance)
!! calculation of \frac{\partial r_{ij}}{\partial \alpha}
                drijdxyz(i3,i4)=0.0d0 ! initialization
                if(i1.eq.i3)then
                  drijdxyz(i3,i4)=(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
                elseif(i2.eq.i3)then
                  drijdxyz(i3,i4)=-1.d0*(xyzstruct(i4,i1)-xyzstruct(i4,i2))/distance
                endif
!!
!! final force calculation
!!------------------------
                if(lscreen) then
                  call getscreenfunctionforelectrostatics(&
                         distance,fscreen,&
                         fscreenderiv,drijdxyz(i3,i4))
!!
                  nnelecforce(i4,i3)=nnelecforce(i4,i3) &
                   - invrij2*((dchargedxyz(i1,i3,i4)*nnatomcharge(i2) &
                              + nnatomcharge(i1)*dchargedxyz(i2,neighboridx_elec(i2,i3),i4))*&
                              distance*fscreen - nnatomcharge(i1)* &
                              nnatomcharge(i2)*drijdxyz(i3,i4)* &
                              fscreen+fscreenderiv*nnatomcharge(i1)* &
                              nnatomcharge(i2)*distance)
                else
                  nnelecforce(i4,i3)=nnelecforce(i4,i3) &
                    - invrij2*((dchargedxyz(i1,neighboridx_elec(i1,i3),i4)*nnatomcharge(i2) &
                        + nnatomcharge(i1)*dchargedxyz(i2,neighboridx_elec(i2,i3),i4))*distance&
                    -nnatomcharge(i1)*nnatomcharge(i2)*drijdxyz(i3,i4))
                endif
!!
              enddo ! i2
            enddo ! i1
          enddo ! i4
        enddo ! i3
!! 
      endif ! num_atoms.gt.1
!!
      return
      end
